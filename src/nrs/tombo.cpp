#include "nrs.hpp"

namespace tombo
{
occa::memory pressureSolve(nrs_t *nrs, double time, int stage)
{
  auto mesh = nrs->meshV;

  double flopCount = 0.0;
  platform->timer.tic("pressure rhs", 1);

  const auto o_stressTerm = [&]()
  {
    auto o_curl = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);

    nrs->curlKernel(mesh->Nelements, 1, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, nrs->o_Ue, o_curl);
    flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Np * mesh->Nq + 36 * mesh->Np);

    oogs::startFinish(o_curl, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

    platform->linAlg->axmyVector(mesh->Nlocal, nrs->fieldOffset, 0, 1.0, nrs->meshV->o_invLMM, o_curl);
    flopCount += mesh->Nlocal;

    auto o_stressTerm = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    nrs->curlKernel(mesh->Nelements, 1, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, o_curl, o_stressTerm);
    flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Np * mesh->Nq + 36 * mesh->Np);
 
    if (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) {
      nrs->pressureStressKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              nrs->o_mue,
                              nrs->o_Ue,
                              nrs->o_div,
                              o_stressTerm);
      flopCount += static_cast<double>(mesh->Nelements) * (18 * mesh->Nq * mesh->Np + 100 * mesh->Np);
    }
    return o_stressTerm;
  }();


  const auto o_gradDiv = [&]()
  { 
    auto o_gradDiv = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    nrs->gradientVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              nrs->o_div,
                              o_gradDiv);
    flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);
    return o_gradDiv;
  }();

  const auto o_lambda0 = [&]()
  {
    auto o_lambda0 = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
    platform->linAlg->adyz(mesh->Nlocal, 1.0, nrs->o_rho, o_lambda0);
    return o_lambda0;
  }();

  const auto o_rhs = [&]()
  { 
    auto o_rhs = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
 
    if (platform->options.compareArgs("PRESSURE VISCOUS TERMS", "TRUE")) {
      nrs->pressureRhsKernel(mesh->Nlocal,
                             nrs->fieldOffset,
                             nrs->o_mue,
                             o_lambda0,
                             nrs->o_BF,
                             o_stressTerm,
                             o_gradDiv,
                             o_rhs);
      flopCount += 12 * static_cast<double>(mesh->Nlocal);
    } else {
  
      auto o_BFx = nrs->o_BF.slice(0 * nrs->fieldOffset, mesh->Nlocal);
      auto o_BFy = nrs->o_BF.slice(1 * nrs->fieldOffset, mesh->Nlocal);
      auto o_BFz = nrs->o_BF.slice(2 * nrs->fieldOffset, mesh->Nlocal);
 
      platform->linAlg->axmy(mesh->Nlocal, 1.0, o_lambda0, o_BFx);
      platform->linAlg->axmy(mesh->Nlocal, 1.0, o_lambda0, o_BFy);
      platform->linAlg->axmy(mesh->Nlocal, 1.0, o_lambda0, o_BFz);

      o_rhs.copyFrom(nrs->o_BF);
    }

    oogs::startFinish(o_rhs, nrs->NVfields, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
    platform->linAlg->axmyVector(mesh->Nlocal, nrs->fieldOffset, 0, 1.0, nrs->meshV->o_invLMM, o_rhs);

    return o_rhs;
  }();

  const auto o_pRhs = [&]()
  { 
    auto o_pRhs = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
 
    nrs->wDivergenceVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, o_rhs, o_pRhs);
    flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);
 
    nrs->pressureAddQtlKernel(mesh->Nlocal, mesh->o_LMM, nrs->g0 * 1 / nrs->dt[0], nrs->o_div, o_pRhs);
    flopCount += 3 * mesh->Nlocal;
 
    nrs->divergenceSurfaceKernel(mesh->Nelements,
                                 mesh->o_sgeo,
                                 mesh->o_vmapM,
                                 nrs->o_EToB,
                                 nrs->g0 * 1 / nrs->dt[0],
                                 nrs->fieldOffset,
                                 o_rhs,
                                 nrs->o_U,
                                 o_pRhs);
    flopCount += 25 * static_cast<double>(mesh->Nelements) * mesh->Nq * mesh->Nq;

    return o_pRhs;
  }();

  platform->timer.toc("pressure rhs");
  platform->flopCounter->add("pressure RHS", flopCount);

  occa::memory o_p = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
  o_p.copyFrom(nrs->o_P);

  nrs->pSolver->solve(o_lambda0, o_NULL, o_pRhs, o_p);

  if (platform->verbose) {
    const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                 1,
                                                                 nrs->fieldOffset,
                                                                 mesh->ogs->o_invDegree,
                                                                 o_p,
                                                                 platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("p norm: %.15e\n", debugNorm);
    }
  }

  return o_p;
}

occa::memory velocitySolve(nrs_t *nrs, double time, int stage)
{
  auto mesh = nrs->meshV;

  double flopCount = 0.0;
  platform->timer.tic("velocity rhs", 1);

  const auto o_gradMueDiv = [&]()
  {
    dfloat scale = 1./3;
    if (platform->options.compareArgs("VELOCITY STRESSFORMULATION", "TRUE")) scale = -2*scale;

    auto o_mueDiv = platform->o_memPool.reserve<dfloat>(nrs->fieldOffset);
    platform->linAlg->axmyz(mesh->Nlocal,
                            scale,
                            nrs->o_mue,
                            nrs->o_div,
                            o_mueDiv);

    auto o_gradMueDiv = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    nrs->gradientVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              o_mueDiv,
                              o_gradMueDiv);
    flopCount += static_cast<double>(mesh->Nelements) * (6 * mesh->Np * mesh->Nq + 18 * mesh->Np);

    return o_gradMueDiv;
  }();

  const auto o_gradP = [&]()
  {
    occa::memory o_gradP = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    nrs->wgradientVolumeKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, nrs->fieldOffset, nrs->o_P, o_gradP);
    flopCount += static_cast<double>(mesh->Nelements) * 18 * (mesh->Np * mesh->Nq + mesh->Np);

    return o_gradP;
  }();

  const auto o_rhs = [&]()
  {
    auto o_rhs = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    nrs->velocityRhsKernel(mesh->Nlocal, nrs->fieldOffset, nrs->o_BF, o_gradMueDiv, o_gradP, o_rhs);
    flopCount += 9 * mesh->Nlocal;

    nrs->velocityNeumannBCKernel(mesh->Nelements,
                                 nrs->fieldOffset,
                                 mesh->o_sgeo,
                                 mesh->o_vmapM,
                                 mesh->o_EToB,
                                 nrs->o_EToB,
                                 time,
                                 mesh->o_x,
                                 mesh->o_y,
                                 mesh->o_z,
                                 nrs->o_rho,
                                 nrs->o_mue,
                                 nrs->o_usrwrk,
                                 nrs->o_Ue,
                                 o_rhs);
 
    flopCount += static_cast<double>(mesh->Nelements) * (3 * mesh->Np + 36 * mesh->Nq * mesh->Nq);

    return o_rhs;
  }();

  platform->timer.toc("velocity rhs");
  platform->flopCounter->add("velocity RHS", flopCount);

  const auto o_lambda0 = nrs->o_mue;
  const auto o_lambda1 = [&]()
  {
    auto o_lambda1 = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
    if (nrs->userVelocityImplicitLinearTerm) {
      auto o_implicitLT = nrs->userVelocityImplicitLinearTerm(time);
      platform->linAlg
          ->axpbyz(mesh->Nlocal, nrs->g0 / nrs->dt[0], nrs->o_rho, 1.0, o_implicitLT, o_lambda1);
    } else {
      platform->linAlg->axpby(mesh->Nlocal, nrs->g0 / nrs->dt[0], nrs->o_rho, 0.0, o_lambda1);
    }
    return o_lambda1;
  }();

  const auto o_U = [&]()
  {
    auto o_U = platform->o_memPool.reserve<dfloat>(nrs->NVfields * nrs->fieldOffset);
    o_U.copyFrom(platform->options.compareArgs("VELOCITY INITIAL GUESS", "EXTRAPOLATION") && stage == 1
                 ? nrs->o_Ue
                 : nrs->o_U);
 
    if (nrs->uvwSolver) {
      nrs->uvwSolver->solve(o_lambda0, o_lambda1, o_rhs, o_U);
    } else {
      const auto o_rhsX = o_rhs.slice(0 * nrs->fieldOffset);
      const auto o_rhsY = o_rhs.slice(1 * nrs->fieldOffset);
      const auto o_rhsZ = o_rhs.slice(2 * nrs->fieldOffset);
      nrs->uSolver->solve(o_lambda0, o_lambda1, o_rhsX, o_U.slice(0 * nrs->fieldOffset));
      nrs->vSolver->solve(o_lambda0, o_lambda1, o_rhsY, o_U.slice(1 * nrs->fieldOffset));
      nrs->wSolver->solve(o_lambda0, o_lambda1, o_rhsZ, o_U.slice(2 * nrs->fieldOffset));
    }

    return o_U;
  }();

  if (platform->verbose) {
    const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                 mesh->dim,
                                                                 nrs->fieldOffset,
                                                                 mesh->ogs->o_invDegree,
                                                                 o_U,
                                                                 platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("U norm: %.15e\n", debugNorm);
    }
  }

  return o_U;
}

} // namespace tombo