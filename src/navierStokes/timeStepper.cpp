#include "nrs.hpp"
#include "neknek.hpp"
#include "avm.hpp"
#include "cfl.hpp"
#include "constantFlowRate.hpp"
#include "nekInterfaceAdapter.hpp"
#include "timeStepper.hpp"
#include "tombo.hpp"
#include "subCycling.hpp"
#include "udf.hpp"
#include "bcMap.hpp"
#include "bdry.hpp"
#include "Urst.hpp"

static void lagFields(nrs_t *nrs)
{
  // lag velocity
  for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
    const auto Nbyte = (nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset;
    nrs->o_U.copyFrom(nrs->o_U, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
  }

  // lag scalars
  if (nrs->Nscalar) {
    auto cds = nrs->cds;
    if(cds->anyEllipticSolver){
      for (int s = std::max(cds->nBDF, cds->nEXT); s > 1; s--) {
        const auto Nbyte = cds->fieldOffsetSum * sizeof(dfloat);
        cds->o_S.copyFrom(cds->o_S, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
      }
    }
  }

  // lag mesh velocity
  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  if (movingMesh) {
    auto mesh = nrs->_mesh;
    for (int s = std::max(nrs->nEXT, mesh->nAB); s > 1; s--) {
      const auto Nbyte = (nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset;
      mesh->o_U.copyFrom(mesh->o_U, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
    }
  }
}

static void extrapolate(nrs_t *nrs)
{
  mesh_t *mesh = nrs->meshV;

  {
    nrs->extrapolateKernel(mesh->Nlocal,
                           nrs->NVfields,
                           nrs->nEXT,
                           nrs->fieldOffset,
                           nrs->o_coeffEXT,
                           nrs->o_U,
                           nrs->o_Ue);

    if (nrs->flow && platform->options.compareArgs("LOWMACH", "TRUE")) {
      if (nrs->pSolver->allNeumann) {
        nrs->p0the = 0.0;
        for (int ext = 0; ext < nrs->nEXT; ++ext) {
          nrs->p0the += nrs->coeffEXT[ext] * nrs->p0th[ext];
        }
      }
    }
  }

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    if (nrs->cht)
      mesh = nrs->_mesh;
    nrs->extrapolateKernel(mesh->Nlocal,
                           nrs->NVfields,
                           nrs->nEXT,
                           nrs->fieldOffset,
                           nrs->o_coeffEXT,
                           mesh->o_U,
                           mesh->o_Ue);
  }

  if(nrs->Nscalar > 0) {
    auto cds = nrs->cds;
    if(!cds->o_Se.isInitialized()) return;
    nrs->extrapolateKernel(cds->mesh[0]->Nlocal,
                           cds->NSfields,
                           cds->nEXT,
                           cds->fieldOffset[0],
                           cds->o_coeffEXT,
                           cds->o_S,
                           cds->o_Se);
  }
}

static void computeDivUErr(nrs_t *nrs, dfloat &divUErrVolAvg, dfloat &divUErrL2)
{
  mesh_t *mesh = nrs->meshV;
  nrs->divergenceVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              nrs->fieldOffset,
                              nrs->o_U,
                              platform->o_mempool.slice0);

  double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
  flops *= static_cast<double>(mesh->Nelements);

  platform->flopCounter->add("divergenceVolumeKernel", flops);

  oogs::startFinish(platform->o_mempool.slice0, 1, nrs->fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, platform->o_mempool.slice0);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, nrs->o_div, -1.0, platform->o_mempool.slice0);
  divUErrL2 = platform->linAlg->weightedNorm2(mesh->Nlocal,
                                              mesh->o_LMM,
                                              platform->o_mempool.slice0,
                                              platform->comm.mpiComm) /
              sqrt(mesh->volume);

  divUErrVolAvg = platform->linAlg->innerProd(mesh->Nlocal,
                                              mesh->o_LMM,
                                              platform->o_mempool.slice0,
                                              platform->comm.mpiComm) /
                  mesh->volume;
  divUErrVolAvg = std::abs(divUErrVolAvg);
}

namespace timeStepper {

void advectionFlops(mesh_t *mesh, int Nfields)
{
  const auto cubNq = mesh->cubNq;
  const auto cubNp = mesh->cubNp;
  const auto Nq = mesh->Nq;
  const auto Np = mesh->Np;
  const auto Nelements = mesh->Nelements;
  double flopCount = 0.0; // per elem basis
  if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
    flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq); // interpolation
    flopCount += 6. * cubNp * cubNq;                                       // apply Dcub
    flopCount += 5 * cubNp; // compute advection term on cubature mesh
    flopCount += mesh->Np;  // weight by inv. mass matrix
  }
  else {
    flopCount += 8 * (Np * Nq + Np);
  }

  flopCount *= Nelements;
  flopCount *= Nfields;

  platform->flopCounter->add("advection", flopCount);
}

void adjustDt(nrs_t* nrs, int tstep)
{
  const double TOLToZero = 1e-12;
  bool initialTimeStepProvided = true;
  if (nrs->dt[0] < TOLToZero && tstep == 1) {
    nrs->dt[0] = 1.0; // startup without any initial timestep guess
    initialTimeStepProvided = false;
  }

  double targetCFL;
  platform->options.getArgs("TARGET CFL", targetCFL);

  const double CFLmax = 1.2 * targetCFL;
  const double CFLmin = 0.8 * targetCFL;

  const double CFL = computeCFL(nrs);

  if (!initialTimeStepProvided) {
    if (CFL > TOLToZero) {
      nrs->dt[0] = targetCFL / CFL * nrs->dt[0];
      nrs->unitTimeCFL = CFL / nrs->dt[0];
    }
    else {
      // estimate from userf
      if (udf.uEqnSource) {
        platform->linAlg->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
        platform->timer.tic("udfUEqnSource", 1);
        double startTime;
        platform->options.getArgs("START TIME", startTime);
        udf.uEqnSource(nrs, startTime, nrs->o_U, nrs->o_FU);
        platform->timer.toc("udfUEqnSource");

        occa::memory o_FUx = nrs->o_FU + (0 * sizeof(dfloat)) * nrs->fieldOffset;
        occa::memory o_FUy = nrs->o_FU + (1 * sizeof(dfloat)) * nrs->fieldOffset;
        occa::memory o_FUz = nrs->o_FU + (2 * sizeof(dfloat)) * nrs->fieldOffset;

        platform->linAlg->abs(3 * nrs->fieldOffset, nrs->o_FU);

        const double maxFUx = platform->linAlg->max(nrs->meshV->Nlocal, o_FUx, platform->comm.mpiComm);
        const double maxFUy = platform->linAlg->max(nrs->meshV->Nlocal, o_FUy, platform->comm.mpiComm);
        const double maxFUz = platform->linAlg->max(nrs->meshV->Nlocal, o_FUz, platform->comm.mpiComm);
        const double maxFU = std::max({maxFUx, maxFUy, maxFUz});

        const double minRho = platform->linAlg->min(nrs->meshV->Nlocal, nrs->o_rho, platform->comm.mpiComm);
        const double maxU = maxFU / minRho;
        const double *x = nrs->meshV->x;
        const double *y = nrs->meshV->y;
        const double *z = nrs->meshV->z;
        double lengthScale = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]) +
                                  (z[0] - z[1]) * (z[0] - z[1]));

        MPI_Allreduce(MPI_IN_PLACE, &lengthScale, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);
        if (maxU > TOLToZero) {
          nrs->dt[0] = sqrt(targetCFL * lengthScale / maxU);
        }
        else {
          nrsAbort(platform->comm.mpiComm,
                   EXIT_FAILURE,
                   "%s\n",
                   "Zero velocity and body force! Please specify an initial timestep!");
        }
      }
    }
    nrs->CFL = CFL;
    return;
  }

  const double unitTimeCFLold = (tstep == 1) ? CFL / nrs->dt[0] : nrs->unitTimeCFL;

  const double CFLold = (tstep == 1) ? CFL : nrs->CFL;

  nrs->CFL = CFL;
  nrs->unitTimeCFL = CFL / nrs->dt[0];

  const double CFLpred = 2.0 * nrs->CFL - CFLold;

  const double TOL = 0.001;

  if (nrs->CFL > CFLmax || CFLpred > CFLmax || nrs->CFL < CFLmin) {
    const double A = (nrs->unitTimeCFL - unitTimeCFLold) / nrs->dt[0];
    const double B = nrs->unitTimeCFL;
    const double C = -targetCFL;
    const double descriminant = B * B - 4 * A * C;
    nrs->dt[1] = nrs->dt[0];
    if (descriminant <= 0.0) {
      nrs->dt[0] = nrs->dt[0] * (targetCFL / nrs->CFL);
    }
    else if (std::abs((nrs->unitTimeCFL - unitTimeCFLold) / nrs->unitTimeCFL) < TOL) {
      nrs->dt[0] = nrs->dt[0] * (targetCFL / nrs->CFL);
    }
    else {
      const double dtLow = (-B + sqrt(descriminant)) / (2.0 * A);
      const double dtHigh = (-B - sqrt(descriminant)) / (2.0 * A);
      if (dtHigh > 0.0 && dtLow > 0.0) {
        nrs->dt[0] = std::min(dtLow, dtHigh);
      }
      else if (dtHigh <= 0.0 && dtLow <= 0.0) {
        nrs->dt[0] = nrs->dt[0] * targetCFL / nrs->CFL;
      }
      else {
        nrs->dt[0] = std::max(dtHigh, dtLow);
      }
    }

    // limit dt change
    dfloat maxAdjustDtRatio = 1;
    dfloat minAdjustDtRatio = 1;
    platform->options.getArgs("MAX ADJUST DT RATIO", maxAdjustDtRatio);
    platform->options.getArgs("MIN ADJUST DT RATIO", minAdjustDtRatio);
    if (tstep > 1)
      nrs->dt[0] = std::max(nrs->dt[0], minAdjustDtRatio * nrs->dt[1]);
    if (tstep > 1)
      nrs->dt[0] = std::min(nrs->dt[0], maxAdjustDtRatio * nrs->dt[1]);
  }
}

void initStep(nrs_t *nrs, dfloat time, dfloat dt, int tstep)
{
  nrs->timePrevious = time;
  nrs->tstep = tstep;

  cds_t *cds = nrs->cds;

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  setDt(nrs, dt, tstep);

  extrapolate(nrs);

  if (nrs->Nsubsteps) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    const auto NbyteCubature = (nrs->NVfields * sizeof(dfloat)) * nrs->cubatureOffset;
    for (int s = nrs->nEXT; s > 1; s--) {
      const auto Nbyte = nrs->fieldOffset * sizeof(dfloat);
      if (movingMesh) {
        mesh->o_divU.copyFrom(mesh->o_divU, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
        nrs->o_relUrst.copyFrom(nrs->o_relUrst,
                                NbyteCubature,
                                (s - 1) * NbyteCubature,
                                (s - 2) * NbyteCubature);
      }
      else {
        nrs->o_Urst.copyFrom(nrs->o_Urst, NbyteCubature, (s - 1) * NbyteCubature, (s - 2) * NbyteCubature);
      }
    }
    if (movingMesh) {
      double flops = 18 * (mesh->Np * mesh->Nq + mesh->Np);
      flops *= static_cast<double>(mesh->Nelements);
      nrs->divergenceVolumeKernel(mesh->Nelements,
                                  mesh->o_vgeo,
                                  mesh->o_D,
                                  nrs->fieldOffset,
                                  mesh->o_U,
                                  mesh->o_divU);
      platform->flopCounter->add("divergenceVolumeKernel", flops);
    }
  }

  computeUrst(nrs);

  if (nrs->Nscalar) {
    if(cds->anyEllipticSolver) {
      platform->timer.tic("makeq", 1);
      platform->linAlg->fillKernel(cds->fieldOffsetSum, 0.0, cds->o_FS);
      makeq(nrs, time, tstep, cds->o_FS, cds->o_BF);
      platform->timer.toc("makeq");
    }
  }

  if (nrs->flow) {
    platform->timer.tic("makef", 1);
    platform->linAlg->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
    makef(nrs, time, tstep, nrs->o_FU, nrs->o_BF);
    platform->timer.toc("makef");
  }

  if (movingMesh) {
    mesh_t *mesh = nrs->_mesh;
    for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
      const auto NbyteScalar = nrs->fieldOffset * sizeof(dfloat);
      mesh->o_LMM.copyFrom(mesh->o_LMM, NbyteScalar, (s - 1) * NbyteScalar, (s - 2) * NbyteScalar);
      mesh->o_invLMM.copyFrom(mesh->o_invLMM, NbyteScalar, (s - 1) * NbyteScalar, (s - 2) * NbyteScalar);
    }

    mesh->move();

    if (nrs->flow) {
      if (bcMap::unalignedMixedBoundary("velocity")) {
        createZeroNormalMask(nrs,
                             nrs->meshV,
                             nrs->uvwSolver->o_EToB,
                             nrs->o_EToBVVelocity,
                             nrs->o_zeroNormalMaskVelocity);
      }
    }

    if (!platform->options.compareArgs("MESH SOLVER", "NONE")) {
      if (bcMap::unalignedMixedBoundary("mesh")) {
        createZeroNormalMask(nrs,
                             mesh,
                             nrs->meshSolver->o_EToB,
                             nrs->o_EToBVMeshVelocity,
                             nrs->o_zeroNormalMaskMeshVelocity);
      }
    }

    if (nrs->cht)
      nrs->meshV->computeInvLMM();
  }

}

void finishStep(nrs_t *nrs)
{
  nrs->dt[2] = nrs->dt[1];
  nrs->dt[1] = nrs->dt[0];

  if (nrs->flow) {
    if (nrs->pSolver->allNeumann && platform->options.compareArgs("LOWMACH", "TRUE"))
      nrs->p0the = nrs->p0th[0];
  }

  platform->device.finish();
}

bool runStep(nrs_t *nrs, std::function<bool(int)> convergenceCheck, int stage)
{
  mesh_t *mesh = nrs->meshV;
  cds_t *cds = nrs->cds;

  const auto tstep = nrs->tstep;
  const auto timeNew = nrs->timePrevious + nrs->dt[0];

  const int isOutputStep = nrs->isOutputStep;

  if (nrs->neknek)
    nrs->neknek->updateBoundary(nrs, tstep, stage);
  
  if (nrs->cvode)
    scalarSolveCvode(nrs, nrs->timePrevious, timeNew, cds->o_S, stage, tstep);
  
  if (stage == 1)
    lagFields(nrs);

  applyDirichlet(nrs, timeNew);

  if (nrs->Nscalar)
    scalarSolve(nrs, timeNew, cds->o_S, stage);
  
  if(udf.postScalar)
    udf.postScalar(nrs, timeNew, tstep);

  evaluateProperties(nrs, timeNew);

  if (udf.div) {
    platform->timer.tic("udfDiv", 1);
    platform->linAlg->fill(mesh->Nlocal, 0.0, nrs->o_div);
    udf.div(nrs, timeNew, nrs->o_div);
    platform->timer.toc("udfDiv");
  }
  
  if(udf.preFluid)
    udf.preFluid(nrs, timeNew, tstep);

  if (nrs->flow)
    fluidSolve(nrs, timeNew, nrs->o_P, nrs->o_U, stage, tstep);

  if (!platform->options.compareArgs("MESH SOLVER", "NONE"))
    meshSolve(nrs, timeNew, nrs->meshV->o_U, stage);

  nrs->timeStepConverged = convergenceCheck(stage);

  platform->timer.tic("udfExecuteStep", 1);
  nek::ifoutfld(0);
  nrs->isOutputStep = 0;
  if (isOutputStep && nrs->timeStepConverged) {
    nek::ifoutfld(1);
    nrs->isOutputStep = 1;
  }
  if (udf.executeStep)
    udf.executeStep(nrs, timeNew, tstep);
  platform->timer.toc("udfExecuteStep");

  if (!nrs->timeStepConverged)
    printInfo(nrs, timeNew, tstep, false, true);

  return nrs->timeStepConverged;
}

void setDt(nrs_t *nrs, double dt, int tstep)
{
  nrs->dt[0] = dt;
  nrsCheck(std::isnan(nrs->dt[0]) || std::isinf(nrs->dt[0]),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s",
           "Unreasonable dt! Dying ...\n");

  nrs->idt = 1 / nrs->dt[0];

  const int bdfOrder = std::min(tstep, nrs->nBDF);
  const int extOrder = std::min(tstep, nrs->nEXT);

  nek::bdfCoeff(&nrs->g0, nrs->coeffBDF, nrs->dt, bdfOrder);
  nek::extCoeff(nrs->coeffEXT, nrs->dt, extOrder, bdfOrder);

  for (int i = nrs->nEXT; i > extOrder; i--)
    nrs->coeffEXT[i - 1] = 0.0;
  for (int i = nrs->nBDF; i > bdfOrder; i--)
    nrs->coeffBDF[i - 1] = 0.0;

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t *mesh = nrs->meshV;
    if (nrs->cht)
      mesh = nrs->cds->mesh[0];
    const int meshOrder = std::min(tstep, mesh->nAB);
    nek::coeffAB(mesh->coeffAB, nrs->dt, meshOrder);
    for (int i = 0; i < meshOrder; ++i)
      mesh->coeffAB[i] *= nrs->dt[0];
    for (int i = mesh->nAB; i > meshOrder; i--)
      mesh->coeffAB[i - 1] = 0.0;
    mesh->o_coeffAB.copyFrom(mesh->coeffAB, mesh->nAB * sizeof(dfloat));
  }

  nrs->ig0 = 1.0 / nrs->g0;
  nrs->o_coeffBDF.copyFrom(nrs->coeffBDF);
  nrs->o_coeffEXT.copyFrom(nrs->coeffEXT);

  if (nrs->Nscalar) {
    nrs->cds->idt = 1 / nrs->cds->dt[0];
    nrs->cds->g0 = nrs->g0;
    nrs->cds->ig0 = nrs->ig0;
  }
}

void makeq(nrs_t *nrs, dfloat time, int tstep, occa::memory o_FS, occa::memory o_BF)
{
  cds_t *cds = nrs->cds;

  if (udf.sEqnSource) {
    platform->timer.tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_FS);
    platform->timer.toc("udfSEqnSource");
  }

  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is] || cds->cvodeSolve[is])
      continue;

    std::string sid = scalarDigitStr(is);

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];
    const dlong isOffset = cds->fieldOffsetScan[is];

    if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "HPFRT")) {
      cds->filterRTKernel(cds->meshV->Nelements,
                          is,
                          1,
                          nrs->fieldOffset,
                          cds->o_applyFilterRT,
                          cds->o_filterMT,
                          cds->o_filterS,
                          cds->o_rho,
                          cds->o_S,
                          o_FS);

      double flops = 6 * mesh->Np * mesh->Nq + 4 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalarFilterRT", flops);
    }
    const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    if (movingMesh && !cds->Nsubsteps) {
      cds->advectMeshVelocityKernel(cds->meshV->Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_D,
                                    isOffset,
                                    cds->vFieldOffset,
                                    cds->o_rho,
                                    mesh->o_U,
                                    cds->o_S,
                                    o_FS);
      double flops = 18 * mesh->Np * mesh->Nq + 21 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalar advectMeshVelocity", flops);
    }

    occa::memory o_Usubcycling = platform->o_mempool.slice0;
    if (platform->options.compareArgs("ADVECTION", "TRUE")) {
      if (cds->Nsubsteps) {
        if (movingMesh)
          o_Usubcycling =
              scalarSubCycleMovingMesh(cds, std::min(tstep, cds->nEXT), time, is, cds->o_U, cds->o_S);
        else
          o_Usubcycling = scalarSubCycle(cds, std::min(tstep, cds->nEXT), time, is, cds->o_U, cds->o_S);
      } else {
        if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
          cds->strongAdvectionCubatureVolumeKernel(cds->meshV->Nelements,
                                                   1,
                                                   0,
                                                   mesh->o_LMM,
                                                   mesh->o_vgeo,
                                                   mesh->o_cubDiffInterpT,
                                                   mesh->o_cubInterpT,
                                                   mesh->o_cubProjectT,
                                                   cds->o_compute + is * sizeof(dlong),
                                                   cds->o_fieldOffsetScan + is * sizeof(dlong),
                                                   cds->vFieldOffset,
                                                   cds->vCubatureOffset,
                                                   cds->o_S,
                                                   cds->o_Urst,
                                                   cds->o_rho,
                                                   o_FS);
        }
        else {
          cds->strongAdvectionVolumeKernel(cds->meshV->Nelements,
                                           1,
                                           0,
                                           mesh->o_LMM,
                                           mesh->o_vgeo,
                                           mesh->o_D,
                                           cds->o_compute + is * sizeof(dlong),
                                           cds->o_fieldOffsetScan + is * sizeof(dlong),
                                           cds->vFieldOffset,
                                           cds->o_S,
                                           cds->o_Urst,
                                           cds->o_rho,
                                           o_FS);
        }
        advectionFlops(cds->mesh[0], 1);
      }
    }
    else {
      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, o_Usubcycling);
    }

    cds->sumMakefKernel(mesh->Nlocal,
                        mesh->o_LMM,
                        cds->idt,
                        cds->o_coeffEXT,
                        cds->o_coeffBDF,
                        cds->fieldOffsetSum,
                        cds->fieldOffset[is],
                        isOffset,
                        cds->o_S,
                        o_Usubcycling,
                        o_FS,
                        cds->o_rho,
                        o_BF);

    dfloat scalarSumMakef = (3 * cds->nEXT + 3) * static_cast<double>(mesh->Nlocal);
    scalarSumMakef += (cds->Nsubsteps) ? mesh->Nlocal : 3 * cds->nBDF * static_cast<double>(mesh->Nlocal);
    platform->flopCounter->add("scalarSumMakef", scalarSumMakef);
  }

  for (int s = std::max(cds->nBDF, cds->nEXT); s > 1; s--) {
    const auto Nbyte = cds->fieldOffsetSum * sizeof(dfloat);
    o_FS.copyFrom(o_FS, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
  }
}

void scalarSolve(nrs_t *nrs, dfloat time, occa::memory o_S, int stage)
{
  cds_t *cds = nrs->cds;

  // avoid invoking the scalarSolve tic in the case there are no elliptic scalars
  if (!cds->anyEllipticSolver)
    return;

  platform->timer.tic("scalarSolve", 1);
  for (int is = 0; is < cds->NSfields; is++) {
    if (!cds->compute[is] || cds->cvodeSolve[is])
      continue;

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];

    cds->setEllipticCoeffKernel(mesh->Nlocal,
                                cds->g0 * cds->idt,
                                cds->fieldOffsetScan[is],
                                nrs->fieldOffset,
                                (cds->o_BFDiag.size()) ? 1 : 0,
                                cds->o_diff,
                                cds->o_rho,
                                cds->o_BFDiag,
                                cds->o_ellipticCoeff);

    occa::memory o_Snew = cdsSolve(is, cds, time, stage);
    o_Snew.copyTo(o_S, cds->fieldOffset[is] * sizeof(dfloat), cds->fieldOffsetScan[is] * sizeof(dfloat));
  }
  platform->timer.toc("scalarSolve");
}

void scalarSolveCvode(nrs_t *nrs, dfloat tn, dfloat time, occa::memory o_S, int stage, int tstep)
{
  // CVODE is not supported with sub-stepping
  if (!nrs->cvode || stage > 1)
    return;

  nrs->cvode->solve(tn, time, tstep);
}

void makef(
    nrs_t *nrs, dfloat time, int tstep, occa::memory o_FU, occa::memory o_BF) {
  mesh_t *mesh = nrs->meshV;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  if (udf.uEqnSource) {
    platform->timer.tic("udfUEqnSource", 1);
    udf.uEqnSource(nrs, time, nrs->o_U, o_FU);
    platform->timer.toc("udfUEqnSource");
  }

  if (platform->options.compareArgs("VELOCITY REGULARIZATION METHOD", "HPFRT")) {
    nrs->filterRTKernel(mesh->Nelements, nrs->o_filterMT, nrs->filterS, nrs->fieldOffset, nrs->o_U, o_FU);
    double flops = 24 * mesh->Np * mesh->Nq + 3 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocityFilterRT", flops);
  }

  if (movingMesh && !nrs->Nsubsteps) {
    nrs->advectMeshVelocityKernel(mesh->Nelements,
                                  mesh->o_vgeo,
                                  mesh->o_D,
                                  nrs->fieldOffset,
                                  mesh->o_U,
                                  nrs->o_U,
                                  o_FU);
    double flops = 54 * mesh->Np * mesh->Nq + 63 * mesh->Np;
    flops *= static_cast<double>(mesh->Nelements);
    platform->flopCounter->add("velocity advectMeshVelocity", flops);
  }

  occa::memory o_Usubcycling = platform->o_mempool.slice0;
  if (platform->options.compareArgs("ADVECTION", "TRUE")) {
    if (nrs->Nsubsteps) {
      if (movingMesh)
        o_Usubcycling = velocitySubCycleMovingMesh(nrs, std::min(tstep, nrs->nEXT), time, nrs->o_U);
      else
        o_Usubcycling = velocitySubCycle(nrs, std::min(tstep, nrs->nEXT), time, nrs->o_U);
    }
    else {
      if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
        nrs->strongAdvectionCubatureVolumeKernel(mesh->Nelements,
                                                 mesh->o_vgeo,
                                                 mesh->o_cubDiffInterpT,
                                                 mesh->o_cubInterpT,
                                                 mesh->o_cubProjectT,
                                                 nrs->fieldOffset,
                                                 nrs->cubatureOffset,
                                                 nrs->o_U,
                                                 nrs->o_Urst,
                                                 platform->o_mempool.slice0);
      else
        nrs->strongAdvectionVolumeKernel(mesh->Nelements,
                                         mesh->o_vgeo,
                                         mesh->o_D,
                                         nrs->fieldOffset,
                                         nrs->o_U,
                                         nrs->o_Urst,
                                         platform->o_mempool.slice0);

      platform->linAlg->axpby(nrs->NVfields * nrs->fieldOffset, -1.0, platform->o_mempool.slice0, 1.0, o_FU);

      advectionFlops(nrs->meshV, nrs->NVfields);
    }
  }
  else {
    if (nrs->Nsubsteps)
      platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Usubcycling);
  }

  nrs->sumMakefKernel(mesh->Nlocal,
                      mesh->o_LMM,
                      nrs->idt,
                      nrs->o_coeffEXT,
                      nrs->o_coeffBDF,
                      nrs->fieldOffset,
                      nrs->o_U,
                      o_Usubcycling,
                      o_FU,
                      o_BF);

  dfloat sumMakefFlops = 0.0;
  if (nrs->Nsubsteps) {
    sumMakefFlops += (6 + 6 * nrs->nEXT) * static_cast<double>(mesh->Nlocal);
  }
  else {
    sumMakefFlops += (6 * nrs->nEXT + 12 * nrs->nBDF) * static_cast<double>(mesh->Nlocal);
  }

  platform->flopCounter->add("sumMakef", sumMakefFlops);

  if (verbose) {
    const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                 nrs->NVfields,
                                                                 nrs->fieldOffset,
                                                                 mesh->ogs->o_invDegree,
                                                                 o_BF,
                                                                 platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0)
      printf("BF norm: %.15e\n", debugNorm);
  }

  for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
    const auto Nbyte = (nrs->NVfields * sizeof(dfloat)) * nrs->fieldOffset;
    o_FU.copyFrom(o_FU, Nbyte, (s - 1) * Nbyte, (s - 2) * Nbyte);
  }
}

void fluidSolve(nrs_t *nrs, dfloat time, occa::memory o_P, occa::memory o_U, int stage, int tstep)
{
  mesh_t *mesh = nrs->meshV;

  platform->timer.tic("pressureSolve", 1);
  nrs->setEllipticCoeffPressureKernel(mesh->Nlocal, nrs->fieldOffset, nrs->o_rho, nrs->o_ellipticCoeff);
  occa::memory o_Pnew = tombo::pressureSolve(nrs, time, stage);
  o_P.copyFrom(o_Pnew, nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("pressureSolve");

  platform->timer.tic("velocitySolve", 1);
  nrs->setEllipticCoeffKernel(mesh->Nlocal,
                              nrs->g0 * nrs->idt,
                              0 * nrs->fieldOffset,
                              nrs->fieldOffset,
                              0,
                              nrs->o_mue,
                              nrs->o_rho,
                              o_NULL,
                              nrs->o_ellipticCoeff);

  occa::memory o_Unew = tombo::velocitySolve(nrs, time, stage);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

  platform->timer.toc("velocitySolve");

  if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
    ConstantFlowRate::adjust(nrs, tstep, time);
  }
}

void printInfo(nrs_t *nrs, dfloat time, int tstep, bool printStepInfo, bool printVerboseInfo)
{
  cds_t *cds = nrs->cds;

  const double elapsedStep = platform->timer.query("elapsedStep", "DEVICE:MAX");
  const double elapsedStepSum = platform->timer.query("elapsedStepSum", "DEVICE:MAX");
  bool verboseInfo = platform->options.compareArgs("VERBOSE SOLVER INFO", "TRUE");
  const dfloat cfl = computeCFL(nrs);
  dfloat divUErrVolAvg, divUErrL2;

  if (verboseInfo) {
    computeDivUErr(nrs, divUErrVolAvg, divUErrL2);
  }
  if (platform->comm.mpiRank == 0) {
    if (verboseInfo && printVerboseInfo) {
      bool cvodePrinted = false;
      for (int is = 0; is < nrs->Nscalar; is++) {
        if (cds->compute[is] && !cds->cvodeSolve[is]) {
          elliptic_t *solver = cds->solver[is];
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projS%02d  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     is,
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("S%02d      : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 is,
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
        } else if (cds->cvodeSolve[is] && !cvodePrinted) {
          nrs->cvode->printInfo(true);
          cvodePrinted = true;
        }
      }

      if (nrs->neknek) {
        printf("neknek   : sync %.2e  exchange %.2e\n", nrs->neknek->tSync, nrs->neknek->tExch);
      }

      if (nrs->flow) {
        elliptic_t *solver = nrs->pSolver;
        if (solver->solutionProjection) {
          const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
          if (prevVecs > 0) {
            printf("projP    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                   solver->res00Norm,
                   solver->res0Norm,
                   solver->res00Norm / solver->res0Norm,
                   prevVecs,
                   solver->solutionProjection->getMaxNumVecsProjection());
          }
        }
        printf("P        : iter %03d  resNorm0 %.2e  resNorm %.2e\n",
               solver->Niter,
               solver->res0Norm,
               solver->resNorm);

        if (nrs->uvwSolver) {
          solver = nrs->uvwSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projUVW  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("UVW      : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e  divErrNorms %.2e %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm,
                 divUErrVolAvg,
                 divUErrL2);
        }
        else {
          solver = nrs->uSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projU    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("U        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e  divErrNorms %.2e %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm,
                 divUErrVolAvg,
                 divUErrL2);
          solver = nrs->vSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projV    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("V        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
          solver = nrs->wSolver;
          if (solver->solutionProjection) {
            const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
            if (prevVecs > 0) {
              printf("projW    : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                     solver->res00Norm,
                     solver->res0Norm,
                     solver->res00Norm / solver->res0Norm,
                     prevVecs,
                     solver->solutionProjection->getMaxNumVecsProjection());
            }
          }
          printf("W        : iter %03d  resNorm0 %.2e  "
                 "resNorm %.2e\n",
                 solver->Niter,
                 solver->res0Norm,
                 solver->resNorm);
        }
      }

      if (nrs->meshSolver) {
        elliptic_t *solver = nrs->meshSolver;
        if (solver->solutionProjection) {
          const int prevVecs = solver->solutionProjection->getPrevNumVecsProjection();
          if (prevVecs > 0) {
            printf("projMSH  : resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
                   solver->res00Norm,
                   solver->res0Norm,
                   solver->res00Norm / solver->res0Norm,
                   prevVecs,
                   solver->solutionProjection->getMaxNumVecsProjection());
          }
        }
        printf("MSH      : iter %03d  resNorm0 %.2e  resNorm %.2e\n",
               solver->Niter,
               solver->res0Norm,
               solver->resNorm);
      }
    }

    if (platform->options.compareArgs("CONSTANT FLOW RATE", "TRUE")) {
      ConstantFlowRate::printInfo(nrs->meshV, verboseInfo && printVerboseInfo);
    }

    if (printStepInfo)
      printf("step= %d  t= %.8e  dt=%.1e  C= %.2f", tstep, time, nrs->dt[0], cfl);

    if (!verboseInfo) {
      bool cvodePrinted = false;
      for (int is = 0; is < nrs->Nscalar; is++)
        if (cds->compute[is] && !cds->cvodeSolve[is]){
          printf("  S: %d", cds->solver[is]->Niter);
        } else if (cds->cvodeSolve[is] && !cvodePrinted){
          nrs->cvode->printInfo(false);
          cvodePrinted = true;
        }

      if (nrs->flow) {
        printf("  P: %d", nrs->pSolver->Niter);
        if (nrs->uvwSolver)
          printf("  UVW: %d", nrs->uvwSolver->Niter);
        else
          printf("  U: %d  V: %d  W: %d", nrs->uSolver->Niter, nrs->vSolver->Niter, nrs->wSolver->Niter);
      }
      if (nrs->meshSolver)
        printf("  MSH: %d", nrs->meshSolver->Niter);
    }

    if (nrs->timeStepConverged && printStepInfo)
      printf("  elapsedStep= %.2es  elapsedStepSum= %.5es\n", elapsedStep, elapsedStepSum);
  }

  bool largeCFLCheck = (cfl > 30) && numberActiveFields(nrs);

  nrsCheck(largeCFLCheck || std::isnan(cfl) || std::isinf(cfl),
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "Unreasonable CFL!");
}


} // namespace timeStepper
