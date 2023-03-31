#include "cvode.hpp"
#include "elliptic.h"
#include "inipp.hpp"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "Urst.hpp"
#include <limits>
#include <array>
#include <numeric>
#include "udf.hpp"

#include "timeStepper.hpp"
#include "plugins/lowMach.hpp"
#include "nekrs.hpp"
#include "bdry.hpp"

#ifdef ENABLE_CVODE
// cvode includes
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_mpiplusx.h>
#ifdef ENABLE_CUDA
#include <nvector/nvector_cuda.h>
#endif
#ifdef ENABLE_HIP
#include <nvector/nvector_hip.h>
#endif
#endif

//#define USE_E_VECTOR_LAYOUT 1

namespace {

#ifdef ENABLE_CVODE
sunrealtype *__N_VGetDeviceArrayPointer(N_Vector u)
{
  bool useDevice = false;
  useDevice |= platform->device.mode() == "CUDA";
  useDevice |= platform->device.mode() == "HIP";
  useDevice |= platform->device.mode() == "OPENCL";

  if (useDevice) {
    return N_VGetDeviceArrayPointer(u);
  }
  else {
    return N_VGetArrayPointer_Serial(u);
  }
}
#endif

int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *)returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  return 0;
}

} // namespace

cvode_t::cvode_t(nrs_t *nrs)
{
}

void cvode_t::initialize(nrs_t *nrs)
{
}

cvode_t::~cvode_t()
{
}

void cvode_t::setupEToLMapping(nrs_t *nrs)
{
}

void cvode_t::setupDirichletMask(nrs_t *nrs)
{
}

void cvode_t::extrapolateDirichlet(nrs_t *nrs, dfloat time, int tstep)
{
}

void cvode_t::rhs(nrs_t *nrs, dfloat time, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::jtvRHS(nrs_t *nrs, dfloat time, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::defaultRHS(nrs_t *nrs, int tstep, dfloat time, dfloat t0, occa::memory o_y, occa::memory o_ydot)
{
}

void cvode_t::makeq(nrs_t *nrs, dfloat time)
{
}

void cvode_t::nrsToCv(nrs_t *nrs, occa::memory o_EField, occa::memory o_LField)
{
}

void cvode_t::cvToNrs(nrs_t *nrs, occa::memory o_LField, occa::memory o_EField)
{
}

void cvode_t::solve(nrs_t *nrs, double t0, double t1, int tstep)
{
}

#ifdef ENABLE_CVODE
long cvode_t::numSteps() const
{
  long int nsteps;
  auto retval = CVodeGetNumSteps(cvodeMem, &nsteps);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  return nsteps - prevNsteps;
}

long cvode_t::numRHSEvals() const
{
  long int nrhs;
  auto retval = CVodeGetNumRhsEvals(cvodeMem, &nrhs);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  return nrhs - prevNrhs;
}

long cvode_t::numNonlinSolveIters() const
{
  long int nni;
  auto retval = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  return nni - prevNni;
}

long cvode_t::numLinIters() const
{
  long int nli;
  auto retval = CVodeGetNumLinIters(cvodeMem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1);
  return nli - prevNli;
}
#else
long cvode_t::numSteps() const { return 0; }

long cvode_t::numRHSEvals() const { return 0; }

long cvode_t::numNonlinSolveIters() const { return 0; }

long cvode_t::numLinIters() const { return 0; }
#endif

void cvode_t::printInfo(bool printVerboseInfo) const
{
}