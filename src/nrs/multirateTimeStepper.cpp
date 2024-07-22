#include "nrs.hpp"
#include "udf.hpp"

void nrs_t::initOuterStep(double time, dfloat _dt, int tstep)
{
  saveSolutionState();
  tStepOuterStart = tstep;
  timeOuterStart = time;

  if (tstep == 1) {
    const bool exchangeAllTimes = false;
    const bool lagState = true;
    neknek->exchange(exchangeAllTimes, lagState);
  }

  neknek->exchangeTimes(time);
}

void nrs_t::finishOuterStep() {}

bool nrs_t::runOuterStep(std::function<bool(int)> convergenceCheck, int stage)
{
  int innerSteps = 1;
  platform->options.getArgs("MULTIRATE STEPS", innerSteps);

  auto tstep = tStepOuterStart;
  auto time = timeOuterStart;

  int requiredCorrectorSteps = 0;
  platform->options.getArgs("MULTIRATE CORRECTOR STEPS", requiredCorrectorSteps);

  const auto correctorStep = stage - 1;
  const bool predictorStep = correctorStep == 0;
  neknek->setPredictor(predictorStep);

  const bool outerConverged = correctorStep >= requiredCorrectorSteps;

  if (stage > 1) {
    restoreSolutionState();
  }

  for (int step = 0; step < innerSteps; ++step) {
    bool last = step == innerSteps - 1;
    initInnerStep(time, dt[0], tstep);
    time += setPrecision(dt[0], 5);

    int innerStage = 1;
    bool converged = false;
    do {
      converged = runInnerStep(convergenceCheck, innerStage++);
    } while (!converged);

    finishInnerStep();

    if (!(last && outerConverged)) {
      printInfo(time, tstep, true, true);
    }

    tstep++;
  }

  const bool exchangeAllTimeStates = outerConverged;
  const bool lagState = outerConverged;
  neknek->exchange(exchangeAllTimeStates, lagState);
  if (!outerConverged) {
    neknek->setCorrectorTime(time);
  }

  return outerConverged;
}