#include "variational/SimpleAdam.h"

#include <cmath>

SimpleAdam::SimpleAdam(double alpha, double beta1, double beta2)
    : mFirstMoment(0.0)
    , mSecondMoment(0.0)
    , mAlpha(alpha)
    , mBeta1(beta1)
    , mBeta2(beta2)
    , mSteps(0)
{
}

void SimpleAdam::reset(double alpha, double beta1, double beta2)
{
    mFirstMoment  = 0.0;
    mSecondMoment = 0.0;
    mAlpha        = alpha;
    mBeta1        = beta1;
    mBeta2        = beta2;
    mSteps        = 0;
}

double SimpleAdam::direction(double grad)
{
    int t         = mSteps + 1;
    mFirstMoment  = mBeta1 * mFirstMoment + (1.0 - mBeta1) * grad;
    mSecondMoment = mBeta2 * mSecondMoment + (1.0 - mBeta2) * grad * grad;

    // bias correction

    double beta1T = std::pow(mBeta1, t);
    double beta2T = std::pow(mBeta2, t);

    double alphaT = mAlpha * std::sqrt(1.0 - beta2T) / (1.0 - beta1T);

    mSteps++;

    return alphaT * mFirstMoment / (std::sqrt(mSecondMoment) + 1e-8);
}

