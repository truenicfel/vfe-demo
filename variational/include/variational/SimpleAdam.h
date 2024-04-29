#pragma once

/**
 * @brief (Simple) Adam implementation provided by https://github.com/Paul-Hi
 */
class SimpleAdam
{
public:
    explicit SimpleAdam(double alpha = 0.001, double beta1 = 0.9, double beta2 = 0.999);

    void reset(double alpha = 0.001, double beta1 = 0.9, double beta2 = 0.999);

    double direction(double grad);

private:
    double mFirstMoment, mSecondMoment, mAlpha, mBeta1, mBeta2;
    int mSteps;
};
