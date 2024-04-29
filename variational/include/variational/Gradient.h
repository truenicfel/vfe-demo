#pragma once

#include <Eigen/Eigen>

#include "variational/SimpleAdam.h"
#include "variational/DoubleBufferedLine.h"

// Helper class that stores optimization parameters such as adam usage, front size, direction and step size.
template<int Dimensions>
class Gradient {

public:

    using Vertex = Eigen::Vector<double, Dimensions>;
    using Line = DoubleBufferedLine<Dimensions>;
    using SampleGradientFunction = std::function<Vertex(const Line&, const int& index)>;

    /**
     * Initialize.
     *
     * @param sampleGradientFunction
     * @param frontSize
     * @param forward
     * @param stepSize (interpreted as learning rate if use adam is true)
     */
    Gradient(const SampleGradientFunction& sampleGradientFunction, const size_t& frontSize, const bool& forward, const bool& useAdam, const double& stepSize)
    : mSampleGradient(sampleGradientFunction)
    , mFrontSize(frontSize)
    , mForward(forward)
    , mUseAdam(useAdam)
    , mAdams(mFrontSize)
    , mStepSize(stepSize)
    {
        if (mUseAdam) {
            for (auto& adamGroup: mAdams) {
                adamGroup.fill(SimpleAdam(stepSize));
            }
        }

    }

    // optimize one step on the given line with the given parameters in ctor
    // useful when conditions need to be checked between optimization steps
    void OptimizeOnce(Line& line) {
        typename Line::IndexRange indexRange = line.GetLastN(mFrontSize);
        if (!mForward) {
            indexRange = line.GetFirstN(mFrontSize);
        }

        // loop over indices that need optimization
#ifdef NDEBUG
#pragma omp parallel for
#endif
        for (int index = indexRange.first; index < indexRange.second; index++) {

            // sample gradient
            Eigen::Vector<double, Dimensions> gradient = mSampleGradient(line, index);
            size_t adamsIndex = index - indexRange.first;
            // pass through adam if necessary
            if (mUseAdam) {
                for (int dimensionIndex = 0; dimensionIndex < Dimensions; ++dimensionIndex) {
                    gradient(dimensionIndex) = mAdams.at(adamsIndex).at(dimensionIndex).direction(gradient(dimensionIndex));
                }
            } else {
                if (gradient.stableNorm() > 1) {
                    gradient.stableNormalize();
                }
                gradient *= mStepSize;
            }

            // "euler step" in gradient direction
            line.Write(index, line.Read(index) - gradient);
        }

        // the line we read from next is the one we optimized last
        //line.Swap();
    }

    // optimize n steps on the given line with the given parameters in ctor
    // useful when no conditions need to be checked between optimization steps
    void OptimizeN(Line& line, int n) {
        typename Line::IndexRange indexRange = line.GetLastN(mFrontSize);
        if (!mForward) {
            indexRange = line.GetFirstN(mFrontSize);
        }

        for (int step = 0; step < n; step++) {

            // loop over indices that need optimization and optimize them once
//#ifdef NDEBUG
//#pragma omp parallel for
//#endif
            for (int index = indexRange.first; index < indexRange.second; index++) {

                // sample gradient
                Eigen::Vector<double, Dimensions> gradient = mSampleGradient(line, index);
                // pass through adam if necessary
                if (mUseAdam) {
                    for (int dimensionIndex = 0; dimensionIndex < Dimensions; ++dimensionIndex) {
                        gradient(dimensionIndex) = mAdams.at(index).at(dimensionIndex).direction(
                                gradient(dimensionIndex));
                    }
                } else {
                    gradient *= mStepSize;
                }

                // "euler step" in gradient direction
                line.Write(index, line.Read(index) - gradient);
            }

            // the line we read from next is the one we optimized last
            line.Swap();

        }
    }


private:

    SampleGradientFunction mSampleGradient;

    size_t mFrontSize;

    bool mForward;

    bool mUseAdam;

    std::vector<std::array<SimpleAdam, Dimensions>> mAdams;

    double mStepSize;
};
