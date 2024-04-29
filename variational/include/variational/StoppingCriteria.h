#pragma once

#include <Eigen/Eigen>

#include "variational/DoubleBufferedLine.h"

template<int Dimensions>
bool leftDomain(const Eigen::Vector<double, Dimensions>& vertex, const Eigen::AlignedBox<double, Dimensions>& domain) {
    return !domain.contains(vertex);
}

template<int Dimensions>
bool lineClosed(const DoubleBufferedLine<Dimensions>& line, const double& epsilon) {
    return line.IsClosed(epsilon);
}