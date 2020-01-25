#pragma once

#include <Eigen/Dense>
#include <vector>

namespace lambdatwist {

    struct CameraPose {
        Eigen::Matrix3d R;
        Eigen::Vector3d t;
    };

    // Solves for camera pose such that: lambda*x = R*X+t  with positive lambda.
    // Note: this code assumes that x is normalized.
    int p3p(const std::vector<Eigen::Vector3d> &x, const std::vector<Eigen::Vector3d> &X, std::vector<CameraPose> *output);

}