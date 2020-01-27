# lambdatwist
Re-implementation of Lambda Twist P3P solver from the paper

*M. Persson, K. Nordberg, Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver, ECCV 2018*

The original implementation from the authors can be found [here](https://github.com/midjji/lambdatwist-p3p). The main difference in this implementation is that it uses [Eigen](http://eigen.tuxfamily.org) and that it solves the cubic with a closed form solution. In our experiments this implementation has similar accuracy and speed as the original code.

The interface for the solver is (see [p3p.h](lambdatwist/p3p.h))
```C++
int p3p(const std::vector<Eigen::Vector3d> &x,
        const std::vector<Eigen::Vector3d> &X,
        std::vector<CameraPose> *output);
```
 where the camera is modeled as
```C++
struct CameraPose {
      Eigen::Matrix3d R;
      Eigen::Vector3d t;
};
```
**Note:** The solver assumes that the vectors in ```x``` are normalized. 
