//  my_function.cpp
#include <cmath>
#include <iostream>
// #include "ellint.h"
extern "C"{
    void ellint(double *k,double *F,double *E){
        // double hpi = std::acos(-1)/2;
        // std::cout << *input << '\n';
        // std::cout << "K(0) = " << std::comp_ellint_1f(0) << '\n'
        //         << "π/2 = " << hpi << '\n'
        //         << "K(0.5) = " << std::comp_ellint_1f(*input) << '\n'
        //         << "F(0.5, π/2) = " << std::ellint_1f(*input, hpi) << '\n';
        // std::cout << "Period of a pendulum length 1 m at 90 degree initial angle is "
        //         << 4*std::sqrt(1/9.80665)*
        //             std::comp_ellint_1f(std::pow(std::sin(hpi/2),2)) << " s\n";
        *F = std::comp_ellint_1f(*k);
        *E = std::comp_ellint_2f(*k);
    };
}
