//
// Created by David Helekal on 2019-01-03.
//

#include "rbf_joint_discrete_kernel.hpp"
#include <cmath>
#include <iostream>

#ifdef _OPENMP

#endif


namespace mmdiff3 {

    rbf_joint_discrete_kernel::~rbf_joint_discrete_kernel() {
    }

    rbf_joint_discrete_kernel::rbf_joint_discrete_kernel(double *LUT, size_t maxval) {
        this->max_dist = maxval;
        this->lookup = LUT;
    }

    void rbf_joint_discrete_kernel::compute_LUT(size_t maxval, double sigma, double *lut) {

        size_t i = 0;
        double val;

#if defined(_OPENMP)
#pragma omp parallel for \
        default(shared) \
        private(i, val) \
        schedule(dynamic, 100)

        for (i = 0; i <= maxval; ++i) {
            assert(!std::isnan(sigma));
            val = exp((-1 / (2 * pow(sigma, 2))) * pow(i, 2));
            lut[i] = val;
        }
#else
        for(i=0; i <= maxval; ++i) {
            assert(!std::isnan(sigma));
            val = exp((-1/(2*pow(sigma,2))) * pow(i, 2));
            lut[i]=val;
        }
#endif
    }
}