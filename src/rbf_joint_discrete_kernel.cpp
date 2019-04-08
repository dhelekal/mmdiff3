//
// Created by David Helekal on 2019-01-03.
//

#include "rbf_joint_discrete_kernel.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace mmdiff3 {

    rbf_joint_discrete_kernel::~rbf_joint_discrete_kernel() {
        free(this->lookup);
    }

    rbf_joint_discrete_kernel::rbf_joint_discrete_kernel(double sigma, size_t min_b, size_t max_b) {
        this->sigma = sigma;
        this->max_dist = abs((long)(max_b-min_b));
        this->lookup = (double*) malloc(sizeof(double[max_dist+1]));

        size_t i = 0;
        double val;

#if defined(_OPENMP)
#pragma omp parallel for \
        default(shared) \
        private(i, val) \
        schedule(dynamic, 100)

        for(i=0; i <= max_dist; ++i) {
            assert(!std::isnan(sigma));
            val = exp((-1/(2*pow(sigma,2))) * pow(i, 2));
            this->lookup[i]=val;
        }
#else
        for(i=0; i <= max_dist; ++i) {
            assert(!std::isnan(sigma));
            auto val = exp((-1/(2*pow(sigma,2))) * pow(i, 2));
            this->lookup[i]=val;
        }
#endif

    }

    inline double rbf_joint_discrete_kernel::compute_kernel(std::tuple<int, int> &a, std::tuple<int, int> &b) const {
        int pos_a;
        int cat_a;

        int pos_b;
        int cat_b;

        double result = 0.0;

        std::tie(pos_a, cat_a) = a;
        std::tie(pos_b, cat_b) = b;

        size_t loc = abs(pos_a - pos_b);
        assert(loc <= this->max_dist);

        if (cat_a == cat_b) result = this->lookup[loc];

        return result;
    }
}