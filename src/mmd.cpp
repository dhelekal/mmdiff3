//
// Created by David Helekal on 2019-01-07.
//

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mmd.hpp"
#include <tuple>
#include <cmath>
#include <iostream>
#include <stdlib.h>

namespace mmdiff3 {
    template<class T>
    double mmd<T>::compute_mmd(std::vector<T> &x, std::vector<T> &y, kernel_function<T> &k) {
        double p = kernel_sum(x, x, k, true), q = kernel_sum(y, y, k, true), pq = kernel_sum(x, y, k, false);
        size_t m = x.size();
        size_t n = y.size();
        return (1. / (m * (m-1))) * p - (2. / (m * n)) * pq + (1. / (n * (n-1))) * q;
    }

    template<class T>
    mmd<T>::mmd() {

    }

    template<class T>
    mmd<T>::~mmd() {
    }

    template<class T>
    double mmd<T>::kernel_sum(std::vector<T> &x,
            std::vector<T> &y,
            kernel_function<T> &k,
            bool no_diag) {
        size_t m = x.size();
        size_t n = y.size();

        double result = 0;
        size_t i = 0;
        size_t j = 0;

#if defined(_OPENMP)
#pragma omp parallel for collapse(2) \
        default(shared) \
        private(i,j) \
        schedule(dynamic, 100) \
        reduction(+:result)

        for(i=0; i < m; ++i){
            for(j=0; j < n; ++j){
                double k_res;
                (no_diag && i==j) ? (k_res = 0) : (k_res = k.compute_kernel(x[i], y[j]));
                result = result + k_res;
            }
        }
#else
        for(i=0; i < m; ++i){
            for(j=0; j < n; ++j){
                double k_res;
                (no_diag && i==j) ? (k_res = 0) : (k_res = k.compute_kernel(x[i], y[j]));
                result = result + k_res;
            }
        }
#endif
        if(std::isnan(result)){
            std::cout << "NaN Encountered";
        }

        return result;
    }

    template class mmd<std::tuple<double, int>>;
}

