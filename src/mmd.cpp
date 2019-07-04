//
// Created by David Helekal on 2019-01-07.
//

#if defined(_OPENMP)

#include <omp.h>

#endif

#include "mmd.hpp"
#include <cmath>
#include <iostream>
#include <cassert>

namespace mmdiff3 {
    template<class T>
    double mmd<T>::compute_mmd(std::vector<T> &x, std::vector<T> &y, kernel_function<T> &k) {
        double p = kernel_sum(x, x, k, false), q = kernel_sum(y, y, k, false), pq = kernel_sum(x, y, k, false);
        size_t m = x.size();
        size_t n = y.size();
        double res = sqrt(((1. / (m * (m))) * p) - ((2. / (m * n)) * pq) + ((1. / (n * (n))) * q));

        assert(!std::isnan(res));
        assert(std::isfinite(res));
        assert(res >= 0);

        return res;
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
        double k_res = 0;
        size_t i = 0;
        size_t j = 0;

#if defined(_OPENMP)
#pragma omp parallel for collapse(2) \
        default(shared) \
        private(i, j, k_res) \
        schedule(dynamic, 1000) \
        reduction(+:result)

        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
                if (no_diag && i == j) {
                    k_res = 0.0;
                } else {
                    k_res = k.compute_kernel(x[i], y[j]);
                }

                assert(!std::isnan(k_res));
                assert(k_res <= 1);

                result += k_res;
            }
        }
#else
        for(i=0; i < m; ++i){
            for(j=0; j < n; ++j){
                if (no_diag && i==j) {
                    k_res = 0.0;
                } else {
                    k_res = k.compute_kernel(x[i], y[j]);
                }

                assert(!std::isnan(k_res));
                result = result + k_res;
            }
        }
#endif
        if (std::isnan(result)) {
            std::cout << "NaN Encountered";
        }

        return result;
    }

    template<class T>
    double mmd<T>::kernel_sum_symmetric(std::vector<T> &x,
                                        kernel_function<T> &k,
                                        bool no_diag) {

        size_t m = x.size();

        double diag = 0;
        for (size_t l = 0; l < m; ++l) {
            diag += k.compute_kernel(x[l], x[l]);
        }

        double result = 0;
        double partial_result = 0;
        double k_res = 0;
        size_t i = 0;
        size_t j = 0;

#if defined(_OPENMP)
#pragma omp parallel for \
        default(shared) \
        private(i, j, k_res, partial_result) \
        schedule(dynamic, 30) \
        reduction(+:result)

        for (i = 0; i < m; ++i) {

            partial_result = 0;
            for (j = 0; j < i; ++j) {
                if (no_diag && i == j) {
                    k_res = 0.0;
                } else {
                    k_res = k.compute_kernel(x[i], x[j]);
                }

                assert(!std::isnan(k_res));
                assert(k_res <= 1);
                partial_result += k_res;
            }

            assert(!std::isnan(partial_result));
            result += partial_result;
        }
#else
        for(i=0; i < m; ++i){

            partial_result = 0;
            for(j=0; j < m; ++j){
                if (no_diag && i==j) {
                    k_res = 0.0;
                } else {
                    k_res = k.compute_kernel(x[i], x[j]);
                }

                assert(!std::isnan(k_res));
                assert(k_res <= 1);
                partial_result += k_res;
            }
            result += partial_result;
        }
#endif
        if (std::isnan(result)) {
            std::cout << "NaN Encountered";
        }

        return 2 * result + diag;
    }

    template
    class mmd<std::tuple<double, int>>;

    template
    class mmd<std::tuple<int, int>>;
}

