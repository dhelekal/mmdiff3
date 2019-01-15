//
// Created by David Helekal on 2019-01-07.
//

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mmd.hpp"
#include <tuple>

namespace mmdiff3 {
    template<class T>
    double mmd<T>::compute_mmd(std::vector<T> &x, std::vector<T> &y, kernel_function<T> &k) {
        double p = kernel_sum(x, x, k), q = kernel_sum(y, y, k), pq = kernel_sum(x, y, k);
        size_t m = x.size();
        size_t n = y.size();
        return (1. / (m * m)) * p - (2. / (m * n)) * pq + (1. / (n * n)) * q;
    }

    template<class T>
    mmd<T>::mmd() {

    }

    template<class T>
    mmd<T>::~mmd() {
    }

    template<class T>
    double mmd<T>::kernel_sum(std::vector<T> &x, std::vector<T> &y, kernel_function<T> &k) {
        size_t m = x.size();
        size_t n = y.size();

        double result = 0;
        int i = 0;

#if defined(_OPENMP)
#pragma omp parallel for \
        default(shared) private(i) \
        schedule(static, 100) \
        reduction(+:result)

        for (i = 0; i < (m * n); ++i) {
            result = result + k.compute_kernel(x[i / n], y[i % m]);
        }
#else
        for (i = 0; i < (m * n); ++i) {
                result = result + k.compute_kernel(x[i / n], y[i % m]);
            }
#endif
        return result;
    }

    template class mmd<std::tuple<double, int>>;
}

