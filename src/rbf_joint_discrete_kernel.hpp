//
// Created by David Helekal on 2019-01-03.
//

#ifndef MMDIFF3_rbf_joint_discrete_kernel_HPP
#define MMDIFF3_rbf_joint_discrete_kernel_HPP

#include <tuple>
#include <vector>
#include <cassert>
#include "kernel_function.hpp"

namespace mmdiff3 {
    class rbf_joint_discrete_kernel : public kernel_function<std::tuple<int, int>> {
    public:

        static void compute_LUT(size_t maxval, double sigma, double *lut);

        explicit rbf_joint_discrete_kernel(double *LUT, size_t maxval);

        ~rbf_joint_discrete_kernel();

        inline double compute_kernel(std::tuple<int, int> &a, std::tuple<int, int> &b) const {
            int pos_a;
            int cat_a;

            int pos_b;
            int cat_b;

            double result = 0.0;

            std::tie(pos_a, cat_a) = a;
            std::tie(pos_b, cat_b) = b;

            if (cat_a == cat_b) {
                size_t loc = abs(pos_a - pos_b);
                assert(loc <= this->max_dist);
                result = this->lookup[loc];
            }

            return result;
        }

    private:
        size_t max_dist;
        double *lookup;
    };
}


#endif //MMDIFF3_KERNEL_HPP
