//
// Created by David Helekal on 2019-01-16.
//
#include "R_init.hpp"

#include "mmdiff3.hpp"

static const R_CallMethodDef callMethods[] = {
        {"jmmd",                 (DL_FUNC) &compute_jmmd,                       6},
        {"kernel_lut",           (DL_FUNC) &compute_LUT,                        2},
        {"kernel_sum",           (DL_FUNC) &compute_joint_kernel_sum,           6},
        {"kernel_sum_symmetric", (DL_FUNC) &compute_joint_kernel_sum_symmetric, 4},
        {NULL, NULL,                                                            0}
};

void R_init_mmdiff3(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}