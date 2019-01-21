//
// Created by David Helekal on 23/10/2018.
//

#ifndef MMDIFF3_MMDIFF_3_CPP
#define MMDIFF3_MMDIFF_3_CPP

#import "mmdiff3.hpp"

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <tuple>
#include <vector>
#include "rbf_joint_kernel.hpp"
#include "mmd.hpp"

using namespace mmdiff3;

SEXP compute_jmmd(SEXP a1, SEXP a2, SEXP b1, SEXP b2, SEXP sigma) {
        double *ra1 = REAL(a1);
        int *ia2 = INTEGER(a2);

        double *rb1 = REAL(b1);
        int *ib2 = INTEGER(b2);

        const double rsigma = REAL(sigma)[0];
        auto ker = rbf_joint_kernel(rsigma);

        mmd<std::tuple<double, int> > run_mmd;
        std::vector<std::tuple<double, int>> vec1, vec2;

        for (int i = 0; i < length(a1); ++i) {
            vec1.emplace_back(std::make_tuple(ra1[i], ia2[i]));
        }

        for (int i = 0; i < length(b1); ++i) {
            vec2.emplace_back(std::make_tuple(rb1[i], ib2[i]));
        }

        SEXP ans = PROTECT(allocVector(REALSXP, 1));

        REAL(ans)[0] = run_mmd.compute_mmd(vec1, vec2, ker);
        UNPROTECT(1);
        return ans;
}
#endif //MMDIFF3_MMDIFF_3_HPP
