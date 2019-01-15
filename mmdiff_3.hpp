//
// Created by David Helekal on 23/10/2018.
//

#ifndef MMDIFF3_MMDIFF_3_HPP
#define MMDIFF3_MMDIFF_3_HPP

#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <tuple>

namespace mmdiff3 {
    extern "C" {
    class mmdiff3 {
    public:
        SEXP compute_jmmd(SEXP a1, SEXP a2, SEXP b1, SEXP b2, SEXP sigma);
    };

    }
}
#endif //MMDIFF3_MMDIFF_3_HPP
