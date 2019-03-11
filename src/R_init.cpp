//
// Created by David Helekal on 2019-01-16.
//
#include "R_init.hpp"

#include <R.h>
#include <R_ext/Rdynload.h>
#include "mmdiff3.hpp"

static const R_CallMethodDef callMethods[]  = {
        {"jmmd", (DL_FUNC) &compute_jmmd, 7},
        {NULL, NULL, 0}
};

void R_init_mmdiff3(DllInfo *info){
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}