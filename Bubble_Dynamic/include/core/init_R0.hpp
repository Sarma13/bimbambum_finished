/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * -----------------------------------------------------------------------------
 * Copyright (C) 2023 Armand Sieber
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file init_R0.hpp
    \brief  Estimates the initial bubble radius, R0, in the case of the bubble dynamics
           initiated with the Rayleigh-Plesset model.

   The value of R0 is found using a root-finding algorithm based on the inputs epsilon and k,
   provided by the user. The root-finding algorithm is provided by the GNU scientific library.
*/

#ifndef INIT_R0_HPP
#define INIT_R0_HPP

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


struct root_params {
    double epsilon;
    double k;
};

inline double root_f(double x, void *p) {
    struct root_params *params = (struct root_params *) p;
    double epsilon = (params->epsilon);
    double k = (params->k);
    double func = (epsilon / (k - 1.0)) * (pow(x, 3.0 * k) - x * x * x) + (1.0 - x * x * x);
    return func;
}

inline double set_R0(double epsilon, double k) {

    int status;
    int iter = 0;
    int max_iter = 100;
    double R_init;
    double x_low = 0.0;
    double x_high = 0.8;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;


    gsl_function F;
    struct root_params params = {epsilon, k};
    F.function = &root_f;
    F.params = &params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, x_low, x_high);

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        R_init = gsl_root_fsolver_root(s);
        x_low = gsl_root_fsolver_x_lower(s);
        x_high = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_low, x_high, 0, 1.0e-8);

    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    return R_init;
}

#endif // INIT_R0_HPP