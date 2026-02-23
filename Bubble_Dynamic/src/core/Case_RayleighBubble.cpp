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

/*! \file Case_RayleighBubble.cpp
    \brief Derived class: handles the dynamics of a bubble that has been initiated with the Rayleigh model.

    This derived class inherits members from BubbleData.cpp and BubbleData.cpp
*/

#include "Case_RayleighBubble.hpp"

#include <armadillo>
#include <cmath>
#include <vector>


using namespace std;
using namespace arma;

/*! Initialize the position of the bubble nodes and the value of the potentials at those nodes based on the Rayleigh model.*/
void Case_RayleighBubble::initialize() {

    const double pi = 3.14159265358979323846264338328;

    t0 = 0.0015527;
    R0 = 1;
    V0 = 0.0; // imposed at zero so that effect of non-condensible gas neglected in time_stepper.hpp
    V = 4.0 / 3.0 * pi * pow(R0, 3.0);

    for (int i = 0; i < Nb + 1; ++i) {
        double alpha_i = i * pi / Nb;
        r_nodes[i] = R0 * sin(alpha_i);
        z_nodes[i] = R0 * cos(alpha_i) - gamma;
        ur_nodes[i] = 0.0; // initialized at 0, not needed for first step computation
        uz_nodes[i] = 0.0; // initialized at 0, not needed for first step computation
        phi_nodes[i] = -R0 * sqrt((2.0 / 3.0 * (pow(1.0 / R0, 3.0) - 1))); // for R0 = 0.1, phi_0 = -2.5807
    }

    // impose first and last node to be exactly on the axis of symmetry
    r_nodes.front() = 0.0;
    r_nodes.back() = 0.0;

}

