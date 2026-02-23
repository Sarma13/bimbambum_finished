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

/*! \file BoundaryData.cpp
    \brief Base class containing function definition, properties and quantities associated with the fluid-fluid interface
    dynamics.
*/

#include "BoundaryData.hpp"

/*! Estimation of most appropriate time step based on the fluid-fluid interface dynamics.*/
double BoundaryData::time_step_boundary() {

    double u_nodes_max = u_nodes1[0];
    for (int i = 0; i < Ns + 1; ++i) {
        u_nodes_max = std::max(u_nodes_max, u_nodes1[i]);
    }

    return delta_phi / (1.0 + u_nodes_max * u_nodes_max);
}
