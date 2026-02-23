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

/*! \file time_stepper.hpp
    \brief  Conducts one time step.

    Updates in time the position (r, z) and potentials (phi, F) of the nodes
    on the bubble and fluid-fluid interfaces.
    The RK1 and RK2 time stepping schemes are implemented to allow the user
    to choose between fast computations (RK1, one stage) and computations with
    an increased stability and accuracy (RK2, two stages).
*/

#ifndef TIME_STEPPER_HPP
#define TIME_STEPPER_HPP

#include <algorithm>
#include <cmath>
#include "Case_SphericalBoundary.hpp"

template<typename BubbleData, typename BoundaryData, typename Inputs, typename BIM_solver>

void time_integration(BubbleData &bubble, BoundaryData &boundary, Inputs &data, BIM_solver &step) {

    // ---------- Heun's method, 2 stages ----------
    if (data.temporal_solver == "RK2") {

        double dt_b; // adaptive time step based on bubble dynamics
        double dt_s; // adaptive time step based on fluid-fluid interface dynamics
        double dt; // adaptive time step considering both dynamics

        //! Adding driving force frequency
        
        const double Omega = data.Omega; // driving force frequency
        const double A = data.A; // driving force amplitude
        

        // first stage  

        boundary->boundary_endpoints_derivatives(); // computes the spatial derivatives
        // at the fluid-fluid interface endpoints

        if (data.surface_elasticity == true) {
            boundary->boundary_curvature();
        }

        bubble->V = bubble->compute_volume();
        
        // Compute spherical volume if using a spherical boundary
        if (data.boundary == "spherical") {
            auto* sphericalBoundary = dynamic_cast<Case_SphericalBoundary*>(boundary.get());
            if (sphericalBoundary) {
                sphericalBoundary->V_spherical = sphericalBoundary->compute_spherical_volume();
            }
        }

        step.compute_un(bubble, boundary, data);
        step.compute_ut(bubble, boundary, data);

        dt_b = bubble->time_step_bubble(data.epsilon, data.k);
        dt_s = boundary->time_step_boundary();
        dt = fmin(dt_b, dt_s); // takes minimum between bubble and fluid-fluid interface time step
        dt = fmax(dt, 0.00001); // impose a minimum time step value to avoid over-filtering
        dt = fmin(dt, 0.01); // impose a maximum time step value for stability
        step.dt = dt;
        // **ADDED**: driving for Stage 1 (time = step.time)
        
        
        for (int i = 0; i < data.Nb + 1; ++i) {
            
            // const double A = 1;
            // const double Omega = 10;
            // const double driving = A * std::sin(Omega * step.time);
            double drivingStage1 = A * std::sin(Omega * step.time);
            bubble->dr1[i] = dt * bubble->ur_nodes[i];
            bubble->dz1[i] = dt * bubble->uz_nodes[i];
            bubble->dphi1[i] = dt * (1 + drivingStage1 + 0.5 * (bubble->u_nodes[i] * bubble->u_nodes[i]) -
                                     data.epsilon * pow(bubble->V0 / bubble->V, data.k) -
                                     data.zeta * (bubble->z_nodes[i] + data.gamma));
            bubble->r_nodes1[i] = bubble->r_nodes[i];
            bubble->z_nodes1[i] = bubble->z_nodes[i];
            bubble->phi_nodes1[i] = bubble->phi_nodes[i];
            bubble->r_nodes[i] = bubble->r_nodes[i] + bubble->dr1[i];
            bubble->z_nodes[i] = bubble->z_nodes[i] + bubble->dz1[i];
            bubble->phi_nodes[i] = bubble->phi_nodes[i] + bubble->dphi1[i];
        }
        for (int i = 0; i < data.Ns + 1; ++i) {

            boundary->dr1[i] = dt * boundary->ur_nodes[i];
            boundary->dz1[i] = dt * boundary->uz_nodes[i];
            boundary->dF1[i] = dt * (boundary->u_vec_prod[i] - 0.5 * (boundary->u_nodes2[i] * boundary->u_nodes2[i] +
                                                                      data.alpha * boundary->u_nodes1[i] *
                                                                      boundary->u_nodes1[i]) +
                                     data.sigma_s * boundary->curv_nodes[i] -
                                     data.zeta * (1.0 - data.alpha) * boundary->z_nodes[i]);
            boundary->r_nodes1[i] = boundary->r_nodes[i];
            boundary->z_nodes1[i] = boundary->z_nodes[i];
            boundary->F_nodes1[i] = boundary->F_nodes[i];
            boundary->r_nodes[i] = boundary->r_nodes[i] + boundary->dr1[i];
            boundary->z_nodes[i] = boundary->z_nodes[i] + boundary->dz1[i];
            boundary->F_nodes[i] = boundary->F_nodes[i] + boundary->dF1[i];
        }



        //Second stage
        boundary->r_nodes[0]      = 0.0;
        boundary->r_nodes[data.Ns] = 0.0;
        boundary->boundary_endpoints_derivatives();
        if (data.surface_elasticity == true) {
            boundary->boundary_curvature();
        }
        bubble->V = bubble->compute_volume();
        
        // Compute spherical volume if using a spherical boundary
        if (data.boundary == "spherical") {
            auto* sphericalBoundary = dynamic_cast<Case_SphericalBoundary*>(boundary.get());
            if (sphericalBoundary) {
                sphericalBoundary->V_spherical = sphericalBoundary->compute_spherical_volume();
            }
        }
        
        step.compute_un(bubble, boundary, data);
        step.compute_ut(bubble, boundary, data);

        for (int i = 0; i < data.Nb + 1; ++i) {
            
            
            double drivingStage2 = data.A * std::sin(data.Omega * (step.time + dt));

            bubble->dr2[i] = dt * bubble->ur_nodes[i];
            bubble->dz2[i] = dt * bubble->uz_nodes[i];
            bubble->dphi2[i] = dt * (1.0 + drivingStage2 + 0.5 * (bubble->u_nodes[i] * bubble->u_nodes[i]) -
                                     data.epsilon * pow(bubble->V0 / bubble->V, data.k) -
                                     data.zeta * (bubble->z_nodes[i] + data.gamma));
            bubble->r_nodes[i] = bubble->r_nodes1[i] + 0.5 * (bubble->dr1[i] + bubble->dr2[i]);
            bubble->z_nodes[i] = bubble->z_nodes1[i] + 0.5 * (bubble->dz1[i] + bubble->dz2[i]);
            bubble->phi_nodes[i] = bubble->phi_nodes1[i] + 0.5 * (bubble->dphi1[i] + bubble->dphi2[i]);
        }

        for (int i = 0; i < data.Ns + 1; ++i) {
            
            boundary->dr2[i] = dt * boundary->ur_nodes[i];
            boundary->dz2[i] = dt * boundary->uz_nodes[i];
            boundary->dF2[i] = dt * (boundary->u_vec_prod[i] - 0.5 * (boundary->u_nodes2[i] * boundary->u_nodes2[i] +
                                                                      data.alpha * boundary->u_nodes1[i] *
                                                                      boundary->u_nodes1[i]) +
                                     data.sigma_s * boundary->curv_nodes[i] -
                                     data.zeta * (1.0 - data.alpha) * boundary->z_nodes[i]);
            boundary->r_nodes[i] = boundary->r_nodes1[i] + 0.5 * (boundary->dr1[i] + boundary->dr2[i]);
            boundary->z_nodes[i] = boundary->z_nodes1[i] + 0.5 * (boundary->dz1[i] + boundary->dz2[i]);
            boundary->F_nodes[i] = boundary->F_nodes1[i] + 0.5 * (boundary->dF1[i] + boundary->dF2[i]);
        }
    }

        // ---------- Euler's method, 1stage ----------
    else if (data.temporal_solver == "RK1") {

        double dt_b;
        double dt_s;
        double dt;

        // first stage
        boundary->boundary_endpoints_derivatives(); // computes the spatial derivatives
        // at the fluid-fluid interface endpoints

        if (data.surface_elasticity == true) {
            boundary->boundary_curvature();
        }

        bubble->V = bubble->compute_volume();
        
        // Compute spherical volume if using a spherical boundary
        if (data.boundary == "spherical") {
            auto* sphericalBoundary = dynamic_cast<Case_SphericalBoundary*>(boundary.get());
            if (sphericalBoundary) {
                sphericalBoundary->V_spherical = sphericalBoundary->compute_spherical_volume();
            }
        }
        
        step.compute_un(bubble, boundary, data);
        step.compute_ut(bubble, boundary, data);

        dt_b = bubble->time_step_bubble(data.epsilon, data.k);
        dt_s = boundary->time_step_boundary();
        dt = fmin(dt_b, dt_s); // takes minimum between bubble and fluid-fluid interface time step
        dt = fmax(dt, 0.00001); // impose a minimum time step value to avoid over-filtering
        dt = fmin(dt, 0.01); // impose a maximum time step value for stability
        step.dt = dt;

        for (int i = 0; i < data.Nb + 1; ++i) {
            const double A = data.A;
            const double Omega = data.Omega;
            const double driving = A * std::sin(Omega * step.time);
            bubble->dr1[i] = dt * bubble->ur_nodes[i];
            bubble->dz1[i] = dt * bubble->uz_nodes[i];
            bubble->dphi1[i] = dt * (1.0 + driving + 0.5 * (bubble->u_nodes[i] * bubble->u_nodes[i]) -
                                     data.epsilon * pow(bubble->V0 / bubble->V, data.k) -
                                     data.zeta * (bubble->z_nodes[i] + data.gamma));
            bubble->r_nodes[i] = bubble->r_nodes[i] + bubble->dr1[i];
            bubble->z_nodes[i] = bubble->z_nodes[i] + bubble->dz1[i];
            bubble->phi_nodes[i] = bubble->phi_nodes[i] + bubble->dphi1[i];
        }
        for (int i = 0; i < data.Ns + 1; ++i) {
            boundary->dr1[i] = dt * boundary->ur_nodes[i];
            boundary->dz1[i] = dt * boundary->uz_nodes[i];
            boundary->dF1[i] = dt * (boundary->u_vec_prod[i] - 0.5 * (boundary->u_nodes2[i] * boundary->u_nodes2[i] +
                                                                      data.alpha * boundary->u_nodes1[i] *
                                                                      boundary->u_nodes1[i]) +
                                     data.sigma_s * boundary->curv_nodes[i] -
                                     data.zeta * (1.0 - data.alpha) * boundary->z_nodes[i]);
            boundary->r_nodes[i] = boundary->r_nodes[i] + boundary->dr1[i];
            boundary->z_nodes[i] = boundary->z_nodes[i] + boundary->dz1[i];
            boundary->F_nodes[i] = boundary->F_nodes[i] + boundary->dF1[i];
        }

    } else {
        std::cerr << "Time stepper unavailable: choose between RK1 and RK2" << std::endl;
        exit(EXIT_FAILURE);
    }

};

#endif // TIME_STEPPER_HPP
