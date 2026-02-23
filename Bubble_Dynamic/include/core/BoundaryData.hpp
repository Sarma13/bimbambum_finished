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

/*! \file BoundaryData.hpp
    \brief Header for BoundaryData.cpp
*/

/*! \class BoundaryData
    \brief Base class containing function declaration, properties and quantities associated with the fluid-fluid
    interface.
*/

#ifndef BOUNDARYDATA_HPP
#define BOUNDARYDATA_HPP

#include <vector>
#include <string>

class BoundaryData {
public:
    template<typename Input>
    BoundaryData(Input &data) {

        Ns = data.Ns;
        delta_phi = data.delta_phi;

        gamma = data.gamma;
        alpha = data.alpha;
        xi = data.xi;

        r_nodes.resize(Ns + 1);
        z_nodes.resize(Ns + 1);
        ur_nodes.resize(Ns + 1);
        uz_nodes.resize(Ns + 1);
        u_nodes1.resize(Ns + 1);
        u_nodes2.resize(Ns + 1);
        u_vec_prod.resize(Ns + 1);
        F_nodes.resize(Ns + 1);
        curv_nodes.resize(Ns + 1);

        r_nodes1.resize(Ns + 1);
        z_nodes1.resize(Ns + 1);
        F_nodes1.resize(Ns + 1);
        dr1.resize(Ns + 1);
        dz1.resize(Ns + 1);
        dF1.resize(Ns + 1);
        dr2.resize(Ns + 1);
        dz2.resize(Ns + 1);
        dF2.resize(Ns + 1);


    }

    virtual ~BoundaryData() = default;


    //Input data;
    int Ns; /*!< \brief The number of elements for the liquid-liquid interface discretization.  */
    double delta_phi; /*!< \brief Constant for adaptive time stepping.*/

    // Bubble and interface properties
    double R0; /*!< \brief Initial bubble radius.*/
    double gamma; /*!< \brief stand-off distance.*/
    double alpha; /*!< \brief Density ratio.*/
    double xi; /*!< \brief Strength parameter.*/
    // Solution vectors
    std::vector<double> r_nodes; /*!< \brief r-coordinate of the nodes on the fluid-fluid interface.*/
    std::vector<double> z_nodes; /*!< \brief z-coordinate of the nodes on the fluid-fluid interface.*/
    std::vector<double> ur_nodes; /*!< \brief r-component of the velocity of the fluid-fluid interface nodes.*/
    std::vector<double> uz_nodes; /*!< \brief z-component of the velocity of the fluid-fluid interface nodes.*/
    std::vector<double> u_nodes1; /*!< \brief Velocity magnitude of the fluid-fluid interface nodes (fluid 1).*/
    std::vector<double> u_nodes2; /*!< \brief Velocity magnitude of the fluid-fluid interface nodes (fluid 2).*/
    std::vector<double> u_vec_prod; /*!< \brief Dot product of the velocity vectors on each side of the fluid-fluid interface.*/
    std::vector<double> F_nodes; /*!< \brief Potential difference across the fluid-fluid interface.*/
    std::vector<double> curv_nodes; /*!< \brief Curvature on the fluid-fluid interface.*/

    // Endpoints derivatives (for cubic spline interpolation)
    double drds1; /*!< \brief Endpoint derivative (the extremum of the domain): dr/ds .*/
    double drds2; /*!< \brief Endpoint derivative (at the axis of symmetry): dr/ds.*/
    double dzds1; /*!< \brief Endpoint derivative (the extremum of the domain): dz/ds.*/
    double dzds2; /*!< \brief Endpoint derivative (at the axis of symmetry): dz/ds.*/

    double dphi1ds1; /*!< \brief Endpoint derivative (the extremum of the domain): dph1/ds.*/
    double dphi1ds2; /*!< \brief Endpoint derivative (at the axis of symmetry): dph1/ds.*/
    double dphi2ds1; /*!< \brief Endpoint derivative (the extremum of the domain): dph2/ds.*/
    double dphi2ds2; /*!< \brief Endpoint derivative (at the axis of symmetry): dph2/ds.*/

    // Time integration supplement storage
    std::vector<double> r_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> z_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> F_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dr1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dz1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dF1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dr2; /*!< \brief Additional vector for storage.*/
    std::vector<double> dz2; /*!< \brief Additional vector for storage.*/
    std::vector<double> dF2; /*!< \brief Additional vector for storage.*/

    /*! \brief Initializes fluid-fluid interface node points.*/
    virtual void initialize() = 0;

    /*! \brief Redistributes the node points position along fluid-fluid interface. */
    virtual void remesh_boundary() = 0;

    /*! \brief Filters the fluid-fluid interface. */
    virtual void filter_boundary() = 0;

    /*! \brief Estimates the derivatives at the endpoints on the fluid-fluid interface. */
    virtual void boundary_endpoints_derivatives() = 0;

    /*! \brief Computes the curvature on the fluid-fluid interface. */
    virtual void boundary_curvature() = 0;

    /*! \brief Computes the time step value based on the fluid-fluid interface dynamics*/
    double time_step_boundary();

};


#endif //BOUNDARYDATA_H
