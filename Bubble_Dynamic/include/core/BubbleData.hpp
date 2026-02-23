/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
 *
 * BubbleData.hpp
 *
 * Description:
 * Base class containing function declaration, properties and quantities
 * associated with the bubble dynamics.
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

/*! \file BubbleData.hpp
    \brief Header for BubbleData.cpp
*/

/*! \class BubbleData
    \brief Base class containing function declaration, properties and quantities associated with the bubble dynamics.
*/

#ifndef BUBBLEDATA_HPP
#define BUBBLEDATA_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "init_R0.hpp"

using namespace std;


class BubbleData {
public://!CONSTRUCTOR
    template<typename Input>
    BubbleData(Input &data) {

        Nb = data.Nb;
        r_nodes.resize(Nb + 1);
        z_nodes.resize(Nb + 1);
        ur_nodes.resize(Nb + 1);
        uz_nodes.resize(Nb + 1);
        u_nodes.resize(Nb + 1);
        phi_nodes.resize(Nb + 1);

        r_nodes1.resize(Nb + 1);
        z_nodes1.resize(Nb + 1);
        phi_nodes1.resize(Nb + 1);
        dr1.resize(Nb + 1);
        dz1.resize(Nb + 1);
        dphi1.resize(Nb + 1);
        dr2.resize(Nb + 1);
        dz2.resize(Nb + 1);
        dphi2.resize(Nb + 1);

        zeta = data.zeta;
        gamma = data.gamma;
        epsilon = data.epsilon;
        k = data.k;
        Rmax = data.Rmax;
        e = data.e;
      

        bubble_dynamics = data.bubble_dynamics;
        delta_phi = data.delta_phi;

    }

    virtual ~BubbleData() = default;


    //Input data;
    int Nb; /*!< \brief The number of elements for the bubble surface discretization.  */
    std::vector<double> r_nodes; /*!< \brief r-coordinate of the nodes on the bubble surface.*/
    std::vector<double> z_nodes; /*!< \brief z-coordinate of the nodes on the bubble surface.*/
    std::vector<double> ur_nodes; /*!< \brief r-component of the velocity of the bubble surface nodes.*/
    std::vector<double> uz_nodes; /*!< \brief z-component of the velocity of the bubble surface nodes.*/
    std::vector<double> u_nodes; /*!< \brief Velocity magnitude of the bubble surface nodes.*/
    std::vector<double> phi_nodes; /*!< \brief Potential values at the bubble surface nodes.*/

    // Time integration supplement storage
    std::vector<double> r_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> z_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> phi_nodes1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dr1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dz1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dphi1; /*!< \brief Additional vector for storage.*/
    std::vector<double> dr2; /*!< \brief Additional vector for storage.*/
    std::vector<double> dz2; /*!< \brief Additional vector for storage.*/
    std::vector<double> dphi2; /*!< \brief Additional vector for storage.*/

    //! MEMBERS VARIABLES
    // Bubble properties
    double R0; /*!< \brief Initial bubble radius.*/
    double V0; /*!< \brief Initial bubble volume.*/
    double V; /*!< \brief Current bubble volume.*/
    double V_spherical; /*!< \brief The spherical droplet volume.*/
    double zeta; /*!< \brief Buoyancy parameter.*/
    double gamma; /*!< \brief Stand-off disatnce.*/
    double epsilon; /*!< \brief Strength parameter.*/
    double k; /*!< \brief Specific heat ratio.*/
    double Rmax; /*!< \brief Maximum bubble radius.*/
    double e; /*!< \brief eccentric.*/
    // Simulation case properties
    std::string bubble_dynamics; /*!< \brief Physics governing the bubble behaviour (Rayleigh or Rayleigh-Plesset).*/
    double delta_phi; /*!< \brief Constant for adaptive time stepping.*/
    double t0; /*!< \brief Initial simulation time.*/

    /*! \brief Initializes the node points on the bubble surface.*/
    virtual void initialize() = 0;

    /*! \brief Filters the bubble surface. */
    void filter_bubble();

    /*! \brief Redistributes the node points position along the bubble surface. */
    void remesh_bubble();

    /*! \brief Computes the bubble volume. */
    double compute_volume();

    /*! \brief Computes the volume of a spherical droplet. */
    double compute_spherical_volume();

    /*! \brief Checks the bubble intersection (transition form singly-connected problem to doubly connected problem --> end of simulation). */
    int intersect();

    /*! \brief Computes the time step value based on the bubble's dynamics*/
    double time_step_bubble(double epsilon, double k);

};

#endif // BUBBLEDATA_HPP
