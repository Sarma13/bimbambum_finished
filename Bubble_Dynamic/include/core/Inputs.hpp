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

/*! \file Inputs.hpp
    \brief Class handling input from configuration file.
*/
/*! \class Input
    \brief Class handling input from configuration file.

    Input interface from the configuration file parser to
    the program.
*/

#ifndef INPUTS_HPP
#define INPUTS_HPP

#include <iostream>
#include <ConfigFileParser.hpp>


class Input {
    public:

        // ---------- Simulation Parameters ----------
        // -------------------------------------------
        /*! \brief The number of elements for the bubble discretization. */
        int Nb;

        /*! \brief The number of elements for the liquid-liquid interface discretization. */
        int Ns;

        /*! \brief The physics governing the bubble.
        The physics takes the values either 'Rayleigh_Bubble' or 'Rayleigh_Plesset_Bubble'.
        */
        std::string bubble_dynamics;

        /*! \brief Buoyancy parameter */
        double zeta;

        /*! \brief Stand-off distance */
        double gamma;

        /*! \brief Density ratio: \f$ \alpha = \rho_1/\rho_2 \f$, with
        \f$ \rho_1\f$ the density of the liquid in which the bubble is located.
        */
        double alpha;

        /*! \brief The strength parameter
        The strength parameter is defined as \f$\epsilon = p_0/ \Delta p\f$ with \f$p_0\f$ the gas pressure
        within the bubble and \f$ \Delta p = p_{inf}-p_{v}\f$.
        */
        double epsilon;

        /*! \brief Specific heat ratio (for air \f$k = 1.4 \f$) */
        double k;

        /*! \brief Whether or not the interface has elasticity */
        bool surface_elasticity;

        /*! \brief The interface surface tension */
        double sigma_s;

        /*! \brief The nature of the initial geometry of the nearby boundary
         N.B. In the current implementation the boundary is initialized as flat and the keyword 'from_code' must be used.
         */
        std::string boundary; //boundary type: "from_code" or "from_file"

        // ---------- Solver Parameters ----------
        // ---------------------------------------

        /*! \brief The order of the temporal discretization.
        The values are either 'RK1' for 1st order or 'RK2' for second order.
         */
        std::string temporal_solver; // 'RK1' or 'RK2'

        /*! \brief Name of the dumper file (stores the solution results). A 'txt' extension is recommended.*/
        std::string dumper_filename;

        /*! \brief Constant for adaptive time stepping.*/
        double delta_phi;

        /*! \brief Filtering frequency (i.e. filtering every nth time steps).*/
        int filtering_freq;

        /*! \brief Number of threads used for parallel processing.*/
        int n_threads;


        //*---------- NEW PARAMETERS ----------
        double Omega;
        double A;
        double xi;
        double Rmax;
        double e;






        /*! \brief Initialize the input data based on the values provided by the user in the config file.*/
        void initialize(std::string input_conf, int argc, const char **argv) {

            ConfigFileParser parser(input_conf);
            parser.PrintConfigFile();

           

            //Simulation parameters
            Nb = parser.GetConfigValueFromList<int>("Bubble data", "Nb", 50 );
            Ns = parser.GetConfigValueFromList<int>("Bubble data", "Ns", 50 );
            bubble_dynamics = parser.GetConfigValueFromList<std::string>("Bubble data", "bubble_dynamics", "Rayleigh_Bubble" );
            zeta = parser.GetConfigValueFromList<double>("Bubble data", "zeta", 0.0 );
            gamma = parser.GetConfigValueFromList<double>("Bubble data", "gamma", 1.0 );
            alpha = parser.GetConfigValueFromList<double>("Bubble data", "alpha", 1.0 );
            epsilon = parser.GetConfigValueFromList<double>("Bubble data", "epsilon", 100.0 );
            k = parser.GetConfigValueFromList<double>("Bubble data", "k", 1.4 );
            surface_elasticity = parser.GetConfigValueFromList<bool>("Bubble data", "surface_elasticity", false );
            sigma_s = parser.GetConfigValueFromList<double>("Bubble data", "sigma_s", 0.1 );
            boundary = parser.GetConfigValueFromList<std::string>("Bubble data", "boundary", "spherical" );
            //* NEW
            Omega = parser.GetConfigValueFromList<double>("Bubble data", "Omega", /*your default*/ 0.0);
            A = parser.GetConfigValueFromList<double>("Bubble data", "A", /*your default*/ 0.0);
            xi = parser.GetConfigValueFromList<double>("Bubble data", "xi", /*your default*/ 2);
            Rmax = parser.GetConfigValueFromList<double>("Bubble data", "Rmax", /*your default*/ 1.0);
            e = parser.GetConfigValueFromList<double>("Bubble data", "e", /*your default*/ 0.0);
            //Solver parameters
            temporal_solver = parser.GetConfigValueFromList<std::string>("Solver parameters", "temporal_solver", "RK1" );
            dumper_filename = parser.GetConfigValueFromList<std::string>("Solver parameters", "dumper_filename", "nodes_position.txt" );
            delta_phi = parser.GetConfigValueFromList<double>("Solver parameters", "delta_phi", 2.0e-2);
            filtering_freq = parser.GetConfigValueFromList<int>("Solver parameters", "filtering_freq", 10);
            n_threads = parser.GetConfigValueFromList<int>("Solver parameters", "n_threads", 1);
            

        }

};


#endif // INPUTS_HPP
