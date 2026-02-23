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

/*! \file main.cpp
    \brief main function for solver execution, called at program startup.
*/

#include <iostream>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "time_stepper.hpp"
#include "Inputs.hpp"
#include "BubbleData.hpp"
#include "BoundaryData.hpp"
#include "BIM_solver.hpp"
#include "Case_SphericalBoundary.hpp"

#include "Case_FlatBoundary.hpp"
#include "Case_RayleighBubble.hpp"
#include "Case_RayleighPlessetBubble.hpp"

int main(int argc, const char **argv) {

    // ---------- User-defined inputs processing ----------
    if (argc != 2) {
        std::cerr << "Too many or too few input arguments: please only provide the JSON configuration file. "
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    const std::string input_conf = argv[1];
    Input data; //reads inputs
    data.initialize(input_conf, argc, argv); // apply inputs to simulation

    // ---------- Parallel processing ----------
    #ifdef _OPENMP
        const int max_threads = omp_get_max_threads();
        int usr_defined_threads = data.n_threads;

        if (usr_defined_threads <= max_threads){
            omp_set_num_threads(usr_defined_threads);
        }
        else{
            std::cerr << "Selected amount of threads 'n_threads' larger than available amount." << std::endl;
            std::cerr << "Maximum number of threads that can be used: " << max_threads << std::endl;
            exit(EXIT_FAILURE);
        }
    #endif

    // ---------- Simulation initialization ----------
    BIM_solver step(data.Nb, data.Ns, data.dumper_filename); //initialize BIM solver

    std::unique_ptr <BubbleData> bubble = nullptr;// bubble surface properties
    std::unique_ptr <BoundaryData> boundary = nullptr; // fluid-fluid interface properties

    // Bubble initial conditions
    if (data.bubble_dynamics == "Rayleigh_Bubble") { // initial bubble properties taken from Rayleigh model
        bubble = std::make_unique<Case_RayleighBubble>(data);
        bubble->initialize();
    } else if (data.bubble_dynamics ==
               "Rayleigh_Plesset_Bubble") { // initial bubble properties taken from Rayleigh-Plesset model
        bubble = std::make_unique<Case_RayleighPlessetBubble>(data);
        bubble->initialize();
    } else {
        std::cerr << "Simulation case unknown" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Fluid-fluid interface initial conditions
    if (data.boundary == "from_code") { // initially flat surface at z = 0 (surface elasticity can be included)
        boundary = std::make_unique<Case_FlatBoundary>(data);
        boundary->initialize();
    } 
    
    else if(data.boundary == "spherical") { // your new droplet/spherical interface
        boundary = std::make_unique<Case_SphericalBoundary>(data);
        boundary->initialize();
    } 
    
    else if (data.boundary == "from_file") { // initial conditions taken from file
        std::cerr << "Boundary from file not implemented yet" << std::endl;
        exit(EXIT_FAILURE);

    }
    
    else {
        std::cerr << "Simulation case unknown" << std::endl;
        exit(EXIT_FAILURE);
    }

    // ---------- COMPUTE LOOP ----------
    // loops until intersection condition reached (usually jet crossing opposite interface)
    bool compute = true;
    int intersection = 0;
    step.time_step = 0;
    step.time = bubble->t0;

    while (compute) {

        step.write_solution(bubble, boundary, data); // write to file simulation status

        // ---------- TIME INTEGRATION ----------
        time_integration(bubble, boundary, data, step); //perform one time step

        // ---------- SURFACES FILTERING ----------
        if ((step.time_step + 1) % data.filtering_freq ==
            0.0) { // bubble and boundary surfaces filtered every N time steps
            bubble->filter_bubble();
            boundary->filter_boundary(); // TODO: implement boundary filtering for Case_BoundaryFromFile
        }

        // ---------- SURFACES REMESHING ----------
        bubble->remesh_bubble(); // bubble and surface re-gridded at every time step
        boundary->remesh_boundary(); // TODO: implement boundary remeshing for Case_BoundaryFromFile

        // ---------- BUBBLE INTERSECTION CHECK ----------
        intersection = bubble->intersect(); // check for time loop stop condition

        if (intersection == 1) {
            compute = false;
        }

        //imposed values at boundaries edges for stability
        bubble->r_nodes[0] = 0.0;
        bubble->r_nodes[data.Nb] = 0.0;
        if (data.boundary == "from_code") {
            std::cout<< "BOUNDARY FROM CODE" << std::endl;
            std::cout<< data.boundary << std::endl;
            boundary->r_nodes[data.Ns] = 0.0;
        }

        if (data.boundary == "spherical") {

            
            std::cout<< data.boundary << std::endl;
            boundary->r_nodes[data.Ns] = 0.0;
            boundary->r_nodes[0] = 0.0;

            

    
        }


        
        // ---------- UPDATE SOLUTION ----------
        step.time_step = step.time_step + 1; // time step count
        step.time = step.time + step.dt; // simulation time
        std::cout << "Time step: " << step.time_step << "\t" << "Simulation time: " << step.time << "\t"
                  << "Time step dt: " << step.dt << std::endl;

    }
    step.nodes_position.close(); // closes output file at the end of the simulation

    return EXIT_SUCCESS;
}
