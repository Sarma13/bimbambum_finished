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

/*! \file integrands_volume.hpp
    \brief  Integrands for computing the bubble volume.
*/

#ifndef INTEGRANDS_VOLUME_HPP
#define INTEGRANDS_VOLUME_HPP

struct volume_params {
    double coeff_r0;
    double coeff_r1;
    double coeff_r2;
    double coeff_r3;
    double coeff_z0;
    double coeff_z1;
    double coeff_z2;
    double coeff_z3;
};


inline double volume_b(double xi, void *p) {

    double pi = 3.14159265358979323846264338328;
    struct volume_params *params = (struct volume_params *) p;
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    double coeff_z0 = (params->coeff_z0);
    double coeff_z1 = (params->coeff_z1);
    double coeff_z2 = (params->coeff_z2);
    double coeff_z3 = (params->coeff_z3);

    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    double dzdxi = coeff_z2 + 2.0 * coeff_z1 * xi + 3.0 * coeff_z0 * xi * xi;

    double volume_b = -r * r * dzdxi * pi;

    return volume_b;
}

/**
 * Calculates the volume integrand for a spherical droplet.
 * 
 * @param xi The parametric coordinate.
 * @param p Pointer to the volume_params structure.
 * @return The volume integrand at the specified parametric coordinate.
 */
inline double volume_spherical_droplet(double xi, void *p) {
    double pi = 3.14159265358979323846264338328;
    struct volume_params *params = (struct volume_params *) p;
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    
    // For a spherical droplet, we compute the radius as a function of xi
    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    
    // For a perfect sphere, you would typically use r^3, but this allows for deformations
    // represented by the polynomial coefficients
    double volume_element = r * r * r * pi * (4.0/3.0);
    
    return volume_element;
}

#endif // INTEGRANDS_VOLUME_HPP