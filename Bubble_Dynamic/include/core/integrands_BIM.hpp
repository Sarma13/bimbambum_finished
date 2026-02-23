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
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file integrands_BIM.hpp
    \brief  Integrands of the boundary integrals equations.

    The discretized integrals consider a linear evolution of the potential and normal velocities
    along the segments separating two adjacent nodes. The positions of the bubble and fluid-fluid interface (r and z)
    are extracted from the cubic splines. The function G1 and G2 contains the discretized integral of the Green's
    function and the functions H1 and H2 the discretized integral of the normal derivative of the Green's function.
*/

#ifndef INTEGRANDS_BIM_HPP
#define INTEGRANDS_BIM_HPP

struct G_H_params {
    double coeff_r0;
    double coeff_r1;
    double coeff_r2;
    double coeff_r3;
    double coeff_z0;
    double coeff_z1;
    double coeff_z2;
    double coeff_z3;
    double r0;
    double z0;
    double a;
    double b;
    int iter;
    int N_b;
    int N_s;
};

/*! \brief  Discretized integral of the Green function (contribution of node i).*/
double G1(double xi, void *p) {
    struct G_H_params *params = (struct G_H_params *) p;

    // Coefficients for the cubic splines
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    double coeff_z0 = (params->coeff_z0);
    double coeff_z1 = (params->coeff_z1);
    double coeff_z2 = (params->coeff_z2);
    double coeff_z3 = (params->coeff_z3);

    // Position of the nodes
    double r0 = (params->r0);
    double z0 = (params->z0);

    // Integration intervals
    double a = (params->a);
    double b = (params->b);


    int Nb = (params->N_b); // Bubble divided in Nb elements
    int Ns = (params->N_s); // Fluid-fluid interface divided in Nb elements
    int j = (params->iter); // jth element

    double M = 0.0; // pre-factor for linear interpolation between adjacent nodes

    // (r, z) coordinates and their spatial derivatives
    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    double z = coeff_z3 + coeff_z2 * xi + coeff_z1 * xi * xi + coeff_z0 * xi * xi * xi;
    double drdxi = coeff_r2 + 2.0 * coeff_r1 * xi + 3.0 * coeff_r0 * xi * xi;
    double dzdxi = coeff_z2 + 2.0 * coeff_z1 * xi + 3.0 * coeff_z0 * xi * xi;


    // Special treatment for M for node points on the symmetry axis.
    // See Curtiss (PhD thesis, University of Birmingham, 2009)
    if (j == 0) {
        M = 1.0 - xi * xi / (b * b);
    } else if (j == Nb - 1 or j == Ns - 1) {
        M = (b - xi) * (b - xi) / ((b - a) * (b - a));
    } else {
        M = (b - xi) / (b - a);
    }

    double s = (r + r0) * (r + r0) + (z - z0) * (z - z0);
    double k = 1.0 - (4.0 * r * r0 / s);
    if (std::abs(k) < 2.22e-16) {
        k = 2.22e-16;
    }

    // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
    double PK[13] = {1.38629436111989061883446424, 0.09657359027997265470861606,
                     0.03088514453248461827359656, 0.01493760036978098687568492,
                     0.00876631219862835129486730, 0.00575489991651211831713086,
                     0.00406819648916235957842217, 0.00316713448114840176286619,
                     0.00385918735043451810914414, 0.00697248927202287553710545,
                     0.00700030498423661873526199, 0.00235535576237663133325157,
                     0.00016175003824586587091022};
    double QK[13] = {0.50000000000000000000000000, 0.12500000000000000000000000,
                     0.07031250000000000000000000, 0.04882812499999999987824278,
                     0.03738403320299965249042380, 0.03028106526770420433989236,
                     0.02544378896278751497219371, 0.02189639358590439516170295,
                     0.01859695172048566289195740, 0.01326644642298080552433290,
                     0.00572150665129845121056799, 0.00098749488654029748460148,
                     0.00003519107157048046293917};

    double P_val =
            PK[0] + k * PK[1] + k * k * PK[2] + k * k * k * PK[3] +
            k * k * k * k * PK[4] + k * k * k * k * k * PK[5] +
            k * k * k * k * k * k * PK[6] +
            k * k * k * k * k * k * k * PK[7] +
            k * k * k * k * k * k * k * k * PK[8] +
            k * k * k * k * k * k * k * k * k * PK[9] +
            k * k * k * k * k * k * k * k * k * k * PK[10] +
            k * k * k * k * k * k * k * k * k * k * k * PK[11] +
            k * k * k * k * k * k * k * k * k * k * k * k * PK[12];
    double Q_val =
            QK[0] + k * QK[1] + k * k * QK[2] + k * k * k * QK[3] +
            k * k * k * k * QK[4] + k * k * k * k * k * QK[5] +
            k * k * k * k * k * k * QK[6] +
            k * k * k * k * k * k * k * QK[7] +
            k * k * k * k * k * k * k * k * QK[8] +
            k * k * k * k * k * k * k * k * k * QK[9] +
            k * k * k * k * k * k * k * k * k * k * QK[10] +
            k * k * k * k * k * k * k * k * k * k * k * QK[11] +
            k * k * k * k * k * k * k * k * k * k * k * k * QK[12];

    double pi = 3.14159265358979323846264338328;
    double K_val = P_val - Q_val * std::log(k);
    if (r0 == 0.0) {
        K_val = pi / 2.0;
    }
    double G1 = 4.0 * r * K_val * std::sqrt(drdxi * drdxi + dzdxi * dzdxi) / std::sqrt(s) * M;

    return G1;
}

/*! \brief  Discretized integral of the Green function (contribution of node i + 1).*/
double G2(double xi, void *p) {
    struct G_H_params *params = (struct G_H_params *) p;

    // Coefficients for the cubic splines
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    double coeff_z0 = (params->coeff_z0);
    double coeff_z1 = (params->coeff_z1);
    double coeff_z2 = (params->coeff_z2);
    double coeff_z3 = (params->coeff_z3);

    // Position of the nodes
    double r0 = (params->r0);
    double z0 = (params->z0);

    // Integration intervals
    double a = (params->a);
    double b = (params->b);


    int Nb = (params->N_b); // Bubble divided in Nb elements
    int Ns = (params->N_s); // Fluid-fluid interface divided in Nb elements
    int j = (params->iter); // jth element

    double M = 0.0; // pre-factor for linear interpolation between adjacent nodes

    // (r, z) coordinates and their spatial derivatives
    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    double z = coeff_z3 + coeff_z2 * xi + coeff_z1 * xi * xi + coeff_z0 * xi * xi * xi;
    double drdxi = coeff_r2 + 2.0 * coeff_r1 * xi + 3.0 * coeff_r0 * xi * xi;
    double dzdxi = coeff_z2 + 2.0 * coeff_z1 * xi + 3.0 * coeff_z0 * xi * xi;

    // Special treatment for M for node points on the symmetry axis.
    // See Curtiss (PhD thesis, University of Birmingham, 2009)
    if (j == 1) {
        M = xi * xi / (b * b);
    } else if (j == Nb or j == Ns) {
        M = 1.0 - (b - xi) * (b - xi) / ((b - a) * (b - a));
    } else {
        M = (xi - a) / (b - a);
    }

    double s = (r + r0) * (r + r0) + (z - z0) * (z - z0);
    double k = 1.0 - (4.0 * r * r0 / s);
    if (std::abs(k) < 2.22e-16) {
        k = 2.22e-16;
    }

    // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
    double PK[13] = {1.38629436111989061883446424, 0.09657359027997265470861606,
                     0.03088514453248461827359656, 0.01493760036978098687568492,
                     0.00876631219862835129486730, 0.00575489991651211831713086,
                     0.00406819648916235957842217, 0.00316713448114840176286619,
                     0.00385918735043451810914414, 0.00697248927202287553710545,
                     0.00700030498423661873526199, 0.00235535576237663133325157,
                     0.00016175003824586587091022};
    double QK[13] = {0.50000000000000000000000000, 0.12500000000000000000000000,
                     0.07031250000000000000000000, 0.04882812499999999987824278,
                     0.03738403320299965249042380, 0.03028106526770420433989236,
                     0.02544378896278751497219371, 0.02189639358590439516170295,
                     0.01859695172048566289195740, 0.01326644642298080552433290,
                     0.00572150665129845121056799, 0.00098749488654029748460148,
                     0.00003519107157048046293917};

    double P_val =
            PK[0] + k * PK[1] + k * k * PK[2] + k * k * k * PK[3] +
            k * k * k * k * PK[4] + k * k * k * k * k * PK[5] +
            k * k * k * k * k * k * PK[6] +
            k * k * k * k * k * k * k * PK[7] +
            k * k * k * k * k * k * k * k * PK[8] +
            k * k * k * k * k * k * k * k * k * PK[9] +
            k * k * k * k * k * k * k * k * k * k * PK[10] +
            k * k * k * k * k * k * k * k * k * k * k * PK[11] +
            k * k * k * k * k * k * k * k * k * k * k * k * PK[12];
    double Q_val =
            QK[0] + k * QK[1] + k * k * QK[2] + k * k * k * QK[3] +
            k * k * k * k * QK[4] + k * k * k * k * k * QK[5] +
            k * k * k * k * k * k * QK[6] +
            k * k * k * k * k * k * k * QK[7] +
            k * k * k * k * k * k * k * k * QK[8] +
            k * k * k * k * k * k * k * k * k * QK[9] +
            k * k * k * k * k * k * k * k * k * k * QK[10] +
            k * k * k * k * k * k * k * k * k * k * k * QK[11] +
            k * k * k * k * k * k * k * k * k * k * k * k * QK[12];

    double pi = 3.14159265358979323846264338328;
    double K_val = P_val - Q_val * std::log(k);
    if (r0 == 0.0) {
        K_val = pi / 2.0;
    }
    double G2 = 4.0 * r * K_val * std::sqrt(drdxi * drdxi + dzdxi * dzdxi) / std::sqrt(s) * M;

    return G2;
}

/*! \brief  Discretized integral of the normal derivative Green function (contribution of node i).*/
double H1(double xi, void *p) {
    struct G_H_params *params = (struct G_H_params *) p;

    // Coefficients for the cubic splines
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    double coeff_z0 = (params->coeff_z0);
    double coeff_z1 = (params->coeff_z1);
    double coeff_z2 = (params->coeff_z2);
    double coeff_z3 = (params->coeff_z3);

    // Position of the nodes
    double r0 = (params->r0);
    double z0 = (params->z0);

    // Integration intervals
    double a = (params->a);
    double b = (params->b);


    int Nb = (params->N_b); // Bubble divided in Nb elements
    int Ns = (params->N_s); // Fluid-fluid interface divided in Nb elements
    int j = (params->iter); // jth element

    double M = 0.0; // pre-factor for linear interpolation between adjacent nodes

    // (r, z) coordinates and their spatial derivatives
    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    double z = coeff_z3 + coeff_z2 * xi + coeff_z1 * xi * xi + coeff_z0 * xi * xi * xi;
    double drdxi = coeff_r2 + 2.0 * coeff_r1 * xi + 3.0 * coeff_r0 * xi * xi;
    double dzdxi = coeff_z2 + 2.0 * coeff_z1 * xi + 3.0 * coeff_z0 * xi * xi;

    // Special treatment for M for node points on the symmetry axis.
    // See Curtiss (PhD thesis, University of Birmingham, 2009)
    if (j == 0) {
        M = 1.0 - xi * xi / (b * b);
    } else if (j == Nb - 1 or j == Ns - 1) {
        M = (b - xi) * (b - xi) / ((b - a) * (b - a));
    } else {
        M = (b - xi) / (b - a);
    }

    double pi = 3.14159265358979323846264338328;

    double H1 = 0.0;
    if (r0 == 0.0) {
        double s = r * r + (z - z0) * (z - z0);
        double E_fac = dzdxi * r - drdxi * (z - z0);
        H1 = -2.0 * pi * r * E_fac / std::pow(s, 1.5) * M;
    } else {
        double s = (r + r0) * (r + r0) + (z - z0) * (z - z0);
        double k = 1.0 - (4.0 * r * r0 / s);
        double k1 = std::sqrt(4.0 * r * r0 / s);

        if (std::abs(k) < 2.22e-16) {
            k = 2.22e-16;
        }

        // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
        double PK[13] = {1.38629436111989061883446424, 0.09657359027997265470861606,
                         0.03088514453248461827359656, 0.01493760036978098687568492,
                         0.00876631219862835129486730, 0.00575489991651211831713086,
                         0.00406819648916235957842217, 0.00316713448114840176286619,
                         0.00385918735043451810914414, 0.00697248927202287553710545,
                         0.00700030498423661873526199, 0.00235535576237663133325157,
                         0.00016175003824586587091022};
        double QK[13] = {0.50000000000000000000000000, 0.12500000000000000000000000,
                         0.07031250000000000000000000, 0.04882812499999999987824278,
                         0.03738403320299965249042380, 0.03028106526770420433989236,
                         0.02544378896278751497219371, 0.02189639358590439516170295,
                         0.01859695172048566289195740, 0.01326644642298080552433290,
                         0.00572150665129845121056799, 0.00098749488654029748460148,
                         0.00003519107157048046293917};

        double P_val =
                PK[0] + k * PK[1] + k * k * PK[2] + k * k * k * PK[3] +
                k * k * k * k * PK[4] + k * k * k * k * k * PK[5] +
                k * k * k * k * k * k * PK[6] +
                k * k * k * k * k * k * k * PK[7] +
                k * k * k * k * k * k * k * k * PK[8] +
                k * k * k * k * k * k * k * k * k * PK[9] +
                k * k * k * k * k * k * k * k * k * k * PK[10] +
                k * k * k * k * k * k * k * k * k * k * k * PK[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * PK[12];
        double Q_val =
                QK[0] + k * QK[1] + k * k * QK[2] + k * k * k * QK[3] +
                k * k * k * k * QK[4] + k * k * k * k * k * QK[5] +
                k * k * k * k * k * k * QK[6] +
                k * k * k * k * k * k * k * QK[7] +
                k * k * k * k * k * k * k * k * QK[8] +
                k * k * k * k * k * k * k * k * k * QK[9] +
                k * k * k * k * k * k * k * k * k * k * QK[10] +
                k * k * k * k * k * k * k * k * k * k * k * QK[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * QK[12];

        // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
        double PE[13] = {1.00000000000000000000000000, 0.44314718055994530941723212,
                         0.05680519270997949103146207, 0.02183137044373718396138156,
                         0.01154452141883701103542361, 0.00714201318820502987066619,
                         0.00485846659881274463594893, 0.00366680346394393045387665,
                         0.00426469424891906813517382, 0.00745727014212456596918847,
                         0.00741871341163044927753980, 0.00248933547336496339904368,
                         0.00017076513539687204438478};
        double QE[13] = {0.00000000000000000000000000, 0.25000000000000000000000000,
                         0.09375000000000000000000000, 0.05859374999999999987183993,
                         0.04272460937486806132127659, 0.03364562817049392175150879,
                         0.02775688834606027631579899, 0.02358179637126350018588892,
                         0.01984699815591322481619125, 0.01407783193717112862114295,
                         0.00605330008329266855825149, 0.00104321202098794509265719,
                         0.00003714782778910401536553};

        double R_val =
                PE[0] + k * PE[1] + k * k * PE[2] + k * k * k * PE[3] +
                k * k * k * k * PE[4] + k * k * k * k * k * PE[5] +
                k * k * k * k * k * k * PE[6] +
                k * k * k * k * k * k * k * PE[7] +
                k * k * k * k * k * k * k * k * PE[8] +
                k * k * k * k * k * k * k * k * k * PE[9] +
                k * k * k * k * k * k * k * k * k * k * PE[10] +
                k * k * k * k * k * k * k * k * k * k * k * PE[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * PE[12];
        double S_val =
                QE[0] + k * QE[1] + k * k * QE[2] + k * k * k * QE[3] +
                k * k * k * k * QE[4] + k * k * k * k * k * QE[5] +
                k * k * k * k * k * k * QE[6] +
                k * k * k * k * k * k * k * QE[7] +
                k * k * k * k * k * k * k * k * QE[8] +
                k * k * k * k * k * k * k * k * k * QE[9] +
                k * k * k * k * k * k * k * k * k * k * QE[10] +
                k * k * k * k * k * k * k * k * k * k * k * QE[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * QE[12];


        double K_val = P_val - Q_val * std::log(k);
        double E_val = R_val - S_val * std::log(k);


        double K_fac = 2.0 * dzdxi * r0 / (k1 * k1);
        double E_fac = (dzdxi * (r + r0) - drdxi * (z - z0) - 2.0 * dzdxi * r0 / (k1 * k1)) / (1.0 - k1 * k1);
        H1 = -4.0 * (r * (E_fac * E_val + K_fac * K_val) / std::pow(s, 1.5)) * M;
    }

    return H1;
}

/*! \brief  Discretized integral of the normal derivative Green function (contribution of node i + 1).*/
double H2(double xi, void *p) {
    struct G_H_params *params = (struct G_H_params *) p;

    // Coefficients for the cubic splines
    double coeff_r0 = (params->coeff_r0);
    double coeff_r1 = (params->coeff_r1);
    double coeff_r2 = (params->coeff_r2);
    double coeff_r3 = (params->coeff_r3);
    double coeff_z0 = (params->coeff_z0);
    double coeff_z1 = (params->coeff_z1);
    double coeff_z2 = (params->coeff_z2);
    double coeff_z3 = (params->coeff_z3);

    // Position of the nodes
    double r0 = (params->r0);
    double z0 = (params->z0);

    // Integration intervals
    double a = (params->a);
    double b = (params->b);


    int Nb = (params->N_b); // Bubble divided in Nb elements
    int Ns = (params->N_s); // Fluid-fluid interface divided in Nb elements
    int j = (params->iter); // jth element

    double M = 0.0; // pre-factor for linear interpolation between adjacent nodes

    // (r, z) coordinates and their spatial derivatives
    double r = coeff_r3 + coeff_r2 * xi + coeff_r1 * xi * xi + coeff_r0 * xi * xi * xi;
    double z = coeff_z3 + coeff_z2 * xi + coeff_z1 * xi * xi + coeff_z0 * xi * xi * xi;
    double drdxi = coeff_r2 + 2.0 * coeff_r1 * xi + 3.0 * coeff_r0 * xi * xi;
    double dzdxi = coeff_z2 + 2.0 * coeff_z1 * xi + 3.0 * coeff_z0 * xi * xi;

    // Special treatment for M for node points on the symmetry axis.
    // See Curtiss (PhD thesis, University of Birmingham, 2009)
    if (j == 1) {
        M = xi * xi / (b * b);
    } else if (j == Nb or j == Ns) {
        M = 1.0 - (b - xi) * (b - xi) / ((b - a) * (b - a));
    } else {
        M = (xi - a) / (b - a);
    }

    double pi = 3.14159265358979323846264338328;

    double H2 = 0.0;
    if (r0 == 0.0) {
        double s = r * r + (z - z0) * (z - z0);
        double E_fac = dzdxi * r - drdxi * (z - z0);
        H2 = -2.0 * pi * r * E_fac / std::pow(s, 1.5) * M;
    } else {
        double s = (r + r0) * (r + r0) + (z - z0) * (z - z0);
        double k = 1.0 - (4.0 * r * r0 / s);
        double k1 = std::sqrt(4.0 * r * r0 / s);

        if (std::abs(k) < 2.22e-16) {
            k = 2.22e-16;
        }

        // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
        double PK[13] = {1.38629436111989061883446424, 0.09657359027997265470861606,
                         0.03088514453248461827359656, 0.01493760036978098687568492,
                         0.00876631219862835129486730, 0.00575489991651211831713086,
                         0.00406819648916235957842217, 0.00316713448114840176286619,
                         0.00385918735043451810914414, 0.00697248927202287553710545,
                         0.00700030498423661873526199, 0.00235535576237663133325157,
                         0.00016175003824586587091022};
        double QK[13] = {0.50000000000000000000000000, 0.12500000000000000000000000,
                         0.07031250000000000000000000, 0.04882812499999999987824278,
                         0.03738403320299965249042380, 0.03028106526770420433989236,
                         0.02544378896278751497219371, 0.02189639358590439516170295,
                         0.01859695172048566289195740, 0.01326644642298080552433290,
                         0.00572150665129845121056799, 0.00098749488654029748460148,
                         0.00003519107157048046293917};

        double P_val =
                PK[0] + k * PK[1] + k * k * PK[2] + k * k * k * PK[3] +
                k * k * k * k * PK[4] + k * k * k * k * k * PK[5] +
                k * k * k * k * k * k * PK[6] +
                k * k * k * k * k * k * k * PK[7] +
                k * k * k * k * k * k * k * k * PK[8] +
                k * k * k * k * k * k * k * k * k * PK[9] +
                k * k * k * k * k * k * k * k * k * k * PK[10] +
                k * k * k * k * k * k * k * k * k * k * k * PK[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * PK[12];
        double Q_val =
                QK[0] + k * QK[1] + k * k * QK[2] + k * k * k * QK[3] +
                k * k * k * k * QK[4] + k * k * k * k * k * QK[5] +
                k * k * k * k * k * k * QK[6] +
                k * k * k * k * k * k * k * QK[7] +
                k * k * k * k * k * k * k * k * QK[8] +
                k * k * k * k * k * k * k * k * k * QK[9] +
                k * k * k * k * k * k * k * k * k * k * QK[10] +
                k * k * k * k * k * k * k * k * k * k * k * QK[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * QK[12];

        // Tabulated polynomials for complete elliptic integrals approximation (1st and 2n kind)
        double PE[13] = {1.00000000000000000000000000, 0.44314718055994530941723212,
                         0.05680519270997949103146207, 0.02183137044373718396138156,
                         0.01154452141883701103542361, 0.00714201318820502987066619,
                         0.00485846659881274463594893, 0.00366680346394393045387665,
                         0.00426469424891906813517382, 0.00745727014212456596918847,
                         0.00741871341163044927753980, 0.00248933547336496339904368,
                         0.00017076513539687204438478};
        double QE[13] = {0.00000000000000000000000000, 0.25000000000000000000000000,
                         0.09375000000000000000000000, 0.05859374999999999987183993,
                         0.04272460937486806132127659, 0.03364562817049392175150879,
                         0.02775688834606027631579899, 0.02358179637126350018588892,
                         0.01984699815591322481619125, 0.01407783193717112862114295,
                         0.00605330008329266855825149, 0.00104321202098794509265719,
                         0.00003714782778910401536553};

        double R_val =
                PE[0] + k * PE[1] + k * k * PE[2] + k * k * k * PE[3] +
                k * k * k * k * PE[4] + k * k * k * k * k * PE[5] +
                k * k * k * k * k * k * PE[6] +
                k * k * k * k * k * k * k * PE[7] +
                k * k * k * k * k * k * k * k * PE[8] +
                k * k * k * k * k * k * k * k * k * PE[9] +
                k * k * k * k * k * k * k * k * k * k * PE[10] +
                k * k * k * k * k * k * k * k * k * k * k * PE[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * PE[12];
        double S_val =
                QE[0] + k * QE[1] + k * k * QE[2] + k * k * k * QE[3] +
                k * k * k * k * QE[4] + k * k * k * k * k * QE[5] +
                k * k * k * k * k * k * QE[6] +
                k * k * k * k * k * k * k * QE[7] +
                k * k * k * k * k * k * k * k * QE[8] +
                k * k * k * k * k * k * k * k * k * QE[9] +
                k * k * k * k * k * k * k * k * k * k * QE[10] +
                k * k * k * k * k * k * k * k * k * k * k * QE[11] +
                k * k * k * k * k * k * k * k * k * k * k * k * QE[12];


        double K_val = P_val - Q_val * std::log(k);
        double E_val = R_val - S_val * std::log(k);


        double K_fac = 2.0 * dzdxi * r0 / (k1 * k1);
        double E_fac = (dzdxi * (r + r0) - drdxi * (z - z0) - 2.0 * dzdxi * r0 / (k1 * k1)) / (1.0 - k1 * k1);
        H2 = -4.0 * (r * (E_fac * E_val + K_fac * K_val) / std::pow(s, 1.5)) * M;
    }

    return H2;
}

#endif // INTEGRANGS_BIM_HPP
