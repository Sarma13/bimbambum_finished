// File: Bubble_Dynamics/src/core/Case_SphericalBoundary.cpp

#include "Case_SphericalBoundary.hpp"
#include "cubic_spline.hpp"
#include "integrands_volume.hpp"
#include <armadillo>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>


using namespace std;
using namespace arma;

static constexpr double pi = 3.14159265358979323846;

Case_SphericalBoundary::Case_SphericalBoundary(const Input &data)
  : BoundaryData(data)
{}

void Case_SphericalBoundary::initialize() {
  // Ns+1 nodes from south to north pole (counterclockwise) on a perfect sphere of radius gamma
  for (int i = 0; i <= Ns; ++i) {
    // Start from south pole (pi) and go to north pole (0)
    double alpha = pi - double(i) * pi / double(Ns);
    r_nodes[i]    = xi * sin(alpha);
    z_nodes[i]    = xi * cos(alpha) - gamma;
    F_nodes[i]    = 0.0;
    curv_nodes[i] = 1.0 / xi;
  }

  double V0;
  V0 = 4.0 / 3.0 * pi * xi * xi * xi;
  V_spherical= V0;
  // enforce r=0 at the poles
  r_nodes.front() = r_nodes.back() = 0.0;
}

void Case_SphericalBoundary::boundary_endpoints_derivatives() {
  // at south pole (now front), tangent is horizontal to +r
  drds1 = +1.0; dzds1 = 0.0;
  // at north pole (now back), tangent is horizontal to -r
  drds2 = -1.0; dzds2 = 0.0;
  dphi1ds1 = dphi1ds2 = dphi2ds1 = dphi2ds2 = 0.0;
}

/*! \brief Remesh the spherical boundary using geometric progression along arc-length.
    Similar to BubbleData::remesh_bubble, but applied to the fluid-fluid interface.
*/
void Case_SphericalBoundary::remesh_boundary() {
  // generate current arc-length distribution for creating later evenly spaced nodes
  vec rs = conv_to<vec>::from(r_nodes);
  vec zs = conv_to<vec>::from(z_nodes);

  vec drs = diff(rs);
  vec dzs = diff(zs);
  vec seg = sqrt(drs % drs + dzs % dzs);
  vec distance_s = join_cols(vec{0.0}, cumsum(seg));
  vector<double> distance = conv_to<vector<double>>::from(distance_s);

  // build clamped-end splines for r(s), z(s), and F(s)
  boundary_endpoints_derivatives();

  cubic_spline spR, spZ, spF;
  spR.set_spline(distance, r_nodes, drds1, drds2);
  spZ.set_spline(distance, z_nodes, dzds1, dzds2);
  spF.set_spline(distance, F_nodes, 0.0, 0.0);

  // geometric progression spacing along total length L
  // Using a ratio that clusters points near the south pole (first node) and sparser near north pole
  int n = Ns + 1;
  double r = 1.001;  // r > 1 gives more points near south pole (index 0)
  double h = (r - 1.0) / (pow(r, n - 1) - 1.0);

  vec e(Ns);

  for (int i = 0; i < Ns; ++i) {
      e(i) = pow(r, 1.0*i);
  }

  vec s_new = h * cumsum(e) * distance[Ns];
  s_new.insert_rows(0, vec(1, fill::zeros));      // s_new size = Ns+1
  s_new(Ns) = distance[Ns];              // ensure last point

  // interpolate to new nodes
  for (int i = 0; i <= Ns; ++i) {
      double s = s_new(i);
      r_nodes[i] = spR.interpolate(s);
      z_nodes[i] = spZ.interpolate(s);
      F_nodes[i] = spF.interpolate(s);
  }
  // enforce axisymmetry
  r_nodes.front() = r_nodes.back() = 0.0;
}

/*! \brief Compute axisymmetric curvature at each boundary node.
  Fits a polynomial in a local 9-point stencil, computes meridional and azimuthal terms.
*/
void Case_SphericalBoundary::boundary_curvature() {
  vector<double> r_copy = r_nodes;
  vector<double> z_copy = z_nodes;
  vector<double> F_copy = F_nodes;

  r_copy.insert(r_copy.begin(), -r_nodes[1]);
  r_copy.insert(r_copy.begin(), -r_nodes[2]);
  r_copy.insert(r_copy.begin(), -r_nodes[3]);
  r_copy.insert(r_copy.begin(), -r_nodes[4]);
  r_copy.push_back(-r_nodes[Ns - 1]);
  r_copy.push_back(-r_nodes[Ns - 2]);
  r_copy.push_back(-r_nodes[Ns - 3]);
  r_copy.push_back(-r_nodes[Ns - 4]);

  z_copy.insert(z_copy.begin(), z_nodes[1]);
  z_copy.insert(z_copy.begin(), z_nodes[2]);
  z_copy.insert(z_copy.begin(), z_nodes[3]);
  z_copy.insert(z_copy.begin(), z_nodes[4]);
  z_copy.push_back(z_nodes[Ns - 1]);
  z_copy.push_back(z_nodes[Ns - 2]);
  z_copy.push_back(z_nodes[Ns - 3]);
  z_copy.push_back(z_nodes[Ns - 4]);

  vec rs = conv_to<vec>::from(r_copy);
  vec zs = conv_to<vec>::from(z_copy);
  vec drs = diff(rs);
  vec dzs = diff(zs);
  vec pointss = sqrt(drs % drs + dzs % dzs);
  vec distance_s = cumsum(pointss);
  distance_s.insert_rows(0, vec(1, fill::zeros));

  int degree = 4;
  int n_elem = 9;

  for (int i = 0; i <= Ns; ++i) {
      vec dist_i = distance_s.subvec(i, i + n_elem - 1);
      vec r_i    = vec(&r_copy[i], n_elem);
      vec z_i    = vec(&z_copy[i], n_elem);
      double s0  = distance_s(i + 4);
      vec r_coeff = polyfit(dist_i, r_i, degree);
      vec z_coeff = polyfit(dist_i, z_i, degree);

      double drds = 4*s0*s0*s0*r_coeff(0) + 3*s0*s0*r_coeff(1) + 2*s0*r_coeff(2) + r_coeff(3);
      double d2r  = 12*s0*s0*r_coeff(0)       + 6*s0*r_coeff(1)       + 2*r_coeff(2);
      double dzds = 4*s0*s0*s0*z_coeff(0) + 3*s0*s0*z_coeff(1) + 2*s0*z_coeff(2) + z_coeff(3);
      double d2z  = 12*s0*s0*z_coeff(0)       + 6*s0*z_coeff(1)       + 2*z_coeff(2);

      double num = drds * d2z - dzds * d2r;
      double den = pow(drds*drds + dzds*dzds, 1.5);
      curv_nodes[i] = num/den + dzds/(r_nodes[i]*sqrt(drds*drds + dzds*dzds));
      if (i==0 || i==Ns) {
          curv_nodes[i] = 2*d2z/(drds*drds*drds);
      }
  }
}

void Case_SphericalBoundary::filter_boundary() {
  vector<double> r_copy(r_nodes.begin(), r_nodes.end());
  vector<double> z_copy(z_nodes.begin(), z_nodes.end());
  vector<double> f_copy(F_nodes.begin(), F_nodes.end());

  r_copy.insert(r_copy.begin(), -r_nodes[1]);
  r_copy.insert(r_copy.begin(), -r_nodes[2]);
  r_copy.push_back(-r_nodes[Ns - 1]);
  r_copy.push_back(-r_nodes[Ns - 2]);

  z_copy.insert(z_copy.begin(), z_nodes[1]);
  z_copy.insert(z_copy.begin(), z_nodes[2]);
  z_copy.push_back(z_nodes[Ns - 1]);
  z_copy.push_back(z_nodes[Ns - 2]);

  f_copy.insert(f_copy.begin(), F_nodes[1]);
  f_copy.insert(f_copy.begin(), F_nodes[2]);
  f_copy.push_back(F_nodes[Ns - 1]);
  f_copy.push_back(F_nodes[Ns - 2]);

  vec rb = conv_to<vec>::from(r_copy);
  vec zb = conv_to<vec>::from(z_copy);
  vec drb = diff(rb), dzb = diff(zb);
  vec seg = sqrt(drb % drb + dzb % dzb);
  vec scum = join_cols(vec{0.0}, cumsum(seg));
  vector<double> dist = conv_to<vector<double>>::from(scum);

  for (int i = 0; i <= Ns; ++i) {
      double s1 = dist[i],   s2 = dist[i+1];
      double s3 = dist[i+2], s4 = dist[i+3];
      double s5 = dist[i+4];
      double c1 = 0.5 * ((s3 - s2) * (s3 - s4)) / ((s1 - s3) * (s1 - s5));
      double c2 = 0.5 * ((s4 - s3)) /       ((s4 - s2));
      double c3 = 0.5 * (1.0 + ((s3 - s2)*(s3 - s4))/((s3 - s1)*(s3 - s5)));
      double c4 = 0.5 * ((s2 - s3)) /       ((s2 - s4));
      double c5 = 0.5 * ((s3 - s2)*(s3 - s4))/((s5 - s1)*(s5 - s3));

      r_nodes[i] = c1*r_copy[i]   + c2*r_copy[i+1] + c3*r_copy[i+2]
                 + c4*r_copy[i+3] + c5*r_copy[i+4];
      z_nodes[i] = c1*z_copy[i]   + c2*z_copy[i+1] + c3*z_copy[i+2]
                 + c4*z_copy[i+3] + c5*z_copy[i+4];
      F_nodes[i] = c1*f_copy[i]   + c2*f_copy[i+1] + c3*f_copy[i+2]
                 + c4*f_copy[i+3] + c5*f_copy[i+4];
  }
  // enforce axisymmetry at poles
  r_nodes.front() = r_nodes.back() = 0.0;
}


/**
 * Computes the volume enclosed by the spherical boundary.
 * 
 * @return The volume of the spherical droplet.
 */
double Case_SphericalBoundary::compute_spherical_volume() {
    double volume = 0.0;

    // Generate surface vectors
    vec rs = conv_to<vec>::from(r_nodes);
    vec zs = conv_to<vec>::from(z_nodes);

    // Compute the arclength
    vec drs = diff(rs);
    vec dzs = diff(zs);
    vec points = sqrt(drs % drs + dzs % dzs);
    vec distance_s = cumsum(points);
    vec begin(1, fill::zeros);
    distance_s = join_cols(begin, distance_s);
    vector<double> distance = conv_to<vector<double>>::from(distance_s);

    // Create cubic splines for the boundary
    cubic_spline spline_rs;
    spline_rs.set_spline(distance, r_nodes, drds1, drds2);

    cubic_spline spline_zs;
    spline_zs.set_spline(distance, z_nodes, dzds1, dzds2);

    // Perform the integration
    double a, b;
    double m_ar, m_br, m_cr, m_dr, m_az, m_bz, m_cz, m_dz;

    for (int i = 0; i < Ns; ++i) {
        m_ar = spline_rs.m_a[i];
        m_br = spline_rs.m_b[i];
        m_cr = spline_rs.m_c[i];
        m_dr = rs(i);

        m_az = spline_zs.m_a[i];
        m_bz = spline_zs.m_b[i];
        m_cz = spline_zs.m_c[i];
        m_dz = zs(i);

        gsl_integration_workspace *w_volume = gsl_integration_workspace_alloc(2000);
        double result, error;

        volume_params params;
        params.coeff_r0 = m_ar;
        params.coeff_r1 = m_br;
        params.coeff_r2 = m_cr;
        params.coeff_r3 = m_dr;
        params.coeff_z0 = m_az;
        params.coeff_z1 = m_bz;
        params.coeff_z2 = m_cz;
        params.coeff_z3 = m_dz;

        gsl_function vol;
        vol.function = &volume_spherical_droplet;  // Use the spherical droplet volume function
        vol.params = &params;

        a = 0.0;
        b = distance_s[1 + i] - distance_s[i];

        gsl_integration_qags(&vol, a, b, 0.0, 1e-10, 2000, w_volume, &result, &error);

        volume = volume + result;
        gsl_integration_workspace_free(w_volume);
    }
    
    return volume;
}





