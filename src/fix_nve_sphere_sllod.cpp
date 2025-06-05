// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nve_sphere_sllod.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "modify.h"
#include "math_extra.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

FixNVESphereSllod::FixNVESphereSllod(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphere(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere/sllod") == 0 && narg < 3)
    utils::missing_cmd_args(FLERR, "fix nve/sphere/sllod", error);
}

/* ---------------------------------------------------------------------- */

void FixNVESphereSllod::init()
{
  FixNVESphere::init();

  // check fix deform remap settings
  auto deform = modify->get_fix_by_style("^deform");
  if (deform.size() < 1)
    error->all(FLERR,"Using fix nve/sphere/sllod with no fix deform defined");

  for (auto &ifix : deform) {
    auto f = dynamic_cast<FixDeform *>(ifix);
    if (f && (f->remapflag != Domain::V_REMAP))
      error->all(FLERR,"Using fix nve/sphere/sllod with inconsistent fix deform remap option");
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereSllod::initial_integrate(int /*vflag*/)
{
  double dtfm,dtirotate,msq,scale,s2,inv_len_mu;
  double g[3], w[3], w_temp[3], a[3];
  double Q[3][3], Q_temp[3][3], R[3][3];
  double h_two[6],vdelu[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfrotate = dtf / inertia;

  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      v[i][0] -= dtf*vdelu[0];
      v[i][1] -= dtf*vdelu[1];
      v[i][2] -= dtf*vdelu[2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  if (extra == 1) {
    double **mu = atom->mu;
    if (dlm == 0) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          if (mu[i][3] > 0.0) {
            g[0] = mu[i][0] + dtv*(omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
            g[1] = mu[i][1] + dtv*(omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
            g[2] = mu[i][2] + dtv*(omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
            msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
            scale = mu[i][3]/sqrt(msq);
            mu[i][0] = g[0]*scale;
            mu[i][1] = g[1]*scale;
            mu[i][2] = g[2]*scale;
          }
    } else {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit && mu[i][3] > 0.0) {
          inv_len_mu = 1.0/mu[i][3];
          a[0] = mu[i][0]*inv_len_mu;
          a[1] = mu[i][1]*inv_len_mu;
          a[2] = mu[i][2]*inv_len_mu;
          s2 = a[0]*a[0] + a[1]*a[1];
          if (s2 != 0.0) {
            scale = (1.0 - a[2])/s2;
            Q[0][0] = 1.0 - scale*a[0]*a[0];
            Q[0][1] = -scale*a[0]*a[1];
            Q[0][2] = -a[0];
            Q[1][0] = -scale*a[0]*a[1];
            Q[1][1] = 1.0 - scale*a[1]*a[1];
            Q[1][2] = -a[1];
            Q[2][0] = a[0];
            Q[2][1] = a[1];
            Q[2][2] = 1.0 - scale*(a[0]*a[0] + a[1]*a[1]);
          } else {
            Q[0][0] = 1.0/a[2];  Q[0][1] = 0.0;       Q[0][2] = 0.0;
            Q[1][0] = 0.0;       Q[1][1] = 1.0/a[2];  Q[1][2] = 0.0;
            Q[2][0] = 0.0;       Q[2][1] = 0.0;       Q[2][2] = 1.0/a[2];
          }
          w[0] = omega[i][0]; w[1] = omega[i][1]; w[2] = omega[i][2];
          matvec(Q,w,w_temp);
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);
          matvec(R,w_temp,w);
          transpose_times3(R,Q,Q_temp);
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);
          matvec(R,w,w_temp);
          transpose_times3(R,Q_temp,Q);
          BuildRzMatrix(R, 2.0*dtf/force->ftm2v*w_temp[2]);
          matvec(R,w_temp,w);
          transpose_times3(R,Q,Q_temp);
          BuildRyMatrix(R, dtf/force->ftm2v*w[1]);
          matvec(R,w,w_temp);
          transpose_times3(R,Q_temp,Q);
          BuildRxMatrix(R, dtf/force->ftm2v*w_temp[0]);
          matvec(R,w_temp,w);
          transpose_times3(R,Q,Q_temp);
          transpose_matvec(Q_temp,w,w_temp);
          omega[i][0] = w_temp[0];
          omega[i][1] = w_temp[1];
          omega[i][2] = w_temp[2];
          mu[i][0] = Q_temp[2][0] * mu[i][3];
          mu[i][1] = Q_temp[2][1] * mu[i][3];
          mu[i][2] = Q_temp[2][2] * mu[i][3];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereSllod::final_integrate()
{
  double dtfm,dtirotate;
  double h_two[6],vdelu[3];

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfrotate = dtf / inertia;

  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      v[i][0] -= dtf*vdelu[0];
      v[i][1] -= dtf*vdelu[1];
      v[i][2] -= dtf*vdelu[2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}

