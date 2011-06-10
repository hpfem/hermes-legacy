#include "../hermes2d.h"

namespace WeakFormsMaxwell {
  DefaultJacobianMagnetostatics::DefaultJacobianMagnetostatics(int i, int j, std::string area,
							       HermesFunction* coeff,
                                                               SymFlag sym, GeomType gt,
                                                               int order_increase)
    : WeakForm::MatrixFormVol(i, j, area, sym), idx_j(j), coeff(coeff), gt(gt),
    order_increase(order_increase) 
  { 
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }
  DefaultJacobianMagnetostatics::DefaultJacobianMagnetostatics(int i, int j, Hermes::vector<std::string> areas,
							       HermesFunction* coeff,
                                                               SymFlag sym, GeomType gt,
                                                               int order_increase)
    : WeakForm::MatrixFormVol(i, j, areas, sym), idx_j(j), coeff(coeff), gt(gt),
    order_increase(order_increase) 
  { 
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  scalar DefaultJacobianMagnetostatics::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
      scalar planar_part = 0;
      scalar axisym_part = 0;
      for (int i = 0; i < n; i++) 
      {
        scalar B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
        if (std::abs(B_i) > 1e-12) 
        {
          planar_part += wt[i] * coeff->derivative(B_i) / B_i
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
            * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) 
          {
            axisym_part += wt[i] * coeff->derivative(B_i) / B_i / e->y[i]
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * u_ext[idx_j]->val[i] * v->dy[i];
          }
          else if (gt == HERMES_AXISYM_Y) 
          {
            axisym_part += wt[i] * coeff->derivative(B_i) / B_i / e->x[i]
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * u_ext[idx_j]->val[i] * v->dx[i];
          }
        }
        planar_part += wt[i] * coeff->value(B_i)
          * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        if (gt == HERMES_AXISYM_X) 
        {
          axisym_part += wt[i] * coeff->value(B_i) / e->y[i]
          * u->val[i] * v->dy[i];
        }
        else if (gt == HERMES_AXISYM_Y) 
        {
          axisym_part += wt[i] * coeff->value(B_i) / e->x[i]
          * u->val[i] * v->dx[i];
        }
      }

      return planar_part + axisym_part;
  }

  Ord DefaultJacobianMagnetostatics::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
      Ord planar_part = 0;
      for (int i = 0; i < n; i++) 
      {
        Ord B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
        planar_part += wt[i] * coeff->derivative(B_i) / B_i
          * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
          * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
        planar_part += wt[i] * coeff->value(B_i)
          * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }

      // This increase is for the axisymmetric part. We are not letting the
      // Ord class do it since it would automatically choose the highest order
      // due to the nonpolynomial 1/r term.
      return planar_part * Ord(order_increase);
  }

  WeakForm::MatrixFormVol* DefaultJacobianMagnetostatics::clone() 
  {
    return new DefaultJacobianMagnetostatics(*this);
  }


  DefaultResidualMagnetostatics::DefaultResidualMagnetostatics(int i, std::string area,
                                                               HermesFunction* coeff,
                                                               GeomType gt, int order_increase)
    : WeakForm::VectorFormVol(i, area), idx_i(i), coeff(coeff),
    gt(gt), order_increase(order_increase) 
  { 
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultResidualMagnetostatics::DefaultResidualMagnetostatics(int i, Hermes::vector<std::string> areas,
                                                               HermesFunction* coeff,
                                                               GeomType gt, int order_increase)
    : WeakForm::VectorFormVol(i, areas), idx_i(i), coeff(coeff), gt(gt),
    order_increase(order_increase) 
  { 
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  scalar DefaultResidualMagnetostatics::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<scalar> *ext) const {
      scalar planar_part = 0;
      scalar axisym_part = 0;
      for (int i = 0; i < n; i++) 
      {
        scalar B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
        planar_part += wt[i] * coeff->value(B_i) *
          (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        if (gt == HERMES_AXISYM_X) axisym_part += wt[i] * coeff->value(B_i) / e->y[i]
          * u_ext[idx_i]->val[i] * v->dy[i];
        else if (gt == HERMES_AXISYM_Y) axisym_part += wt[i] * coeff->value(B_i) / e->x[i]
          * u_ext[idx_i]->val[i] * v->dx[i];
      }
      return planar_part + axisym_part;
  }

  Ord DefaultResidualMagnetostatics::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord planar_part = 0;
    for (int i = 0; i < n; i++) 
    {
      Ord B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
      planar_part += wt[i] * coeff->value(B_i) *
        (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
    }
    return planar_part * Ord(order_increase);
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::VectorFormVol* DefaultResidualMagnetostatics::clone() 
  {
    return new DefaultResidualMagnetostatics(*this);
  }
};

