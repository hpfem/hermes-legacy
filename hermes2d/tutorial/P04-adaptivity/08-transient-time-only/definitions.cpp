#include "definitions.h"

WeakFormHeatTransferNewtonTimedep::WeakFormHeatTransferNewtonTimedep(double alpha, double tau, Solution* sln_prev_time) : WeakForm(1) 
{ 
  add_matrix_form(new MatrixFormVolHeatTransfer(0, 0, alpha, tau));
  add_vector_form(new VectorFormVolHeatTransfer(0, alpha, tau));
}

template<typename Real, typename Scalar>
Scalar WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                                                 Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result -= wt[i] * (dlam_du<Real>(u_prev_newton->val[i]) * u->val[i] *
                       (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + lam<Real>(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

scalar WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                                           Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                                      Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real>
Real WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::lam(Real u) const 
{ 
  return 1 + pow(u, alpha); 
}

template<typename Real>
Real WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::dlam_du(Real u) const 
{ 
  return alpha*pow(u, alpha - 1); 
}
    
WeakForm::MatrixFormVol* WeakFormHeatTransferNewtonTimedep::MatrixFormVolHeatTransfer::clone() 
{
  return new MatrixFormVolHeatTransfer(*this);
}


template<typename Real, typename Scalar>
Scalar WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                                                 Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result -= wt[i] * (lam<Real>(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] 
                       + u_prev_newton->dy[i] * v->dy[i])
	           - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

scalar WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                        Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                   ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

template<typename Real>
Real WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::heat_src(Real x, Real y) const 
{
  return 1.0;
}

template<typename Real>
Real WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::lam(Real u) const 
{ 
  return 1 + pow(u, alpha); 
}

WeakForm::VectorFormVol* WeakFormHeatTransferNewtonTimedep::VectorFormVolHeatTransfer::clone() 
{
  return new VectorFormVolHeatTransfer(*this);
}


EssentialBCNonConst::EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
{
  markers.push_back(marker);
}

EssentialBoundaryCondition::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
{
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

scalar EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
{
  return (x+10)*(y+10)/100.;
}

scalar InitialSolutionHeatTransfer::value (double x, double y) const 
{
  return (x+10)*(y+10)/100.;
} 

void InitialSolutionHeatTransfer::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (y+10)/10.;
  dy = (x+10)/10.;
}

Ord InitialSolutionHeatTransfer::ord(Ord x, Ord y) const 
{
  return (x+10)*(y+10)/100.;
}
