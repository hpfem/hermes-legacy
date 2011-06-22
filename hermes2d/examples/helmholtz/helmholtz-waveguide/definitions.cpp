#include "definitions.h"

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
  return 100*cos(y*M_PI/0.1);
}


WeakFormHelmholtz::WeakFormHelmholtz(double eps, double mu, double omega, double sigma, double beta, double E0, double h) : WeakForm(2)
    {
        add_matrix_form(new MatrixFormHelmholtzEquation_real_real(0, 0, eps, omega, mu));
        add_matrix_form(new MatrixFormHelmholtzEquation_real_imag(0, 1, mu, omega, sigma));
        add_matrix_form(new MatrixFormHelmholtzEquation_imag_real(1, 0, mu, omega, sigma));
        add_matrix_form(new MatrixFormHelmholtzEquation_imag_imag(1, 1, eps, mu, omega));

        add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_imag(0, 1, "Bdy_impedance", beta));
        add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_real(1, 0, "Bdy_impedance", beta));
    }

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
            return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
        }

scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }


        template<typename Real, typename Scalar>
        Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            return -omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
        }

  scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

 Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
 

 

        template<typename Real, typename Scalar>
        Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            return  omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
        }

   scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

 Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
 

 

        template<typename Real, typename Scalar>
        Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
        }
   scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

 Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }

 


 

        template<typename Real, typename Scalar>
        Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            return beta * int_u_v<Real, Scalar>(n, wt, u, v);
        }

  scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form_surf<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

     Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
 



 

        template<typename Real, typename Scalar>
        Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            return - beta*int_u_v<Real, Scalar>(n, wt, u, v);
        }

 scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form_surf<double, scalar>(n, wt, u_ext, u, v, e, ext);
        }

  Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }

