#include "hermes2d.h"

/* Weak forms */

class WeakFormHeatTransferNewtonTimedep : public WeakForm
{
public:
  WeakFormHeatTransferNewtonTimedep(double alpha, double tau, Solution* sln_prev_time);

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, double alpha, double tau) 
          : WeakForm::MatrixFormVol(i, j), alpha(alpha), tau(tau) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const;

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) const;
    
    double alpha;
    double tau;
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i, double alpha, double tau) 
          : WeakForm::VectorFormVol(i), alpha(alpha), tau(tau) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) const;

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const;

    double alpha;
    double tau;
  };
};

/* Essential boundary conditions */

class EssentialBCNonConst : public EssentialBoundaryCondition {
public:
  EssentialBCNonConst(std::string marker);

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

/* Initial condition */

class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};
