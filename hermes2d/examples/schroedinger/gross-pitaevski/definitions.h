#include "hermes2d.h"
#include "runge_kutta.h"

/* Initial solution */

/*
class CustomInitialSolution : public ExactSolutionScalar
{
public:
  CustomInitialSolution(Mesh* mesh)
               : ExactSolutionScalar(mesh) {};

  virtual double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};
*/

/* Weak forms */

class CustomWeakFormGPRK : public WeakForm
{
public:
  CustomWeakFormGPRK(double H, double M, double G, double OMEGA);

private:

  class CustomFormMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomFormMatrixFormVol(int i, int j, double H, double M, double G, double OMEGA) 
          : WeakForm::MatrixFormVol(i, j), H(H), M(M), G(G), OMEGA(OMEGA) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_rk(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                    Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    // Members.
    double H, M, G, OMEGA;
  };


  class CustomFormVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomFormVectorFormVol(int i, double H, double M, double G, double OMEGA)
          : WeakForm::VectorFormVol(i), H(H), M(M), G(G), OMEGA(OMEGA) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    // Members.
    double H, M, G, OMEGA;
  };
  double H, M, G, OMEGA;
};

