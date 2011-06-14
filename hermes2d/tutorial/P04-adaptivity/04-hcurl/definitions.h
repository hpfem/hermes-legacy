#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/hcurl.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/weakforms_hcurl.h"


class CustomExactSolution : public ExactSolutionVector
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionVector(mesh) {};
  ~CustomExactSolution() {};

  virtual scalar2 value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar2& dx, scalar2& dy) const;
  
  virtual Ord ord(Ord x, Ord y) const;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double mu_r, double kappa);

  class CustomVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurf() : WeakForm::VectorFormSurf(0) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  };
};

