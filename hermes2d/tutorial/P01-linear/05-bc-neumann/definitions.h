#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoissonNeumann : public WeakForm
{
public:
  CustomWeakFormPoissonNeumann(std::string mat_al, HermesFunction* lambda_al,
                               std::string mat_cu, HermesFunction* lambda_cu,
                               HermesFunction* vol_src_term, std::string bdy_heat_flux,
                               HermesFunction* surf_src_term);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition 
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  protected:
    double A, B, C;
};

