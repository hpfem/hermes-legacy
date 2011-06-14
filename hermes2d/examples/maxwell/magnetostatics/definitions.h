#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm
{
public:
  CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2,
                               CubicSpline* mu_inv_iron, std::string material_air,
                               std::string material_copper, double mu_vacuum,
                               double current_density, int order_inc = 3);
};

class HERMES_API FilterVectorPotencial : public MagFilter
{
public:
  FilterVectorPotencial(Hermes::vector<MeshFunction*> solutions, Hermes::vector<int> items) 
        : MagFilter(solutions, items) {};

protected:
  void filter_fn(int n, Hermes::vector<scalar*> values, scalar* result, Geom<double> *e);
};

class HERMES_API FilterFluxDensity : public Filter
{
public:
  FilterFluxDensity(Hermes::vector<MeshFunction*> solutions)
        : Filter(solutions) {};

  virtual scalar get_pt_value(double x, double y, int item = H2D_FN_VAL);

protected:
  void precalculate(int order, int mask);
};
