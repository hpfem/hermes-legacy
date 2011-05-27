////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

using namespace WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion; 

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   Hermes::vector<MeshFunction*>& iterates,
                   double init_keff, std::string bdy_vacuum);
};

// Integral over the active core.
double integrate(MeshFunction* sln, std::string area);