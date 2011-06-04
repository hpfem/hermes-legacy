#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_motor, double eps_motor, 
                        std::string mat_air, double eps_air);
};
