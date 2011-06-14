#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormAcoustics : public WeakForm
{
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                          double sound_speed, double omega);
};
