#include "hermes2d.h"

/* Weak forms */

class WeakFormEigenLeft : public WeakForm
{
public:
  WeakFormEigenLeft();
};

class WeakFormEigenRight : public WeakForm
{
public:
  WeakFormEigenRight();
};
