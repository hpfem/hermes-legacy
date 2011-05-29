#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_al, HermesFunction* lambda_al,
                        std::string mat_cu, HermesFunction* lambda_cu,
                        HermesFunction* src_term);
};
