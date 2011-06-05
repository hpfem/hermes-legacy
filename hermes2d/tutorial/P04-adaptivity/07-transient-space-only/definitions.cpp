#include "definitions.h"

CustomNonlinearity::CustomNonlinearity(double alpha): HermesFunction()
  {
    this->is_const = false;
    this->alpha = alpha;
  }

scalar CustomNonlinearity::value(double u) const
  {
    return -1 - pow(u, alpha);
  }

Ord CustomNonlinearity::value(Ord u) const
  { 
    return Ord(10);
  }

scalar CustomNonlinearity::derivative(double u) const
  {
    return -alpha * pow(u, alpha - 1.0);
  }

Ord CustomNonlinearity::derivative(Ord u) const
  {
    return Ord(10);
  }


EssentialBCNonConst::EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

EssentialBoundaryCondition::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
  { 
    return EssentialBoundaryCondition::BC_FUNCTION; 
  }

scalar EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return (x + 10) * (y + 10) / 100.;
  }


void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = (y+10)/100.;
    dy = (x+10)/100.;
  };

scalar CustomInitialCondition::value (double x, double y) const 
  {
    return (x+10)*(y+10)/100.;
  }

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
  {
    return (x+10)*(y+10)/100.;
  }
