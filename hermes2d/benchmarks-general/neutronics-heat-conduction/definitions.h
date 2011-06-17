#ifndef __DEFINITIONS_H__INCLUDED__
#define __DEFINITIONS_H__INCLUDED__

#include "hermes2d.h"

struct FuelProperties
{
  static const double Tref;                   // Reference temperature for defining temperature dependence of parameters.
  static const double neutron_velocity;       // Averaged neutron speed.
  static const double D;                      // Neutron diffusion coefficient.
  static const double nu;                     // Number of neutrons emitted per fission event.
  static const double Sigma_f;                // Fission cross section.
  static const double kappa;                  // Energy per fission.
  static const double rho;                    // Fuel density.
  static const double cp;                     // Fuel heat capacity.
  
  // Thermal conductivity dependence on temperature.
  static const double k0;  // Constant coefficient of heat conduction dependence on temperature.
  static const double k1;  // Linear coefficient of heat conduction dependence on temperature.                    
  
  class NegativeConductivityTemperatureProfile : public HermesFunction
  {                                           
    public:                                   
      virtual double value(double T) const    
      {                                       
        return -(k0 + k1 * (T - Tref))/(rho*cp);
      }                                       
      
      virtual Ord value(Ord T) const 
      {
        return T;
      }
      
      virtual double derivative(double T) const
      {
        return -k1/(rho*cp);
      }
      
      virtual Ord derivative(Ord T) const 
      {
        return Ord(1);
      }
  };
  
  // Removal cross section dependence on temperature.
  static const double Sigma_r_ref;            // Absorption cross section for reference temperature Tref.
  static const double doppler_coeff;          // Coefficient for modelling the dependence of absorption cross section
                                              // on temperature (Doppler effect).
  class RemovalCrossSectionTemperatureProfile : public HermesFunction
  {
    public:      
      virtual double value(double T) const
      {
        return neutron_velocity * (Sigma_r_ref + doppler_coeff * (sqrt(T) - sqrt(Tref)));
      }
      
      virtual Ord value(Ord T) const 
      {
        return sqrt(T);
      }
      
      virtual double derivative(double T) const
      {
        return neutron_velocity * doppler_coeff / (2*sqrt(T));
      }
      
      virtual Ord derivative(Ord T) const 
      {
        return Ord(10);
      }
  };
};

struct TemperatureField
{ 
  static const double CT = 1.0;
  static const double rT = 1.0;
  
  // Time dependence of the temperature.
  static double transient_profile(double time) 
  {
    //  return 1.0;
    return 1+tanh(rT*time);
  }
  
  static double transient_profile_derivative(double time) 
  {
    //  return 0.0;
    return rT*(1-pow(tanh(rT*time),2));
  }
  
  
  class ExactDistribution : public ExactSolutionScalar
  {
    double LX, LY;  // Domain extents in the x and y directions.
    double time;    // Current time.
        
    public:
      ExactDistribution(Mesh* mesh, double lx, double ly) : ExactSolutionScalar(mesh), LX(lx), LY(ly), time(0) {};
      
      virtual scalar value(double x, double y) const 
      {
        return CT*transient_profile(time)*sin(M_PI*x/LX)*sin(M_PI*y/LY);
      }
      
      virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
      {
        dx = CT*transient_profile(time)*M_PI/LX*cos(M_PI*x/LX)*sin(M_PI*y/LY);
        dy = CT*transient_profile(time)*sin(M_PI*x/LX)*M_PI/LY*cos(M_PI*y/LY);
      }
      
      virtual Ord ord(Ord x, Ord y) const 
      {
        return Ord(Ord::get_max_order());
      }
      
      void update(double new_time) { time = new_time; reinit(); }
  };
  
  class SourceTerm : public HermesFunction
  {
    double LX, LY; // Domain extents in the x and y directions.
    
    FuelProperties fuel;
    
    public:
      SourceTerm(double lx, double ly) : HermesFunction(), LX(lx), LY(ly) {};
      
      template<typename Real> Real val(Real x, Real y, double t) const;
      
      virtual double value(double x, double y, double t) const
      {
        return val<double>(x, y, t);
      }
      
      virtual Ord value(Ord x, Ord y, double t) const 
      {
        //return val<Ord>(x, y, t);
        return Ord(10);
      }
  };
    
  class Residual : public WeakForm::VectorFormVol
  {
    FuelProperties::NegativeConductivityTemperatureProfile minus_lambda;
    SourceTerm qT;
    double fission_source_coeff;
    
    public:
      Residual(double lx, double ly) : WeakForm::VectorFormVol(0), qT(lx, ly)
      {
        fission_source_coeff = FuelProperties::kappa * FuelProperties::Sigma_f / (FuelProperties::rho * FuelProperties::cp);
      };
      
      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext) const
      {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }                    

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }
      
      virtual WeakForm::VectorFormVol* clone() {
        return new Residual(*this);
      }
  };
};

struct NeutronField
{
  static const double CF = 1.0;
  static const double rF = 0.25;
  
  // Time dependence of the neutron flux.
  static double transient_profile(double time) 
  {
    //  return Temperature::transient_profile(time);
    return 1+exp(rF*time);
  }
  
  static double transient_profile_derivative(double time) 
  {
    //  return Temperature::transient_profile_derivative(time);
    return rF*exp(rF*time);
  }
  
  class ExactDistribution : public ExactSolutionScalar
  {
    double LX, LY;  // Domain extents in the x and y directions.
    double time;    // Current time.
    
    public:
      ExactDistribution(Mesh* mesh, double lx, double ly) : ExactSolutionScalar(mesh), LX(lx), LY(ly), time(0) {};
      
      virtual scalar value(double x, double y) const 
      {
        return CF*transient_profile(time)*sin(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX*y/LY;
      }
      
      virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
      {
        dx = CF*transient_profile(time)*(sin(M_PI*x/LX)*sin(M_PI*y/LY)/LX*y/LY
              + M_PI/LX*cos(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX*y/LY);
        dy = CF*transient_profile(time)*(sin(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX/LY
              + sin(M_PI*x/LX)*M_PI/LY*cos(M_PI*y/LY)*x/LX*y/LY);
      }
      
      virtual Ord ord(Ord x, Ord y) const 
      {
        return Ord(Ord::get_max_order());
      }
      
      void update(double new_time) { time = new_time; reinit(); }
  };
  
  class SourceTerm : public HermesFunction
  {
    double LX, LY; // Domain extents in the x and y directions.
    
    FuelProperties fuel;
    
    public:
      SourceTerm(double lx, double ly) : HermesFunction(), LX(lx), LY(ly) {};
      
      template<typename Real> Real val(Real x, Real y, double t) const;
      
      virtual double value(double x, double y, double t) const
      {
        return val<double>(x, y, t);
      }
      
      virtual Ord value(Ord x, Ord y, double t) const 
      {
        //return val<Ord>(x, y, t);
        return Ord(10);
      }
  };
  
  class TemperatureDerivative : public WeakForm::MatrixFormVol
  {
    FuelProperties::RemovalCrossSectionTemperatureProfile Sigma_r;
    
    public:
      TemperatureDerivative() : WeakForm::MatrixFormVol(0, 0) {};
        
      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
      {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }                      

      virtual WeakForm::MatrixFormVol* clone() {
        return new TemperatureDerivative(*this);
      }
  };
  
  class NeutronFluxDerivative : public WeakForm::MatrixFormVol
  {  
    FuelProperties::RemovalCrossSectionTemperatureProfile Sigma_r;
    double diffusion_coeff, fission_yield_coeff;
    
    public:
      NeutronFluxDerivative() : WeakForm::MatrixFormVol(1, 1, HERMES_ANY, HERMES_NONSYM) 
      {
        diffusion_coeff = FuelProperties::neutron_velocity * FuelProperties::D;
        fission_yield_coeff = FuelProperties::neutron_velocity * FuelProperties::nu * FuelProperties::Sigma_f;
      };
        
      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
      {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }                      

      virtual WeakForm::MatrixFormVol* clone() {
        return new NeutronFluxDerivative(*this);
      }
  };
  
  class Residual : public WeakForm::VectorFormVol
  {
    FuelProperties::RemovalCrossSectionTemperatureProfile Sigma_r;
    SourceTerm qF;
    double diffusion_coeff, fission_yield_coeff;
    
    public:
      Residual(double lx, double ly) : WeakForm::VectorFormVol(1), qF(lx, ly)
      {
        diffusion_coeff = FuelProperties::neutron_velocity * FuelProperties::D;
        fission_yield_coeff = FuelProperties::neutron_velocity * FuelProperties::nu * FuelProperties::Sigma_f;
      };
      
      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext) const
      {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }                    

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }
      
      virtual WeakForm::VectorFormVol* clone() {
        return new Residual(*this);
      }
  };
};

class CustomWeakForm : public WeakForm
{
  public:
    CustomWeakForm(double lx, double ly);
};

class Views
{
  protected:
    char title[100]; // Character array to store the title for an actual view and time step.
    
    ScalarView *sview_T;
    ScalarView *sview_phi;
    ScalarView *sview_T_exact;
    ScalarView *sview_phi_exact;
  
  public:
    Views(unsigned int width, unsigned int height);
    virtual ~Views();
    
    virtual void show_solutions(double current_time, Hermes::vector< Solution* > solutions);
    virtual void show_exact(double current_time, Hermes::vector< Solution* > exact);
};

#endif