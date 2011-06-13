#include "hermes2d.h"

enum GalerkinMethod
{
  CG,
  CG_STAB_SUPG,    // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  CG_STAB_GLS_1,   // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  CG_STAB_GLS_2,   // assumes H2D_SECOND_DERIVATIVES_ENABLED, quadratic elements
  CG_STAB_SGS,     // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  CG_STAB_SGS_ALT, // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  DG
};

const std::string method_names[7] = 
{
  "unstabilized continuous Galerkin",
  "streamline upwind Petrov-Galerkin",
  "Galerkin least squares for linear elements",
  "Galerkin least squares for quadratic elements",
  "subgrid scale stabilized Galerkin",
  "alternative subgrid scale stabilized Galerkin",
  "discontinuous Galerkin"
};

//////////////////////////////////////////////////// PROBLEM DATA /////////////////////////////////////////////////////

// Diffusivity.
const double EPSILON = 0.01; 

// Constant flow field (\beta = (0,1), div(\beta) = 0).

// Inner product between two 2D vectors.
template<typename Real>
inline Real dot2(Real x1, Real y1, Real x2, Real y2) {
  return x1*x2 + y1*y2;
}

struct ConstFlowField
{  
  static const double b1;
  static const double b2;
  
  // Inner product between the flow field and a general vector at given point.
  template<typename Real>
  static Real dot(Real vx, Real vy) {
    return dot2<Real>(b1, b2, vx, vy);
  }
  
  // Inner product between the flow field and normal at given quadrature point.
  template<typename Real>
  static Real dot_n(const Geom<Real>* e, unsigned int pt) {
    return dot<Real>(e->nx[pt], e->ny[pt]);
  }
  
  // Inner product between the flow field and the gradient of a given function 
  // at given quadrature point.
  template<typename Real>
  static Real dot_grad(Func<Real> *u, const Geom<Real>* e, unsigned int pt) {
    return dot<Real>(u->dx[pt], u->dy[pt]);
  }
  
  static double norm() {
    return dot<double>(b1, b2);
  }
  
};

// Boundary condition represented by a non-vanishing function.
class NonzeroBoundaryValues : public ExactSolutionScalar 
{
  public:
    NonzeroBoundaryValues(Mesh* mesh) : ExactSolutionScalar(mesh) {};
    
    virtual scalar value(double x, double y) const {
      return sin(M_PI*x);
    }
    virtual Ord value(Ord x, Ord y) const {
      return ord(x, y);
    }
    virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
      dx = M_PI * cos(M_PI*x); dy = 0.0;      
    }
    virtual Ord ord(Ord x, Ord y) const {
      return Ord(10);
    }
};

// This class allows Hermes to treat essential boundary conditions either in a standard way
// (automatically during construction of the space) or weakly (by integrating their values 
// in a user-defined integral form). The container EssentialBCs may be used to pass all
// weakly imposable essential b.c.'s of the problem to the space constructor, or to obtain
// the b.c. associated with currently assembled boundary edge inside the integral form.
// In the latter case, the retrieved pointer to EssentialBoundaryCondition must be 
// static-cast to WeaklyImposableBC in order to correctly handle calls to 'value' with both
// double and Ord arguments.
class WeaklyImposableBC : public EssentialBoundaryCondition
{  
  EssentialBCValueType type;
  ExactSolutionScalar* exact_solution;
  
  public:
    WeaklyImposableBC(const Hermes::vector<std::string>& markers, scalar value) 
      : EssentialBoundaryCondition(markers), type(EssentialBoundaryCondition::BC_CONST)
    {
      this->value_const = value;  
    };
    WeaklyImposableBC(const Hermes::vector<std::string>& markers, ExactSolutionScalar* exact_solution) 
      : EssentialBoundaryCondition(markers), type(EssentialBoundaryCondition::BC_FUNCTION)
    {
      this->exact_solution = exact_solution;
    }
    
    inline EssentialBCValueType get_value_type() const { return type; }
    
    virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
      return value(x, y);
    }
    
    virtual scalar value(double x, double y) const {
      return (type == BC_CONST) ? value_const : exact_solution->value(x, y);
    }
    virtual Ord value(Ord x, Ord y) const {
      return (type == BC_CONST) ? Ord(1) : exact_solution->ord(x, y);
    }
};

//////////////////////////////////////////////////// EXACT SOLUTION ///////////////////////////////////////////////////

class CustomExactSolution : public ExactSolutionScalar
{
  double D, m1, m2, ie1me2; // Auxiliary variables for calculating the value/derivatives.
  
  public:
    CustomExactSolution(Mesh* mesh, double epsilon) : ExactSolutionScalar(mesh) {
      D = sqrt(1.+4.*sqr(epsilon*M_PI));
      m1 = (1.-D) / (2.*epsilon);
      m2 = (1.+D) / (2.*epsilon);
      ie1me2 = 1. / (exp(m1) - exp(m2));
    };
    
    virtual scalar value(double x, double y) const {
      return ie1me2 * (exp(m1+m2*y) - exp(m2+m1*y)) * sin(M_PI*x);
    }
    
    virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
      double e12y = exp(m1+m2*y);
      double e21y = exp(m2+m1*y);
      dx = ie1me2 * (e12y - e21y) * M_PI * cos(M_PI*x);
      dy = ie1me2 * (e12y*m2 - e21y*m1) * sin(M_PI*x);
    }
    
    virtual Ord ord(Ord x, Ord y) const {
      return Ord(Ord::get_max_order());
    }
};

////////////////////////////////////////////////////// WEAK FORMS /////////////////////////////////////////////////////

/////////////////////////////////////////////// CONTINUOUS APPROXIMATION //////////////////////////////////////////////
//
// References:
//
// R. Codina:
// Comparison of some finite element methods for solving the diffusion-convection-reaction equation.
// Comput. Meth. Appl. Mech. Engrg. 156 (1998) 185-210.
//
// F. Shakib, T.J.R. Hughes, Z. Johan: 
// A new finite element formulation for computational fluid dynamics: X. The compressible Euler and Navier-Stokes 
// equations. Comput. Methods Appl. Mech. Engrg. 89 (1991) 141-219.
//

class CustomWeakFormContinuousGalerkin : public WeakForm
{
  HermesFunction *fn_epsilon, *fn_b1, *fn_b2;
  
  public:
    CustomWeakFormContinuousGalerkin(GalerkinMethod method, double epsilon);
    ~CustomWeakFormContinuousGalerkin() {
      delete fn_epsilon; delete fn_b1; delete fn_b2;
    }
                                     
  private:
    
    class StabilizationJacobian : public WeakForm::MatrixFormVol
    {
      GalerkinMethod method;
      double epsilon, b_norm;
      
      double element_Peclet_number(double h_e) const { 
        return b_norm * h_e / (2*epsilon);
      }
      Ord element_Peclet_number(Ord h_e) const { 
        return h_e;
      }
      
      public:
        StabilizationJacobian(GalerkinMethod method, double epsilon) 
          : WeakForm::MatrixFormVol(0, 0, HERMES_ANY, 
                                    (method==CG_STAB_GLS_1 || method==CG_STAB_GLS_2) ? HERMES_SYM : HERMES_NONSYM),
            method(method), epsilon(epsilon), b_norm(ConstFlowField::norm())
        {
        #ifndef H2D_SECOND_DERIVATIVES_ENABLED   
            if (method != CG && method != DG)
              error("Stabilization forms require H2D_SECOND_DERIVATIVES_ENABLED defined in h2d_common.h.");
        #endif    
        };

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                            Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

        virtual WeakForm::MatrixFormVol* clone() {
          return new StabilizationJacobian(*this);
        }
    };
    
    class StabilizationResidual : public WeakForm::VectorFormVol
    {
      GalerkinMethod method;
      double epsilon, b_norm;
      
      double element_Peclet_number(double h_e) const { 
        return b_norm * h_e / (2*epsilon);
      }
      Ord element_Peclet_number(Ord h_e) const { 
        return h_e;
      }
      
      public:
        StabilizationResidual(GalerkinMethod method, double epsilon)
          : WeakForm::VectorFormVol(0), method(method), epsilon(epsilon), b_norm(ConstFlowField::norm())
        {};
        
        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                            Geom<double> *e, ExtData<scalar> *ext) const;
        
        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

        virtual WeakForm::VectorFormVol* clone() {
          return new StabilizationResidual(*this);
        }
    };
};

///////////////////////////////////////////// DISCONTINUOUS APPROXIMATION /////////////////////////////////////////////

// Advection : upwind flux
// Diffusion : penalized interface jumps
//              - SIPG (Symmetric Interior Penalty Galerkin): theta = -1, C_W nontrivial, optimal convergence rate
//              - IIPG (Incomplete Interior Penalty Galerkin): theta = 0, C_W nontrivial, sub-optimal convergence rate
//              - NIPG (Nonsymmetric Interior Penalty Galerkin): theta = 1, C_W arbitrary, sub-optimal convergence rate
// Ess. b.c.:  weakly imposed (as natural conditions). 
//

class CustomWeakFormDiscontinuousGalerkin : public WeakForm
{
  int theta, C_W;
  
  HermesFunction *fn_epsilon;
  
  public:
    CustomWeakFormDiscontinuousGalerkin(const EssentialBCs& boundary_values,
                                        double epsilon, int theta = 1, int C_W = 1);
    ~CustomWeakFormDiscontinuousGalerkin() { delete fn_epsilon; }
      
  private:
    
    static scalar upwind_flux(double u_cent, double u_neib, double beta_dot_n) {
      scalar eff_val;
      if (beta_dot_n > 0) eff_val = u_cent;
      else if (beta_dot_n < 0) eff_val = u_neib;
      else eff_val = 0.5*(u_cent + u_neib);
      return beta_dot_n * eff_val;
    }
    
    static Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) {
      return a_dot_n * (u_cent + u_neib); 
    }
    
    struct Advection
    {
      class VolumetricJacobian : public WeakForm::MatrixFormVol
      {
        public:
          VolumetricJacobian() : WeakForm::MatrixFormVol(0,0) {};
          
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
            return new VolumetricJacobian(*this);
          }
      };
      
      class VolumetricResidual : public WeakForm::VectorFormVol
      {
        public:
          VolumetricResidual() : WeakForm::VectorFormVol(0) {};
          
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
            return new VolumetricResidual(*this);
          }
      };
      
      class BoundaryJacobian : public WeakForm::MatrixFormSurf
      {
        public:
          BoundaryJacobian() : WeakForm::MatrixFormSurf(0,0) {};
          
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

          virtual WeakForm::MatrixFormSurf* clone() {
            return new BoundaryJacobian(*this);
          }
      };
      
      class BoundaryResidual : public WeakForm::VectorFormSurf
      {
        const EssentialBCs& boundary_values;
        
        public:
          BoundaryResidual(const EssentialBCs& boundary_values) : WeakForm::VectorFormSurf(0), 
            boundary_values(boundary_values) 
          {};
          
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

          virtual WeakForm::VectorFormSurf* clone() {
            return new BoundaryResidual(*this);
          }
      };
      
      class InterfaceJacobian : public WeakForm::MatrixFormSurf
      {
        public:
          InterfaceJacobian() : WeakForm::MatrixFormSurf(0,0,H2D_DG_INNER_EDGE) {};
          
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

          virtual WeakForm::MatrixFormSurf* clone() {
            return new InterfaceJacobian(*this);
          }
      };
      
      class InterfaceResidual : public WeakForm::VectorFormSurf
      {
        public:
          InterfaceResidual() : WeakForm::VectorFormSurf(0,H2D_DG_INNER_EDGE) {};
          
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

          virtual WeakForm::VectorFormSurf* clone() {
            return new InterfaceResidual(*this);
          }
      };
    };
    
    struct Diffusion
    {
      class BoundaryJacobian : public WeakForm::MatrixFormSurf
      {
        double epsilon;
        int theta, C_W;
        
        public:
          BoundaryJacobian(double epsilon, int theta, int C_W) : WeakForm::MatrixFormSurf(0,0),
            epsilon(epsilon), theta(theta), C_W(C_W)
          {};
          
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

          virtual WeakForm::MatrixFormSurf* clone() {
            return new BoundaryJacobian(*this);
          }
      };
      
      class BoundaryResidual : public WeakForm::VectorFormSurf
      {
        const EssentialBCs& boundary_values;
        
        double epsilon;
        int theta, C_W;
        
        public:
          BoundaryResidual(const EssentialBCs& boundary_values, double epsilon, int theta, int C_W) 
            : WeakForm::VectorFormSurf(0), boundary_values(boundary_values), epsilon(epsilon), theta(theta), C_W(C_W)
          {};
          
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

          virtual WeakForm::VectorFormSurf* clone() {
            return new BoundaryResidual(*this);
          }
      };
      
      class InterfaceJacobian : public WeakForm::MatrixFormSurf
      {
        double epsilon;
        int theta, C_W;
        
        public:
          InterfaceJacobian(double epsilon, int theta, int C_W) 
            : WeakForm::MatrixFormSurf(0,0,H2D_DG_INNER_EDGE), epsilon(epsilon), theta(theta), C_W(C_W) 
          {};
          
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

          virtual WeakForm::MatrixFormSurf* clone() {
            return new InterfaceJacobian(*this);
          }
      };
      
      class InterfaceResidual : public WeakForm::VectorFormSurf
      {
        double epsilon;
        int theta, C_W;
        
        public:
          InterfaceResidual(double epsilon, int theta, int C_W) 
            : WeakForm::VectorFormSurf(0,H2D_DG_INNER_EDGE), epsilon(epsilon), theta(theta), C_W(C_W)
          {};
          
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

          virtual WeakForm::VectorFormSurf* clone() {
            return new InterfaceResidual(*this);
          }
      };
    };
};
