#include "hermes2d.h"

enum DiscretizationMethod
{
  CG,
  SUPG,
  DG
};

const std::string method_names[3] = 
{
  "unstabilized continuous Galerkin",
  "streamline upwind Petrov-Galerkin",
  "discontinuous Galerkin",
};

class InflowVariation : public HermesFunction
{
  public:
    InflowVariation() : HermesFunction() {
      is_const = false;
    };
    
    virtual scalar value(double x, double y) const {
      return sqr(sin(M_PI*y));
    }
    virtual Ord value(Ord x, Ord y) const {
      //return sqr(sin(M_PI*y));
      // Hack to remove the warning about overcoming the maximal order:
      return Ord(15);
    }
};

class ReactionTerm : public HermesFunction
{
  public:
    ReactionTerm(DiscretizationMethod method) : HermesFunction() {
      const_value = (method == DG) ? 11. : 0.; // DG: 11 = -div(beta)
    };
};

// Inner product between two 2D vectors.
template<typename Real>
inline Real dot2(Real x1, Real y1, Real x2, Real y2) {
  return x1*x2 + y1*y2;
}

// Class that represents the flow field and reaction variables and associated operations.
template<typename Real> struct FlowField
{
  static Real fn_a(Real x, Real y) {
    return 10.*sqr(y) - 12.*x + 1.;
  }
  static Real fn_b(Real x, Real y) {
    return 1. + y;
  }
  
  // Inner product between the flow field and a general vector at given point.
  static Real dot(Real x, Real y, Real vx, Real vy) {
    return dot2<Real>(fn_a(x,y), fn_b(x,y), vx, vy);
  }
  
  // Inner product between the flow field and normal at given quadrature point.
  static Real dot_n(const Geom<Real>* e, unsigned int pt) {
    return dot(e->x[pt], e->y[pt], e->nx[pt], e->ny[pt]);
  }
  
  // Inner product between the flow field and the gradient of a given function 
  // at given quadrature point.
  static Real dot_grad(Func<Real> *u, const Geom<Real>* e, unsigned int pt) {
    return dot(e->x[pt], e->y[pt], u->dx[pt], u->dy[pt]);
  }
};

// Weighted integral over the specified boundary edge. A semi-analytic value is used 
// as the exact value for benchmarking purposes.
class BoundaryIntegral
{
  std::string marker;
  FlowField<double> beta;
  
  // Weighting function.
  double weight_fn(double x, double y) const {
    return sin(M_PI*x*0.5);
  }
  
  public:
    BoundaryIntegral(std::string bnd_marker) : marker(bnd_marker) {};  
    double value(MeshFunction* sln) const;
    static double exact() { return 0.246500343856481; }
};

// Class that represents the semi-analytic solution on the rectangle [0,1]x[0,1].
// Contains coordinates of the Gauss-Kronrod quadrature nodes and the corresponding
// function values and quadrature weights, as read from the file supplied to the
// constructor.
class SemiAnalyticSolution
{  
  std::vector<long double> x; //  x coordinate.
  std::vector<long double> y; //  y coordinate.
  std::vector<long double> u; //  u(x,y).
  std::vector<long double> w; //  Gauss-Kronrod quadrature weights.
  
  unsigned long int n;  // Number of quadrature points (size of vectors x,y,u,w).
  
  bool NA;  // Indicates that the exact solution could not be read from the input file.
  
  public:
    SemiAnalyticSolution(std::string file);
    
    // The following two functions use the quadrature specified by the input file to
    // calculate integrals over the rectangle [0,1]x[0,1].
    //
    // Calculates L2-norm of the exact solution.
    double get_l2_norm();     
    //
    // Calculates L2-norm of the relative difference between
    // the exact solution and the one computed by Hermes.  
    double get_l2_rel_err(Solution *approx_sln);  
};

class CustomWeakFormContinuousGalerkin : public WeakForm
{
  public:
    CustomWeakFormContinuousGalerkin(bool supg_stabilization = true);
    
  private:
    class VolumetricJacobian : public WeakForm::MatrixFormVol
    {
      ReactionTerm *mu;
      
      public:
        VolumetricJacobian(DiscretizationMethod disc_method) : WeakForm::MatrixFormVol(0,0), 
          mu(new ReactionTerm(disc_method)) 
        {};
        
        ~VolumetricJacobian() { delete mu; }
          
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
      ReactionTerm *mu;
      
      public:
        VolumetricResidual(DiscretizationMethod disc_method) : WeakForm::VectorFormVol(0), 
          mu(new ReactionTerm(disc_method)) 
        {};
        
        ~VolumetricResidual() { delete mu; }
        
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
    
    class InflowBoundaryJacobian : public WeakForm::MatrixFormSurf
    {
      const HermesFunction *inflow;
      
      public:
        InflowBoundaryJacobian(const std::string& boundary_marker, const HermesFunction *boundary_function)
          : WeakForm::MatrixFormSurf(0,0,boundary_marker), inflow(boundary_function)
        {};
        
        ~InflowBoundaryJacobian() { delete inflow; }
        
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
          return new InflowBoundaryJacobian(*this);
        }
    };
    
    class InflowBoundaryResidual : public WeakForm::VectorFormSurf
    {
      const HermesFunction *inflow;
      
      public:
        InflowBoundaryResidual(const std::string& boundary_marker, const HermesFunction *boundary_function)
          : WeakForm::VectorFormSurf(0,boundary_marker), inflow(boundary_function)
        {};
        
        ~InflowBoundaryResidual() { delete inflow; }
        
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
          return new InflowBoundaryResidual(*this);
        }
    };
    
    class StabilizationJacobian : public WeakForm::MatrixFormVol
    {
      public:
        StabilizationJacobian() : WeakForm::MatrixFormVol(0,0) {};
        
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
          // return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
          // Hack to remove the warning about overcoming the maximal order:
          return Ord(15);
        } 

        virtual WeakForm::MatrixFormVol* clone() {
          return new StabilizationJacobian(*this);
        }
    };
    
    class StabilizationResidual : public WeakForm::VectorFormVol
    {
      public:
        StabilizationResidual() : WeakForm::VectorFormVol(0) {};
        
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
          // vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
          // Hack to remove the warning about overcoming the maximal order:
          return Ord(15); 
        }

        virtual WeakForm::VectorFormVol* clone() {
          return new StabilizationResidual(*this);
        }
    };
};

class CustomWeakFormDiscontinuousGalerkin : public WeakForm
{
  public:
    // Theta = 0.5 for upwinding
    CustomWeakFormDiscontinuousGalerkin(double theta = 0.5);
    
  private :
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
      const HermesFunction *inflow;
      
      public:
        BoundaryJacobian(const std::string& boundary_marker, const HermesFunction *boundary_function = NULL)
          : WeakForm::MatrixFormSurf(0,0,boundary_marker), inflow(boundary_function)
        {};
        
        ~BoundaryJacobian() { delete inflow; }
        
        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                            Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

        virtual WeakForm::MatrixFormSurf* clone() {
          return new BoundaryJacobian(*this);
        }
    };
    
    class BoundaryResidual : public WeakForm::VectorFormSurf
    {
      const HermesFunction *inflow;
      
      public:
        BoundaryResidual(const std::string& boundary_marker, const HermesFunction *boundary_function = NULL)
          : WeakForm::VectorFormSurf(0,boundary_marker), inflow(boundary_function) 
        {};
        
        ~BoundaryResidual() { delete inflow; }
        
        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                            Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;

        virtual WeakForm::VectorFormSurf* clone() {
          return new BoundaryResidual(*this);
        }
    };
    
    class InterfaceJacobian : public WeakForm::MatrixFormSurf
    {
      double theta;   // Stabilization parameter. Standard upwind scheme is obtained for theta = 0.5.
      
      public:
        InterfaceJacobian(double theta = 0.5) : WeakForm::MatrixFormSurf(0,0,H2D_DG_INNER_EDGE),
          theta(theta) 
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
      double theta;   // Stabilization parameter. Standard upwind scheme is obtained for theta = 0.5.
      
      public:
        InterfaceResidual(double theta = 0.5) : WeakForm::VectorFormSurf(0,H2D_DG_INNER_EDGE),
          theta(theta) 
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