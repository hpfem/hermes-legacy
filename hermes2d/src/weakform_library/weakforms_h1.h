// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_H1_WEAK_FORMS_H
#define __H2D_H1_WEAK_FORMS_H

#include "../integrals/h1.h"

namespace WeakFormsH1 
{
  /* Default volumetric matrix form \int_{area} coeff(x, y) * u * v \bfx
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    DefaultMatrixFormVol(int i = 0, int j = 0, std::string area = HERMES_ANY,
                         HermesFunction* coeff = HERMES_ONE,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
                         HermesFunction* coeff = HERMES_ONE,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultMatrixFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default volumetric matrix form \int_{area} coeff'(u_ext[0]) u \nabla u_ext[0] \cdot \nabla v
     + coeff(u_ext[0]) * \nabla u \cdot \nabla v d\bfx
     coeff... (generally nonconstant) function of the solution
  */

  class HERMES_API DefaultJacobianDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianDiffusion(int i = 0, int j = 0, std::string area = HERMES_ANY, 
                             HermesFunction* coeff = HERMES_ONE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas,
                             HermesFunction* coeff = HERMES_ONE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultJacobianDiffusion();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      int idx_j;
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default volumetric matrix form
     \int_{area} coeff`(u_ext[0]) * u * u_ext[0]->dx * v
     + coeff1(u_ext[0]) * u->dx * v
     + coeff2`(u_ext[0]) * u * u_ext[0]->dy * v
     + coeff2(u_ext[0]) * u->dy * v d\bfx.
     coeff1, coeff2... generally non-constant functions of the solution
  */

  class HERMES_API DefaultJacobianAdvection : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianAdvection(int i = 0, int j = 0, std::string area = HERMES_ANY, 
                             HermesFunction* coeff1 = HERMES_ONE,
                             HermesFunction* coeff2 = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

    DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas,
                             HermesFunction* coeff1 = HERMES_ONE,
                             HermesFunction* coeff2 = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

   ~DefaultJacobianAdvection();

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
      int idx_j;
      HermesFunction* coeff1, *coeff2;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} coeff(x, y) * v d\bfx
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultVectorFormVol(int i = 0, std::string area = HERMES_ANY,
                         HermesFunction* coeff = HERMES_ONE,
                         GeomType gt = HERMES_PLANAR);

    DefaultVectorFormVol(int i, Hermes::vector<std::string> areas,
                         HermesFunction* coeff = HERMES_ONE,
                         GeomType gt = HERMES_PLANAR);

    ~DefaultVectorFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} coeff(x, y) * u_ext[0] * v d\bfx
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultResidualVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualVol(int i = 0, std::string area = HERMES_ANY,
                       HermesFunction* coeff = HERMES_ONE,
                       GeomType gt = HERMES_PLANAR);
    DefaultResidualVol(int i, Hermes::vector<std::string> areas,
                       HermesFunction* coeff = HERMES_ONE,
                       GeomType gt = HERMES_PLANAR);

    ~DefaultResidualVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} coeff(u_ext[0]) *
     \nabla u_ext[0] \cdot \nabla v d\bfx
     coeff... (generally nonconstant) function of the solution
  */

  class HERMES_API DefaultResidualDiffusion : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualDiffusion(int i = 0, std::string area = HERMES_ANY,
                             HermesFunction* coeff = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

    DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas,
                             HermesFunction* coeff = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

    ~DefaultResidualDiffusion();

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default volumetric vector form \int_{area} coeff1(u_ext[0]) * u->dx * v->val
     + coeff2(u_ext[0]) * u->dy * v->val d\bfx
     coeff1, coeff2... (generally nonconstant) functions of the solution
  */

  class HERMES_API DefaultResidualAdvection : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualAdvection(int i = 0, std::string area = HERMES_ANY, 
                             HermesFunction* coeff1 = HERMES_ONE,
                             HermesFunction* coeff2 = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);
    DefaultResidualAdvection(int i, Hermes::vector<std::string> areas,
                             HermesFunction* coeff1 = HERMES_ONE,
                             HermesFunction* coeff2 = HERMES_ONE,
                             GeomType gt = HERMES_PLANAR);

    ~DefaultResidualAdvection();

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    private:
      int idx_i;
      HermesFunction* coeff1, *coeff2;
      GeomType gt;
  };

  /* Default surface matrix form \int_{area} coeff(x, y) * u * v dS
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultMatrixFormSurf(int i = 0, int j = 0, std::string area = HERMES_ANY,
                          HermesFunction* coeff = HERMES_ONE,
                          GeomType gt = HERMES_PLANAR);

    DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                          HermesFunction* coeff = HERMES_ONE,
                          GeomType gt = HERMES_PLANAR);

    ~DefaultMatrixFormSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;
      
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

    private:
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default surface matrix form \int_{area} coeff'(u_ext[0]) * u_ext[0] * u * v
     + coeff(u_ext[0]) * u * v dS
     coeff... (generally nonconstant) function of the solution
  */

  class HERMES_API DefaultJacobianFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultJacobianFormSurf(int i = 0, int j = 0, std::string area = HERMES_ANY,
                            HermesFunction* coeff = HERMES_ONE,
                            GeomType gt = HERMES_PLANAR);
    DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas,
                            HermesFunction* coeff = HERMES_ONE,
                            GeomType gt = HERMES_PLANAR);

    ~DefaultJacobianFormSurf();
      
    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

    private:
      int idx_j;
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default surface vector form \int_{area} coeff(x, y) * v dS
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultVectorFormSurf(int i = 0, std::string area = HERMES_ANY,
                          HermesFunction* coeff = HERMES_ONE,
                          GeomType gt = HERMES_PLANAR);
    DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas,
                          HermesFunction* coeff = HERMES_ONE,
                          GeomType gt = HERMES_PLANAR);

    ~DefaultVectorFormSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

    private:
      HermesFunction* coeff;
      GeomType gt;
  };

  class HERMES_API DefaultMultiComponentVectorFormSurf : public WeakForm::MultiComponentVectorFormSurf
  {
  public:
    DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                        std::string area = HERMES_ANY,
                                        Hermes::vector<scalar> coeffs = Hermes::vector<scalar>(1.0),
                                        GeomType gt = HERMES_PLANAR);
    DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                        Hermes::vector<std::string> areas,
                                        Hermes::vector<scalar> coeffs, GeomType gt = HERMES_PLANAR);

    virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    /* FIXME
    virtual WeakForm::VectorFormSurf* clone() {
      return new DefaultMultiComponentVectorFormSurf(*this);
    }
    */

    private:
      Hermes::vector<scalar> coeffs;
      GeomType gt;
  };

  /* Default surface vector form \int_{area} coeff(x, y) * u_ext[0] v dS
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultResidualSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultResidualSurf(int i = 0, std::string area = HERMES_ANY,
                        HermesFunction* coeff = HERMES_ONE,
                        GeomType gt = HERMES_PLANAR);
    DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
                        HermesFunction* coeff = HERMES_ONE,
                        GeomType gt = HERMES_PLANAR);

    ~DefaultResidualSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

    private:
      int idx_i;
      HermesFunction* coeff;
      GeomType gt;
  };

  /* Default weak form for the Laplace equation -div(coeff(u) grad u) = 0. */

  class HERMES_API DefaultWeakFormLaplace : public WeakForm
  {
  public:
    DefaultWeakFormLaplace(std::string area = HERMES_ANY, 
                           HermesFunction* coeff = HERMES_ONE,
                           GeomType gt = HERMES_PLANAR);
  };


  /* Default weak form for the Poisson equation -div(coeff(u) grad u) + f(x, y) = 0. */

  class HERMES_API DefaultWeakFormPoisson : public WeakForm
  {
  public:
    // For not having to indirectly call the WeakForm constructor.
    DefaultWeakFormPoisson();
    DefaultWeakFormPoisson(std::string area, 
                           HermesFunction* coeff = HERMES_ONE,
                           HermesFunction* f = HERMES_ONE,
                           GeomType gt = HERMES_PLANAR);
  };
};

#endif
