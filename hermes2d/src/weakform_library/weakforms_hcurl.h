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

#ifndef __H2D_HCURL_WEAK_FORMS_H
#define __H2D_HCURL_WEAK_FORMS_H

#include "../integrals/hcurl.h"

namespace WeakFormsHcurl 
{
  /* Default volumetric matrix form \int_{area} coeff(x, y) * E \cdot F d\bfx
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    DefaultMatrixFormVol(int i, int j, std::string area = HERMES_ANY,
                         HermesFunction* coeff = HERMES_ONE,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
                         HermesFunction* coeff = HERMES_ONE,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultMatrixFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
        HermesFunction* coeff;
        GeomType gt;
  };

  /* FIXME
     Default volumetric matrix form \int_{area} coeff(x, y) \curl E \curl F d\bfx
     coeff... so far only constant number possible
  */

  class HERMES_API DefaultJacobianCurlCurl : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianCurlCurl(int i, int j, std::string area = HERMES_ANY,
                            HermesFunction* coeff = HERMES_ONE,
                            SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas,
                            HermesFunction* coeff = HERMES_ONE,
                            SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    ~DefaultJacobianCurlCurl();

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

  /* FIXME
     Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx
     coeff0, coeff1... constant numbers
  */

  class HERMES_API DefaultVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultVectorFormVol(int i, std::string area = HERMES_ANY, 
                         HermesFunction* coeff0 = HERMES_ONE,
                         HermesFunction* coeff1 = HERMES_ONE,
                         GeomType gt = HERMES_PLANAR);

    DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, 
                         HermesFunction* coeff0 = HERMES_ONE,
                         HermesFunction* coeff1 = HERMES_ONE,
                         GeomType gt = HERMES_PLANAR);

    ~DefaultVectorFormVol();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

  private:
    HermesFunction* coeff0, *coeff1;
    GeomType gt;
  };

  /* Default volumetric vector form \int_{area} coeff(x, y) * u_ext[0] * v d\bfx
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultResidualVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualVol(int i, std::string area = HERMES_ANY,
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

  /*
     Default volumetric vector form \int_{area} coeff(u_ext[0]) *
     \nabla u_ext[0] \cdot \nabla v d\bfx
     coeff... (generally nonconstant) function of the solution
  */

  class HERMES_API DefaultResidualCurlCurl : public WeakForm::VectorFormVol
  {
  public:
    DefaultResidualCurlCurl(int i, std::string area = HERMES_ANY,
                            HermesFunction* coeff = HERMES_ONE,
                            GeomType gt = HERMES_PLANAR);

    DefaultResidualCurlCurl(int i, Hermes::vector<std::string> areas,
                            HermesFunction* coeff = HERMES_ONE,
                            GeomType gt = HERMES_PLANAR);

    ~DefaultResidualCurlCurl();

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

  /* FIXME
     Default surface matrix form \int_{area} coeff e tau f tau dS
     coeff... constant number
  */

  class HERMES_API DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultMatrixFormSurf(int i, int j, std::string area = HERMES_ANY,
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

  /* FIXME
     Default surface vector form \int_{area} coeff(x, y) * v dS
     coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultVectorFormSurf(int i, std::string area = HERMES_ANY,
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

  /* FIXME
     Default surface residual form \int_{area} coeff u_ext[0] tau f tau dS
     coeff... constant number
  */

  class HERMES_API DefaultResidualSurf : public WeakForm::VectorFormSurf
  {
  public:
    DefaultResidualSurf(int i, std::string area = HERMES_ANY,
                        HermesFunction* coeff = HERMES_ONE,
                        GeomType gt = HERMES_PLANAR);

    DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
                        HermesFunction* coeff = HERMES_ONE,
                        GeomType gt = HERMES_PLANAR);

    ~DefaultResidualSurf();

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

  private:
      HermesFunction* coeff;
      GeomType gt;
  };
}

#endif
