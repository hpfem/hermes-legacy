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
  /* Default volumetric matrix form \int_{area} const_coeff * function_coeff(x, y) * E \cdot F d\bfx
     const_coeff... constant number
     function_coeff... (generally nonconstant) function of x, y
  */

  class HERMES_API DefaultMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    DefaultMatrixFormVol(int i, int j, std::string area = HERMES_ANY,
                         scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);
    DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
                         scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                         SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;

  };

  /* Default volumetric matrix form \int_{area} const_coeff \curl E \curl F d\bfx
     coeff... constant number
  */

  class HERMES_API DefaultJacobianCurlCurl : public WeakForm::MatrixFormVol
  {
  public:
    DefaultJacobianCurlCurl(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                            CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                            SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);
    DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                            CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                            SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

  private:
      scalar coeff;
  };

  /* Default surface matrix form \int_{area} coeff e tau f tau dS
      coeff... constant number
  */

  class HERMES_API DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    DefaultMatrixFormSurf(int i, int j, std::string area = HERMES_ANY,
                          scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);
    DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                          scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR);

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormSurf* clone();

  private:
      scalar coeff;
  };

  /* Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx
     coeff0, coeff1... constant numbers
  */

  class HERMES_API DefaultVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    DefaultVectorFormVol(int i, scalar coeff0, scalar coeff1);
    DefaultVectorFormVol(int i, std::string area, scalar coeff0, scalar coeff1);

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

  private:
    scalar coeff0, coeff1;
  };
};

#endif
