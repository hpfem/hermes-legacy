#include "../hermes2d.h"

namespace WeakFormsH1 
{
  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, std::string area, HermesFunction* coeff, SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, Hermes::vector<std::string> areas,
    HermesFunction* coeff, SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultMatrixFormVol::~DefaultMatrixFormVol() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  Ord DefaultMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::MatrixFormVol* DefaultMatrixFormVol::clone() 
  {
    return new DefaultMatrixFormVol(*this);
  }


  DefaultJacobianDiffusion::DefaultJacobianDiffusion(int i, int j, std::string area,
                                                     HermesFunction* coeff,
                                                     SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), idx_j(j), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  };

  DefaultJacobianDiffusion::DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas, 
                                                     HermesFunction* coeff, SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym), idx_j(j), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultJacobianDiffusion::~DefaultJacobianDiffusion() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultJacobianDiffusion::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                  (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                  + coeff->value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + coeff->value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + coeff->value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
    }

    return result;
  }

  Ord DefaultJacobianDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                  (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                  + coeff->value(u_ext[idx_j]->val[i])
                                  * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
                result += wt[i] * e->y[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                          (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                          + coeff->value(u_ext[idx_j]->val[i])
              * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + coeff->value(u_ext[idx_j]->val[i])
                      * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
    }

    return result;
  }

  WeakForm::MatrixFormVol* DefaultJacobianDiffusion::clone() 
  {
    return new DefaultJacobianDiffusion(*this);
  }
  

  DefaultJacobianAdvection::DefaultJacobianAdvection(int i, int j, std::string area, 
                                                     HermesFunction* coeff1,
                                                     HermesFunction* coeff2,
                                                     GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM),
      idx_j(j), coeff1(coeff1), coeff2(coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    if (coeff2 == HERMES_ONE) this->coeff2 = new HermesFunction(1.0);
  }

  DefaultJacobianAdvection::DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas, 
                                                     HermesFunction* coeff1,
                                                     HermesFunction* coeff2,
                                                     GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, HERMES_NONSYM),
      idx_j(j), coeff1(coeff1), coeff2(coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    if (coeff2 == HERMES_ONE) this->coeff2 = new HermesFunction(1.0);
  }

  DefaultJacobianAdvection::~DefaultJacobianAdvection() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff1 != HERMES_ONE) delete coeff1;
    //if (coeff2 != HERMES_ONE) delete coeff2;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianAdvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                               Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (  coeff1->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dx[i] * v->val[i]
                             + coeff1->value(u_ext[idx_j]->val[i]) * u->dx[i] * v->val[i]
                             + coeff2->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dy[i] * v->val[i]
                             + coeff2->value(u_ext[idx_j]->val[i]) * u->dy[i] * v->val[i]);
        }
        return result;
      }

  scalar DefaultJacobianAdvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianAdvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::MatrixFormVol* DefaultJacobianAdvection::clone() 
  {
    return new DefaultJacobianAdvection(*this);
  }


  DefaultVectorFormVol::DefaultVectorFormVol(int i, std::string area,
                                             HermesFunction* coeff,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultVectorFormVol::DefaultVectorFormVol(int i, Hermes::vector<std::string> areas,
                                             HermesFunction* coeff,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultVectorFormVol::~DefaultVectorFormVol() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }
    return result;
  }

  Ord DefaultVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormVol* DefaultVectorFormVol::clone() 
  {
    return new DefaultVectorFormVol(*this);
  }


  DefaultResidualVol::DefaultResidualVol(int i, std::string area,
                                         HermesFunction* coeff,
                                         GeomType gt)
    : WeakForm::VectorFormVol(i, area), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultResidualVol::DefaultResidualVol(int i, Hermes::vector<std::string> areas,
                                         HermesFunction* coeff,
                                         GeomType gt)
    : WeakForm::VectorFormVol(i, areas), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultResidualVol::~DefaultResidualVol() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultResidualVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                   Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }
    return result;
  }

  Ord DefaultResidualVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                              Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualVol::clone() 
  {
    return new DefaultResidualVol(*this);
  }

    
  DefaultResidualDiffusion::DefaultResidualDiffusion(int i, std::string area,
                                                     HermesFunction* coeff, GeomType gt)
    : WeakForm::VectorFormVol(i, area), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  };

  DefaultResidualDiffusion::DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas,
                                                     HermesFunction* coeff, GeomType gt)
    : WeakForm::VectorFormVol(i, areas), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultResidualDiffusion::~DefaultResidualDiffusion() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultResidualDiffusion::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                         Geom<double> *e, ExtData<scalar> *ext) const
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(u_ext[idx_i]->val[i])
                        * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(u_ext[idx_i]->val[i])
                          * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(u_ext[idx_i]->val[i])
                          * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
      }
    }

    return result;
  }

  Ord DefaultResidualDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * coeff->value(u_ext[idx_i]->val[i])
                      * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
    }
    if (gt != HERMES_PLANAR) result = result * Ord(1);

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualDiffusion::clone() 
  {
    return new DefaultResidualDiffusion(*this);
  }
 
  
  DefaultResidualAdvection::DefaultResidualAdvection(int i, std::string area, 
                                                     HermesFunction* coeff1,
                                                     HermesFunction* coeff2,
                                                     GeomType gt)
    : WeakForm::VectorFormVol(i, area), idx_i(i), coeff1(coeff1), coeff2(coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    if (coeff2 == HERMES_ONE) this->coeff2 = new HermesFunction(1.0);
  }
  
  DefaultResidualAdvection::DefaultResidualAdvection(int i, Hermes::vector<std::string> areas,\
                                                     HermesFunction* coeff1,
                                                     HermesFunction* coeff2,
                                                     GeomType gt)
    : WeakForm::VectorFormVol(i, areas),
      idx_i(i), coeff1(coeff1), coeff2(coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    if (coeff2 == HERMES_ONE) this->coeff2 = new HermesFunction(1.0);
  }

  DefaultResidualAdvection::~DefaultResidualAdvection() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff1 != HERMES_ONE) delete coeff1;
    //if (coeff2 != HERMES_ONE) delete coeff2;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultResidualAdvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                               Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    Func<Scalar>* u_prev = u_ext[idx_i];
    for (int i = 0; i < n; i++) {
      result += wt[i] * (coeff1->value(u_prev->val[i]) * (u_prev->dx[i] * v->val[i])
                          + coeff2->value(u_prev->val[i]) * (u_prev->dy[i] * v->val[i]));
    }
    return result;
  }

  scalar DefaultResidualAdvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                         Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualAdvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

  WeakForm::VectorFormVol* DefaultResidualAdvection::clone() 
  {
    return new DefaultResidualAdvection(*this);
  }
 
  
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, std::string area,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }
  
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultMatrixFormSurf::~DefaultMatrixFormSurf() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultMatrixFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  Ord DefaultMatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                 Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::MatrixFormSurf* DefaultMatrixFormSurf::clone() 
  {
    return new DefaultMatrixFormSurf(*this);
  }
 
  
  DefaultJacobianFormSurf::DefaultJacobianFormSurf(int i, int j, std::string area,
                                                   HermesFunction* coeff,
                                                   GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, area),
      idx_j(j), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }
  
  DefaultJacobianFormSurf::DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas,
                                                   HermesFunction* coeff,
                                                   GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultJacobianFormSurf::~DefaultJacobianFormSurf() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianFormSurf::matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                   Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u_ext[idx_j]->val[i]
                          + coeff->value(u_ext[idx_j]->val[i]))
                * u->val[i] * v->val[i];
    }
    return result;
  }

  scalar DefaultJacobianFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                        Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form_surf<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                   Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  WeakForm::MatrixFormSurf* DefaultJacobianFormSurf::clone() 
  {
    return new DefaultJacobianFormSurf(*this);
  }


  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, std::string area,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }
  
  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultVectorFormSurf::~DefaultVectorFormSurf() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return result;
  }

  Ord DefaultVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormSurf* DefaultVectorFormSurf::clone() 
  {
    return new DefaultVectorFormSurf(*this);
  }
 
  
  DefaultMultiComponentVectorFormSurf::DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                                                           std::string area,
                                                                           Hermes::vector<scalar> coeffs,
                                                                           GeomType gt)
  : WeakForm::MultiComponentVectorFormSurf(coordinates, area), coeffs(coeffs), gt(gt) 
  { 
  }
  
  DefaultMultiComponentVectorFormSurf::DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                                                           Hermes::vector<std::string> areas,
                                                                           Hermes::vector<scalar> coeffs, GeomType gt)
  : WeakForm::MultiComponentVectorFormSurf(coordinates, areas), coeffs(coeffs), gt(gt) 
  { 
  }

  void DefaultMultiComponentVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                                  Geom<double> *e, ExtData<scalar> *ext, 
                                                  Hermes::vector<scalar>& result) const 
  {
    scalar result_base = 0;
    if (gt == HERMES_PLANAR)
      result_base = int_v<double>(n, wt, v);
    else
      if (gt == HERMES_AXISYM_X)
        result_base = int_y_v<double>(n, wt, v, e);
      else
        result_base = int_x_v<double>(n, wt, v, e);

    for(unsigned int result_i = 0; result_i < this->coordinates.size(); result_i++)
      result.push_back(result_base * coeffs[result_i]);
  }

  Ord DefaultMultiComponentVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                               Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    if (gt == HERMES_PLANAR)
      return int_v<Ord>(n, wt, v);
    else
      if (gt == HERMES_AXISYM_X)
        return int_y_v<Ord>(n, wt, v, e);
      else
        return int_x_v<Ord>(n, wt, v, e);
  }

  /*
  WeakForm::VectorFormSurf* DefaultMultiComponentVectorFormSurf::clone() {
    return new DefaultMultiComponentVectorFormSurf(*this);
  }
  */

  DefaultResidualSurf::DefaultResidualSurf(int i, std::string area,
                                           HermesFunction* coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, area), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }
  
  DefaultResidualSurf::DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
                                           HermesFunction* coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, areas), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
  }

  DefaultResidualSurf::~DefaultResidualSurf() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultResidualSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  Ord DefaultResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                               Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormSurf* DefaultResidualSurf::clone() 
  {
    return new DefaultResidualSurf(*this);
  }

  
  DefaultWeakFormLaplace::DefaultWeakFormLaplace(std::string area, 
                                                 HermesFunction* coeff,
                                                 GeomType gt) : WeakForm()
  {
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, area, coeff, HERMES_NONSYM, gt));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, area, coeff, gt));
  };
  
  DefaultWeakFormPoisson::DefaultWeakFormPoisson() : WeakForm()
  {
  }
  
  DefaultWeakFormPoisson::DefaultWeakFormPoisson(std::string area,
                                                 HermesFunction* coeff,
                                                 HermesFunction* f,
                                                 GeomType gt) : WeakForm()
  {
    // Jacobian.
    // NOTE: The flag HERMES_NONSYM is important here.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, area, coeff, HERMES_NONSYM, gt));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, area, coeff, gt));
    add_vector_form(new DefaultVectorFormVol(0, area, f, gt));
  };
};
