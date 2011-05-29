#include "../hermes2d.h"

namespace WeakFormsHcurl 
{
  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, std::string area,
    HermesFunction* coeff, SymFlag sym, 
    GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }
 
  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, Hermes::vector<std::string> areas,
    HermesFunction* coeff, SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultMatrixFormVol::~DefaultMatrixFormVol() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    if (coeff->is_constant() != true) error("Nonconstant coeff in curl curl Jacobian not implemented yet.");
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) 
                        * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    if (coeff->is_constant() != true) error("Nonconstant coeff in curl curl Jacobian not implemented yet.");
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) 
                        * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  WeakForm::MatrixFormVol* DefaultMatrixFormVol::clone() 
  {
    return new DefaultMatrixFormVol(*this);
  }
     

  DefaultJacobianCurlCurl::DefaultJacobianCurlCurl(int i, int j, std::string area,
                                                   HermesFunction* coeff,
                                                   SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), idx_j(j), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  };

  DefaultJacobianCurlCurl::DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas, 
                                                   HermesFunction* coeff,
                                                   SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym),
      idx_j(j), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultJacobianCurlCurl::~DefaultJacobianCurlCurl() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultJacobianCurlCurl::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                        Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    if (coeff->is_constant() != true) error("Nonconstant coeff in curl curl Jacobian not implemented yet.");
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * (u->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultJacobianCurlCurl::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    if (coeff->is_constant() != true) error("Nonconstant coeff in curl curl Jacobian not implemented yet.");
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * (u->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  WeakForm::MatrixFormVol* DefaultJacobianCurlCurl::clone() 
  {
    return new DefaultJacobianCurlCurl(*this);
  }
     
      
  DefaultVectorFormVol::DefaultVectorFormVol(int i, std::string area,
                                             HermesFunction* coeff0, HermesFunction* coeff1,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, area), coeff0(coeff0), coeff1(coeff1), gt(gt)
  { 
    // If coeff0 is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff0 == HERMES_ONE) this->coeff0 = new HermesFunction(1.0);
    else if(!coeff0->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
    // If coeff1 is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    else if(!coeff1->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormVol::DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, 
                                             HermesFunction* coeff0, HermesFunction* coeff1,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, areas), coeff0(coeff0), coeff1(coeff1), gt(gt)
  { 
    // If coeff0 is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff0 == HERMES_ONE) this->coeff0 = new HermesFunction(1.0);
    else if(!coeff0->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
    // If coeff1 is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff1 == HERMES_ONE) this->coeff1 = new HermesFunction(1.0);
    else if(!coeff1->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormVol::~DefaultVectorFormVol() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff0 != HERMES_ONE) delete coeff0;
    //if (coeff1 != HERMES_ONE) delete coeff1;
  };

  scalar DefaultVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * coeff0->value(e->x[i], e->y[i]) * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * coeff1->value(e->x[i], e->y[i]) * v->val1[i];
    return int_v0 + int_v1;
  }

  Ord DefaultVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * coeff0->value(e->x[i], e->y[i]) * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * coeff1->value(e->x[i], e->y[i]) * v->val1[i];
    return int_v0 + int_v1;
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
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualVol::DefaultResidualVol(int i, Hermes::vector<std::string> areas,
                                         HermesFunction* coeff,
                                         GeomType gt)
    : WeakForm::VectorFormVol(i, areas), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
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
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * (u_ext[idx_i]->val0[i] * v->val0[i] +
                                                            u_ext[idx_i]->val1[i] * v->val1[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

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
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualVol::clone() 
  {
    return new DefaultResidualVol(*this);
  }


  DefaultResidualCurlCurl::DefaultResidualCurlCurl(int i, std::string area,
                                                   HermesFunction* coeff,
                                                   GeomType gt)
    : WeakForm::VectorFormVol(i, area), idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  };

  DefaultResidualCurlCurl::DefaultResidualCurlCurl(int i, Hermes::vector<std::string> areas,
                                                   HermesFunction* coeff,
                                                   GeomType gt)
    : WeakForm::VectorFormVol(i, areas),
      idx_i(i), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultResidualCurlCurl::~DefaultResidualCurlCurl() 
  {
    // FIXME: Should be deleted here only if it was created here.
    //if (coeff != HERMES_ONE) delete coeff;
  };

  scalar DefaultResidualCurlCurl::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                        Geom<double> *e, ExtData<scalar> *ext) const
  {
    Func<scalar>* u_prev = u_ext[idx_i];
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        double mag0_i = std::abs(u_prev->val0[i]);
        double mag1_i = std::abs(u_prev->val1[i]);
        double mag_i = sqrt(sqr(mag0_i) + sqr(mag1_i));
        result += wt[i] * coeff->value(mag_i) 
	                * (u_prev->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultResidualCurlCurl::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Func<Ord>* u_prev = u_ext[idx_i];
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        Ord mag0_i = u_prev->val0[i];
        Ord mag1_i = u_prev->val1[i];
        Ord mag_i = sqrt(sqr(mag0_i) + sqr(mag1_i));
        result += wt[i] * coeff->value(mag_i) 
	                * (u_prev->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualCurlCurl::clone() 
  {
    return new DefaultResidualCurlCurl(*this);
  }
 

  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, std::string area,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
  }
  
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
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
      for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * 
                          ((u->val0[i] * e->tx[i] + u->val1[i] * e->ty[i]) *
                           conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  Ord DefaultMatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
       for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) * 
                          ((u->val0[i] * e->tx[i] + u->val1[i] * e->ty[i]) *
                           conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));     
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  WeakForm::MatrixFormSurf* DefaultMatrixFormSurf::clone() 
  {
    return new DefaultMatrixFormSurf(*this);
  }


  DefaultResidualSurf::DefaultResidualSurf(int i, std::string area,
                                           HermesFunction* coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualSurf::DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
                                           HermesFunction* coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant functions in Hcurl forms not implemented yet.");
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
      for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) 
                        * (    (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
                           conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  Ord DefaultResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[],
                               Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++)
        result += wt[i] * coeff->value(e->x[i], e->y[i]) 
                        * (    (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
                           conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
     }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  WeakForm::VectorFormSurf* DefaultResidualSurf::clone()
  {
    return new DefaultResidualSurf(*this);
  }


  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, std::string area,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }
  
  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas,
                                               HermesFunction* coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, areas), coeff(coeff), gt(gt)
  {
    // If coeff is HERMES_ONE, initialize it to be constant 1.0.
    if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
    else if(!coeff->is_constant()) error("Nonconstant coefficients in Hcurl forms not implemented yet.");
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
        result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  WeakForm::VectorFormSurf* DefaultVectorFormSurf::clone() 
  {
    return new DefaultVectorFormSurf(*this);
  }
 
    
};
