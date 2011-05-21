/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double K_squared) : WeakForm(1) {
    // Jacobian.
    add_matrix_form(new CustomMatrixFormVol(0, 0, K_squared));
    // Residual.
    add_vector_form(new CustomVectorFormVol(0, K_squared));
  };

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j, double K_squared)
          : WeakForm::MatrixFormVol(i, j), K_squared(K_squared) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        val = val + wt[i] * K_squared * u->val[i] * v->val[i];
      }

      return val;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double K_squared;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, double K_squared)
          : WeakForm::VectorFormVol(i), K_squared(K_squared) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        val = val + wt[i] * K_squared * u_ext[0]->val[i] * v->val[i];
        val -= wt[i] * K_squared * v->val[i]; 
      }

      return val;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        val = val + wt[i] * K_squared * u_ext[0]->val[i] * v->val[i];
        val -= wt[i] * K_squared * v->val[i]; 
      }

      return val;
    }

    double K_squared;
  };
};

