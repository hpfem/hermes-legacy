#include "hermes2d.h"
#include "weakform/weakform.h"
#include <integrals/h1.h>
#include "boundaryconditions/essential_bcs.h"

class EssentialBCNonConst : public EssentialBoundaryCondition 
{
public:
    EssentialBCNonConst(std::string marker);

    ~EssentialBCNonConst() {};

    virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

    virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

class WeakFormHelmholtz : public WeakForm
{
public:
    // Problems parameters
    WeakFormHelmholtz(double eps, double mu, double omega, double sigma, double beta, double E0, double h);

private:
    class MatrixFormHelmholtzEquation_real_real : public WeakForm::MatrixFormVol
    {
    public:
        MatrixFormHelmholtzEquation_real_real(int i, int j, double eps, double omega, double mu)
              : WeakForm::MatrixFormVol(i, j), eps(eps), omega(omega), mu(mu) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    private:
        // Memebers.
        double eps;
        double omega;
        double mu;
    };

    class MatrixFormHelmholtzEquation_real_imag : public WeakForm::MatrixFormVol
    {
    private:
        // Memebers.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_real_imag(int i, int j, double mu, double omega, double sigma) 
              : WeakForm::MatrixFormVol(i, j), mu(mu), omega(omega), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    };

    class MatrixFormHelmholtzEquation_imag_real : public WeakForm::MatrixFormVol
    {
    private:
        // Memebers.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_imag_real(int i, int j, double mu, double omega, double sigma) 
              : WeakForm::MatrixFormVol(i, j), mu(mu), omega(omega), sigma(sigma) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    };

    class MatrixFormHelmholtzEquation_imag_imag : public WeakForm::MatrixFormVol
    {
    public:
        MatrixFormHelmholtzEquation_imag_imag(int i, int j, double eps, double mu, double omega) 
              : WeakForm::MatrixFormVol(i, j), eps(eps), mu(mu), omega(omega) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

        // Memebers.
        double eps;
        double mu;
        double omega;
    };


    class MatrixFormSurfHelmholtz_real_imag : public WeakForm::MatrixFormSurf
    {
    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_real_imag(int i, int j, std::string area, double beta)
              : WeakForm::MatrixFormSurf(i, j, area), beta(beta){};

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    };



    class MatrixFormSurfHelmholtz_imag_real : public WeakForm::MatrixFormSurf
    {

    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_imag_real(int i, int j, std::string area, double beta)
              : WeakForm::MatrixFormSurf(i, j, area), beta(beta){};

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    };
};
