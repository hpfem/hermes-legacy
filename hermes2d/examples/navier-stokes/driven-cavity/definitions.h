#include "hermes2d.h"

using namespace WeakFormsH1;

class WeakFormNSNewton : public WeakForm
{
public:
  WeakFormNSNewton(bool Stokes, double Reynolds, double time_step, Solution* x_vel_previous_time, 
                   Solution* y_vel_previous_time);

  class BilinearFormSymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), 
                        Reynolds(Reynolds), time_step(time_step) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormNonsymVel_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j, bool Stokes) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const; 

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j, bool Stokes) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j, bool Stokes) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j, bool Stokes) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
  };

  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

  };

  class BilinearFormNonsymYVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
  };

  class VectorFormNS_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step) : WeakForm::VectorFormVol(i), 
                   Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step) 
          : WeakForm::VectorFormVol(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_2 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_2(int i) : WeakForm::VectorFormVol(i) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution* x_vel_previous_time;
  Solution* y_vel_previous_time;
};

/* Essential boundary conditions */

// Time-dependent surface x-velocity of inner circle.
class EssentialBCNonConstX : public EssentialBoundaryCondition
{
public:
  EssentialBCNonConstX(Hermes::vector<std::string> markers, double vel, double startup_time) 
             : EssentialBoundaryCondition(markers), startup_time(startup_time), vel(vel)  {};

  EssentialBCNonConstX(std::string marker, double vel, double startup_time) 
             : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)), startup_time(startup_time), vel(vel)  {};
  
  ~EssentialBCNonConstX() {};

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

protected:
  double startup_time;
  double vel;
};

// Time-dependent surface y-velocity of inner circle.
class EssentialBCNonConstY : public EssentialBoundaryCondition
{
public:
  EssentialBCNonConstY(Hermes::vector<std::string> markers, double vel, double startup_time) 
             : EssentialBoundaryCondition(markers), startup_time(startup_time), vel(vel)  {};

  EssentialBCNonConstY(std::string marker, double vel, double startup_time) 
             : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)), startup_time(startup_time), vel(vel)  {};
  
  ~EssentialBCNonConstY() {};

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

protected:
  double startup_time;
  double vel;
};

