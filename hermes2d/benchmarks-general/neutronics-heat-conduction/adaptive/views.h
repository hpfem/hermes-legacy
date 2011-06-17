#include "../definitions.h"

class ViewsAdapt : public Views
{
  OrderView *oview_T;
  OrderView *oview_phi;
  
  public:
    ViewsAdapt(unsigned int width, unsigned int height);
    virtual ~ViewsAdapt();

    virtual void show_solutions(double current_time, Hermes::vector< Solution* > solutions) {
      Views::show_solutions(current_time, solutions);
    }
    
    void show_solutions(double current_time, int adaptivity_step, Hermes::vector< Solution* > solutions);
    
    void show_orders(double current_time, Hermes::vector< Space* > spaces);
    void show_orders(double current_time, int adaptivity_step, Hermes::vector< Space* > spaces);
};