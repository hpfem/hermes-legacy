#include "views.h"

ViewsAdapt::ViewsAdapt(unsigned int width, unsigned int height) : Views(width, height)
{
  oview_T = new OrderView("", new WinGeom(2*width, 0, width, height));
  oview_T->fix_scale_width(80);
  oview_phi = new OrderView("", new WinGeom(2*width, height+50, width, height));
  oview_phi->fix_scale_width(80);
}

ViewsAdapt::~ViewsAdapt()
{
  delete oview_T;
  delete oview_phi;
}

void ViewsAdapt::show_solutions(double current_time, int adaptivity_step, Hermes::vector< Solution* > solutions)
{
  // Show the new time level solution.
  sprintf(title, "T (coarse mesh), t = %g s, adapt step %d", current_time, adaptivity_step);
  sview_T->set_title(title); 
  sview_T->show(solutions[0]);
  
  sprintf(title, "phi (coarse mesh), t = %g s, adapt step %d", current_time, adaptivity_step);
  sview_phi->set_title(title);
  sview_phi->show(solutions[1]);
}

void ViewsAdapt::show_orders(double current_time, int adaptivity_step, Hermes::vector< Space* > spaces)
{
  sprintf(title, "T mesh (coarse), t = %g, adapt step %d", current_time, adaptivity_step);
  oview_T->set_title(title);
  oview_T->show(spaces[0]);
  
  sprintf(title, "phi mesh (coarse), t = %g, adapt step %d", current_time, adaptivity_step);
  oview_phi->set_title(title);
  oview_phi->show(spaces[1]);
}

void ViewsAdapt::show_orders(double current_time, Hermes::vector< Space* > spaces)
{
  sprintf(title, "T mesh (coarse), t = %g, final mesh", current_time);
  oview_T->set_title(title);
  oview_T->show(spaces[0]);
  
  sprintf(title, "phi mesh (coarse), t = %g, final mesh", current_time);
  oview_phi->set_title(title);
  oview_phi->show(spaces[1]);
}