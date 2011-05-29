#include "definitions.h"

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string omega_1, std::string omega_2, 
                                             std::string omega_3, std::string omega_4, 
                                             std::string omega_5, std::string bdy_left, 
                                             std::string bdy_top, std::string bdy_right, 
                                             std::string bdy_bottom) : WeakForm(1),
  
  omega_1(omega_1), omega_2(omega_2), omega_3(omega_3), 
  omega_4(omega_4), omega_5(omega_5), 

  p_1(25.0),
  p_2(7.0),
  p_3(5.0),
  p_4(0.2),
  p_5(0.05),

  q_1(25.0),
  q_2(0.8),
  q_3(0.0001),
  q_4(0.2),
  q_5(0.05),

  f_1(0.0),
  f_2(1.0),
  f_3(1.0),
  f_4(0.0),
  f_5(0.0),

  bdy_left(bdy_left), 
  bdy_top(bdy_top), 
  bdy_right(bdy_right), 
  bdy_bottom(bdy_bottom),

  c_left(0.0),
  c_top(1.0),
  c_right(2.0),
  c_bottom(3.0),

  g_n_left(0.0),
  g_n_top(3.0),
  g_n_right(2.0),
  g_n_bottom(1.0)

{
    add_matrix_form(new CustomMatrixFormVol(0, 0));
    add_vector_form(new CustomVectorFormVol(0));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_bottom));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_right));
    add_matrix_form_surf(new CustomMatrixFormSurf(0, 0, bdy_top));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_bottom));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_top));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_left));
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_right));
}
