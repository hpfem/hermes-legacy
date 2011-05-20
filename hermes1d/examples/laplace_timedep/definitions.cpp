// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution (all solution components)
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += weights[i] * (u[i] * v[i] / TAU + dudx[i] * dvdx[i]);
  }
  return val;
};

// (nonlinear) form for the residual vector
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...approximate solution
// v...test function
// u_prev...previous solution (all solution components)
double residual(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  // FIXME.
  double u_prev_time[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM];
  error("Registration of external functions in weak forms not implemented yet.");

  double val = 0;
  int comp = 0;    // Solution component.
  int si = 0;      // Solution index (only 0 is relevant for this example).
  for(int i = 0; i<num; i++) {
    val += weights[i] * (u_prev[si][comp][i] - u_prev_time[si][comp][i]) * v[i] / TAU 
                         + du_prevdx[si][comp][i] * dvdx[i] - f(x[i] * v[i]);
  }

  return val;
};
