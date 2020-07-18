#ifndef GLOBALS_H
#define GLOBALS_H

typedef struct SimulationGrid {
  double vel_u;         // x-component of the velocity
  double vel_v;         // y-component of the velocity
  double pressure;      // pressure
  double F_term;        // F-term eq. (23)
  double G_term;        // G-term eq. (24)
  double poisson_rhs;   // RHS of the poisson eq. (40)
  double residual;      // residual eq. (43)
} SimulationGrid;

typedef struct Gravity {
  double g_x;           // gravity x-component
  double g_y;           // gravity y-component
} Gravity;

#endif /* GLOABALS_H */ 