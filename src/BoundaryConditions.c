#include "Globals.h"
#include "BoundaryConditions.h"
#include <math.h>
#include <stdio.h>
// No-Slip condition (essentially inflow condition with both velocity components set to zero)
void no_slip_condition(SimulationGrid **cell, int i_max, int j_max, char boundary) {
  inflow_condition(cell, i_max, j_max, 0, 0, boundary);
}

// Inflow condition: velocities are set to a given value / averaged; r = right, l = left, b = bottom, t = top
void inflow_condition(SimulationGrid **cell, int i_max, int j_max, double inflow_vel_x, double inflow_vel_y, char boundary) {
  switch(boundary) {
    case 'r':
      for (int j = 1; j <= j_max; j++) {
        cell[i_max][j].vel_u = inflow_vel_x;
        cell[i_max + 1][j].vel_v = 2 * inflow_vel_y - cell[i_max][j].vel_v;
      }
      break;
    case 'l':
      for (int j = 1; j <= j_max; j++) {
        cell[0][j].vel_u = inflow_vel_x;
        cell[0][j].vel_v = 2 * inflow_vel_y - cell[1][j].vel_v;
      }
      break;
    case 'b':
      for (int i = 1; i <= i_max; i++) {
        cell[i][0].vel_u = 2 * inflow_vel_x - cell[i][1].vel_u;
        cell[i][0].vel_v = inflow_vel_y;
      }
      break;
    case 't':
      for (int i = 1; i <= i_max; i++) {
        cell[i][j_max + 1].vel_u = 2 * inflow_vel_x - cell[i][j_max].vel_u;
        cell[i][j_max].vel_v = inflow_vel_y;
      }
      printf("Nice");
      for (int i = 1; i <= i_max; i++){
        // Choice one
        // cell[i][j_max].pressure = 1.0;
        // cell[i][1].pressure = 1.0;
        // Choice two
        cell[i][j_max].pressure=sin(M_PI*i/i_max)*sinh(M_PI);
      }
      break;
  }
}
