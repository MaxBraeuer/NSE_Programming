#include "Globals.h"
#include "MainCalculations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double dp_dx(SimulationGrid **cell, int i, int j, double delta_x);

static double dp_dy(SimulationGrid **cell, int i, int j, double delta_y);

static double d2u_dx2(SimulationGrid **cell, int i, int j, double delta_x);

static double d2u_dy2(SimulationGrid **cell, int i, int j, double delta_y);

static double d2v_dx2(SimulationGrid **cell, int i, int j, double delta_x);

static double d2v_dy2(SimulationGrid **cell, int i, int j, double delta_y);

static double du2_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight);

static double duv_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight);

static double dv2_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight);

static double duv_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight);

static void calculate_pressures(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y, double omega_relax);

static void calculate_residuals(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y);

static double L2_norm(SimulationGrid **cell, int i_max, int j_max, char mode);


void calculate_F_G(SimulationGrid **cell, int i_max, int j_max, double delta_x,
                  double delta_y, double delta_time, double reynold, double gamma_weight, Gravity g_accel) {

  for (int i = 1; i <= i_max - 1; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].F_term = cell[i][j].vel_u + delta_time *
                          ((1 / reynold) * (d2u_dx2(cell, i, j, delta_x) + d2u_dy2(cell, i, j, delta_y)) -
                          du2_dx(cell, i, j, delta_x, gamma_weight) - duv_dy(cell, i, j, delta_y, gamma_weight) + g_accel.g_x);
	}
  }

  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max - 1; j++) {
      cell[i][j].G_term = cell[i][j].vel_v + delta_time *
                          ((1 / reynold) * (d2v_dx2(cell, i, j, delta_x) + d2v_dy2(cell, i, j, delta_y)) -
                          duv_dx(cell, i, j, delta_x, gamma_weight) - dv2_dy(cell, i, j, delta_y, gamma_weight) + g_accel.g_y);
    }
  }
}

void rhs_pressure(SimulationGrid **cell, int i_max, int j_max, double delta_time, double delta_x, double delta_y) {
  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].poisson_rhs = (1.0 / delta_time) * ((cell[i][j].F_term - cell[i-1][j].F_term) / delta_x + (cell[i][j].G_term - cell[i][j-1].G_term) / delta_y);
    }
  }
}

void sor_loop(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y,
             double omega_relax, double epsilon_tolerance, int max_iterations) {
  int count = 0;
  int i, j;
  double pressure_L2_norm_before_sor_step = L2_norm(cell, i_max, j_max, 'p');
  while (count < max_iterations) {
    for (i = 1; i <= i_max; i++) {
      cell[i][0].pressure = cell[i][1].pressure;
      cell[i][j_max + 1].pressure = cell[i][j_max].pressure;
    }
    for (j = 1; j <= j_max; j++) {
      cell[0][j].pressure = cell[1][j].pressure;
      cell[i_max + 1][j].pressure = cell[i_max][j].pressure;
    }

    calculate_pressures(cell, i_max, j_max, delta_x, delta_y, omega_relax);
    calculate_residuals(cell, i_max, j_max, delta_x, delta_y);

    if (L2_norm(cell, i_max, j_max, 'r') < epsilon_tolerance * pressure_L2_norm_before_sor_step) {
      printf("SOR converged after %d iterations.\n", count + 1);
      return;
    }

    count++;
  }
}


void calculate_velocities(SimulationGrid **cell, int i_max, int j_max, double delta_time, double delta_x, double delta_y) {

  for (int i = 1; i <= i_max - 1; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].vel_u = cell[i][j].F_term - delta_time * dp_dx(cell, i, j, delta_x);

	  }
	}
  }

  for (int i = 1; i <= i_max; i++) {
	for (int j = 1; j <= j_max - 1; j++) {
      cell[i][j].vel_v = cell[i][j].G_term - delta_time * dp_dy(cell, i, j, delta_y);

	}
  }
}

// Stencils are verified
static double dp_dx(SimulationGrid **cell, int i, int j, double delta_x) {
  return ((cell[i+1][j].pressure - cell[i][j].pressure) / delta_x);
}

static double dp_dy(SimulationGrid **cell, int i, int j, double delta_y) {
  return ((cell[i][j+1].pressure - cell[i][j].pressure) / delta_y);
}

static double d2u_dx2(SimulationGrid **cell, int i, int j, double delta_x) {
  return ((cell[i+1][j].vel_u - 2 * cell[i][j].vel_u + cell[i-1][j].vel_u) / (delta_x * delta_x));
}

static double d2u_dy2(SimulationGrid **cell, int i, int j, double delta_y) {
  return ((cell[i][j+1].vel_u - 2 * cell[i][j].vel_u + cell[i][j-1].vel_u) / (delta_y * delta_y));
}

static double d2v_dx2(SimulationGrid **cell, int i, int j, double delta_x) {
  return ((cell[i+1][j].vel_v - 2 * cell[i][j].vel_v + cell[i-1][j].vel_v) / (delta_x * delta_x));
}

static double d2v_dy2(SimulationGrid **cell, int i, int j, double delta_y) {
  return ((cell[i][j+1].vel_v - 2 * cell[i][j].vel_v + cell[i][j-1].vel_v) / (delta_y * delta_y));
}

static double du2_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight) {
  return ((1 / (4 * delta_x)) * ((pow(cell[i][j].vel_u + cell[i+1][j].vel_u, 2) - pow(cell[i-1][j].vel_u + cell[i][j].vel_u, 2)) +
    gamma_weight * (fabs(cell[i][j].vel_u + cell[i+1][j].vel_u) * (cell[i][j].vel_u - cell[i+1][j].vel_u) -
                    fabs(cell[i-1][j].vel_u + cell[i][j].vel_u) * (cell[i-1][j].vel_u - cell[i][j].vel_u))));
}

static double duv_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight) {
  return ((1 / (4 * delta_y)) * (((cell[i][j].vel_v + cell[i+1][j].vel_v) * (cell[i][j].vel_u + cell[i][j+1].vel_u) -
                                  (cell[i][j-1].vel_v + cell[i+1][j-1].vel_v) * (cell[i][j-1].vel_u + cell[i][j].vel_u)) +
                                gamma_weight * (fabs(cell[i][j].vel_v + cell[i+1][j].vel_v) * (cell[i][j].vel_u - cell[i][j+1].vel_u) -
                                fabs(cell[i][j-1].vel_v + cell[i+1][j-1].vel_v) * (cell[i][j-1].vel_u - cell[i][j].vel_u))));
}

static double dv2_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight) {
  return ((1 / (4 * delta_y)) * ((pow(cell[i][j].vel_v + cell[i][j+1].vel_v, 2) - pow(cell[i][j-1].vel_v + cell[i][j].vel_v, 2)) +
                                gamma_weight * (fabs(cell[i][j].vel_v + cell[i][j+1].vel_v) * (cell[i][j].vel_v - cell[i][j+1].vel_v)
                                - fabs(cell[i][j-1].vel_v + cell[i][j].vel_v) * (cell[i][j-1].vel_v - cell[i][j].vel_v))));
}

static double duv_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight) {
  return ((1 / (4 * delta_x)) * (((cell[i][j].vel_u + cell[i][j+1].vel_u) * (cell[i][j].vel_v + cell[i+1][j].vel_v) -
                                  (cell[i-1][j].vel_u + cell[i-1][j+1].vel_u) * (cell[i-1][j].vel_v + cell[i][j].vel_v)) +
                                gamma_weight * (fabs(cell[i][j].vel_u + cell[i][j+1].vel_u) * (cell[i][j].vel_v - cell[i+1][j].vel_v) -
                                fabs(cell[i-1][j].vel_u + cell[i-1][j+1].vel_u) * (cell[i-1][j].vel_v - cell[i][j].vel_v))));
}

static void calculate_pressures(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y, double omega_relax) {
  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].pressure = (1 - omega_relax) * cell[i][j].pressure + omega_relax / (2 * ((1 / (delta_x * delta_x)) + (1 / (delta_y * delta_y)))) *
            ((cell[i+1][j].pressure + cell[i-1][j].pressure) / (delta_x * delta_x) +
            (cell[i][j+1].pressure + cell[i][j-1].pressure) / (delta_y * delta_y)  - cell[i][j].poisson_rhs);
    }
  }
}

static void calculate_residuals(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y) {
  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].residual = ((cell[i+1][j].pressure - 2 * cell[i][j].pressure + cell[i-1][j].pressure) / (delta_x * delta_x) +
          (cell[i][j+1].pressure - 2 * cell[i][j].pressure + cell[i][j-1].pressure) / (delta_y * delta_y) - cell[i][j].poisson_rhs);
    }
  }
}

static double L2_norm(SimulationGrid **cell, int i_max, int j_max, char mode) {
  double square_sum = 0.0;
  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      if (mode == 'r') {
          square_sum += cell[i][j].residual * cell[i][j].residual;
      } else if (mode == 'p') {
          square_sum += cell[i][j].pressure * cell[i][j].pressure;
      }
    }
  }
  return (sqrt(square_sum / (i_max * j_max)));
}
