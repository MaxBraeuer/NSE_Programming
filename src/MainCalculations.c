#include "Globals.h"
#include "MainCalculations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// linear derivative of the pressure with respect to x
static double dp_dx(SimulationGrid **cell, int i, int j, double delta_x);

// linear derivative of the pressure with respect to y
static double dp_dy(SimulationGrid **cell, int i, int j, double delta_y);

// second linear derivative of u (velocity x-component) with respect to x, calculated by standard central difference
static double d2u_dx2(SimulationGrid **cell, int i, int j, double delta_x);

// second linear derivative of u (velocity x-component) with respect to y, calculated by standard central difference
static double d2u_dy2(SimulationGrid **cell, int i, int j, double delta_y);

// second linear derivative of v (velocity y-component) with respect to x, calculated by standard central difference
static double d2v_dx2(SimulationGrid **cell, int i, int j, double delta_x);

// second linear derivative of v (velocity y-component) with respect to y, calculated by standard central difference
static double d2v_dy2(SimulationGrid **cell, int i, int j, double delta_y);

// nonlinear derivative of u^2 (velocity x-component) with respect to x, calculated by a linear combination of central differences and so called donor-cell stencils
static double du2_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight);

// nonlinear derivative of u*v with respect to y, calculated by a linear combination of central differences and so called donor-cell stencils
static double duv_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight);

// nonlinear derivative of v^2 (velocity y-component) with respect to y, calculated by a linear combination of central differences and so called donor-cell stencils
static double dv2_dy(SimulationGrid **cell, int i, int j, double delta_y, double gamma_weight);

// nonlinear derivative of u*v with respect to x, calculated by a linear combination of central differences and so called donor-cell stencils
static double duv_dx(SimulationGrid **cell, int i, int j, double delta_x, double gamma_weight);

// calculates the pressure according to eq. (42)
static void calculate_pressures(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y, double omega_relax);

// Calculates the residual according to eq. (43)
static void calculate_residuals(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y);

// calculates the discrete L2-norm of the residuals or the pressures (mode 'r' for the residual, mode 'p' for pressure)
static double L2_norm(SimulationGrid **cell, int i_max, int j_max, char mode);


void calculate_F_G(SimulationGrid **cell, int i_max, int j_max, double delta_x,
                  double delta_y, double delta_time, double reynold, double gamma_weight, Gravity g_accel) {

  // F:
  for (int i = 1; i <= i_max - 1; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].F_term = cell[i][j].vel_u + delta_time *
                          ((1 / reynold) * (d2u_dx2(cell, i, j, delta_x) + d2u_dy2(cell, i, j, delta_y)) -
                          du2_dx(cell, i, j, delta_x, gamma_weight) - duv_dy(cell, i, j, delta_y, gamma_weight) + g_accel.g_x);
	}
  }

  // G:
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
    // fill the boundary cells by copying the old pressure values from the neighbouring cells
    for (i = 1; i <= i_max; i++) {
      cell[i][0].pressure = cell[i][1].pressure;
      cell[i][j_max + 1].pressure = cell[i][j_max].pressure;
    }
    for (j = 1; j <= j_max; j++) {
      cell[0][j].pressure = cell[1][j].pressure;
      cell[i_max + 1][j].pressure = cell[i_max][j].pressure;
    }

    // calculate the pressure and residuals for all the other cells
    calculate_pressures(cell, i_max, j_max, delta_x, delta_y, omega_relax);
    calculate_residuals(cell, i_max, j_max, delta_x, delta_y);

    // stop the iteration if the norm of the residual falls below a specified error tolerance
    if (L2_norm(cell, i_max, j_max, 'r') < epsilon_tolerance * pressure_L2_norm_before_sor_step) {
      printf("SOR converged after %d iterations.\n", count + 1);
      return;
    }

    count++;
  }

  // if the maximum number of iterations gets exceeded throw a warning
  printf("SOR did not converge! M%f >= %f\n",L2_norm(cell, i_max, j_max, 'r'),epsilon_tolerance * pressure_L2_norm_before_sor_step);
}


void calculate_velocities(SimulationGrid **cell, int i_max, int j_max, double delta_time, double delta_x, double delta_y) {

  // vel_u
  for (int i = 1; i <= i_max - 1; i++) {
    for (int j = 1; j <= j_max; j++) {
      cell[i][j].vel_u = cell[i][j].F_term - delta_time * dp_dx(cell, i, j, delta_x);

	  // if the velocity component is nan or +/-inf raise an error
	  if (!isfinite(cell[i][j].vel_u)) {
        fprintf(stderr, "Error while updating the velocites. Something went wrong. \n");
        exit(EXIT_FAILURE);
	  }
	}
  }

  // vel_v
  for (int i = 1; i <= i_max; i++) {
	for (int j = 1; j <= j_max - 1; j++) {
      cell[i][j].vel_v = cell[i][j].G_term - delta_time * dp_dy(cell, i, j, delta_y);

	  // if the velocity component is nan or +/-inf raise an error
	  if (!isfinite(cell[i][j].vel_v)) {
        fprintf(stderr, "Error while updating the velocites. Something went wrong. \n");
        exit(EXIT_FAILURE);
	  }
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
