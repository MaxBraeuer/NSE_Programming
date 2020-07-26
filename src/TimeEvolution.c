#include "Globals.h"
#include "TimeEvolution.h"
#include "BoundaryConditions.h"
#include "MainCalculations.h"
#include "InputOutput.h"

#include <stdarg.h>
#include <math.h>
#include <stdio.h>

/* declare static functions */

// Search the maximum (absolute value) velocity component in x-direction
static double search_max_vel_u(SimulationGrid **cell, int i_max, int j_max);
// Search the maximum (absolute value) velocity component in y-direction
static double search_max_vel_v(SimulationGrid **cell, int i_max, int j_max);
// Variadic function that returns the smallest double
static double search_min_double(int n_args, ...);


// Main time evolution function
void time_evolution(int problem_type, double inflow_vel, double frequency, int i_max, int j_max,
                   double reynold, double a_size, double b_size, double tau_safety, double omega_relax,
                   double epsilon_tolerance, int max_iterations, Gravity g_accel,
                   double max_int_time, int max_int_steps, SimulationGrid **cell) {

  // Declare some new needed variables
  double vel_u_max_abs, vel_v_max_abs;  // maximum velocity components (absolute value)
  double gamma_weight;                  // factor for donor-cell-stencil weighing
  double delta_time;                    // size of the timestep
  double int_time = 0.0;                // initialize the integration time
  int time_step_number = 1;				// initialize the time step number
  double delta_x = a_size / i_max;      // cell size in x-direction
  double delta_y = b_size / j_max;      // cell size in y-direction
  double reynold_stability_condition = reynold / (2.0 * ((1.0 / (delta_x * delta_x)) + (1.0 / (delta_y * delta_y))));

  no_slip_condition(cell, i_max, j_max, 'r');
  no_slip_condition(cell, i_max, j_max, 'l');
  no_slip_condition(cell, i_max, j_max, 'b');

  // the top boundary is set using the inflow condition;
  // depending on the problem type the velocity x-component is either constant or varies
  if (problem_type == 0) {
    inflow_condition(cell, i_max, j_max, inflow_vel, 0.0, 't');
  } else if (problem_type == 1) {
    inflow_condition(cell, i_max, j_max, inflow_vel * sin(frequency * int_time), 0.0, 't');
  }
  while (int_time < max_int_time && time_step_number <= max_int_steps) {
    // calculate adaptive size of the timestep to ensure a stable simulation: eq. (26)
    vel_u_max_abs = search_max_vel_u(cell, i_max, j_max);
    vel_v_max_abs = search_max_vel_v(cell, i_max, j_max);

    // Catch a possible division by zero. C doesn't really care but just to be on the safe side
    if (vel_u_max_abs == 0 || vel_v_max_abs == 0) {
      delta_time = tau_safety * reynold_stability_condition;
    } else {
      delta_time = tau_safety * search_min_double(3, reynold_stability_condition, delta_x / vel_u_max_abs, delta_y / vel_v_max_abs);
    }

    // implement boundary conditions: eq. (8-10)
    // the right, left and bottom boundaries are set using the no-slip condition
    export_cells(cell, i_max, j_max, time_step_number-1, int_time);

    // calculate F and G: eq. (29-30)
    gamma_weight = fmax(vel_u_max_abs * delta_time / delta_x, vel_v_max_abs * delta_time / delta_y);

    // calculate_F_G(cell, i_max, j_max, delta_x, delta_y, delta_time, reynold, gamma_weight, g_accel);

    // Calculate the RHS of the pressure eq. (40)
    // rhs_pressure(cell, i_max, j_max, delta_time, delta_x, delta_y);

    // begin SOR-loop
    sor_loop(cell, i_max, j_max, delta_x, delta_y, omega_relax, epsilon_tolerance, max_iterations);

    // calculate new velocity components u_vel and v_vel: eq. (21-22)
    // calculate_velocities(cell, i_max, j_max, delta_time, delta_x, delta_y);
    // printf("Velocities updated!\n");

    // increment the integration time by the size of the adaptive timestep
    int_time += delta_time;
    printf("Integration time elapsed: %lg / %lg \n", int_time, max_int_time);
    printf("Time steps done: %d / %d \n", time_step_number, max_int_steps);

    // export the cells
    export_cells(cell, i_max, j_max, time_step_number, int_time);

    time_step_number += 1;

  }
}

static double search_max_vel_u(SimulationGrid **cell, int i_max, int j_max) {
  double max_val = 0.0;

  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      if (max_val < fabs(cell[i][j].vel_u)) {
        max_val = fabs(cell[i][j].vel_u);
      }
    }
  }
  return max_val;
}

static double search_max_vel_v(SimulationGrid **cell, int i_max, int j_max) {
  double max_val = 0.0;

  for (int i = 1; i <= i_max; i++) {
    for (int j = 1; j <= j_max; j++) {
      if (max_val < fabs(cell[i][j].vel_v)) {
        max_val = fabs(cell[i][j].vel_v);
      }
    }
  }
  return max_val;
}

static double search_min_double(int n_args, ...) {
  va_list num_list;
  va_start(num_list, n_args);

  double temp_val;
  double min_val = va_arg(num_list, double);

  for (int i = 1; i < n_args; i++) {
    temp_val = va_arg(num_list, double);
    if (min_val > temp_val) {
      min_val = temp_val;
    }
  }
  va_end(num_list);

  return min_val;
}
