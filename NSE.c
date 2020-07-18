//
// ********************************************************************
// * Main function of the navier-stokes equation project
//
// ********************************************************************
//
#include "Globals.h"
#include "MemoryManager.h"
#include <stdio.h>

int main(void) {

  /* Declare all parameters */

  // user parameters (editable through UserInput.txt)
  int problem_type;             // constant (0) or periodic (1) boundary conditions
  double inflow_vel;            // magnitude of the inflow velocity
  double frequency;             // frequency of periodic boundary conditions
  int i_max, j_max;             // # of grid points in x- and y direction respectively
  double reynold;               // Reynolds number
  double a_size, b_size;        // sizes of the grid
  double tau_safety;            // safety factor tau
  double omega_relax;           // relaxation factor omega
  double epsilon_tolerance;     // relative tolerance for the residual epsilon
  int max_iterations;           // maximum # of iterations for the SOR-loop
  double max_int_time;          // maximum allowed integration time
  Gravity g_accel;              // gravitational acceleration
  // simulation parameters
  SimulationGrid ** cell;       // 2D-Array simulation parameters

  // Initialize the user parameters by parsing their values
  initialize_user_parameters(&problem_type, &inflow_vel, &frequency, &i_max, &j_max, &a_size, &b_size,
                            &g_accel, &reynold, &tau_safety, &omega_relax, &epsilon_tolerance,
                            &max_iterations, &max_int_time);

  // Allocate memory for all fields
  allocate_memory(&cell, i_max, j_max);
  printf("Memory allocated! \n");

  initialize_cells(cell, i_max, j_max);
  printf("Cells initialized to zero! \n");

  // Start of the time evolution loop

  // Deallocate the memory for the fields
  free_memory(&cell);
  printf("Memory deallocated. Goodbye :) \n");
}
