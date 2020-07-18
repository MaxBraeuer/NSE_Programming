#include "Globals.h"
#include "TimeEvolution.h"
#include "BoundaryConditions.h"
#include "MainCalculations.h"
#include "InputOutput.h"

#include <math.h>
#include <stdio.h>


void time_evolution(int problem_type, double inflow_vel, double frequency, int i_max, int j_max, 
                   double reynold, double a_size, double b_size, double tau_safety, double omega_relax, 
                   double epsilon_tolerance, int max_iterations, Gravity g_accel, 
                   double max_int_time, SimulationGrid **cell) {
					 
  // Declare some new needed variables
  double vel_u_max_abs, vel_v_max_abs;  // maximum velocity components (absolute value)
  double gamma_weight;                  // factor for donor-cell-stencil weighing
  double delta_time;                    // size of the timestep
  double int_time = 0.0;                // initialize the integration time  
  double delta_x = a_size / i_max;      // cell size in x-direction
  double delta_y = b_size / j_max;      // cell size in y-direction
  
  while (int_time < max_int_time) {
    // calculate adaptive size of the timestep to ensure a stable simulation: eq. (26)
    vel_u_max_abs = ;
    vel_v_max_abs = ;

	delta_time = ;

    // implement boundary conditions: eq. (8-10)
    // the right, left and bottom boundaries are set using the no-slip condition
	//

    // the top boundary is set using the inflow condition;
    // depending on the problem type the velocity x-component is either constant or varies 
	//

    // calculate F and G: eq. (29-30)
    gamma_weight = fmax(vel_u_max_abs * delta_time / delta_x, vel_v_max_abs * delta_time / delta_y);
	
    calculate_F_G(cell, i_max, j_max, delta_x, delta_y, delta_time, reynold, gamma_weight, g_accel);

    // Calculate the RHS of the pressure eq. (40)
    rhs_pressure(cell, i_max, j_max, delta_time, delta_x, delta_y);

    // begin SOR-loop
    sor_loop(cell, i_max, j_max, delta_x, delta_y, omega_relax, epsilon_tolerance, max_iterations);

    // calculate new velocity components u_vel and v_vel: eq. (21-22)
    calculate_velocities(cell, i_max, j_max, delta_time, delta_x, delta_y);
    printf("Velocities updated!\n");

    // increment the integration time by the size of the adaptive timestep
    int_time += delta_time;
    printf("Integration time elapsed: %lg / %lg \n", int_time, max_int_time);

    // export the cells
	//

  }
}