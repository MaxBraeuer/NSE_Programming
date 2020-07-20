#include "Globals.h"
#include "InputOutput.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
void initialize_user_parameters(int *problem_type,double *inflow_vel, double *frequency, int *i_max, int *j_max, 
                                double *a_size, double *b_size, Gravity *g_accel, double *reynold,
                                double *tau_safety, double *omega_relax, double *epsilon_tolerance,
                                int *max_iterations, double *max_int_time, int *max_int_steps) {
  char buffer[256];

  FILE* input_file;
  input_file = fopen("UserInput.txt", "r");

  // Error if the file is not found or corrupt
  if (input_file == NULL) {
    perror("Error while opening the file 'UserInput.txt'");
    exit(EXIT_FAILURE);
  }

  // Read the User Input file line by line with fgets to the buffer 
  // and assign the corresponding values to their variables

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%d", problem_type);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", inflow_vel);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", frequency);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%d", i_max);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%d", j_max);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", a_size);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", b_size);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", &g_accel->g_x);
  
  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", &g_accel->g_y);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", reynold);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", tau_safety);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", omega_relax);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", epsilon_tolerance);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%d", max_iterations);

  fgets(buffer, 256, input_file);
  sscanf(buffer, "%lg", max_int_time);
  
  fgets(buffer, 256, input_file);
  sscanf(buffer, "%d", max_int_steps);
  
  fclose(input_file);

}

/*
void export_cells() {
}
*/

