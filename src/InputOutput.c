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
 int error_count = 0;

 FILE* input_file;
 input_file = fopen("UserInput.txt", "r");

 // Error if the file is not found or corrupt
 if (input_file == NULL) {
   perror("Error while opening the file 'UserInput.txt'");
   exit(EXIT_FAILURE);
 }

 // Read the User Input file line by line with fgets to the buffer
 // and assign the corresponding values to their variables

 // In C you can't extract the type of a variable at runtime (only at compilation-time), so if the user
 // types in a float where the value is supposed to be an integer, we're kind of screwed (potentially).

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%d", problem_type);
 if (*problem_type != 0 && *problem_type != 1) {
   fprintf(stderr, "Error in 'UserInput.txt': The problem type has to be 0 or 1.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", inflow_vel);

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", frequency);

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%d", i_max);
 if (*i_max <= 0) {
   fprintf(stderr, "Error in 'UserInput.txt': i_max has to be an integer greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%d", j_max);
 if (*j_max <= 0) {
   fprintf(stderr, "Error in 'UserInput.txt': j_max has to be an integer greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", a_size);
 if (*a_size <= 0.0) {
   fprintf(stderr, "Error in 'UserInput.txt': side a has to have a length greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", b_size);
 if (*b_size <= 0.0) {
   fprintf(stderr, "Error in 'UserInput.txt': side b has to have a length greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", &g_accel->g_x);
 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", &g_accel->g_y);

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", reynold);
 if (*reynold <= 0.0) {
   fprintf(stderr, "Error in 'UserInput.txt': Reynolds number has to be greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", tau_safety);
 if (*tau_safety <= 0.0 || *tau_safety > 1.0) {
   fprintf(stderr, "Error in 'UserInput.txt': The safety factor tau has to be in the interval (0,1].\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", omega_relax);
 if (*omega_relax <= 0.0) {
   fprintf(stderr, "Error in 'UserInput.txt': The relaxation factor omega has to be greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", epsilon_tolerance);
 if (*epsilon_tolerance <= 0.0 || *epsilon_tolerance > 1.0) {
   fprintf(stderr, "Error in 'UserInput.txt': The tolerance factor epsilon has to be in the interval (0,1].\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%d", max_iterations);
 if (*max_iterations <= 0) {
   fprintf(stderr, "Error in 'UserInput.txt': The number of maximum SOR-iterations has to be an integer greater than zero.\n");
   error_count += 1;
 }

 fgets(buffer, 256, input_file);
 sscanf(buffer, "%lg", max_int_time);
 if (*max_int_time <= 0.0) {
   fprintf(stderr, "Error in 'UserInput.txt': The maximum integration time has to be greater than zero.\n");
   error_count += 1;
 }
 
 fgets(buffer, 256, input_file);
 sscanf(buffer, "%d", max_int_steps);
 if (*max_int_steps <= 0) {
   fprintf(stderr, "Error in 'UserInput.txt': The number of maximum integration steps has to be an integer greater than zero.\n");
   error_count += 1;
 }

 fclose(input_file);

 if (error_count > 0) {
   exit(EXIT_FAILURE);
 }
}

void export_cells(SimulationGrid **cell, int i_max, int j_max, int current_time_step_number, double current_int_time) {

 char filename_u[32];
 char filename_v[32];
 char filename_p[32];

 sprintf(filename_u, "output/u_%d.txt", current_time_step_number);
 sprintf(filename_v, "output/v_%d.txt", current_time_step_number);
 sprintf(filename_p, "output/p_%d.txt", current_time_step_number);

 FILE *out_vel_u, *out_vel_v, *out_pressure;

 out_vel_u = fopen(filename_u, "w");
 out_vel_v = fopen(filename_v, "w");
 out_pressure = fopen(filename_p, "w");

 // error if the folder output/ does not exist
 if (out_vel_u == NULL || out_vel_v == NULL || out_pressure == NULL) {
   fprintf(stderr, "Error while creating the output. Make sure the folder 'output' exists within the build directory.");
   exit(EXIT_FAILURE);
 }

 // write the array contents row by row to the text files
 for (int i = 0; i <= i_max + 1; i++) {
   for (int j = 0; j <= j_max + 1; j++) {
     fprintf(out_vel_u, "%lg\n", cell[i][j].vel_u);
     fprintf(out_vel_v, "%lg\n", cell[i][j].vel_v);
     fprintf(out_pressure, "%lg\n", cell[i][j].pressure);
   }
 }

 fclose(out_vel_u);
 fclose(out_vel_v);
 fclose(out_pressure);
}
