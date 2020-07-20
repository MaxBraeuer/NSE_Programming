void initialize_user_parameters(int *problem_type, double *inflow_vel, double *frequency, int *i_max, 
                                int *j_max, double *a_size, double *b_size, Gravity *g_accel, double *reynold,
                                double *tau_safety, double *omega_relax, double *epsilon_tolerance,
                                int *max_iterations, double *max_int_time, int *max_int_steps);
								
void export_cells(SimulationGrid **cell, int i_max, int j_max, int current_time_step_number, double current_int_time);