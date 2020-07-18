void time_evolution(int problem_type, double inflow_vel, double frequency, int i_max, int j_max, 
                   double reynold, double a_size, double b_size, double tau_safety, double omega_relax, 
                   double epsilon_tolerance, int max_iterations, Gravity g_accel, 
                   double max_int_time, SimulationGrid **cell);