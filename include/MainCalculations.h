void calculate_F_G(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y,
                  double delta_time, double reynold, double gamma_weight, Gravity g_accel);

void rhs_pressure(SimulationGrid **cell, int i_max, int j_max, double delta_time, double delta_x, double delta_y);

void sor_loop(SimulationGrid **cell, int i_max, int j_max, double delta_x, double delta_y,
             double omega_relax, double epsilon_tolerance, int max_iterations);

void calculate_velocities(SimulationGrid **cell, int i_max, int j_max, double delta_time, double delta_x, double delta_y);
