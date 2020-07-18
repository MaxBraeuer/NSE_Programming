#include "MemoryManager.h"

#include <stdlib.h>
#include <stdio.h>

// we do _not_ use calloc because calloc sets all bits to zero, which does not necessarily result in a float value of 0.0
void allocate_memory(SimulationGrid ***cell, int i_max, int j_max) {
  *cell = malloc((i_max + 2) * sizeof(SimulationGrid));
  for (int i = 0; i < i_max + 2; i++) {
    (*cell)[i] = malloc((j_max + 2) * sizeof(SimulationGrid));
  }

  if (*cell) {
    return;
  } else {
    fprintf(stderr, "Error while allocating memory.");
    exit(EXIT_FAILURE);
  }
}

// Initialize all cells to 0.0 (float)
void initialize_cells(SimulationGrid **cell, int i_max, int j_max) {
  int i, j;
  for (i = 0; i < i_max + 2; i++) {
    for (j = 0; j < j_max + 2; j++) {
      cell[i][j].vel_u = 0.0;
      cell[i][j].vel_v = 0.0;
      cell[i][j].pressure = 0.0;
      cell[i][j].F_term = 0.0;
      cell[i][j].G_term = 0.0;
      cell[i][j].poisson_rhs = 0.0;
      cell[i][j].residual = 0.0;
    }
  }
}

void free_memory(SimulationGrid ***cell) {
  free(*cell);
}
