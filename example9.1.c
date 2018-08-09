#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  // Constants
  int M = 100;           // Grid square on a side
  double V = 1.0;        // Voltage at top wall
  double target = 1e-6;  // Target accuracy

  // calloc() initialize with 0's
  double **phi = calloc((M+1), sizeof(*phi));
  for (int i = 0; i < M+1; i++) {
    phi[i] = calloc((M+1), sizeof(*phi));
  }

  double **phiprime = malloc((M+1) * sizeof(*phiprime));
  for (int i = 0; i < M+1; i++) {
    phiprime[i] = malloc((M+1) * sizeof(*phiprime));
  }

  for (int i = 0; i < M+1; i++) {
    phi[0][i] = V;
  }

  double delta = 1.0;
  while (delta > target) {
    for (int i = 0; i < M+1; i++) {
      for (int j = 0; j < M+1; j++) {
        if (i == 0 || i == M || j == 0 || j == M) {
          phiprime[i][j] = phi[i][j];
        }

        else {
          phiprime[i][j] = phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1];
          phiprime /= 4; // a = 1cm
        }
      }
    }

    delta = MaxAbs(phi, phiprime);
    swap(phi, phiprime);

  }

  // Plot the phi matrix
  return 0;
}
