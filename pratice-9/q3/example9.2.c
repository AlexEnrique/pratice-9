#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

void createGnuplotScript();
void swap(double ***phi, double ***newphi, int M);
double calcDelta(double **phi, double **newphi, int M);
double getMaxAbs(double **phi, int M);

int main() {
  // Constants
  int M = 100;           // Grid
  // double V = 0.0;     // Voltage at walls (unuseful - see line 17)
  double target = 1e-6;  // Target accuracy

  // calloc() initialize with 0's (we don't need to fill the walls's elements with V = 0.0)
  double **phi = calloc((M+1), sizeof(*phi));
  for (int i = 0; i < M+1; i++) {
    phi[i] = calloc((M+1), sizeof(*phi));
  }

  double **newphi = malloc((M+1) * sizeof(*newphi));
  for (int i = 0; i < M+1; i++) {
    newphi[i] = malloc((M+1) * sizeof(*newphi));
  }

  double delta = 1.0;
  while (delta > target) {

    for (int i = 0; i < M+1; i++) {
      for (int j = 0; j < M+1; j++) {
        if (i == 0 || i == M || j == 0 || j == M) { // at the walls V == 0.0 (original phi[i][j] - line 17)
          newphi[i][j] = phi[i][j];
        } // end if

        else {
          // a = 1cm
          newphi[i][j] = phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1];

          // rho (x, y) <=> rho(i, j); q_i = i * deltaQ (q == x or q == y)
          if (i >= 20 && i <= 40 && j >= 20 && j <= 40)
            newphi[i][j] += -1; // eps_0 == 1
          else if (i >= 60 && i <= 80 && j >= 60 && j <= 80)
            newphi[i][j] += 1; // eps_0 == 1

          newphi[i][j] /= 4;
        } // end else
      } // end 'j' for
    } // end 'i' for

    delta = calcDelta(phi, newphi, M);
    printf("target: %.2e -- delta: %.2e\n", target, delta);

    swap(&phi, &newphi, M);
  } // end while

  // to normalize the data
  double norm = getMaxAbs(phi, M);

  // Plot the phi matrix
  FILE *fp = fopen("ex3.dat", "w");

  for (int i = 0; i < M+1; i++) {
    for (int j = 0; j < M+1; j++) {
      fprintf(fp, "%d  %d  %lf\n", i, j, phi[i][j] / norm);
    }
  }

  fclose(fp);

  createGnuplotScript();
  system("gnuplot < script.gnu --persist");

  return 0;
}

void swap(double ***phi, double ***newphi, int M) {
  double aux;

  for (int i = 0; i < M+1; i++) {
    for (int j = 0; j < M+1; j++) {
      aux = (*phi)[i][j];
      (*phi)[i][j] = (*newphi)[i][j];
      (*newphi)[i][j] = aux;
    }
  }
}

double calcDelta(double **phi, double **newphi, int M) {
  double delta = -DBL_MAX;

  for (int i = 0; i < M+1; i++) {
    for (int j = 0; j < M+1; j++) {
      if ( fabs(phi[i][j] - newphi[i][j]) > delta )
        delta = fabs(phi[i][j] - newphi[i][j]);
    }
  }

  return delta;
}

double getMaxAbs(double **phi, int M) {
  double max = -DBL_MAX;

  for (int i = 0; i < M+1; i++) {
    for (int j = 0; j < M+1; j++) {
      if (phi[i][j] > max)
        max = phi[i][j];
    }
  }

  return max;
}

void createGnuplotScript() {
  FILE *script = fopen("script.gnu", "w");

  fprintf(script, "set palette gray\n");

  fprintf(script, "unset border\n");
  fprintf(script, "unset xtics\n");
  fprintf(script, "unset ytics\n");
  fprintf(script, "set pm3d\n");

  fprintf(script, "set xrange[0:100]\n");
  fprintf(script, "set yrange[0:100]\n");

  fprintf(script, "plot 'ex3.dat' using 1:2:3 w p lw 7 palette t ''\n");

  fclose(script);
}
