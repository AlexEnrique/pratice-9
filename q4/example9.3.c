#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define OUTPUT_FNAME "ex4.dat"
#define SCRIPT_NAME "script.gnu"
#define BUFF_SIZE 70

// Functions
void swap(double **Temp, double **__Temp, int Ngrid);
void createGnuplotScript(double tEnd);

int main(int argc, char *argv[]) {
  if (argc == 1) {
    printf( "========================================================================\n"   );
    printf( "= ERROR: you have to pass the end time via command line                =\n"   );
    printf( "= Call the programa and write the end time after the name              =\n"   );
    printf( "= \t(example: $ ./program.x 3.1415)                                =\n"       );
    printf( "= The above example calls the program passing \"3.1415\" as the end time =\n" );
    printf( "========================================================================\n"   );
    printf("\a");

    return 0;
  }

  // Constants
  int Ngrid = 100;                  // Grid
  double L = .01;                   // Thickness
  double a = L / Ngrid;             // spacing
  double h = 1e-4;                  // time-step
  double D = 4.25e-6;               // Therm. diffsvty.
  double tEnd = atof(*(argv + 1));  // Final time (to plot)

  double Thig = 50.0;               // Higher Temperature
  double Tlow = 0.0;                // Lower Temperature
  double Tintermd = 20.0;           // Intermediate (initial) Temperature

  double t = 0.0;                   // initial time
  double c = h * D / pow(a,2);      // c := hD/a^2; for simplicity at calculating Temp[]

  double *Temp = malloc(Ngrid * sizeof(*Temp));      // Temperature array
  double *__Temp = malloc(Ngrid * sizeof(*__Temp));  // Temporary Temperature array

  // Initializing Temperature array
  __Temp[0] = Thig;                      // First element of Temp
  __Temp[Ngrid - 1] = Tlow;              // Last element of __Temp
  Temp[0] = Thig;                        // First element of Temp
  Temp[Ngrid - 1] = Tlow;                // Last element of Temp
  for (int i = 1; i < Ngrid - 1; i++) {  // Intermediates elements of Temp
    Temp[i] = Tintermd;
  }
  // The Intermediates elements of __Temp are always calculated

  int count = 0; // this counter is used below to avoid swapping Temp and __Temp unnecessarily for each iteration
  do { // using do {...} while (...) to accept tEnd = 0.0 as possible argument
    if (count % 2 == 0) {
      for (int i = (0 + 1); i < (Ngrid - 1); i++) {
        __Temp[i] = Temp[i] + c * (Temp[i+1] + Temp[i-1] - 2 * Temp[i]);
      }
    } // end if

    else {
      for (int i = (0 + 1); i < (Ngrid - 1); i++) {
        Temp[i] = __Temp[i] + c * (__Temp[i+1] + __Temp[i-1] - 2 * __Temp[i]);
      }
    } // end else

    t += h;
    count++;
  } while (t < tEnd); // end do

  if (count % 2 == 1) // if count % 2 == 1, the last ocupied array was __Temp, not Temp, because there is a "count++" at the end of the do loop
    swap(&Temp, &__Temp, Ngrid);

  free(__Temp);

  // Output the calculated data to a file
  FILE *fp = fopen(OUTPUT_FNAME, "w");

  for (int i = 0; i < Ngrid; i++) {
    fprintf(fp, "%d  %lf\n", i, Temp[i]);
  }

  fclose(fp);
  free(Temp);

  createGnuplotScript(tEnd);
  system("gnuplot < script.gnu --persist");

  return 0;
}

void swap(double **Temp, double **__Temp, int Ngrid) {
  for (int i = 0; i < Ngrid; i++) {
    (*Temp)[i] = (*__Temp)[i];
  }
}

void createGnuplotScript(double tEnd) {
  char *buffer = malloc(BUFF_SIZE * sizeof(*buffer));
  FILE *script = fopen(SCRIPT_NAME, "w");

  fprintf(script, "set xlabel \"x\"\n");
  fprintf(script, "set ylabel \"T\"\n");

  fprintf(script, "set xtics auto\n");
  fprintf(script, "set ytics auto\n");

  snprintf(buffer, BUFF_SIZE, "plot '%s' title 'T(x,t), for t = %.3lf' with line lt 7 lw 1.7", OUTPUT_FNAME, tEnd);
  fprintf(script, "%s\n", buffer);

  fclose(script);
}
