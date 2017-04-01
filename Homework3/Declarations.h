// This file declares various functions that are useful over the course
// of the semester.  It must be "included" above the main program.
// Version Spring 2016.


// These are some standard C function libraries.
#include <stdlib.h>  // "standard library"
#include <math.h>    // various math functions, such as exp()
#include <stdio.h>   // various "input/output" functions
#include <time.h>    // functions for timing computations

// The numerical equivalent of zero.
#define ZERO 0.00000001



// These functions are found in "Definitions.h".

// - Random number generator.
double MTUniform (unsigned int);

// - Cumulative normal functions.
double Psi (double);
double PsiInv (double);


// - Histogram functions.
void Histogram (double, double, double, int, int);
void NormalHistogram (double, int, int);
void UniformHistogram (double, int, int);


// - Miscellanceous functions.
void   Pause (void);
void   Exit (void);
int GetInteger (char *);
double GetDouble (char *);

// P-value functions.
double PValueForT (int, double);
double PValueForF (int, int, double);


// - Linear algebra functions.
double **Array (int, int);
void Show (double **);
void Write (double **, FILE *);
double **Transpose (double **);
double **Multiply (double **, double **);
double **Invert (double **);
double **Copy (double **);
double **Identity (int);
double Det(double **);
double **Cholesky (double **);
double **ScalarMultiple (double, double **);
double **Add (double **, double **);
int Rows (double **);
int Columns (double **);
void Free (double **);
double **Evalues (double **);
double **Evector (double **, double **);
double **MeanZero (double **);
double **Covariance (double **);
double **Correlation (double **);
