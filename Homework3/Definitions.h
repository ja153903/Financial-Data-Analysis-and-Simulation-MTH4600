// This file contains various functions that are useful over the course
//   of the semester.  They must be declared above the main program and this
//   file must be "included" following the main program.
// Version Spring 2016.


// First are some random number generators.

////////////////////////////////////////////////////////////////////////////////
// MERSENNE TWISTER
// By M. Matsumoto and T. Nishimura (1998).
// "Mersenne Twister: a 623-dimensionally equidistributed uniform pseudo-random
//   number generator".
// ACM Transactions of Modeling and Computer Simulation 8(1):3-30.
// Any coding errors introduced are my own (C.D. Howard).

// An unsigned integer is represented by 32 bits in base 2.  The largest unsigned integer is:
// 2^32 - 1 = 4,294,967,295 (base 10);
//          = ffffffff (base 16);
//          = 1111 1111 1111 1111 1111 1111 1111 1111 (base 2).
// The digits in hexadecimal (base 16) are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f.
double MTUniform (unsigned int seed) {

   static unsigned int X[1248], m[2], initialized = 0, k;
   unsigned int N, Y;
   unsigned int Temper (unsigned int);

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      X[0] = (seed ? seed : 1);
      // Seed X[1],X[2],...,X[623] with your favorite LCG.
      for (k = 1; k < 624; k++) {
         X[k] = 22695477 * X[k-1] + 1;
      }
      m[0] = 0; m[1] = 0x9908b0df;
      // The counter "k" is now 624.
      initialized = 1;
   }

   // Concatenate the first bit of X[k-624] with the last 31 bits of X[k-623],
   //    "|" is "bit-wise or".
   Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

   // Now apply the invertible linear transformation A to Y and bit-wise add X[k-227].
   X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

   // Re-load X[0],X[1],...,X[623] as you go.
   X[k-624] = X[k];

   // Apply the tempering function.
   N = Temper (X[k]);

   // Increment the counter; shift vectors when appropriate.
   k ++;
   if (k == 1248) k = 624;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}




////////////////////////////////////////////////////////////////////////////////
// Tempering function used by the Mersenne Twister.
unsigned int Temper (unsigned int N) {

   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   return N;

}






////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then CONTINUES the program.
void Pause () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to continue program... ");
   fgets (input, 9, stdin);

   return;

}



////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then EXITS the program.
// It prevents the window from closing up before the output can be viewed.
void Exit () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to exit program... ");
   fgets (input, 9, stdin);

   exit (0);

}



////////////////////////////////////////////////////////////////////////////////
// Read in an integer from the screen in response to a question.
////////////////////////////////////////////////////////////////////////////////
int GetInteger (char *question) {

   int i;
   char input[100];

   printf (question);
   fgets (input, 99, stdin);
   sscanf (input, "%d", &i);
   return i;

}



////////////////////////////////////////////////////////////////////////////////
// Read in a double from the screen in response to a question.
////////////////////////////////////////////////////////////////////////////////
double GetDouble (char *question) {

   double x;
   char input[100];

   printf (question);
   fgets (input, 99, stdin);
   sscanf (input, "%lf", &x);
   return x;

}


////////////////////////////////////////////////////////////////////////////////
// Allocate array space for an m x n array.
////////////////////////////////////////////////////////////////////////////////
double **Array (int m, int n) {

   int i;
   double **A;

   A = (double **) calloc (m+1, sizeof (double *));
   for (i = 0; i <= m; i++) {
      A[i] = (double *) calloc (n+1, sizeof (double));
   }

   // Record dimensions in the 0^th row, which is unused in matrix operations.
   A[0][0] = m;
   A[0][1] = n;

   return (A);

}


////////////////////////////////////////////////////////////////////////////////
// Multiply the matrix A by the scalar c.
////////////////////////////////////////////////////////////////////////////////
double **ScalarMultiple (double c, double **A) {

   int n, m, i, j;
   double **cA;

   m = Rows (A);
   n = Columns (A);

   cA = Array (m, n);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         cA[i][j] = c * A[i][j];
      }
   }

   return cA;

}


////////////////////////////////////////////////////////////////////////////////
// Compute the transpose of the matrix A.
////////////////////////////////////////////////////////////////////////////////
double **Transpose (double **A) {

   int i, j, n, m;
   double **At;

   m = Rows (A);
   n = Columns (A);

   At = Array (n,m);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
         At[i][j] = A[j][i];
      }
   }

   return At;

}


////////////////////////////////////////////////////////////////////////////////
// Use Gaussian elimination to invert the square matrix A0.
////////////////////////////////////////////////////////////////////////////////
double **Invert (double **A0) {

   int i, j, k, n, rmax;
   double **A, **Ainv, c;

   // Make a copy so that the original matrix is not changed during G.E.
   A = Copy (A0);

   n = Rows (A);

   if (n != Columns (A)) {
      printf ("Trying to invert a non-square matrix.\n");
      Pause ();
      exit(1);
   }

   // Start with the n x n identity matrix.  This matrix will eventually hold
   //    the inverse of A.
   Ainv = Identity (n);

   // Work on column j of A.
   for (j = 1; j <= n; j++) {

      // Find the largest non-zero entry in the column on or below the diagonal.
      c = 0;
      rmax = 0;
      for (i = j; i <= n; i++) {
         if (fabs(A[i][j]) > c) {
            c = fabs (A[i][j]);
            rmax = i;
         }
      }

      // If they are all 0 the matrix is singular.
      if (c < ZERO) {
         printf ("Trying to invert a singular matrix.\n");
         Pause ();
         exit(1);
      }

      // Swap rows j and rmax in both A and Ainv.
      i = j;
      for (k = 1; k <= n; k++) {
         c = A[i][k];
         A[i][k] = A[rmax][k];
         A[rmax][k] = c;
         c = Ainv[i][k];
         Ainv[i][k] = Ainv[rmax][k];
         Ainv[rmax][k] = c;
      }

      // Scale so the pivot is 1.
      c = A[i][i];
      for (k = 1; k <= n; k++) {
         A[i][k] /= c;
         Ainv[i][k] /= c;
      }

      // Make rest of column j equal to 0. Apply same row operations to Ainv.
      for (i = 1; i <= n; i++) if (i != j) {
         c = A[i][j];
         for (k = 1; k <= n; k++) {
            A[i][k] += -c * A[j][k];
            Ainv[i][k] += -c * Ainv[j][k];
         }
      }

   }

   // We're done with A.
   Free (A);

   return Ainv;

}

////////////////////////////////////////////////////////////////////////////////
// Use Gaussian elimination to compute the determinant of A0.
////////////////////////////////////////////////////////////////////////////////
double Det (double **A0) {

   int i, j, k, n, rmax;
   double **A, c, detA;

   // Make a copy so the original matrix is unchanged.
   A = Copy (A0);

   // Make sure the matrix is square.
   n = Rows (A);
   if (n != Columns (A)) {
      printf ("Trying to get determinant of a non-square matrix.\n");
      Pause ();
      exit(1);
   }

   // Initialize the determinant value.
   detA = 1;

   // Work on column j of A.
   for (j = 1; j <= n; j++) {

      // Find the largest non-zero entry in the column on or below the diagonal.
      c = 0;
      rmax = 0;
      for (i = j; i <= n; i++) {
         if (fabs(A[i][j]) > c) {
            c = fabs (A[i][j]);
            rmax = i;
         }
      }

      // If all such entries are 0 the matrix is singular, so Det(A) = 0.
      if (c < ZERO) {
         return 0;
      }

      i = j;
      // Swap row i and row rmax if different. Here det E = -1.
      if (rmax != i) {
         for (k = 1; k <= n; k++) {
            c = A[i][k];
            A[i][k] = A[rmax][k];
            A[rmax][k] = c;
         }
         // Update the determinant.
         detA *= -1;
      }

      // Scale so the pivot is 1. Here det E = 1/c.
      c = A[i][i];
      for (k = 1; k <= n; k++) {
         A[i][k] /= c;
      }
      // Update the determinant.
      detA *= c;

      // Make rest of column j equal to 0. For these row operations det E = 1.
      for (i = 1; i <= n; i++) if (i != j) {
         c = A[i][j];
         for (k = 1; k <= n; k++) {
            A[i][k] += -c * A[j][k];
         }
      }

   }

   // We're done with A.
   Free (A);

   return detA;

}


////////////////////////////////////////////////////////////////////////////////
// Compute the eigenvectors of V0 corresponding to the eigenvalues that are held
//   in E[1][1], E[1][2], ... , E[1][n]. This is done by solving the system
//   of equations V0*x = lambda*x, where lambda is the eigenvalue in question.
//   The eigenvectors are put into the columns of the matrix Q.
////////////////////////////////////////////////////////////////////////////////
double **Evector (double **V0, double **E) {

   int i, j, k, n, r, rmax;
   double **V, **Q, c;

   // Make sure that A0 is square and symmetric.
   n = Rows (V0);
   if (n != Columns(V0)) {
      printf ("In Eigenvalue function V0 is not square.\n");
      Pause ();
      exit (1);
   }

   for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
         if (V0[i][j] != V0[j][i]) {
            printf ("In Eigenvalue function V0 is not symmetric.\n");
            Pause ();
            exit (1);
         }
      }
   }

   // Allocate working array space.
   V = Array (n,n);
   Q = Array (n,n);


   // Work on eigenvalue number r, which is E[1][r].
   for (r = 1; r <= n; r++) {

      // Subtract the eigenvalue from the diagonal of the matrix V0, put the
      //    result in V.
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= n; j++) {
            V[i][j] = (i == j ? V0[i][j] - E[1][r] : V0[i][j]);
         }
      }


      // Use Gaussian elimination to solve the system V*x=0.
      // Work on column j.
      for (j = 1; j < n; j++) {

      // Find the largest non-zero entry in the column on or below the diagonal.
         c = 0;
         rmax = 0;
         for (i = j; i <= n; i++) {
            if (fabs(V[i][j]) > c) {
               c = fabs (V[i][j]);
               rmax = i;
            }
         }

         if (c < ZERO) {
            printf ("Eigenspace has dimension >= 2\n");
            Pause ();
            exit (1);
         }

         i = j;
         // Swap rows i and rmax if different.
         if (rmax != i) {
            for (k = 1; k <= n; k++) {
               c = V[i][k];
               V[i][k] = V[rmax][k];
               V[rmax][k] = c;
            }
         }


         // Scale so the pivot is 1.
         c = V[i][i];
         for (k = 1; k <= n; k++) {
            V[i][k] /= c;
         }

         // Make rest of column j equal to 0.
         for (i = 1; i <= n; i++) if (i != j) {
            c = V[i][j];
            for (k = 1; k <= n; k++) {
               V[i][k] += -c * V[j][k];
            }
         }

      }

      // The eigenvector is now (0,0,...,0,1)^T minus the last column of V.
      // Compute the norm, c, of this vector so that we can re-scale it to 1.
      c = 0;
      for (i = 1; i <= n - 1; i++) {
         c += V[i][n] * V[i][n];
      }
      c += 1;
      c = sqrt(c);

      // Re-scale it and put it into column r of the matrix Q.
      for (i = 1; i <= n - 1; i++) {
         Q[i][r] = -V[i][n] / c;
      }
      Q[n][r] = 1.0 / c;

   } // Next r.

   // We're done with V.
   Free (V);

   return Q;

}


////////////////////////////////////////////////////////////////////////////////
// Compute the matrix product AB.
////////////////////////////////////////////////////////////////////////////////
double **Multiply (double **A, double **B) {

   int i, j, k, nA, mA, nB, mB;
   double **AB;

   mA = Rows (A);
   nA = Columns (A);

   mB = Rows (B);
   nB = Columns (B);

   if (nA != mB) {
      printf ("Dimensions don't match in matrix multiplication.\n");
      Pause ();
      exit(1);
   }

   // Allocate space for the product.
   AB = Array (mA, nB);

   // Compute (AB)_ij.
   for (i = 1; i <= mA; i++) {
      for (j = 1; j <= nB; j++) {
         for (k = 1; k <= nA; k++) {
            AB[i][j] += A[i][k] * B[k][j];
         }
      }
   }

   return AB;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the matrix sum A+B.
////////////////////////////////////////////////////////////////////////////////
double **Add (double **A, double **B) {

   int i, j, n, m;
   double **S;

   m = Rows (A);
   n = Columns (A);

   // Check that the dimensions match.
   if ( (m != Rows (B)) || (n != Columns (B)) ) {
      printf ("Dimensions don't match in matrix addition.\n");
      Pause ();
      exit(1);
   }

   // Allocate space for the sum.
   S = Array (m, n);

   // Compute (A+B)_ij
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         S[i][j] = A[i][j] + B[i][j];
      }
   }

   return S;

}


////////////////////////////////////////////////////////////////////////////////
// Compute the Cholesky decomposition V = LL^T for the symmetric PD matrix V.
////////////////////////////////////////////////////////////////////////////////
double **Cholesky (double **V) {

   int n, i, j, k;
   double **L, sum;

   // Check that the matrix V is square and symmetric.
   n = Rows (V);
    if (n != Columns(V)) {
      printf ("In Cholesky function V is not square.\n");
      Pause ();
      exit (1);
   }

   for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
         if (V[i][j] != V[j][i]) {
            printf ("In Cholesky function V is not symmetric.\n");
            Pause ();
            exit (1);
         }
      }
   }

   // Allocate space for L.
   L = Array (n, n);

   // Work on row k+1 of L.
   for (k = 0; k < n; k++) {

      for (j = 1; j <= k; j++) {
         sum = 0;
         for (i = 1; i < j; i++) {
            sum += L[k+1][i] * L[j][i];
         }
         L[k+1][j] = (V[k+1][j] - sum) / L[j][j];
      }

      sum = V[k+1][k+1];
      for (i = 1; i < k+1; i++) {
         sum -= L[k+1][i] * L[k+1][i];
      }

      if (sum <= ZERO) {
         printf ("In Cholesky function V is not positive definite.\n");
         Pause ();
         exit (1);
      }

      L[k+1][k+1] = sqrt (sum);

   }

   return L;

}


////////////////////////////////////////////////////////////////////////////////
// Print the matrix A to the screen.
////////////////////////////////////////////////////////////////////////////////
void Show (double **A) {

   int i, j, m, n;

   m = Rows (A);
   n = Columns (A);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         printf ("%10.4f", A[i][j]);
      }
      printf ("\n");
   }

   printf ("\n");
   Pause ();
   printf ("\n");

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Print the matrix A to a file.
////////////////////////////////////////////////////////////////////////////////
void Write (double **A, FILE *fp) {

   int i, j, m, n;

   m = Rows (A);
   n = Columns (A);

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         fprintf (fp, "%10.4f", A[i][j]);
      }
      fprintf (fp, "\n");
   }

   return;

}


////////////////////////////////////////////////////////////////////////////////
// Copy the matrix A into another matrix (C).
////////////////////////////////////////////////////////////////////////////////
double **Copy (double **A) {

   double **C;
   int i, j, m, n;

   // Allocate space for C.
   m = Rows (A);
   n = Columns (A);
   C = Array (m, n);

   // Copy A_ij into C_ij.
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         C[i][j] = A[i][j];
      }
   }

   return C;

}


////////////////////////////////////////////////////////////////////////////////
// Generate the n x n identity matrix.
////////////////////////////////////////////////////////////////////////////////
double **Identity (int n) {

   double **I;
   int i, j;

   I = Array (n, n);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
         I[i][j] = (i == j);
      }
   }

   return I;

}


////////////////////////////////////////////////////////////////////////////////
// Report the number of rows in the matrix A.  This is held in A[0][0].  Note
//   that row 0 of A is not used in any matrix computations.
////////////////////////////////////////////////////////////////////////////////
int Rows (double **A) {

   return (int) A[0][0];

}

////////////////////////////////////////////////////////////////////////////////
// Report the number of columns in the matrix A. Held in A[0][1].
////////////////////////////////////////////////////////////////////////////////
int Columns (double **A) {

   return (int) A[0][1];

}



////////////////////////////////////////////////////////////////////////////////
// Release allocated space for the matrix A.
////////////////////////////////////////////////////////////////////////////////
void Free (double **A) {

   int i, m;

   m = Rows (A);
   for (i = 0; i <= m; i++) {
      free (A[i]);
   }
   free (A);

   return;

}


////////////////////////////////////////////////////////////////////////////////
// Compute the eigenvalues for the PD matrix V. Put the result in
//    E[1][1], E[1][2], ... , E[1][n].  Do this by finding the roots of the
//    charateristic polynomial det(V - x*I). Order them from largest to smallest.
////////////////////////////////////////////////////////////////////////////////
double **Evalues (double **V) {

   int n, i, j, k, n0, N = 100000;
   double **tildeV, **E;
   double trace, x, y, xlast, ylast, dx, root;
   FILE *fp;

   // Make sure V is square...
   n = Rows (V);
   if (n != Columns(V)) {
      printf ("In Eigenvalue function V is not square.\n");
      Pause ();
      exit (1);
   }

   // ... and symmetric.
   for (i = 1; i <= n; i++) {
      for (j = 1; j < i; j++) {
         if (V[i][j] != V[j][i]) {
            printf ("In Eigenvalue function V is not symmetric.\n");
            Pause ();
            exit (1);
         }
      }
   }

   // Allocate array space for the eigenvalues.
   E = Array (1,n);

   // Open an output file to report the characteristic polynomial (c.p.) for
   //    viewing using TeX software.
   fp = fopen ("CharPolyData.txt", "w");

   // Make a copy of V.
   tildeV = Copy (V);

   // Compute the trace of V.  Since all the eigenvalues are positive and
   //   sum to trace(V), all the roots of the c.p. must lie between 0 and
   //   trace(V).
   trace = 0;
   for (i = 1; i <= n; i++) {
      trace += V[i][i];
   }

   // Compute c.p.(x) for N values of x between 0 and trace(V).
   dx = trace / N;

   // n0 is the number of roots found so far.
   n0 = 0;

   // Compute the y-intercept of the c.p.
   xlast = 0;
   ylast = Det (V);

   // Search for the roots.
   for (i = 1; i <= N; i++) {

      // Update the x value.
      x = i * dx;

      // Compute Det (V - x*I). Put it into y.
      for (j = 1; j <= n; j++) {
         tildeV[j][j] = V[j][j] - x;
      }
      y = Det (tildeV);

      // If reasonably small, report it to the output file.
      if (fabs(y) <= 5) {
         fprintf (fp, "%8.5f %8.5f\n", x, y);
      }

      // Look for a zero or a change in sign --- either indicates the presence
      //    of a root.
      if (y == 0.0 || ylast * y < 0.0) {

         // Augment the number of roots found.
         n0 ++;

         // Interpolate between (xlast, ylast) and (x,y), which are on opposite
         //   sides of the x-axis, to approximate the location of the root.
         root = xlast - ylast * (x-xlast) / (y-ylast);

         // Record the eigenvalues in the E[1][*] array from largest to smallest.
         k = n + 1 - n0;
         E[1][k] = root;

      }

      // Quit when all n roots have been found.
      if (n0 == n) {
         break;
      }

      // Otherwise, update (xlast,ylast) to be the current value of (x,y) and
      //    continue the search.
      xlast = x;
      ylast = y;

   }

   // Close up the output file.
   fclose (fp);

   // We're done with tildeV.
   Free (tildeV);

   return E;

}


////////////////////////////////////////////////////////////////////////////////
// Make the columns of an array have mean 0.
////////////////////////////////////////////////////////////////////////////////
double **MeanZero (double **A) {

   double **C, **mu;
   int i, j, m, n;

   // Allocate space for C.
   m = Rows (A);
   n = Columns (A);
   C = Array (m, n);
   mu = Array (1, n);

   // Compute the means of each column.
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         mu[1][j] += A[i][j] / m;
      }
   }

   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         C[i][j] = A[i][j] - mu[1][j];
      }
   }

   Free (mu);

   return C;

}

double **Covariance (double **X0) {

   int n;
   double **X, **XT, **C, **Cov;

   n = Rows (X0);

   X = MeanZero (X0);

   XT = Transpose(X);

   C = Multiply (XT, X);

   Cov = ScalarMultiple (1.0/(n-1.0), C);

   Free (X);
   Free (XT);
   Free (C);

   return Cov;

}


double **Correlation (double **X0) {

   int m, n, i, j;
   double s_ii;
   double **X, **XT, **C, **Cor;

   m = Rows (X0);
   n = Columns (X0);

   X = MeanZero (X0);

   XT = Transpose(X);

   C = Multiply (XT, X);

   Cor = ScalarMultiple (1.0/(m-1.0), C);

   for (i = 1; i <= n; i++) {
      s_ii = sqrt(Cor[i][i]);
      for (j = 1; j <= n; j++) {
         Cor[i][j] /= s_ii;
         Cor[j][i] /= s_ii;
      }
   }

   Free (X);
   Free (XT);
   Free (C);

   return Cor;

}


////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of randomly generated numbers. It is
// typically used with continuous (not discrete)  data.
////////////////////////////////////////////////////////////////////////////////
void   Histogram (double x,          // Item of data to be added to the histogram.
                  double a,          // Lower bound of histogram.
                  double b,          // Upper bound of histogram.
                  int n,             // Number of "buckets".
                  int done)          // 0 if not done, 1 if done.

{

   int k;
   double x0, x1, biggest;
   FILE *fp;

   static int initialized=0;
   static double dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

   }


   // If finished, create the TeX files that generate the plot. This is the
   // final call to the histogram function.
    else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "      From         To  Frequency\n");
      fprintf (fp, "========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f\n", x0, x1, freq[k]);
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Generate TeX output files.
      // Scale so that the biggest bucket has 1 (for graphing convenience).

      // First find the largest bucket value.
      biggest = 0;
      for (k = 0; k < n; k++) {
         if (freq[k] > biggest) {
            biggest = freq[k];
         }
      }

      // Now re-scale all the bucket values.
      for (k = 0; k < n; k++) {
         freq[k] /= biggest;
      }

      // Report data to the output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1) * dx, freq[k]);
      }
      fclose (fp);

      // Free up the freq[] array.
      free (freq);

      // No longer initialized.
      initialized = 0;

      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <%8.3f truein, 2.5 truein>\n", 4.0 / (b-a));
      fprintf (fp, "\\setplotarea x from %8.3f to %8.3f, y from  0 to 1.0\n", a, b);
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from %8.2f to %8.2f by %8.2f\n", a, b, (b-a)/10);
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot %8.3f 0  %8.3f 0  %8.3f 1  %8.3f 1 /\n",
                a-.045*(b-a), a-.03*(b-a), a-.03*(b-a), a-.045*(b-a) );
      fprintf (fp, "\\put {0} [cr] at %8.3f  0\n", a-.055*(b-a));
      fprintf (fp, "\\put {%d} [cr] at %8.3f  1\n", (int) biggest, a-.055*(b-a));
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\put {\\sl Histogram for a Continuous Random Variable} at %8.3f 1.2\n", (a+b)/2.0);
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);


   }

   return;

}



////////////////////////////////////////////////////////////////////////////////
// This histogram function is specially designed for data that should be
// Normal(0,1) in distribution.  The function computes a properly scaled
// standard normal density function for comparison to the data.
////////////////////////////////////////////////////////////////////////////////
void NormalHistogram (double x,          // Data to be added to the histogram.
                      int n,             // Number of "buckets".
                      int done)          // 0 if not done, 1 if done.

{

   int k;
   double z, x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the standard normal.
      a = -5; b = 5;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // ...otherwise, when finished, create the TeX output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n",
                       x0, x1, freq[k], N * (Psi(x1) - Psi(x0)));
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <0.45 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from -5 to 5, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from -5 to 5 by 1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot \"Normal.txt\"\n");
      fprintf (fp, "\\plot -5.2 0 -5.1 0 -5.1 1  -5.2 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -5.3 0\n");
      fprintf (fp, "\\put {$\\frac{1}{\\sqrt{2\\pi}}$} [cr] at -5.3 1\n");
      fprintf (fp, "\\put {\\sl Standard Normal Histogram} at 0 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is sqrt(2pi) (for comparing to the Normal(0,1) density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      area /= sqrt(2*3.14159);
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram to an output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);

      // Report the appropriately scaled normal density function for comparison.
      fp = fopen ("Normal.txt", "w");
      for (z = -5; z <= 5; z += .01) {
         fprintf (fp, "%8.3f %8.4f\n", z, exp(-z*z/2));
      }
      fclose (fp);

   }

   return;

}




////////////////////////////////////////////////////////////////////////////////
// This program creates a histogram of the data x, and plots the Uniform[0,1]
// density function for comparison.
////////////////////////////////////////////////////////////////////////////////
void UniformHistogram (double x,         // data to be added to the histogram
                      int n,             // number of "buckets"
                      int done)          // 1 when done, 0 otherwise

{

   int k;
   double x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the Uniform[0,1].
      a = 0; b = 1;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // When finished, create the output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1) * dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n", x0, x1, freq[k], N*dx);
      }
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <4.5 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from 0 to 1, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from 0 to 1 by .1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot 0 1  1 1 /\n");    // the Uniform[0,1] density function
      fprintf (fp, "\\plot -.06 0 -.05 0 -.05 1  -.06 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -.07 0\n");
      fprintf (fp, "\\put {1} [cr] at -.07 1\n");
      fprintf (fp, "\\put {\\sl Uniform (0,1) Histogram} at 0.5 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is 1 (for comparing to the uniform[0,1] density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram data.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);


   }

   return;

}


////////////////////////////////////////////////////////////////////////////////
// Normal (cumulative) distribution function evaluated at "x".
// It is accurate to 7 decimal places for all values of x.
////////////////////////////////////////////////////////////////////////////////
double Psi (double x) {

   static double    A =    0.0293868870,
                    B =    0.2161934875,
                    C =    0.6503029656,
                    D =    0.7978845608,
                    E =    0.0594864800,
                    F =    0.4160822924;

   int sign=1;
   double R;

   if (x < 0) {
      sign = -1;
      x *= -1;
   }

   x = fabs(x);
   R = (((A*x+B)*x+C)*x+D)*x / ((E*x+F)*x+1);

   // Return Psi (x)
   return (0.5 + sign * 0.5 * (1.0 - exp (-R)));

}


////////////////////////////////////////////////////////////////////////////////
// This algorithm is due to Peter John Acklam. Any coding errors are my own.
// It is a very good approximation of Psi^(-1).
////////////////////////////////////////////////////////////////////////////////
double PsiInv (double u) {

   static double
    A1 =  -3.969683028665376e+01,
    A2 =   2.209460984245205e+02,
    A3 =  -2.759285104469687e+02,
    A4 =   1.383577518672690e+02,
    A5 =  -3.066479806614716e+01,
    A6 =   2.506628277459239e+00,
    B1 =  -5.447609879822406e+01,
    B2 =   1.615858368580409e+02,
    B3 =  -1.556989798598866e+02,
    B4 =   6.680131188771972e+01,
    B5 =  -1.328068155288572e+01,
    C1 =  -7.784894002430293e-03,
    C2 =  -3.223964580411365e-01,
    C3 =  -2.400758277161838e+00,
    C4 =  -2.549732539343734e+00,
    C5 =   4.374664141464968e+00,
    C6 =   2.938163982698783e+00,
    D1 =   7.784695709041462e-03,
    D2 =   3.224671290700398e-01,
    D3 =   2.445134137142996e+00,
    D4 =   3.754408661907416e+00,
    P0 =  0.02425,
    P1 =  0.97575;

   double N, q, r;

   // Left tail.
   if (u < P0) {
      q = sqrt(-2*log(u));
      N = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   // Middle, general case.
   else if (u <= P1) {
      q = u - 0.5;
      r = q*q;
      N = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
   }

   // Right tail.
   else {
      q = sqrt(-2*log(1.0-u));
      N = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   return (N);

}

double PValueForT (int n, double x) {

   int i, N = 10000;

   double sum, c, P, L, dx, y, epsilon = 0.00001;
   double PdfT (int, double, double);

   // The pdf for the t distribugtion with n degrees of freedom is
   //  pdf(x) = c * (1 + x^2/n)^(-(n+1)/2) for the normalizing factor c.


   // First use Simpson's rule to estimate c.
   L = y = 1.0;
   while (y > epsilon) {
      L ++;
      y = PdfT (n, 1.0, L);
   }
   dx = L / N;

   // Endpoints.
   sum = (PdfT(n, 1.0, 0.0) + PdfT(n, 1.0, L)) * dx / 3.0;

   for (i = 1; i < N; i++) {
      sum += (i % 2 ? 4 : 2) * PdfT(n, 1.0, i*dx) * dx / 3.0;
   }

   c = 1.0 / (2.0 * sum);

   // Now use Simpson's rule to estimate the integral from 0 to x.
   dx = x / N;

   // Endpoints.
   sum = (PdfT(n, c, 0.0) + PdfT(n, c, x)) * dx / 3.0;

   for (i = 1; i < N; i++) {
      sum += (i % 2 ? 4 : 2) * PdfT(n, c, i*dx) * dx / 3.0;
   }

   // Return p-value for a two-sided test.
   P = 1.0 - 2.0*sum;

   return P;

}

double PdfT (int n, double c, double x) {

   return c * pow (1.0 + x*x/n, -(n+1)/2.0);

}

double PValueForF (int n, int d, double x) {

   int i, N = 100000;

   double sum, c, P, L, dx, y, epsilon = 0.00000001;
   double PdfF (int, int, double, double);

   // The pdf for the t distribugtion with n degrees of freedom is
   //  pdf(x) = c * (1 + x^2/n)^(-(n+1)/2) for the normalizing factor c.


   // First use Simpson's rule to estimate c.
   L = y = 1.0;
   while (y > epsilon) {
      L ++;
      y = PdfF (n, d, 1.0, L);
   }
   dx = L / N;

   // Endpoints.
   sum = (PdfF(n, d, 1.0, epsilon) + PdfF(n, d, 1.0, L)) * dx / 3.0;

   for (i = 1; i < N; i++) {
      sum += (i % 2 ? 4 : 2) * PdfF(n, d, 1.0, i*dx) * dx / 3.0;
   }

   c = 1.0 / sum;

   // Now use Simpson's rule to estimate the integral from 0 to x.
   dx = x / N;

   // Endpoints.
   sum = (PdfF(n, d, c, epsilon) + PdfF(n, d, c, x)) * dx / 3.0;

   for (i = 1; i < N; i++) {
      sum += (i % 2 ? 4 : 2) * PdfF(n, d, c, i*dx) * dx / 3.0;
      //printf ("%8.4f  %8.4f\n", i*dx, PdfF (n, d, c, i*dx));
   }

   // Return p-value for a one-sided test.
   P = 1.0 - sum;

   return P;

}

double PdfF (int n, int d, double c, double x) {

   return c * pow (x, (n/2.0) - 1.0) * pow (1.0 + (n*x)/d, -(n+d)/2.0);

}