/*
**     Project1: d)
**     The algorithm for solving the tridiagonal matrix
**     equation is implemented (requiering O(8n) FLOPS) and
**     compared to the method of LU decomposition. The time
**     usage results are written to a txt file as a table.
*/
# include <iostream>
# include <fstream>
# include <iomanip>
# include <time.h>
# include <new>
# include <string>
# include <cstdio>
# include <cstdlib>
# include <cmath>
# include <cstring>
# include "lib.h"

using namespace std;
ofstream ofile;

// Declaration of functions:
void initialize (int *, int *, int *);
double f(double);
double Solution(double);
double eps(double, double);
double time_of_tridiagonal_method(int);
double time_of_LU_method(int);

// Functions:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double f(double x) {return 100*exp(-10*x);}

double eps(double u, double v) {return log10(fabs(v-u));}

void initialize (int *initial_n, int *number_of_steps, int *step_increase)
{
  printf("Read in from screen initial_n, number_of_steps and step_increase \n");
  scanf("%d %d %d",initial_n, number_of_steps, step_increase);
  return;
}

double time_of_LU_method(int n) {
  double *b_twidd = new double[n];
  double *x = new double[n+2];
  double **A;
  double calculation_time;
  double h = 1.0/(n+1.0);

  // Constructing x-array:
  for (int i=0; i<=n+1; i++) {x[i] = i*h;}

  // Constructing right hand side of equation:
  for (int i=0; i<n; i++) {b_twidd[i] = h*h*f(x[i+1]);}

  // Constructing matrix A the brute force way
  // without use of Armadillo here.
  A = new double*[n];
  for (int i=0; i<n; i++) {
    A[i] = new double[n];
  }
  int flag;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
       flag = 0;
       if (i == j) {
        A[i][j] = 2;
        flag = 1;
       }
       if (i == j+1 || j == i+1 ) { 
        A[i][j] = -1;
        flag = 1;
       }
       if (flag == 0) {
        A[i][j] = 0;
       }
    // Could print results to check that A comes out right:
    //cout << " " << A[i][j] << " ";
    }
  //cout << endl;
  }
  
  // LU decomposition of A:
  // allocate space in memory:
  int *indx;
  double d;
  indx = new int[n];

  // Time of calculation starts here:
  clock_t start, finish;
  start = clock();
  
  ludcmp(A, n, indx, &d);   // LU decompose  A[][] 
  lubksb(A, n, indx, b_twidd); // Solve linear equation

  // The solution v is now in b_twidd!
  finish = clock();
  calculation_time = (finish-start)/(double)CLOCKS_PER_SEC;

  // Print to compare results:
  //for (int i = 0; i<n; i++) {
  //  cout << " LU-method: v =  " << b_twidd[i] << endl;
  //}
  
  delete [] x;
  delete [] b_twidd;
  for (int i=0; i<n; i++) {
    delete [] A[i];
  }
  delete [] A;

  return calculation_time;
}

double time_of_tridiagonal_method(int n) {
    // Constants of the problem:
    double calculation_time;
    double h = 1.0/(n+1.0);
    double *x = new double[n+2];
    double *b_twidd = new double[n+1];
    b_twidd[0] = 0;

    // The constituents of the tridiagonal matrix A:
    // Zeroth element not needed, but included to make indexing easy:
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    // Temporal variabel in Gaussian elimination:
    double *diag_temp = new double[n+1];

    // Real solution and approximated one:
    double *u = new double[n+2]; // Analytical solution
    double *v = new double[n+2]; // Numerical solution
    // Including extra points to make the indexing easy:
    u[0] = 0;
    v[0] = 0;

    // Filling up x-array:
    for (int i=0; i<=n+1; i++) {
        // Making x[0] = 0 and x[n+1] = 1:
        x[i] = i*h;
    }

    // Filling up b_twiddle array, i.e. right hand side of equation:
    for (int i=1; i<=n; i++) {
        b_twidd[i] = h*h*f(x[i]);
        u[i] = Solution(x[i]);
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;

    // Algorithm for finding v:
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = b_twidd(i)
    // Row reduction; forward substitution:
    // Time of calculation starts here:
    clock_t start, finish;
    start = clock();

    double b_temp = b[1];
    v[1] = b_twidd[1]/b_temp;
    for (int i=2;i<=n;i++) {
       // Temporary value needed also in next loop:
       diag_temp[i] = c[i-1]/b_temp;
       // Temporary diagonal element:
       b_temp = b[i] - a[i]*diag_temp[i];
       // Updating right hand side of matrix equation:
       v[i] = (b_twidd[i]-v[i-1]*a[i])/b_temp;
    }

    // Row reduction; backward substition:
    for (int i=n-1;i>=1;i--) {
        // Can use diag_temp here since a[i] = c[i] = -1:
        v[i] -= diag_temp[i+1]*v[i+1];
    }

    finish = clock();
    calculation_time = (finish-start)/(double)CLOCKS_PER_SEC;

    // Could print to compare results:
    //for (int i=1;i<n+1;i++) {
    //  cout << " Tridiagonal method: v =  " << v[i] << endl;
    //}

    // Delete arrays, free memory:
    delete [] diag_temp;
    delete [] x;
    delete [] b_twidd;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] v;

  return calculation_time;
}

int main(int argc, char* argv[]) {
  // Declaration of initial variables:
  int initial_n, number_of_steps, step_increase;
  double max_error, calculation_time_used_trid, calculation_time_used_LU;

  // Want to store table of result in a .txt file:
  string outfilename;
  outfilename = "Time_comparison.txt";

  initialize (&initial_n, &number_of_steps, &step_increase);
  
  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "       n:       time (LU):     time (tridiag): " << endl;
  int n = initial_n;
  // Loop over all n producing table of time used:
  for (int i=0; i<number_of_steps; i++) {
      ofile << setw(10) << setprecision(5) << n;

      // Calculate only with LU decompostition if n is not "too big":
      if ( n <= 3000) {
        ofile << setw(15) << setprecision(8) << time_of_LU_method(n);
       }
      if (n > 3000) {
        ofile << setw(15) << setprecision(8) << " - - ";
      }
      ofile << setw(15) << setprecision(8) << time_of_tridiagonal_method(n) << endl;

      n *= step_increase;
  }
  ofile.close();

  return 0;
}
