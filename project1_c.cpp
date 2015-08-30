/*
**     Project1: c)
**     The algorithm for solving the tridiagonal matrix
**     equation is implemented (requiering O(8n) FLOPS).
**     Log of the absolute error is written to file for
**     several different values of n. initial n, number of steps
**     and step_increase is read from screen. Name of output file
**     is read from command line. Results read and plotted in python
**     script project1_c_plot.py
*/
# include  <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>

using namespace std;
ofstream ofile;

// Declaration of functions:
void initialize (int *, int *, int *);
double f(double);
double Solution(double);
double eps(double, double);
double calculate_max_eps(int);

// Functions:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double f(double x) {return 100*exp(-10*x);}

double eps(double u, double v) {return log10(fabs(v-u));}

void initialize (int *initial_n, int *number_of_steps, int *step_increase)
{
  printf("Read in from screen initial_n, number_of_steps and step_increase \n");
  scanf("%d %d %d",initial_n, number_of_steps , step_increase);
  return;
}

double calculate_max_eps(int n) {
    // Constants of the problem:
    double max_eps;
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

    // Filling up x-array, making x[0] = 0 and x[n+1] = 1:
    for (int i=0; i<=n+1; i++) {x[i] = i*h;}

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
        v[i] -= diag_temp[i+1]*v[i+1];
    }
    
    // Finding largest relative error and returning it:
    max_eps = fabs(eps(u[1],v[1]));
    for (int i=1; i<=n; i++) {
      // Recall that eps() calculates log_10 of error. Want this to
      // be as small as possible in aboslute value, meaning
      // the largest possible error:
      if (fabs(eps(u[i],v[i])) < fabs(max_eps)) {
         // found new largest value:
         max_eps = eps(u[i],v[i]);
      }
    }

    // Delete arrays, free memory:
    delete [] diag_temp;
    delete [] x;
    delete [] b_twidd;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] v;

  return max_eps;
}

int main(int argc, char* argv[]) {
  // Declaration of initial variables:
  int initial_n, number_of_steps, step_increase;
  double max_error;
  char *outfilename;

  if( argc <= 1 ){
      cout << "Bad Usage: " << argv[0] <<
          " read also output file on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1];
    }

  initialize (&initial_n, &number_of_steps, &step_increase);
  // More variables to be declared:
  double *h_log10 = new double[number_of_steps];
  double *eps_max_log10 = new double[number_of_steps];
  
  // Write results to txt file:
  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << " n:            log10(h):      log10(max(eps)): " << endl;
  int n = initial_n;
  // Loop over all n producing log_10(h) and log_10(max_eps):
  for (int i=0; i<number_of_steps; i++) {
      h_log10[i] = log10(1.0/(n+1.0));
      eps_max_log10[i] = calculate_max_eps(n);
      ofile << setw(10) << setprecision(5) << n;
      ofile << setw(15) << setprecision(8) << h_log10[i];
      ofile << setw(15) << setprecision(8) << eps_max_log10[i] << endl;

      n *= step_increase;
  }
  ofile.close();
  return 0;
}
