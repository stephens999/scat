#include <string>
#include <map>
#include <vector>

using namespace std;

double ranf();
// generate a random integer according to a user-defined density
int rint2 ( const std::vector<double> & , double psum = -1.0 );
double rgamma(double n,double lambda);
void rdirichlet(const double * a, const int k, double * b);
void rdirichlet(const vector<double> & a, const int k, vector<double> & b);
void rperm(vector<int> & perm,int n);
double rnorm(double,double);
double dnorm(double);
