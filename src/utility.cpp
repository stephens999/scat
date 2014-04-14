#include "utility.hpp"

#include <numeric>
#include <cstdlib>
#include <cmath>

using namespace std;


double ranf(){
  return random()/(RAND_MAX+1.0);
}

// generate a random integer according to a user-defined density
int rint2 ( const vector<double> & prob, double psum )
{
    double csum = prob[0];
    double u = ranf();
    
    if(psum == 0.0){ // return a uniform random number if all zeros
      return (int) floor(prob.size() * u);
    }
    else
      {
	if ( psum > 0.0 ) {
	  u *= psum;
	  for (int i = 0; i < prob.size() - 1; ++i) {
            if ( u < csum ) return i;
            csum += prob[i+1];
	  }
	} else {
	  vector<double> cumprob ( prob.size(), 0.0 );
	  // Calculate cdf
	  std::partial_sum ( prob.begin(), prob.end(), cumprob.begin());
	  u *= cumprob[prob.size()-1];
	  for (int i = 0; i < cumprob.size() - 1; ++i) {
            if ( u < cumprob[i] ) return i;
	  }
	}
	return prob.size() - 1;
      }
   
}

//
// Random permutation of 0 to n-1, in perm
//
void rperm(vector<int> & perm,int n)
{
  int i,s,temp,t;
  for(i=0;i<n;i++)
    perm[i]=i;
  for(i=0;i<n;i++)
    {
      t=n-i-1; // t runs from n-1 down to 0
      s=(int) floor((t+1)*ranf()); // s unif on 0 to t
      temp=perm[s]; // swap s and t
      perm[s]=perm[t];
      perm[t]=temp;
    }
}

// normal density
double dnorm(double x)
{
  double PI = 3.141592;
  return (1.0/sqrt(2*PI)) * exp(-0.5*x*x);
}
//
// normal random generator
// mean mu and variance 1/t
// (Ripley (1987), alg 3.17, P82)
//
double rnorm(double mu,double sigma)
{
	double u=0,v=0,x=0,z=0;

	loopstart:
	u=ranf();
	v=0.8578*(2*ranf()-1);
	x=v/u;
	z=0.25*x*x;
	if(z<(1-u)) goto loopend;
	if(z>(0.259/u+0.35)) goto loopstart;
	if(z>(-log(u))) goto loopstart;
	loopend:
	return mu+sigma*x;
}
//
// gamma random generator
// from Ripley, 1987, P230
//
double rgamma(double n,double lambda)
{

	double x=0.0;
	if(n<1)
	{
		const double E=2.71828182;
		const double b=(n+E)/E;
		double p=0.0;
		one: 
		p=b*ranf();
		if(p>1) goto two;
		x=exp(log(p)/n);
		if(x>-log(ranf())) goto one;
		goto three;
		two: 
		x=-log((b-p)/n);
		if (((n-1)*log(x))<log(ranf())) goto one;
		three:;	
	}
	else if(n==1.0)
//
// exponential random variable, from Ripley, 1987, P230
//	
	{
		double a=0.0;
		double u,u0,ustar;
	ten:
		u=ranf();
		u0=u;
	twenty:
		ustar=ranf();
		if(u<ustar) goto thirty;
		u=ranf();
		if(u<ustar) goto twenty;
		a++;
		goto ten;
	thirty:
		return (a+u0)/lambda;
	}
	else
	{
		double static nprev=0.0;
		double static c1=0.0;
		double static c2=0.0;
		double static c3=0.0;
		double static c4=0.0;
		double static c5=0.0;
		double u1;
		double u2;
		if(n!=nprev)
		{
			c1=n-1.0;
			double aa=1.0/c1;
			c2=aa*(n-1/(6*n));
			c3=2*aa;
			c4=c3+2;
			if(n>2.5) c5=1/sqrt(n);
		}
		four:
		u1=ranf();
		u2=ranf();
		if(n<=2.5) goto five;
		u1=u2+c5*(1-1.86*u1);
		if ((u1<=0) || (u1>=1)) goto four;
		five:
		double w=c2*u2/u1;
		if(c3*u1+w+1.0/w < c4) goto six;
		if(c3*log(u1)-log(w)+w >=1) goto four;
		six:
		x=c1*w;		
		nprev=n;
	}	

	return x/lambda;
}

//
// dirichlet random generator
// set b to be ~Dirichlet(a)
//
void rdirichlet(const double * a, const int k, double * b)
{
  int i;
	double sum=0.0;
	for(i=0;i<k;i++)
	{
		b[i]=rgamma(a[i],1);
		sum += b[i];
	}
	for(i=0;i<k;i++)
	{
		b[i] /= sum;
	}
}

//
// dirichlet random generator
// set b to be ~Dirichlet(a)
//
void rdirichlet(const vector<double> & a, const int k, vector<double> & b)
{
  int i;
	double sum=0.0;
	for(i=0;i<k;i++)
	{
		b[i]=rgamma(a[i],1);
		sum += b[i];
	}
	for(i=0;i<k;i++)
	{
		b[i] /= sum;
	}
}

