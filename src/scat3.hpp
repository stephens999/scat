#ifndef SCAT2_HPP
#define SCAT2_HPP

// SCAT version 3.0.3

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cmath>
#include "string.h"
#include <map>

extern "C" void dpotrf_(
	const char &uplo,		// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int &info			// (output)
	);

using namespace std; 
const double PI = 3.141592; 
const string VERSION="3.0.3";

const int FLAGINT = -1;
    
const int MAXSPECIES = 2;
// This line is 120 for microsats and 2 for SNPs
// const int MAXNALLELE = 120; 
const int MARGIN = 7;  // margin around grid

typedef std::vector<double> DoubleVec1d;
typedef std::vector<DoubleVec1d> DoubleVec2d;
typedef std::vector<DoubleVec2d> DoubleVec3d;
typedef std::vector<DoubleVec3d> DoubleVec4d;
typedef std::vector<int> IntVec1d;
typedef std::vector<IntVec1d> IntVec2d;
typedef std::vector<IntVec2d> IntVec3d;
typedef std::vector<IntVec3d> IntVec4d;

class Mapgrid;

const int ValidateAssumptions = 1;

const double EPSILON = 1e-100;

const int DELTALENGTH = 3;

void error_and_exit(const string& msg);

bool CheckThetaValues(const DoubleVec4d& Theta, const DoubleVec4d& ExpTheta,
  const DoubleVec3d& SumExpTheta);

bool CheckPsiValues(const DoubleVec2d& Psi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi);

bool CheckCountValues(const IntVec4d& Count, const IntVec3d& SumCount);

bool CheckConsistency(const DoubleVec4d& Theta, const DoubleVec4d& ExpTheta,
  const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount,
  const DoubleVec2d& Psi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi);

bool CompareLogLiks(const DoubleVec3d& newLL, const DoubleVec3d& oldLL);

double to_degrees(double radianvalue);

bool InRange(double x, double y, const DoubleVec1d& BoundaryX, const DoubleVec1d& BoundaryY, const Mapgrid& mymapgrid);

double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2);

int IsInsideBoundary( double x, double y, const vector<double> & BoundaryX, const vector<double> & BoundaryY);

void InitialiseXY(vector<double> & BoundaryX, vector<double> & BoundaryY, vector<double> & Xcoord, vector<double> & Ycoord, const Mapgrid& mymapgrid);

void ReadInBoundary(ifstream & bfile, vector<double> & BoundaryX, vector<double> & BoundaryY);
	
int GetLocationNumber(vector<int> & RegionsPresent, int r);

int GetLocationNumberAdd(vector<int> & RegionsPresent, int r);

void permute_regions (vector<int> & Region, vector<int> & Perm);

string getline(streambuf * pbuf);

void input_genotype_data( ifstream & input, vector<int>  & Region,
  vector<int> & Species, vector<vector<vector<int> > > & Genotype, vector<int> & NMissing, vector<string> & Id, vector<int> & RegionsPresent,bool useregion);

void OutputLatLongs(ostream & locatefile, double x, double y, double loglik);

void OutputRegionNames(ostream & freqfile, const vector<string> & RegionName, const vector<int> & Perm);

void OutputAcceptRates(ostream & ostr);

void output_positions_data(const vector<string> & RegionName, const vector<int> & Region, const vector<double> & x, const vector<double> & y, const vector<string> & Id);

void input_positions_data( ifstream & input, vector<double> & x, vector<double> & y, vector<string> & RegionName, vector<int> & SubRegion, vector<int> & Region, vector<int> & Perm, vector<int> & RegionsPresent);

void output_genotypes(const vector<vector<vector<int> > > & Genotype, const vector<string> & Id);

void recode_genotypes(vector<vector<vector<int> > > & OriginalGenotype, vector<vector<vector<int> > > & RecodedGenotype, vector<map<int,int> >& Coding, vector<int> & Nallele);

void SubtractFromCount(int ind, IntVec4d& Count, IntVec3d& SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector<int> > > & Genotype);

void AddToCount(int ind, IntVec4d& Count, IntVec3d& SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype);

void count_up_alleles(IntVec4d& Count, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype);

void calc_SumCount(const IntVec4d& Count, IntVec3d& SumCount);

double simple_Prob(const IntVec4d& Count, const IntVec3d& SumCount,vector<vector<vector<int> > > &  Genotype, int ind, int r);

double log_hybrid_Prob(const IntVec4d& Count, const IntVec3d& SumCount,vector<vector<vector< int> > > &  Genotype, int ind, int r, int s);

double log_hybrid_Prob(const DoubleVec4d& Freq,vector<vector<vector<int> > > & Genotype, int ind, int r, int s);

void calc_ExpTheta_and_SumExpTheta(DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& SumExpTheta);

void output_counts(const IntVec4d& Count, vector<int> & Perm);

double FittedCovariance(const vector<double> & Alpha, double d);
  
double Covariance(int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu);

double Correlation(int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu);

double Covariance(int locus, int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu);

double Correlation(int l, int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu);

// compute euclidean distance from positions in radians
//from http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/TSPFAQ.html
double Distance(double x1, double y1, double x2, double y2);

double Distance(int r, int s, vector<double> & Xcoord, vector<double> & Ycoord);

void calc_L(DoubleVec1d& L,const vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord);

double CurrentLogLik(const DoubleVec3d& LogLik, int r);

void calc_LogLik(DoubleVec3d& LogLik, const DoubleVec4d& Theta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount);

double calc_LogLikWithoutPseudoCounts(const DoubleVec4d& Theta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount);

double SumLogLik(const DoubleVec3d& LogLik);

void NormaliseMeanCov(DoubleVec4d& MeanCov, DoubleVec4d& MeanFittedCov, double totaliter);

void CheckSumExpTheta(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta);

void UpdateMeans(const DoubleVec4d& ExpTheta, DoubleVec4d& X, vector<vector<double> > & Pi, const DoubleVec3d& SumExpTheta, DoubleVec4d& MeanFreq, DoubleVec4d& MeanX, DoubleVec4d& MeanX2, vector<vector<double> > & MeanPi, DoubleVec4d& MeanCov, DoubleVec4d& MeanFittedCov, DoubleVec4d& MeanCor, DoubleVec4d& MeanFittedCor, vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu);

void UpdateSubRegionProb(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, double *** SubRegionProb, vector<vector<vector<int> > > & Genotype, vector<int> & SubRegion);

void UpdateLocusMeanProb(int ind, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, DoubleVec4d& LocusMeanProb, vector<vector<vector<int> > > & Genotype);

void UpdateLocusMeanProb(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, DoubleVec4d& LocusMeanProb, vector<vector<vector<int> > > & Genotype);

void NormaliseMeanPi(vector<vector<double> > & MeanPi);

void NormaliseMeanFreq(DoubleVec4d& MeanFreq,double totaliter);

void ComputeLogMeanProb(DoubleVec3d& MeanProb, const DoubleVec4d& LocusMeanProb, double totaliter);

void NormaliseSubRegionProb(DoubleVec3d& SubRegionProb);

void OutputPi(vector<vector<double> > & Pi,vector<string> & RegionName,ostream & ostr, vector<int> & Perm);
  
void OutputLogMeanProb(const DoubleVec3d& MeanProb, ofstream & output, vector<int> & NMissing, vector<int> & Region, vector<string> & Id, vector<string> & RegionName, vector<int> & Perm, int first, int last, vector<int> & RegionsPresent);

void OutputLogMeanProb2(const DoubleVec4d& MeanFreq, vector<vector<vector<int> > > & Genotype, ostream & output, vector<int> & NMissing, vector<int> & Region, vector<int> & Perm);

void OutputParameters(ofstream & paramfile, vector<double> & Alpha, double & Beta, vector<double> & Gamma, vector<double> & Delta, double & Eta, const DoubleVec1d& Lambda, const DoubleVec3d& LogLik);


// ------------------------- The MCMC update routines ---------------------

// The next 3 functions I've just added in response to the interest D+J
// expressed; so pretty much untested (previously Species set to 
// 0 for all individuals)
// 
// update vector of allocation variables (Species) which holds
// which species (= subpopulation) each individual is currently
// assigned to
void update_Species(vector<int> & Species, vector< vector<double> > & Pi, vector<int> & Region, vector<vector<vector<int> > > & Genotype, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta);

void compute_Pi(DoubleVec2d& Pi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi);

void update_Pi(vector< vector<double> > & Pi, vector<int> & Species, vector<int> & Region);

void update_Nu(DoubleVec3d& Nu, vector<double> & Gamma, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount);

void update_Mu(DoubleVec2d& Mu, double Beta, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount);

void update_Lambda(DoubleVec1d& Lambda, double Eta, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Region, vector<int> & Species);

void update_Alpha0(vector<double> & Alpha, const DoubleVec4d& X);

void update_Delta0(vector<double> & Delta, const DoubleVec2d& Y);

void update_Beta(double & Beta, const DoubleVec2d& Mu);

void update_Eta(double & Eta, const DoubleVec1d& Lambda);

bool InElephantRange(double x, double y);

bool InForest(double x,double y);

void update_Location(vector<double> & Alpha, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& SumExpTheta, DoubleVec3d& LogLik, DoubleVec1d& L, const IntVec4d& Count, const IntVec3d& SumCount, vector<double> & Xcoord, vector<double> & Ycoord, vector<int> & Species, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<double> & BoundaryX, vector<double> & BoundaryY, const Mapgrid& mymapgrid);

// update the parameters in the covariance matrix
void update_Alpha(vector<double> & Alpha, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, DoubleVec1d& L, vector<double> & Xcoord, vector<double> & Ycoord);

// update the parameters in the covariance matrix M 
void update_Delta(vector<double> & Delta, const DoubleVec2d& Y, const DoubleVec1d& Lambda, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Species, vector<int> & Region, DoubleVec1d& M, vector<double> & Xcoord, vector<double> & Ycoord);

double divLogLikValue(int r, int k, int l, int j, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, DoubleVec1d& L);

double calcNewdivLogLik(int r, int k, int l, int j, const DoubleVec1d& ExpTheta, const DoubleVec1d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L);

void update_XJoint(vector<double> & Alpha, DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, DoubleVec1d& L);

void update_XSingle(vector<double> & Alpha, DoubleVec4d& X, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, DoubleVec1d& L);

void update_YSingle(vector<double> & Delta, DoubleVec2d& Y, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Region, vector<int> & Species, DoubleVec1d& M);

void DoAllUpdates(DoubleVec4d& X, double & Beta,  DoubleVec1d& Gamma, DoubleVec1d& Alpha, DoubleVec2d& Mu, DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, IntVec4d& Count, IntVec3d& SumCount, DoubleVec1d& L, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<int> & Species, vector<vector<double> > & Pi, vector<double> & Xcoord, vector<double> & Ycoord, DoubleVec2d& Y, double & Eta, vector<double> & Delta, DoubleVec1d& Lambda, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, DoubleVec1d& M, vector<double> & BoundaryX, vector<double> & BoundaryY, const Mapgrid& mymapgrid );

void InitialiseTheta(DoubleVec4d& Theta, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec1d& L);

void Initialise(DoubleVec4d& X, double & Beta,  vector<double> & Gamma, vector<double> & Alpha, DoubleVec2d& Mu, DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec1d& L, vector<double> & Xcoord, vector<double> & Ycoord);

void output_empirical_freqs(vector<string> & RegionName,vector<int> & Perm, vector<map<int,int> >& Coding, const IntVec4d& Count, const IntVec3d& SumCount);
	
void OutputMeanFreq(ofstream & freqfile, vector<string> & RegionName, vector<int> & Perm, vector<map<int,int> >& Coding, const DoubleVec4d& MeanFreq );

void OutputEstimatedFreqs(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, ostream& output, vector<map<int,int> >& Coding, vector<int>& Perm);

std::string ToString(int);
std::string ToString(unsigned long);
std::string ToString(unsigned long long);
std::string ToString(long);
std::string ToString(long long);
std::string ToString(double);

void cholesky_in_place(DoubleVec1d& myL, int nregion);

int main ( int argc, char** argv);
#endif  // SCAT2_HPP

