// scat2.cpp : Defines the entry point for the console application.
//

// SCAT version 3.0.2

#include "scat2.hpp"
#include "utility.hpp"
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
#include "readboundary.hpp"

extern "C" void dpotrf_(
	const char &uplo,		// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int &info			// (output)
	);

using namespace std; 
    
string tracefilename("tracefile.txt");
ofstream TRACEFILE(tracefilename.c_str());

void error_and_exit(const string& msg) {
  cerr << msg << endl;
  exit(-1);
}


bool CheckThetaValues(const DoubleVec4d& Theta, const DoubleVec4d& ExpTheta,
  const DoubleVec3d& SumExpTheta) {
  bool correct = true;
  double epsilon = 0.0001;
  for(int r=0; r<NREGION; r++) {
    for(int k=0; k<NSPECIES; k++) {
      for(int l=0; l<NLOCI; l++) {
        double sumExpTheta = 0.0;
        for(int j=0; j<Nallele[l]; j++) {
          // check expTheta
          double expTheta = exp(Theta[r][k][l][j]);
          if(ExpTheta[r][k][l][j] != expTheta) {
            cerr << "ExpTheta incorrect for " << r << " ";
            cerr << k << " " << l << " " << j << endl;
            correct = false;
          }
          sumExpTheta += expTheta;
        }
        // check SumExpTheta
        // we do not check for exact equality here as apparently this variable
        // accumulates a small rounding error, as reflected in original author's
        // CheckSumExpTheta function.
        if(fabs(sumExpTheta - SumExpTheta[r][k][l]) > epsilon) {
          cerr << "SumExpTheta incorrect for " << r << " ";
          cerr << k << " " << l << endl;
          cerr << "Expected " << SumExpTheta[r][k][l] << " calculated " << sumExpTheta << endl;
          correct = false;
        }
      }
    }
  }
  return correct;
}

bool CheckPsiValues(const DoubleVec2d& Psi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi) {
  bool correct = true;
  for(int r=0; r<NREGION; r++) {
    double sumExpPsi = 0.0;
    for(int k=0; k<NSPECIES; k++) {
      // check expPsi
      double expPsi = exp(Psi[r][k]);
      if(ExpPsi[r][k] != expPsi) {
        cerr << "ExpPsi incorrect for " << r << " " << k << endl;
        cerr << "Expected " << ExpPsi[r][k] << " " << " calculated " << expPsi << endl;
        correct = false;
      }
      sumExpPsi += expPsi;
    }
    if(SumExpPsi[r] != sumExpPsi) {
      cerr << "SumExpPsi incorrect for " << r << endl;
      cerr << "Expected " << SumExpPsi[r] << " calculated " << sumExpPsi << endl;
      correct = false;
    }
  }
  return correct;
}

bool CheckCountValues(const IntVec4d& Count, const IntVec3d& SumCount)
{
  bool correct = true;
  for(int r=0; r<NREGION; r++) {
    for(int k=0; k<NSPECIES; k++) {
      for(int l=0; l<NLOCI; l++) {
        double sumCount = 0.0;
        for(int j=0; j<Nallele[l]; j++) {
	  sumCount += Count[r][k][l][j];
        }
        // Check SumCount
        if(sumCount != SumCount[r][k][l]) {
          cerr << "SumCount incorrect for " << r << " " << k << " " << l << endl;
          cerr << "Expected " << SumCount[r][k][l] << " calculated " << sumCount << endl;
          correct = false;
        }
      }
    }
  }
  return correct;
}


bool CheckConsistency(const DoubleVec4d& Theta, const DoubleVec4d& ExpTheta,
  const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount,
  const DoubleVec2d& Psi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi) {

  bool correct = CheckThetaValues(Theta,ExpTheta,SumExpTheta);
  correct = correct && CheckPsiValues(Psi,ExpPsi,SumExpPsi);
  correct = correct && CheckCountValues(Count, SumCount);
  return correct;
}

bool CompareLogLiks(const DoubleVec3d& newLL, const DoubleVec3d& oldLL) {
  // region species locus allele
  double discrepancy = 0.0;
  double epsilon = 0.0001;
  for(int r = 0; r < NREGION; ++r) {
    for(int k=0; k < NSPECIES; ++k) {
      for(int l=0; l < NLOCI; ++l) {
        discrepancy += fabs(newLL[r][k][l] - oldLL[r][k][l]);
      }
    }
  }
  if (discrepancy > epsilon) {
    cerr << "Likelihood has changed:  difference between old and new " << discrepancy << endl;
    return false;
  }
  return true;
}

double to_degrees(double radianvalue) {
  return (radianvalue * 180.0/PI);
}


bool InRange(double x, double y, const DoubleVec1d& BoundaryX, const DoubleVec1d& BoundaryY, const Mapgrid& mymapgrid) {
  bool isinrange;
  if(READGRID) {
    // REVERSAL here because mymapgrid uses (lat,long) whereas SCAT2 normally uses (long,lat)
    // CONVERSION here because mymapgrid uses degrees and SCAT2 normally uses radians
    isinrange = mymapgrid.in_range(y*180.0/PI,x*180.0/PI);
  } else if (FORESTONLY) {
    isinrange = (InElephantRange(x,y) && InForest(x,y));
  } else if (SAVANNAHONLY) {
    isinrange = (InElephantRange(x,y) && !InForest(x,y));
  } else {
    assert(READBOUNDARY);
    isinrange = IsInsideBoundary(x,y,BoundaryX,BoundaryY);
  }
  return isinrange;
}

  
double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2){
  return( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
}

//  winding number test for a point in a polygon
//  modelled on code from softsurfer.com, by Dan Sunday
//      Input:   x,y = a point,
//               BoundaryX and BoundaryY = points of a polygon with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if (x,y) is outside polygon)
int IsInsideBoundary( double x, double y, const DoubleVec1d& BoundaryX, const DoubleVec1d& BoundaryY)
{
  if(BoundaryX.size() == 0) // if no boundary, just return 1
    return 1;

  int    wn = 0;    // the winding number counter
  
  // loop through all edges of the polygon
  for (int i=0; i < (BoundaryX.size()-1); i++) {   // edge from V[i] to V[i+1]
    if (BoundaryY[i] <= y) {         // start y <= P.y
      if (BoundaryY[i+1] > y)      // an upward crossing
	if (isLeft( BoundaryX[i], BoundaryY[i], BoundaryX[i+1], BoundaryY[i+1], x, y) > 0)  // P left of edge
	  ++wn;            // have a valid up intersect
    }
    else {                       // start y > P.y (no test needed)
      if (BoundaryY[i+1] <= y)     // a downward crossing
	if (isLeft( BoundaryX[i], BoundaryY[i], BoundaryX[i+1] , BoundaryY[i+1], x, y) < 0)  // P right of edge
	  --wn;            // have a valid down intersect
    }
  }
  return wn;
}
//===================================================================

void InitialiseXY(vector<double> & BoundaryX, vector<double> & BoundaryY, vector<double> & Xcoord, vector<double> & Ycoord, const Mapgrid& mymapgrid)
{
  if(FORESTONLY){
    Xcoord[NREGION-1] = 0.2;
    Ycoord[NREGION-1] = 0.1;
  } else if(SAVANNAHONLY){
    Xcoord[NREGION-1] = 0.52;
    Ycoord[NREGION-1] = -0.17;
  } else if(CHEAT){
    Xcoord[NREGION-1] = Xcoord[TRUEREGION]+rnorm(0,0.001); // these lines "cheat" by starting with 
    Ycoord[NREGION-1] = Ycoord[TRUEREGION]+rnorm(0,0.001); // initial guess close to true location
  } else {
    double Xcenter = 0;
    double Ycenter = 0;
    for(int r = 0; r<(NREGION-1); r++){
      Xcenter += Xcoord[r];
      Ycenter += Ycoord[r];
    }
    Xcenter /= (NREGION-1);
    Ycenter /= (NREGION-1);
    int ntries = 0;
    do{
      if (ntries == 100) {
        cerr << "SCAT InitialiseXY(), tried 100 times and fails" << endl;
        exit(-1);
      }
      double lambda = ranf();
      int r = (int) (ranf() * (NREGION-1));
      Xcoord[NREGION-1] = lambda * Xcenter + (1-lambda) * Xcoord[r];
      Ycoord[NREGION-1] = lambda * Ycenter + (1-lambda) * Ycoord[r];
      ntries++;
    } while(InRange(Xcoord[NREGION-1],Ycoord[NREGION-1],BoundaryX,BoundaryY,mymapgrid)==0);
  }
}


void ReadInBoundary(ifstream & bfile, vector<double> & BoundaryX, vector<double> & BoundaryY)
{
   double x,y;
   do{
      bfile >> y;
      bfile >> x;
      if (ECHOINPUTS)
        cout << y << "," << x << endl;
      x = PI*x/180;
      y = PI*y/180;
      BoundaryX.push_back(x);
      BoundaryY.push_back(y);
   } while(BoundaryX.size()<2 || x!=BoundaryX[0] || y!=BoundaryY[0]);

}
	
int GetLocationNumber(vector<int> & RegionsPresent, int r)
{
  if(r <0)
    return r;
  
  int s;
  for(s = 0; s<RegionsPresent.size(); s++){
    if(RegionsPresent[s] == r)
      break;
  }
  return s;
}

int GetLocationNumberAdd(vector<int> & RegionsPresent, int r)
{
	  if(r <0)
		      return r;

	    int s;
	      for(s = 0; s<RegionsPresent.size(); s++){
		          if(RegionsPresent[s] == r)
				        break;
			    }
	        if(s == RegionsPresent.size())
			    RegionsPresent.push_back(r);

		  return s;
}



void permute_regions (vector<int> & Region, vector<int> & Perm){
	for(int r = 0; r< Region.size(); r++){
		if(Region[r]>=0){
			Region[r] = Perm[Region[r]];
		}
	}
}


string getline(streambuf * pbuf)
{
    char ch;
    string str;
    size_t pos;
    while((ch = pbuf->sgetc()) != EOF)
    {
        if(ch != '\n' && ch != '\r')
        {
            str.push_back(ch);
            ch = pbuf->snextc();
        }
        else {
            pbuf->sbumpc();  //chomp;
            pos = str.find_first_not_of(";, \t", 0);
            if(str.empty() || pos  == string::npos || str.at(pos) == '#')
            {
                str.clear();
                continue;
            }
            else
                break;
        }
    }                                                                  
    return str;
}   //this getline use ;, \t as delimit, and ignore the lines either full of delimit or starting with #. 


void input_genotype_data( ifstream & input, vector<int>  & Region,
  vector<int> & Species, vector<vector<vector<int> > > & Genotype, vector<int> & NMissing, vector<string> & Id, vector<int> & RegionsPresent,bool useregion)
{
  string delimit(" \t");
  streambuf * pbuf = input.rdbuf();
  int badallele = -1000;
  
  string line;
  line.assign(getline(pbuf));
  // cout << "READING GENOTYPE DATA" << endl;
  cout << "Reading Genotype Data" << endl;
  while(line.size()>0){
	for(int chrom = 0; chrom<2; chrom++){
  		int beg = line.find_first_not_of(delimit, 0);
    	int end = line.find_first_of(delimit, beg);
    	if(end == beg) break;
    	string sv(line, beg, end-beg);
    	if(chrom==0){
  	  Id.push_back(sv);
	  NMissing.push_back(0);
	  vector<vector<int> > blank(2,vector<int>(NLOCI,badallele));
	  Genotype.push_back(blank);
        } else {
          if (sv != Id.back()) {
            string msg = "Second entry for " + Id.back() + "not found";
            error_and_exit(msg);
          }
        }
		
	beg = line.find_first_not_of(delimit, end);
	end = line.find_first_of(delimit, beg);
	sv = line.substr(beg, end-beg);
	int r = atoi(sv.data());	
        if(!useregion)
			r= -1;
		
		r = GetLocationNumberAdd(RegionsPresent, r);
     
    	if(r <(-1)) { 
			cerr << "Error in numbering of locations in genotype file" << endl;
			cerr << "Individal ID = " << Id[Id.size()-1] << endl;
			exit(1);
     	}
    
        if(chrom==1) {
          if(r>=0) {
            if(r != Region.back()) {
              string msg = "Two haplotypes of the same individual in different regions";
              error_and_exit(msg);
            }
          }
        }
     	if(chrom==0){
			if(r>=0)
	  			Region.push_back(r);
			else{
	  			if(NSPECIES > 1){
	    			cerr << "Error: mustn't have individuals of unknown location with more than 1 species" << endl;
	    		exit(1);
	  			} else {
	    			Region.push_back(r);
	  			}		
			}
	
			Species.push_back((int) (NSPECIES * (1-ranf())));
		}
		
		int dummy;
    	for(int skip = 0; skip<SKIPCOL; skip++){
			beg = line.find_first_not_of(delimit, end);
			end = line.find_first_of(delimit, beg);
	  	}

	int individual = Genotype.size()-1;
      	for(int j=0; j<NLOCI; j++){
      		int allele;
		beg = line.find_first_not_of(delimit, end);
  		end = line.find_first_of(delimit, beg);
  		string sv(line, beg, end-beg);
  		allele = atoi(sv.data());  			  
	 
  		Genotype[individual][chrom][j] = allele;
  		if((allele<0) && (chrom ==0))
     		NMissing[individual]++;
	  	}
        // Check that all NLOCI alleles have been found; if we detect
        // 'badallele' they have not.
        for(int j=0; j<NLOCI; ++j) {
          if (Genotype[individual][chrom][j] == badallele) {
             string msg = "Too few markers found for individual " + to_string(individual);
             error_and_exit(msg);
          }
        }
        line.assign(getline(pbuf));
	}
  }
}

void OutputLatLongs(ostream & locatefile, double x, double y, double loglik){
  // convert x and y into lat and long
  locatefile << 180 * y/PI << " " << 180 * x/PI <<  " " << loglik << endl; 
}

void OutputRegionNames(ostream & freqfile, const vector<string> & RegionName, const vector<int> & Perm)
{
  for(int r=0;r<NREGION;r++){     
    freqfile << std::fixed << setw(9-RegionName[Perm[r]].size()) << RegionName[Perm[r]] << " ";
  }
  freqfile << endl;
}


void OutputAcceptRates(ostream & ostr)
{
  for(int alphaparam =1; alphaparam < ALPHALENGTH; alphaparam++){
    if(ALPHAATTEMPT[alphaparam] >0)
      ostr << (1.0*ALPHAACCEPT[alphaparam])/ALPHAATTEMPT[alphaparam] << " ";
  }
  if(XATTEMPT >0)
    ostr << " " << (1.0*XACCEPT)/XATTEMPT;
  if(MUATTEMPT > 0)
    ostr << " " << (1.0*MUACCEPT)/MUATTEMPT;
  if(NUATTEMPT>0)
    ostr << " " << (1.0*NUACCEPT)/NUATTEMPT;
  
  for(int deltaparam =0; deltaparam < DELTALENGTH; deltaparam++){
    if(DELTAATTEMPT[deltaparam] > 0)
      ostr << (1.0*DELTAACCEPT[deltaparam])/DELTAATTEMPT[deltaparam] << " ";
  }
  if(YATTEMPT>0)
    ostr << " " <<  (1.0*YACCEPT)/YATTEMPT;
  if(LAMBDAATTEMPT>0)
    ostr << " " << (1.0*LAMBDAACCEPT)/LAMBDAATTEMPT;
  ostr << endl;
}    

void output_positions_data(const vector<string> & RegionName, const vector<int> & Region, const vector<double> & x, const vector<double> & y, const vector<string> & Id){

	for(int i=0; i<Region.size(); i++){
		if(Region[i]>=0)
		cout << Id[i] << " : " << RegionName[Region[i]] << "," << x[Region[i]] << "," << y[Region[i]] << endl;
		else
			cout << Id[i] << " : " << Region[i] << endl;
	}
}


void input_positions_data( ifstream & input, vector<double> & x, vector<double> & y, vector<string> & RegionName, vector<int> & SubRegion, vector<int> & Region, vector<int> & Perm, vector<int> & RegionsPresent)
{
  int region,s;
  string regionname;
  int subregion=0;
  double tempx,tempy;
  cout << "Inputting location information" << endl;
  int r=0;
  while(r< RegionsPresent.size()){
    input >> regionname;    
    input >> region;
    if (ECHOINPUTS) {
      cout << "Location " << region << ", Name " << regionname << endl;
    }
    
    region = GetLocationNumber(RegionsPresent, region);
    
    if(region <(-1) ){ 
      cerr << "Error in numbering of locations in location file" << endl;
      cerr << "Line of file = " << (r+1) << endl;
      exit(1);
    }
    
    if(USESUBREGION == 1)
      input >> subregion;

    if(region < RegionsPresent.size()){
       RegionName[Perm[region]] = regionname;
       SubRegion[Perm[region]] = subregion;
       if (ECHOINPUTS) {
         cout << "Region " << region << ", Name " << regionname << endl;
       }
       input >> tempy;
       input >> tempx;

       // convert x and y into radians (x and y should be decimal lat and long)
       x[Perm[region]] = PI * tempx / 180.0; 
       y[Perm[region]] = PI * tempy / 180.0;
       r++;
    } else {
      input >> tempy;
      input >> tempx;
    }
  } 
}

void output_genotypes(const vector<vector<vector<int> > > & Genotype, const vector<string> & Id)
{
  for(int i=0;i<NIND;i++){
    for(int chrom=0; chrom<2; chrom++){
      cout << "Individual " << Id[i]  <<  " : ";
      for(int locus=0; locus<NLOCI; locus++){
	cout << Genotype[i][chrom][locus] << " ";
      }
      cout << endl;
    }   
  }
}


void recode_genotypes(vector<vector<vector<int> > > & OriginalGenotype, vector<vector<vector<int> > > & RecodedGenotype, IntVec2d& Coding, vector<int> & Nallele)
{
	RecodedGenotype = OriginalGenotype; 
  for(int locus=0; locus<NLOCI; locus++){
    vector<int> AllelesPresent;
    for(int ind =0; ind<NIND; ind++){
      for(int chrom=0; chrom<2; chrom++){
	int allele = OriginalGenotype[ind][chrom][locus];
	if(allele>0){
	  if(find(AllelesPresent.begin(),AllelesPresent.end(),allele)==AllelesPresent.end()) // if allele is not in current list of AllelesPresent, add it
	    AllelesPresent.push_back(allele);	
	}
	else
	  RecodedGenotype[ind][chrom][locus] = allele;
      }
    }
    sort(AllelesPresent.begin(),AllelesPresent.end());
    Nallele.push_back(AllelesPresent.size());
	
	if( AllelesPresent.size() > MAXNALLELE ) {
		cerr << "MAXNALLELE is insufficient: " << AllelesPresent.size()  << " needed." << endl;
		exit(1);
	}

    for(int allele =0; allele<AllelesPresent.size(); allele++){
      for(int ind =0; ind<NIND; ind++){
	for(int chrom=0; chrom<2; chrom++){
	  if(OriginalGenotype[ind][chrom][locus]==AllelesPresent[allele]){
	    RecodedGenotype[ind][chrom][locus] = allele;
	    Coding[locus][allele] = AllelesPresent[allele];
	  }
	}
      }
    }
  }
}

void SubtractFromCount(int ind, IntVec4d& Count, IntVec3d& SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector<int> > > & Genotype)
{
  if(Region[ind]>=0){
    for(int chrom=0; chrom<2; chrom++){
      for(int locus=0; locus<NLOCI; locus++){
	if(Genotype[ind][chrom][locus]>=0){
	  Count[Region[ind]][Species[ind]][locus][Genotype[ind][chrom][locus]]--;
	  SumCount[Region[ind]][Species[ind]][locus]--;
	}
      }
    }
  }
}


void AddToCount(int ind, IntVec4d& Count, IntVec3d& SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype)
{
  if(Region[ind]>=0){
    for(int chrom=0; chrom<2; chrom++){
      for(int locus=0; locus<NLOCI; locus++){
	if(Genotype[ind][chrom][locus]>=0){
	  Count[Region[ind]][Species[ind]][locus][Genotype[ind][chrom][locus]]++;
	  SumCount[Region[ind]][Species[ind]][locus]++;
	}
      }
    }
  }
}


void count_up_alleles(IntVec4d& Count, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype)
{
  for(int r=0; r<NREGION; r++){
    for(int k=0; k< NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	for(int j=0; j<Nallele[l]; j++){
	  Count[r][k][l][j] = PSEUDOCOUNT; // start by adding 1 to every allele as a "pseudo-count" 
	  //to avoid possible problems with non-observed alleles having very low freq estimates
	}
      }
    }
  }

  for(int i=0;i<NIND;i++){
    for(int chrom=0; chrom<2; chrom++){
      for(int locus=0; locus<NLOCI; locus++){
	if(Region[i]>=0){
	  if(Genotype[i][chrom][locus]>=0){
	    Count[Region[i]][Species[i]][locus][Genotype[i][chrom][locus]]++;
	  }
	}
      }
    }
  }
}

void calc_SumCount(const IntVec4d& Count, IntVec3d& SumCount)
{
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	SumCount[r][k][l] = 0;
	for(int j=0;j<Nallele[l];j++){
	  SumCount[r][k][l]+=Count[r][k][l][j];
	}
      }
    }
  }
}

// compute "whichrun" type probability for
// individual
// WARNING NEVER CALLED
double simple_Prob(const IntVec4d& Count, const IntVec3d& SumCount,vector<vector<vector<int> > > &  Genotype, int ind, int r)
{
  int k=0;
  double prob = 1;
  for(int l=0;l<NLOCI; l++){
    for(int chr =0; chr <2; chr++){
      if(Genotype[ind][chr][l]>=0){
	 double p=0;
	 if(SumCount[r][k][l]>0)
	   p = (Count[r][k][l][Genotype[ind][chr][l]]*1.0+1.0)/(SumCount[r][k][l]+MAXNALLELE);
	 if(p==0)
	   p = 1.0/ (SumCount[r][k][l]+1);
	 cout << p << endl;
	 prob *= p;
      }	
    }	
  }
  return prob;  
}


// compute prob of individual's genotype as hybrid of region r and s
double log_hybrid_Prob(const IntVec4d& Count, const IntVec3d& SumCount,vector<vector<vector< int> > > &  Genotype, int ind, int r, int s)
{
  int k=0;

  double lprob = 0;
  for(int l=0;l<NLOCI; l++){
    double llocusprob = 0;
    if(Genotype[ind][0][l]>=0 && Genotype[ind][1][l]>=0){
      double pr0, pr1, ps0, ps1;
      pr0 = (1-DELTA) * (Count[r][k][l][Genotype[ind][0][l]] + 1.0)/(SumCount[r][k][l]+Nallele[l]) + DELTA/Nallele[l];
      pr1 = (1-DELTA) * (Count[r][k][l][Genotype[ind][1][l]] + 1.0)/(SumCount[r][k][l]+Nallele[l]) + DELTA/Nallele[l];
      ps0 = (1-DELTA) * (Count[s][k][l][Genotype[ind][0][l]] + 1.0)/(SumCount[s][k][l]+Nallele[l]) + DELTA/Nallele[l];
      ps1 = (1-DELTA) * (Count[s][k][l][Genotype[ind][1][l]] + 1.0)/(SumCount[s][k][l]+Nallele[l]) + DELTA/Nallele[l];
      if(Genotype[ind][0][l] == Genotype[ind][1][l]) // if homozygote
	llocusprob += log( NULLPROB * 0.5 * (pr0 + ps0) + (1-NULLPROB) * (pr0 * ps0));
      else{
	double mult = 1;
	llocusprob += log((1-NULLPROB) * mult * (pr0 * ps1 + pr1 * ps0));
      }
    }
    lprob += llocusprob;
  } 
  return lprob;

}


double log_hybrid_Prob(const DoubleVec4d& Freq,vector<vector<vector<int> > > & Genotype, int ind, int r, int s)
{
  int k=0;

  double lprob = 0;
  for(int l=0;l<NLOCI; l++){
    double llocusprob = 0;
    if(Genotype[ind][0][l]>=0){
      double pr0, pr1, ps0, ps1;
      pr0 = (1-DELTA) * Freq[r][k][l][Genotype[ind][0][l]] + DELTA/Nallele[l];
      pr1 = (1-DELTA) * Freq[r][k][l][Genotype[ind][1][l]] + DELTA/Nallele[l];
      ps0 = (1-DELTA) * Freq[s][k][l][Genotype[ind][0][l]] + DELTA/Nallele[l];
      ps1 = (1-DELTA) * Freq[s][k][l][Genotype[ind][1][l]] + DELTA/Nallele[l];
      if(Genotype[ind][0][l] == Genotype[ind][1][l]) // if homozygote
	llocusprob += log( NULLPROB * 0.5 * (pr0 + ps0) + (1-NULLPROB) * (pr0 * ps0));
      else{
	double mult = 1;
	llocusprob += log((1-NULLPROB) * mult * (pr0 * ps1 + pr1 * ps0));
      }
    }
    lprob += llocusprob;
  } 
  return lprob;

}


void calc_ExpTheta_and_SumExpTheta(DoubleVec4d& Theta, DoubleVec4d& ExpTheta, 
  DoubleVec3d& SumExpTheta)
{
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	SumExpTheta[r][k][l] = 0;
	for(int j=0;j<Nallele[l];j++){
	  ExpTheta[r][k][l][j] = exp(Theta[r][k][l][j]);
	  SumExpTheta[r][k][l] += ExpTheta[r][k][l][j];
	}
      }
    }
  }
}

void output_counts(const IntVec4d& Count, vector<int> & Perm)
{
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	for(int j=0;j<Nallele[l];j++){
	  cout << Perm[r] << "," << k  << "," << l  << "," << j  << ":" << Count[r][k][l][j] << endl;
	}
      }
    }
  }
}

double FittedCovariance(vector<double> & Alpha, double d){
  
  if(d>0){
    if(USESPATIAL)
      return exp( - exp(Alpha[2]* log(d / Alpha[1])) );
    else
      return 0;
  }
  else if(INCLUDENUGGET)
    return 1.0+Alpha[3]; // Alpha[3] is "nugget" effect
  
  return 1.0;
}


// (sample) covariance of thetas in region r0 and r1; species k0 and k1
double Covariance(int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu){

  double Er0k0 = 0; //Expectation of Theta-mu-Nu in r0,k0
  double Er1k1 = 0;
  double Er0k0r1k1 = 0; // Cross term
  int total=0;

  for(int l=0; l<NLOCI; l++){
    for(int j=0;j<Nallele[l];j++){
      Er0k0 += Theta[r0][k0][l][j]-Mu[l][j]-Nu[k0][l][j];
      Er1k1 +=  Theta[r1][k1][l][j]-Mu[l][j]-Nu[k1][l][j];
      Er0k0r1k1 += (Theta[r0][k0][l][j]-Mu[l][j]-Nu[k0][l][j])*(Theta[r1][k1][l][j]-Mu[l][j]-Nu[k1][l][j]);
      total += 1;
    }
  }

  Er0k0 /= total;
  Er1k1 /= total;
  Er0k0r1k1 /= total;

  return (Er0k0r1k1 - Er0k0*Er1k1);

}

// as above, but correlation
double Correlation(int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu){

  double Er0k0 = 0; //Expectation of Theta-Mu-Nu in r0,k0
  double Er1k1 = 0;
  double E2r0k0 = 0; // Expectation of (Theta-Mu-Nu)^2
  double E2r1k1 = 0;
  double Er0k0r1k1 = 0; // Cross term
  int total=0;

  for(int l=0; l<NLOCI; l++){
    for(int j=0;j<Nallele[l];j++){
      Er0k0 += Theta[r0][k0][l][j] - Mu[l][j] - Nu[k0][l][j];     
      E2r0k0 += (Theta[r0][k0][l][j]-Mu[l][j]-Nu[k0][l][j])*(Theta[r0][k0][l][j]-Mu[l][j]-Nu[k0][l][j]);
      Er1k1 +=  Theta[r1][k1][l][j]-Mu[l][j] -Nu[k1][l][j];
      E2r1k1 += (Theta[r1][k1][l][j]-Mu[l][j] -Nu[k1][l][j] )*(Theta[r1][k1][l][j]-Mu[l][j] -Nu[k1][l][j] );
      Er0k0r1k1 += (Theta[r0][k0][l][j]-Mu[l][j] -Nu[k0][l][j] )*(Theta[r1][k1][l][j]-Mu[l][j] -Nu[k1][l][j]);
      total += 1;
    }
  }

  Er0k0 /= total;
  E2r0k0 /= total;
  Er1k1 /= total;
  E2r1k1 /=total;
  Er0k0r1k1 /= total;

  return (Er0k0r1k1 - Er0k0*Er1k1)/sqrt((E2r0k0-Er0k0*Er0k0)*(E2r1k1-Er1k1*Er1k1));

}

// (sample) covariance of thetas at locus "locus" in region r0 and r1; species k0 and k1
double Covariance(int locus, int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu){

  double Er0k0 = 0; //Expectation of Theta-mu in r0,k0
  double Er1k1 = 0;
  double Er0k0r1k1 = 0; // Cross term
  int total=0;

  for(int j=0;j<Nallele[locus];j++){
    Er0k0 += Theta[r0][k0][locus][j]-Mu[locus][j] - Nu[k0][locus][j] ;
    Er1k1 +=  Theta[r1][k1][locus][j]-Mu[locus][j] - Nu[k1][locus][j] ;
    Er0k0r1k1 += (Theta[r0][k0][locus][j]-Mu[locus][j]- Nu[k0][locus][j] )*(Theta[r1][k1][locus][j]-Mu[locus][j]- Nu[k1][locus][j] );
    total += 1;
  }
  
  Er0k0 /= total;
  Er1k1 /= total;
  Er0k0r1k1 /= total;
  return (Er0k0r1k1 - Er0k0*Er1k1);

}

// as above, but correlation
double Correlation(int l, int r0, int k0, int r1, int k1, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu){

  double Er0k0 = 0; //Expectation of Theta-Mu in r0,k0
  double Er1k1 = 0;
  double E2r0k0 = 0; // Expectation of (Theta-Mu)^2
  double E2r1k1 = 0;
  double Er0k0r1k1 = 0; // Cross term
  int total=0;

  for(int j=0;j<Nallele[l];j++){
    Er0k0 += Theta[r0][k0][l][j]-Mu[l][j]- Nu[k0][l][j];
    E2r0k0 += (Theta[r0][k0][l][j]-Mu[l][j]- Nu[k0][l][j] )*(Theta[r0][k0][l][j]-Mu[l][j]- Nu[k0][l][j]);
    Er1k1 +=  Theta[r1][k1][l][j]-Mu[l][j]- Nu[k1][l][j];
    E2r1k1 += (Theta[r1][k1][l][j]-Mu[l][j]- Nu[k1][l][j])*(Theta[r1][k1][l][j]-Mu[l][j]- Nu[k1][l][j]);
    Er0k0r1k1 += (Theta[r0][k0][l][j]-Mu[l][j]- Nu[k0][l][j])*(Theta[r1][k1][l][j]-Mu[l][j]- Nu[k1][l][j]);
    total += 1;
  }

  Er0k0 /= total;
  E2r0k0 /= total;
  Er1k1 /= total;
  E2r1k1 /=total;
  Er0k0r1k1 /= total;

  return (Er0k0r1k1 - Er0k0*Er1k1)/sqrt((E2r0k0-Er0k0*Er0k0)*(E2r1k1-Er1k1*Er1k1));

}

// compute euclidean distance from positions in radians
//from http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/TSPFAQ.html
double Distance(double x1, double y1, double x2, double y2){
  double RRR = 6378.388; 
  double q1 = cos( x1 - x2 ); 
  double q2 = cos( y1 - y2 ); 
  double q3 = cos( y1 + y2 ); 
  return ( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) ); 
}


// distance between region r and region s
double Distance(int r, int s, vector<double> & Xcoord, vector<double> & Ycoord){
  return Distance(Xcoord[r], Ycoord[r], Xcoord[s], Ycoord[s]);
}


// compute the matrix L (uses the Lapack routine dpotrf)
void calc_L(double * L,vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord)
{
  for(int r=0; r<NREGION; r++){
    for(int s=0; s <NREGION; s++){
      L[r+ s*NREGION] = FittedCovariance(Alpha,Distance(r,s,Xcoord,Ycoord));
    }
  }

  int INFO = 0;
  char UPLO = 'L';
  dpotrf_(UPLO,NREGION,L,NREGION,INFO);

  // INFO returns something about whether the computation was successful
  if(INFO>0) {
    cerr << "Warning: INFO=" << INFO << endl;
    cerr << "This probably means likelihood calculation has failed due to variable " << INFO << endl;
    cerr << "If the message continues to appear after burn-in the results cannot be trusted." << endl;
  }
  if(INFO<0) {
    cerr << "Warning: INFO=" << INFO << endl;
    cerr << "This probably means an illegal value in variable " << INFO << endl;
    cerr << "If the message continues to appear after burn-in the results cannot be trusted." << endl;
  }
}
 

double CurrentLogLik(const DoubleVec3d& LogLik, int r){
  double sum = 0.0;
  for(int l=0; l<NLOCI; l++){
	  sum += LogLik[r][0][l];
  }
  return sum;


}
// compute the Loglikelihood
// LogLik[r][k][l] holds the loglikelihood for individuals in region r, 
// species k, at locus l
void calc_LogLik(DoubleVec3d& LogLik, const DoubleVec4d& Theta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount){
  double loglik = 0;
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	LogLik[r][k][l] = 0;
	LogLik[r][k][l] -= SumCount[r][k][l] * log(SumExpTheta[r][k][l]);
	for(int j=0; j<Nallele[l]; j++){
	  LogLik[r][k][l] += Count[r][k][l][j] * Theta[r][k][l][j];	  
	}
      }
    }
  }
}


//used in testing: compute the Loglikelihood but with the pseudocounts, if any, removed
double calc_LogLikWithoutPseudoCounts(const DoubleVec4d& Theta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount){
  double loglik = 0;
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	loglik -= (SumCount[r][k][l]-PSEUDOCOUNT * Nallele[l]) * log(SumExpTheta[r][k][l]);
	for(int j=0; j<Nallele[l]; j++){
	  loglik += (Count[r][k][l][j]-PSEUDOCOUNT) * Theta[r][k][l][j];	  
	}
      }
    }
  }
  return loglik;
}



// Compute overall loglikelihood by summing the array LogLik
double SumLogLik(const DoubleVec3d& LogLik){
  double loglik = 0;
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	loglik += LogLik[r][k][l];
      }
    }
  }
  return loglik;
}

void NormaliseMeanCov(DoubleVec4d& MeanCov, DoubleVec4d& MeanFittedCov, double totaliter)
{
  for(int r=0;r<NREGION;r++){
    for(int k=0; k<NSPECIES; k++){
      for(int r1 = 0; r1<NREGION; r1++){
	for(int k1 = 0; k1 < NSPECIES; k1++){
	  MeanCov[r][k][r1][k1] /= totaliter;
	  MeanFittedCov[r][k][r1][k1] /= totaliter;
	}
      }
    }
  }
}


void CheckSumExpTheta(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta)
{
  for(int r=0; r<NREGION; r++){
	 for(int k = 0; k<NSPECIES; k++){
		for(int l=0;l<NLOCI; l++){
		   double sum = 0;
           for(int allele=0; allele<Nallele[l]; allele++){
		   	  sum += ExpTheta[r][k][l][allele]; 				
           }
		   
		   cout << r << "," << k << "," << l << ", SumExpTheta = " << SumExpTheta[r][k][l] << ", Sum = " << sum << endl;
	 	   if(abs(SumExpTheta[r][k][l]- sum) > 0.1)
			   	exit(1);
		}
	 }
  }
}
//
// keep a running count of the posterior means of some of the quantities
//
void UpdateMeans(const DoubleVec4d& ExpTheta, DoubleVec4d& X, vector<vector<double> > & Pi, const DoubleVec3d& SumExpTheta, DoubleVec4d& MeanFreq, DoubleVec4d& MeanX, DoubleVec4d& MeanX2, vector<vector<double> > & MeanPi, DoubleVec4d& MeanCov, DoubleVec4d& MeanFittedCov, DoubleVec4d& MeanCor, DoubleVec4d& MeanFittedCor, vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord, const DoubleVec4d& Theta, const DoubleVec2d& Mu, const DoubleVec3d& Nu)
{
  for(int r=0;r<NREGION;r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0;l<NLOCI; l++){
	for(int allele=0; allele<Nallele[l]; allele++){
      
	  MeanFreq[r][k][l][allele] += ExpTheta[r][k][l][allele]/SumExpTheta[r][k][l];
	  MeanX[r][k][l][allele] += X[r][k][l][allele];
	  MeanX2[r][k][l][allele] += X[r][k][l][allele] * X[r][k][l][allele];
	  
	}
      }
      MeanPi[r][k] += Pi[r][k];

      for(int r1 = 0; r1<NREGION; r1++){
	for(int k1 = 0; k1 < NSPECIES; k1++){
	  MeanCov[r][k][r1][k1] += Covariance(r,k,r1,k1,Theta,Mu,Nu);
	  MeanFittedCov[r][k][r1][k1] += FittedCovariance(Alpha,Distance(r,r1,Xcoord,Ycoord))/Alpha[0]; // should really have separate Alpha for each k!?
	  MeanCor[r][k][r1][k1] += Covariance(r,k,r1,k1,Theta,Mu,Nu) * Alpha[0];
	  MeanFittedCor[r][k][r1][k1] += FittedCovariance(Alpha,Distance(r,r1,Xcoord,Ycoord)); // should really have separate Alpha for each k!?
	}
      }
    }
  }
}

// next few functions are updates for "assignment" tests on each individual
// (trying to assign each individual to one or other region/species on the basis
// of its genotype)
// Compute prob of individual ind's genotype, for each subregion, species combination
void UpdateSubRegionProb(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, double *** SubRegionProb, vector<vector<vector<int> > > & Genotype, vector<int> & SubRegion)
{
  for(int ind = 0; ind < NIND; ind++){
    for(int r=0;r<NREGION;r++){
      for(int k=0; k<NSPECIES; k++){
	double prob = 1;
	for(int l=0;l<NLOCI; l++){
	  for(int chr =0; chr <2; chr++){
	    if(Genotype[ind][chr][l]>=0){
	      prob *= ExpTheta[r][k][l][Genotype[ind][chr][l]]/SumExpTheta[r][k][l];
	    }	
  	  }	
	}
	SubRegionProb[ind][SubRegion[r]][k] += prob;
      }
    }
  }
}



// Compute prob of individual ind's genotype, at each locus, for each region, species combination.
void UpdateLocusMeanProb(int ind, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, DoubleVec4d& LocusMeanProb, vector<vector<vector<int> > > & Genotype)
{ 
  //  cout << "Ind = " << ind << endl;
  double llocusprob;

  for(int r=0;r<NREGION;r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0;l<NLOCI; l++){
	llocusprob = 0;
	int allele1 = Genotype[ind][0][l];
	int allele2 = Genotype[ind][1][l];
	double p1 = 1.0;
	double p2 = 1.0;
	if(allele1>=0){
	  p1 = (1-DELTA) * ExpTheta[r][k][l][allele1]/SumExpTheta[r][k][l] + DELTA/Nallele[l];
	}
	if(allele2>=0){
	  p2 = (1-DELTA) * ExpTheta[r][k][l][allele2]/SumExpTheta[r][k][l] + DELTA/Nallele[l];
	}
	if(allele1==allele2){
	  llocusprob += log( NULLPROB*p1 + (1-NULLPROB) * p1 * p1);
	} else {
	  if(allele1>=0){
	    llocusprob += log ( (1-NULLPROB) * p1 * p2);
	  }
	}
	LocusMeanProb[ind][r][k][l] += exp(llocusprob);
      }
    }
  }
}

// compute probs for all individuals
void UpdateLocusMeanProb(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, DoubleVec4d& LocusMeanProb, vector<vector<vector<int> > > & Genotype)
{ 
  for(int ind =0; ind< NIND; ind++){
    UpdateLocusMeanProb(ind, ExpTheta, SumExpTheta, LocusMeanProb, Genotype);
  }
}



void NormaliseMeanPi(vector<vector<double> > & MeanPi)
{

  for(int r=0; r<NREGION; r++){
    double sum = 0;
    for(int k=0; k<NSPECIES; k++){
      sum+= MeanPi[r][k];
    }
    for(int k=0; k<NSPECIES; k++){
      MeanPi[r][k] /= sum;
    }
  }
}

void NormaliseMeanFreq(DoubleVec4d& MeanFreq,double totaliter){
  for(int l=0;l<NLOCI; l++){
    for(int r=0;r<NREGION;r++){     
      for(int allele=0; allele<Nallele[l]; allele++){
	MeanFreq[r][0][l][allele]/=totaliter;
      }
    }
  } 
}


// modified to avoid underflow; it's slower but this runs only once
void ComputeLogMeanProb(DoubleVec3d& MeanProb, const DoubleVec4d& LocusMeanProb, double totaliter){
  double ltiter = log(totaliter);
  for(int ind = 0; ind < NIND; ind++){
    for(int r=0;r<NREGION;r++){
      for(int k=0; k<NSPECIES; k++){
        double logmp = 0;
	for(int l =0; l<NLOCI; l++){
          logmp += log(LocusMeanProb[ind][r][k][l]) - ltiter;
	}
        MeanProb[ind][r][k] = logmp;
      }
    }
  }
}

void NormaliseSubRegionProb(DoubleVec3d& SubRegionProb){
  for(int ind = 0; ind < NIND; ind++){
    double sum = 0;
    for(int r=0;r<NSUBREGION;r++){
      for(int k=0; k<NSPECIES; k++){
	sum += SubRegionProb[ind][r][k];
      }
    }
    for(int r=0;r<NSUBREGION;r++){
      for(int k=0; k<NSPECIES; k++){
	SubRegionProb[ind][r][k] /= sum; 
      }
    }
  }
}

// some output routines
void OutputPi(vector<vector<double> > & Pi,vector<string> & RegionName,ostream & ostr, vector<int> & Perm){
  for(int k = 0; k< NSPECIES; k++){
    for(int r = 0; r<Pi.size(); r++){
      //ostr.setf(ios::fixed);
      ostr.setf(ios::showpoint);
      ostr.precision(2); 
      ostr << Pi[Perm[r]][k] << "  ";
    }
  }
  ostr << endl;
}
  
// WARNING NEVER CALLED
void OutputEstimatedFreqs(const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, ostream & output, const IntVec2d& Coding, vector<int> & Perm)
{
  for(int l=0;l<NLOCI; l++){
    for(int allele=0; allele<Nallele[l]; allele++){
      output << setw(5) << Coding[l][allele] << " : ";
      for(int r=0;r<NREGION;r++){     
	output << std::fixed << setprecision(3) << setw(5) << ExpTheta[Perm[r]][0][l][allele]/SumExpTheta[Perm[r]][0][l] << " ";
      }
      output << endl;
    }
  }
}


void OutputLogMeanProb(const DoubleVec3d& MeanProb, ofstream & output, vector<int> & NMissing, vector<int> & Region, vector<string> & Id, vector<string> & RegionName, vector<int> & Perm, int first, int last, vector<int> & RegionsPresent)
{
  if(!LOCATE)
    output << "Warning: these results not produced via the -A option, so not cross-validated" << endl;

  output << "Id " << " " << "#Genotypes" << " Location" << " " << "Assigned" << " ";
  if(NSPECIES>1)
    output << "AssignedK ";

  OutputRegionNames(output,RegionName,Perm);
  
  for(int ind=first; ind<=last; ind++){
    int r;
    if(Region[ind] == -1)
      r = -1;
    else {
      for(r=0; r< NREGION; r++){
	if(Perm[r] == Region[ind])
	  break;
      }
    }
    
    output << Id[ind] << " " << (NLOCI - NMissing[ind]) << " " << ((r>=0) ? RegionsPresent[r]:r) << " ";

    int assignedr = 0;
    int assignedk = 0;
    double maxlogprob = MeanProb[ind][Perm[0]][0];
    
    for(int k=0; k<NSPECIES; k++){
      for(int r=0;r<(NREGION-LOCATE );r++){    
	if(MeanProb[ind][Perm[r]][k] > maxlogprob){
	  assignedr = r; assignedk = k; maxlogprob = MeanProb[ind][Perm[r]][k];
	}
      }
    }
			    
    output << RegionsPresent[assignedr] << " ";
    if(NSPECIES>1)
      output << assignedk << " ";

    for(int k=0; k<NSPECIES; k++){
      for(int r=0;r<(NREGION-LOCATE);r++){     
	output << std::scientific <<  MeanProb[ind][Perm[r]][k] << " ";
      }
    }
    output << endl;
  }
}


void OutputLogMeanProb2(const DoubleVec4d& MeanFreq, vector<vector<vector<int> > > & Genotype, ostream & output, vector<int> & NMissing, vector<int> & Region, vector<int> & Perm)
{
  for(int ind=0; ind<NIND; ind++){
    output << NMissing[ind] << " " << Region[ind] << " ";
    for(int r=0;r<NREGION;r++){     
      for(int k=0; k<NSPECIES; k++){
	output << std::scientific << log_hybrid_Prob(MeanFreq,Genotype, ind, Perm[r], Perm[r]) << " ";
      }
    }
    output << endl;
  }
}


void OutputParameters(ofstream & paramfile, vector<double> & Alpha, double & Beta, vector<double> & Gamma, vector<double> & Delta, double & Eta, const DoubleVec1d& Lambda, const DoubleVec3d& LogLik)
{
  for(int a = 0; a < Alpha.size(); a++)
    paramfile << Alpha[a] << " ";
  paramfile <<  Beta <<  " ";
  if(NSPECIES > 1){
    for(int k = 0; k< NSPECIES; k++)
      paramfile << Gamma[k] << " ";
    for(int a = 0; a < Delta.size(); a++)
      paramfile << Delta[a] << " ";
    paramfile <<  Eta <<  " ";
	for(int k = 0; k< NSPECIES; k++)
	  paramfile << Lambda[k] << " ";
  }
  paramfile <<  SumLogLik(LogLik);
  paramfile << endl;     
           
}     

void OutputTheta(double **** Theta, ostream & output, int ** Coding,vector<int> & Perm)
{
  for(int l=0;l<NLOCI; l++){
    for(int allele=0; allele<Nallele[l]; allele++){
      output << setw(5) << Coding[l][allele] << " : ";
      for(int r=0;r<NREGION;r++){     
	output << std::fixed << setprecision(3) << setw(5) << Theta[Perm[r]][0][l][allele] << " ";
      }
      output << endl;
    }
  }

}


// ------------------------- The MCMC update routines ---------------------

// The next 3 functions I've just added in response to the interest D+J
// expressed; so pretty much untested (previously Species set to 
// 0 for all individuals)
// 
// update vector of allocation variables (Species) which holds
// which species (= subpopulation) each individual is currently
// assigned to
void update_Species(vector<int> & Species, vector< vector<double> > & Pi, vector<int> & Region, vector<vector<vector<int> > > & Genotype, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta)
{
  vector<double> Prob(NSPECIES,0);
  for(int ind=0; ind<NIND; ind++){
    int r = Region[ind]; 
    if(r>=0){
      double sumprob = 0;
      for(int k=0; k<NSPECIES; k++){
	//compute the vector of  probs of individual is genotypes, 
	//for each species possible    
	Prob[k] = Pi[r][k];
	for(int l=0; l<NLOCI; l++){
	  for(int chr = 0; chr < 1; chr ++){
	    int allele = Genotype[ind][chr][l];
	    if(allele>=0) // if not missing data
	      Prob[k] *= ExpTheta[r][k][l][allele]/SumExpTheta[r][k][l];
	  }
	}
	sumprob += Prob[k];
      }
      Species[ind] = rint2(Prob,sumprob);
    }
  }
}


void compute_Pi(DoubleVec2d& Pi, const DoubleVec2d& ExpPsi, const DoubleVec1d& SumExpPsi )
{
  for(int r =0; r<NREGION; r++){
    for(int k = 0; k<NSPECIES; k++){
      Pi[r][k] = ExpPsi[r][k]/SumExpPsi[r];
    }
  }
}

void update_Pi(vector< vector<double> > & Pi, vector<int> & Species, vector<int> & Region)
{
  for(int r =0; r<NREGION; r++){
    vector<double> DirichletParam(NSPECIES,1); // set prior on Pi to be Di(0.1,0.1) [temporary to get something going]
    for(int ind =0; ind < NIND; ind++){
      if(Region[ind] == r)
         DirichletParam[Species[ind]] ++;
    }
    rdirichlet(DirichletParam,Pi[r].size(),Pi[r] );
  }
}



// update Nu (the species-specific adjustment to the background "ancestral" allele freqs)
void update_Nu(DoubleVec3d& Nu, vector<double> & Gamma, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount)
{

  double NewNu;
  // do not need initialization as initialized inline
  static DoubleVec2d NewTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewExpTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewSumExpTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewLogLik(NREGION,DoubleVec1d(NSPECIES,0.0));

  for(int k=0; k<NSPECIES; k++){
    for(int l=0; l<NLOCI; l++){
      for(int j=0;j<Nallele[l];j++){
	double LogLikRatio = 0;
	NewNu = rnorm(Nu[k][l][j],1);
	LogLikRatio = Gamma[k] * 0.5 * (Nu[k][l][j] * Nu[k][l][j] - NewNu * NewNu); //prior on Nu[k][l][j] is N(0,1/Gamma[k]); 
	for(int r=0; r<NREGION; r++){
	  NewTheta[r][k] = Theta[r][k][l][j] - Nu[k][l][j] + NewNu;
	  NewExpTheta[r][k] = exp(NewTheta[r][k]);
	  NewSumExpTheta[r][k] = SumExpTheta[r][k][l] - ExpTheta[r][k][l][j] + NewExpTheta[r][k];
	  NewLogLik[r][k] = LogLik[r][k][l] 
	    + Count[r][k][l][j] * (NewTheta[r][k] - Theta[r][k][l][j]) 
	    - SumCount[r][k][l] * (log(NewSumExpTheta[r][k]) - log(SumExpTheta[r][k][l]));
	  LogLikRatio += NewLogLik[r][k] - LogLik[r][k][l];
	}
	
	double A = exp(LogLikRatio); // acceptance prob
	NUATTEMPT +=1;
	if( ranf()<A ){ //accept move 
	  NUACCEPT +=1;
	  for(int r=0; r<NREGION; r++){
	    Theta[r][k][l][j] = NewTheta[r][k];
	    ExpTheta[r][k][l][j] = NewExpTheta[r][k];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r][k];
	    LogLik[r][k][l] = NewLogLik[r][k];
	  }
	  Nu[k][l][j] = NewNu;
	}
      }
    }
  }
}


// update Mu (the background "ancestral" allele freqs)
void update_Mu(DoubleVec2d& Mu, double Beta, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount)
{
  static DoubleVec2d NewTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewExpTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewSumExpTheta(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewLogLik(NREGION,DoubleVec1d(NSPECIES,0.0));
  double NewMu;

  for(int l=0; l<NLOCI; l++){
    for(int j=0;j<Nallele[l];j++){
      double LogLikRatio = 0;
      NewMu = rnorm(Mu[l][j],sqrt(1.0/Beta));
      LogLikRatio = Beta * 0.5 * (Mu[l][j] * Mu[l][j] - NewMu * NewMu); //prior on Mu is N(0,1/Beta) 
	
      for(int r=0; r<NREGION; r++){
	for(int k=0; k<NSPECIES; k++){
	  NewTheta[r][k] = Theta[r][k][l][j] - Mu[l][j] + NewMu;
	  NewExpTheta[r][k] = exp(NewTheta[r][k]);
	  NewSumExpTheta[r][k] = SumExpTheta[r][k][l] - ExpTheta[r][k][l][j] + NewExpTheta[r][k];
	  NewLogLik[r][k] = LogLik[r][k][l] 
	    + Count[r][k][l][j] * (NewTheta[r][k] - Theta[r][k][l][j]) 
	    - SumCount[r][k][l] * (log(NewSumExpTheta[r][k]) - log(SumExpTheta[r][k][l]));
	  LogLikRatio += NewLogLik[r][k] - LogLik[r][k][l];
	}
      }
	
      double A = exp(LogLikRatio); // acceptance prob

      MUATTEMPT +=1;
      if( ranf()<A ){ //accept move 
	MUACCEPT +=1;
	for(int r=0; r<NREGION; r++){
	  for(int k=0; k<NSPECIES; k++){
	    Theta[r][k][l][j] = NewTheta[r][k];
	    ExpTheta[r][k][l][j] = NewExpTheta[r][k];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r][k];
	    LogLik[r][k][l] = NewLogLik[r][k];
	  }
	}
	Mu[l][j] = NewMu;
      }
    }
  }
}

// update Lambda (the background species abundance)
void update_Lambda(DoubleVec1d& Lambda, double Eta, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Region, vector<int> & Species)
{
  static DoubleVec1d NewPsi(NREGION,0.0);
  static DoubleVec1d NewExpPsi(NREGION,0.0);
  static DoubleVec1d NewSumExpPsi(NREGION,0.0);

  for(int k=0;k<NSPECIES;k++){
    double NewLambda = rnorm(Lambda[k],1);
    double LogLikRatio = Eta * 0.5 * (Lambda[k] * Lambda[k] - NewLambda * NewLambda); //prior on Lambda[k] is N(0,1/eta) 
    
    for(int r=0; r<NREGION; r++){
      NewPsi[r] = Psi[r][k] - Lambda[k] + NewLambda;
      NewExpPsi[r] = exp(NewPsi[r]);
      NewSumExpPsi[r] = SumExpPsi[r] - ExpPsi[r][k] + NewExpPsi[r];
    }
    
    // maybe worth checking this?
    for(int ind=0; ind< NIND; ind++){
      int r0 = Region[ind];
      if(r0>=0){
	if(Species[ind] == k)
	  LogLikRatio += NewPsi[r0] - Psi[r0][k];
	LogLikRatio += log(SumExpPsi[r0]) - log(NewSumExpPsi[r0]);
      }
    }
      
    double A = exp(LogLikRatio); // acceptance prob
        
    LAMBDAATTEMPT +=1;
    if( ranf()<A ){ //accept move 
      LAMBDAACCEPT +=1;
      for(int r=0; r<NREGION; r++){
	Psi[r][k] = NewPsi[r];
	ExpPsi[r][k] = NewExpPsi[r];
	SumExpPsi[r] = NewSumExpPsi[r];
      }
      Lambda[k] = NewLambda;
    }
  }
}

void update_Alpha0(vector<double> & Alpha, const DoubleVec4d& X)
{
  double sum = 0;
  double sumsq = 0;
  int total =0;

  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0;l<NLOCI; l++){
	for(int allele=0; allele<Nallele[l]; allele++){
	  sumsq += X[r][k][l][allele] * X[r][k][l][allele];
          total += 1;
	}
      }
    }
  }
  Alpha[0] = rgamma(0.5*total+NA[0],0.5*sumsq+LA[0]);
}


void update_Delta0(vector<double> & Delta, const DoubleVec2d& Y)
{
  double sum = 0;
  int total =0;
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      sum += Y[r][k] * Y[r][k];
      total += 1;
    }
  }
  Delta[0] = rgamma(0.5*total+ND[0],0.5*sum+LD[0]);
}


// beta is the prior precision for Mu (ie Mu is N(0,1/Beta)
// prior on beta is Gamma(NBETA,LBETA)
void update_Beta(double & Beta, const DoubleVec2d& Mu)
{
  double sumsq = 0;
   double sum = 0;
  int total =0;
  for(int l=0;l<NLOCI; l++){
    for(int allele=0; allele<Nallele[l]; allele++){
      sumsq += Mu[l][allele] * Mu[l][allele];
      total += 1;
    }
  }
  Beta = rgamma(0.5*total+NBETA,0.5*sumsq+LBETA);
}



// Eta is the prior precision for Lambda (ie Lambda is N(0,1/Eta)
// prior on Eta is Gamma(NETA,LETA)
void update_Eta(double & Eta, const DoubleVec1d& Lambda)
{
  double sum = 0;
  int total =0;
  for(int k=0;k<NSPECIES; k++){
    sum += Lambda[k] * Lambda[k];
    total += 1;
  }
  Eta = rgamma(0.5*total+NETA,0.5*sum+LETA);
}


bool InElephantRange(double x, double y){
  double x0,y0,x1,y1,a,b;

// check (crudely) for falling in sea or sahara
  if(y>0.29)
    return false; // Sahara region (N of Timbuktoo)
  
  x0 = 0.746; // exclude NE ethiopia
  y0 = -0.007;
  x1 = 0.648;
  y1 = 0.140;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  x0 = 0.572; 
  y0 = 0.182;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  x1 = 0.445;
  y1 = 0.135;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b)){
    x0 = 0.391;
    y0 = 0.253;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y > (a*x + b))
      return false;
  }
  x0 = 0.391;
  y0 = 0.253;
  x1 = -0.022;
  y1 = 0.293;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if(y > (a*x + b))
    return false;
  
  // exclude south of africa
  
  if(y < -0.368){
    x0 = 0.559;
    y0 = -0.501;
    x1 = 0.449;
    y1 = -0.440;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y < (a*x + b))
      return false;
    x0 = 0.465;
    y0 = -0.368;
    a = (y0-y1)/(x0-x1);
    b = y0 - a* x0;
    if(y > (a*x + b))
      return false;
  }

  if(  (x < 0.157 ) 
       && ( y < (0.25*x+ PI*5.59/180) ) 
       && ( y < (-0.56*x+ PI*8.38/180) ) )
    return false; // bump in Lagos/Accra region
  
  if((x < 0.157) && (y<0.073) ) // SW coast
    return false;
  if((x < -0.2)) //W coast
    return false;
  if(y < -0.6) // S coast
    return false;
  if((x > 0.70) && (y<(-0.05)) )
    return false; // SE coast
  if((y < (-1.38 + 1.59 * x)))
    return false; // SEcoast, Durbin/Cidade de Nacala
  
  
  if((x > 0.70) && (y< -0.66 + (13.0/15)* x))
    return false; // Ecoast
 
  return true;
}

bool InForest(double x,double y){
  double x0,y0,x1,y1,a,b,a0,b0,a1,b1,a2,b2;
  bool forest = false;
  x0 = -0.003; // the part to the west of Mole
  y0 = 0.349;
  x1 = 0.038;
  y1 = 0.143;
  a = (y0-y1)/(x0-x1);
  b = y0 - a* x0;
  if( y < a*x + b){
    forest = true;
    //cout << "Is to the west of Mole" << endl;
  }
  
  // now check the rest of the forest areas
  x0 = 0.038;
  y0 = 0.143;
  x1 = 0.524;
  y1 = 0.087;
  a0 = (y0-y1)/(x0-x1);
  b0 = y0 - a0* x0;
  
  x0 = 0.489;
  y0 = -0.165;
  a1 = (y0-y1)/(x0-x1);
  b1 = y0 - a1* x0;
  
  x1 = 0.241;
  y1 = -0.157;
  a2 = (y0-y1)/(x0-x1);
  b2 = y0 - a2* x0;
  
  if( ( y < a0*x + b0) &&
      ( y > a1*x + b1) &&
      ( y > a2*x + b2))
    forest = true;


  // rule out anything north of 10deg

  if(y> (PI*10/180))
	  forest = false;
	  
  return forest;  
}


// This routine used to take LogLik by reference.  However, its internal
// calculations are a different formula than the LogLik of the rest of
// the program.  This required an immediate corrective call to calc_Loglik
// after each call to this function, which was an accident waiting to
// happen.  It now uses internal storage only.  -- Mary 12/30/2020

void update_Location(vector<double> & Alpha, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& SumExpTheta, DoubleVec3d& LogLik, double * L, const IntVec4d& Count, const IntVec3d& SumCount, vector<double> & Xcoord, vector<double> & Ycoord, vector<int> & Species, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<double> & BoundaryX, vector<double> & BoundaryY, const Mapgrid& mymapgrid)
{
  // do not need initialization, first use is assignment
  static DoubleVec2d NewTheta(NLOCI,DoubleVec1d(MAXNALLELE,0.0));
  static DoubleVec2d NewExpTheta(NLOCI,DoubleVec1d(MAXNALLELE,1.0));
  static DoubleVec1d NewSumExpTheta(NLOCI,0.0);
  static DoubleVec1d NewLogLik(NLOCI,0.0);
  static DoubleVec1d OldLogLik(NLOCI,0.0);

  // this one has to be an array for use in calc_L, alas.
  // It does not need initialization
  double * NewL = new double [NREGION * NREGION];

  for(int rep = 0; rep < 10; rep++){
    vector<double> NewXcoord(Xcoord);
    vector<double> NewYcoord(Ycoord);
    double h; // proposal sd
    double LogLikRatio = 0;
    
    if(ranf()<0.9){ // propose small move from current position
      h=0.02; //0.07;
      NewXcoord[NREGION-1] += rnorm(0,h);
      NewYcoord[NREGION-1] += rnorm(0,h);
      // do not bother with further details if proposed location is out of range
      if (!InRange(NewXcoord[NREGION-1],NewYcoord[NREGION-1],BoundaryX,BoundaryY,mymapgrid)) {
        // DEBUG
        TRACEFILE << to_degrees(NewXcoord[NREGION-1]) << " " << to_degrees(NewYcoord[NREGION-1]) << " Reject out of range" << endl;
        continue;
      }
      LogLikRatio = 0; // Hastings ratio
    } else { // propose jump to new randomly-chosen location
      h=0.04; //0.02; //sd of proposal
      int newregion = (int) (ranf() * (NREGION - 1));
      NewXcoord[NREGION-1] = Xcoord[newregion] + rnorm(0,h);
      NewYcoord[NREGION-1] = Ycoord[newregion] + rnorm(0,h);
      // do not bother with further details if proposed location is out of range
      if (!InRange(NewXcoord[NREGION-1],NewYcoord[NREGION-1],BoundaryX,BoundaryY,mymapgrid)) {
        // DEBUG
        TRACEFILE << to_degrees(NewXcoord[NREGION-1]) << " " << to_degrees(NewYcoord[NREGION-1]) << " Reject out of range" << endl;
        continue;
      }
      double forwardsprob = 0; // prob density of proposed move
      double backwardsprob = 0; // prob density of backwards move
      for(int r = 0; r < (NREGION-1); r++){
	forwardsprob += dnorm((NewXcoord[NREGION-1]-Xcoord[r])/h) * dnorm((NewYcoord[NREGION-1]-Ycoord[r])/h);
	backwardsprob += dnorm((Xcoord[NREGION-1]-Xcoord[r])/h) * dnorm((Ycoord[NREGION-1]-Ycoord[r])/h);
      }
      LogLikRatio = log(backwardsprob)-log(forwardsprob);
    }


    if(NONUNIFORMPRIOR){ // put prior that favours sample coming from close to 
      double h=0.000004; // one of the sampling locations
      double newprior = 0;
      double oldprior = 0;
      for(int r = 0; r < (NREGION-1); r++){
	 newprior += dnorm((NewXcoord[NREGION-1]-Xcoord[r])/sqrt(h)) * dnorm((NewYcoord[NREGION-1]-Ycoord[r])/sqrt(h));
	 oldprior += dnorm((Xcoord[NREGION-1]-Xcoord[r])/sqrt(h)) * dnorm((Ycoord[NREGION-1]-Ycoord[r])/sqrt(h));
      }
      LogLikRatio += log(newprior)-log(oldprior);
    }
      
    calc_L(NewL,Alpha,NewXcoord,NewYcoord);
    
    int r = NREGION - 1;
    int k = Species[SAMPLETOLOCATE];
    
    for(int l=0; l<NLOCI; l++){
      NewSumExpTheta[l] = 0; 
      for(int j=0;j<Nallele[l];j++){  
	NewTheta[l][j] = Mu[l][j]+Nu[k][l][j];
	for(int s = 0; s<=r; s++){
	  NewTheta[l][j] += NewL[r+s*NREGION] * X[s][k][l][j];
	}
	NewExpTheta[l][j] = exp(NewTheta[l][j]);
	NewSumExpTheta[l] += NewExpTheta[l][j];
      }
    }

    double tempprob = 1;
    double newtempprob = 1;

    
    for(int l=0; l<NLOCI; l++){
      NewLogLik[l] = 0; 
      OldLogLik[l] = 0;
      for(int ind =0; ind < NIND; ind++){
	if(Region[ind] == r){
	  int allele1 = Genotype[ind][0][l];
	  int allele2 = Genotype[ind][1][l];
	  double Newp1=1.0;
	  double Newp2 = 1.0;
	  double p1 = 1.0;
	  double p2 = 1.0;
	  if(allele1>=0){
	    Newp1 = (1-DELTA) * NewExpTheta[l][allele1]/NewSumExpTheta[l] + DELTA/Nallele[l];
	    p1 = (1-DELTA) * ExpTheta[r][k][l][allele1]/SumExpTheta[r][k][l] + DELTA/Nallele[l];
	  }
	  if(allele2>=0){
	    Newp2 = (1-DELTA) * NewExpTheta[l][allele2]/NewSumExpTheta[l] + DELTA/Nallele[l];
	    p2 = (1-DELTA) * ExpTheta[r][k][l][allele2]/SumExpTheta[r][k][l] + DELTA/Nallele[l];
	  }
	  
	  if(allele1 == allele2){
	    NewLogLik[l] += log(NULLPROB*Newp1 + (1-NULLPROB) * Newp1 * Newp1);
	    OldLogLik[l] += log( NULLPROB*p1 + (1-NULLPROB) * p1 * p1);
	  }
	  else{
	    NewLogLik[l] += log((1-NULLPROB) *  Newp1 * Newp2);
	    OldLogLik[l] += log((1-NULLPROB) *  p1 * p2 );
	  }
	  
	}
      }
      LogLikRatio += NewLogLik[l] - OldLogLik[l];
      tempprob *= exp(OldLogLik[l]); 

      newtempprob *= exp(NewLogLik[l]);
    }

    LOCATTEMPT +=1;
    double A = exp(LogLikRatio); // acceptance prob
    
    // test
    if( ranf()<A ){ //accept move
      // DEBUG
      TRACEFILE << to_degrees(NewXcoord[NREGION-1]) << " " << to_degrees(NewYcoord[NREGION-1]) << " Accept" << endl;
      LOCACCEPT +=1;
       Xcoord = NewXcoord;
       Ycoord = NewYcoord;
       for(int r0=0; r0<NREGION; r0++){
	 for(int s=0; s <NREGION; s++){
	   L[r0+ s*NREGION] = NewL[r0 + s*NREGION];
	 }
       }
       for(int l=0; l<NLOCI; l++){	
         // following line no longer needed
	 // OldLogLik[l] = NewLogLik[l];
	 SumExpTheta[r][k][l] = NewSumExpTheta[l];
	 for(int j=0;j<Nallele[l];j++){  
	   Theta[r][k][l][j] = NewTheta[l][j];
	   ExpTheta[r][k][l][j] = exp(NewTheta[l][j]);
	 }
       }
    } else {
      TRACEFILE << to_degrees(NewXcoord[NREGION-1]) << " " << to_degrees(NewYcoord[NREGION-1]) << " Reject" << endl;
    }
  }

  delete [] NewL;
  // we have modified things on which the log likelihood depends, so update it here
  calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
}


// update the parameters in the covariance matrix 
void update_Alpha(vector<double> & Alpha, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L, vector<double> & Xcoord, vector<double> & Ycoord)
{

static DoubleVec4d NewTheta(NREGION,DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,
  DoubleVec1d(MAXNALLELE,0.0))));
static DoubleVec4d NewExpTheta(NREGION,DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,
  DoubleVec1d(MAXNALLELE,0.0))));
static DoubleVec3d NewSumExpTheta(NREGION,DoubleVec2d(NSPECIES,
  DoubleVec1d(NLOCI,0.0)));
static DoubleVec3d NewLogLik(NREGION,DoubleVec2d(NSPECIES,
  DoubleVec1d(NLOCI,0.0)));

// This one has to be an array
static double* NewL;

static int FirstTime = 1;
static int FirstNREGION;
static int RecursionCheck = 0;

	if( FirstTime == 1 ) {
 		NewL = new double[NREGION*NREGION];
 		FirstNREGION = NREGION;
 		FirstTime = 0;
 	}
 	if( ValidateAssumptions ) {
 		if( RecursionCheck++ > 1 || NREGION !=FirstNREGION ) {
 			cerr << "Assumption violated" << endl;
 			exit(1);
 		}
 	}

  update_Alpha0(Alpha, X);

  for(int alphaparam = 1; alphaparam < ALPHALENGTH; alphaparam++){
    ALPHAATTEMPT[alphaparam] +=1;
    vector<double> NewAlpha(Alpha);
  
    double LogLikRatio = 0;
    
    // random walk on log alphas
    
    vector<double> AlphaUpdateSD = vector<double>(4,ALPHAUPDATESD);
    double h=rnorm(0,AlphaUpdateSD[alphaparam]);
    
    if(alphaparam != 2){
      NewAlpha[alphaparam] *= exp(h); // add h to log alpha
    }
    else{
      NewAlpha[alphaparam] += h;
    }

    
// reject new choice if proposal falls outside valid range    
    if(NewAlpha[alphaparam]<ALPHAMIN[alphaparam])
      continue;
   if(NewAlpha[alphaparam]>ALPHAMAX[alphaparam])
      continue;
    
    calc_L(NewL,NewAlpha,Xcoord,Ycoord);
    
    for(int r=0; r<NREGION; r++){
      for(int k=0; k<NSPECIES; k++){
	for(int l=0; l<NLOCI; l++){
	  NewSumExpTheta[r][k][l] = 0;
	  for(int j=0;j<Nallele[l];j++){  
	    NewTheta[r][k][l][j] = Mu[l][j]+Nu[k][l][j];
	    for(int s = 0; s<=r; s++){
	      NewTheta[r][k][l][j] += NewL[r+s*NREGION] * X[s][k][l][j];
	    }
            NewExpTheta[r][k][l][j] = exp(NewTheta[r][k][l][j]);
            NewSumExpTheta[r][k][l] += NewExpTheta[r][k][l][j];
	  }
        }
      }
    }
     
    // This one is definitely needed
    calc_LogLik(NewLogLik,NewTheta,NewSumExpTheta,Count,SumCount);

    for(int r=0; r<NREGION; r++){ 
      for(int k=0; k<NSPECIES; k++){
	for(int l=0; l<NLOCI; l++){
	  LogLikRatio += NewLogLik[r][k][l] - LogLik[r][k][l];
	}
      }
    }
    
    double A = exp(LogLikRatio); // acceptance prob
    
    if( ranf()<A ){ //accept move 
      ALPHAACCEPT[alphaparam] +=1;
      Alpha[alphaparam] = NewAlpha[alphaparam];
      for(int r=0; r<NREGION; r++){
	for(int s=0; s <NREGION; s++){
	  L[r+ s*NREGION] = NewL[r + s*NREGION];
	}
      }
      for(int r=0; r<NREGION; r++){ 
	for(int k=0; k<NSPECIES; k++){
	  for(int l=0; l<NLOCI; l++){
	    for(int j=0;j<Nallele[l];j++){  
	      Theta[r][k][l][j] = NewTheta[r][k][l][j];
	      ExpTheta[r][k][l][j] = NewExpTheta[r][k][l][j];
	    }
	    SumExpTheta[r][k][l] = NewSumExpTheta[r][k][l];
	    LogLik[r][k][l] = NewLogLik[r][k][l];
	  }
	}
      }
    } 
  }

  if( ValidateAssumptions ) RecursionCheck--;

}
// update the parameters in the covariance matrix M 
void update_Delta(vector<double> & Delta, const DoubleVec2d& Y, const DoubleVec1d& Lambda, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Species, vector<int> & Region, double * M, vector<double> & Xcoord, vector<double> & Ycoord)
{
  // These don't need initialization as it's handled inline
  static DoubleVec2d NewPsi(NREGION,DoubleVec1d(NSPECIES,0.0));
  static DoubleVec2d NewExpPsi(NREGION,DoubleVec1d(NSPECIES,1.0));
  static DoubleVec1d NewSumExpPsi(NREGION,0.0);

  // This one has to be an array for calc_L use
  double * NewM = new double [NREGION * NREGION];

  
  update_Delta0(Delta, Y);

  for(int deltaparam = 1; deltaparam < DELTALENGTH; deltaparam++){
    DELTAATTEMPT[deltaparam] +=1;
    vector<double> NewDelta(Delta);
  
    double LogLikRatio = 0;
    
    // random walk on log deltas
    
    vector<double> DeltaUpdateSD = vector<double>(4,0.4);
    double h=rnorm(0,DeltaUpdateSD[deltaparam]);
    
    if(deltaparam != 2){
      NewDelta[deltaparam] *= exp(h); // add h to log alpha
    }
    else{
      NewDelta[deltaparam] += h;
    }

    
// reject new choice if proposal falls outside valid range    
    if(NewDelta[deltaparam]<DELTAMIN[deltaparam])
      continue;
   if(NewDelta[deltaparam]>DELTAMAX[deltaparam])
      continue;
    
    calc_L(NewM,NewDelta,Xcoord,Ycoord);
    
    for(int r=0; r<NREGION; r++){
      NewSumExpPsi[r] = 0;
      for(int k=0; k<NSPECIES; k++){
	NewPsi[r][k] = Lambda[k];
	for(int s = 0; s<=r; s++){
	  NewPsi[r][k] += NewM[r+s*NREGION] * Y[s][k];
	}
	NewExpPsi[r][k] = exp(NewPsi[r][k]);
	NewSumExpPsi[r] += NewExpPsi[r][k];
      }
    }
    
    // maybe worth improving efficiency here?
    for(int ind=0; ind< NIND; ind++){
      int r = Region[ind];
      int k = Species[ind];
      LogLikRatio += log( (NewExpPsi[r][k]/NewSumExpPsi[r]) / (ExpPsi[r][k]/SumExpPsi[r]) );
    }

    double A = exp(LogLikRatio); // acceptance prob
    
    if( ranf()<A ){ //accept move 
      DELTAACCEPT[deltaparam] +=1;
      Delta[deltaparam] = NewDelta[deltaparam];
      for(int r=0; r<NREGION; r++){
	for(int s=0; s <NREGION; s++){
	  M[r+ s*NREGION] = NewM[r + s*NREGION];
	}
      }
      for(int r=0; r<NREGION; r++){ 
	for(int k=0; k<NSPECIES; k++){
	  Psi[r][k] = NewPsi[r][k];
	  ExpPsi[r][k] = NewExpPsi[r][k];
	}
	SumExpPsi[r] = NewSumExpPsi[r];
      }
    }
    
  }

 delete [] NewM;

}

// derivative of LogLik with respect to x(r,k,l,j) 
double divLogLikValue(int r, int k, int l, int j, const DoubleVec4d& ExpTheta, const DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L)
{
  double sum = 0;
  for(int s=r; s<NREGION; s++){
    sum += (Count[s][k][l][j] - SumCount[s][k][l] * ExpTheta[s][k][l][j]/SumExpTheta[s][k][l]) * L[s + r*NREGION];
  }
  return sum;
}

double calcNewdivLogLik(int r, int k, int l, int j, const DoubleVec1d& ExpTheta, const DoubleVec1d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L)
{
  double sum = 0;
  for(int s=r; s<NREGION; s++){
    sum += (Count[s][k][l][j] - SumCount[s][k][l] * ExpTheta[s]/SumExpTheta[s]) * L[s + r*NREGION];
  }
  return sum;
}


// update the Xs for a particular allele and locus, in a particular species, across all regions at once
void update_XJoint(vector<double> & Alpha, DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L)
{

  // These do not need initialization as are initialized at use
  static DoubleVec1d NewTheta(NREGION,0.0);
  static DoubleVec1d NewExpTheta(NREGION,0.0);
  static DoubleVec1d NewSumExpTheta(NREGION,0.0);

  static DoubleVec1d NewLogLik(NREGION,0.0);
  static DoubleVec1d NewX(NREGION,0.0);
  double NewdivLogLik;
  double h = XPROPOSALFACTOR * sqrt(1.0/Alpha[0]);

  for(int k=0; k<NSPECIES; k++){
    for(int l=0; l<NLOCI; l++){
      for(int j=0;j<Nallele[l];j++){
	double LogLikRatio = 0;
	for(int r=0; r<NREGION; r++){
	  double propmean = X[r][k][l][j]; // proposal mean
	  if(USELANGEVIN == 1)
	    propmean += (h*h/2)*(divLogLikValue(r,k,l,j,ExpTheta,SumExpTheta,Count,SumCount,L)- Alpha[0]* X[r][k][l][j]);
	  NewX[r] = rnorm(propmean,h);
	  LogLikRatio += 0.5*Alpha[0]*((X[r][k][l][j]*X[r][k][l][j]) - (NewX[r]*NewX[r])); // this is the "prior" part of the acceptance ratio
	  if(USELANGEVIN == 1)
	    LogLikRatio += 0.5* (NewX[r] - propmean)*  (NewX[r] - propmean)/(h*h); // bottom part of Hastings ratio X -> NewX
	} 
	
	for(int r=0; r<NREGION; r++){
	  NewTheta[r] = Mu[l][j] + Nu[k][l][j];
	  NewSumExpTheta[r] = 0;
	  for(int s = 0; s<=r; s++){
	    NewTheta[r] += L[r+s*NREGION] * NewX[s];
	  }
	  NewExpTheta[r] = exp(NewTheta[r]);
	  NewSumExpTheta[r] = SumExpTheta[r][k][l] - ExpTheta[r][k][l][j] + NewExpTheta[r];
	  
	  NewLogLik[r] = LogLik[r][k][l] 
	    + Count[r][k][l][j] * (NewTheta[r] - Theta[r][k][l][j]) 
	    - SumCount[r][k][l] * 
	    (log(NewSumExpTheta[r]) - log(SumExpTheta[r][k][l]));
	  LogLikRatio += (NewLogLik[r] - LogLik[r][k][l])/TEMPERATURE; // Likelihood ratio part of acceptance ratio
	}
	
	if(USELANGEVIN == 1){ // compute top part of Hastings ratio (H ratio is 1 if Langevin updates not used)
	  for(int r=0; r<NREGION; r++){ 
	    // backpropmean is the proposal mean when going from NewX
	    NewdivLogLik = calcNewdivLogLik(r,k,l,j,NewExpTheta,NewSumExpTheta,Count,SumCount,L);
	    double backpropmean = NewX[r] + (h/2) * (NewdivLogLik - Alpha[0] * NewX[r]);
	    LogLikRatio -= 0.5 * (X[r][k][l][j] - backpropmean)*(X[r][k][l][j] - backpropmean)/(h*h); // top part of Hastings ratio
	  }
	}
	double A = exp(LogLikRatio); // acceptance prob
	XATTEMPT +=1;
	
	if( ranf()<A ){ //accept move 
	  XACCEPT +=1;
	  for(int r=0; r<NREGION; r++){
	    Theta[r][k][l][j] = NewTheta[r];
	    ExpTheta[r][k][l][j] = NewExpTheta[r];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r];
	    X[r][k][l][j] = NewX[r];
	    LogLik[r][k][l] = NewLogLik[r];    
	  }
	}
      }
    }
  }
}

// update each X individually (better acceptance rate/ larger proposal variance,// but more likelihood evaluations!)
void update_XSingle(vector<double> & Alpha, DoubleVec4d& X, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, const IntVec4d& Count, const IntVec3d& SumCount, double * L)
{
   static DoubleVec1d NewTheta(NREGION,0.0);
   static DoubleVec1d NewExpTheta(NREGION,0.0);
   static DoubleVec1d NewSumExpTheta(NREGION,0.0);
   static DoubleVec1d NewLogLik(NREGION,0.0);
 
   static int FirstTime = 1;
   static int FirstNREGION;
   static int RecursionCheck = 0;

   double NewX;
  double NewdivLogLik;
  double h = XPROPOSALFACTOR * sqrt(1.0/Alpha[0]);

     if( FirstTime == 1 ) {
 	  FirstNREGION = NREGION;
 	  FirstTime = 0;
   }
   if( ValidateAssumptions ) {
 	  if( RecursionCheck++ > 1 || NREGION != FirstNREGION ) {
 		  cerr << "Assumption violated" << endl;
 		  exit(1);
 	  }
   }

  for(int k=0; k<NSPECIES; k++){
    for(int l=0; l<NLOCI; l++){
      for(int j=0;j<Nallele[l];j++){	
	for(int r=0; r<NREGION; r++){
	  double LogLikRatio = 0;
	  double propmean = X[r][k][l][j]; // proposal mean
	  
	  if(USELANGEVIN == 1) 
	    propmean += (h*h/2)*(divLogLikValue(r,k,l,j,ExpTheta,SumExpTheta,Count,SumCount,L) - Alpha[0] * X[r][k][l][j]);
	  
	  NewX = rnorm(propmean,h);
	  LogLikRatio += 0.5* Alpha[0] * ( X[r][k][l][j] * X[r][k][l][j] - NewX * NewX); // this is the "prior" part of the acceptance ratio
	  if(USELANGEVIN == 1)
	    LogLikRatio += 0.5* (NewX - propmean)*  (NewX - propmean)/(h*h); // bottom part of Hastings ratio X -> NewX
	
	
	  for(int s=r; s<NREGION; s++){
	    NewTheta[s] = Theta[s][k][l][j] + L[s+r*NREGION] * (NewX - X[r][k][l][j]);
	    NewExpTheta[s] = exp(NewTheta[s]);
	    NewSumExpTheta[s] = SumExpTheta[s][k][l] - ExpTheta[s][k][l][j] + NewExpTheta[s];
	    NewLogLik[s] = LogLik[s][k][l] 
	      + Count[s][k][l][j] * (NewTheta[s] - Theta[s][k][l][j]) 
	      - SumCount[s][k][l] * 
	      (log(NewSumExpTheta[s]) - log(SumExpTheta[s][k][l]));
	    LogLikRatio += (NewLogLik[s] - LogLik[s][k][l])/TEMPERATURE; // Likelihood ratio part of acceptance ratio
	  }
	  
	  if(USELANGEVIN == 1){ // compute top part of Hastings ratio (H ratio is 1 if Langevin updates not used)
	      // backpropmean is the proposal mean when going from NewX
	    NewdivLogLik = calcNewdivLogLik(r,k,l,j,NewExpTheta,NewSumExpTheta,Count,SumCount,L);
	    double backpropmean = NewX + (h*h/2) * (NewdivLogLik - Alpha[0] * NewX);
	    LogLikRatio -= 0.5 * (X[r][k][l][j] - backpropmean)*(X[r][k][l][j] - backpropmean)/(h*h); // top part of Hastings ratio
	  }
	  
	  double A = exp(LogLikRatio); // acceptance prob
	  XATTEMPT +=1;
	  
	  if( ranf()<A ){ //accept move 
	    XACCEPT +=1;
	    X[r][k][l][j] = NewX;
	    for(int s=r; s<NREGION; s++){
	      Theta[s][k][l][j] = NewTheta[s];
	      ExpTheta[s][k][l][j] = NewExpTheta[s];
	      SumExpTheta[s][k][l] = NewSumExpTheta[s];  
	      LogLik[s][k][l] = NewLogLik[s];
	    }
	  }
	}
      }
    }
  }
  if( ValidateAssumptions ) RecursionCheck--;
}

// update each Y individually (better acceptance rate/ larger proposal variance,// but more likelihood evaluations!)
void update_YSingle(vector<double> & Delta, DoubleVec2d& Y, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, vector<int> & Region, vector<int> & Species, double * M)
{
  // These do not need initialization as are initialized inline
  static DoubleVec1d NewPsi(NREGION,0.0);
  static DoubleVec1d NewExpPsi(NREGION,0.0);
  static DoubleVec1d NewSumExpPsi(NREGION,0.0);
 
  double NewY;
  double h = YPROPOSALFACTOR * sqrt(1.0/Delta[0]);

  for(int k=0; k<NSPECIES; k++){   
    for(int r=0; r<NREGION; r++){
      double LogLikRatio = 0;
      double propmean = Y[r][k]; // proposal mean
      NewY = rnorm(propmean,h);
      
      LogLikRatio += 0.5* Delta[0] * ( Y[r][k] * Y[r][k] - NewY * NewY); // this is the "prior" part of the acceptance ratio
      	
      for(int s=r; s<NREGION; s++){
	NewPsi[s] = Psi[s][k] + M[s+r*NREGION] * (NewY - Y[r][k]);
	NewExpPsi[s] = exp(NewPsi[s]);
	NewSumExpPsi[s] = SumExpPsi[s] - ExpPsi[s][k] + NewExpPsi[s];	
      }
      
      // maybe worth checking this?
      for(int ind=0; ind< NIND; ind++){
	int r0 = Region[ind];
	if(r0 >= r){
	  if(Species[ind] == k)
	    LogLikRatio += NewPsi[r0] - Psi[r0][k];
	  LogLikRatio += log(SumExpPsi[r0]) - log(NewSumExpPsi[r0]);
	}
      }
      
      double A = exp(LogLikRatio); // acceptance prob
      
      YATTEMPT +=1;
      
      if( ranf()<A ){ //accept move 
	YACCEPT +=1;
	Y[r][k] = NewY;
	for(int s=r; s<NREGION; s++){
	  Psi[s][k] = NewPsi[s];
	  ExpPsi[s][k] = NewExpPsi[s];
	  SumExpPsi[s] = NewSumExpPsi[s];  
	}
      }
    }
  }
}


void DoAllUpdates(DoubleVec4d& X, double & Beta,  vector<double> & Gamma, vector<double> & Alpha, DoubleVec2d& Mu, DoubleVec3d& Nu, DoubleVec4d& Theta, DoubleVec4d& ExpTheta, DoubleVec3d& LogLik, DoubleVec3d& SumExpTheta, IntVec4d& Count, IntVec3d& SumCount, double * L, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<int> & Species, vector<vector<double> > & Pi, vector<double> & Xcoord, vector<double> & Ycoord, DoubleVec2d& Y, double & Eta, vector<double> & Delta, DoubleVec1d& Lambda, DoubleVec2d& Psi, DoubleVec2d& ExpPsi, DoubleVec1d& SumExpPsi, double * M, vector<double> & BoundaryX, vector<double> & BoundaryY, const Mapgrid& mymapgrid )
{
  if(UPDATEBETA ==1)
    update_Beta(Beta,Mu);
  if(UPDATEX){
    if(UPDATEJOINT ==1)
      update_XJoint(Alpha,X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L);
    else
      update_XSingle(Alpha, X,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L);
  }
  if(UPDATEALPHA ==1)
    update_Alpha(Alpha,X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Xcoord,Ycoord);

  if(UPDATEMU == 1)
    update_Mu(Mu,Beta,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount);

  if(UPDATENU == 1){
    update_Nu(Nu,Gamma,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount);
    for(int k=0; k<NSPECIES; k++)
      update_Beta(Gamma[k],Nu[k]);
  }

  if(LOCATE)
    update_Location(Alpha,X,Mu,Nu,Theta,ExpTheta,SumExpTheta,LogLik,L,Count,SumCount,Xcoord,Ycoord,Species,Genotype,Region,BoundaryX, BoundaryY, mymapgrid);

  if(NSPECIES>1){
    update_Species(Species, Pi, Region, Genotype, ExpTheta, SumExpTheta);
    count_up_alleles(Count,Region,Species,Genotype);
    calc_SumCount(Count, SumCount);
    update_Lambda(Lambda,  Eta, Psi, ExpPsi, SumExpPsi, Region, Species);
    update_Eta(Eta,Lambda);
    update_YSingle(Delta,Y,Psi,ExpPsi,SumExpPsi,Region,Species,M);
    update_Delta(Delta,Y,Lambda,Psi,ExpPsi,SumExpPsi,Species,Region,M,Xcoord,Ycoord);
    compute_Pi(Pi, ExpPsi, SumExpPsi);

    // This one is probably needed -- not checking at this time
    calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
  }
#ifndef NDEBUG
  if (!CheckConsistency(Theta,ExpTheta,SumExpTheta,Count,SumCount,Psi,
     ExpPsi,SumExpPsi)) exit(-1);
#endif
 
}

void InitialiseTheta(DoubleVec4d& Theta, const DoubleVec4d& X, const DoubleVec2d& Mu, const DoubleVec3d& Nu, double * L){
  
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	for(int j=0;j<MAXNALLELE;j++){	 
	  Theta[r][k][l][j] = Mu[l][j] + Nu[k][l][j];
	  for(int s = 0; s<=r; s++){
	    Theta[r][k][l][j] += L[r+s*NREGION] * X[s][k][l][j];
	  }
	}
      }
    } 
  }
}

void Initialise(DoubleVec4d& X, double & Beta,  vector<double> & Gamma, vector<double> & Alpha, DoubleVec2d& Mu, DoubleVec3d& Nu, DoubleVec4d& Theta, double * L, vector<double> & Xcoord, vector<double> & Ycoord){ 

  if(INCLUDENUGGET)
    Alpha[3] = 1.0;

  Alpha[0] = ALPHA0;
  Alpha[1] = ALPHA1;
  Alpha[2] = ALPHA2;
  Beta = BETA;

  calc_L(L,Alpha,Xcoord,Ycoord);
  
  Gamma = vector<double>(NSPECIES,0.3); // Gamma[k] is 1/variance of Nu[k]
  
  for(int l=0; l<NLOCI; l++){
    for(int j=0; j<MAXNALLELE; j++){
      if(UPDATEMU==1){
	Mu[l][j] = rnorm(0,sqrt(1.0/Beta));
      }
      else
	Mu[l][j] = 0;
    }
  }
  for(int k=0; k<NSPECIES; k++){
    for(int l=0; l<NLOCI; l++){
      for(int j=0; j<MAXNALLELE; j++){
	if(UPDATENU==1)
	  Nu[k][l][j] = rnorm(0,1.0/sqrt(Gamma[k]));
	else
	  Nu[k][l][j] = 0;
      }
    }
  }
  
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	for(int j=0;j<MAXNALLELE;j++){	 
	  if(UPDATEX==1)
	    X[r][k][l][j] = rnorm(0,sqrt(1.0/Alpha[0])); //rnorm(0,1);
	  else
	    X[r][k][l][j] = 0;
	}
      }
    }
  }
  InitialiseTheta(Theta,X,Mu,Nu,L);
}


void output_empirical_freqs(vector<string> & RegionName,vector<int> & Perm,const IntVec2d& Coding, const IntVec4d& Count, const IntVec3d& SumCount){
	
 	cout << "Empirical Frequencies:" << endl;
    for(int l=0;l<NLOCI; l++){
		cout << "Locus " << (l+1) << endl;
		cout << setw(5) << " " << " : ";
		OutputRegionNames(cout, RegionName, Perm);
		for(int allele=0; allele<Nallele[l]; allele++){
		    cout << setw(5) << Coding[l][allele] << " : ";
		    for(int r=0;r<NREGION;r++){
		        cout << std::fixed << setprecision(3) << setw(10) << (1.0*Count[Perm[r]][0][l][allele]) << "/" << SumCount[Perm[r]][0][l] << " ";
		    }
		    cout << endl;
		}
	}
}

void OutputMeanFreq(ofstream & freqfile, vector<string> & RegionName, vector<int> & Perm, const IntVec2d& Coding, const DoubleVec4d& MeanFreq )
{
   freqfile << "Posterior Mean Freqs:" << endl;
   for(int l=0;l<NLOCI; l++){
      freqfile << "Locus " << (l+1) << endl;
       freqfile << setw(5) << " " << " : ";
  	    OutputRegionNames(freqfile, RegionName, Perm);
  	    for(int allele=0; allele<Nallele[l]; allele++){
   	      freqfile << setw(5) << Coding[l][allele] << " : ";
          for(int r=0;r<NREGION-LOCATE;r++){
   		    double f=MeanFreq[Perm[r]][0][l][allele];
  		    if(f<1e-4) f=0;
   		    freqfile << std::fixed << setprecision(3) << setw(10) << f << " ";
  	      }
   	      freqfile << endl;
        }
   }
}

int main ( int argc, char** argv)
{
  // MENU PROCESSING PHASE
  cout << "SCAT version " << VERSION << endl;
  int SEED = 0;
  map<string, string> filenames; 
  while( ( argc > 1 ) && ( argv[1][0] == '-' ) ) {
    switch(argv[1][1]) {
 
    case 'a':
      ++argv; --argc; ALPHAUPDATESD = atof(&argv[1][0]);
      break;

    case 'A': // estimate the location of a sample
      LOCATE = 1;
      ++argv;
      --argc;
      FIRSTSAMPLETOLOCATE = atoi(&argv[1][0])-1;
      ++argv;
      --argc;
      LASTSAMPLETOLOCATE = atoi(&argv[1][0])-1;
      cout << "Locating Samples " << FIRSTSAMPLETOLOCATE+1 << " to " <<  LASTSAMPLETOLOCATE+1 << endl;
      break;
 
    case 'b': // don't update beta
      UPDATEBETA = 0;
      break;

    case 'B': // read in boundary file 
      filenames["boundaryfile"] = argv[1]+2;
      READBOUNDARY = 1; 
      break;

    case 'C':
      ++argv; --argc;  SKIPCOL =  atoi(&argv[1][0]);
      break;

    case 'D':
      FORESTONLY = 1;
      break;

    case 'd':
      SAVANNAHONLY = 1;
      break;

    case 'e':
      ++argv; --argc; DELTA = atof(&argv[1][0]);
      ++argv; --argc; NULLPROB = atof(&argv[1][0]);
      break;

    case 'E':
      ECHOINPUTS = 1;
      break;

    case 'f': // fix alpha and beta, to values given in subsequent arguments
      UPDATEALPHA = 0;
      UPDATEBETA = 0;
      ++argv; --argc; ALPHA0 = atof(&argv[1][0]);
      ++argv;--argc; ALPHA1 = atof(&argv[1][0]);    
      ++argv;--argc;  ALPHA2 = atof(&argv[1][0]);
      ++argv;--argc; BETA = atof(&argv[1][0]);
      break;

    case 'F':
      OUTPUTSAMPLEFREQ = 1;
      break;
    
    case 'g': // read in grid file, replacement for a boundary file
      ++argv; --argc;
      filenames["gridfile"] = argv[1];
      READGRID = 1;
      break;

    case 'h': // set proposal variance for X update
      ++argv; --argc; XPROPOSALFACTOR = atof(&argv[1][0]);
      cout << "X proposal factor = " << XPROPOSALFACTOR << endl;
      break;
    
    case 'H':
      HYBRIDCHECK = atoi(&argv[1][2]);
	  break;
     
    case 'i' :
      ++argv; --argc; ALPHA0 = atof(&argv[1][0]);
      ++argv;--argc; ALPHA1 = atof(&argv[1][0]);    
      ++argv;--argc;  ALPHA2 = atof(&argv[1][0]);
      ++argv;--argc; BETA = atof(&argv[1][0]);
      break;
      
    case 'I': // don't permute the regions - use input order
      PERMUTE = 0;
      break;
   
    case 'j': // update x jointly
      UPDATEJOINT = 1;
      break;

    case 'M': // filename of samples to be assigned
      ASSIGNFILE = 1;
      filenames["assignfile"] = argv[1]+2;
      break;
		
    case 'm': // fix mu to be 0 (mimic "indep frequencies" model)
      UPDATEMU = 0;
      UPDATEBETA = 0;
      break;

    case 'n':
      UPDATENU =0;
      break;

    case 'N':
      INCLUDENUGGET = 1;
      break;

    case 'p':
      PSEUDOCOUNT = atoi(&argv[1][2]);
      break;

    case 'r': // not random walk
      USELANGEVIN = 1;
      break;

    case 'R': // remove all samples from a region when doing location
      REMOVEREGION = 1;
      break;

    case 'S': // seed
      ++argv;
      --argc;
      SEED = atoi(&argv[1][0]);
      break;

    case 'T': // number regions from 0
       ++argv;
      --argc;
      NUMBERREGIONSFROM = atoi(&argv[1][0]);
      break;

    case 'v':
      VERBOSE = 1;
      OUTPUTX = 1;
      filenames["Xfile"] = argv[1]+2;
      break;

    case 'w': // don't use spatial smoothing
      USESPATIAL = 0;
      break;

    case 'W': // locate whole region
      LOCATEWHOLEREGION = 1;
      break;

    case 'x': // nonuniformprior
      NONUNIFORMPRIOR = 1;
      break;
    
    case 'X': // start location close to true location
      CHEAT = 1;
      break;
      
    case 'Z' : // include subregion info in location file
      USESUBREGION = 1;
      break;

    default: 
      cerr << "Error: option " << argv[1] << "unrecognized" << endl;
      return 1;
    
    }
    ++argv;
    --argc;
  }

  // error checking
  if (SAVANNAHONLY && FORESTONLY) {
    string msg = "Asked for savannah only AND forest only!?";
    error_and_exit(msg);
  }

  if (SKIPCOL < 0) {
    error_and_exit("Asked to skip negative columns?!");
  }
  
  
  if(ASSIGNFILE==1)
	  if(HYBRIDCHECK==0)
		  LOCATE = 1;
  
  if(LOCATE) {
    if(LASTSAMPLETOLOCATE < 0)
      LASTSAMPLETOLOCATE = FIRSTSAMPLETOLOCATE;
    if(LASTSAMPLETOLOCATE < FIRSTSAMPLETOLOCATE) {
      string msg = "First and last samples to locate are incoherent";
      error_and_exit(msg);
    }
  }

  if(argc<5){
    cerr << "Usage is ./SCAT genotypefile locationfile outputdir NLOCI Niter Nthin Nburn " << endl;
    cerr << "Niter, Nthin and Nburn are optional" << endl;
    cerr << "See instructions for all additional options" << endl;
    exit(1);
  }

 

  cout << "MAXNALLELE = " << MAXNALLELE <<endl;
  init_genrand(SEED);

  filenames["input"]  = argv[1];
  filenames["regions"] = argv[2];
  filenames["outputdir"] = argv[3];
  
  NLOCI = atoi(argv[4]);

  int Niter = 100;
  int Nthin = 10;
  int Nburn = 100;
 
  if(argc > 5) 
    Niter = atoi(argv[5]);
  if(argc > 6)
    Nthin = atoi(argv[6]);
  if(argc > 7)
    Nburn = atoi(argv[7]);
  
  if(argc > 8)
    NSPECIES = atoi(argv[8]);

  if( NSPECIES > MAXSPECIES ) {
    cerr << "MAXSPECIES is insufficient: " << NSPECIES << " needed." << endl;
    exit(1);
  }
  
  if(NSPECIES == 1)
    UPDATENU = 0;

  ifstream input (filenames["input"].c_str());
  if(!input.is_open()) {
    error_and_exit("Failed to open " + filenames["input"]);
  }
  ifstream regionfile (filenames["regions"].c_str());
  if(!regionfile.is_open()) {
    error_and_exit("Failed to open " + filenames["regions"]);
  }

  // This one cannot be checked as the stream is used even if it
  // did not open!
  ifstream assignfile (filenames["assignfile"].c_str());

  string outputfilename = filenames["outputdir"] + "/Output_probs";
  ofstream output (outputfilename.c_str());

  string freqfilename = filenames["outputdir"] + "/Output_freqs";
  ofstream freqfile (freqfilename.c_str());

  string acceptfilename = filenames["outputdir"] + "/Output_accept";
  ofstream acceptfile (acceptfilename.c_str());

  string mapinfofilename = filenames["outputdir"] + "/Output_mapinfo";
  ofstream mapinfofile (mapinfofilename.c_str());

  string pifilename = filenames["outputdir"] + "/Output_pi";
  ofstream pifile;
  if(NSPECIES>1)
    pifile.open(pifilename.c_str());

  string corrfilename = filenames["outputdir"] + "/Output_corr";
  ofstream corrfile (corrfilename.c_str());

  string paramfilename = filenames["outputdir"] + "/Output_params";
  ofstream paramfile (paramfilename.c_str());
  ofstream Xfile (filenames["Xfile"].c_str());
  
  //assure ( input, "inputfile" );

  vector<int> Region;

  // END OF MENU PROCESSING PHASE

  // input genotype data from inputfile

  vector<int> RegionsPresent(0);
  
  vector<string> Id(0);

    // declare memory for Genotype and OriginalGenotype
  vector<vector<vector<int> > > Genotype;
  vector<vector<vector<int> > > OriginalGenotype;
   
   // Coding[i] is the actual allele that i codes for
  IntVec2d Coding(NLOCI,IntVec1d(MAXNALLELE,0.0));

  vector<int> NMissing(NIND,0);
  vector<int> Species;

  input_genotype_data(input,Region,Species,OriginalGenotype,NMissing,Id, RegionsPresent,true);
  if(ASSIGNFILE){
	      FIRSTSAMPLETOLOCATE = OriginalGenotype.size();
		  input_genotype_data(assignfile,Region,Species,OriginalGenotype,NMissing,Id,RegionsPresent,false);
    	  LASTSAMPLETOLOCATE = OriginalGenotype.size()-1; 
  }
   NIND = OriginalGenotype.size();
  if(!LOCATE){
	   FIRSTSAMPLETOLOCATE = 0;
	   LASTSAMPLETOLOCATE = (NIND-1);
  }
  
  if (ECHOINPUTS) {
    output_genotypes(OriginalGenotype,Id);
  }
  recode_genotypes(OriginalGenotype,Genotype,Coding,Nallele);
  if (ECHOINPUTS) {
    cout << "Number of Alleles at each locus:" << endl;
    int j=1;
    cout << "Locus : #alleles" << endl;
    for(vector<int>::iterator i = Nallele.begin(); i!=Nallele.end(); i++)
      cout << j++ << " : " << *i << endl;
  } 
  
  NREGION = RegionsPresent.size()+LOCATE;
  
  if (ECHOINPUTS) {
    cout << "Number of Regions: " << NREGION << endl;
  }
  
  vector<int> Perm(NREGION,0);

  for(int r = 0; r< NREGION; r++)
    Perm[r] = r;
  if(PERMUTE){
    rperm(Perm,NREGION-1);
    permute_regions(Region,Perm);
  }
  
  if(VERBOSE){
    cout << "Permutation:" << endl;
    for(int r = 0; r< NREGION; r++)
      cout << Perm[r] << " ";
    cout << endl;
  }

  vector<string> RegionName(NREGION);
  vector< vector<double> > Pi(NREGION,vector<double>(NSPECIES,1.0/NSPECIES)); // Pi[r][k] is proportion
  // of inds at location r that are in species k
  
  // declare memory for X and Y coords of region
  DoubleVec1d Xcoord(NREGION,0);
  DoubleVec1d Ycoord(NREGION,0);
  DoubleVec1d BoundaryX;
  DoubleVec1d BoundaryY;
  if(READBOUNDARY){
    cout << "Reading in Boundary data from boundary file " << filenames["boundaryfile"] << endl;
    ifstream bfile (filenames["boundaryfile"].c_str());
    if (!bfile.is_open()) {
      error_and_exit("Could not open boundary file");
    }
    ReadInBoundary(bfile,BoundaryX,BoundaryY);
  }

  
  Mapgrid mymapgrid;
  if(READGRID){
    cout << "Reading in gridfile data from file " << filenames["gridfile"] << endl;
    ifstream gridfile (filenames["gridfile"].c_str());
    if (!gridfile.is_open()) {
      error_and_exit("Could not open grid file" + filenames["gridfile"]);
    }
    mymapgrid.Initialize(gridfile, MARGIN);
  }

  // we presumably now know our map boundaries, so report on them
  if (READGRID) {
    mymapgrid.WriteMapInfo(mapinfofile);
  } else {
    mapinfofile << "Using default map information for ";
    if (SAVANNAHONLY) {
      mapinfofile << "Savannah elephants " << endl;
    } 
    if (FORESTONLY) {
      mapinfofile << "Forest elephants " << endl;
    }
    if (READBOUNDARY) {
      mapinfofile << "Custom boundary file " << endl;
    }
  }

  if(READGRID) mymapgrid.PrintGrid(mapinfofile);

  vector<int> SubRegion(NREGION,0);
  

  //declare memory for theta, and SumExpTheta
  // Count[r][k][l][j] is number of allele j in region r, species k, locus l
  // SumCount[r][k][l] is the sum of Count over j

  IntVec4d Count(NREGION,IntVec3d(NSPECIES,IntVec2d(NLOCI,IntVec1d(MAXNALLELE,0))));
  IntVec3d SumCount(NREGION,IntVec2d(NSPECIES,IntVec1d(NLOCI,0)));

  // Theta[r][k][l][j] is log relative freq of allele j in region r,
  // species k, locus l
  // ExpTheta is exponential of theta (stored to save computational time) 
  // SumExpTheta is sum of ExpTheta over j
  // So Freq[r][k][l][j] = ExpTheta[r][k][l][j]/SumExpTheta[r][k][l]
  // mu, X and L are such that Theta = mu + LX
  // where X are iid N(0,Alpha[0])
  // L is the (lower-triangular) Cholesky decomposition of the covariance of the thetas
  // mu is a priori N(0,1/beta), where beta is a priori Gamma(NBETA,LBETA)

   DoubleVec4d Theta(NREGION,DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,DoubleVec1d(MAXNALLELE,0.0)))); 
   DoubleVec4d ExpTheta(NREGION,DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,DoubleVec1d(MAXNALLELE,0.0)))); 
   DoubleVec3d SumExpTheta(NREGION,DoubleVec2d(NSPECIES,DoubleVec1d(NLOCI,0.0)));
   DoubleVec4d X(NREGION,DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,DoubleVec1d(MAXNALLELE,0.0))));

  DoubleVec2d Mu(NLOCI,DoubleVec1d(MAXNALLELE,0.0));
  DoubleVec3d Nu(NSPECIES,DoubleVec2d(NLOCI,DoubleVec1d(MAXNALLELE,0.0)));
  
  // Psi[r][k] is log relative abundance of species k in region r
  // ExpPsi and SumExpPsi as in Theta
  // So Pi[r][k] = ExpPsi[r][k]/SumExpPsi[r]
  // Lambda, Y and M are such that Psi = Lambda + MY
  // with Y being iid N(0,Delta[0])
  // and M being lower-triangular Cholesky decomp of covariance of Psi
  // which itself is determined by Delta
  // Lambda[k] are apriori iid N(1/eta), eta is a priori Gamma(NETA,LETA)

  DoubleVec2d Psi(NREGION,DoubleVec1d(NSPECIES,0.0));
  DoubleVec2d ExpPsi(NREGION,DoubleVec1d(NSPECIES,1.0));
  // The following line represents a bug fix. In the original code
  // this variable was initialized to NREGION, but this seems
  // incorrect.
  DoubleVec1d SumExpPsi(NREGION,NSPECIES);
  DoubleVec2d Y(NREGION,DoubleVec1d(NSPECIES,0.0));
  DoubleVec1d Lambda(NSPECIES,0.0);

  // these hold posterior means of Frequencies, X and X^2 in each region
  DoubleVec4d MeanFreq(NREGION, DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,
    DoubleVec1d(MAXNALLELE,0.0))));
  DoubleVec4d MeanX(NREGION, DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,
    DoubleVec1d(MAXNALLELE,0.0))));
  DoubleVec4d MeanX2(NREGION, DoubleVec3d(NSPECIES,DoubleVec2d(NLOCI,
    DoubleVec1d(MAXNALLELE,0.0))));

  // These are dimensioned by region1, species1, region2, species2
  // and represent covariances and correlations among regions
  DoubleVec4d MeanCov(NREGION, DoubleVec3d(NSPECIES, DoubleVec2d(NREGION,
    DoubleVec1d(NSPECIES,0.0))));
  DoubleVec4d MeanFittedCov(NREGION, DoubleVec3d(NSPECIES, DoubleVec2d(NREGION,
    DoubleVec1d(NSPECIES,0.0))));
  DoubleVec4d MeanCor(NREGION, DoubleVec3d(NSPECIES, DoubleVec2d(NREGION,
    DoubleVec1d(NSPECIES,0.0))));
  DoubleVec4d MeanFittedCor(NREGION, DoubleVec3d(NSPECIES, DoubleVec2d(NREGION,
    DoubleVec1d(NSPECIES,0.0))));

  vector<vector<double> > MeanPi(NREGION, vector<double>(NSPECIES,0.0)); // mean value of Pi

  // LogLik[r][k][l] holds the log-likelihood for individuals in region r, 
  // species k, at locus l 
  DoubleVec3d LogLik(NREGION,DoubleVec2d(NSPECIES,DoubleVec1d(NLOCI,0.0)));
 
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
        for(int j=0;j<MAXNALLELE;j++){	 
          Theta[r][k][l][j] = Mu[l][j] + Nu[k][l][j];
          ExpTheta[r][k][l][j] = exp(Theta[r][k][l][j]);
        }
      }
    }
  }
 
  DoubleVec1d Alpha(ALPHALENGTH,1); // parameters in covariance structure
  double Beta; // 1/variance of mus

  DoubleVec1d Delta(DELTALENGTH,1);
  double Eta;
  
  ALPHAMIN[0] = 0.2; // 1/ALPHA[0]
  ALPHAMAX[0] = 100; // controls amount of variation of frequencies about mean
  ALPHAMIN[1] = 1; // scale on which frequency correlations drop off
  ALPHAMAX[1] = 100000;
  ALPHAMIN[2] = 0.1;
  ALPHAMAX[2] = 2;
  if(ALPHALENGTH>3){
    ALPHAMIN[3] = 0;
    ALPHAMAX[3] = 100;
  }

  DELTAMIN[0] = 0.1;
  DELTAMAX[0] = 100;
  DELTAMIN[1] = 1;
  DELTAMAX[1] = 10000;
  DELTAMIN[2] = 0.1;
  DELTAMAX[2] = 2;

  DoubleVec1d Gamma (NSPECIES,0.3); // Gamma[k] is 1/variance of Nu[k]

  // input data from region file (Xcoords and Ycoords; regions and subregions)
  input_positions_data(regionfile,Xcoord,Ycoord,RegionName,SubRegion,Region,Perm,RegionsPresent);
  if (ECHOINPUTS) {
    output_positions_data(RegionName, Region, Xcoord, Ycoord, Id);
  }

  //declare memory for L, the matrix in the Cholesky decomp of Sigma
  // these are arrays for use in a C interface
  double * L = new double [NREGION * NREGION];
  double * M = new double [NREGION * NREGION];

  Initialise(X, Beta, Gamma, Alpha, Mu, Nu, Theta, L, Xcoord, Ycoord);

  // MeanProb[ind][r][k] is Pr of ind's genotype in region r, species k
  // LocusMeanProb [ind][r][k][l] is Pr of ind's genotype at locus l in region r, species k
  DoubleVec3d MeanProb(NIND, DoubleVec2d(NREGION, DoubleVec1d(NSPECIES, 0.0)));
  DoubleVec4d LocusMeanProb(NIND, DoubleVec3d(NREGION, DoubleVec2d(NSPECIES,
    DoubleVec1d(NLOCI, 0.0))));
  DoubleVec3d SubRegionProb(NIND, DoubleVec2d(NSUBREGION,
    DoubleVec1d(NSPECIES, 0.0)));

  // count up the number in each region
  count_up_alleles(Count,Region,Species,Genotype);
  if(VERBOSE)
    output_counts(Count,Perm);
  calc_SumCount(Count, SumCount);

  if(VERBOSE)
    output_empirical_freqs(RegionName,Perm,Coding,Count,SumCount);

  calc_ExpTheta_and_SumExpTheta(Theta,ExpTheta,SumExpTheta);
  // This one is definitely needed -- initialization
  calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
 
  double totaliter = 0; // number of iterations used in computing mean frequencies
 
  // start by doing burn-in (whether or not cross-validating)
  int templocate = LOCATE;
  LOCATE = 0; // don't locate during burnin

  int printdot = max(Nburn * SCREENPROGRESS_INTERVAL,1.0);
  int printincr = printdot;

  if(!LOCATEWHOLEREGION){
    // cerr << "Performing Global Burn-in iterations" << endl;
    cout << "Performing Global Burn-in iterations" << flush;
    for(int iter=0; iter< Nburn; iter++){
      for(int nthin=0; nthin < Nthin; nthin++){
	DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, BoundaryX, BoundaryY, mymapgrid);
      }
      
      OutputAcceptRates(acceptfile);

      // testing says that a call to calc_LogLik is not needed here
      // as DoAllUpdates adequately updates it
      //calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
      
      if(NSPECIES > 1)
	OutputPi( Pi, RegionName, pifile, Perm );
       
      // cerr << "Iteration:" << (iter+1) << "\033[A" << endl;
      if (iter == printdot) {
        printdot += printincr;
        cout << "." << flush;
      }

      OutputParameters(paramfile,Alpha,Beta,Gamma,Delta,Eta,Lambda,LogLik);

      if(OUTPUTX==1){
	for(int r=0; r < NREGION; r++){
	  for(int a =0; a<Nallele[0]; a++){
	    Xfile << X[r][0][0][a]  << " ";
	  }
	}
	Xfile << endl;
      }  
    }
  }
  cout << "." << endl;

  LOCATE = templocate;
   
  // different ways of running program: Locate (Continuous and smooth
  // assignment test, cross-validated); HYBRIDCHECK (undocumented);
  // Non-cross val, assignment test without cross-validation
  totaliter = 0;
	
  if(LOCATE){
    for(SAMPLETOLOCATE = FIRSTSAMPLETOLOCATE; SAMPLETOLOCATE <= LASTSAMPLETOLOCATE; SAMPLETOLOCATE++){
      
      string LOCATEFILE(filenames["outputdir"]);
	  LOCATEFILE.append("/");
      LOCATEFILE.append( Id[SAMPLETOLOCATE] );
      LOCACCEPT = 0;
      LOCATTEMPT = 0;

      ofstream locatefile (LOCATEFILE.c_str());
      
      // cerr << "Individual:" << (SAMPLETOLOCATE+1) << endl;      
      string outname = Id[SAMPLETOLOCATE];
      if (outname.length() > MAXOUTCHARS_INNAME) {
        int diff = outname.length() - MAXOUTCHARS_INNAME;
        outname.erase(MAXOUTCHARS_INNAME,diff);
      }
      // cout << "Individual: " << outname << endl;
      cout << "Individual: " << outname << flush;

      // reset theta, counts, etc, ignoring individual ind
      TRUEREGION = Region[SAMPLETOLOCATE];      
      Region[SAMPLETOLOCATE] = NREGION-1;
      if(LOCATEWHOLEREGION){ // locate all the inds from that region
	     for(int i = 0; i < NIND; i++){
	        if(Region[i] == TRUEREGION)
	           Region[i] = NREGION - 1;    
	     }
      }
      
      //cout << "Initialising XY position... ";
      //cerr << "Initialising XY position... ";
      InitialiseXY(BoundaryX,BoundaryY,Xcoord,Ycoord,mymapgrid);
      //cerr << "Done" << endl;

      calc_L(L,Alpha,Xcoord,Ycoord);
      count_up_alleles(Count,Region,Species,Genotype);
		
      if(REMOVEREGION){ // remove all individuals from the true region
	     for(int i = 0; i < NIND; i++){
	        if(Region[i] == TRUEREGION)
	            SubtractFromCount(i, Count, SumCount, Region, Species, Genotype);    
	     }
      }
  	  
      calc_SumCount(Count, SumCount);
      
      if(VERBOSE) output_empirical_freqs(RegionName,Perm,Coding,Count,SumCount);
    
      for(int r=0; r<NREGION; r++){
	for(int k=0; k<NSPECIES; k++){
	  for(int l=0; l<NLOCI; l++){
	    for(int j=0;j<Nallele[l];j++){  
	      X[r][k][l][j] = 0;
	    }
	  }
	}
      }

      InitialiseTheta(Theta,X,Mu,Nu,L);
      calc_ExpTheta_and_SumExpTheta(Theta,ExpTheta,SumExpTheta);
      // This one is definitely needed
      calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
       
      printdot = max(Nburn * SCREENPROGRESS_INTERVAL,1.0);
      printincr = printdot;
      cout << "   BurnIn" << flush;

      for(int iter=0; iter< Nburn; iter++){
	// cout << "Burnin Iteration:" << (iter+1) << "\033[A" << endl;
        if (iter == printdot) {
          cout << "." << flush;
          printdot += printincr;
        }
	for(int nthin=0; nthin < Nthin; nthin++){
	  DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, BoundaryX, BoundaryY, mymapgrid);
	}
	// this should not be needed due to call in update_Locations
        // calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);

 	//CheckSumExpTheta(ExpTheta,SumExpTheta);
	double sumloglik = SumLogLik(LogLik);
	OutputLatLongs(locatefile,Xcoord[NREGION-1],Ycoord[NREGION-1],sumloglik);

      }
      
      // cout << endl;
      cout << "." << "Searching" << flush;

      printdot = max(Niter * SCREENPROGRESS_INTERVAL,1.0);
      printincr = printdot;

      for(int iter=0; iter< Niter; iter++){
	// cout << "Iteration:" << (iter+1) << "\033[A" << endl;
        if (iter == printdot) {
          cout << "." << flush;
          printdot += printincr;
        }
	for(int nthin=0; nthin < Nthin; nthin++){
	  DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, BoundaryX, BoundaryY, mymapgrid);
	}
        // This should not be needed due to call in update_Locations
        // calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
	
	double sumloglik = SumLogLik(LogLik);
        assert(InRange(Xcoord[NREGION-1],Ycoord[NREGION-1],BoundaryX,BoundaryY,mymapgrid));
	OutputLatLongs(locatefile,Xcoord[NREGION-1],Ycoord[NREGION-1],sumloglik);

	UpdateMeans(ExpTheta, X, Pi, SumExpTheta, MeanFreq, MeanX, MeanX2, MeanPi, MeanCov, MeanFittedCov, MeanCor, MeanFittedCor, Alpha, Xcoord, Ycoord,Theta,Mu,Nu);	
	UpdateLocusMeanProb(SAMPLETOLOCATE, ExpTheta, SumExpTheta, LocusMeanProb, Genotype);
	
	totaliter++;
      }	
      cout << "." << endl;
      
      locatefile << "Acceptance rate: " << LOCACCEPT*1.0/LOCATTEMPT  << endl;

      locatefile.close();
      
      if(REMOVEREGION){ // add back all individuals from the true region
	for(int i = 0; i < NIND; i++){
	  if(Region[i] == TRUEREGION)
	    AddToCount(i, Count, SumCount, Region, Species, Genotype);    
	}
      }
      Region[SAMPLETOLOCATE] = TRUEREGION; // reset to original region
      if(LOCATEWHOLEREGION){
	for(int i = 0; i < NIND; i++){
	  if(Region[i] == NREGION-1)
	    Region[i] = TRUEREGION;
	}
      }
      // cerr << endl;

    }
  } else if(HYBRIDCHECK) {
    
    string outputfilename = filenames["outputdir"] + "/Output_hybrid";
    ofstream hybridout (outputfilename.c_str());

    double loglik = 0;
    
    hybridout << "Hybridcheck =" << HYBRIDCHECK << "\n";

    if(HYBRIDCHECK <10){ // set regions to be subregion
      for(int ind = 0 ; ind<NIND; ind++){
	if(Region[ind]>=0)
	  Region [ind] = SubRegion[Region[ind]];
      }
      NREGION = NSUBREGION;
    }

    // column headers for tab delimited output file
    hybridout << "Individual";
    for(int r = 0; r < NREGION; r++){
      for(int s = 0; s< NREGION; s++){
        if(HYBRIDCHECK==1 || HYBRIDCHECK == 11 || r==s) {
          hybridout << "\t" << r << "/" << s;
        }
      }
    }

    hybridout << "\n";
    count_up_alleles(Count,Region,Species,Genotype);

    calc_SumCount(Count, SumCount); 
    if (VERBOSE)
      output_counts(Count, Perm);
    hybridout.precision(4);

    for(int ind = 0; ind < NIND; ind++){ 
      
      SubtractFromCount(ind, Count, SumCount, Region, Species, Genotype); 
      calc_SumCount(Count, SumCount);

      hybridout << Id[ind];

      for(int r = 0; r < NREGION; r++){
	for(int s = 0; s< NREGION; s++){
	  if(HYBRIDCHECK==1 || HYBRIDCHECK == 11 || r==s) {
            double lhprob = log_hybrid_Prob(Count,SumCount,Genotype,ind,r,s);
            hybridout << "\t" << lhprob;
          }
	  if((r == Region[ind]) && (r==s))
	    loglik += log_hybrid_Prob(Count,SumCount,Genotype,ind,r,s);
	}
      }
      hybridout << "\n";

      AddToCount(ind, Count, SumCount, Region, Species, Genotype);   
      calc_SumCount(Count, SumCount);
    }
    hybridout.close();

  } else { // if not cross-validating

    XACCEPT = 0;
    XATTEMPT = 0;
    // cerr << "Performing Main Iterations " << endl;
    cout << "Performing Main Iterations" << flush;
    printdot = max(Niter * SCREENPROGRESS_INTERVAL,1.0);
    printincr = printdot;

    for(int iter=0; iter< Niter; iter++){     
      for(int nthin=0; nthin < Nthin; nthin++){
	DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, BoundaryX, BoundaryY, mymapgrid);
      }

      // Test if this one is needed  -- not tested as not in execution path for test data
      // It is probably not needed but retaining it to be safe
#ifndef NDEBUG
      DoubleVec3d debug_LogLik(LogLik);
#endif  
      calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
#ifndef NDEBUG
      if(!CompareLogLiks(LogLik,debug_LogLik)) cerr << "At line 3204" << endl;
#endif

      if(NSPECIES > 1)
	OutputPi( Pi, RegionName, pifile, Perm );
      
      UpdateMeans(ExpTheta, X, Pi, SumExpTheta, MeanFreq, MeanX, MeanX2, MeanPi, MeanCov, MeanFittedCov, MeanCor, MeanFittedCor, Alpha, Xcoord, Ycoord, Theta, Mu, Nu);
      totaliter++;
      UpdateLocusMeanProb(ExpTheta, SumExpTheta, LocusMeanProb, Genotype);
      // DEBUG  I am removing this as I believe it is lava flow
      // if(USESUBREGION==1){
	// UpdateSubRegionProb(ExpTheta, SumExpTheta, SubRegionProb, Genotype, SubRegion);
      // }
      //OutputEstimatedFreqs( ExpTheta, SumExpTheta, output, Coding, Perm);   
      //OutputTheta( MeanX, output, Coding, Perm);
      //cerr << "Iteration:" << (iter+1) << "\033[A" << endl;
      if (iter == printdot) {
        cout << "." << flush;
        printdot += printincr;
      }

      OutputParameters(paramfile, Alpha, Beta, Gamma, Delta, Eta, Lambda, LogLik);


      if(OUTPUTX==1){
	for(int r=0; r < NREGION; r++){
	  for(int a =0; a<Nallele[0]; a++){
	    Xfile << X[r][0][0][a] << " ";
	  }
	}
	Xfile << endl;
      }
    }
    // cerr << endl;
    cout << "." << endl;
  }

  NormaliseMeanFreq(MeanFreq,totaliter);
  ComputeLogMeanProb(MeanProb,LocusMeanProb,totaliter); 
  NormaliseMeanPi(MeanPi);
  NormaliseMeanCov(MeanCov, MeanFittedCov, totaliter);
  NormaliseMeanCov(MeanCor, MeanFittedCor, totaliter);

   // output probs of each individual being in each region
  OutputLogMeanProb(MeanProb,output,NMissing,Region,Id, RegionName, Perm, FIRSTSAMPLETOLOCATE, LASTSAMPLETOLOCATE, RegionsPresent);
  
  // output probs of each individual being in each region, according to mean freqs
  OutputLogMeanProb2(MeanFreq,Genotype,output,NMissing,Region,Perm);

  // DEBUG I am removing this because I believe it is lava flow
  //if(USESUBREGION==1){
   // NormaliseSubRegionProb(SubRegionProb);
  //}

  // output correlations and covariances (fitted and empirical) to corrfile

  corrfile << "location1 location2 distance cov fittedcov corr fittedcorr" << endl;

     for(int r0=0; r0<NREGION; r0++){
      for(int k0=0; k0<NSPECIES; k0++){
	for(int r1=0; r1<=r0; r1++){
	  for(int k1=0; k1<=k0; k1++){
	    corrfile << RegionName[r0] << " " << RegionName[r1] << " " << Distance(r0,r1,Xcoord,Ycoord) << " " << MeanCov[r0][k0][r1][k1] << " " << MeanFittedCov[r0][k0][r1][k1]  << " " << MeanCor[r0][k0][r1][k1] << " " << MeanFittedCor[r0][k0][r1][k1]<< " " << endl;
	  }
	}
      }
    }
    
  if(NSPECIES > 1){
    OutputPi(MeanPi, RegionName, pifile, Perm );
  }

  // 
  // Output Posterior Mean Frequency Estimates
  //

  OutputMeanFreq(freqfile,RegionName,Perm,Coding,MeanFreq);

  freqfile << "Empirical Freqs:" << endl;
  for(int l=0;l<NLOCI; l++){
    freqfile << "Locus " << (l+1) << endl;
    freqfile << setw(5) << " " << " : ";
    OutputRegionNames(freqfile, RegionName, Perm);
    for(int allele=0; allele<Nallele[l]; allele++){
      freqfile << setw(5) << Coding[l][allele] << " : ";
      for(int r=0;r<NREGION- LOCATE;r++){     
	freqfile << std::fixed << setprecision(3) << setw(10) << (1.0*Count[Perm[r]][0][l][allele])/SumCount[Perm[r]][0][l] << " ";
      }
      freqfile << endl;
    }
  }
 
cout << "Program finished" << endl;
}
