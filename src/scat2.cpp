//#include "input.hpp"
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

extern "C" void dpotrf_(
	const char &uplo,		// (input)
	const int &n,			// (input)
	double *a,			// a[n][lda] (input/output)
	const int &lda,			// (input)
	int &info			// (output)
	);

using namespace std; 
const double PI = 3.141592; 
    
int USESUBREGION = 0; // indicates whether region file also holds subregion data
int NSUBREGION = 6;

int DCSPECIAL = 0;
int PSEUDOCOUNT = 0;
int LOCATE = 0; // whether to try to estimate the position of a particular sample
int ASSIGNFILE = 0;
int SAMPLETOLOCATE = 0;
int FIRSTSAMPLETOLOCATE=0;
int LASTSAMPLETOLOCATE=-1;
int LOCATEWHOLEREGION = 0;
int INCLUDENUGGET = 0;
int SKIPCOL = 0;

int NUMBERREGIONSFROM = 1; // what's the first region numbered? see T flag
int READBOUNDARY = 0;

int PERMUTE = 1;
int VERBOSE = 0;

char TAG = 'a'; // tag added to output files
int OFFSET = 1; // start numbering output files from this number
int OUTPUTSAMPLEFREQ = 0;
double TEMPERATURE = 1; //0.0001;
double ALPHAUPDATESD = 0.4;

double DELTA = 0.05; //0.05; //0.1; // prob of genotyping error in this one sample (set to 0 for no error)
double NULLPROB = 0; // 0.1;// prob of null allele
    
int TRUEREGION; // stores true region of the sample to be located
int CHEAT = 0; // stores whether to cheat by starting estimates at the true location
int NONUNIFORMPRIOR = 0; // non-uniform prior weights sample as being more 
// likely to be near one of the sampling locations
int REMOVEREGION = 0;

int USELANGEVIN = 0;
int USESPATIAL =1; 
int CROSSVAL = 0;
int HYBRIDCHECK = 0;

int STARTCROSSVAL = 0; // ind to start CV at
int ENDCROSSVAL = 0; // ind to end CV at

int UPDATEALPHA =1;
int INITALPHA = 0;
int UPDATEMU = 1;
int UPDATENU = 1;
int UPDATEBETA = 1;
int UPDATEX = 1;

int FORESTONLY = 0;
int SAVANNAHONLY = 0;

int OUTPUTX = 0;
int UPDATEJOINT = 0;

double XPROPOSALFACTOR = 0.5; // factor by which to multiply sd
double YPROPOSALFACTOR = 0.5;

int NIND;
const int MAXNALLELE = 60; 
int NREGION;
int NSPECIES = 1;

int NLOCI = 16;
const double EPSILON = 1e-100;

double NBETA = 0.001; // parameters on gamma prior on beta 
double LBETA = 0.001; // (beta is prior precision of mu)

double NETA = 0.001; // params of gamma prior on eta,
double LETA = 0.001; // (eta is prior precision of lambda) 

double ALPHA0 =1; // values at which to initialise Alpha
double ALPHA1 = 1;
double ALPHA2 = 1;
double BETA = 1;

const int ALPHALENGTH = 3+INCLUDENUGGET;
vector<double> NA(ALPHALENGTH-1,0.001); // gamma prior on Alpha[0] and Alpha[1]
vector<double> LA(ALPHALENGTH-1,0.001); // (Gamma(NA,LA)) 
vector<double> ALPHAMAX(ALPHALENGTH,0);
vector<double> ALPHAMIN(ALPHALENGTH,0);
vector<int> ALPHAATTEMPT(ALPHALENGTH,0);
vector<int> ALPHAACCEPT(ALPHALENGTH,0);

const int DELTALENGTH = 3;
vector<double> ND(DELTALENGTH-1,0.001); // gamma prior on Delta[0] and Delta[1]
vector<double> LD(DELTALENGTH-1,0.001); // (Delta(ND,LD)) 
vector<double> DELTAMAX(DELTALENGTH,0);
vector<double> DELTAMIN(DELTALENGTH,0);
vector<int> DELTAATTEMPT(DELTALENGTH,0);
vector<int> DELTAACCEPT(DELTALENGTH,0);


int XACCEPT =0;
int XATTEMPT = 0;

int MUACCEPT = 0;
int MUATTEMPT =0;

int NUACCEPT = 0;
int NUATTEMPT =0;

int YACCEPT =0;
int YATTEMPT = 0;

int LAMBDAACCEPT = 0;
int LAMBDAATTEMPT =0;

int LOCATTEMPT = 0;
int LOCACCEPT = 0;

vector<int> Nallele;


double isLeft(double x0, double y0, double x1, double y1,  double x2, double y2){
  return( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
}

//  winding number test for a point in a polygon
//  modelled on code from softsurfer.com, by Dan Sunday
//      Input:   x,y = a point,
//               BoundaryX and BoundaryY = points of a polygon with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if (x,y) is outside polygon)
int IsInsideBoundary( double x, double y, vector<double> & BoundaryX, vector<double> & BoundaryY)
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

void InitialiseXY(vector<double> & BoundaryX, vector<double> & BoundaryY, vector<double> & Xcoord, vector<double> & Ycoord)
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
    do{
      double lambda = ranf();
      int r = (int) (ranf() * (NREGION-1));
      Xcoord[NREGION-1] = lambda * Xcenter + (1-lambda) * Xcoord[r];
      Ycoord[NREGION-1] = lambda * Ycenter + (1-lambda) * Ycoord[r];
    } while(IsInsideBoundary(Xcoord[NREGION-1],Ycoord[NREGION-1],BoundaryX, BoundaryY)==0);
  }
}


void ReadInBoundary(ifstream & bfile, vector<double> & BoundaryX, vector<double> & BoundaryY)
{
   double x,y;
   do{
      bfile >> y;
      bfile >> x;
      cerr << y << "," << x << endl;
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
  
  string line;
  line.assign(getline(pbuf));
  cout << "READING GENOTYPE DATA" << endl;
  while(line.size()>0){
	for(int chrom = 0; chrom<2; chrom++){
  		int beg = line.find_first_not_of(delimit, 0);
    	int end = line.find_first_of(delimit, beg);
    	if(end == beg) break;
    	string sv(line, beg, end-beg);
    	if(chrom==0){
			Id.push_back(sv.data());
			NMissing.push_back(0);
			vector<vector<int> > blank(2,vector<int>(NLOCI,0));
			Genotype.push_back(blank);
		   	cout << "Inputting data for " << Id[Id.size()-1] << endl;
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
	
			//input >> s;
			//Species.push_back(s);
			Species.push_back((int) (NSPECIES * (1-ranf())));
		}
		
		int dummy;
    	for(int skip = 0; skip<SKIPCOL; skip++){
			beg = line.find_first_not_of(delimit, end);
			end = line.find_first_of(delimit, beg);
	  	}

      	for(int j=0; j<NLOCI; j++){
      		int allele;
      		int i = Genotype.size()-1;
			beg = line.find_first_not_of(delimit, end);
	  		end = line.find_first_of(delimit, beg);
	  		string sv(line, beg, end-beg);
	  		allele = atoi(sv.data());  			  
	 
	  		Genotype[i][chrom][j] = allele;
	  		if((allele<0) && (chrom ==0))
	     		NMissing[i]++;
	  
	  	}
        line.assign(getline(pbuf));
	}
  }
}

void OutputLatLongs(ostream & locatefile, double x, double y, double loglik){
  // convert x and y into lat and long
  locatefile << 180 * y/PI << " " << 180 * x/PI <<  " " << loglik << endl; 
}

void OutputRegionNames(ostream & freqfile, vector<string> & RegionName, vector<int> & Perm)
{
  for(int r=0;r<NREGION;r++){     
    freqfile << setiosflags(ios::fixed) << setw(9-RegionName[Perm[r]].size()) << RegionName[Perm[r]] << " ";
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

void output_positions_data( vector<string> & RegionName, vector<int> & Region, vector<double> & x, vector<double> & y, vector<string> & Id){

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
    cout << "Location " << region << ", Name " << regionname << endl;
    
    region = GetLocationNumber(RegionsPresent, region);
    
    //    region -=NUMBERREGIONSFROM;

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
       cout << "Region " << region << ", Name " << regionname << endl;
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

void output_genotypes(vector<vector<vector<int> > > & Genotype,vector<string> & Id)
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


void recode_genotypes(vector<vector<vector<int> > > & OriginalGenotype, vector<vector<vector<int> > > & RecodedGenotype, int ** Coding, vector<int> & Nallele)
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

    for(int allele =0; allele<AllelesPresent.size(); allele++){
      //cout << "Locus:" << locus << ",Allele:" << AllelesPresent[allele] << endl;
      for(int ind =0; ind<NIND; ind++){
	for(int chrom=0; chrom<2; chrom++){
	  if(OriginalGenotype[ind][chrom][locus]==AllelesPresent[allele]){
	    RecodedGenotype[ind][chrom][locus] = allele;
	    Coding[locus][allele] = AllelesPresent[allele];
	    //cout << "allele:" << allele << endl;
	  }
	}
      }
    }
    
  }


}

void SubtractFromCount(int ind, int **** Count, int *** SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector<int> > > & Genotype)
{
  //cout << "Subtracting Ind " << ind << " from region " << Region[ind] << "," << Species[ind] << endl;
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


void AddToCount(int ind, int **** Count, int *** SumCount, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype)
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


void count_up_alleles(int **** Count, vector<int> & Region, vector<int> & Species, vector<vector<vector< int> > > & Genotype)
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

void calc_SumCount(int **** Count, int *** SumCount)
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
double simple_Prob(int **** Count, int *** SumCount,vector<vector<vector<int> > > &  Genotype, int ** Coding, int ind, int r)
{
  int k=0;
  double prob = 1;
  for(int l=0;l<NLOCI; l++){
    for(int chr =0; chr <2; chr++){
      if(Genotype[ind][chr][l]>=0){
	//cout << "allele: " << Coding[l][Genotype[ind][chr][l]] << ",";
	 double p=0;
	 if(SumCount[r][k][l]>0)
	   p = (Count[r][k][l][Genotype[ind][chr][l]]*1.0+1.0)/(SumCount[r][k][l]+MAXNALLELE);
	 if(p==0)
	   p = 1.0/ (SumCount[r][k][l]+1);
	 cout << p << endl;
	 prob *= p;
      }	
    }	

    //    if(Genotype[ind][0][l]>=0)
    //  if(Genotype[ind][0][l] != Genotype[ind][1][l])
    //	prob*=2;

  }
  return prob;  
}


// compute prob of individual's genotype as hybrid of region r and s
double log_hybrid_Prob(int **** Count, int *** SumCount,vector<vector<vector< int> > > &  Genotype, int ** Coding, int ind, int r, int s)
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
    //cout << llocusprob << " ";
  } 
  //cout << endl;
  return lprob;

}


double log_hybrid_Prob(double **** Freq,vector<vector<vector<int> > > & Genotype, int ind, int r, int s)
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
    //cout << llocusprob << " ";
  } 
  //cout << endl;
  return lprob;

}

// void calc_SumExpTheta(double **** Theta, double *** SumExpTheta)
// {
//   for(int r=0; r<NREGION; r++){
//     for(int k=0; k<NSPECIES; k++){
//       for(int l=0; l<NLOCI; l++){
// 	SumExpTheta[r][k][l] = 0;
// 	for(int j=0;j<Nallele[l];j++){
// 	  SumExpTheta[r][k][l]+=exp(Theta[r][k][l][j]);
// 	}
//       }
//     }
//   }
// }

void calc_ExpTheta_and_SumExpTheta(double **** Theta, double **** ExpTheta, double *** SumExpTheta)
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

void output_counts(int **** Count, vector<int> & Perm)
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
double Covariance(int r0, int k0, int r1, int k1, double **** Theta, double ** Mu, double *** Nu){

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
double Correlation(int r0, int k0, int r1, int k1, double **** Theta, double ** Mu, double *** Nu){

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
double Covariance(int locus, int r0, int k0, int r1, int k1, double **** Theta, double ** Mu, double *** Nu){

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
double Correlation(int l, int r0, int k0, int r1, int k1, double **** Theta, double ** Mu, double *** Nu){

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

// Cartesian distance
//  double Distance(double x1, double y1, double x2, double y2){
//    return sqrt( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
//  }

// distance between region r and region s
double Distance(int r, int s, vector<double> & Xcoord, vector<double> & Ycoord){
  return Distance(Xcoord[r], Ycoord[r], Xcoord[s], Ycoord[s]);
}


// compute the matrix L (uses the Lapack routine dpotrf)
void calc_L(double * L,vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord)
{
  //cout << "Computing L" << endl;
  //cout << "Sigma value" << endl;
  for(int r=0; r<NREGION; r++){
    for(int s=0; s <NREGION; s++){
      L[r+ s*NREGION] = FittedCovariance(Alpha,Distance(r,s,Xcoord,Ycoord));
      //cout << setprecision(10) << L[r+s*NREGION] << " ";
    }
    //cout << endl;
  }

  int INFO = 0;
  char UPLO = 'L';
  dpotrf_(UPLO,NREGION,L,NREGION,INFO);

  //cout << "Computed value of L" << endl;
  // for(int r=0; r<NREGION; r++){
  //   for(int s=0; s <NREGION; s++){
       //cout << L[r+s*NREGION] << " ";
  //   }
     //cout << endl;
  // }




  // INFO returns something about whether the computation was successful
  if(INFO>0)
    cerr << "Warning: INFO=" << INFO << endl;
}
 

double CurrentLogLik(double *** LogLik, int r){
  double sum;
  for(int l=0; l<NLOCI; l++){
	  sum += LogLik[r][0][l];
  }
  return sum;


}
// compute the Loglikelihood
// LogLik[r][k][l] holds the loglikelihood for individuals in region r, 
// species k, at locus l
void calc_LogLik(double *** LogLik,double **** Theta, double *** SumExpTheta, int **** Count, int *** SumCount){
  double loglik = 0;
  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	//cout << r << " " << k << " " << l << " " << endl;
	//cout << "SumCount:" << SumCount[r][k][l] << endl;
	//cout << "SumExpTheta:" << SumExpTheta[r][k][l] << endl;
	LogLik[r][k][l] = 0;
	//if(SumCount[r][k][l] > 0)
	LogLik[r][k][l] -= SumCount[r][k][l] * log(SumExpTheta[r][k][l]);
	  //else{
	  //cout << "HERE:" << SumCount[r][k][l] << endl;
	  //}
	for(int j=0; j<Nallele[l]; j++){
	  //  cout << Count[r][k][l][j] << endl;
	  LogLik[r][k][l] += Count[r][k][l][j] * Theta[r][k][l][j];	  
	}
	//cout << "LogLik:" <<  LogLik[r][k][l] << endl;
      }
    }
  }
}


//used in testing: compute the Loglikelihood but with the pseudocounts, if any, removed
double calc_LogLikWithoutPseudoCounts(double **** Theta, double *** SumExpTheta, int **** Count, int *** SumCount){
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
double SumLogLik(double *** LogLik){
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

void NormaliseMeanCov(double **** MeanCov, double **** MeanFittedCov, double totaliter)
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


void CheckSumExpTheta(double **** ExpTheta, double *** SumExpTheta)
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
void UpdateMeans(double **** ExpTheta,double **** X, vector<vector<double> > & Pi, double *** SumExpTheta, double **** MeanFreq, double **** MeanX, double **** MeanX2, vector<vector<double> > & MeanPi, double **** MeanCov, double **** MeanFittedCov, double **** MeanCor, double **** MeanFittedCor, vector<double> & Alpha, vector<double> & Xcoord, vector<double> & Ycoord, double **** Theta, double ** Mu, double *** Nu)
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
void UpdateSubRegionProb(double **** ExpTheta, double *** SumExpTheta, double *** SubRegionProb, vector<vector<vector<int> > > & Genotype, vector<int> & SubRegion)
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
void UpdateLocusMeanProb(int ind, double **** ExpTheta, double *** SumExpTheta, double **** LocusMeanProb, vector<vector<vector<int> > > & Genotype)
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
 	  
	// if(Genotype[ind][chr][l]>=0){
// 	    double p = (1-DELTA) * ExpTheta[r][k][l][Genotype[ind][chr][l]]/SumExpTheta[r][k][l] + DELTA/Nallele[l];
// 	    if(p>1) p=1; // allow for rounding errors
// 	    prob *= p;
// 	    llocusprob += log(p);
// 	  }	
// 	}	
	LocusMeanProb[ind][r][k][l] += exp(llocusprob);
      }
    }
  }
}

// compute probs for all individuals
void UpdateLocusMeanProb(double **** ExpTheta, double *** SumExpTheta, double **** LocusMeanProb, vector<vector<vector<int> > > & Genotype)
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

void NormaliseMeanFreq(double **** MeanFreq,double totaliter){
  for(int l=0;l<NLOCI; l++){
    for(int r=0;r<NREGION;r++){     
      for(int allele=0; allele<Nallele[l]; allele++){
	MeanFreq[r][0][l][allele]/=totaliter;
      }
    }
  } 
}

// void NormaliseMeanProb(double *** MeanProb){
//   for(int ind = 0; ind < NIND; ind++){
//     double sum = 0;
//     int numregion = NREGION;
//     if(LOCATE) numregion --;
//     for(int r=0;r<numregion;r++){
//       for(int k=0; k<NSPECIES; k++){
// 	sum += MeanProb[ind][r][k];
//       }
//     }
//     for(int r=0;r<NREGION;r++){
//       for(int k=0; k<NSPECIES; k++){
// 	MeanProb[ind][r][k] /= sum; 
//       }
//     }
//   }
// }

// void CorrectMeanProb(double **** MeanProb,double totaliter){
//   for(int ind = 0; ind < NIND; ind++){
//     for(int r=0;r<NREGION;r++){
//       for(int k=0; k<NSPECIES; k++){
// 	for(int l =0; l<NLOCI; l++){
// 	  MeanProb[ind][r][k][l] /= totaliter;
// 	}
//       }
//     }
//   }
// }

void ComputeMeanProb(double *** MeanProb,double **** LocusMeanProb, double totaliter){
  for(int ind = 0; ind < NIND; ind++){
    for(int r=0;r<NREGION;r++){
      for(int k=0; k<NSPECIES; k++){
	MeanProb[ind][r][k] = 1;
	for(int l =0; l<NLOCI; l++){
	  MeanProb[ind][r][k] *= (LocusMeanProb[ind][r][k][l] / totaliter);
	}
      }
    }
  }
}

void NormaliseSubRegionProb(double *** SubRegionProb){
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
  
//   for(int k = 0; k<NSPECIES; k++){
//     for(int r = 0; r<Pi.size(); r++){
//       ostr << RegionName[r] << "(" << k << ") ";
//     }
//   }
//   ostr << endl;
   
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
  
void OutputEstimatedFreqs(double **** ExpTheta, double *** SumExpTheta, ostream & output, int ** Coding, vector<int> & Perm)
{
  for(int l=0;l<NLOCI; l++){
    for(int allele=0; allele<Nallele[l]; allele++){
      output << setw(5) << Coding[l][allele] << " : ";
      for(int r=0;r<NREGION;r++){     
	//output << Theta[r][0][0][0] << " " << Theta[r][0][0][1] << " ";
	output << setiosflags(ios::fixed) << setprecision(3) << setw(5) << ExpTheta[Perm[r]][0][l][allele]/SumExpTheta[Perm[r]][0][l] << " ";
      }
      output << endl;
    }
  }

}

void OutputMeanPos(double *** MeanProb, vector<double> & Xcoord, vector<double> & Ycoord, vector<int> & Region, vector<string> & RegionName, ostream & output, vector<int> & NMissing, vector<int> & Perm)
{
  double meanx, meany, truex, truey, sumprob;
  
  for(int ind=0; ind<NIND; ind++){
    if(Region[ind]>=0){
      truex = Xcoord[Region[ind]];
      truey = Ycoord[Region[ind]];
   
      meanx=0;
      meany=0;
      sumprob=0;
      int numregion = NREGION;
      if(LOCATE) numregion --; // if LOCATE is true, ignore last region (which contains only the ind to be located)
      
      
      for(int r=0;r<numregion;r++){
	meanx+= MeanProb[ind][r][0] * Xcoord[r];
	meany+= MeanProb[ind][r][0] * Ycoord[r];
	sumprob += MeanProb[ind][r][0];
    }    
      meanx /= sumprob;
      meany /= sumprob;
      output << ind << " " << NLOCI - NMissing[ind] << " " << RegionName[Region[ind]] << " " << setiosflags(ios::fixed) << setprecision(3) << setw(5) << meanx << " " << meany << " " << truex << " " << truey << " " <<  Distance(meanx,meany,truex,truey) << endl;
      
    }
  }
}

void OutputMeanProb(double *** MeanProb, ostream & output, vector<int> & NMissing, vector<int> & Region, vector<int> & Perm)
{
  for(int ind=0; ind<NIND; ind++){
    output << NMissing[ind] << " " << Region[ind] << " ";
    for(int r=0;r<NREGION;r++){     
      for(int k=0; k<NSPECIES; k++){
	//output << setiosflags(ios::fixed) << setprecision(3) << setw(6) << MeanProb[ind][r][k];
	output << setiosflags(ios::scientific) <<  MeanProb[ind][Perm[r]][k] << " ";
      }
    }
    output << endl;   
  }
}

void OutputLogMeanProb(double *** MeanProb, ofstream & output, vector<int> & NMissing, vector<int> & Region, vector<string> & Id, vector<string> & RegionName, vector<int> & Perm, int first, int last, vector<int> & RegionsPresent)
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
      //r += NUMBERREGIONSFROM;
    }
    
    output << Id[ind] << " " << (NLOCI - NMissing[ind]) << " " << ((r>=0) ? RegionsPresent[r]:r) << " ";

    int assignedr = 0;
    int assignedk = 0;
    double maxlogprob = log(MeanProb[ind][Perm[0]][0]);
    
    for(int k=0; k<NSPECIES; k++){
      for(int r=0;r<(NREGION-LOCATE );r++){    
	if(log(MeanProb[ind][Perm[r]][k]) > maxlogprob){
	  assignedr = r; assignedk = k; maxlogprob = log(MeanProb[ind][Perm[r]][k]);
	}
      }
    }
			    
    output << RegionsPresent[assignedr] << " ";
    if(NSPECIES>1)
      output << assignedk << " ";

    for(int k=0; k<NSPECIES; k++){
      for(int r=0;r<(NREGION-LOCATE);r++){     
	output << setiosflags(ios::scientific) <<  log(MeanProb[ind][Perm[r]][k]) << " ";
      }
    }
    output << endl;
  }
}


void OutputLogMeanProb2(double **** MeanFreq, vector<vector<vector<int> > > & Genotype, ostream & output, vector<int> & NMissing, vector<int> & Region, vector<int> & Perm)
{
  for(int ind=0; ind<NIND; ind++){
    output << NMissing[ind] << " " << Region[ind] << " ";
    for(int r=0;r<NREGION;r++){     
      for(int k=0; k<NSPECIES; k++){
	//output << setiosflags(ios::fixed) << setprecision(3) << setw(6) << MeanProb[ind][r][k];
	output << setiosflags(ios::scientific) << log_hybrid_Prob(MeanFreq,Genotype, ind, Perm[r], Perm[r]) << " ";
      }
    }
    output << endl;
  }
}


void OutputParameters(ofstream & paramfile, vector<double> & Alpha, double & Beta, vector<double> & Gamma, vector<double> & Delta, double & Eta,double * Lambda, double *** LogLik)
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
	//output << Theta[r][0][0][0] << " " << Theta[r][0][0][1] << " ";
	output << setiosflags(ios::fixed) << setprecision(3) << setw(5) << Theta[Perm[r]][0][l][allele] << " ";
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
void update_Species(vector<int> & Species, vector< vector<double> > & Pi, vector<int> & Region, vector<vector<vector<int> > > & Genotype, double **** ExpTheta, double *** SumExpTheta)
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


void compute_Pi(vector< vector<double> > & Pi, double ** ExpPsi, double * SumExpPsi )
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
void update_Nu(double *** Nu, vector<double> & Gamma, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount,
vector <vector<double> > & NewTheta,vector <vector<double> > & NewExpTheta,vector <vector<double> > & NewSumExpTheta,vector <vector<double> > & NewLogLik )
{

  double NewNu;

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
	
	//cout << "A=" << A << endl;
	
	NUATTEMPT +=1;
	if( ranf()<A ){ //accept move 
	  NUACCEPT +=1;
	  for(int r=0; r<NREGION; r++){
	    Theta[r][k][l][j] = NewTheta[r][k];
	    ExpTheta[r][k][l][j] = NewExpTheta[r][k];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r][k];
	    LogLik[r][k][l] = NewLogLik[r][k];
	    //cout << "Accept!" << endl;	    
	  }
	  Nu[k][l][j] = NewNu;
	}
      }
    }
  }
}


// update Mu (the background "ancestral" allele freqs)
void update_Mu(double ** Mu, double Beta, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount,
vector <vector<double> > & NewTheta,vector <vector<double> > & NewExpTheta,vector <vector<double> > & NewSumExpTheta,vector <vector<double> > & NewLogLik )
{

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

	//cout << "A=" << A << endl;
	
      MUATTEMPT +=1;
      if( ranf()<A ){ //accept move 
	MUACCEPT +=1;
	for(int r=0; r<NREGION; r++){
	  for(int k=0; k<NSPECIES; k++){
	    Theta[r][k][l][j] = NewTheta[r][k];
	    ExpTheta[r][k][l][j] = NewExpTheta[r][k];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r][k];
	    LogLik[r][k][l] = NewLogLik[r][k];
	    //cout << "Accept!" << endl;
	  }
	}
	Mu[l][j] = NewMu;
      }
    }
  }
}

// update Lambda (the background species abundance)
void update_Lambda(double * Lambda, double Eta, double ** Psi, double ** ExpPsi, double * SumExpPsi, vector<int> & Region, vector<int> & Species, vector<double> & NewPsi,vector<double>  & NewExpPsi,vector <double > & NewSumExpPsi )
{

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
	if((Species[ind] == k))
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

void update_Alpha0(vector<double> & Alpha, double **** X)
{
  double sum = 0;
  double sumsq = 0;
  int total =0;

  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0;l<NLOCI; l++){
	for(int allele=0; allele<Nallele[l]; allele++){
	  sumsq += X[r][k][l][allele] * X[r][k][l][allele];
	  //sum += X[r][k][l][allele];
          total += 1;
	}
      }
    }
  }
  //cout << "Mean X = " << sum/total << endl;
  //cout << "MeanSqX = " << sumsq/total << endl;
  Alpha[0] = rgamma(0.5*total+NA[0],0.5*sumsq+LA[0]);
}


void update_Delta0(vector<double> & Delta, double ** Y)
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
void update_Beta(double & Beta, double ** Mu)
{
  double sumsq = 0;
   double sum = 0;
  int total =0;
  for(int l=0;l<NLOCI; l++){
    for(int allele=0; allele<Nallele[l]; allele++){
      sumsq += Mu[l][allele] * Mu[l][allele];
      //sum += Mu[l][allele];
      total += 1;
    }
  }
  //cout << "Mean Mu = " << sum/total << endl;
  //cout << "Sum Sq Mu = " << sumsq/total << endl;

  Beta = rgamma(0.5*total+NBETA,0.5*sumsq+LBETA);
}



// Eta is the prior precision for Lambda (ie Lambda is N(0,1/Eta)
// prior on Eta is Gamma(NETA,LETA)
void update_Eta(double & Eta, double * Lambda)
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


// update location of the sample whose location is being estimated
void update_Location(vector<double> & Alpha, double **** X, double ** Mu, double *** Nu, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount, double * L, vector<double> & Xcoord, vector<double> & Ycoord, vector<int> & Species, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<double> & BoundaryX, vector<double> & BoundaryY)
{
  double ** NewTheta = new double * [NLOCI];
  double ** NewExpTheta = new double * [NLOCI];
  double * NewLogLik = new double [NLOCI];
  double * NewSumExpTheta = new double [NLOCI];
  double * NewL = new double [NREGION * NREGION];

 
  for(int l=0; l<NLOCI; l++){
    NewTheta[l] = new double [MAXNALLELE];
    NewExpTheta[l] = new double [MAXNALLELE];	
    for(int j=0;j<MAXNALLELE;j++){	 
      NewTheta[l][j] = 0;
      NewExpTheta[l][j] = 1;
    }
    NewLogLik[l] = 0;
    NewSumExpTheta[l] = 0;	
  }


  for(int rep = 0; rep < 10; rep++){
    vector<double> NewXcoord(Xcoord);
    vector<double> NewYcoord(Ycoord);
    double h; // proposal sd
    double LogLikRatio = 0;
    
    if(ranf()<0.9){ // propose small move from current position
      h=0.02; //0.07;
      NewXcoord[NREGION-1] += rnorm(0,h);
      NewYcoord[NREGION-1] += rnorm(0,h);

      LogLikRatio = 0; // Hastings ratio
    } else { // propose jump to new randomly-chosen location
      h=0.04; //0.02; //sd of proposal
      //int newregion = 11;
      int newregion = (int) (ranf() * (NREGION - 1));
      //cout << "proposing move near region " << newregion << endl;
      NewXcoord[NREGION-1] = Xcoord[newregion] + rnorm(0,h);
      NewYcoord[NREGION-1] = Ycoord[newregion] + rnorm(0,h);
      double forwardsprob = 0; // prob density of proposed move
      double backwardsprob = 0; // prob density of backwards move
      for(int r = 0; r < (NREGION-1); r++){
	forwardsprob += dnorm((NewXcoord[NREGION-1]-Xcoord[r])/h) * dnorm((NewYcoord[NREGION-1]-Ycoord[r])/h);
	backwardsprob += dnorm((Xcoord[NREGION-1]-Xcoord[r])/h) * dnorm((Ycoord[NREGION-1]-Ycoord[r])/h);
      }
      //cout << "forwards prob = " << forwardsprob << endl;
      //cout << "backwards prob = " << backwardsprob << endl;
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
      
    bool isinrange;

    if(FORESTONLY){
      isinrange = InElephantRange(NewXcoord[NREGION-1],NewYcoord[NREGION-1]) && InForest(NewXcoord[NREGION-1],NewYcoord[NREGION-1]);
    } else if(SAVANNAHONLY){
      isinrange = InElephantRange(NewXcoord[NREGION-1],NewYcoord[NREGION-1]) && !InForest(NewXcoord[NREGION-1],NewYcoord[NREGION-1]);
    } else
      isinrange = (IsInsideBoundary(NewXcoord[NREGION-1],NewYcoord[NREGION-1], BoundaryX, BoundaryY) != 0);
  
    if(!isinrange)
      continue;
        
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

    
    // for(int l=0; l<NLOCI; l++){
//       NewLogLik[l] = 0;
//       NewLogLik[l] -= (SumCount[r][k][l]) * log(NewSumExpTheta[l]+Nallele[l]*DELTA);
//       LogLikRatio -= SumCount[r][k][l] * (log(NewSumExpTheta[l]+Nallele[l]*DELTA)-log(SumExpTheta[r][k][l]+Nallele[l]*DELTA));
//       for(int j=0; j<Nallele[l]; j++){
// 	NewLogLik[l] += Count[r][k][l][j] * log(NewExpTheta[l][j]+DELTA) ;
// 	LogLikRatio += Count[r][k][l][j] * (log(NewExpTheta[l][j]+DELTA) - log( ExpTheta[r][k][l][j] + DELTA));	  
//       }    
    //}
   

    double tempprob = 1;
    double newtempprob = 1;

    
    for(int l=0; l<NLOCI; l++){
      NewLogLik[l] = 0;
      LogLik[r][k][l] = 0;
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
	  
	  if((allele1 == allele2)){
	    NewLogLik[l] += log(NULLPROB*Newp1 + (1-NULLPROB) * Newp1 * Newp1);
	    LogLik[r][k][l] += log( NULLPROB*p1 + (1-NULLPROB) * p1 * p1);
	  }
	  else{
	    NewLogLik[l] += log((1-NULLPROB) *  Newp1 * Newp2);
	    LogLik[r][k][l] += log((1-NULLPROB) *  p1 * p2 );
	  }
	  
	}
      }
      LogLikRatio += NewLogLik[l] - LogLik[r][k][l];
      tempprob *= exp(LogLik[r][k][l]); 

      newtempprob *= exp(NewLogLik[l]);
    }

    //    cout << setiosflags(ios::scientific) << "Tempprobs:" << tempprob << " , " << newtempprob << endl;
    
    LOCATTEMPT +=1;
    double A = exp(LogLikRatio); // acceptance prob
    
    // test
    if( ranf()<A ){ //accept move
      LOCACCEPT +=1;
       Xcoord = NewXcoord;
       Ycoord = NewYcoord;
       for(int r0=0; r0<NREGION; r0++){
	 for(int s=0; s <NREGION; s++){
	   L[r0+ s*NREGION] = NewL[r0 + s*NREGION];
	 }
       }
       for(int l=0; l<NLOCI; l++){	
	 LogLik[r][k][l] = NewLogLik[l];
	 SumExpTheta[r][k][l] = NewSumExpTheta[l];
	 for(int j=0;j<Nallele[l];j++){  
	   Theta[r][k][l][j] = NewTheta[l][j];
	   ExpTheta[r][k][l][j] = exp(NewTheta[l][j]);
	 }
       }
    }
 
    //    compute the different region probs for this individual's genotype
    //cout << "Region probs:" << endl;
    //for(int r=0;r<NREGION;r++){
    //  for(int k=0; k<NSPECIES; k++){
//	double prob = 1;
//	for(int l=0;l<NLOCI; l++){
//	  for(int chr =0; chr <2; chr++){
//	    if(Genotype[SAMPLETOLOCATE][chr][l]>=0){
//	      prob *= ExpTheta[r][k][l][Genotype[SAMPLETOLOCATE][chr][l]]/SumExpTheta[r][k][l];
//	    }	
//	  }	
//	}
	//cout << setiosflags(ios::scientific) << r << ":" << prob << endl;
 //     }
//    }
  
  }


  for(int l=0; l<NLOCI; l++){
    delete [] NewTheta[l];
    delete [] NewExpTheta[l];
  }

  delete [] NewTheta;
  delete [] NewExpTheta;
  delete [] NewSumExpTheta;
  delete [] NewLogLik;
  delete [] NewL;
}


// update the parameters in the covariance matrix (would be neater
// not to declare and delete memory every time!)
void update_Alpha(vector<double> & Alpha, double **** X, double ** Mu, double *** Nu, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount, double * L, vector<double> & Xcoord, vector<double> & Ycoord)
{
  double **** NewTheta = new double *** [NREGION];
  double **** NewExpTheta = new double *** [NREGION];
  double *** NewLogLik = new double ** [NREGION];
  double *** NewSumExpTheta = new double ** [NREGION];
  
  update_Alpha0(Alpha, X);

  for(int r=0; r<NREGION; r++){
    NewTheta[r] = new double ** [NSPECIES];
    NewExpTheta[r] = new double ** [NSPECIES];    
    NewLogLik[r] = new double * [NSPECIES];
    NewSumExpTheta[r] = new double * [NSPECIES];
    for(int k=0; k<NSPECIES; k++){
      NewTheta[r][k] = new double * [NLOCI];
      NewExpTheta[r][k] = new double * [NLOCI]; 
      NewLogLik[r][k] = new double [NLOCI];
      NewSumExpTheta[r][k] = new double  [NLOCI];
      for(int l=0; l<NLOCI; l++){
	NewTheta[r][k][l] = new double [MAXNALLELE];
	NewExpTheta[r][k][l] = new double [MAXNALLELE];	
	for(int j=0;j<MAXNALLELE;j++){	 
	  NewTheta[r][k][l][j] = 0;
	  NewExpTheta[r][k][l][j] = 1;
	}
	NewLogLik[r][k][l] = 0;
	NewSumExpTheta[r][k][l] = 0;	
      }
    }
  }
    
  double * NewL = new double [NREGION * NREGION];
  
  
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
     
    calc_LogLik(NewLogLik,NewTheta,NewSumExpTheta,Count,SumCount);

    for(int r=0; r<NREGION; r++){ 
      for(int k=0; k<NSPECIES; k++){
	for(int l=0; l<NLOCI; l++){
	  //cout << NewLogLik[r][k][l] << "-" << LogLik[r][k][l] << endl;
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
      //cout << "Accept!" << endl; 
    } else {
      //cout << "Reject!" << endl;
    }
  }

  for(int r=0; r<NREGION; r++){
    for(int k=0; k<NSPECIES; k++){
      for(int l=0; l<NLOCI; l++){
	delete [] NewTheta[r][k][l];
	delete [] NewExpTheta[r][k][l];
      }
      delete [] NewTheta[r][k];
      delete [] NewExpTheta[r][k];
      delete [] NewLogLik[r][k];
      delete [] NewSumExpTheta[r][k];
    }
    delete [] NewTheta[r];
    delete [] NewExpTheta[r];
    delete [] NewLogLik[r];
    delete [] NewSumExpTheta[r];
  }
	
  delete [] NewTheta;
  delete [] NewExpTheta;
  delete [] NewSumExpTheta;
  delete [] NewLogLik;
  delete [] NewL;

}
// update the parameters in the covariance matrix M (would be neater
// not to declare and delete memory every time!)
void update_Delta(vector<double> & Delta, double ** Y, double * Lambda, double ** Psi, double ** ExpPsi, double * SumExpPsi, vector<int> & Species, vector<int> & Region, double * M, vector<double> & Xcoord, vector<double> & Ycoord)
{
  double ** NewPsi = new double * [NREGION];
  double ** NewExpPsi = new double * [NREGION];
  double * NewSumExpPsi = new double [NREGION];
  
  update_Delta0(Delta, Y);

  for(int r=0; r<NREGION; r++){
    NewPsi[r] = new double [NSPECIES];
    NewExpPsi[r] = new double [NSPECIES];    
    for(int k=0; k<NSPECIES; k++){
      NewPsi[r][k] = 0;
      NewExpPsi[r][k] = 1;
    }
    NewSumExpPsi[r] = 0;	
  }
    
  double * NewM = new double [NREGION * NREGION];

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

  for(int r=0; r<NREGION; r++){
    delete [] NewPsi[r];
    delete [] NewExpPsi[r];
  }
  	
  delete [] NewPsi;
  delete [] NewExpPsi;
  delete [] NewSumExpPsi;
  delete [] NewM;

}

// derivative of LogLik with respect to x(r,k,l,j) 
double divLogLikValue(int r, int k, int l, int j, double **** ExpTheta, double *** SumExpTheta, int **** Count, int *** SumCount, double * L)
{
  double sum = 0;
  for(int s=r; s<NREGION; s++){
    sum += (Count[s][k][l][j] - SumCount[s][k][l] * ExpTheta[s][k][l][j]/SumExpTheta[s][k][l]) * L[s + r*NREGION];
  }
  return sum;
}

double calcNewdivLogLik(int r, int k, int l, int j, double * ExpTheta, double * SumExpTheta, int **** Count, int *** SumCount, double * L)
{
  double sum = 0;
  for(int s=r; s<NREGION; s++){
    sum += (Count[s][k][l][j] - SumCount[s][k][l] * ExpTheta[s]/SumExpTheta[s]) * L[s + r*NREGION];
  }
  return sum;
}


// update the Xs for a particular allele and locus, in a particular species, across all regions at once
void update_XJoint(vector<double> & Alpha, double **** X, double ** Mu, double *** Nu,  double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount, double * L)
{

  double * NewTheta = new double [NREGION];
  double * NewExpTheta = new double [NREGION];
  double * NewSumExpTheta = new double [NREGION];
 double * NewLogLik = new double [NREGION];
  double * NewX = new double [NREGION];
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
	  //cout << X[r][k][l][j] << " " << divLogLik[r][k][l][j] << endl;
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
	//cout << "proposing move " << currenttheta << "->" << newtheta << endl;
	double A = exp(LogLikRatio); // acceptance prob

	
	//cout << "A=" << A << endl;
	XATTEMPT +=1;
	
	if( ranf()<A ){ //accept move 
	  XACCEPT +=1;
	  for(int r=0; r<NREGION; r++){
	    Theta[r][k][l][j] = NewTheta[r];
	    ExpTheta[r][k][l][j] = NewExpTheta[r];
	    SumExpTheta[r][k][l] = NewSumExpTheta[r];
	    X[r][k][l][j] = NewX[r];
	    LogLik[r][k][l] = NewLogLik[r];    
	    //cout << "Accept!" << endl;
	  }
	}
      }
    }
  }

  delete [] NewTheta;
  delete [] NewExpTheta;
  delete [] NewSumExpTheta;
  delete [] NewX;
  delete [] NewLogLik;

}

// update each X individually (better acceptance rate/ larger proposal variance,// but more likelihood evaluations!)
void update_XSingle(vector<double> & Alpha, double **** X, double ** Mu, double *** Nu, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount, double * L)
{

  double * NewTheta = new double [NREGION];
  double * NewExpTheta =  new double [NREGION];
  double * NewSumExpTheta = new double [NREGION];
  double * NewLogLik = new double [NREGION];
  double NewX;
  double NewdivLogLik;
  double h = XPROPOSALFACTOR * sqrt(1.0/Alpha[0]);


  for(int k=0; k<NSPECIES; k++){
    for(int l=0; l<NLOCI; l++){
      for(int j=0;j<Nallele[l];j++){	
	for(int r=0; r<NREGION; r++){
	  double LogLikRatio = 0;
	  double propmean = X[r][k][l][j]; // proposal mean
	  
	  if(USELANGEVIN == 1) 
	    propmean += (h*h/2)*(divLogLikValue(r,k,l,j,ExpTheta,SumExpTheta,Count,SumCount,L) - Alpha[0] * X[r][k][l][j]);
	  
	  //cout << X[r][k][l][j] << " " << divLogLik[r][k][l][j] << endl;
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
	    //cout << "HERE!" << r << k << l << j << endl;
	    NewdivLogLik = calcNewdivLogLik(r,k,l,j,NewExpTheta,NewSumExpTheta,Count,SumCount,L);
	    double backpropmean = NewX + (h*h/2) * (NewdivLogLik - Alpha[0] * NewX);
	    LogLikRatio -= 0.5 * (X[r][k][l][j] - backpropmean)*(X[r][k][l][j] - backpropmean)/(h*h); // top part of Hastings ratio
	  }
	  
	  //cout << "proposing move " << currenttheta << "->" << newtheta << endl;
	  double A = exp(LogLikRatio); // acceptance prob
	
	//cout << "A=" << A << endl;
	  XATTEMPT +=1;
	  
	  if( ranf()<A ){ //accept move 
	    XACCEPT +=1;
	    X[r][k][l][j] = NewX;
	    for(int s=r; s<NREGION; s++){
	      Theta[s][k][l][j] = NewTheta[s];
	      ExpTheta[s][k][l][j] = NewExpTheta[s];
	      SumExpTheta[s][k][l] = NewSumExpTheta[s];  
	      LogLik[s][k][l] = NewLogLik[s];
	      //cout << "Accept!" << endl;
	    }
	  }

	}
      }
    }
  }

  delete [] NewTheta;
  delete [] NewExpTheta;
  delete [] NewSumExpTheta;
  delete [] NewLogLik;

}

// update each Y individually (better acceptance rate/ larger proposal variance,// but more likelihood evaluations!)
void update_YSingle(vector<double> & Delta, double ** Y, double ** Psi, double ** ExpPsi, double * SumExpPsi, vector<int> & Region, vector<int> & Species, double * M)
{

  double * NewPsi = new double [NREGION];
  double * NewExpPsi =  new double [NREGION];
  double * NewSumExpPsi = new double [NREGION];
 
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
	  if((Species[ind] == k))
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

  delete [] NewPsi;
  delete [] NewExpPsi;
  delete [] NewSumExpPsi;


}


void DoAllUpdates(double **** X, double & Beta,  vector<double> & Gamma, vector<double> & Alpha,  double ** Mu, double *** Nu, double **** Theta, double **** ExpTheta, double *** LogLik, double *** SumExpTheta, int **** Count, int *** SumCount, double * L, vector<vector<vector<int> > > & Genotype, vector<int> & Region, vector<int> & Species, vector<vector<double> > & Pi, vector<double> & Xcoord, vector<double> & Ycoord,vector <vector<double> > & NewTheta,vector <vector<double> > & NewExpTheta,vector <vector<double> > & NewSumExpTheta,vector <vector<double> > & NewLogLik, double ** Y, double & Eta, vector<double> & Delta, double * Lambda, double ** Psi, double ** ExpPsi, double * SumExpPsi, double * M,  vector<double>  & NewPsi,vector<double>  & NewExpPsi,vector<double> & NewSumExpPsi,vector<double> & BoundaryX, vector<double> & BoundaryY )
{
  if(UPDATEBETA ==1)
    update_Beta(Beta,Mu);
  if(UPDATEX){
    if(UPDATEJOINT ==1)
      update_XJoint(Alpha,X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L);
    else
      update_XSingle(Alpha, X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L);
  }
  if(UPDATEALPHA ==1)
    update_Alpha(Alpha,X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Xcoord,Ycoord);

  if(UPDATEMU == 1)
    update_Mu(Mu,Beta,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount, NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik);

  if(UPDATENU == 1){
    update_Nu(Nu,Gamma,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount, NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik);
    for(int k=0; k<NSPECIES; k++)
      update_Beta(Gamma[k],Nu[k]);
  }

  if(LOCATE)
    update_Location(Alpha,X,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Xcoord,Ycoord,Species,Genotype,Region,BoundaryX, BoundaryY);

  if(NSPECIES>1){
    update_Species(Species, Pi, Region, Genotype, ExpTheta, SumExpTheta);
    count_up_alleles(Count,Region,Species,Genotype);
    calc_SumCount(Count, SumCount);

    //update_Pi(Pi, Species, Region);
    
    update_Lambda(Lambda,  Eta, Psi, ExpPsi, SumExpPsi, Region, Species, NewPsi, NewExpPsi,NewSumExpPsi );
    update_Eta(Eta,Lambda);
    update_YSingle(Delta,Y,Psi,ExpPsi,SumExpPsi,Region,Species,M);
    update_Delta(Delta,Y,Lambda,Psi,ExpPsi,SumExpPsi,Species,Region,M,Xcoord,Ycoord);
    compute_Pi(Pi, ExpPsi, SumExpPsi);

    calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
  }
 
}

void InitialiseTheta(double **** Theta, double **** X,double ** Mu, double *** Nu, double * L){
  
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

void Initialise(double **** X, double & Beta,  vector<double> & Gamma, vector<double> & Alpha,  double ** Mu, double *** Nu, double **** Theta, double * L, vector<double> & Xcoord, vector<double> & Ycoord){ 

  if(INCLUDENUGGET)
    Alpha[3] = 1.0;

  Alpha[0] = ALPHA0;
  Alpha[1] = ALPHA1;
  Alpha[2] = ALPHA2;
  Beta = BETA;

//   cerr << "Initialising: " << endl;
//   cerr << "Alpha[0] = " << ALPHA0 << endl;
//   cerr << "Alpha[1] = " << ALPHA1 << endl;
//   cerr << "Alpha[2] = " << ALPHA2 << endl;
//   cerr << "Beta     = " << BETA << endl;  

//   if(INITALPHA == 0){
//     Alpha[0] = 0.28; //10.0;
//     Alpha[1] = 1390; //10; //1.0/5000;
//     Alpha[2] = 1.24; //1.0;
//     Beta = 0.42;
//   }

//   if(INITALPHA == 1){ // point values for ak.highqual.neworder.inp
//     Alpha[0] = 0.28; //10.0;
//     Alpha[1] = 1390; //10; //1.0/5000;
//     Alpha[2] = 1.24; //1.0;Alpha[0] = 0.27; 
//     Beta = 0.42;
//   }

//   if(INITALPHA == 2){ // for ak.wcf.highqual.inp
//     Alpha[0] = 0.43; 
//     Alpha[1] = 5300; 
//     Alpha[2] = 0.32; 
//     Beta = 2.3;
//   }

//   if(INITALPHA == 3){ // for ak.nse.highqual.inp
//     Alpha[0] = 0.6; 
//     Alpha[1] = 5200; 
//     Alpha[2] = 0.84; 
//     Beta = 0.5;
//     Alpha[3] = 0;
//   }


//   if(INITALPHA == 4){ // for ak.wcf.highqual.inp, with no nugget.
//     Alpha[0] = 0.42; 
//     Alpha[1] = 5200; 
//     Alpha[2] = 0.32; 
//     Beta = 2.3;
//     Alpha[3] = 0;
//   }

//   if(INITALPHA == 5){ // for ak.nse.highqual.inp
//     Alpha[0] = 0.445; 
//     Alpha[1] = 6757; 
//     Alpha[2] = 0.855; 
//     Beta = 0.576;
//     Alpha[3] = 0;
//   }


// //   if(INITALPHA == 5){ //
// //     Alpha[0] = 1e-10; // 1e10 for testing, to have no variation about mean 
// //     Alpha[1] = 5200; 
// //     Alpha[2] = 0.84; 
// //     Beta = 1e100; // no variability in the mus
// //   }

//   if(INITALPHA == 6){ // for forest, after adding nugget
//     Alpha[0] = 1.184; 
//     Alpha[1] = 2269; 
//     Alpha[2] = 1.45; 
//     Alpha[3] = 0.89;
//     Beta = 1.18; 
//   }
 
//   if(INITALPHA == 7){ // for forest, after adding nugget
//     Alpha[0] = 0.01; 
//     Alpha[1] = 2269; 
//     Alpha[2] = 1.45; 
//     Alpha[3] = 0.89;
//     Beta = 1.18; 
//   }

//    if(INITALPHA == 8){ // for forest, after adding nugget
//     Alpha[0] = 100; 
//     Alpha[1] = 2269; 
//     Alpha[2] = 1.45; 
//     Alpha[3] = 0.89;
//     Beta = 1.18; 
//   }

  



  //for(int r =0; r<NREGION; r++)
  //  cout << "HERE: " << Xcoord[r] << "," << Ycoord[r] << endl;
  
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

  InitialiseTheta(Theta, X,Mu,Nu,L);

}


void output_empirical_freqs(vector<string> & RegionName,vector<int> & Perm,int ** Coding,int **** Count,int *** SumCount){
	
 	cout << "Empirical Frequencies:" << endl;
    for(int l=0;l<NLOCI; l++){
		cout << "Locus " << (l+1) << endl;
		cout << setw(5) << " " << " : ";
		OutputRegionNames(cout, RegionName, Perm);
		for(int allele=0; allele<Nallele[l]; allele++){
		    cout << setw(5) << Coding[l][allele] << " : ";
		    for(int r=0;r<NREGION;r++){
		        cout << setiosflags(ios::fixed) << setprecision(3) << setw(10) << (1.0*Count[Perm[r]][0][l][allele]) << "/" << SumCount[Perm[r]][0][l] << " ";
		    }
		    cout << endl;
		}
	}
}

void OutputMeanFreq(ofstream & freqfile, vector<string> & RegionName, vector<int> & Perm, int ** Coding, double **** MeanFreq )
{
   freqfile << "Posterior Mean Freqs:" << endl;
   for(int l=0;l<NLOCI; l++){
      freqfile << "Locus " << (l+1) << endl;
       freqfile << setw(5) << " " << " : ";
  	    OutputRegionNames(freqfile, RegionName, Perm);
  	    for(int allele=0; allele<Nallele[l]; allele++){
   	      freqfile << setw(5) << Coding[l][allele] << " : ";
          for(int r=0;r<NREGION-LOCATE;r++){
	     // for(int r=0;r<NREGION;r++){ 
   		    double f=MeanFreq[Perm[r]][0][l][allele];
  		    if(f<1e-4) f=0;
   		    freqfile << setiosflags(ios::fixed) << setprecision(3) << setw(10) << f << " ";
  	      }
   	      freqfile << endl;
        }
   }
}

int main ( int argc, char** argv)
{
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

    case 'c': // cross-validate
      CROSSVAL = 1;
      STARTCROSSVAL = atoi(&argv[1][2]);
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
	  cout << "Assigning!" << endl;
	  
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

    case 'o':
      ++argv; --argc; OFFSET = atoi(&argv[1][0]);
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

	case 's': // directory for output
	  filenames["assigndir"] = argv[1]+2;
	  break;
	  
    case 'S': // seed
      ++argv;
      --argc;
      SEED = atoi(&argv[1][0]);
      break;

      
    case 't': //tag to use for location output files
      ++argv;
      --argc;
      TAG = argv[1][0];
      cerr << "tag = " << TAG << endl;
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

    case 'w': // don't use spatial smooting
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
      
    case 'z' : // DC Special
      DCSPECIAL = 1;
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
  
  if(ASSIGNFILE==1)
	  if(HYBRIDCHECK==0)
		  LOCATE = 1;
  
  if(LOCATE)
    if(LASTSAMPLETOLOCATE < 0)
      LASTSAMPLETOLOCATE = FIRSTSAMPLETOLOCATE;

  if(argc<5){
    cerr << "Usage is ./SCAT genotypefile locationfile outputfile NLOCI Niter Nthin Nburn " << endl;
    exit(1);
  }

  srandom(SEED);

  filenames["input"]  = argv[1];
  filenames["regions"] = argv[2];
  filenames["output"] = argv[3];
  
  //NIND = atoi(argv[4]);
  NLOCI = atoi(argv[4]);
  //NREGION = atoi(argv[6]);

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

 
  if(NSPECIES == 1)
    UPDATENU = 0;

  ifstream input (filenames["input"].c_str());
  ifstream regionfile (filenames["regions"].c_str());
  string outputfilename = filenames["output"] + "_probs";
  ofstream output (outputfilename.c_str());

  ifstream assignfile (filenames["assignfile"].c_str());
  string freqfilename = filenames["output"]+"_freqs";
  ofstream freqfile (freqfilename.c_str());

  string acceptfilename = filenames["output"] + "_accept";
  ofstream acceptfile (acceptfilename.c_str());

  string pifilename = filenames["output"] + "_pi";
  ofstream pifile;
  if(NSPECIES>1)
    pifile.open(pifilename.c_str());

  string corrfilename = filenames["output"] + "_corr";
  ofstream corrfile (corrfilename.c_str());

  string paramfilename = filenames["output"] + "_params";
  ofstream paramfile (paramfilename.c_str());
  ofstream Xfile (filenames["Xfile"].c_str());
  
  //assure ( input, "inputfile" );

  vector<int> Region;

  // input genotype data from inputfile

  vector<int> RegionsPresent(0);
  
  vector<string> Id(0);

    // declare memory for Genotype and OriginalGenotype
  vector<vector<vector<int> > > Genotype;
  vector<vector<vector<int> > > OriginalGenotype;
   
                                 // Coding[i] is the actual allele that i codes for
  int ** Coding = new int * [NLOCI];
  for(int l=0; l<NLOCI; l++){
     Coding[l] = new int [MAXNALLELE];
  }
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
  
  output_genotypes(OriginalGenotype,Id);
  recode_genotypes(OriginalGenotype,Genotype,Coding,Nallele);
  //
  //           //  output_genotypes(Genotype,Id);
  //
  cout << "Number of Alleles at each locus:" << endl;
  int j=1;
  cout << "Locus : #alleles" << endl;
  for(vector<int>::iterator i = Nallele.begin(); i!=Nallele.end(); i++)
      cout << j++ << " : " << *i << endl;
  
  
  NREGION = RegionsPresent.size()+LOCATE;
  
  cout << "Number of Regions: " << NREGION << endl;
  
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

  vector<double> Xcoord(NREGION,0);
  vector<double> Ycoord(NREGION,0);
  vector<double> BoundaryX;
  vector<double> BoundaryY;
  if(READBOUNDARY){
    cerr << "Reading in Boundary data from boundary file " << filenames["boundaryfile"] << endl;
    ifstream bfile (filenames["boundaryfile"].c_str());
    ReadInBoundary(bfile,BoundaryX,BoundaryY);
  }

  vector<int> SubRegion(NREGION,0);
  

  //declare memory for theta, and SumExpTheta
  //theta[r][k][l][j] is log relative freq of allele j at locus l
  // in region r, species k, relative to allele j=0.
  // SumExpTheta[r][k][l] is the sum of exp(theta) over j

  // declare memory
  // Count[r][k][l][j] is number of allele j in region r, species k, locus l
  // SumCount[r][k][l] is the sum of Count over j
  int **** Count = new int *** [NREGION];
  int *** SumCount = new int ** [NREGION];

  // Theta[r][k][l][j] is log relative freq of allele j in region r,
  // species k, locus l
  // ExpTheta is exponential of theta (stored to save computational time) 
  // SumExpTheta is sum of ExpTheta over j
  // So Freq[r][k][l][j] = ExpTheta[r][k][l][j]/SumExpTheta[r][k][l]
  // mu, X and L are such that Theta = mu + LX
  // where X are iid N(0,Alpha[0])
  // L is the (lower-triangular) Cholesky decomposition of the covariance of the thetas
  // mu is a priori N(0,1/beta), where beta is a priori Gamma(NBETA,LBETA)

  double **** Theta = new double *** [NREGION];
  double **** ExpTheta = new double *** [NREGION];
  double *** SumExpTheta = new double ** [NREGION];
  double **** X = new double *** [NREGION]; // Theta = mu + Nu + LX
  double ** Mu = new double * [NLOCI]; 
  double *** Nu = new double ** [NSPECIES];
  
  // these vectors used to hold proposed new values
  vector< vector<double> > NewTheta(NREGION, vector<double>(NSPECIES));
  vector< vector<double> > NewExpTheta(NREGION, vector<double>(NSPECIES));
  vector< vector<double> > NewSumExpTheta(NREGION, vector<double>(NSPECIES));
  vector< vector<double> > NewLogLik(NREGION, vector<double>(NSPECIES));

  // Psi[r][k] is log relative abundance of species k in region r
  // ExpPsi and SumExpPsi as in Theta
  // So Pi[r][k] = ExpPsi[r][k]/SumExpPsi[r]
  // Lambda, Y and M are such that Psi = Lambda + MY
  // with Y being iid N(0,Delta[0])
  // and M being lower-triangular Cholesky decomp of covariance of Psi
  // which itself is determined by Delta
  // Lambda[k] are apriori iid N(1/eta), eta is a priori Gamma(NETA,LETA)

  double ** Psi = new double * [NREGION];
  double ** ExpPsi = new double * [NREGION];
  double * SumExpPsi = new double [NREGION];
  double ** Y = new double * [NREGION]; // Theta = mu + Nu + LX
  double * Lambda = new double [NSPECIES];

  vector<double> NewPsi(NREGION);
  vector<double> NewExpPsi(NREGION);
  vector<double> NewSumExpPsi(NREGION);
  
  // these hold posterior means of Frequencies, X and X^2 in each region
  double **** MeanFreq = new double *** [NREGION];
  double **** MeanX = new double *** [NREGION]; // mean value of X
  double **** MeanX2 = new double *** [NREGION]; // mean squared value of X

  double **** MeanCov = new double *** [NREGION];
  double **** MeanFittedCov = new double *** [NREGION];
  double **** MeanCor = new double *** [NREGION];
  double **** MeanFittedCor = new double *** [NREGION];

  vector<vector<double> > MeanPi(NREGION, vector<double>(NSPECIES,0.0)); // mean value of Pi

  // LogLik[r][k][l] holds the log-likelihood for individuals in region r, 
  // species k, at locus l 
  double *** LogLik = new double ** [NREGION];
 
  for(int l=0; l<NLOCI; l++){
    Mu[l] = new double [MAXNALLELE];
  }
  for(int k=0; k<NSPECIES; k++){
    Nu[k] = new double * [NLOCI];
    for(int l=0; l<NLOCI; l++){
      Nu[k][l] = new double [MAXNALLELE];
    }
  }

  for(int r=0; r<NREGION; r++){
    Psi[r] =  new double [NSPECIES];
    ExpPsi[r] =  new double [NSPECIES];
    Y[r] = new double [NSPECIES];
    for(int k=0; k < NSPECIES; k++){
      Psi[r][k] =0;
      ExpPsi[r][k] = 1;
      Y[r][k] = 0;
    }
    SumExpPsi[r] = NREGION;
  }

  
  for(int r=0; r<NREGION; r++){
    Count[r] = new int ** [NSPECIES];
    Theta[r] = new double ** [NSPECIES];
    ExpTheta[r] =  new double ** [NSPECIES];
    X[r] = new double ** [NSPECIES];   
    LogLik[r] = new double * [NSPECIES];
    SumCount[r] = new int * [NSPECIES];
    SumExpTheta[r] = new double * [NSPECIES];
    MeanFreq[r] = new double ** [NSPECIES];
    MeanX[r] = new double ** [NSPECIES];
    MeanX2[r] = new double ** [NSPECIES];
    MeanCov[r] = new double ** [NSPECIES];
    MeanFittedCov[r] = new double ** [NSPECIES];
    MeanCor[r] = new double ** [NSPECIES];
    MeanFittedCor[r] = new double ** [NSPECIES];
    for(int k=0; k<NSPECIES; k++){
      Count[r][k] = new int * [NLOCI];
      Theta[r][k] = new double * [NLOCI];
      ExpTheta[r][k] = new double * [NLOCI];
      X[r][k] = new double * [NLOCI];
      LogLik[r][k] = new double [NLOCI];
      SumCount[r][k] = new int  [NLOCI];
      SumExpTheta[r][k] = new double  [NLOCI];
      MeanFreq[r][k] = new double * [NLOCI];
      MeanX[r][k] = new double * [NLOCI];
      MeanX2[r][k] = new double * [NLOCI];
      MeanCov[r][k] = new double * [NREGION];
      MeanFittedCov[r][k] = new double * [NREGION]; 
      MeanCor[r][k] = new double * [NREGION];
      MeanFittedCor[r][k] = new double * [NREGION];
      for(int r1=0; r1<NREGION; r1++){
	MeanCov[r][k][r1] = new double [NSPECIES];
	MeanFittedCov[r][k][r1] = new double [NSPECIES];
	MeanCor[r][k][r1] = new double [NSPECIES];
	MeanFittedCor[r][k][r1] = new double [NSPECIES];
	for (int k1 = 0; k1<NSPECIES; k1++){
	  MeanCov[r][k][r1][k1] = 0;
	  MeanFittedCov[r][k][r1][k1] = 0;
	  MeanCor[r][k][r1][k1] = 0;
	  MeanFittedCor[r][k][r1][k1] = 0;
	}
      }
      for(int l=0; l<NLOCI; l++){
	Count[r][k][l] = new int [MAXNALLELE];
	Theta[r][k][l] = new double [MAXNALLELE];
	ExpTheta[r][k][l] = new double [MAXNALLELE];
	X[r][k][l] = new double [MAXNALLELE];
	MeanFreq[r][k][l] = new double [MAXNALLELE];
	MeanX[r][k][l] = new double [MAXNALLELE];
	MeanX2[r][k][l] = new double [MAXNALLELE];
	for(int j=0;j<MAXNALLELE;j++){	 
	  Count[r][k][l][j] = 0;
	  Theta[r][k][l][j] = Mu[l][j] + Nu[k][l][j];
	  ExpTheta[r][k][l][j] = exp(Theta[r][k][l][j]);
	  X[r][k][l][j] = 0;
	  MeanFreq[r][k][l][j] = 0;
	  MeanX[r][k][l][j] = 0;
	  MeanX2[r][k][l][j] = 0;
	}
	LogLik[r][k][l] = 0;
	SumCount[r][k][l] = 0;
	SumExpTheta[r][k][l] = 0;	
      }
    }
  }

  
 
  vector<double> Alpha(ALPHALENGTH,1); // parameters in covariance structure
  double Beta; // 1/variance of mus

  vector<double> Delta(DELTALENGTH,1);
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

 
  vector<double> Gamma (NSPECIES,0.3); // Gamma[k] is 1/variance of Nu[k]

  // input data from region file (Xcoords and Ycoords; regions and subregions)
  input_positions_data(regionfile,Xcoord,Ycoord,RegionName,SubRegion,Region,Perm,RegionsPresent);
  output_positions_data(RegionName, Region, Xcoord, Ycoord, Id);

  //declare memory for L, the matrix in the Cholesky decomp of Sigma
  double * L = new double [NREGION * NREGION];
  double * M = new double [NREGION * NREGION];

  Initialise(X, Beta, Gamma, Alpha, Mu, Nu, Theta, L, Xcoord, Ycoord);

  // MeanProb[ind][r][k] is Pr of ind's genotype in region r, species k
  // LocusMeanProb [ind][r][k][l] is Pr of ind's genotype at locus l in region r, species k
  double *** MeanProb = new double ** [NIND];
  double **** LocusMeanProb = new double *** [NIND];
  double *** SubRegionProb = new double ** [NIND];

  for(int ind = 0; ind< NIND; ind ++){
    MeanProb[ind] = new double * [NREGION];
    LocusMeanProb[ind] = new double ** [NREGION];
    SubRegionProb[ind] = new double * [NSUBREGION];
    for(int r=0;  r<NREGION; r++){
      MeanProb[ind][r] = new double [NSPECIES];
      LocusMeanProb[ind][r] = new double * [NSPECIES];
      for(int k =0; k<NSPECIES; k++){
	MeanProb[ind][r][k] = 0;
	LocusMeanProb[ind][r][k] = new double [NLOCI];
	for(int l=0; l<NLOCI; l++){
	  LocusMeanProb[ind][r][k][l] = 0;
	}
      }
    }
    for(int r=0;  r<NSUBREGION; r++){
      SubRegionProb[ind][r] = new double [NSPECIES];
      for(int k =0; k<NSPECIES; k++){
	SubRegionProb[ind][r][k] = 0;
      }
    }
  }

  // count up the number in each region
  count_up_alleles(Count,Region,Species,Genotype);
  if(VERBOSE)
    output_counts(Count,Perm);
  calc_SumCount(Count, SumCount);

  if(VERBOSE)
    output_empirical_freqs(RegionName,Perm,Coding,Count,SumCount);


  calc_ExpTheta_and_SumExpTheta(Theta,ExpTheta,SumExpTheta);
  calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
 
  double totaliter = 0; // number of iterations used in computing mean frequencies
  
 
  // start by doing burn-in (whether or not cross-validating)
  int templocate = LOCATE;
  LOCATE = 0; // don't locate during burnini


  if(!LOCATEWHOLEREGION){
    cerr << "Performing Burn-in iterations" << endl;
    for(int iter=0; iter< Nburn; iter++){
      for(int nthin=0; nthin < Nthin; nthin++){
	DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik, Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, NewPsi, NewExpPsi, NewSumExpPsi, BoundaryX, BoundaryY);
      }
      
      OutputAcceptRates(acceptfile);

      calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
      
      //OutputEstimatedFreqs( ExpTheta, SumExpTheta, output, Coding);
      //OutputTheta( Theta, output, Coding, Perm);
      if(NSPECIES > 1)
	OutputPi( Pi, RegionName, pifile, Perm );
       
      cerr << "Iteration:" << (iter+1) << "\033[A" << endl;

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
    cerr << endl;
  }

  
  LOCATE = templocate;
   
  // different ways of running program: Locate (Continuous and smooth
  // assignment test, cross-validated); HYBRIDCHECK (undocumented);
  // Non-cross val, assignment test without cross-validation
  totaliter = 0;
	
  if(LOCATE){
    for(SAMPLETOLOCATE = FIRSTSAMPLETOLOCATE; SAMPLETOLOCATE <= LASTSAMPLETOLOCATE; SAMPLETOLOCATE++){
      char firstdigit = (char) ('0' + (SAMPLETOLOCATE+OFFSET-FIRSTSAMPLETOLOCATE) / 100);
      char seconddigit = (char) ('0' +((SAMPLETOLOCATE+OFFSET-FIRSTSAMPLETOLOCATE) % 100) / 10);
      char thirddigit = (char)  ('0' +((SAMPLETOLOCATE+OFFSET-FIRSTSAMPLETOLOCATE) % 10));    
      
      string LOCATEFILE(filenames["assigndir"]);
      LOCATEFILE.push_back(firstdigit);
      LOCATEFILE.push_back(seconddigit);
      LOCATEFILE.push_back(thirddigit);
      LOCATEFILE.push_back(TAG);
      LOCACCEPT = 0;
      LOCATTEMPT = 0;

      ofstream locatefile (LOCATEFILE.c_str());
      
      cerr << "Individual:" << (SAMPLETOLOCATE+1) << endl;      
      //paramfile << "Individual:" << SAMPLETOLOCATE << endl;

      // reset theta, counts, etc, ignoring individual ind
      
      TRUEREGION = Region[SAMPLETOLOCATE];      
      Region[SAMPLETOLOCATE] = NREGION-1;
      if(LOCATEWHOLEREGION){ // locate all the inds from that region
	     for(int i = 0; i < NIND; i++){
	        if(Region[i] == TRUEREGION)
	           Region[i] = NREGION - 1;    
	     }
      }
      
	  cout << "Initialising XY position... ";
      cerr << "Initialising XY position... ";
      InitialiseXY(BoundaryX,BoundaryY,Xcoord,Ycoord);
      cerr << "Done" << endl;

      calc_L(L,Alpha,Xcoord,Ycoord);
      count_up_alleles(Count,Region,Species,Genotype);
		
      if(DCSPECIAL){ // case where we remove all but 3 of the DC samples
      // DC samples are 43-96
         for(int ind = 46; ind < 97; ind++) 
	        if(ind != SAMPLETOLOCATE)
	           SubtractFromCount(ind,Count,SumCount,Region,Species,Genotype);
      }
    
      if(REMOVEREGION){ // remove all individuals from the true region
	     for(int i = 0; i < NIND; i++){
	        if(Region[i] == TRUEREGION)
	            SubtractFromCount(i, Count, SumCount, Region, Species, Genotype);    
	     }
      }
  	  
      calc_SumCount(Count, SumCount);
      
	  if(VERBOSE)
	      output_empirical_freqs(RegionName,Perm,Coding,Count,SumCount);

		 
	  //output_counts(Count, Perm);
       //cout << "Simple Probabilities: " << endl;
       //for(int r=0; r<NREGION; r++)
 	//cout << r << ":" << simple_Prob(Count,SumCount,Genotype,Coding,SAMPLETOLOCATE,r) << endl;
    
      for(int r=0; r<NREGION; r++){
	for(int k=0; k<NSPECIES; k++){
	  for(int l=0; l<NLOCI; l++){
	    for(int j=0;j<Nallele[l];j++){  
	      X[r][k][l][j] = 0;
	      //Mu[l][j] = rnorm(0,1);
	    }
	  }
	}
      }
      InitialiseTheta(Theta,X,Mu,Nu,L);
      calc_ExpTheta_and_SumExpTheta(Theta,ExpTheta,SumExpTheta);
      calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
       
      for(int iter=0; iter< Nburn; iter++){
	cout << "Iteration:" << (iter+1) << "\033[A" << endl;
	for(int nthin=0; nthin < Nthin; nthin++){
	  DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, NewPsi, NewExpPsi, NewSumExpPsi, BoundaryX, BoundaryY);
	}
	// prevent rounding errors
	calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
 	//CheckSumExpTheta(ExpTheta,SumExpTheta);
	double sumloglik = SumLogLik(LogLik);
	OutputLatLongs(locatefile,Xcoord[NREGION-1],Ycoord[NREGION-1],sumloglik);
      	//locatefile << Xcoord[NREGION-1] << " " << Ycoord[NREGION-1] <<  " " << sumloglik << endl; 

      }
      
      for(int iter=0; iter< Niter; iter++){
	cout << "Iteration:" << (iter+1) << "\033[A" << endl;
	for(int nthin=0; nthin < Nthin; nthin++){
	  DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, NewPsi, NewExpPsi, NewSumExpPsi, BoundaryX, BoundaryY);
	}
	calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);
	
	double sumloglik = SumLogLik(LogLik);
	OutputLatLongs(locatefile,Xcoord[NREGION-1],Ycoord[NREGION-1],sumloglik);

	
//    cout << 180*Ycoord[NREGION-1]/PI << "," << 180*Xcoord[NREGION-1]/PI << "," << CurrentLogLik(LogLik,NREGION-1) << ":";
//    cout << 180*Ycoord[NREGION-2]/PI << "," << 180*Xcoord[NREGION-2]/PI << "," << endl;
//	cout << endl;
//    for(int r=0;r<NREGION; r++){
//		cout << Covariance(r,0,NREGION-1,0,Theta,Mu,Nu) << ",";
//	}
//	cout << endl;
//	for(int r=0;r<NREGION; r++){
//		cout << FittedCovariance(Alpha,Distance(r,NREGION-1,Xcoord,Ycoord)) << ",";
//	}
//	cout << endl;
//	
//     for(int r=0;r<NREGION; r++){
//		         cout << Distance(r,NREGION-1,Xcoord,Ycoord) << ",";
//	 }
//	 cout << endl;

		 
	

	UpdateMeans(ExpTheta, X, Pi, SumExpTheta, MeanFreq, MeanX, MeanX2, MeanPi, MeanCov, MeanFittedCov, MeanCor, MeanFittedCor, Alpha, Xcoord, Ycoord,Theta,Mu,Nu);	
	UpdateLocusMeanProb(SAMPLETOLOCATE, ExpTheta, SumExpTheta, LocusMeanProb, Genotype);
	
	totaliter++;
      }	
      
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
      cerr << endl;

    }
  } else if(HYBRIDCHECK) {
    
    double loglik = 0;
    
    cout << "Hybridcheck =" << HYBRIDCHECK << endl;
    // 1 or 2 for region-based (1 does hybrid, 2 does pure only); 11 or 12 for subregion-based

    if(HYBRIDCHECK <10){ // set regions to be subregion
      for(int ind = 0 ; ind<NIND; ind++){
	if(Region[ind]>=0)
	  Region [ind] = SubRegion[Region[ind]];
	cout << ind << "," << Region[ind] << endl;
      }
      NREGION = NSUBREGION;
    }
    count_up_alleles(Count,Region,Species,Genotype);

    if(DCSPECIAL){ // case where we remove all but 3 of the DC samples
      // DC samples are 43-96
      for(int ind = 46; ind < 97; ind++) 
	SubtractFromCount(ind,Count,SumCount,Region,Species,Genotype);
    }
    
    calc_SumCount(Count, SumCount); 
    output_counts(Count, Perm);
    cout << "Hybrid probs" << endl;
    cout.precision(4); 

    for(int ind = 0; ind < NIND; ind++){ 
      
      if(!DCSPECIAL)
	SubtractFromCount(ind, Count, SumCount, Region, Species, Genotype); 
      calc_SumCount(Count, SumCount);

      for(int r = 0; r < NREGION; r++){
	for(int s = 0; s< NREGION; s++){
	  if(HYBRIDCHECK==1 || HYBRIDCHECK == 11 || r==s)
	    cout << log_hybrid_Prob(Count,SumCount,Genotype,Coding,ind,r,s) << " ";
	  if((r == Region[ind]) && (r==s))
	    loglik += log_hybrid_Prob(Count,SumCount,Genotype,Coding,ind,r,s);
	}
      }
      cout << endl;

      if(!DCSPECIAL)
	AddToCount(ind, Count, SumCount, Region, Species, Genotype);   
      calc_SumCount(Count, SumCount);
    }
    
    cout.precision(6);
    cout << "Loglik = " << loglik << endl;

  } else { // if not cross-validating

    XACCEPT = 0;
    XATTEMPT = 0;
    cerr << "Performing Main Iterations " << endl;

    for(int iter=0; iter< Niter; iter++){     
      for(int nthin=0; nthin < Nthin; nthin++){
	DoAllUpdates(X,Beta,Gamma,Alpha,Mu,Nu,Theta,ExpTheta,LogLik,SumExpTheta,Count,SumCount,L,Genotype,Region,Species,Pi,Xcoord,Ycoord,NewTheta,NewExpTheta,NewSumExpTheta,NewLogLik,Y,Eta, Delta,Lambda,  Psi,ExpPsi, SumExpPsi, M, NewPsi, NewExpPsi, NewSumExpPsi, BoundaryX, BoundaryY);
      }

      //      cout << "LogLik Before =" << setprecision(10) << SumLogLik(LogLik) << endl;

      calc_LogLik(LogLik,Theta,SumExpTheta,Count,SumCount);     
      
      //cout << "LogLik After =" << setprecision(10) << SumLogLik(LogLik) << endl;
      //cout << "Acceptance Rate = " << (1.0*XACCEPT)/XATTEMPT << endl;
      //cout << "LogLik without Pseudocounts = " << calc_LogLikWithoutPseudoCounts(Theta,SumExpTheta,Count,SumCount);  

      if(NSPECIES > 1)
	OutputPi( Pi, RegionName, pifile, Perm );
      
      UpdateMeans(ExpTheta, X, Pi, SumExpTheta, MeanFreq, MeanX, MeanX2, MeanPi, MeanCov, MeanFittedCov, MeanCor, MeanFittedCor, Alpha, Xcoord, Ycoord, Theta, Mu, Nu);
      totaliter++;
      UpdateLocusMeanProb(ExpTheta, SumExpTheta, LocusMeanProb, Genotype);
      if(USESUBREGION==1){
	UpdateSubRegionProb(ExpTheta, SumExpTheta, SubRegionProb, Genotype, SubRegion);
      }
      //OutputEstimatedFreqs( ExpTheta, SumExpTheta, output, Coding, Perm);   
      //OutputTheta( MeanX, output, Coding, Perm);
      cerr << "Iteration:" << (iter+1) << "\033[A" << endl;

      OutputParameters(paramfile, Alpha, Beta, Gamma, Delta, Eta, Lambda, LogLik);


      if(OUTPUTX==1){
	for(int r=0; r < NREGION; r++){
	  for(int a =0; a<Nallele[0]; a++){
	    Xfile << X[r][0][0][a] << " ";
	  }
	}
	Xfile << endl;
      }
      
      // corrfile << "locus region1 region2 subregion1 subregion2 distance corr covar fittedcorr fittedcovar" << endl;
//       for(int l=0; l< NLOCI; l++){
// 	for(int r0=0; r0<NREGION; r0++){
// 	  for(int k0=0; k0<NSPECIES; k0++){
// 	    for(int r1=0; r1<=r0; r1++){
// 	      for(int k1=0; k1<=k0; k1++){
// 		corrfile << l << " " << r0 << " " << r1 << " " << SubRegion[r0] << " " << SubRegion[r1] << " " << Distance(r0,r1,Xcoord,Ycoord) << " " << Correlation(r0,k0,r1,k1,Theta,Mu, Nu) << " " << Covariance(r0,k0,r1,k1,Theta,Mu,Nu) << " " << FittedCovariance(Alpha,Distance(r0,r1,Xcoord,Ycoord))/FittedCovariance(Alpha,0) << " " << FittedCovariance(Alpha,Distance(r0,r1,Xcoord,Ycoord)) << endl;
// 	      }
// 	    }
// 	  }
// 	}
//       }
       
    }
    cerr << endl;
  }

 
  //NormaliseMeanProb(MeanProb); //CorrectMeanProb(MeanProb, totaliter); // divide MeanProb elements by totaliter
  
  NormaliseMeanFreq(MeanFreq,totaliter);
  ComputeMeanProb(MeanProb,LocusMeanProb,totaliter); 
  NormaliseMeanPi(MeanPi);
  NormaliseMeanCov(MeanCov, MeanFittedCov, totaliter);
  NormaliseMeanCov(MeanCor, MeanFittedCor, totaliter);

  cerr << "Creating output files" << endl;

  // output estimated position of each individual
  //OutputMeanPos(MeanProb,Xcoord,Ycoord,Region,RegionName,output,NMissing);

   // output probs of each individual being in each region
  OutputLogMeanProb(MeanProb,output,NMissing,Region,Id, RegionName, Perm, FIRSTSAMPLETOLOCATE, LASTSAMPLETOLOCATE, RegionsPresent);
  
  // output probs of each individual being in each region, according to mean freqs
  OutputLogMeanProb2(MeanFreq,Genotype,output,NMissing,Region,Perm);

  if(USESUBREGION==1){
    NormaliseSubRegionProb(SubRegionProb);
  }

  //
  // output misclassification matrices
  //

 //  double ** RegionMisclass = new double * [NREGION];
//   for(int r=0; r<NREGION; r++){
//     RegionMisclass[r] = new double [NREGION];
//     for(int s=0; s<NREGION; s++){
//       RegionMisclass[r][s] = 0;
//     }
//   }
  
//   double ** SubRegionMisclass = new double * [NSUBREGION];
//   for(int r=0; r<NSUBREGION; r++){
//     SubRegionMisclass[r] = new double [NSUBREGION];
//     for(int s=0;s<NSUBREGION; s++){
//       SubRegionMisclass[r][s] = 0;
//     }
//   }
  
  
//   for(int maxmissing =0; maxmissing < NLOCI; maxmissing++){

//     for(int r=0; r<NREGION; r++){
//       for(int s=0; s<NREGION; s++){
// 	RegionMisclass[r][s] = 0;
//       }
//     }
//     for(int r=0; r<NSUBREGION; r++){
//       for(int s=0; s<NSUBREGION; s++){
// 	SubRegionMisclass[r][s] = 0;
//       }
//     }

    
//     for(int i=0; i<NIND; i++){
//       if(NMissing[i]<=maxmissing){
// 	int bestclass = 0;
// 	for(int r=0; r<NREGION; r++){
// 	  if(MeanProb[i][r][0]>MeanProb[i][bestclass][0]){
// 	    bestclass = r;
// 	  }
// 	}
// 	RegionMisclass[Region[i]][bestclass]++;
// 	int bestsubregion = 0;
// 	for(int r=0; r<NSUBREGION; r++){
// 	  if(SubRegionProb[i][r][0] > SubRegionProb[i][bestsubregion][0]){
// 	    bestsubregion = r;
// 	  }
// 	}
      
// 	SubRegionMisclass[SubRegion[Region[i]]][bestsubregion]++;
//       }
//     }

//     for(int r=0; r< NREGION; r++){
//       for(int s=0; s<NREGION; s++){
// 	output << RegionMisclass[r][s] << " ";
//       }
//       output << endl;
//     }
    
//     for(int r=0; r< NSUBREGION; r++){
//       for(int s=0; s<NSUBREGION; s++){
// 	output << SubRegionMisclass[r][s] << " ";
//       }
//       output << endl;
//     }
//   }

  //  int ind =2;
//    double LR = 0;
//    for(int l = 0; l < NLOCI; l++){
//      cout << "Locus " << l << endl;
//      for(int r=0;r<NREGION;r++){
//        double prob = 1;	
//        for(int k=0; k<NSPECIES; k++){
//  	for(int chr =0; chr <2; chr++){
//  	  if(Genotype[ind][chr][l]>=0){
//  	    prob *= exp(Theta[r][k][l][Genotype[ind][chr][l]])/SumExpTheta[r][k][l];
//  	  }	
//  	}	
//        }
//        output << r << ":" << prob << endl;
//        if(r==1)
//  	LR += log(prob);
//        if(r==10)
//  	LR -= log(prob);
//      }
//      output << "LR=" << LR << endl;
//    }
  


//   output << "Variance of Xs" << endl;
//   for(int l=0;l<NLOCI; l++){
//     output << "Locus " << l << endl;
//     for(int allele=0; allele<Nallele[l]; allele++){
//       output << setw(5) << Coding[l][allele] << " : ";
//       for(int r=0;r<NREGION;r++){     
// 	output << setprecision(3) << setw(5) << 
// 	  MeanX2[r][0][l][allele]/totaliter - (MeanX[r][0][l][allele]/totaliter)*(MeanX[r][0][l][allele]/totaliter) << " ";
//       }
//       output << endl;
//     }
//   }

//   output << "MeanX values" << endl;
//   for(int l=0;l<NLOCI; l++){
//     output << "Locus " << l << endl;
//     for(int allele=0; allele<Nallele[l]; allele++){
//       output << setw(5) << Coding[l][allele] << " : ";
//       for(int r=0;r<NREGION;r++){     
// 	output  << setprecision(3) << setw(5) << MeanX[r][0][l][allele]/totaliter <<  " ";
//       }
//       output << endl;
//     }
//   }

//   output << "Mu Values" << endl;
//   for(int l=0;l<NLOCI; l++){
//     output << "Locus " << l << endl;
//     for(int allele=0; allele<Nallele[l]; allele++){
//       output << setw(5) << Coding[l][allele] << " : ";     
//       output << setprecision(3) << setw(5) << Mu[l][allele] <<  " ";
//       output << endl;
//     }
//   }
  

  //
  // output correlations and covariances (fitted and empirical) to corrfile
  //

  //corrfile << "locus location1 location2 distance corr covar fittedcorr fittedcovar" << endl;
  

  corrfile << "location1 location2 distance cov fittedcov corr fittedcorr" << endl;

   //for(int l=0; l< NLOCI; l++){
     for(int r0=0; r0<NREGION; r0++){
      for(int k0=0; k0<NSPECIES; k0++){
	for(int r1=0; r1<=r0; r1++){
	  for(int k1=0; k1<=k0; k1++){
	    //corrfile << l << " " << r0 << " " << r1 << " " << Distance(r0,r1,Xcoord,Ycoord) << " " << Correlation(r0,k0,r1,k1,Theta,Mu, Nu) << " " << Covariance(r0,k0,r1,k1,Theta,Mu,Nu) << " " << FittedCovariance(Alpha,Distance(r0,r1,Xcoord,Ycoord))/FittedCovariance(Alpha,0) << " " << FittedCovariance(Alpha,Distance(r0,r1,Xcoord,Ycoord)) << endl;
	    corrfile << RegionName[r0] << " " << RegionName[r1] << " " << Distance(r0,r1,Xcoord,Ycoord) << " " << MeanCov[r0][k0][r1][k1] << " " << MeanFittedCov[r0][k0][r1][k1]  << " " << MeanCor[r0][k0][r1][k1] << " " << MeanFittedCor[r0][k0][r1][k1]<< " " << endl;
	  }
	}
      }
    }
     //}
    
  if(NSPECIES > 1){
    //    cout << "Posterior Mean Pi:" << endl;
    OutputPi(MeanPi, RegionName, pifile, Perm );
  }

  // 
  // Output Final Frequency Estimates
  //
  // output << "Final Freqs:" << endl;
//   for(int l=0;l<NLOCI; l++){
//     for(int allele=0; allele<Nallele[l]; allele++){
//       output << setw(5) << l << "   " << Coding[l][allele] << "  ";
//       for(int r=0;r<NREGION;r++){     
// 	double p = ExpTheta[r][0][l][allele]/SumExpTheta[r][0][l];
// 	if(p<1e-4) 
// 	  p=0;
// 	output << setiosflags(ios::fixed) << setprecision(3) << setw(10) << p << " ";
//       }
//       output << endl;
//     }
//   }


  
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
	freqfile << setiosflags(ios::fixed) << setprecision(3) << setw(10) << (1.0*Count[Perm[r]][0][l][allele])/SumCount[Perm[r]][0][l] << " ";
      }
      freqfile << endl;
    }
  }

  cerr << "Finished!" << endl;

 
}
