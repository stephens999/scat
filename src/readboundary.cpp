#include "readboundary.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>    // for trunc

// SCAT version 3.0.1
// VORONOI version 2.0.0

// This file contains routines to handle grid-style species range maps
// It is used by both SCAT and VORONOI but exists separately in both
// projects for technical reasons; the two copies should be identical.
// Written by Mary Kuhner and Jon Yamato in 2020

using namespace std;

Mapgrid::Mapgrid() {
}

// Assumes ifstream is attached to a file giving the 
// in-range grid squares as latitude and longitude

void Mapgrid::Initialize(ifstream& gridfile,int margin) {
   // read the map data
   assert(gridfile.is_open());
   string line;
   double latitude, longitude;
   DoubleVec1d latitudes;
   DoubleVec1d longitudes;
   vector<pair<double, double> > gridvals; 
   while (getline(gridfile,line)) {
     istringstream iss(line);
     int intest = iss.peek();
     if (intest == EOF) {
       cerr << "  gridfile contains an empty line" << endl;
       exit(-1);
     }
     iss >> latitude;
     if (iss.fail()) {
       cerr << "  gridfile line '" << line << "' does not start with a decimal latitude" << endl;
       exit(-1);
     }
     intest = iss.peek();
     if (intest == EOF) {
       cerr << "  gridfile line '" << line << "' lacks second entry (a longitude)" << endl;
       exit(-1);
     }
     // is it either a "E","N","S" or "W"
     if (intest == 69 || intest == 78 || intest == 83 || intest == 87) { 
       cerr << "  gridfile line '" << line << "' latitudes and longitudes must be decimal" << endl;
       exit(-1);
     }
     iss >> longitude;
     if (iss.fail()) {
       cerr << "  gridfile line '" << line << "' failed to parse a decimal longitude" << endl;
       exit(-1);
     }
     intest = iss.peek();
     // is it either a "E","N","S" or "W"
     if (intest == 69 || intest == 78 || intest == 83 || intest == 87) { 
       cerr << "  gridfile line '" << line << "' latitudes and longitudes must be decimal" << endl;
       exit(-1);
     }
     if (intest != EOF) {
       cerr << "  gridfile had more than two entries on line: '" << line << "'" << endl;
       exit(-1);
     }
     latitudes.push_back(latitude);
     longitudes.push_back(longitude);
     gridvals.push_back(make_pair(latitude,longitude));
   }
   gridfile.close();
  
   // establish boundaries and size of grid 
   double minlat = *min_element(latitudes.begin(),latitudes.end());
   double maxlat = *max_element(latitudes.begin(),latitudes.end());
   double minlong = *min_element(longitudes.begin(),longitudes.end());
   double maxlong = *max_element(longitudes.begin(),longitudes.end());
   // note use of "floor" here so that a value of -1.5 ends up as -2, not -1.
   // margin is a boundary around the grid, which Voronoi needs
   m_gridminlat = floor(minlat) - margin;
   // +1 in following line because we need the square associated with
   // the maximum value to be LIVE
   m_gridmaxlat = floor(maxlat) + margin + 1;
   m_gridminlong = floor(minlong) - margin;
   // +1 here for same reason
   m_gridmaxlong = floor(maxlong) + margin + 1;
   m_gridsize_lat = m_gridmaxlat - m_gridminlat;
   m_gridsize_long = m_gridmaxlong - m_gridminlong;

   // temporary code to make the grid square, as VORONOI is not prepared
   // for a non-square grid. In the long run, should allow this.
   if (m_gridsize_lat != m_gridsize_long) {
     int maxgrid = max(m_gridsize_lat,m_gridsize_long);
     int mingrid = min(m_gridsize_lat,m_gridsize_long);
     int diff = maxgrid - mingrid;
     int leftadd = diff/2;
     int rightadd = diff - leftadd;
     if (maxgrid == m_gridsize_lat) {   // latitude is bigger
       m_gridminlong -= leftadd;
       m_gridmaxlong += rightadd;
       m_gridsize_long = m_gridmaxlong - m_gridminlong;
       assert(m_gridsize_long == m_gridsize_lat);
     } else {  // longitude is bigger
       m_gridminlat -= leftadd;
       m_gridmaxlat += rightadd;
       m_gridsize_lat = m_gridmaxlat - m_gridminlat;
       assert(m_gridsize_long == m_gridsize_lat);
     }
   }
   // DEBUG
   cout << "Grid sizes " << m_gridsize_lat << " " << m_gridsize_long << endl;

   IntVec2d tempvec(m_gridsize_lat,IntVec1d(m_gridsize_long,0));
 
   vector<pair<double, double>>::iterator latlong;
   for(latlong = gridvals.begin(); latlong != gridvals.end(); ++latlong) {
      int latitude(floor(latlong->first));
      int longitude(floor(latlong->second));
      latitude -= m_gridminlat;
      longitude -= m_gridminlong;
      if (tempvec[latitude][longitude] == 1) {
        cerr << "Found two elements in the grid file that both map to grid square";
        cerr << latitude << ", " << longitude << endl;
        cerr << "One of them is latitude " << latlong->first << " and longitude " << latlong->second << endl;
        exit(-1);
      }
      tempvec[latitude][longitude] = 1;
   }
   m_grid = tempvec;

   m_initialized = true;
}

Mapgrid::Mapgrid(ifstream& gridfile, int margin) {
  Initialize(gridfile, margin);
}  

void Mapgrid::WriteMapInfo(std::ofstream& outfile) const {
  if(!m_initialized) {
    cerr << "Attempted to write map info for uninitialized object" << endl;
    exit(-1);
  }
  outfile << "Lower left grid square at " << m_gridminlat << ", " << m_gridminlong << endl;
  outfile << "Upper right grid square at " << m_gridmaxlat << ", " << m_gridmaxlong << endl;
  outfile << "North-South extent of grid:  " << m_gridsize_lat << " squares" << endl;
  outfile << "East-West extent of grid:  " << m_gridsize_long << " squares" << endl;
}

void Mapgrid::PrintGrid(ofstream& output) const {
  // The grid is stored with 0,0 being the SW corner, but we print it with north
  // at the top, hence the reversals.
  if(!m_initialized) {
    cerr << "Attempted to write map info for uninitialized object" << endl;
    exit(-1);
  }
  for (int i = m_gridsize_lat - 1; i >= 0; --i) {
    for (int j = 0; j < m_gridsize_long; ++j) {
      output << m_grid[i][j];
    }
    output << endl;
  }
}

int Mapgrid::get_gridsize() const {
  assert(m_gridsize_lat == m_gridsize_long);
  return m_gridsize_lat;
}

pair<int,int> Mapgrid::get_latbounds() const {
  return make_pair(m_gridminlat,m_gridmaxlat);
}

pair<int,int> Mapgrid::get_longbounds() const {
  return make_pair(m_gridminlong,m_gridmaxlong);
}

int Mapgrid::lat_to_glat(double latitude) const {
  if (!m_initialized) {
    cerr << "Error: Using an uninitialized Mapgrid object" << endl;
    exit(-1);
  }
  return (int) trunc(latitude - m_gridminlat);
}

int Mapgrid::long_to_glong(double longitude) const {
  if (!m_initialized) {
    cerr << "Error: Using an uninitialized Mapgrid object" << endl;
    exit(-1);
  }
  return (int) trunc(longitude - m_gridminlong);
}

bool Mapgrid::in_range(double latitude, double longitude) const {
  if (!m_initialized) {
    cerr << "Error: Using an uninitialized Mapgrid object" << endl;
    exit(-1);
  }
  // if we are outside the grid we are not "in range"
  int glat = lat_to_glat(latitude);
  if (glat < 0 || glat >= m_gridsize_lat) return false;
  int glong = long_to_glong(longitude);
  if (glong < 0 || glong >= m_gridsize_long) return false; 
  if (m_grid[glat][glong] == 1) return true;
  else return false;
}
