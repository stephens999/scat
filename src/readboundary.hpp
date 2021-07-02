#ifndef READBOUNDARY_HPP
#define READBOUNDARY_HPP
#include <iosfwd>
#include <vector>
#include <iostream>
#include <fstream>

// SCAT version 3.0.1
// VORONOI version 2.0.0

typedef std::vector<double> DoubleVec1d;
typedef std::vector<DoubleVec1d> DoubleVec2d;
typedef std::vector<int> IntVec1d;
typedef std::vector<IntVec1d> IntVec2d;

class Mapgrid
{
public:
  Mapgrid();
  Mapgrid(std::ifstream& gridfile, int margin);

  void Initialize(std::ifstream& gridfile, int margin);
  bool in_range(double latitude, double longitude) const;
  int get_gridsize_lat() const { return m_gridsize_lat; };
  int get_gridsize_long() const { return m_gridsize_long; };
  int get_gridsize() const;  // assumes it's square, others do not!
  std::pair<int,int> get_latbounds() const;
  std::pair<int,int> get_longbounds() const;
  void WriteMapInfo(std::ofstream& outfile) const;
  void PrintGrid(std::ofstream& output) const;   

private:
  Mapgrid(const Mapgrid&);              // deliberately disabled
  Mapgrid& operator=(Mapgrid const&);   // deliberately disabled

  IntVec2d m_grid;
  int lat_to_glat(double latitude) const;
  int long_to_glong(double longitude) const;

  // minimum and maximum entries in data
  int m_gridminlat = -999;
  int m_gridminlong = -999;
  int m_gridmaxlat = -999;
  int m_gridmaxlong = -999;
  int m_gridsize_lat = -999;
  int m_gridsize_long = -999;

  bool m_initialized = false;
  
};
#endif // READBOUNDARY_HPP
