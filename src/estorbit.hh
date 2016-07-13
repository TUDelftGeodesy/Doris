#ifndef ESTORB_H
#define ESTORB_H

using namespace std;

#include "constants.hh"		
#include "readinput.hh"		
#include "orbitbk.hh"		
#include "slcimage.hh"          
#include "productinfo.hh"       
#include "bk_baseline.hh"

void estorbits(const input_estorbits  &estorbitinput,
	       const input_gen        &generalinput,
	       const input_ell        &ellips,
	       const slcimage         &master,         
	       const slcimage         &slave,           
                     orbit            &masterorbit,      
	             orbit            &slaveorbit,      
	       const productinfo      &interferogram, 
	       const productinfo      &coherence,
	       const BASELINE         &baseline);

void modifyOrbits(const input_morbits &morbitsinput,
		  const input_gen     &generalinput,
		  const int16         FILEID,
		  const input_ell     &ellips,
		  const slcimage      &image,
                  const slcimage      &masterimage,
                        orbit         &thisorbit,
                        orbit         &masterorbit);

void disableOldOrbits(const char resfile[]); 

void getCoordinateFrame(cn              &ea,
			cn              &er,
			cn              &ex, 
			const real8     tAzi,
			orbit           &masterorbit);


// the following definitions do not necessarily need to be in the headerfile, 
// since they are only used in estorbit.cc and nowhere else. - at least for now. [HB]

class sigmaTable
{
private:
  real4 table[101];
public:
  sigmaTable(uint multilook);
  real4 getSigma(real4 coherence) const;
};

class observationdata
{
private:
  matrix<uint> line;      // multilooked line numbers of observations
  matrix<uint> pixel;     // multilooked pixel number os observations
  matrix<real4> phi;      // observed interferometric phases
  matrix<real4> weight;   
  matrix<real4> height;   // interpolated DEM height 
  matrix<bool> use;       // if false, observation has been rejected by outlier test (LSQ only)
  matrix<real8> residual; 
  matrix<real4> test;     // outlier test (t-test)
  uint nobs;              // current number of used observations, not counting rejected ones
  uint maxobs;            // theoretical number of observations, ignoring thresholding by EO_THRESHOLD
  const uint degree;      // degree of orbit error polynomial
  const uint npar;        // number of parameters, is: 2*(degree+1)
  uint gridsize;          // size of cells of the grid defined for selection of observagtions
  uint gridlines;         //   total number of rows of this grid in azimuth
  uint gridpixels;        //   total number of columns of this grid in range
  const bool useheights;  // true if DEM is specified and used, otherwise height is not initialised
  real8 azoffset;         // time offset of master w. r. t. reference acquisition (EO_REFORBIT)
  const real8 tmin;       // acquisition start time of reference acquisition (used for normalisation)
  const real8 tmax;       // acquisition end time of reference acquisition (used for normalisation)
  real8 t_start;          // acquisition start time w. r. t. start time of reference acquisition
  const input_ell ellips;
  const slcimage *master;
  const slcimage *slave;
  orbit *masterorbit;
  orbit *slaveorbit;
  orbit *referenceorbit;
  const productinfo interferogram;
public:
  observationdata(const input_estorbits &estorbitinput,
		  const input_gen       &generalinput,
		  const input_ell       &ellips,
		  const slcimage        *master,         
		  const slcimage        *slave,
		  const slcimage        *reference,
		        orbit           *masterorbit,      
                        orbit           *slaveorbit,      
                        orbit           *referenceorbit,
		  const productinfo     &interferogram,
		  const productinfo     &coherence);

  // constants needed for lineOfDesignMatrix
  static const int16 DESIGN_NORMAL=0;   
  static const int16 DESIGN_EXTENDED=1; 
  static const int16 DESIGN_DUMP=2;     

  matrix<real8> lineOfDesignMatrix(const uint index, int16 version);
  uint nRequested() const {return maxobs;};
  uint nInit() const {return phi.lines();};
  uint nUsed() const {return nobs;};
  real8 operator () (const uint i) const {return phi(i,0);};
  real8 getWeight(const uint i) const {return weight(i,0);};    // real8 necessary?
  void reweight(const uint i, const real4 newWeight) {weight(i,0) = newWeight;};
  void setRes(const uint i, const real8 v) {residual(i,0) = v;};
  real8 getRes(const uint i) const {return residual(i,0);};
  uint getLine(const uint i) const {return line(i,0);};
  uint getPixel(const uint i) const {return pixel(i,0);};
  void reject(const uint i) {use(i,0) = false; nobs--;};
  bool isUsed(const uint i) const {return use(i,0);};
  void setT(const uint i, const real4 T) {test(i,0) = T;};
  real4 getT(const uint i) const {return test(i,0);};
  real4 getHeight(const uint i) const {return useheights ? height(i,0) : 0;};
  uint getGridSize() const {return gridsize;};
  uint getGridLines() const {return gridlines;};
  uint getGridPixels() const {return gridpixels;};
};

void repairMatrix(matrix<real8> &M);

real8 getEigenTheta(const matrix<real8> &Q);

void constrainSolution(matrix<real8> &Qxx, matrix<real8> &rhs, const matrix<int16> &constraints);

void readOrbitDataPoints(const char resfile[], 
			 matrix<real8> &data);

real8 getAzOffset(const input_ell &ellips,
		  const slcimage  &master,
		  const slcimage  &slave,
                        orbit     &masterorbit,
                        orbit     &slaveorbit);

void writeResiduals(const input_estorbits &estorbitsinput,
		    const observationdata &obs);



#endif // ESTORB_H
