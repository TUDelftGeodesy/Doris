#include "constants.hh"       
#include "ioroutines.hh"      
#include "orbitbk.hh"         
#include "slcimage.hh"        
#include "productinfo.hh"     
#include "exceptions.hh"      
#include "bk_baseline.hh"     
#include "coregistration.hh"  
#include "estorbit.hh"

#include <iomanip>      
#include <cstdlib>      
#include <cmath>        
#include <algorithm>    
#include <cstdio>        



/****************************************************************
 *    estorbits                                                 *
 *                                                              *
 * estimate baseline error polynomials in horizontal (Bh) and   *
 * vertical (Bv) baseline. If d is the polynomial degree, the   *
 * output (coeff) is a 2*(d+1) x 1 vector containing the        *
 * following baseline error coefficients:                       *
 *                                                              *
 * (dBh0;dBh1;...;dBv0;dBv1;...)                                *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void estorbits(const input_estorbits &estorbitsinput,
	       const input_gen       &generalinput,
	       const input_ell       &ellips,
	       const slcimage        &master,         
	       const slcimage        &slave,           
                     orbit           &masterorbit,      
	             orbit           &slaveorbit,      
	       const productinfo     &interferogram, 
	       const productinfo     &coherence,
	       const BASELINE        &baseline)
{
  TRACE_FUNCTION("estorbits (HB 17-Mar-2011)");

  // ______ check if cutout windows of interferogram and coherence image are congruent ______
  if (!(specified(estorbitsinput.ifpositions) && estorbitsinput.weighting==eo_noweighting && estorbitsinput.threshold==0))
    if (interferogram.win.linelo != coherence.win.linelo ||
	interferogram.win.pixlo != coherence.win.pixlo ||
	int(interferogram.win.lines()/interferogram.multilookL) != 
	int(coherence.win.lines()/coherence.multilookL) ||
	int(interferogram.win.pixels()/interferogram.multilookP) != 
	int(coherence.win.pixels()/coherence.multilookP) ||
	interferogram.multilookL != coherence.multilookL ||
	interferogram.multilookP != coherence.multilookP)
      // direct comparison of the two windows is not possible, since linehi or pixhi
      // might differ although lines() and pixels() are the same.
      { 
	ERROR << "Interferogram and coherence image are not congruent "
	      << "(same cutout window, same multilook factors), which is required "
	      << "for step ESTORBITS." << endl
	      << "Int: ln(" << interferogram.win.linelo << "," << interferogram.win.linehi 
	      << "), px(" << interferogram.win.pixlo << "," << interferogram.win.pixhi 
	      << "), size(" << interferogram.multilookL << "," << interferogram.multilookP 
	      << ")" << endl
	      << "Coh : ln(" << coherence.win.linelo << "," << coherence.win.linehi 
	      << "), px(" << coherence.win.pixlo << "," << coherence.win.pixhi 
	      << "), size(" << coherence.multilookL << "," << coherence.multilookP 
	      << ")" << endl;
	PRINT_ERROR(ERROR.get_str())
	  throw(some_error);
      }

  const real4 meanLine  = (master.currentwindow.linelo+master.currentwindow.linehi)/2;
  const real4 meanPix   = (master.currentwindow.pixlo+master.currentwindow.pixhi)/2;
  const int16 degree    = estorbitsinput.poldegree; // degree of polynomial to be estimated
  const int16 npar      = 2*(degree+1);             // number of parameters in polynomial (*2 for Bh and Bv)
  real8 theta;                                      // mean look angle
  matrix<real8> Bh, Bv, Bpar, Bperp;                // baseline parameters

  // ______ choose reference ______
  const slcimage *reference;
  orbit *referenceorbit;
  if (specified(estorbitsinput.reforbitfile))
    {
      slcimage *tempptr = new slcimage();
      (*tempptr).fillslcimage(estorbitsinput.reforbitfile);
      reference = tempptr;
      referenceorbit = new orbit();
      (*referenceorbit).set_interp_method(generalinput.orb_interp);
      (*referenceorbit).initialize(estorbitsinput.reforbitfile);
      INFO << "Using reference orbit from file: " << estorbitsinput.reforbitfile;
      INFO.print();
    }
  else
    {
      INFO.print("Using master orbit as reference orbit.");
      reference = &master;
      referenceorbit = &masterorbit;
    }

  // ______ read observation data ______
  observationdata obs(estorbitsinput,generalinput,ellips,&master,&slave,reference,
		      &masterorbit,&slaveorbit,referenceorbit,interferogram,coherence);

  // ______ check estimability and compute redundancy ______
  uint minimum_number_of_obs = npar+1;
  if (estorbitsinput.method==eo_gridsearch) minimum_number_of_obs -= estorbitsinput.constraints.lines();
  // For the sake of an easier implementation, estimability of the unconstrained solution 
  // is required for method LSQ, even if the solution is constrained.
  if (obs.nUsed()<minimum_number_of_obs)
    {
      ERROR << "Only " << obs.nUsed() << " pixels found with coherence > " << estorbitsinput.threshold << ".";
      ERROR.print();
      ERROR << "At least " << minimum_number_of_obs << " observations are required.";
      ERROR.print();
      throw(unhandled_case_error);
    }
  int32 redundancy = obs.nUsed()-npar-1+estorbitsinput.constraints.lines();  
  // "minus one" for the eliminated constant phase shift parameter
  DEBUG << "Initial redundancy: " << redundancy;
  DEBUG.print();
  real4 averageredundancy = real4(redundancy)/obs.nUsed();
  DEBUG << "Average redundancy number per observation (0...1): " 
	<< setprecision(8) << averageredundancy;
  DEBUG.print();
  if (averageredundancy<0.5)  
    WARNING.print("Reliability of estimates is questionable due to poor redundancy.");
  

  // ====== frist output to logfile ======
  ofstream scratchlogfile("scratchlogestorbit", ios::out | ios::trunc);
  bk_assert(scratchlogfile,"estorbit: scratchlogestorbit",__FILE__,__LINE__);
  scratchlogfile.setf(ios_base::fixed);
  scratchlogfile 
    << "\n\n*******************************************************************"
    << "\n* ORBIT ERROR ESTIMATION"
    << "\n*******************************************************************"
    << "\nMethod:                                 ";
  switch (estorbitsinput.method)
    {
    case eo_lsq: 
      scratchlogfile << "Least Squares Adjustment";
      break;
    case eo_gridsearch: 
      scratchlogfile << "Gridsearch";
      break;
    }
  scratchlogfile
    << "\nPolynomial degree:                      " << estorbitsinput.poldegree
    << "\nConstraints:                            ";
  if (estorbitsinput.constrained)
      for (uint i=0; i<estorbitsinput.constraints.lines(); i++)
	{
	  if (estorbitsinput.constraints(i,0)==0) scratchlogfile << "BPAR";
	  else if (estorbitsinput.constraints(i,0)==1) scratchlogfile << "BPERP";
	  scratchlogfile << setw(1) << estorbitsinput.constraints(i,1) << " ";
	}
  else
    scratchlogfile << "none"; 
  scratchlogfile << "\nInterpolated DEM:                       ";
  if (specified(estorbitsinput.fiheightmap))
    scratchlogfile << estorbitsinput.fiheightmap;
  else
    scratchlogfile << "not used";
  scratchlogfile  << "\nReference orbit:                        ";
  if (specified(estorbitsinput.reforbitfile))
    scratchlogfile << estorbitsinput.reforbitfile;
  else 
    scratchlogfile << "master orbit";
  scratchlogfile << "\nWeighting scheme                        ";
  switch (estorbitsinput.weighting)
    {
    case eo_noweighting: scratchlogfile << "none";                    break;
    case eo_coh:         scratchlogfile << "coherence";               break;
    case eo_coh2:        scratchlogfile << "coherence^2";             break;
    case eo_pdf:         scratchlogfile << "theoretical probability"; break;
    }
  scratchlogfile << "\nCoherence threshold:                    " 
		 << setprecision(3) << estorbitsinput.threshold;
  if (specified(estorbitsinput.ifpositions))
    scratchlogfile << "\nObservation locations from file:        " 
		   << estorbitsinput.ifpositions;
  else
    {
      scratchlogfile
	<< "\nGridwidth:                              " << obs.getGridSize()
	<< "\nGridrows:                               " << obs.getGridLines()
	<< "\nGridcolumns:                            " << obs.getGridPixels()
	<< "\nObservation pixels requested            " << estorbitsinput.nobs;
    }
  scratchlogfile
    << "\nMaximum number of available pixels:     " << obs.nRequested()
    << "\nCoherent observation pixels found:      " << obs.nInit();

  
  // ====== compute conversion factors ======
  cn M, P;
  M = masterorbit.getxyz(master.line2ta(meanLine));
  lp2xyz(meanLine,master.currentwindow.pixhi,ellips,master,masterorbit,P);
  real8 deltaTheta = acos((-M.normalize()).in((P-M).normalize()));
  lp2xyz(meanLine,master.currentwindow.pixlo,ellips,master,masterorbit,P);
  deltaTheta -= acos((-M.normalize()).in((P-M).normalize()));

  const real8 rateConversion = 4. / (reference->originalwindow.lines()-1) * reference->prf;
  const real8 dBperp_one_fringe = -(master.wavelength+slave.wavelength)/4/deltaTheta;
  const real8 timeOfAcquisition = (master.currentwindow.lines()-1)/master.prf;;
  const real8 dBparRate_one_fringe = -(master.wavelength+slave.wavelength)/4/timeOfAcquisition;


  // ====== setup normal equations, compute mean theta  ======
  matrix<real8> NEQ(npar+1,npar+1);      // bordered normal equation matrix  [A'*P*A , A'*P*y ; y'*P*A , y'*P*y]
  matrix<real8> NEQc(npar+1,npar+1);     // bordered normal equation matrix after Cholesky decomposition
  matrix<real8> Qxx(npar,npar);          // cofactor matrix of estimated polynomial coefficients
  matrix<real8> schreibersums(1,npar+1); // columnsums of [A,l]
  matrix<real8> rhs(npar,1);             // estimated polynomial coefficients ("right hand side")
  real8 weightsum=0;                     // sum of al weights

  NEQ = 0;   
  schreibersums = 0;

  // ______ stack normal equations ______
  for (uint i=0; i<obs.nInit(); i++)
    {
      matrix<real8> Ai = obs.lineOfDesignMatrix(i,observationdata::DESIGN_EXTENDED);
      NEQ += obs.getWeight(i) * matTxmat(Ai,Ai);
      schreibersums += obs.getWeight(i) * Ai;
      weightsum += obs.getWeight(i);
    }

  NEQ -= matTxmat(schreibersums,schreibersums)/weightsum;  // apply Schreiber's sum equation
  // This accounts for constant phase shift parameter phi_0.

  // ______ compute mean theta ______
  if (redundancy-int32(estorbitsinput.constraints.lines()) >= 0)  
    // estimability of unconstrained solution required for determination of theta
    { 
      // invert normal equation matrix and compute mean theta from the eigenspaces
      rhs = NEQ.getdata(window(0,npar-1,npar,npar));
      NEQc = NEQ;
      choles(NEQc);
      Qxx = NEQc.getdata(window(0,npar-1,0,npar-1));
      solvechol(Qxx, rhs);
      invertchol(Qxx);
      repairMatrix(Qxx);
      theta = getEigenTheta(Qxx);
    }
  else
    {
      // define theta as the look angle to the scene centre
      cn M, P;
      M = masterorbit.getxyz(master.line2ta(meanLine));
      lp2xyz(meanLine,meanPix,ellips,master,masterorbit,P);
      theta = acos((-M.normalize()).in((P-M).normalize()));
    }


  // ====== method GRIDSEARCH ======
  if (estorbitsinput.method==eo_gridsearch)
    {
      PROGRESS.print("Selected method: GRIDSEARCH.");

      // some empirical parameters for iterative refinement of search space
      const uint numberOfIterations = 6;  
      real4 step[numberOfIterations] = {0.8,0.2,0.05,0.02,0.008,0.003}; // gridwidth [fringes]
      uint  zoom[numberOfIterations-1] = {7,10,2,2,2}; 

      // define initial searchspace
      real4 dBparRateUpperBound = abs(estorbitsinput.maxfringesaz*dBparRate_one_fringe)/rateConversion;
      real4 dBperpUpperBound = abs(estorbitsinput.maxfringesrg*dBperp_one_fringe);
      real4 dBparRateLowerBound = -dBparRateUpperBound;
      real4 dBperpLowerBound = -dBperpUpperBound;

      // some parameters for reliability assessment
      real4 maxGamma;                    // maximum coherence value
      real4 ratioOfTwoHighestMaxima;     
      real4 distanceOfTwoHighestMaxima;  //  [fringes]

      // ______ initialisation ______
      Bpar = matrix<real8>(2,1);
      Bpar = 0;
      Bperp = matrix<real8>(2,1);
      Bperp = 0;
      matrix<real8> T(4,2);
      T = 0;
      T(0,1) = cos(theta);
      T(1,0) = sin(theta);
      T(2,1) = sin(theta);
      T(3,0) = -cos(theta);


      // ====== iteration loop ======
      for (uint it=0;;it++)
	{
	  PROGRESS << "Iteration " << it+1 << " of 6.";
	  PROGRESS.print();
	  real4 dBparRateStep = step[it]*abs(dBparRate_one_fringe)/rateConversion;
	  real4 dBperpStep = step[it]*abs(dBperp_one_fringe);
	  uint n1 = ceil((dBparRateUpperBound-dBparRateLowerBound)/dBparRateStep)+1;
	  uint n2 = ceil((dBperpUpperBound-dBperpLowerBound)/dBperpStep)+1;
	  real4 dBparRateGrid[n1], dBperpGrid[n2];
	  matrix<complr4> GAMMA(n1,n2);

	  // ______ initialise search grid ______
	  GAMMA=0;
	  dBparRateGrid[0] = dBparRateLowerBound;
	  dBperpGrid[0]    = dBperpLowerBound;
	  for (uint i=1;i<n1;i++)
	    dBparRateGrid[i] = dBparRateGrid[i-1] + dBparRateStep;
	  for (uint i=1;i<n2;i++)
	    dBperpGrid[i] = dBperpGrid[i-1] + dBperpStep;

	  // ______ compute gamma for the whole grid ______
	  for (uint i=0; i<obs.nInit(); i++)
	    {
	      matrix<real8> Ai = obs.lineOfDesignMatrix(i,observationdata::DESIGN_NORMAL);
	      Ai -= schreibersums.getdata(window(0,0,0,npar-1)) / weightsum;
	      Ai = Ai*T;
	      
	      for (uint i1=0;i1<n1;i1++)
		{
		  matrix<real8> x(2,1);
		  x(0,0) = dBparRateGrid[i1];
		  for (uint i2=0;i2<n2;i2++)
		    {
		      x(1,0) = dBperpGrid[i2];
		      real8 phase = obs(i)-(Ai*x)(0,0);
		      GAMMA(i1,i2) += complr4(obs.getWeight(i)*cos(phase),obs.getWeight(i)*sin(phase));  
		    }
		}
	    }

	  // ______ find maximum ______
	  uint max1, max2;
	  matrix<real4> GAMMAabs(n1,n2);
 	  for (uint i1=0;i1<n1;i1++)
 	    for (uint i2=0;i2<n2;i2++)
 	      GAMMAabs(i1,i2) = abs(GAMMA(i1,i2));
	  max(GAMMAabs,max1,max2);
	  maxGamma = GAMMAabs(max1,max2)/obs.nUsed();

	  if (max1<n1 && max2<n2)
	    {
	      Bpar(1,0) = dBparRateGrid[max1];
	      Bperp(0,0) = dBperpGrid[max2];
	    }
	  else 
	    {
	      ERROR.print("Method GRIDSEARCH: No maximum found.");
	      throw(some_error);
	    }

	  // ______ compute reliability indicators ______
	  if (it==1)
	    {
	      // find second-highest local maximum by setting the environments of (max1,max2)
	      // to -1 as long as the gradient is negative. Then apply max-operator again.
	      real4 maxGAMMA = GAMMAabs(max1,max2);;
	      GAMMAabs(max1,max2) = -1;
	      while (true)
		{
		  uint max1b, max2b;
		  max(GAMMAabs,max1b,max2b);
		  if(max1b>=n1 || max2b>=n2 || GAMMAabs(max1b,max2b)==-1)
		    {
		      ratioOfTwoHighestMaxima = NaN;
		      distanceOfTwoHighestMaxima = NaN;
		      break;
		    }
		  if ( (max1b>0 && max2b>0 && GAMMAabs(max1b-1,max2b-1)==-1) ||
		       (max1b>0 && GAMMAabs(max1b-1,max2b)==-1) ||
		       (max1b>0 && max2b<(n2-1) && GAMMAabs(max1b-1,max2b+1)==-1) ||
		       (max2b>0 && GAMMAabs(max1b,max2b-1)==-1) ||
		       (max2b<(n2-1) && GAMMAabs(max1b,max2b+1)==-1) ||
		       (max1b<(n1-1) && max2b>0 && GAMMAabs(max1b+1,max2b-1)==-1) ||
		       (max1b<(n1-1) && GAMMAabs(max1b+1,max2b)==-1) ||
		       (max1b<(n1-1) && max2b<(n2-1) && GAMMAabs(max1b+1,max2b+1)==-1) )
		    GAMMAabs(max1b,max2b) = -1;
		  else
		    {
		      ratioOfTwoHighestMaxima = maxGAMMA/GAMMAabs(max1b,max2b);
		      distanceOfTwoHighestMaxima = (abs(int(max1-max1b))+abs(int(max2-max2b)))*step[it];  // [fringes]
		      INFO << "Ratio of highest and second-highest local maximum: " << ratioOfTwoHighestMaxima;
		      INFO.print();
		      if (ratioOfTwoHighestMaxima<1.5)
			WARNING.print("No distinct local maximum found (ratio<1.5). There is a high chance of unreliable estimates.");
		      break;
		    }
		}
	    }
	  
	  // ______ refine search space ______
	  if (it==numberOfIterations-1) break;
	  dBparRateUpperBound = Bpar(1,0)    + zoom[it] * step[it] * abs(dBparRate_one_fringe) / rateConversion;
	  dBperpUpperBound    = Bperp(0,0)   + zoom[it] * step[it] * abs(dBperp_one_fringe);
	  dBparRateLowerBound = 2*Bpar(1,0)  - dBparRateUpperBound;
	  dBperpLowerBound    = 2*Bperp(0,0) - dBperpUpperBound;
	}
      if (strcmp(estorbitsinput.foresiduals," "))  // output if filename has been specified by EO_RESIDUALS
	writeResiduals(estorbitsinput,obs);
	  
      // ______ convert (Bpar1,Bperp) -> (Bh0,Bh1,Bv0,Bv1) ______
      Bh = matrix<real8>(2,1);
      Bh(0,0) = cos(theta)*Bperp(0,0);
      Bh(1,0) = sin(theta)*Bpar(1,0);
      Bv = matrix<real8>(2,1);
      Bv(0,0) = sin(theta)*Bperp(0,0);
      Bv(1,0) = -cos(theta)*Bpar(1,0);


      // ====== second output to logfile (gridsearch) ======
      scratchlogfile 
	<< "\nSearchspace in azimuth [fringes]:       " << estorbitsinput.maxfringesaz
	<< "\nSearchspace in range [fringes]:         " << estorbitsinput.maxfringesrg
	<< "\nOrbit coherence maximum (gridsearch)    " << setprecision(3) << maxGamma
 	<< "\nRatio of two highest local maxima:      " << setprecision(3) << ratioOfTwoHighestMaxima
 	<< "\nDistance of 2 highest maxima [fringes]: " << setprecision(1) << distanceOfTwoHighestMaxima;

      // ______ write baseline error parameters ______
      scratchlogfile.precision(8);
      scratchlogfile.width(13);
      scratchlogfile << "\n\nEstimated coefficients:\n"
		     << "\ndBh_hat\n";
      for (uint i=0; i<=degree; i++)
	{
	  if (Bh(i,0) >= 0.) scratchlogfile << " ";
	  scratchlogfile <<  Bh(i,0) << endl;
	}
      scratchlogfile << "\ndBv_hat\n";
      for (uint i=0; i<=degree; i++)
	{
	  if (Bv(i,0) >= 0.) scratchlogfile << " ";
 	  scratchlogfile <<  Bv(i,0) << endl;
 	}
      scratchlogfile << "\ndBpar_hat\n";
      for (uint i=0; i<=degree; i++)
	{
	  if (Bpar(i,0) >= 0.) scratchlogfile << " ";
	  scratchlogfile <<  Bpar(i,0) << endl;
	}
      scratchlogfile << "\ndBperp_hat\n";
      for (uint i=0; i<=degree; i++)
 	{
	  if (Bperp(i,0) >= 0.) scratchlogfile << " ";
 	  scratchlogfile <<  Bperp(i,0) << endl;
 	}

    } // END iteration loop


  // ====== method LSQ ======
  else 
    {
      PROGRESS.print("Selected method: LSQ.");
      if (estorbitsinput.constrained) 
	constrainSolution(Qxx,rhs,estorbitsinput.constraints);
      
      // ______ compute residuals and iterate for outlier removal  ______
      uint iteration=0, iTmax;
      real8 Omega=0;
      
      // ______ intitial computation of Omega (weighted sum of squared residuals) ______
      for (uint i=0; i<obs.nInit(); i++)
	{
	  matrix<real8> Ai = obs.lineOfDesignMatrix(i,observationdata::DESIGN_NORMAL);
	  Ai -= schreibersums.getdata(window(0,0,0,npar-1)) / weightsum;
	  matrix<real8> vi = Ai*rhs;
	  obs.setRes(i, vi(0,0)-(obs(i)-schreibersums(0,npar)/weightsum));
	  Omega += obs.getRes(i)*obs.getWeight(i)*obs.getRes(i);
	}

      // ______ prepare logfile for data snooping ______
     scratchlogfile << "\n\nData snooping, rejected observations:\n"
		     << "\n__line _pixel residual ____test" << endl;

      
      // ====== data snooping loop (iterative outlier detection and rejection) ======
      real4 Tmax;
      while (true)
	{
	  Tmax = 0;
	  // ______ computation of residuals, redundancy numbers and tests ______
	  for (uint i=0; i<obs.nInit(); i++)
	    if (obs.isUsed(i))
	      {
		matrix<real8> Ai = obs.lineOfDesignMatrix(i,observationdata::DESIGN_NORMAL);								
		Ai -= schreibersums.getdata(window(0,0,0,npar-1)) / weightsum;
		matrix<real8> vi = Ai*rhs;
		obs.setRes(i, vi(0,0)-(obs(i)-schreibersums(0,npar)/weightsum));
		
		real8 qvivi  = 1/obs.getWeight(i)-matxmatT(Ai*Qxx,Ai)[0][0];
		real8 redund = qvivi * obs.getWeight(i);
		
		real8 sigma_test_sqr = (Omega-obs.getRes(i)*obs.getRes(i)/qvivi)/(redundancy-1);
		obs.setT(i, obs.getRes(i)/sqrt(sigma_test_sqr*qvivi));
		
		if (abs(obs.getT(i))>abs(Tmax))
		  {
		    Tmax = obs.getT(i);
		    iTmax = i;
		  }
	      } 

	  // ______ check if residuals have been computed correctly ______
	  if (!estorbitsinput.constrained)
	    {
	      DEBUG << "sigma from residuals: " 
		    << setprecision(15) << sqrt(Omega/redundancy);
	      DEBUG.print();
	      DEBUG << "sigma from cholesky:  " 
		    << setprecision(15) << NEQc(npar,npar)/sqrt(real8(redundancy));
	      DEBUG.print();
	    }
	  
	  // ______ check if more iterations are required ______
	  if (iteration==estorbitsinput.maxiter || redundancy==0)
	    { 
	      PROGRESS << "Maximum number of iterations reached." << ends;
	      PROGRESS.print();
	      break;
	    }
	  else if (abs(obs.getT(iTmax))<estorbitsinput.k_alpha)
	    {
	      PROGRESS << iteration << " iteration(s) done, Tmax=" << obs.getT(iTmax) << ".";
	      PROGRESS.print();
	      break;
	    }
	  
	  // ______ reject observation with highest test ______
	  // remove its contribution from the normal equations
	  obs.reject(iTmax);
	  matrix<real8> Ai = obs.lineOfDesignMatrix(iTmax,observationdata::DESIGN_EXTENDED);
	  NEQ += matTxmat(schreibersums,schreibersums)/weightsum;
	  NEQ -= obs.getWeight(iTmax) * matTxmat(Ai,Ai);
	  schreibersums -= obs.getWeight(iTmax) * Ai;
	  weightsum -= obs.getWeight(iTmax);
	  NEQ -= matTxmat(schreibersums,schreibersums)/weightsum;
	  Omega -= obs.getRes(iTmax)*obs.getWeight(iTmax)*obs.getRes(iTmax);
	  redundancy--;

	  // ______ report rejection to logfile ______	  
 	  scratchlogfile << setw(6) << setprecision(0) << obs.getLine(iTmax)+1 << " "
 			 << setw(6) << setprecision(0) << obs.getPixel(iTmax)+1 << " "
 			 << setw(8) << setprecision(3) << obs.getRes(iTmax) << " "
 			 << setw(8) << setprecision(3) << obs.getT(iTmax) << " "
			 << endl;

 	  DEBUG << "Iteration " << ++iteration << ", rejecting (l,p)=(" << obs.getLine(iTmax)+1 
		<< "," << obs.getPixel(iTmax)+1 << "), T=" << setprecision(4) << obs.getT(iTmax) 
		<< ", v=" << setprecision(4) << obs.getRes(iTmax) << ends;
 	  DEBUG.print();
 	  if (iteration%100==0)
 	    {
 	      PROGRESS << "Outlier rejection, iteration " << setw(3) << iteration
		       << " of max. " << estorbitsinput.maxiter << ".";
 	      PROGRESS.print();
 	    }
	  
	  // ______ re-estimation with updated normal equations ______
	  rhs = NEQ.getdata(window(0,npar-1,npar,npar));
	  NEQc = NEQ;
	  choles(NEQc);
	  Qxx = NEQc.getdata(window(0,npar-1,0,npar-1));
	  solvechol(Qxx, rhs);
	  invertchol(Qxx);
	  repairMatrix(Qxx);
	  theta = getEigenTheta(Qxx);
	  if (estorbitsinput.constrained) 
	    constrainSolution(Qxx,rhs,estorbitsinput.constraints);
	  
	} // END data snooping loop
      if (strcmp(estorbitsinput.foresiduals," "))  // output if filename has been specified by EO_RESIDUALS
	writeResiduals(estorbitsinput,obs);
	  
      // ______ compute baseline components ______
      Bh = rhs.getdata(window(0,degree,0,0));
      Bv = rhs.getdata(window(degree+1,npar-1,0,0));
      Bpar  = Bh * sin(theta) - Bv * cos(theta);
      Bperp = Bh * cos(theta) + Bv * sin(theta);
      
      // ______ error propagation for (par,perp) baseline components ______
      matrix<real8> Rot(npar,npar);
      for (uint i=0; i<=degree; i++)
	{
	  uint i2 = i + degree + 1;
	  Rot(i ,i ) =  sin(theta);
	  Rot(i ,i2) = -cos(theta);
	  Rot(i2, i) =  cos(theta);
	  Rot(i2,i2) =  sin(theta); 
	}
      matrix<real8> Qpp = matxmatT(Rot*Qxx, Rot);
      

      // ====== second output to logfile (lsq) ======     
      real8 sigma = sqrt(Omega/redundancy);
      scratchlogfile 
	<< "\nRejected observation pixels:            " << iteration
	<< "\nObservation pixels finally used:        " << obs.nUsed()
	<< "\nRedundancy:                             " << redundancy
	<< "\nMaximum test:                           " << setprecision(2) << abs(Tmax)
	<< "\nCritical value:                         " << setprecision(2) << estorbitsinput.k_alpha
	<< "\nEstimated phase standard deviation      " << setprecision(3) << sigma;

      // ______ write baseline error parameters and standard deviations ______
      scratchlogfile.precision(8);
      scratchlogfile.width(13);
      scratchlogfile << "\n\nEstimated coefficients:\n"
		     << "\ndBh_hat \tstd.\n";
      for (uint i=0; i<=degree; i++)
	{
	  if (Bh(i,0) >= 0.) scratchlogfile << " ";
	  scratchlogfile <<  Bh(i,0) << " \t" << sigma*sqrt(Qxx(i,i)) << endl;
	}
      scratchlogfile << "\ndBv_hat \tstd.\n";
      for (uint i=0; i<=degree; i++)
	{
 	  uint i2 = i + degree + 1;
	  if (Bv(i,0) >= 0.) scratchlogfile << " ";
 	  scratchlogfile <<  Bv(i,0) << " \t" << sigma*sqrt(Qxx(i2,i2)) << endl;
 	}
      scratchlogfile << "\ndBpar_hat \tstd.\n";
      for (uint i=0; i<=degree; i++)
	{
	  if (Bpar(i,0) >= 0.) scratchlogfile << " ";
	  scratchlogfile <<  Bpar(i,0) << " \t" << sigma*sqrt(Qpp(i,i)) << endl;
	}
      scratchlogfile << "\ndBperp_hat \tstd.\n";
      for (uint i=0; i<=degree; i++)
 	{
 	  uint i2 = i + degree + 1;
	  if (Bperp(i,0) >= 0.) scratchlogfile << " ";
 	  scratchlogfile <<  Bperp(i,0) << " \t" << sigma*sqrt(Qpp(i2,i2)) << endl;
 	}
      
      // ______ write covariance matrix of baseline error parameters ______
      scratchlogfile << "\nCovariance matrix of estimated parameters (h/v):"
		     << "\n------------------------------------------------\n"
		     << setiosflags(ios_base::scientific);
      real8 sigma2 = sigma*sigma;
      for (uint i=0; i<npar; i++)
	{
	  for (uint j=0; j<npar; j++)
	    scratchlogfile 
	      << setw(13) << setprecision(7)
	      << sigma2*Qxx(i,j) << " ";
	  scratchlogfile << endl;
	}
      scratchlogfile << resetiosflags(ios_base::scientific);
    } // END method switch


  // ====== output to resultfile ======
  ofstream scratchresfile("scratchresestorbit", ios::out | ios::trunc);
  bk_assert(scratchresfile,"estorbit: scratchresestorbit",__FILE__,__LINE__);
  scratchresfile.setf(ios::right, ios::adjustfield);
  scratchresfile
    << "\n\n*******************************************************************"
    << "\n*_Start_" << processcontrol[pr_i_estorbits]
    << "\n*******************************************************************"
    << "\nMethod:                                 ";
  switch (estorbitsinput.method)
    {
    case eo_lsq: 
      scratchresfile << "Least Squares Adjustment";
      break;
    case eo_gridsearch: 
      scratchresfile << "Gridsearch";
      break;
    }
  scratchresfile << "\nReference orbit:                        ";
  if (specified(estorbitsinput.reforbitfile))
    scratchresfile << estorbitsinput.reforbitfile;
  else
    scratchresfile << "master orbit";
  scratchresfile
    << "\nMean look angle [deg]:                  " << setprecision(6) << rad2deg(theta)
    << endl << endl;

  scratchresfile << "Degree_eo:      " << estorbitsinput.poldegree;
  scratchresfile.precision(8);
  scratchresfile.flags(ios_base::fixed);
  scratchresfile << "\nEstimated_coefficients_horizontal:\n";
  for (uint i=0; i<=degree; i++)
    {
      if (Bh(i,0) >= 0.) scratchresfile << " ";
      scratchresfile << Bh(i,0) << endl;
    }
  scratchresfile << "\nEstimated_coefficients_vertical:\n";
  for (uint i=0; i<=degree; i++)
    {
      if (Bv(i,0) >= 0.) scratchresfile << " ";
      scratchresfile << Bv(i,0) << endl;
    }
  scratchresfile << "\nEstimated_coefficients_parallel:\n";
  for (uint i=0; i<=degree; i++)
    {
      if (Bpar(i,0) >= 0.) scratchresfile << " ";
      scratchresfile << Bpar(i,0) << endl;        
    }
  scratchresfile << "\nEstimated_coefficients_perpendicular:\n";
  for (uint i=0; i<=degree; i++)
    {
      if (Bperp(i,0) >= 0.) scratchresfile << " ";
      scratchresfile << Bperp(i,0) << endl; 
    }
  
  scratchresfile
    << "*******************************************************************"
    << "\n* End_" << processcontrol[pr_i_estorbits] << "_NORMAL"
    << "\n*******************************************************************\n";
  scratchresfile.close();

      
  // ====== final output to logfile ======
  scratchlogfile 
    << "\nMean look angle [deg]:                  " << setprecision(4) << rad2deg(theta)
    << "\nEstimated error in azimuth [fringes]:   " << setprecision(1) << Bpar(1,0)/dBparRate_one_fringe*rateConversion
    << "\nEstimated error in range [fringes]:     " << setprecision(1) << Bperp(0,0)/dBperp_one_fringe 
    << "\ndBpar/dt conversion factor to m/s:      " << setprecision(8) << rateConversion
    << "\ndBpar/dt equ. one fringe in azimuth:    " << setprecision(5) << dBparRate_one_fringe/rateConversion
    << "\ndBpar/dt equ. one fringe in az. [m/s]:  " << setprecision(5) << dBparRate_one_fringe
    << "\ndBperp equ. one fringe in range [m]:    " << setprecision(3) << dBperp_one_fringe 
    << endl
    << "\nUse_the_following_cards_to_correct_for_the_estimated_orbit_error___"
    << "\nPROCESS            M_MORBITS";
  for (uint i=0; i<=degree; i++)
    scratchlogfile << "\nM_MO_DBH" << setw(1) << i << "          " 
		   << setprecision(8) << -Bh(i,0)/2;
  for (uint i=0; i<=degree; i++)
    scratchlogfile << "\nM_MO_DBV" << setw(1) << i << "          " 
		   << setprecision(8) << -Bv(i,0)/2;
  if (specified(estorbitsinput.reforbitfile))
    scratchlogfile << "\nM_MO_REFORBIT      " << estorbitsinput.reforbitfile;
  scratchlogfile << "\nPROCESS            S_MORBITS";
  for (uint i=0; i<=degree; i++)
    scratchlogfile << "\nS_MO_DBH" << setw(1) << i << "          " 
		   << setprecision(8) << Bh(i,0)/2;
  for (uint i=0; i<=degree; i++)
    scratchlogfile << "\nS_MO_DBV" << setw(1) << i << "          " 
		   << setprecision(8) << Bv(i,0)/2;
  if (specified(estorbitsinput.reforbitfile))
    scratchlogfile << "\nS_MO_REFORBIT      " << estorbitsinput.reforbitfile;
  scratchlogfile
    << "\n*******************************************************************\n";
  scratchlogfile.close();
  
  PROGRESS.print("Orbit estimation finished.");


  // ====== dump data if requested ======
  if (specified(estorbitsinput.foobsdata))
    {
      PROGRESS << "Dumping observation data to file: " << estorbitsinput.foobsdata << ".";
      PROGRESS.print();

      ofstream outfile;
      outfile.open(estorbitsinput.foobsdata, ios::out);
      bk_assert(outfile,"estorbits: dump data file",__FILE__,__LINE__);
      outfile.setf(ios::fixed);

      const uint cen_lin = (master.currentwindow.linelo+master.currentwindow.linehi)/2;
      const uint cen_pix = (master.currentwindow.pixlo +master.currentwindow.pixhi) /2;

      outfile
	<<   "method                   " << setprecision(1) << (estorbitsinput.method==eo_lsq ? "LSQ" : "GRIDSEARCH")
	<< "\ndBpar1                   " << setprecision(5) << Bpar(1,0)
	<< "\ndBperp0                  " << setprecision(4) << Bperp(0,0)
	<< "\nBperp                    " << setprecision(1) << baseline.get_bperp(cen_lin,cen_pix,0)
	<< "\nBtemp                    " << setprecision(0) << Btemp(master.utc1,slave.utc1)
	<< "\nwavelength_master        " << setprecision(7) << master.wavelength
	<< "\nwavelength_slave         " << setprecision(7) << slave.wavelength
	<< "\nrange_bandwidth_master   " << setprecision(0) << master.rbw
	<< "\nrange_bandwidth_slave    " << setprecision(0) << slave.rbw
	<< "\nmean_range               " << setprecision(0) << SOL*master.pix2tr(cen_pix)
	<< "\nmean_incidence_angle     " << setprecision(6) << rad2deg(baseline.get_theta_inc(cen_lin,cen_pix,0))
	<< "\nconv_Bpar1               " << setprecision(8) << rateConversion
	<< "\ndBpar1_onefringe         " << setprecision(6) << dBparRate_one_fringe
	<< "\ndBperp_onefringe         " << setprecision(4) << dBperp_one_fringe
	<< "\nnumber_of_observations   " << setprecision(0) << obs.nUsed()
	<< "\n*******************************************************************"
	<< "****************************************"
	<< "\n__line _pixel ___tAzi__ _____ahm_____ _____avm_____ _____ahs_____"
	<< " _____avs_____ ____weight___ ____phase____";

      for (uint i=0; i<obs.nInit(); i++)
	{
	  if (obs.isUsed(i))
	    {
	      matrix<real8> Ai = obs.lineOfDesignMatrix(i,observationdata::DESIGN_DUMP);
	      outfile 
		<< endl << resetiosflags(ios_base::scientific) << setiosflags(ios_base::fixed) 
		<< setw( 6) << setprecision(0) << obs.getLine(i)+1 << " "
		<< setw( 6) << setprecision(0) << obs.getPixel(i)+1 << " "
		<< setw( 9) << setprecision(6) << Ai(0,0) << " " 
		<< resetiosflags(ios_base::fixed) << setiosflags(ios_base::scientific)
		<< setw(13) << setprecision(6) << Ai(0,1) << " "
		<< setw(13) << setprecision(6) << Ai(0,2) << " "
		<< setw(13) << setprecision(6) << Ai(0,3) << " "
		<< setw(13) << setprecision(6) << Ai(0,4) << " "
		<< setw(13) << setprecision(6) << obs.getWeight(i) << " "
		<< setw(13) << setprecision(6) << obs(i);
	    }
	}
       outfile.close();
    }

  // ______ tidy up ______
  if (specified(estorbitsinput.reforbitfile))
    delete reference, referenceorbit;
} // END estorbits



/****************************************************************
 *    modifyOrbits                                              *
 *                                                              *
 * Modify orbits according to baseline correction coefficients. *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void modifyOrbits(const input_morbits &morbitsinput,
		  const input_gen     &generalinput,
		  const int16         FILEID,
		  const input_ell     &ellips,
		  const slcimage      &image,
                  const slcimage      &masterimage,
                        orbit         &thisorbit,
                        orbit         &masterorbit)
{
  TRACE_FUNCTION("modifyOrbits (HB 17-Mar-2011)");

  // ______ some settings ______
  const int precision_time = 8;    // decimal places in the resultfile
  const int precision_coo  = 4;    // decimal places in the resultfile
  const int precision_vel  = 6;    // probably too many digits ???

  // ______ choose reference ______
  const slcimage *referenceimage;
  orbit *referenceorbit;
  if (specified(morbitsinput.reforbitfile))
    {
      slcimage *tempptr = new slcimage();
      (*tempptr).fillslcimage(morbitsinput.reforbitfile);
      referenceimage = tempptr;
      referenceorbit = new orbit();
      (*referenceorbit).set_interp_method(generalinput.orb_interp);
      (*referenceorbit).initialize(morbitsinput.reforbitfile);
      INFO << "Using reference orbit from file: " << morbitsinput.reforbitfile;
      INFO.print();
    }
  else
    {
      INFO.print("Using master orbit as reference orbit.");
      referenceimage = &masterimage;
      referenceorbit = &masterorbit;
    }

  // ______ apply modification to orbit and dump to outfile ______
  real8 azoffset = getAzOffset(ellips,*referenceimage,image,*referenceorbit,thisorbit);
  const real8 tshift = -image.t_azi1 + referenceimage->t_azi1 - azoffset;
  const real8 tmin = referenceimage->t_azi1;   
  const real8 tmax = referenceimage->t_azi1
    + (referenceimage->originalwindow.lines()-1)/referenceimage->prf;
  matrix<real8> datapoints = thisorbit.modify(morbitsinput.coeff,*referenceorbit,
					      tshift,tmin,tmax);
 
  // ______ output to file ______
  ofstream outfile("scratchdatapoints", ios::out | ios::trunc);	// do replace
  bk_assert(outfile,"modifyOrbits: scratchdatapoints",__FILE__,__LINE__);
  outfile 
    << "\n\n*******************************************************************"
    << "\n*_Start_";
  if (FILEID==MASTERID) outfile << processcontrol[pr_m_morbits];
  else if (FILEID==SLAVEID) outfile << processcontrol[pr_s_morbits];
  outfile
    << "\n*******************************************************************"
    << "\nt(s)\t\tX(m)\t\tY(m)\t\tZ(m)";
  if (datapoints.pixels()==7) outfile << "\t\tX_V(m/s)\t\tY_V(m/s)\t\tZ_V(m/s)";
  outfile << "\nNUMBER_OF_DATAPOINTS: \t\t\t" << thisorbit.npoints() << endl;
  outfile.setf(ios_base::fixed);
  outfile.setf(ios_base::showpoint);
  for (uint i=0; i<thisorbit.npoints(); i++)
    {
      outfile
	<< setw(6+precision_time) << setprecision(precision_time) << datapoints(i,0) << "\t" 
	<< setw(8+precision_coo) << setprecision(precision_coo)
	<< datapoints(i,1) << "\t" << datapoints(i,2) << "\t" << datapoints(i,3);
      if (datapoints.pixels()==7)
	outfile
	  << setw(4+precision_vel) << setprecision(precision_coo)
	  << datapoints(i,4) << "\t" << datapoints(i,5) << "\t" << datapoints(i,6);
      outfile << endl;
    }
  outfile
    << "\n*******************************************************************"
    << "\n* End_";
  if (FILEID==MASTERID) outfile << processcontrol[pr_m_morbits];
  else if (FILEID==SLAVEID) outfile << processcontrol[pr_s_morbits];
  outfile
    << "\n*******************************************************************\n";
  outfile.close();

  // ______ tidy up ______
  if (specified(morbitsinput.reforbitfile))
    delete referenceimage, referenceorbit;
} // END modifyOrbits



/****************************************************************
 *    disableOldOrbits                                          *
 *                                                              *
 * Initial strategy to account for multiple sections with orbit *
 * data in result file: uncomment old section.              .   *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void disableOldOrbits(const char resultfile[])
{
  TRACE_FUNCTION("disableOldOrbits (HB 17-Mar-2011)");

  // ______Copy relevant part to tmp file______
  // ______ then rename file to old name ______
  //ifstream ifile(file, ios::in | ios::nocreate);
  ifstream ifile(resultfile, ios::in);
  bk_assert(ifile,resultfile,__FILE__,__LINE__);
  ofstream tmpfile("scratchtmp", ios::out | ios::trunc);
  bk_assert(tmpfile,"disableOldOrbits: scratchtmp",__FILE__,__LINE__);
  
  char dummyline[ONE27];
  ifile.getline(dummyline,ONE27,'\n');         
  while (ifile && strncmp(dummyline,"NUMBER_OF_DATAPOINTS:",21)) 
    {
      tmpfile << dummyline << endl;
      ifile.getline(dummyline,ONE27,'\n');
    }
  if (ifile)
    { 
    tmpfile << "#" << dummyline << endl;
    INFO << "Disabled " << setw(2) << "one orbit dataset in " << resultfile << ".";
    INFO.print();
    }
  while (ifile)
    {
      ifile.getline(dummyline,ONE27,'\n');
      tmpfile << dummyline << endl;
    }

  ifile.close();
  tmpfile.close();

  if (rename("scratchtmp",resultfile))                        // rename file
    {
      ERROR << "Could not rename file: scratchtmp to: " << resultfile;
    PRINT_ERROR(ERROR.get_str())
      throw(file_error);
    }
} // END disableOldOrbits



/****************************************************************
 *    getCoordinateFrame                                        *
 *                                                              *
 * Returns Frenet frame of three unit vectors at a specified    *
 * position of an orbit                                         *
 *                                                              *
 * ea: unit vector in along-track direction                     *
 * er: unit vector in radial (vertical) direction               *
 * ex: unit vector in across-track (horizontal) direction       *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011
 ****************************************************************/
void getCoordinateFrame(cn &ea,
			cn &er,
			cn &ex, 
			const real8 tAzi,
			orbit &masterorbit)
{
  TRACE_FUNCTION("getCoordinateFrame (HB 17-Mar-2011)");
 
  ea = masterorbit.getxyzdot(tAzi).normalize();
  er = masterorbit.getxyz(tAzi).normalize();
  er = (er - ea * ea.in(er)).normalize();
  ex = (ea.out(er)).normalize();
} // END getCoordinateFrame



/****************************************************************
 *    repairMatrix                                              *
 *                                                              *
 * Copies the lower triangular half of a quadratic matrix to    *
 * its upper triangular half.                                   *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void repairMatrix(matrix<real8> &M)
{
  TRACE_FUNCTION("repairMatrix (HB 17-Mar-2011)");
 
  for (uint i=0; i<M.lines(); i++)
    for (uint j=i+1; j<M.pixels(); j++)
      M(i,j) = M(j,i);
} // END repairMatrix



/****************************************************************
 *    sigmaTable::sigmaTable                                    *
 *                                                              *
 * Fills lookup-table for phase standard deviations. All        *
 * entries are hard-coded, which seemed to be the most          *
 * efficient way of implementation (?), since their computation *
 * takes a lot of time.                                         *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
  ****************************************************************/
sigmaTable::sigmaTable(uint multilook)
{
  TRACE_FUNCTION("sigmaTable::sigmaTable (HB 17-Mar-2011)");
  if (multilook<1)
    {
      PRINT_ERROR("MULTILOOK FACTOR < 1.");
      throw(some_error);// exit
    }
  
  const uint N_MULTI = 32;
  uint mtable[N_MULTI] = {1,2,3,4,5,6,8,10,13,16,20,25,32,40,50,63,79,100,126,158,200,251,316,398,501,631,794,1000,1259,1585,1995,2512};
  real4 data[101][N_MULTI] = {
    {1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111,1.8111},
    {1.8024,1.7981,1.7949,1.7922,1.7898,1.7876,1.7839,1.7806,1.7762,1.7722,1.7676,1.7623,1.7558,1.7492,1.7418,1.7332,1.7238,1.7127,1.7005,1.6714,1.6714,1.6543,1.6349,1.6131,1.5886,1.5610,1.5302,1.4955,1.4567,1.4134,1.3652,1.3115},
    {1.7938,1.7851,1.7786,1.7731,1.7684,1.7641,1.7565,1.7498,1.7409,1.7331,1.7236,1.7131,1.6999,1.6866,1.6716,1.6542,1.6351,1.6127,1.5881,1.5292,1.5292,1.4950,1.4561,1.4126,1.3643,1.3104,1.2511,1.1856,1.1143,1.0373,0.9552,0.8689},
    {1.7851,1.7721,1.7623,1.7541,1.7469,1.7404,1.7289,1.7189,1.7055,1.6936,1.6793,1.6634,1.6436,1.6234,1.6008,1.5746,1.5458,1.5121,1.4751,1.3876,1.3876,1.3371,1.2804,1.2177,1.1493,1.0748,0.9951,0.9105,0.8228,0.7338,0.6461,0.5625},
    {1.7764,1.7590,1.7459,1.7349,1.7253,1.7166,1.7012,1.6877,1.6698,1.6539,1.6348,1.6135,1.5869,1.5599,1.5297,1.4947,1.4565,1.4119,1.3631,1.2491,1.2491,1.1845,1.1130,1.0357,0.9534,0.8670,0.7785,0.6896,0.6035,0.5232,0.4513,0.3894},
    {1.7677,1.7459,1.7294,1.7156,1.7036,1.6927,1.6734,1.6565,1.6340,1.6140,1.5901,1.5633,1.5301,1.4963,1.4587,1.4152,1.3678,1.3130,1.2535,1.1166,1.1166,1.0407,0.9587,0.8724,0.7838,0.6948,0.6087,0.5278,0.4553,0.3927,0.3404,0.2971},
    {1.7590,1.7327,1.7129,1.6963,1.6818,1.6687,1.6455,1.6251,1.5980,1.5740,1.5452,1.5131,1.4732,1.4329,1.3881,1.3365,1.2806,1.2164,1.1474,0.9922,0.9922,0.9087,0.8209,0.7316,0.6439,0.5603,0.4843,0.4173,0.3608,0.3140,0.2752,0.2425},
    {1.7503,1.7195,1.6963,1.6769,1.6599,1.6446,1.6175,1.5936,1.5620,1.5339,1.5003,1.4629,1.4166,1.3699,1.3182,1.2591,1.1954,1.1230,1.0462,0.8779,0.8779,0.7905,0.7015,0.6147,0.5335,0.4601,0.3969,0.3437,0.2999,0.2633,0.2324,0.2057},
    {1.7415,1.7063,1.6797,1.6575,1.6380,1.6204,1.5893,1.5620,1.5258,1.4938,1.4555,1.4129,1.3604,1.3076,1.2495,1.1834,1.1129,1.0337,0.9508,0.7748,0.7748,0.6872,0.6012,0.5208,0.4491,0.3874,0.3360,0.2934,0.2579,0.2277,0.2016,0.1788},
    {1.7327,1.6930,1.6630,1.6380,1.6160,1.5962,1.5612,1.5304,1.4897,1.4537,1.4108,1.3632,1.3047,1.2462,1.1821,1.1099,1.0336,0.9490,0.8620,0.6837,0.6837,0.5989,0.5189,0.4472,0.3859,0.3346,0.2923,0.2569,0.2269,0.2009,0.1782,0.1582},
    {1.7239,1.6796,1.6463,1.6184,1.5939,1.5719,1.5329,1.4988,1.4536,1.4137,1.3663,1.3139,1.2497,1.1859,1.1165,1.0389,0.9580,0.8695,0.7803,0.6045,0.6045,0.5250,0.4527,0.3903,0.3383,0.2953,0.2595,0.2290,0.2027,0.1798,0.1596,0.1419},
    {1.7151,1.6663,1.6295,1.5987,1.5718,1.5475,1.5047,1.4671,1.4176,1.3739,1.3221,1.2651,1.1956,1.1269,1.0528,0.9708,0.8864,0.7957,0.7062,0.5369,0.5369,0.4641,0.4000,0.3462,0.3019,0.2649,0.2338,0.2068,0.1833,0.1628,0.1446,0.1286},
    {1.7062,1.6529,1.6126,1.5790,1.5496,1.5231,1.4764,1.4355,1.3817,1.3343,1.2783,1.2169,1.1425,1.0694,0.9913,0.9059,0.8193,0.7278,0.6396,0.4798,0.4798,0.4143,0.3582,0.3116,0.2731,0.2406,0.2128,0.1886,0.1674,0.1487,0.1322,0.1176},
    {1.6973,1.6394,1.5957,1.5593,1.5274,1.4987,1.4481,1.4039,1.3459,1.2949,1.2350,1.1694,1.0905,1.0137,0.9323,0.8445,0.7567,0.6659,0.5805,0.4321,0.4321,0.3738,0.3245,0.2838,0.2497,0.2206,0.1954,0.1733,0.1539,0.1368,0.1217,0.1083},
    {1.6884,1.6259,1.5788,1.5395,1.5052,1.4743,1.4199,1.3724,1.3103,1.2559,1.1921,1.1228,1.0399,0.9598,0.8760,0.7866,0.6989,0.6100,0.5286,0.3924,0.3924,0.3406,0.2971,0.2609,0.2302,0.2037,0.1806,0.1603,0.1425,0.1267,0.1127,0.1003},
    {1.6794,1.6124,1.5618,1.5197,1.4829,1.4498,1.3916,1.3410,1.2749,1.2173,1.1499,1.0771,0.9906,0.9080,0.8224,0.7325,0.6458,0.5600,0.4832,0.3593,0.3593,0.3131,0.2743,0.2416,0.2136,0.1892,0.1679,0.1492,0.1326,0.1179,0.1050,0.0934},
    {1.6704,1.5988,1.5448,1.4998,1.4606,1.4253,1.3635,1.3098,1.2398,1.1790,1.1084,1.0324,0.9429,0.8583,0.7716,0.6821,0.5974,0.5154,0.4438,0.3315,0.3315,0.2901,0.2549,0.2250,0.1992,0.1767,0.1569,0.1394,0.1240,0.1103,0.0982,0.0874},
    {1.6614,1.5851,1.5277,1.4799,1.4382,1.4009,1.3354,1.2787,1.2050,1.1413,1.0676,0.9888,0.8969,0.8108,0.7239,0.6356,0.5536,0.4760,0.4097,0.3079,0.3079,0.2704,0.2383,0.2107,0.1867,0.1656,0.1472,0.1308,0.1164,0.1036,0.0922,0.0821},
    {1.6524,1.5714,1.5106,1.4600,1.4159,1.3764,1.3074,1.2477,1.1706,1.1041,1.0276,0.9464,0.8525,0.7657,0.6791,0.5927,0.5141,0.4412,0.3802,0.2877,0.2877,0.2534,0.2237,0.1980,0.1756,0.1559,0.1386,0.1232,0.1096,0.0976,0.0869,0.0773},
    {1.6433,1.5577,1.4934,1.4400,1.3935,1.3520,1.2795,1.2170,1.1365,1.0675,0.9885,0.9053,0.8099,0.7229,0.6374,0.5535,0.4787,0.4107,0.3546,0.2702,0.2702,0.2384,0.2108,0.1868,0.1658,0.1472,0.1309,0.1164,0.1036,0.0922,0.0821,0.0731},
    {1.6342,1.5439,1.4762,1.4200,1.3712,1.3276,1.2517,1.1865,1.1028,1.0315,0.9503,0.8654,0.7692,0.6825,0.5986,0.5177,0.4470,0.3838,0.3324,0.2547,0.2547,0.2252,0.1993,0.1767,0.1569,0.1394,0.1240,0.1103,0.0981,0.0874,0.0778,0.0693},
    {1.6250,1.5301,1.4590,1.4000,1.3489,1.3033,1.2241,1.1563,1.0696,0.9961,0.9130,0.8269,0.7304,0.6445,0.5627,0.4853,0.4187,0.3602,0.3129,0.2410,0.2410,0.2134,0.1890,0.1677,0.1489,0.1323,0.1177,0.1047,0.0932,0.0830,0.0739,0.0658},
    {1.6158,1.5163,1.4417,1.3800,1.3266,1.2790,1.1966,1.1263,1.0369,0.9615,0.8768,0.7898,0.6934,0.6089,0.5296,0.4559,0.3935,0.3393,0.2956,0.2287,0.2287,0.2027,0.1797,0.1594,0.1417,0.1259,0.1121,0.0997,0.0888,0.0790,0.0704,0.0627},
    {1.6066,1.5023,1.4244,1.3600,1.3043,1.2548,1.1693,1.0967,1.0047,0.9276,0.8416,0.7542,0.6584,0.5756,0.4992,0.4293,0.3710,0.3208,0.2803,0.2176,0.2176,0.1930,0.1712,0.1520,0.1351,0.1201,0.1069,0.0951,0.0847,0.0754,0.0672,0.0598},
    {1.5973,1.4884,1.4070,1.3399,1.2820,1.2307,1.1422,1.0673,0.9730,0.8945,0.8076,0.7200,0.6253,0.5446,0.4713,0.4053,0.3509,0.3043,0.2665,0.2075,0.2075,0.1841,0.1634,0.1451,0.1290,0.1147,0.1021,0.0909,0.0809,0.0721,0.0642,0.0572},
    {1.5880,1.4743,1.3896,1.3199,1.2598,1.2066,1.1152,1.0383,0.9419,0.8622,0.7747,0.6873,0.5940,0.5158,0.4458,0.3836,0.3329,0.2894,0.2541,0.1982,0.1982,0.1760,0.1563,0.1388,0.1234,0.1098,0.0977,0.0870,0.0774,0.0690,0.0614,0.0547},
    {1.5786,1.4603,1.3722,1.2998,1.2376,1.1827,1.0886,1.0097,0.9115,0.8307,0.7429,0.6561,0.5647,0.4891,0.4225,0.3640,0.3167,0.2760,0.2427,0.1897,0.1897,0.1686,0.1497,0.1330,0.1183,0.1052,0.0937,0.0834,0.0742,0.0661,0.0589,0.0525},
    {1.5692,1.4462,1.3547,1.2798,1.2155,1.1588,1.0621,0.9815,0.8816,0.8002,0.7123,0.6264,0.5371,0.4644,0.4012,0.3463,0.3020,0.2638,0.2323,0.1819,0.1819,0.1617,0.1436,0.1276,0.1135,0.1010,0.0899,0.0800,0.0713,0.0635,0.0566,0.0504},
    {1.5597,1.4320,1.3373,1.2597,1.1934,1.1351,1.0359,0.9536,0.8524,0.7705,0.6828,0.5982,0.5113,0.4416,0.3818,0.3302,0.2886,0.2526,0.2228,0.1746,0.1746,0.1553,0.1379,0.1226,0.1091,0.0971,0.0864,0.0769,0.0685,0.0610,0.0544,0.0484},
    {1.5502,1.4178,1.3197,1.2397,1.1714,1.1115,1.0100,0.9262,0.8239,0.7417,0.6546,0.5714,0.4872,0.4206,0.3640,0.3156,0.2764,0.2424,0.2139,0.1679,0.1679,0.1493,0.1327,0.1179,0.1049,0.0934,0.0832,0.0740,0.0659,0.0587,0.0523,0.0466},
    {1.5406,1.4035,1.3022,1.2197,1.1494,1.0880,0.9844,0.8993,0.7960,0.7138,0.6275,0.5461,0.4647,0.4012,0.3478,0.3021,0.2652,0.2329,0.2057,0.1616,0.1616,0.1437,0.1277,0.1136,0.1011,0.0899,0.0801,0.0713,0.0635,0.0566,0.0504,0.0449},
    {1.5310,1.3892,1.2846,1.1997,1.1276,1.0647,0.9590,0.8728,0.7689,0.6869,0.6017,0.5221,0.4438,0.3833,0.3329,0.2898,0.2548,0.2240,0.1981,0.1557,0.1557,0.1385,0.1231,0.1095,0.0975,0.0867,0.0772,0.0688,0.0613,0.0546,0.0486,0.0433},
    {1.5213,1.3748,1.2670,1.1798,1.1058,1.0415,0.9340,0.8468,0.7425,0.6609,0.5770,0.4996,0.4243,0.3668,0.3191,0.2785,0.2452,0.2158,0.1909,0.1501,0.1501,0.1336,0.1188,0.1057,0.0940,0.0837,0.0746,0.0664,0.0591,0.0527,0.0469,0.0418},
    {1.5116,1.3603,1.2494,1.1598,1.0841,1.0186,0.9093,0.8213,0.7168,0.6359,0.5534,0.4783,0.4061,0.3516,0.3065,0.2679,0.2363,0.2081,0.1842,0.1450,0.1450,0.1290,0.1147,0.1021,0.0908,0.0809,0.0720,0.0641,0.0571,0.0509,0.0453,0.0404},
    {1.5018,1.3458,1.2318,1.1399,1.0626,0.9957,0.8850,0.7963,0.6919,0.6118,0.5311,0.4583,0.3892,0.3375,0.2948,0.2581,0.2279,0.2009,0.1779,0.1401,0.1401,0.1247,0.1109,0.0986,0.0878,0.0782,0.0696,0.0620,0.0552,0.0492,0.0438,0.0391},
    {1.4919,1.3313,1.2141,1.1201,1.0411,0.9731,0.8610,0.7719,0.6678,0.5887,0.5098,0.4395,0.3735,0.3244,0.2839,0.2490,0.2200,0.1941,0.1719,0.1354,0.1354,0.1206,0.1073,0.0954,0.0849,0.0756,0.0674,0.0600,0.0534,0.0476,0.0424,0.0378},
    {1.4820,1.3167,1.1965,1.1002,1.0197,0.9507,0.8374,0.7479,0.6444,0.5665,0.4896,0.4219,0.3589,0.3122,0.2737,0.2404,0.2126,0.1876,0.1663,0.1311,0.1311,0.1167,0.1038,0.0924,0.0822,0.0732,0.0652,0.0581,0.0517,0.0461,0.0411,0.0366},
    {1.4720,1.3020,1.1788,1.0805,0.9985,0.9285,0.8142,0.7245,0.6218,0.5452,0.4704,0.4053,0.3452,0.3009,0.2643,0.2324,0.2057,0.1816,0.1610,0.1269,0.1269,0.1130,0.1006,0.0895,0.0797,0.0709,0.0632,0.0563,0.0501,0.0447,0.0398,0.0355},
    {1.4619,1.2873,1.1611,1.0607,0.9774,0.9065,0.7913,0.7017,0.5999,0.5249,0.4523,0.3898,0.3325,0.2903,0.2554,0.2248,0.1991,0.1758,0.1559,0.1230,0.1230,0.1095,0.0974,0.0867,0.0772,0.0687,0.0612,0.0545,0.0486,0.0433,0.0386,0.0344},
    {1.4518,1.2725,1.1433,1.0411,0.9565,0.8847,0.7689,0.6795,0.5788,0.5054,0.4352,0.3752,0.3206,0.2804,0.2470,0.2176,0.1928,0.1704,0.1511,0.1192,0.1192,0.1062,0.0945,0.0841,0.0749,0.0667,0.0594,0.0529,0.0471,0.0420,0.0374,0.0333},
    {1.4416,1.2576,1.1256,1.0215,0.9356,0.8631,0.7468,0.6578,0.5585,0.4868,0.4189,0.3615,0.3094,0.2711,0.2391,0.2108,0.1869,0.1652,0.1465,0.1156,0.1156,0.1030,0.0917,0.0816,0.0726,0.0647,0.0576,0.0513,0.0457,0.0407,0.0363,0.0323},
    {1.4313,1.2427,1.1079,1.0019,0.9150,0.8418,0.7252,0.6367,0.5389,0.4691,0.4036,0.3486,0.2990,0.2623,0.2316,0.2043,0.1812,0.1602,0.1422,0.1122,0.1122,0.1000,0.0890,0.0792,0.0705,0.0628,0.0559,0.0498,0.0444,0.0396,0.0352,0.0314},
    {1.4209,1.2277,1.0901,0.9825,0.8945,0.8208,0.7041,0.6162,0.5201,0.4522,0.3891,0.3365,0.2891,0.2540,0.2245,0.1982,0.1759,0.1555,0.1380,0.1090,0.1090,0.0971,0.0864,0.0769,0.0685,0.0610,0.0543,0.0484,0.0431,0.0384,0.0342,0.0305},
    {1.4104,1.2126,1.0724,0.9631,0.8741,0.8000,0.6833,0.5963,0.5020,0.4361,0.3754,0.3250,0.2798,0.2462,0.2177,0.1923,0.1707,0.1510,0.1340,0.1058,0.1058,0.0943,0.0839,0.0747,0.0665,0.0593,0.0528,0.0470,0.0419,0.0373,0.0333,0.0296},
    {1.3999,1.1975,1.0546,0.9437,0.8540,0.7794,0.6630,0.5769,0.4846,0.4208,0.3624,0.3143,0.2710,0.2387,0.2113,0.1867,0.1658,0.1467,0.1302,0.1029,0.1029,0.0917,0.0816,0.0726,0.0647,0.0576,0.0513,0.0457,0.0407,0.0363,0.0323,0.0288},
    {1.3893,1.1823,1.0368,0.9245,0.8340,0.7592,0.6432,0.5582,0.4680,0.4062,0.3502,0.3041,0.2626,0.2316,0.2051,0.1814,0.1611,0.1425,0.1265,0.1000,0.1000,0.0891,0.0793,0.0706,0.0629,0.0560,0.0499,0.0444,0.0396,0.0353,0.0314,0.0280},
    {1.3786,1.1671,1.0191,0.9053,0.8142,0.7392,0.6238,0.5400,0.4520,0.3923,0.3386,0.2945,0.2547,0.2248,0.1992,0.1762,0.1565,0.1385,0.1230,0.0972,0.0972,0.0867,0.0771,0.0687,0.0612,0.0545,0.0485,0.0432,0.0385,0.0343,0.0306,0.0272},
    {1.3677,1.1518,1.0013,0.8863,0.7945,0.7195,0.6048,0.5224,0.4367,0.3791,0.3275,0.2853,0.2472,0.2184,0.1936,0.1713,0.1522,0.1347,0.1197,0.0946,0.0946,0.0843,0.0750,0.0668,0.0595,0.0530,0.0472,0.0421,0.0375,0.0334,0.0298,0.0265},
    {1.3568,1.1364,0.9835,0.8673,0.7751,0.7001,0.5864,0.5054,0.4220,0.3666,0.3171,0.2766,0.2400,0.2122,0.1882,0.1665,0.1480,0.1310,0.1164,0.0920,0.0920,0.0820,0.0730,0.0650,0.0579,0.0516,0.0459,0.0409,0.0365,0.0325,0.0290,0.0258},
    {1.3458,1.1209,0.9658,0.8485,0.7559,0.6811,0.5683,0.4889,0.4080,0.3546,0.3072,0.2684,0.2331,0.2062,0.1830,0.1620,0.1440,0.1275,0.1133,0.0895,0.0895,0.0798,0.0711,0.0633,0.0563,0.0502,0.0447,0.0398,0.0355,0.0316,0.0282,0.0251},
    {1.3347,1.1053,0.9480,0.8297,0.7369,0.6623,0.5508,0.4730,0.3946,0.3432,0.2977,0.2605,0.2265,0.2005,0.1779,0.1576,0.1401,0.1241,0.1102,0.0872,0.0872,0.0777,0.0692,0.0616,0.0549,0.0489,0.0435,0.0388,0.0346,0.0308,0.0274,0.0244},
    {1.3234,1.0897,0.9303,0.8111,0.7181,0.6438,0.5337,0.4577,0.3817,0.3323,0.2887,0.2529,0.2201,0.1950,0.1731,0.1533,0.1363,0.1208,0.1073,0.0849,0.0849,0.0756,0.0674,0.0600,0.0534,0.0476,0.0424,0.0378,0.0336,0.0300,0.0267,0.0238},
    {1.3121,1.0740,0.9126,0.7926,0.6996,0.6257,0.5171,0.4429,0.3695,0.3220,0.2801,0.2457,0.2140,0.1897,0.1685,0.1492,0.1327,0.1176,0.1045,0.0826,0.0826,0.0737,0.0656,0.0584,0.0520,0.0463,0.0413,0.0368,0.0328,0.0292,0.0260,0.0232},
    {1.3006,1.0582,0.8948,0.7742,0.6812,0.6079,0.5010,0.4287,0.3577,0.3121,0.2719,0.2388,0.2082,0.1845,0.1639,0.1453,0.1292,0.1145,0.1017,0.0805,0.0805,0.0717,0.0639,0.0569,0.0507,0.0451,0.0402,0.0358,0.0319,0.0284,0.0253,0.0226},
    {1.2890,1.0424,0.8771,0.7559,0.6631,0.5904,0.4853,0.4149,0.3464,0.3026,0.2640,0.2321,0.2025,0.1796,0.1596,0.1414,0.1258,0.1115,0.0991,0.0784,0.0784,0.0699,0.0622,0.0554,0.0493,0.0440,0.0392,0.0349,0.0311,0.0277,0.0247,0.0220},
    {1.2773,1.0264,0.8594,0.7378,0.6453,0.5732,0.4700,0.4017,0.3356,0.2936,0.2565,0.2257,0.1970,0.1748,0.1554,0.1377,0.1225,0.1086,0.0965,0.0763,0.0763,0.0681,0.0606,0.0540,0.0481,0.0428,0.0382,0.0340,0.0303,0.0270,0.0241,0.0214},
    {1.2654,1.0104,0.8418,0.7197,0.6277,0.5564,0.4553,0.3889,0.3253,0.2849,0.2492,0.2195,0.1917,0.1701,0.1513,0.1341,0.1193,0.1057,0.0940,0.0744,0.0744,0.0663,0.0590,0.0526,0.0468,0.0417,0.0372,0.0331,0.0295,0.0263,0.0234,0.0209},
    {1.2534,0.9942,0.8241,0.7019,0.6103,0.5399,0.4409,0.3766,0.3154,0.2766,0.2422,0.2135,0.1866,0.1656,0.1473,0.1306,0.1162,0.1030,0.0916,0.0724,0.0724,0.0646,0.0575,0.0512,0.0456,0.0406,0.0362,0.0323,0.0287,0.0256,0.0228,0.0203},
    {1.2413,0.9780,0.8065,0.6842,0.5932,0.5238,0.4270,0.3648,0.3059,0.2686,0.2355,0.2077,0.1816,0.1612,0.1434,0.1272,0.1132,0.1003,0.0892,0.0706,0.0706,0.0629,0.0560,0.0499,0.0445,0.0396,0.0353,0.0314,0.0280,0.0250,0.0222,0.0198},
    {1.2290,0.9617,0.7889,0.6666,0.5763,0.5080,0.4135,0.3534,0.2967,0.2609,0.2290,0.2021,0.1768,0.1570,0.1397,0.1239,0.1103,0.0977,0.0869,0.0688,0.0688,0.0613,0.0546,0.0486,0.0433,0.0386,0.0344,0.0306,0.0273,0.0243,0.0217,0.0193},
    {1.2165,0.9453,0.7713,0.6492,0.5597,0.4925,0.4005,0.3424,0.2879,0.2534,0.2226,0.1966,0.1721,0.1529,0.1360,0.1206,0.1074,0.0952,0.0846,0.0670,0.0670,0.0597,0.0532,0.0474,0.0422,0.0376,0.0335,0.0298,0.0266,0.0237,0.0211,0.0188},
    {1.2039,0.9288,0.7537,0.6319,0.5434,0.4774,0.3878,0.3318,0.2794,0.2462,0.2165,0.1913,0.1675,0.1488,0.1324,0.1175,0.1046,0.0927,0.0824,0.0653,0.0653,0.0582,0.0518,0.0461,0.0411,0.0366,0.0326,0.0291,0.0259,0.0231,0.0206,0.0183},
    {1.1910,0.9122,0.7362,0.6148,0.5274,0.4626,0.3756,0.3215,0.2712,0.2393,0.2106,0.1862,0.1630,0.1449,0.1289,0.1144,0.1019,0.0903,0.0803,0.0636,0.0636,0.0567,0.0505,0.0450,0.0401,0.0357,0.0318,0.0283,0.0252,0.0225,0.0200,0.0179},
    {1.1781,0.8955,0.7187,0.5978,0.5116,0.4482,0.3637,0.3117,0.2633,0.2326,0.2048,0.1812,0.1587,0.1410,0.1255,0.1114,0.0992,0.0879,0.0782,0.0619,0.0619,0.0552,0.0492,0.0438,0.0390,0.0347,0.0310,0.0276,0.0246,0.0219,0.0195,0.0174},
    {1.1649,0.8786,0.7012,0.5811,0.4960,0.4341,0.3522,0.3021,0.2557,0.2260,0.1992,0.1763,0.1544,0.1373,0.1222,0.1085,0.0966,0.0856,0.0762,0.0603,0.0603,0.0538,0.0479,0.0426,0.0380,0.0338,0.0302,0.0269,0.0239,0.0213,0.0190,0.0169},
    {1.1515,0.8617,0.6838,0.5645,0.4808,0.4203,0.3411,0.2929,0.2482,0.2197,0.1937,0.1715,0.1503,0.1336,0.1190,0.1056,0.0940,0.0834,0.0741,0.0587,0.0587,0.0524,0.0466,0.0415,0.0370,0.0330,0.0294,0.0262,0.0233,0.0208,0.0185,0.0165},
    {1.1379,0.8446,0.6664,0.5481,0.4658,0.4069,0.3303,0.2840,0.2410,0.2135,0.1884,0.1668,0.1462,0.1300,0.1158,0.1028,0.0915,0.0812,0.0722,0.0572,0.0572,0.0510,0.0454,0.0404,0.0360,0.0321,0.0286,0.0255,0.0227,0.0202,0.0180,0.0161},
    {1.1240,0.8275,0.6491,0.5318,0.4511,0.3937,0.3198,0.2753,0.2341,0.2075,0.1831,0.1622,0.1422,0.1265,0.1127,0.1000,0.0891,0.0790,0.0703,0.0556,0.0556,0.0496,0.0442,0.0394,0.0351,0.0312,0.0278,0.0248,0.0221,0.0197,0.0175,0.0156},
    {1.1100,0.8102,0.6318,0.5157,0.4367,0.3810,0.3097,0.2670,0.2272,0.2016,0.1780,0.1577,0.1383,0.1231,0.1096,0.0973,0.0867,0.0769,0.0684,0.0541,0.0541,0.0483,0.0430,0.0383,0.0341,0.0304,0.0271,0.0241,0.0215,0.0192,0.0171,0.0152},
    {1.0957,0.7927,0.6145,0.4999,0.4226,0.3685,0.2998,0.2588,0.2206,0.1958,0.1730,0.1533,0.1345,0.1197,0.1066,0.0946,0.0843,0.0748,0.0665,0.0527,0.0527,0.0470,0.0418,0.0373,0.0332,0.0296,0.0263,0.0235,0.0209,0.0186,0.0166,0.0148},
    {1.0811,0.7751,0.5973,0.4842,0.4087,0.3563,0.2903,0.2509,0.2141,0.1902,0.1681,0.1490,0.1307,0.1163,0.1036,0.0920,0.0820,0.0727,0.0647,0.0512,0.0512,0.0457,0.0407,0.0362,0.0323,0.0288,0.0256,0.0228,0.0203,0.0181,0.0162,0.0144},
    {1.0662,0.7574,0.5801,0.4687,0.3951,0.3444,0.2810,0.2432,0.2078,0.1846,0.1633,0.1448,0.1270,0.1131,0.1007,0.0894,0.0797,0.0707,0.0629,0.0498,0.0498,0.0444,0.0395,0.0352,0.0314,0.0280,0.0249,0.0222,0.0198,0.0176,0.0157,0.0140},
    {1.0511,0.7396,0.5629,0.4534,0.3818,0.3328,0.2719,0.2356,0.2016,0.1792,0.1585,0.1406,0.1234,0.1098,0.0978,0.0869,0.0774,0.0687,0.0611,0.0484,0.0484,0.0431,0.0384,0.0342,0.0305,0.0272,0.0242,0.0216,0.0192,0.0171,0.0153,0.0136},
    {1.0356,0.7215,0.5458,0.4383,0.3687,0.3215,0.2631,0.2283,0.1955,0.1739,0.1539,0.1365,0.1198,0.1066,0.0950,0.0844,0.0752,0.0667,0.0593,0.0470,0.0470,0.0419,0.0373,0.0332,0.0296,0.0264,0.0235,0.0209,0.0187,0.0166,0.0148,0.0132},
    {1.0198,0.7034,0.5288,0.4233,0.3559,0.3105,0.2544,0.2211,0.1895,0.1686,0.1493,0.1324,0.1163,0.1035,0.0922,0.0819,0.0730,0.0647,0.0576,0.0456,0.0456,0.0407,0.0362,0.0323,0.0288,0.0256,0.0228,0.0203,0.0181,0.0162,0.0144,0.0128},
    {1.0037,0.6850,0.5118,0.4086,0.3433,0.2997,0.2460,0.2140,0.1836,0.1634,0.1447,0.1284,0.1128,0.1004,0.0895,0.0795,0.0708,0.0628,0.0559,0.0443,0.0443,0.0395,0.0352,0.0313,0.0279,0.0249,0.0222,0.0197,0.0176,0.0157,0.0140,0.0124},
    {0.9871,0.6665,0.4949,0.3940,0.3310,0.2891,0.2378,0.2071,0.1778,0.1583,0.1402,0.1245,0.1093,0.0973,0.0867,0.0770,0.0686,0.0609,0.0542,0.0429,0.0429,0.0383,0.0341,0.0304,0.0271,0.0241,0.0215,0.0191,0.0171,0.0152,0.0135,0.0121},
    {0.9702,0.6477,0.4780,0.3797,0.3189,0.2788,0.2297,0.2003,0.1721,0.1533,0.1358,0.1205,0.1059,0.0943,0.0840,0.0746,0.0665,0.0590,0.0525,0.0416,0.0416,0.0371,0.0330,0.0294,0.0262,0.0234,0.0208,0.0185,0.0165,0.0147,0.0131,0.0117},
    {0.9528,0.6288,0.4611,0.3655,0.3070,0.2686,0.2217,0.1935,0.1664,0.1483,0.1314,0.1166,0.1025,0.0912,0.0813,0.0722,0.0644,0.0571,0.0508,0.0402,0.0402,0.0359,0.0320,0.0285,0.0254,0.0226,0.0201,0.0180,0.0160,0.0143,0.0127,0.0113},
    {0.9350,0.6097,0.4443,0.3515,0.2953,0.2587,0.2139,0.1869,0.1608,0.1433,0.1270,0.1128,0.0991,0.0882,0.0786,0.0699,0.0623,0.0552,0.0491,0.0389,0.0389,0.0347,0.0309,0.0276,0.0245,0.0219,0.0195,0.0174,0.0155,0.0138,0.0123,0.0110},
    {0.9166,0.5903,0.4275,0.3376,0.2838,0.2489,0.2062,0.1803,0.1552,0.1384,0.1227,0.1089,0.0957,0.0852,0.0760,0.0675,0.0602,0.0534,0.0475,0.0376,0.0376,0.0336,0.0299,0.0266,0.0237,0.0211,0.0188,0.0168,0.0150,0.0133,0.0119,0.0106},
    {0.8977,0.5707,0.4108,0.3240,0.2725,0.2393,0.1986,0.1738,0.1497,0.1335,0.1183,0.1051,0.0924,0.0823,0.0733,0.0652,0.0581,0.0515,0.0458,0.0363,0.0363,0.0324,0.0289,0.0257,0.0229,0.0204,0.0182,0.0162,0.0144,0.0129,0.0115,0.0102},
    {0.8782,0.5508,0.3941,0.3104,0.2614,0.2298,0.1910,0.1673,0.1442,0.1286,0.1140,0.1013,0.0890,0.0793,0.0707,0.0628,0.0560,0.0497,0.0442,0.0350,0.0350,0.0312,0.0278,0.0248,0.0221,0.0197,0.0175,0.0156,0.0139,0.0124,0.0111,0.0099},
    {0.8580,0.5306,0.3774,0.2971,0.2504,0.2204,0.1835,0.1609,0.1387,0.1237,0.1097,0.0975,0.0857,0.0763,0.0680,0.0605,0.0539,0.0478,0.0425,0.0337,0.0337,0.0301,0.0268,0.0239,0.0213,0.0189,0.0169,0.0150,0.0134,0.0119,0.0106,0.0095},
    {0.8371,0.5102,0.3607,0.2838,0.2396,0.2112,0.1761,0.1544,0.1332,0.1188,0.1054,0.0937,0.0823,0.0733,0.0654,0.0581,0.0518,0.0460,0.0409,0.0324,0.0324,0.0289,0.0257,0.0229,0.0204,0.0182,0.0162,0.0145,0.0129,0.0115,0.0102,0.0091},
    {0.8154,0.4894,0.3440,0.2707,0.2288,0.2019,0.1686,0.1480,0.1277,0.1140,0.1011,0.0898,0.0790,0.0704,0.0627,0.0557,0.0497,0.0441,0.0392,0.0311,0.0311,0.0277,0.0247,0.0220,0.0196,0.0175,0.0156,0.0139,0.0124,0.0110,0.0098,0.0087},
    {0.7928,0.4682,0.3273,0.2576,0.2181,0.1928,0.1612,0.1415,0.1221,0.1090,0.0967,0.0860,0.0756,0.0673,0.0600,0.0534,0.0476,0.0422,0.0376,0.0298,0.0298,0.0265,0.0236,0.0211,0.0188,0.0167,0.0149,0.0133,0.0118,0.0105,0.0094,0.0084},
    {0.7693,0.4466,0.3106,0.2446,0.2075,0.1836,0.1537,0.1350,0.1166,0.1041,0.0924,0.0821,0.0722,0.0643,0.0573,0.0510,0.0454,0.0403,0.0359,0.0284,0.0284,0.0254,0.0226,0.0201,0.0179,0.0160,0.0142,0.0127,0.0113,0.0101,0.0090,0.0080},
    {0.7446,0.4246,0.2939,0.2316,0.1968,0.1744,0.1462,0.1284,0.1109,0.0991,0.0879,0.0782,0.0687,0.0612,0.0546,0.0485,0.0433,0.0384,0.0342,0.0271,0.0271,0.0241,0.0215,0.0192,0.0171,0.0152,0.0136,0.0121,0.0108,0.0096,0.0085,0.0076},
    {0.7186,0.4020,0.2770,0.2186,0.1861,0.1651,0.1385,0.1218,0.1052,0.0940,0.0834,0.0742,0.0652,0.0581,0.0518,0.0461,0.0411,0.0364,0.0324,0.0257,0.0257,0.0229,0.0204,0.0182,0.0162,0.0144,0.0129,0.0115,0.0102,0.0091,0.0081,0.0072},
    {0.6911,0.3789,0.2600,0.2056,0.1753,0.1557,0.1308,0.1150,0.0994,0.0888,0.0788,0.0701,0.0616,0.0549,0.0490,0.0435,0.0388,0.0344,0.0306,0.0243,0.0243,0.0217,0.0193,0.0172,0.0153,0.0136,0.0122,0.0108,0.0097,0.0086,0.0077,0.0068},
    {0.6620,0.3551,0.2428,0.1924,0.1644,0.1462,0.1229,0.1081,0.0935,0.0835,0.0741,0.0659,0.0580,0.0517,0.0461,0.0410,0.0365,0.0324,0.0288,0.0228,0.0228,0.0204,0.0182,0.0162,0.0144,0.0128,0.0114,0.0102,0.0091,0.0081,0.0072,0.0064},
    {0.6307,0.3305,0.2253,0.1790,0.1533,0.1364,0.1148,0.1010,0.0873,0.0780,0.0693,0.0616,0.0542,0.0483,0.0431,0.0383,0.0341,0.0303,0.0270,0.0214,0.0214,0.0191,0.0170,0.0151,0.0135,0.0120,0.0107,0.0095,0.0085,0.0076,0.0067,0.0060},
    {0.5970,0.3050,0.2075,0.1653,0.1418,0.1263,0.1063,0.0936,0.0810,0.0724,0.0643,0.0571,0.0503,0.0448,0.0400,0.0355,0.0317,0.0281,0.0250,0.0198,0.0198,0.0177,0.0157,0.0140,0.0125,0.0111,0.0099,0.0088,0.0079,0.0070,0.0063,0.0056},
    {0.5603,0.2783,0.1891,0.1512,0.1299,0.1158,0.0976,0.0859,0.0743,0.0664,0.0590,0.0525,0.0461,0.0411,0.0367,0.0326,0.0291,0.0258,0.0230,0.0182,0.0182,0.0162,0.0145,0.0129,0.0115,0.0102,0.0091,0.0081,0.0072,0.0064,0.0057,0.0051},
    {0.5195,0.2501,0.1699,0.1363,0.1173,0.1047,0.0883,0.0777,0.0673,0.0601,0.0534,0.0475,0.0418,0.0372,0.0332,0.0295,0.0263,0.0234,0.0208,0.0165,0.0165,0.0147,0.0131,0.0117,0.0104,0.0093,0.0083,0.0074,0.0066,0.0058,0.0052,0.0046},
    {0.4735,0.2199,0.1497,0.1205,0.1039,0.0927,0.0782,0.0689,0.0597,0.0533,0.0474,0.0421,0.0371,0.0330,0.0295,0.0262,0.0234,0.0207,0.0184,0.0146,0.0146,0.0130,0.0116,0.0104,0.0092,0.0082,0.0073,0.0065,0.0058,0.0052,0.0046,0.0041},
    {0.4198,0.1870,0.1277,0.1032,0.0891,0.0796,0.0672,0.0592,0.0512,0.0458,0.0407,0.0362,0.0318,0.0284,0.0253,0.0225,0.0201,0.0178,0.0159,0.0126,0.0126,0.0112,0.0100,0.0089,0.0079,0.0071,0.0063,0.0056,0.0050,0.0045,0.0040,0.0036},
    {0.3539,0.1495,0.1028,0.0834,0.0720,0.0644,0.0544,0.0479,0.0415,0.0371,0.0330,0.0293,0.0258,0.0230,0.0205,0.0182,0.0163,0.0144,0.0128,0.0102,0.0102,0.0091,0.0081,0.0072,0.0064,0.0057,0.0051,0.0045,0.0041,0.0037,0.0034,0.0032},
    {0.2633,0.1032,0.0716,0.0583,0.0505,0.0451,0.0381,0.0336,0.0291,0.0260,0.0231,0.0206,0.0181,0.0161,0.0144,0.0128,0.0114,0.0101,0.0090,0.0071,0.0071,0.0064,0.0057,0.0051,0.0045,0.0040,0.0036,0.0034,0.0032,0.0030,0.0029,0.0027},
    {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000}
  };
  
  real4 logmultilook = log(real4(multilook));
  real4 logmulti[N_MULTI];
  logmulti[0] = log(real4(mtable[0]));
  uint l1, l2;
  real4 logmtable[N_MULTI];
  for (l2=1; l2<N_MULTI; l2++)
    {
      logmtable[l2] = log(real4(mtable[l2]));
      if (logmultilook<logmtable[l2]) break;
    }
  if (l2==N_MULTI)
     {
       PRINT_ERROR("MULTILOOK FACTOR TOO HIGH.");
       throw(some_error);// exit
     }
  l1 = l2 - 1;
  real4 p1 = (logmtable[l2]-logmultilook)/(logmtable[l2]-logmtable[l1]);
  real4 p2 = 1 - p1;
  for (uint c=0; c<101; c++)
    table[c] = p1*data[c][l1] + p2 * data[c][l2];
} // END sigmaTable::sigmaTable



/****************************************************************
 *    sigmaTable::getSigma                                      *
 *                                                              *
 * Retrieve phase standard deviation from lookup-table.         *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
real4 sigmaTable::getSigma(real4 coh) const
{
  uint c1 = int(100*coh);
  real4 p2 = 100*coh-c1;
  return((1-p2)*table[c1]+p2*table[c1+1]);
} // END sigmaTable::getSigma



/****************************************************************
 *    getEigenTheta                                             *
 *                                                              *
 * Return orientation angle of the error ellipses of the        *
 * estimated baseline error parameters.                         *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
real8 getEigenTheta(const matrix<real8> &Q)
{
  TRACE_FUNCTION("getEigenTheta (HB 17-Mar-2011)");
  int16 index1 = 0;
  int16 index2 = Q.lines()/2;
  real8 trace = Q(index1,index1) + Q(index2,index2);
  real8 det = Q(index1,index1) * Q(index2,index2) - Q(index1,index2) * Q(index2,index1);
  real8 lambda1 = trace/2 + sqrt((trace*trace)/4-det);
  real8 theta = atan(Q(index1,index2)/(Q(index1,index1)-lambda1));

  // ______ if available average theta derived from constant and linear component ______
  if (Q.lines()>=4)
    {
      index1++; index2++;
      trace = Q(index1,index1) + Q(index2,index2);
      det = Q(index1,index1) * Q(index2,index2) - Q(index1,index2) * Q(index2,index1);
      lambda1 = trace/2 + sqrt((trace*trace)/4-det);
      theta = 0.5*(theta + atan(Q(index1,index2)/(Q(index1,index1)-lambda1))+PI/2);
    }
  return theta;
} // END getEigenTheta



/****************************************************************
 *    constrainSolution                                         *
 *                                                              *
 * apply constraints to estimated orbit error parameters and    *
 * their covariance matrix.                                     *
 *                                                              *
 * parameters:                                                  *
 *   - npar x npar covariance matrix                            *
 *   - npar x 1 parameter vector (Bh0,Bh1,...,Bv0,Bv1,...)      *
 *   - constraints-array (see definition of struct              *
 *     input_estorbits in readinput.hh for explanation)         *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void constrainSolution(matrix<real8> &Qxx, matrix<real8> &rhs, const matrix<int16> &constraints)
{
  uint nc = constraints.lines();
  uint npar = rhs.lines();
  uint degree = npar/2-1;
  matrix<real8> R(nc,npar);       // coefficient matrix for constraints
  matrix<real8> Qww(nc,nc);       // cofactor matrix of misclosures
  matrix<real8> Proj(npar,npar);  // projector applied to solution

  real8 theta = getEigenTheta(Qxx);

  // ______ setup coefficient matrix for constraints ______
  R = 0;
  for (uint i=0; i<nc; i++)
    {
      if (constraints(i,0)==0) // BPAR
	{
	  R(i,constraints(i,1)) = sin(theta);
	  R(i,degree+1+constraints(i,1)) = -cos(theta);
	}
      else if (constraints(i,0)==1) // BPERP
	{
	  R(i,constraints(i,1)) = cos(theta);
	  R(i,degree+1+constraints(i,1)) = sin(theta);
	}
    }
  
  // ______ compute projector (constrained_solution = Proj * unconstrained_solution ______
  Qww = matxmatT(R*Qxx,R);
  choles(Qww);
  invertchol(Qww);
  repairMatrix(Qww);
  
  Proj = - Qxx * matTxmat(R,Qww*R);
  for (uint i=0; i<npar; i++)
    Proj(i,i) += 1;
  
  // ______ apply to solution ______
  rhs = Proj * rhs;
  Qxx = Proj * Qxx;
} // END constrainSolution



/****************************************************************
 *    getAzOffset                                               *
 *                                                              *
 * returns azimuth time offset between two acquisitions, which  *
 * is the same as the COARSEORB-shift in azimuth converted to   *
 * time                                                         *
 *                                                              *
 * Note: The offset is positive if                              *
 * - COARSEORB-offset in lines is positive.                     *
 * - slave line number > corresponding master line number.      *
 * - slave acquisition starts before master acquisition.        *
 * - the first few slave lines do not overlap with master.      *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
real8 getAzOffset(const input_ell &ellips,
                  const slcimage  &master,
		  const slcimage  &slave,
                        orbit     &masterorbit,
                        orbit     &slaveorbit)
{
  const uint cen_pix = (slave.currentwindow.pixlo +slave.currentwindow.pixhi)/2;
  cn P;
  real8 tAzi, tRg;
  lp2xyz(1,cen_pix,ellips,slave,slaveorbit,P);
  xyz2t(tAzi,tRg,master,masterorbit,P);
  return master.t_azi1-tAzi;
} // END getAzOffset



/****************************************************************
 *    observationdata::observationdata                          *
 *                                                              *
 * Reads interferometric phases, coherence estimates and        *
 * and heights from file (bufferwise), selects observations,    *
 * computes weights and stores all necessary information in     *
 * private class members.                                       *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
observationdata::observationdata(const input_estorbits &estorbitsinput,
				 const input_gen       &generalinput,
				 const input_ell       &ellips,
				 const slcimage        *master,         
				 const slcimage        *slave,
				 const slcimage        *reference,
				       orbit           *masterorbit,      
                                       orbit           *slaveorbit,      
                                       orbit           *referenceorbit,
				 const productinfo     &interferogram,
				 const productinfo     &coherence)
  : ellips(ellips),
    master(master),
    slave(slave),
    masterorbit(masterorbit),
    slaveorbit(slaveorbit),
    referenceorbit(referenceorbit),
    useheights(specified(estorbitsinput.fiheightmap)),
    degree(estorbitsinput.poldegree),
    interferogram(interferogram),
    npar(2*(estorbitsinput.poldegree+1)),
    tmin(reference->t_azi1),
    tmax(tmin+(reference->originalwindow.lines()-1)/reference->prf)
{
  TRACE_FUNCTION("observationdata::observationdata (HB 17-Mar-2011)");
  const uint GRIDBORDER = 50;  // define a band of 50 (non-multilooked) pixels at both near
  // range and far range from where no pixels are used as observations. In ERS-SLC are some-
  // times processing-artefacts in these regions. Don't know about other sensors, but anyway,
  // discarding these small borders does not downgrade estimation quality significantly.

  const uint ifgLines   = uint(interferogram.win.lines() /interferogram.multilookL);
  const uint ifgPixels  = uint(interferogram.win.pixels()/interferogram.multilookP);
  const uint border     = ceil(real4(GRIDBORDER)/interferogram.multilookP);

  // ______ compute azimuth timing offset ______
  azoffset = getAzOffset(ellips,*reference,*master,*referenceorbit,*masterorbit);
  t_start = tmin - azoffset;
  DEBUG << "Azimuth timing offset w. r. t. reference orbit: " << azoffset;
  DEBUG.print();

  // ______ reading of coherence image required? ______
  bool readcoherence = true;  
  if (specified(estorbitsinput.ifpositions) && estorbitsinput.weighting==eo_noweighting && estorbitsinput.threshold==0)
    readcoherence = false;

  // ______ determine buffersize, check file format of ifg. ______
  uint numberofbigimages = 2; // 1 for ifg, +1 for the small stuff
  if (interferogram.formatflag == FORMATCR4)
    numberofbigimages += 2;
  else if (!(interferogram.formatflag == FORMATR4))
    { 
      ERROR << interferogram.file << " must be real4 or complr4.";
      ERROR.print();
      throw(unhandled_case_error);
    }
  if (readcoherence)
    {
      numberofbigimages++;
      if (coherence.formatflag == FORMATCR4)
	numberofbigimages += 2;
      else if (!(coherence.formatflag == FORMATR4))
	{ 
	  ERROR << coherence.file << " must be real4 or complr4.";
	  ERROR.print();
	  throw(unhandled_case_error);
	}
    }
  if (useheights)
    numberofbigimages++;
  uint bufferlines = generalinput.memory/ifgPixels/4/numberofbigimages;
    
  // ______ read positions from file or select from the cells of a grid ______
  bool readposfromfile;  
  if (specified(estorbitsinput.ifpositions))
    {
      readposfromfile = true;
      gridsize = 0;   
      gridlines = 0;  
      gridpixels = 0;
      maxobs = estorbitsinput.nobs;
    }
  else
    {
      readposfromfile = false;
      gridsize = floor(sqrt(ifgLines*ifgPixels/estorbitsinput.nobs));
      gridlines  = ceil(ifgLines/real4(gridsize));
      gridpixels = ceil((ifgPixels-border)/real4(gridsize));
      maxobs = gridlines*gridpixels;
      bufferlines = (bufferlines/gridsize)*gridsize;  // integer division!
    }
 
  // ______ check buffersize ______
  if (bufferlines>ifgLines) 
    bufferlines = ifgLines;
  else if (bufferlines<1)
    {
      PRINT_ERROR("Please increase memory (MEMORY card) or decrease multiL (INT_MULTILOOK card).")
	throw(input_error);
    }
  const uint nbuffers = ceil(ifgLines/real4(bufferlines));

  // ______ initialise storage matrices ______
  line = matrix<uint>(maxobs,1);
  pixel = matrix<uint>(maxobs,1);
  phi = matrix<real4>(maxobs,1);
  weight = matrix<real4>(maxobs,1);
  if (useheights) height = matrix<real4>(maxobs,1);

  // ______ read positions from file EO_IN_POS ______
  if (readposfromfile)
    {
      ifstream ifpos;
      char dummyline[ONE27];
      openfstream(ifpos,estorbitsinput.ifpositions);
      bk_assert(ifpos,estorbitsinput.ifpositions,__FILE__,__LINE__);
      PROGRESS << "Reading positions of observations from file: " << estorbitsinput.ifpositions;
      PROGRESS.print();

      uint ll,pp;
      for (uint i=0; i<maxobs; ++i)
	{
	  ifpos >> ll >> pp;
	  if (ll<1 || pp<1 || ll>ifgLines || pp>ifgPixels)
	    {
	      ERROR << "Pixel location (" << ll << "," << pp << ") is outside the multilooked crop. ";
	      ERROR.print();
	      throw(unhandled_case_error);
	    }
	  line(i,0) = uint(ll)-1;               
	  pixel(i,0) = uint(pp)-1;              
	  ifpos.getline(dummyline,ONE27,'\n');  
	}
      ifpos.close();

      // ______ Check last point ivm.EOL after last position in file ______
      if (line(maxobs-1,0) == line(maxobs-2,0) &&
	  pixel(maxobs-1,0) == pixel(maxobs-2,0))
	{
	  maxobs--;
	  WARNING << "EO: there should be no EOL after last point in file: "
		  << estorbitsinput.ifpositions;
	  WARNING.print();
	}
      DEBUG << maxobs << " positions read from file " << estorbitsinput.ifpositions;
      DEBUG.print();
    }
  
  // ______ open files with ifg, coherence and height ______
  ifstream ifint, ifcoh, ifheight;
  openfstream(ifint,interferogram.file);
  bk_assert(ifint,interferogram.file,__FILE__,__LINE__);
  if (readcoherence)
    {
      openfstream(ifcoh,coherence.file);
      bk_assert(ifcoh,coherence.file,__FILE__,__LINE__);
      cerr << coherence.file << endl;
    }
  else
    DEBUG.print("Coherence image will not be read, since values are not required..");
  if (useheights)
    {
      openfstream(ifheight,estorbitsinput.fiheightmap);
      bk_assert(ifheight,estorbitsinput.fiheightmap,__FILE__,__LINE__);
    }
  else 
    {
      INFO.print("No DEM specified with EO_IN_DEM_LP.");
      INFO.print("Assuming flat topography for orbit error estimation.");
    }


  // ====== read data bufferwise ======
  PROGRESS << "Reading image data in " << nbuffers << " buffer(s).";
  PROGRESS.print();
  nobs = 0;
  uint linesdone = 0;
  uint restlines = ifgLines;
  for (uint buf=0; buf<nbuffers; buf++) 
    {
      uint linesInThisBuffer;
      if (restlines>=bufferlines)
	linesInThisBuffer=bufferlines;
      else
	linesInThisBuffer=restlines;
      DEBUG << "Reading buffer " << buf+1 << ", (multilooked) lines " << linesdone+1 
	    << " to " <<  linesdone+linesInThisBuffer+1 << ".";
      DEBUG.print();

      matrix<real4> INT, COH, HEIGHT;
      if (interferogram.formatflag == FORMATR4)
	{
	  INT = matrix<real4>(linesInThisBuffer,ifgPixels);
	  ifint >> INT;
	}
      else
	{
	  matrix<complr4> CINT(linesInThisBuffer,ifgPixels);
	  ifint >> CINT;
	  INT = angle(CINT);
	}
      if (readcoherence)
	{
	  if (coherence.formatflag == FORMATR4)
	    {
	      COH = matrix<real4>(linesInThisBuffer,ifgPixels);
	      ifcoh >> COH;
	    }
	  else
	    {
	      matrix<complr4> CCOH(linesInThisBuffer,ifgPixels);
	      ifcoh >> CCOH;
	      COH = matrix<real4>(linesInThisBuffer,ifgPixels);
	      for (register uint l=0; l<linesInThisBuffer; l++)
		for (register uint p=0; p<ifgPixels; p++)
		  COH(l,p) = abs(CCOH(l,p));
	    }
	}
	  
      if (useheights) 
	{
	  HEIGHT = matrix<real4>(linesInThisBuffer,ifgPixels);
	  ifheight >> HEIGHT;
	}

      // ______ either use positions from file ______
      if (readposfromfile)
	{
	  uint lmax = linesdone + bufferlines;
	  for (uint i=0; i<maxobs; i++)
	    if (line(i,0)>=linesdone && line(i,0)<lmax)
	      {
		uint l = line(i,0)-linesdone;
		if (readcoherence) 
		  {
		    if (COH(l,pixel(i,0))<estorbitsinput.threshold) continue;
		    weight(nobs,0) = COH(l,pixel(i,0));
		  }
		phi(nobs,0) = INT(l,pixel(i,0));
		if (useheights) height(nobs,0) = HEIGHT(l,pixel(i,0));
		nobs++;
	      }
	}

      // ______ or select observations by selecting most coherent pixels from a grid ______
      else 
	{
	  for (int32 l0=0; l0<linesInThisBuffer; l0+=gridsize)
	    {
	      uint lmax = l0 + min(gridsize,linesInThisBuffer-l0);
	      for (int32 p0=border; p0<ifgPixels-border; p0+=gridsize)
		{
		  uint pmax = p0 + min(gridsize,ifgPixels-border-p0);
		  real4 maxcoh = 0;
		  int32 maxcoh_l = -1;
		  int32 maxcoh_p = -1;
		  
		  for (uint l=l0; l<lmax; l++)
		    for (uint p=p0; p<pmax; p++)
		      {
			if (COH(l,p)>maxcoh)
			  {
			    maxcoh = COH(l,p);
			    maxcoh_l = l;
			    maxcoh_p = p;
			  }
		      }
		  if (maxcoh>=estorbitsinput.threshold) 
		    {
		      phi(nobs,0) = INT(maxcoh_l,maxcoh_p);
		      weight(nobs,0) = COH(maxcoh_l,maxcoh_p);
		      if (useheights) height(nobs,0) = HEIGHT(maxcoh_l,maxcoh_p);
		      line(nobs,0) = linesdone+maxcoh_l;
		      pixel(nobs,0) = maxcoh_p;
		      nobs++;
		    }
		}
	    }
	}

      linesdone += bufferlines;
      restlines -= bufferlines;
    } // END buffer-loop

   ifint.close();
   ifcoh.close();
   ifheight.close();
   INFO << nobs << " of " << maxobs << " pixel observations have coherence >= " 
	<< estorbitsinput.threshold << ".";
   INFO.print();
   if (nobs<maxobs)
     {
       PROGRESS << maxobs-nobs << " observations are discarded due to lacking coherence.";
       PROGRESS.print();
     }

   // ______ apply weighting scheme ______
   switch(estorbitsinput.weighting)
     {
     case eo_noweighting:
       weight = 1;
       break;
     case eo_coh:
       // do nothing, keep coherence value
       break;
     case eo_coh2:
       for (uint i=0; i<nobs; i++)
	 weight(i,0) *= weight(i,0);  
       break;
     case eo_pdf:
       sigmaTable LUT(coherence.multilookL*coherence.multilookP);
       for (uint i=0; i<nobs; i++)
	 {
	   weight(i,0) = LUT.getSigma(weight(i,0));
	   weight(i,0) = 1/(weight(i,0)*weight(i,0));
	 }
     }

   // ______ resize storage matrices ______
   if (nobs>0)
     {
       if (nobs != maxobs)
	 {
	   phi = phi.getdata(window(0,nobs-1,0,0));
	   weight = weight.getdata(window(0,nobs-1,0,0));
	   if (useheights) height = height.getdata(window(0,nobs-1,0,0));
	   line = line.getdata(window(0,nobs-1,0,0));
	   pixel = pixel.getdata(window(0,nobs-1,0,0));
	 }
       use = matrix<bool>(nobs,1);
       use = true;
       residual = matrix<real8>(nobs,1);
       test = matrix<real4>(nobs,1);
     }
} // END observationdata::observationdata



/****************************************************************
 *    observationdata::lineOfDesignMatrix                       *
 *                                                              *
 * compute design matrix coefficients for specific observation. *
 *                       .                                      *
 * input:                                                       *
 *   - index of the observation                                 *
 *   - output version, defined by one of three constants:       *
 *                                                              *
 * DESIGN_NORMAL:   (ah0,ah1,...,av0,av1)                       *
 * DESIGN_EXTENDED: (ah0,ah1,...,av0,av1,phi)                   *
 * DESIGN_DUMP:     (tAzi,-ah0m,-av0m,ah0s,av0s)                *
 *                                                              *
 * tAzi: azimuth time normalised to [-2,2]                      *
 * ah0:  coefficient for Bh, averaged for master and slave      *
 * ah1:  tAzi*ah0                                               *
 * av0:  coefficient for Bv, averaged for master and slave      *
 * av1:  tAzi*av0                                               *
 * ah0m: coefficient for Bh, computed for master                *
 * avm0: coefficient for Bv, computed for master                *
 * ah0m: coefficient for Bh, computed for slave                 *
 * avm0: coefficient for Bv, computed for slave                 *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
matrix<real8> observationdata::lineOfDesignMatrix(const uint index,int16 version)
{
  TRACE_FUNCTION("observationdata::lineOfDesignMatrix (HB 17-Mar-2011)");

  // ______ compute non-multilooked line/pixel and azimuth time ______
  cn er, ea, ex;
  real8 ln = interferogram.win.linelo - .5 + (2.*line(index,0)+1)*interferogram.multilookL/2.;
  real8 tAzi = t_start + (ln-1)/master->prf;
  getCoordinateFrame(ea,er,ex,tAzi,*referenceorbit);
  real8 px = interferogram.win.pixlo -0.5 + (2.*pixel(index,0)+1)*interferogram.multilookP/2.;
  tAzi = normalize(tAzi,tmin,tmax);  // normalise tAzi now for later use

  // ______ consider terrain height ______
  input_ell ELLIPS = ellips;
  if (useheights)
    {
      ELLIPS.a += height(index,0);
      ELLIPS.b += height(index,0);
    }

  // ______ find corresponding pos. on master orbit (M), slave orbit (S) and surface (P) ______
  cn M, S, P;
  M = (*masterorbit).getxyz(master->t_azi1+(ln-1)/master->prf); // compute point on master orbit
  lp2xyz(ln, px, ELLIPS, *master, *masterorbit, P);          // find corresponding surface point
  xyz2orb(S, *slave, *slaveorbit, P);                         // find point on slave orbit

  // ______ initialise Ai ______
  matrix<real8> Ai;
  switch (version)
    {
    case DESIGN_NORMAL  : Ai = matrix<real8>(1,npar);   break;
    case DESIGN_EXTENDED: Ai = matrix<real8>(1,npar+1); break;
    case DESIGN_DUMP    : Ai = matrix<real8>(1,5);      break;
    }

  // ______ fill Ai ______
  switch (version)
    {
    case DESIGN_NORMAL: 
    case DESIGN_EXTENDED:
      {
	// compute coefficients
	cn R_scaled = ( (P-M).normalize()/master->wavelength + (P-S).normalize()/slave->wavelength) *2*PI;
	real8 ah = -R_scaled.in(ex);   // coefficient for horizontal baseline component
	real8 av = -R_scaled.in(er);   // coefficient for vertical baseine component
	
	// setup line of design matrix 
	for (uint i=0; i<=degree; i++)
	  {
	    Ai(0,i) = ah;
	    Ai(0,degree+1+i) = av;
	  }
	for (uint i=1; i<=degree; i++)
	  for (uint j=i; j<=degree; j++)
	    {
	      Ai(0,j) *= tAzi;
	      Ai(0,degree+1+j) *= tAzi;
	    }
	
	if (version==DESIGN_EXTENDED)
	  Ai(0,npar) = real8(phi(index,0));
      }
      break;

    case DESIGN_DUMP:
      cn P_M_norm = (M-P).normalize();
      cn P_S_norm = (S-P).normalize();
	    
      Ai(0,0) = tAzi;                                        // normalised azimuth time
      Ai(0,1) = -4*PI/master->wavelength * P_M_norm.in(ex);  // across track coefficient for master
      Ai(0,2) = -4*PI/master->wavelength * P_M_norm.in(er);  // radial coefficient for master
      Ai(0,3) =  4*PI/slave->wavelength  * P_S_norm.in(ex);  // across track coefficient for slave
      Ai(0,4) =  4*PI/slave->wavelength  * P_S_norm.in(er);  // radial coefficient for slave
      break;
    }

  return Ai;
} // END observationdata::lineOfDesignMatrix



/****************************************************************
 *    writeResiduals                                            *
 *                                                              *
 * Write some information on individual observation pixels to   *
 * file.                 .                                      *
 *                                                              *
 * Hermann Baehr, 17-Mar-2011                                   *
 ****************************************************************/
void writeResiduals(const input_estorbits &estorbitsinput, 
		    const observationdata &obs)
{
  ofstream resdata(estorbitsinput.foresiduals, ios::out | ios::trunc);
  bk_assert(resdata,"estorbits: foresiduals",__FILE__,__LINE__);
  if (estorbitsinput.method==eo_lsq)
    resdata << "This file contains residuals from orbit error"
	    << "\nestimation. There are seven columns with: observation"
	    << "\nnumber, line (multilooked), pixel (multilooked)," 
	    << "\nphase [rad], weight, residual [rad], test."
	    << "\n number posLm posPm    phase   weight residual    test"
	    << "\n------------------------------------------------------\n";
  else if (estorbitsinput.method==eo_gridsearch)
    resdata << "This file contains residuals from orbit error"
	    << "\nestimation. There are five columns with: observation"
	    << "\nnumber, line (multilooked), pixel (multilooked)," 
	    << "\nphase [rad], weight."
	    << "\n number posLm posPm    phase   weight"
	    << "\n------------------------------------------------------\n";
  resdata.close();
  FILE *resfile;
  resfile = fopen(estorbitsinput.foresiduals,"a");
  if (estorbitsinput.method==eo_lsq)
    {
      for (uint i=0; i<obs.nInit(); i++)
	if (obs.isUsed(i))
	  fprintf(resfile,"%7i %5i %5i %8.3f %8.3f %8.3f %7.2f\n",
		  i+1, obs.getLine(i)+1, obs.getPixel(i)+1, obs(i),
		  obs.getWeight(i), obs.getRes(i), abs(obs.getT(i)));
    }
  else if (estorbitsinput.method==eo_gridsearch)
    {
      for (uint i=0; i<obs.nInit(); i++)
	if (obs.isUsed(i))
	  fprintf(resfile,"%7i %5i %5i %8.3f %8.3f\n",
		  i+1, obs.getLine(i)+1, obs.getPixel(i)+1, obs(i),
		  obs.getWeight(i));
    }
  fclose(resfile);
} // END writeResiduals
