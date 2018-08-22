// Copyright (c) 2008, 2009, 2010 Regents of the University of California.
//
// ADModelbuilder and associated libraries and documentations are
// provided under the general terms of the "BSD" license.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2.  Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3.  Neither the name of the  University of California, Otter Research,
// nor the ADMB Foundation nor the names of its contributors may be used
// to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
// MultiSpecies Size Structured Assessment Model MS3AM
// 		based on LeMANS as modified by M. Fogarty and S. Gaichas
//		coded by S. Gaichas, help from K. Curti, G. Fay, T. Miller
//		testing version, February 2012
//		initial working version, July 2012
//
// Renamed hydra_sim and modified for use as size structued simulation model
//		operating model for testing produciton, nonlinear time series models
//		May 2013
//
// Dec 29, 2016: Adapted by Andy Beet
//             : Dynamically change F (function of effort) as a response to
//               Guild level biomass levels and/or species biomass
//                dropping below some threshold. Several data lines added to .dat file
//                and additional code added. Assessment module
//=======================================================================================
GLOBALS_SECTION
//=======================================================================================
  //Including C++ libraries
//  #include <morefun.cxx>
  #include "statsLib.h"
  ofstream simout("simout.csv");
  ofstream test("test.csv");

//=======================================================================================
DATA_SECTION
//=======================================================================================

//Debug statements from Kiersten Curti's model
  //  0's = Main body of program
      //  =  1:  Exits after data section
      //  =  2:  Exits after parameter section
      //  =  3:  Exits at end of procedure section
      //  =  4:  Prints checkpoints after each function in the procedure section, except those within the year loop
      //  =  5:  Prints checkpoints after each function within the year loop
      //  =  6:  Prints parameter estimates after each iteration
      //  =  7:  Prints pop dy variables after each iteration
      //  =  8:  Prints trophic matrices at end of year loop, then exits
      //  =  9:  Prints out predicted indices that go into the objective function after each iteration
  //10's = Initial states function
      //  = 10:  Outputs food-selection parameters at the end of the function and then exits
      //  = 11:  Outputs food selection paramters at the end of each iteration
      //  = 12:  Outputs fishery and survey selectivity matrices at the end of the function and then exits
      //  = 13:  Outputs abundance and biomass arrays and then exits
      //  = 14:  Outputs Yr1, Age1 and iFt matrices to ensure 'means + devt'ns' parameterized correctly and then exits
      //  = 15:  Outputs initial N and proportion mature arrays
  //20's = Suitability function
      //  = 20:  Output at end of function
      //  = 21:  Output at end of function, then exits
      //  = 22:  Outputs suitability and scaled suitability for each predator sp and then exits
      //  = 23:  Outputs Eta, Sigma, WtRatio, G for pred, prey combos and then exits
  //30's = Predation mortality function
      //  = 30:  Prints checkpoints throughout the function
      //  = 31:  Prints intermediate arrays (Avail, scAv, Consum) for each predator sp and then exits
      //  = 32:  Prints intermediate arrays for each prey sp (Consum, sumCon) and then exits
      //  = 33:  Prints sumCon and B right before M2 is calculated, and M2 after it is calculated
      //  = 34:  Debug 31 but does not exit
      //  = 35:  Prints nf and Other Food in every year
  //40's = Population dynamics function
      //  = 40:  Outputs N, C_hat and Cprop_hat at end of function
      //  = 41:  Outputs N, C_hat and Cprop_hat at end of function and then exits
      //  = 42:  Outputs mortality components for each species in each year and exits at end of year loop after trophic =1
      //  = 43:  Outputs mortality components at end of function and then exits
  //50's = Survey abundance function
      //  = 50:  Prints intermediate arrays for species where survey data is one contiguous time series and then exits
      //  = 51:  Prints intermediate arrays for species where the survey data is split into multiple segments and then exits
      //  = 52:  Prints predicted q, FICs, FIC_hat and N for each species and then exits
      //  = 53:  Prints estimated q matrix at the end of each iteration
  //60's = Log likelihood function
      //  = 60: Prints checkpoints after each likelihood component
      //  = 61: Prints checkpoints for multinomial components within the year loop
      //  = 62: Prints predicted and log predicted indices for TotC and TotFIC
      //  = 63: Prints predicted and log predicted indices for Cprop
      //  = 64: Prints predicted and log predicted indices for Sprop
      //  = 65: FHprop, when added
      //  = 66: Prints summary of objective function components at end of each iteration
  //70's = Food habits function
      //  = 70:  Prints Avpd (scAv) and FHtmp for each predator, year and then exits
      //  = 71:  Prints Avpd, Avpy for each prey species within Avpd, the colsum of Avpy, and FHtmp; exits after all predator species
      //  = 72:  Prints bin years, FHtmp, FHsum, and average FH, FH_hat, for each FH bin, and exits after all predator species
      //  = 73:  Prints total %W for each pred age, summed across prey sp, for each pred sp and FH bin.  This value should always == 1.
  //80's = Penalty functions
      // = 80: Biomass penalty function: Prints pre- and post- B for every species and year, regardless of whether a penalty is imposed
      // = 81: Biomass penalty function: Prints pre- and post- biomass and assorted arrays when biomass falls below threshold
      // = 82: Yr1 penalty function: Prints avgZ, thYr1, Yr1 and Ypen arrays and then exits at end of function
      // = 83: Recruitment penalty function: Prints Age1 parameter estimates, calculated CV and recruitment penalty

//take random number seeds from command line, 
//code from https://groups.nceas.ucsb.edu/non-linear-modeling/projects/skate/ADMB/skatesim.tpl
//#1: spinydog
//#2: winterskate
//#3: Aherring
//#4: Acod
//#5: haddock
//#6: yellowtailfl
//#7: winterfl
//#8: Amackerel
//#9: silverhake
//#10: goosefish

  int sim;
  int rseed;
 LOC_CALCS
    int on,opt;
    sim=0;
//     rseed=1;
    rseed=123456;
    //the following line checks for the "-sim" command line option
    //if it exists the if statement retreives the random number seed
    //that is required for the simulation model
    if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
    {
    //  cout<<"in if statement \n";
      sim=1;
         rseed=atoi(ad_comm::argv[on+1]);
    //if(SimFlag)exit(1);
    }
 END_CALCS

  init_int debug;                          //Debugger switch, 0 off, other above, now in dat file
 
//read in bounding indices from .dat file
  init_int Nyrs  												//number of years
  init_int Nspecies												//number of species
  init_int Nsizebins											//number of size bins
  init_int Nareas												//number of areas
  init_int Nfleets												//number of fleets
//  init_int Nages												//number of age classes
  int Totsizebins
  !!  Totsizebins = Nspecies*Nsizebins;
  int spp                         //loop counter species
  int size                         //loop counter lengthbins
  int area                        //loop counter areas
  int pred								 //loop counter predators
  int prey								 //loop counter prey
  int t									 //loop counter timestep
  int yr								 //loop counter year
  int fleet              //loop counter fleet
  int Nstepsyr           //model timesteps per year, set using max(growthprob_phi) below
  int Tottimesteps       //total timesteps for dimensioning state variable arrays
  int yrct               //year counter for dynamics
  int iguild             // loop counter in calc_survey_abundance
  int iassess            // loop counter in calc_assessment_strategy
  int ithreshold         // loop counter in calc_assessment_strategy
  ivector maxThreshold(1,Nareas) // stored most severe  threshold detected

  number o 
  !!  o = 0.001;         //small number to add before taking logs in objective function calcs
  // should add 1 since log of 1 = 0. objecive function not changed! Andy Beet
  
  init_number wtconv     //multiply by weight in grams to desired units for biomass, catch
 
//read in length bin sizes, weight at length parameters from .dat file
  init_matrix binwidth(1,Nspecies,1,Nsizebins) 					//length bin width in cm
  //init_3darray binwidth(1,Nareas,1,Nspecies,1,Nsizebins) 		//length bin width in cm, area specific
  init_vector lenwt_a(1,Nspecies)				//weight = lenwt_a * length ^ lenwt_b UNITS cm to g
  init_vector lenwt_b(1,Nspecies)				//weight = lenwt_a * length ^ lenwt_b UNITS cm to g
  //init_matrix lenwt_a(1,Nareas,1,Nspecies) 					//weight = lenwt_a * length ^ lenwt_b, area specific
  //init_matrix lenwt_b(1,Nareas,1,Nspecies) 					//weight = lenwt_a * length ^ lenwt_b, area specific

//calculate length and weight attributes of bins
//dimension all by area and add outside area loop to make these area specific
  matrix lbinmax(1,Nspecies,1,Nsizebins)						//upper end of length bins 
  matrix lbinmin(1,Nspecies,1,Nsizebins)						//lower end of length bins
  matrix lbinmidpt(1,Nspecies,1,Nsizebins)						//midpoint of length bins
  matrix wtbinmax(1,Nspecies,1,Nsizebins)						//max weight of length bins
  matrix wtbinmin(1,Nspecies,1,Nsizebins)						//min weight of length bins
  matrix wtatlbinmidpt(1,Nspecies,1,Nsizebins)					//wt at length midpoint of length bins
  matrix binavgwt(1,Nspecies,1,Nsizebins)						//average weight for length bins
  matrix powlbinmaxb(1,Nspecies,1,Nsizebins)
  !!	for (spp=1; spp<=Nspecies; spp++){
  !!		lbinmax(spp,1)   = binwidth(spp,1);
  !!		lbinmin(spp,1)   = 0.0;						//lowest bin assumed to start at 0!
  !!		lbinmidpt(spp,1) = binwidth(spp,1)/2.0;		
  !!		for (size=2; size<=Nsizebins; size++){
  !!			lbinmax(spp, size)   = binwidth(spp, size) + lbinmax(spp, size-1);
  !!			lbinmin(spp, size)   = lbinmax(spp, size-1);
  !!			lbinmidpt(spp, size) = binwidth(spp, size)/2.0 + lbinmax(spp, size-1);
  !!		}
  !!    	wtbinmax(spp) = lenwt_a(spp)* pow(lbinmax(spp), lenwt_b(spp));
  !!    	wtbinmin(spp) = lenwt_a(spp)* pow(lbinmin(spp), lenwt_b(spp));
  !!    	wtatlbinmidpt(spp) = lenwt_a(spp)* pow(lbinmidpt(spp), lenwt_b(spp));
  !!	}	
  !!	binavgwt = (wtbinmin + wtbinmax)/2.0;

//read in covariate information from .dat file
  init_int Nrecruitment_cov  									//number of recruitment covariates
  init_int Nmaturity_cov  										//number of maturity covariates
  init_int Ngrowth_cov  										//number of growth covariates
  init_matrix recruitment_cov(1,Nrecruitment_cov,1,Nyrs)		//time series of recruitment covariates
  init_matrix maturity_cov(1,Nmaturity_cov,1,Nyrs)				//time series of maturity covariates
  init_matrix growth_cov(1,Ngrowth_cov,1,Nyrs)					//time series of growth covariates
  //init_3darray recruitment_cov(1,Nareas,1,Nrecruitment_cov,1,Nyrs)  //time series of recruitment covariates, area specific
  //init_3darray maturity_cov(1,Nareas,1,Nmaturity_cov,1,Nyrs)  //time series of maturity covariates, area specific
  //init_3darray growth_cov(1,Nareas,1,Ngrowth_cov,1,Nyrs)   	//time series of growth covariates, area specific
 
//read in survey and catch observations from .dat file
  init_3darray obs_survey_biomass(1,Nareas,1,Nspecies,1,Nyrs)  	//spring or fall? units needed
  init_3darray obs_catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  	//total catch in tons
  init_3darray obs_effort(1,Nareas,1,Nfleets,1,Nyrs)  	//standardized effort units needed

  //init_4darray for survey size comp by area, species, year? 
  //init_5darray for catch at size by area, species, fleet, year? 
  
//read in mean stomach content weight time series from .dat file for intake calculation
  init_4darray mean_stomwt(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)     
  
  //want variance for this in this for fitting? 
  
//read in temperature time series from .dat file for intake calculation
  init_matrix obs_temp(1,Nareas,1,Nyrs)       //want variance measure? data source?
  
//read in estimation phases from .dat file
  init_int yr1Nphase            //year 1 N at size estimation phase
  init_int recphase				//recruitment parameter estimation phase
  init_int avg_rec_phase		//average recruitment estimation phase (could make species specific, currently global)
  init_int avg_F_phase			//average fishing mort estimation phase (could make species specific, currently global)
  init_int dev_rec_phase		//recruitment deviation estimation phase (could make species specific, currently global)
  init_int dev_F_phase			//fishing mort deviation estimation phase (could make species specific, currently global)
  init_int fqphase              //fishery q estimation phase
  init_int sqphase              //survey q estimation phase
  init_int ssig_phase           //survey sigma (obs error) phase
  init_int csig_phase           //catch sigma (obs error) phase
  
//read in lists of species names, area names, fleet names, actual years from .dat file

//the following are parameters that will be fixed in initial runs, so read in as "data" from .dat
//to estimate any of them within the model, place in parameter section, initialize from .pin file
//recruitment parameters from .dat file
  init_matrix recGamma_alpha(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model alpha
  init_matrix recGamma_shape(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model shape parameter
  init_matrix recGamma_beta(1,Nareas,1,Nspecies)			//eggprod gamma Ricker model beta
  
  init_matrix recDS_alpha(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model alpha
  init_matrix recDS_shape(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model shape parameter
  init_matrix recDS_beta(1,Nareas,1,Nspecies)			//SSB Deriso-Schnute model beta

  init_matrix recGamSSB_alpha(1,Nareas,1,Nspecies)			//SSB gamma alpha
  init_matrix recGamSSB_shape(1,Nareas,1,Nspecies)			//SSB gamma shape parameter
  init_matrix recGamSSB_beta(1,Nareas,1,Nspecies)			//SSB gamma beta

  init_matrix recRicker_alpha(1,Nareas,1,Nspecies)			//SSB Ricker model alpha
  init_matrix recRicker_shape(1,Nareas,1,Nspecies)			//SSB Ricker model shape parameter=1.0 not used
  init_matrix recRicker_beta(1,Nareas,1,Nspecies)			//SSB Ricker model beta

  init_matrix recBH_alpha(1,Nareas,1,Nspecies)			//SSB Beverton Holt model alpha
  init_matrix recBH_shape(1,Nareas,1,Nspecies)			//SSB Beverton Holt model shape parameter=1.0 not used
  init_matrix recBH_beta(1,Nareas,1,Nspecies)			//SSB Beverton Holt model beta

  init_matrix recShepherd_alpha(1,Nareas,1,Nspecies)			//SSB Shepherd model alpha
  init_matrix recShepherd_shape(1,Nareas,1,Nspecies)			//SSB Shepherd model shape parameter=1.0 not used
  init_matrix recShepherd_beta(1,Nareas,1,Nspecies)			//SSB Shepherd model beta


  init_ivector rectype(1,Nspecies)  //switch for alternate recruitment functions
  
  init_ivector stochrec(1,Nspecies)  //switch for stochastic recruitment
   
  matrix rec_alpha(1,Nareas,1,Nspecies)
  matrix rec_shape(1,Nareas,1,Nspecies)
  matrix rec_beta(1,Nareas,1,Nspecies)
 
  //initialize recruitment parameters for each type (case 9 no functional form, dummy pars)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!	  switch (rectype (spp)){
  !!       case 1:	  				//egg production based recruitment, 3 par gamma (Ricker-ish) 
  !!		  rec_alpha(area,spp) = recGamma_alpha(area,spp);
  !!		  rec_shape(area,spp) = recGamma_shape(area,spp);
  !!		  rec_beta(area,spp) = recGamma_beta(area,spp);
  !!	   break;
  !!	   case 2:                   //SSB based recruitment, 3 par Deriso-Schnute; see Quinn & Deriso 1999 p 95
  !!          rec_alpha(area,spp) = recDS_alpha(area,spp);
  !!          rec_shape(area,spp) = recDS_shape(area,spp);
  !!          rec_beta(area,spp) = recDS_beta(area,spp);
  !!       break;
  !!	   case 3:                   //SSB based recruitment, 3 par gamma
  !!          rec_alpha(area,spp) = recGamSSB_alpha(area,spp);
  !!          rec_shape(area,spp) = recGamSSB_shape(area,spp);
  !!          rec_beta(area,spp) = recGamSSB_beta(area,spp);
  !!       break;
  !!	   case 4:                   //SSB based recruitment, 2 par Ricker
  !!          rec_alpha(area,spp) = recRicker_alpha(area,spp);
  !!          rec_shape(area,spp) = recRicker_shape(area,spp);
  !!          rec_beta(area,spp) = recRicker_beta(area,spp);
  !!       break;
  !!	   case 5:                   //SSB based recruitment, 2 par BevHolt
  !!          rec_alpha(area,spp) = recBH_alpha(area,spp);
  !!          rec_shape(area,spp) = recBH_shape(area,spp);
  !!          rec_beta(area,spp) = recBH_beta(area,spp);
  !!       break;
  !!       case 6:                // SSB based recruitment, 3 parameters Shepherd
  !!          rec_alpha(area,spp) = recShepherd_alpha(area,spp);
  !!          rec_shape(area,spp) = recShepherd_shape(area,spp);
  !!          rec_beta(area,spp) = recShepherd_beta(area,spp);
  !!       break;
  !!	   case 9:                   //no functional form, uses average+devs in .pin file
  !!          rec_alpha(area,spp) = 0;
  !!          rec_shape(area,spp) = 0;
  !!          rec_beta(area,spp) = 0;
  !!       break;
  !!       default:
  !!          cout<<"undefined recruitment type, check .dat file"<<endl;
  !!          exit(1);
  !!		}
  !!    }
  !! }

  init_matrix sexratio(1,Nareas,1,Nspecies)  //this is proportion female
  init_matrix recruitment_covwt(1,Nspecies,1,Nrecruitment_cov)	//recruitment covariate weighting factor
  //init_3darray recruitment_covwt(1,Nareas,1,Nspecies,1,Nrecruitment_cov) //area specific weighting
  
//fecundity parameters from .dat file and calculate fecundity at length
  init_matrix fecund_d(1,Nareas,1,Nspecies)
  init_matrix fecund_h(1,Nareas,1,Nspecies)
  init_3darray fecund_theta(1,Nareas,1,Nspecies,1,Nsizebins)
  3darray fecundity(1,Nareas,1,Nspecies,1,Nsizebins)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!    	fecundity(area, spp) = elem_prod(fecund_theta(area, spp),
  !!											(fecund_d(area, spp)
  !!                            	 			* pow(lbinmidpt(spp),fecund_h(area, spp))));
  !!	}	
  !!  }

//maturity parameters from .dat file 
  init_matrix maturity_nu(1,Nareas,1,Nspecies)
  init_matrix maturity_omega(1,Nareas,1,Nspecies)
  init_matrix maturity_covwt(1,Nspecies,1,Nmaturity_cov)
  matrix covariates_M(1,Nspecies,1,Nyrs) // intermediate calculation to obtain maturity covariates //AndyBeet
  
//growth parameters from .dat file and calculate simple (no cov) prob of growing through length interval
  init_matrix growth_psi(1,Nareas,1,Nspecies)    //power function growth length=psi*age^kappa
  init_matrix growth_kappa(1,Nareas,1,Nspecies)  //power function growth length=psi*age^kappa
  init_matrix growth_covwt(1,Nspecies,1,Ngrowth_cov)
  init_matrix vonB_Linf(1,Nareas,1,Nspecies)    //alternate parameterization, vonB growth
  init_matrix vonB_k(1,Nareas,1,Nspecies)       //alternate parameterization, vonB growth
  init_ivector growthtype(1,Nspecies)                          //switch for alternate growth types
  init_number phimax
  4darray growthprob_phi(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)              
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!     for(yr=1; yr<=Nyrs; yr++){
  !!      switch (growthtype (spp)){
  !!        case 1:	  	 //exponential no covariates 
  !!    	  growthprob_phi(area, spp, yr) = 1/pow((lbinmax(spp)/growth_psi(area, spp)),
  !!													(1.0/growth_kappa(area, spp)));
  !!        break;
  !!        case 2:       //exponential with covariates
  !!    	  growthprob_phi(area, spp, yr) = 1/pow((lbinmax(spp)/
  !!                                     growth_psi(area, spp)* mfexp(growth_covwt(spp)*trans(growth_cov)(yr))),
  !!													(1.0/growth_kappa(area, spp)));
  !!        break;
  !!        case 3:       //VonB no covariates
  !!          growthprob_phi(area, spp, yr) = vonB_k(area, spp)/log(
  !!                                          elem_div((vonB_Linf(area, spp)-lbinmin(spp)),
  !!                                                   (vonB_Linf(area, spp)-lbinmax(spp)))); 
  !!        break;
  !!        case 4:       //VonB with covariates
  !!          growthprob_phi(area, spp, yr) = vonB_k(area, spp)/log(
  !!                                          elem_div((vonB_Linf(area, spp)*mfexp(growth_covwt(spp)*trans(growth_cov)(yr))-lbinmin(spp)),
  !!                                                   (vonB_Linf(area, spp)*mfexp(growth_covwt(spp)*trans(growth_cov)(yr))-lbinmax(spp)))); 
  !!        break;
  !!        default:
  !!          cout<<"undefined growth type, check .dat file"<<endl;
  !!          exit(1);
  !!        }
  !!      growthprob_phi(area, spp, yr)(Nsizebins) = 0.0; //set prob of outgrowing highest bin to 0
  !!      double tempmax =  max(growthprob_phi(area, spp, yr));
  !!      phimax = max(tempmax,phimax);
  !!	  }	
  !!	}	
  !!    growthprob_phi(area) /= phimax;  //rescale so no group has >1 prob growing out
  !!  }
  !!//  growthprob_phi /= phimax;   //rescale so no group has >1 prob growing out--not working on 4d array
  !!  Nstepsyr = round(phimax);            //set model timestep to phimax
  !!  Tottimesteps = Nstepsyr*Nyrs;        //for scaling state variable arrays

//intake parameters from .dat file
  init_matrix intake_alpha(1,Nareas,1,Nspecies)
  init_matrix intake_beta(1,Nareas,1,Nspecies)
  4darray intake(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins)
  !!  for (area=1; area<=Nareas; area++){
  !!	for(spp=1; spp<=Nspecies; spp++){
  !!      for(yr=1; yr<=Nyrs; yr++){
  !!        intake(area, spp, yr) = 24.0 * (intake_alpha (area, spp) * 
  !!                              mfexp(intake_beta (area, spp) * obs_temp (area,yr))) *
  !!                              mean_stomwt(area, spp,yr) * //daily intake in g 
  !!                              365.0 / //annual intake 
  !!                              Nstepsyr;  //intake per model timestep
  !!      }
  !!    }
  !!  }
            
//natural mortality parameters from .dat file and calculate weight ratio, size preference, suitability
  init_3darray M1(1,Nareas,1,Nspecies,1,Nsizebins)
  init_3darray isprey(1,Nareas,1,Nspecies,1,Nspecies)    //preds in columns, prey in rows
  init_matrix preferred_wtratio(1,Nareas,1,Nspecies)     //pred specific, not size
  init_vector sd_sizepref(1,Nspecies)              //pred specific, not size
  4darray wtratio(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins)  //2nd dim pred spp, 3rd all spp as prey lengths, 4th pred lengths
  !!  for (area=1; area<=Nareas; area++){
  !!  	for (pred=1; pred<=Nspecies; pred++){             
  !!    	for(prey=1; prey<=Nspecies; prey++){
  !!			dmatrix wttemp = outer_prod(wtatlbinmidpt(prey), 1.0/wtatlbinmidpt(pred)); 
  !!			wttemp.rowshift(prey*Nsizebins-(Nsizebins-1));
  !!            wtratio(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) = wttemp;
  !!		}
  !!    }
  !!  }
//  old 4darray sizepref(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins) //2nd dim pred spp, 3rd all spp as prey lengths, 4th pred lengths
//  !!  //sizepref = exp(-square(log(wtratio)-preferred_wtratio))/2.0*variance_sizepref
//  !!  for (area=1; area<=Nareas; area++){
//  !!  	for (pred=1; pred<=Nspecies; pred++){ 
//  !!    	for(prey=1; prey<=Nspecies; prey++){
//  !!			  dmatrix logratioprey = log(wtratio(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins));  
//  !!			  sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) = 
//  !!									exp(-square(logratioprey-preferred_wtratio(area, pred))/
//  !!									2.0*var_sizepref(pred));
//  !!		}
//  !!	}
//  !!  } 
  4darray sizepref(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins) //2nd dim pred spp, 3rd all spp as prey lengths, 4th pred lengths
  !!  //sizepref = exp(-square(log(wtratio)-preferred_wtratio))/2.0*variance_sizepref--not full lognormal density!! full below
  !!  for (area=1; area<=Nareas; area++){
  !!  	for (pred=1; pred<=Nspecies; pred++){ 
  !!    	for(prey=1; prey<=Totsizebins; prey++){
  !!     	    for(int isize=1; isize<=Nsizebins; isize++){
  !!              double wtratioprey = wtratio(area, pred, prey, isize);
  !!			  sizepref(area, pred, prey, isize) = 
  !!              1/(wtratioprey*sd_sizepref(pred)*sqrt(2*M_PI))*exp(-square(log(wtratioprey)-preferred_wtratio(area, pred))/(2*square(sd_sizepref(pred))));
  !!           }
  !!		}
  !!	}
  !!  } 
  4darray suitability(1,Nareas,1,Nspecies,1,Totsizebins,1,Nsizebins)           
  !!  //suitability = sizepref * isprey
  !!  for (area=1; area<=Nareas; area++){
  !!  	for (pred=1; pred<=Nspecies; pred++){ 
  !!    	for(prey=1; prey<=Nspecies; prey++){
  !!			suitability(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins) = 
  !!						sizepref(area, pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins)*
  !!						isprey(area, prey, pred);
  !!		}
  !!	}
  !!  } 
  

  //fishery selectivity pars from dat file, for now not area specific
  init_matrix fishsel_c(1,Nspecies,1,Nfleets)  //fishery selectivity c par
  init_matrix fishsel_d(1,Nspecies,1,Nfleets)  //fishery selectivity d par

  // Added by andyBeet
  init_matrix B0(1,Nareas,1,Nspecies) // Equilibrium biomass. Obtained from a baseline run with zero fishing effort and no Errors added(recruitment, survey, catch)
  init_int Nguilds // number of guilds
  imatrix maxGuildThreshold(1,Nareas,1,Nguilds) // most severe exceedence for each guild
  imatrix maxSpeciesThreshold(1,Nareas,1,Nspecies) // most severe exceedence for each species. currently a binary response
  init_vector guildMembers(1,Nspecies) // assign each species to a guild. 1,2,3 etc.
  init_int AssessmentPeriod // time (yrs) when we assess guildlevel biomass levels

  matrix B0_guilds(1,Nareas,1,Nguilds) // equilibrium biomass for guild.
  // calculates the unfished equilibrium biomass of guild
  !! for (area=1; area<=Nareas; area++) {
  !!     for (iguild=1; iguild<=Nguilds; iguild++ ) {
  !!          for (spp=1; spp<=Nspecies; spp++) {
  !!               if (guildMembers(spp)== iguild) {
  !!                  // sum up equilibr biomass for each guild
  !!                  B0_guilds(area,iguild) += B0(area,spp);
  !!
  !!              }
  !!          }
  !!      }
  !! }

  init_int Nthresholds // number of thresholds used for corrective fishing rules
  init_vector threshold_proportion(1,Nthresholds) //levels at which action is taken
  init_vector exploitation_levels(1,Nthresholds) //levels to drop exploitation to if threshold is exceeded
  init_vector threshold_species(1,Nspecies) // individual species thresholds (fraction)
  init_int AssessmentOn // binary, yes or no
  init_int speciesDetection // binary yes or no. Determins if species level should influence exploitation rate change during assessment

  init_int LFI_size // determins the size deemed a large fish. Used in LFI calc_health_indices module


//flag marking end of file for data input          
  init_int eof;

//debugging section, check inputs and initial calculations
	LOCAL_CALCS

  if (debug == 1)
    {
    cout<<"Nyrs\n"<<Nyrs<<endl;    
    cout<<"Nspecies\n"<<Nspecies<<endl;
    cout<<"Nsizebins\n"<<Nsizebins<<endl;
    cout<<"Nareas\n"<<Nareas<<endl;
    cout<<"Nfleets\n"<<Nfleets<<endl;
    cout<<"wtconv\n"<<wtconv<<endl;
    cout<<"Totsizebins\n"<<Totsizebins<<endl;
    cout<<"binwidth\n"<<binwidth<<endl;
    cout<<"lenwt_a\n"<<lenwt_a<<endl;
    cout<<"lenwt_b\n"<<lenwt_b<<endl;
    cout<<"lbinmax\n"<<lbinmax<<endl;
    cout<<"lbinmin\n"<<lbinmin<<endl;
    cout<<"lbinmidpt\n"<<lbinmidpt<<endl;
    cout<<"wtbinmax\n"<<wtbinmax<<endl;
    cout<<"wtbinmin\n"<<wtbinmin<<endl;
    cout<<"wtatlbinmidpt\n"<<wtatlbinmidpt<<endl;
    cout<<"binavgwt\n"<<binavgwt<<endl;
    cout<<"Nrecruitment_cov\n"<<Nrecruitment_cov<<endl;
    cout<<"Nmaturity_cov\n"<<Nmaturity_cov<<endl;
    cout<<"Ngrowth_cov\n"<<Ngrowth_cov<<endl;
    cout<<"recruitment_cov\n"<<recruitment_cov<<endl;
    cout<<"maturity_cov\n"<<maturity_cov<<endl;
    cout<<"growth_cov\n"<<growth_cov<<endl;
    cout<<"obs_survey_biomass\n"<<obs_survey_biomass<<endl;
    cout<<"obs_catch_biomass\n"<<obs_catch_biomass<<endl;
    cout<<"mean_stomwt\n"<<mean_stomwt<<endl;
    cout<<"obs_temp\n"<<obs_temp<<endl;
    cout<<"recruitment_covwt\n"<<recruitment_covwt<<endl;
    cout<<"rectype\n"<<rectype<<endl;
    cout<<"stochrec\n"<<stochrec<<endl;
    cout<<"rec_alpha\n"<<rec_alpha<<endl;
    cout<<"rec_shape\n"<<rec_shape<<endl;
    cout<<"rec_beta\n"<<rec_beta<<endl;
    cout<<"fecund_d\n"<<fecund_d<<endl;
    cout<<"fecund_h\n"<<fecund_h<<endl;
    cout<<"fecund_theta\n"<<fecund_theta<<endl;
    cout<<"fecundity\n"<<fecundity<<endl;
    cout<<"maturity_nu\n"<<maturity_nu<<endl;
    cout<<"maturity_omega\n"<<maturity_omega<<endl;
    cout<<"maturity_covwt\n"<<maturity_covwt<<endl;
    cout<<"growth_psi\n"<<growth_psi<<endl;
    cout<<"growth_kappa\n"<<growth_kappa<<endl;
    cout<<"growth_covwt\n"<<growth_covwt<<endl;
    cout<<"vonB_Linf\n"<<vonB_Linf<<endl;
    cout<<"vonB_k\n"<<vonB_k<<endl;
    cout<<"growthtype (1 power, 2 power/covariates, 3 vonB, 4 vonB covariates)\n"<<growthtype<<endl;
    cout<<"growthprob_phi\n"<<growthprob_phi<<endl;
    cout<<"phimax\n"<<phimax<<endl;
    cout<<"Nstepsyr\n"<<Nstepsyr<<endl;
    cout<<"Tottimesteps\n"<<Tottimesteps<<endl;
    cout<<"intake_alpha\n"<<intake_alpha<<endl;
    cout<<"intake_beta\n"<<intake_beta<<endl;
    cout<<"intake\n"<<intake<<endl;
    cout<<"M1\n"<<M1<<endl;
    cout<<"isprey\n"<<isprey<<endl;
    cout<<"preferred_wtratio\n"<<preferred_wtratio<<endl;
    cout<<"sd_sizepref\n"<<sd_sizepref<<endl;
    cout<<"wtratio\n"<<wtratio<<endl;
    cout<<"sizepref\n"<<sizepref<<endl;
    cout<<"suitability\n"<<suitability<<endl;
    cout<<"B0\n"<<B0<<endl;
    cout<<"Nguilds\n"<<Nguilds<<endl;
    cout<<"guildMembers\n"<<guildMembers<<endl;
    cout<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;

//    cout<<setprecision(10)<<"FH\n"<<FH<<endl;
//    cout<<setprecision(10)<<"FHideal\n"<<FHideal<<endl;
    cout<<"eof\n"<<eof<<endl;
    }

  if(eof != 54321) {cout<<"Stop, data input failed"<<endl<<"eof: "<<eof<<endl; exit(1);}

  if (debug == 1) {cout<<"\nManually exiting at end of data section..."<<endl;  exit(-1);}

	END_CALCS
	
//=======================================================================================
INITIALIZATION_SECTION
//=======================================================================================
 
//=======================================================================================
PARAMETER_SECTION
//=======================================================================================
  matrix effort_updated(1,Nareas,1,Nfleets)     // calculated new Effort for each fleet following threshold exceedence calc_assessment_strategy

  //Initial N year 1
  init_3darray yr1N(1,Nareas,1,Nspecies,1,Nsizebins,yr1Nphase)       //initial year N at size, millions
  
  //recruitment parameters (alts in .dat file read in with switch for rec form by spp)
  init_matrix recruitment_alpha(1,Nareas,1,Nspecies,recphase)			//recruitment model alpha
  init_matrix recruitment_shape(1,Nareas,1,Nspecies,recphase)			//recruitment model shape parameter
  init_matrix recruitment_beta(1,Nareas,1,Nspecies,recphase)			//recruitment model beta
 
  //proportion mature
  4darray propmature(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins) //from maturity pars and covs

  //egg production
  3darray eggprod(1,Nareas,1,Nspecies,1,Nyrs) //from fecundity, propmature, sexratio, N, in millions
  
  //recruitment: average annual, annual devs, actual (avg+dev)
  init_matrix avg_recruitment(1,Nareas,1,Nspecies,avg_rec_phase)  //average annual recruitment by area, species
  init_3darray recruitment_devs(1,Nareas,1,Nspecies,1,Nyrs,dev_rec_phase)  //recruitment deviations by area, species
  3darray recruitment(1,Nareas,1,Nspecies,1,Nyrs)  //by definition into first size bin for each species, millions
  
  //recruitment simulation
  init_matrix recsigma(1,Nareas,1,Nspecies)    //sigma for stochastic recruitment from SR curve
  3darray rec_procError(1,Nareas,1,Nspecies,1,Nyrs)   //to generate deviations from SR curve
  
  //growth options
  //independent of bioenergetics--age based, use growthprob_phi above based on pars of age predicting length
  //4darray length(1,Nareas,1,Nspecies,1,Nages,1,Nyrs) //not needed in basic calculations so leave aside for now
  //bioenergetics based, depends on consumption--to be added
  
  //fishing mort: average annual, annual devs, actual (avg+dev)
  //**needs to be done by fleet, currently assuming each fleet has the same avg_F, F_devs and Ftot**
  //**then they sum to give F by species**
  //init_3darray avg_F(1,Nareas,1,Nspecies,1,Nfleets,avg_F_phase)  //logspace average annual fishing mort by area, species
  //init_3darray F_devs(1,Nspecies,1,Nfleets,1,Nyrs,dev_F_phase)  //logspace F deviations by area, species--NEEDS TO BE 4D, CANT DO?, FIX
  //
  //***********June 2014 replace with F = q*E formulation*****************************
  4darray Fyr(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)  //array to get annual Fs by fleet from either avg/devs or q*effort, logspace
  
  4darray suitpreybio(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins);  //suitable prey for each predator size and year, see weight unit above for units
 
  //N, B, F, Z, M2, C, need total prey consumed per pred, suitability, available prey, food habits?
  4darray N(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //numbers by area, species, size,  timestep , in millions
  4darray B(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //biomass by area, species, size,  timestep , see weight unit above for units
  4darray F(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //fishing mort by area, species, size,  timestep 
  4darray C(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //catch numbers by area, species, size, timestep , in millions
  4darray Z(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //total mort by area, species, size,  timestep 
  4darray M2(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //predation mort by area, species, size,  timestep 
  4darray eatN(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //number eaten by predators by area, species, size,  timestep 
  4darray otherDead(1,Nareas,1,Nspecies,1,Tottimesteps,1,Nsizebins) //number unknown deaths by area, species, size,  timestep
  
  //Fishery selectivities, fleet F and catch, fishery q
  4darray fishsel(1,Nareas,1,Nspecies,1,Nfleets,1,Nsizebins)  //fishery selectivity
  5darray Ffl(1,Nareas,1,Nspecies,1,Nfleets,1,Tottimesteps,1,Nsizebins) //fleet specific Fs
  5darray Cfl(1,Nareas,1,Nspecies,1,Nfleets,1,Tottimesteps,1,Nsizebins) //fleet specific Catch in numbers
  init_3darray fishery_q(1,Nareas,1,Nspecies,1,Nfleets,fqphase)
  3darray  mean_guild_fishery_q(1,Nareas,1,Nguilds,1,Nfleets) // mean q for guild and fleet// andybeet
  matrix  mean_fishery_q(1,Nareas,1,Nfleets) // mean q for fleet. ignore values of zero //andybeet
  // calculates the mean fishery_q for each guild (over fleets)
  !! for (area=1; area<=Nareas; area++) {
  !!     for (iguild=1; iguild<=Nguilds; iguild++ ) {
  !!           for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
  !!             int icount = 0;
  !!               for (spp=1; spp<=Nspecies; spp++) {
  !!                 if (guildMembers(spp)== iguild) {
  !!                    icount++;
  !!                    // sum up q's
  !!                    mean_guild_fishery_q(area,iguild,ifleet) += fishery_q(area,spp,ifleet);
  !!                 }
  !!              }
  !!                    mean_guild_fishery_q(area,iguild,ifleet) =  mean_guild_fishery_q(area,iguild,ifleet)/icount;
  !!          }
  !!      }
  !! }

  // calculates the mean q for each fleet ignoring zero q's
  !! for (area=1; area<=Nareas; area++) {
  !!           for (int ifleet=1;ifleet<=Nfleets;ifleet++) {
  !!             int icount = 0;
  !!               for (spp=1; spp<=Nspecies; spp++) {
  !!                    if (fishery_q(area,spp,ifleet) < 1e-29) {
  !!                      //ignore
  !!                    } else {
  !!                       icount++;
  !!                       mean_fishery_q(area,ifleet) += fishery_q(area,spp,ifleet);
  !!                    }
  !!              }
  !!              mean_fishery_q(area,ifleet) =  mean_fishery_q(area,ifleet)/icount;
  !!          }
  !! }
  
  //Survey qs (add selectivities)
  init_matrix survey_q(1,Nareas,1,Nspecies,sqphase)
  
  //survey obs error
  init_matrix surv_sigma(1,Nareas,1,Nspecies,ssig_phase)
  3darray surv_obsError(1,Nareas,1,Nspecies,1,Nyrs)
  
  //catch obs error
  init_3darray catch_sigma(1,Nareas,1,Nspecies,1,Nfleets,csig_phase)
  4darray catch_obsError(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs)
  
  //annual total B and SSB
  3darray avByr(1,Nareas,1,Nspecies,1,Nyrs) //uses B units--"true" simulated biomass
  3darray SSB(1,Nareas,1,Nspecies,1,Nyrs) //uses B units--"true" SSB
  
  //annual eaten B
  3darray eaten_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units
  3darray otherDead_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units
  
  //estimated fishery catch and survey catch
  3darray catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units--"true" simulated catch
  4darray fleet_catch_biomass(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs) //B units--"true" catch by fleet
  4darray est_fleet_catch_biomass(1,Nareas,1,Nspecies,1,Nfleets,1,Nyrs) //B units, can have obs error and q
  3darray est_catch_biomass(1,Nareas,1,Nspecies,1,Nyrs)  //uses B units, sums over est_fleet_catch_biomass
  4darray est_fleet_catch_guild_biomass(1,Nareas,1,Nguilds,1,Nfleets,1,Nyrs)// uses B units. sums over species. calc_catc_etx
  4darray est_fleet_catch_guild_assessment(1,Nareas,1,Nguilds,1,Nfleets,1,Nyrs)// moving average of catch over  NumAssessment years
  3darray est_catch_guild_biomass(1,Nareas,1,Nguilds,1,Nyrs) // sums over species and fllet for guild
  3darray est_survey_biomass_assessment(1,Nareas,1,Nspecies,1,Nyrs)// moving average for each species over NumAssessment years

  3darray est_survey_biomass(1,Nareas,1,Nspecies,1,Nyrs) //uses B units, can have q, obs error 
  3darray est_survey_guild_biomass(1,Nareas,1,Nguilds,1,Nyrs) //uses B units, sums over species calc_survey_abundance
  3darray est_survey_guild_biomass_assessment(1,Nareas,1,Nguilds,1,Nyrs)   // moving average of guild biomass over NumAssessment years
  
  //objective function, penalties, components of objective function (placeholder, not yet developed)
  3darray resid_catch(1,Nareas,1,Nspecies,1,Nyrs)   //log(obs)-log(est) catch
  3darray resid_bio(1,Nareas,1,Nspecies,1,Nyrs)     //log(obs)-log(est) survey bio
  matrix totcatch_fit(1,Nareas,1,Nspecies)  //fit to total catch in weight by area and species
  matrix catchcomp_fit(1,Nareas,1,Nspecies) //fit to catch at length composition
  matrix totbio_fit(1,Nareas,1,Nspecies)    //fit to total survey biomass by area and species
  matrix biocomp_fit(1,Nareas,1,Nspecies)   //fit to survey catch at length composition
  //matrix agelencomp_fit(1,Nareas,1,Nspecies) //fit to age at length composition, where available

// calc_health_indices variables AndyBeet
  matrix simpsonsIndex(1,Nareas,1,Nyrs); // index values
  3darray LFI_Biomass(1,Nareas,1,Nspecies,1,Nyrs) // large fish defined as # in top size class. For each species
  3darray LFI_Catch(1,Nareas,1,Nspecies,1,Nyrs) // large fish defined as # in top size class. For each species
  3darray LFI_Survey(1,Nareas,1,Nspecies,1,Nyrs) // large fish defined as # in top size class. For each species
  matrix LFI_threshold(1,Nareas,1,Tottimesteps) // large fish index. large fish defined as exceeding x cm (parameter read in)
  vector prob_species(1,Nspecies); //function of N
  number LF_Biomass;
  vector B_total(1,Nspecies); // total B by species at each time t
  vector B_largestClass(1,Nspecies); // biomass of largest size class
  4darray N_tot(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins); // acumulative N over the year
  4darray B_tot(1,Nareas,1,Nspecies,1,Nyrs,1,Nsizebins); // total B over year 

  matrix objfun_areaspp(1,Nareas,1,Nspecies) //sum over components for area and species
  
  objective_function_value objfun
  
  	LOCAL_CALCS
    if (debug == 2)
      {
       cout<<"rectype\n"<<rectype<<endl;
       cout<<"recruitment_alpha\n"<<recruitment_alpha<<endl;
       cout<<"recruitment_shape\n"<<recruitment_shape<<endl;
       cout<<"recruitment_beta\n"<<recruitment_beta<<endl;
      //cout<<"aAge1\n"<<aAge1<<endl<<"aFt\n"<<aFt<<endl;
      //for (i=1; i<=nsp; i++)  {
      //  cout<<"species: "<<i<<endl;
      //  cout<<"idAge1\n"<<idAge1(i)<<endl<<"idFt\n"<<idFt(i)<<endl;
      //  cout<<"iagesel\n"<<iagesel(i)<<endl<<"iFICsel\n"<<iFICsel(i)<<endl;
      //  cout<<"iYr1\n"<<iYr1(i)<<endl<<endl;
      //  }
      //cout<<"iRho\n"<<iRho<<endl;
      cout<<"\nManually exiting at the end of the parameter section...\n"<<endl;
      exit(-1);
      }
	END_CALCS


//=======================================================================================
PRELIMINARY_CALCS_SECTION
//=======================================================================================
  recruitment_alpha = rec_alpha;
  recruitment_shape = rec_shape;
  recruitment_beta = rec_beta;
  
//=======================================================================================
PROCEDURE_SECTION
//=======================================================================================

  calc_initial_states();  if (debug == 4) {cout<<"completed Initial States"<<endl;}
 
  yrct=1;
 
  
  for (t=2; t<=Tottimesteps; t++) 
     {
  //beet                // make N(t) =  N(t-1)
		//if (debug == 3) {cout<<yrct<<" "<<t<<endl;}
                // abeet
                if (t % Nstepsyr == 1) {yrct++;} // first time step in new year = > increment year

                calc_update_N(); // N(t) = N(t-1)
                // add recruits at start of year

                calc_pred_mortality(); if (debug == 4) {cout<<"completed Predation Mortality"<<endl;}
		
                calc_fishing_mortality(); if (debug == 4) {cout<<"completed Fishing Mortality"<<endl;}

                calc_total_mortality(); // We calculate Z(t) = M1 + M2 + F

		calc_catch_etc(); if (debug == 4) {cout<<"completed Catch"<<endl;} // split F among fleets

		calc_pop_dynamics(); if (debug == 4) {cout<<"completed Pop Dynamics"<<endl;} // update N - death + growth
                
            //    calc_recruitment_SSB(); if (debug == 4) {cout<<"completed Recruitment"<<endl;}
               // update N to add recruits to bin 1 at end of year.

		calc_movement(); if (debug == 4) {cout<<"completed Movement"<<endl;}
		
                calc_survey_abundance();  if (debug == 4) {cout<<"completed Survey Abundance"<<endl;}

                calc_health_indices();  if (debug == 4) {cout<<"completed Survey Abundance"<<endl;}

                calc_assessment_strategy(); if (debug == 4) {cout<<"completed calc_assessment_strategy"<<endl;}

	 }
 
  if (debug == 4) {cout<<"completed timestep loop"<<endl;}  

  evaluate_the_objective_function(); if (debug == 4) {cout<<"completed Log Likelihood"<<endl;}
  
  if (debug == 3 && sim == 0) 
    {
      cout<<"rseed\n"<<rseed<<endl;
	  cout<<"eggprod\n"<<eggprod<<endl;
	  cout<<"rectype (1=gamma/'Ricker' eggprod, 2=Deriso-Schnute SSB, 3=SSB gamma, 4=SSB Ricker, 5=SSB Beverton Holt, 9=avg+dev)\n"<<rectype<<endl;
	  cout<<"recruitment_alpha\n"<<recruitment_alpha<<endl;
      cout<<"recruitment_shape\n"<<recruitment_shape<<endl;
      cout<<"recruitment_beta\n"<<recruitment_beta<<endl;
      cout<<"stochrec\n"<<stochrec<<endl;
      cout<<"recsigma\n"<<recsigma<<endl;
      cout<<"recruitment\n"<<recruitment<<endl;
      cout<<"SSB\n"<<SSB<<endl;
      cout<<"avByr\n"<<avByr<<endl;
      cout<<"suitpreybio\n"<<suitpreybio<<endl;
      cout<<"M2\n"<<M2<<endl;
      cout<<"F\n"<<F<<endl;
      cout<<"Z\n"<<Z<<endl;
      cout<<"N\n"<<N<<endl;
      cout<<"B\n"<<B<<endl;
      cout<<"eatN\n"<<eatN<<endl;
      cout<<"otherDead\n"<<otherDead<<endl;
      cout<<"eaten_biomass\n"<<eaten_biomass<<endl;
      cout<<"otherDead_biomass\n"<<otherDead_biomass<<endl;  
      cout<<"fishsel\n"<<fishsel<<endl;
      cout<<"Ffl\n"<<Ffl<<endl;
      cout<<"C\n"<<C<<endl;
      cout<<"Cfl\n"<<Cfl<<endl;
      cout<<"fleet_catch_biomass\n"<<fleet_catch_biomass<<endl; 
      cout<<"catch_biomass\n"<<catch_biomass<<endl; 
      cout<<"est_fleet_catch_biomass\n"<<est_fleet_catch_biomass<<endl;
      cout<<"est_catch_biomass\n"<<est_catch_biomass<<endl;
      cout<<"est_survey_biomass\n"<<est_survey_biomass<<endl;
      cout<<"obs_survey_biomass\n"<<obs_survey_biomass<<endl;
      cout<<"obs_catch_biomass\n"<<obs_catch_biomass<<endl;
      cout<<"est_survey_guild_biomass_assessment\n"<< est_survey_guild_biomass_assessment<<endl;
      cout<<"B0\n"<<B0<<endl;
      cout<<"Nguilds\n"<<Nguilds<<endl;
      cout<<"guildMembers\n"<<guildMembers<<endl;
      cout<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;
      cout<<"\npin file inputs\n"<<endl;
      cout<<"yr1N\n"<<yr1N<<endl;
      cout<<"\navg_recruitment and recruitment_devs used only if rectype=9\n"<<endl;
      cout<<"avg_recruitment\n"<<avg_recruitment<<endl;
      cout<<"recruitment_devs\n"<<recruitment_devs<<endl;
      //cout<<"avg_F\n"<<avg_F<<endl;
      //cout<<"F_devs\n"<<F_devs<<endl;
      cout<<"survey_q\n"<<survey_q<<endl;
      cout<<"surv_sigma\n"<<surv_sigma<<endl;
      cout<<"catch_sigma\n"<<catch_sigma<<endl;
      cout <<"guild biomass\n"<<est_survey_guild_biomass<<endl;
      cout<<endl<<"manually exiting at end of procedure section.....\n"<<endl;
      exit(0);
    }
    
    if (debug == 3 && sim == 1) 
    { 
      write_simout();
      cout<<"rseed\n"<<rseed<<endl;
	  cout<<"rectype (1=gamma/'Ricker' eggprod, 2=Deriso-Schnute SSB, 3=SSB gamma, 4=SSB Ricker, 5=SSB Beverton Holt, 9=avg+dev)\n"<<rectype<<endl;
	  cout<<"recruitment_alpha\n"<<recruitment_alpha<<endl;
      cout<<"recruitment_shape\n"<<recruitment_shape<<endl;
      cout<<"recruitment_beta\n"<<recruitment_beta<<endl;
      cout<<"Nyrs\n"<<Nyrs<<endl;
      cout<<"Nstepsyr\n"<<Nstepsyr<<endl;
      cout<<"stochrec\n"<<stochrec<<endl;
      cout<<"recsigma\n"<<recsigma<<endl;
      cout<<"recruitment\n"<<recruitment<<endl;
      cout<<"SSB\n"<<SSB<<endl;
      cout<<"avByr\n"<<avByr<<endl;
      cout<<"M2\n"<<M2<<endl;
      cout<<"F\n"<<F<<endl;
      cout<<"Z\n"<<Z<<endl;
      cout<<"N\n"<<N<<endl;
      cout<<"eaten_biomass\n"<<eaten_biomass<<endl;
      cout<<"otherDead_biomass\n"<<otherDead_biomass<<endl;  
      cout<<"fleet_catch_biomass\n"<<fleet_catch_biomass<<endl; 
      cout<<"catch_biomass\n"<<catch_biomass<<endl; 
      cout<<"est_fleet_catch_biomass\n"<<est_fleet_catch_biomass<<endl;
      cout<<"est_fleet_catch_guild_biomass\n"<<est_fleet_catch_guild_biomass<<endl;
      cout<<"est_catch_guild_biomass\n"<<est_catch_guild_biomass<<endl;
      cout<<"est_catch_biomass\n"<<est_catch_biomass<<endl;
      cout<<"est_survey_biomass\n"<<est_survey_biomass<<endl;
      cout<<"est_survey_guild_biomass\n"<<est_survey_guild_biomass<<endl;
      cout<<"obs_survey_biomass\n"<<obs_survey_biomass<<endl;
      cout<<"obs_catch_biomass\n"<<obs_catch_biomass<<endl;
      cout<<"est_survey_guild_biomass_assessment\n"<< est_survey_guild_biomass_assessment<<endl;
      cout<<"B0\n"<<B0<<endl;
      cout<<"B0_guilds\n"<<B0_guilds<<endl;
      cout<<"Nguilds\n"<<Nguilds<<endl;
      cout<<"guildMembers\n"<<guildMembers<<endl;
      cout<<"Nthresholds\n"<<Nthresholds<<endl;
      cout<<"threshold_proportion\n"<<threshold_proportion<<endl;
      cout<<"exploitation_levels\n"<<exploitation_levels<<endl;
      cout<<"threshold_species\n"<<threshold_species<<endl;
      cout<<"AssessmentPeriod\n"<<AssessmentPeriod<<endl;
      cout<<"SpeciesDetection\n"<<speciesDetection<<endl;
      cout<<"AssessmentOn\n"<<AssessmentOn<<endl;
      cout<<"simpsonsIndex\n"<<simpsonsIndex<<endl;
      cout<<"LFI_threshold\n"<<LFI_threshold<<endl;
      cout<<"LFI_Biomass\n"<<LFI_Biomass<<endl;
      cout<<"LFI_Catch\n"<<LFI_Catch<<endl;
      cout<<"LFI_Survey\n"<<LFI_Survey<<endl;
      
      cout<<"\npin file inputs\n"<<endl;
      cout<<"yr1N\n"<<yr1N<<endl;
      //cout<<"avg_F\n"<<avg_F<<endl;
      //cout<<"F_devs\n"<<F_devs<<endl;
      cout<<"survey_q\n"<<survey_q<<endl;
      cout<<"surv_sigma\n"<<surv_sigma<<endl;
      cout<<"catch_sigma\n"<<catch_sigma<<endl;
      cout<<endl<<"manually exiting at end of procedure section....\n"<<endl;

      exit(0);
    }


//----------------------------------------------------------------------------------------
FUNCTION calc_initial_states
//----------------------------------------------------------------------------------------

  propmature.initialize();
  eggprod.initialize();
  recruitment.initialize();
  rec_procError.initialize();
  suitpreybio.initialize();
  Fyr.initialize();
  fishsel.initialize(); Ffl.initialize(); Cfl.initialize();
  N.initialize(); B.initialize(); F.initialize(); C.initialize(); 
  Z.initialize(); M2.initialize(); eatN.initialize(); otherDead.initialize();
  avByr.initialize(); SSB.initialize();
  eaten_biomass.initialize();
  otherDead_biomass.initialize();
  surv_obsError.initialize();
  catch_obsError.initialize();
  est_survey_biomass.initialize(); est_catch_biomass.initialize();
   // andy beet
  est_survey_guild_biomass.initialize();
  est_fleet_catch_guild_biomass.initialize();
  est_catch_guild_biomass.initialize();
  est_survey_guild_biomass_assessment.initialize();
  est_survey_biomass_assessment.initialize();
  est_fleet_catch_guild_assessment.initialize();
  covariates_M.initialize();
  //andybeet
  fleet_catch_biomass.initialize(); est_fleet_catch_biomass.initialize();
  catch_biomass.initialize();
  
  totcatch_fit.initialize(); catchcomp_fit.initialize();
  totbio_fit.initialize(); biocomp_fit.initialize();
  
  //need year 1 N to do recruitment, pred mort, N for following years
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
         N(area, spp, 1) = yr1N(area, spp);
         B(area, spp, 1) = wtconv*elem_prod(N(area, spp, 1),binavgwt(spp));
         avByr(area, spp, 1) = sum(B(area,spp,1))/Nstepsyr;
         est_survey_biomass(area,spp,1) = avByr(area, spp, 1);  //perfect surveys as placeholder
      }
  }
 
  
  //propmature do once for whole cov time series, estimate maturity params, covwt, or both in later phases
  //propmature = 1/(1+exp(-(maturity_nu+maturity_omega*lbinmidpt)*sum(maturity_covwt*maturity_cov)))
  // first calculate the covarate part to add. andybeet
  for (spp=1; spp<=Nspecies; spp++) {
      for (yr=1; yr<=Nyrs; yr++) {
          for (int icov=1; icov<=Nmaturity_cov;icov++) {
              covariates_M(spp,yr) += maturity_covwt(spp,icov)*maturity_cov(icov,yr);
          }
      }
  }


  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
		for(yr=1; yr<=Nyrs; yr++){
//			propmature(area, spp)(yr) = 1/(1+exp(-(maturity_nu(area, spp) + 
// Andybeet                                          maturity_omega(area, spp)*lbinmidpt(spp)) +
//                                          maturity_covwt(spp)*trans(maturity_cov)(yr)));
                   for (int isizebin=1; isizebin<=Nsizebins; isizebin++) {
			propmature(area,spp,yr,isizebin) = 1/(1+exp(-(maturity_nu(area, spp) + 
                                          maturity_omega(area, spp)*lbinmidpt(spp,isizebin)) +
                                          covariates_M(spp,yr)));
                   }
		}
    }   
  }
  
  //as long as growth is fit outside the model and transition probs done once at the beginning, could move here
 
  //fill F arrays; either start with avg F and devs by fleet, or calculate from q and effort by fleet
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
	  for(fleet=1; fleet<=Nfleets; fleet++){	
          //fill Fyr array with area, species, fleet, year specific Fs
          // Fyr(area,spp,fleet) = avg_F(area,spp,fleet) + F_devs(spp,fleet); //WARNING ONLY WORKS WITH 1 AREA:REDO 
//            Fyr(area,spp,fleet) = log(fishery_q(area,spp,fleet)*obs_effort(area,fleet)); //Andy Beet 
            Fyr(area,spp,fleet) = fishery_q(area,spp,fleet)*obs_effort(area,fleet); //Andy Beet 
      }		
    }   
  }






///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////   Random Number Generators ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

  // Add simulated observation errors for survey
    random_number_generator rng (rseed);
    dvector obsError(1,Nyrs);
    for (area=1; area<=Nareas; area++){
  	  for(spp=1; spp<=Nspecies; spp++){
         obsError.fill_randn(rng);
         surv_obsError(area,spp) = obsError;
      }   
    }
    
  // Add simulated observation errors for catch (fleet specific)
    random_number_generator rng2 (rseed+10000);
    dvector CobsError(1,Nyrs);
    for (area=1; area<=Nareas; area++){
  	  for(spp=1; spp<=Nspecies; spp++){
		 for(fleet=1; fleet<=Nfleets; fleet++){
           CobsError.fill_randn(rng2);
           catch_obsError(area,spp,fleet) = CobsError;
	     }
      }   
    }
    
    // Add simulated process errors for recruitment
    random_number_generator rng3 (rseed+20000);
    dvector RprocError(1,Nyrs);
    for (area=1; area<=Nareas; area++){
  	  for(spp=1; spp<=Nspecies; spp++){
         RprocError.fill_randn(rng3);
         rec_procError(area,spp) = RprocError;
      }   
    }

    

  if (debug == 15){
    cout<<"Ninit\n"<<N<<endl;
    cout<<"propmature\n"<<propmature<<endl;
    cout<<"Fyr\n"<<Fyr<<endl;
    cout<<"surv_obsError\n"<<surv_obsError<<endl;
    cout<<"catch_obsError\n"<<catch_obsError<<endl;
    cout<<"rec_procError\n"<<rec_procError<<endl;
    cout<<endl<<"manually exiting after calc_initial_states...\n"<<endl;
    exit(-1);
  }
  //other covariate sums here too? time series will be input



//----------------------------------------------------------------------------------------
FUNCTION calc_update_N
//----------------------------------------------------------------------------------------
// simply make N(t) = N(t-1)
 for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
            for(int isize=1; isize<=Nsizebins; isize++){
               N(area,spp,t,isize) = N(area,spp,t-1,isize);
            }
         }
  }


//----------------------------------------------------------------------------------------
FUNCTION calc_pred_mortality
//----------------------------------------------------------------------------------------
  
  //totalconsumedbypred = allmodeledprey(pred,predsize) + otherprey

  for (area=1; area<=Nareas; area++){
  	for(pred=1; pred<=Nspecies; pred++){
	    for(prey=1; prey<=Nspecies; prey++){
		  dmatrix suittemp = suitability(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins);
		  suittemp.rowshift(1); //needed to match up array bounds
	      suitpreybio(area,pred,t) += wtconv*(elem_prod(binavgwt(prey),N(area,prey,t)) * suittemp);
      }
    }
  } 
  
  //M2(area, prey, preysize,t) = sumover_preds_predsizes(intake*N(area,pred,predsize,t)*suitability(area,predpreysize)/
  //								sumover_preds_predsizes(totalconsumedbypred))
  for (area=1; area<=Nareas; area++){
  	for(prey=1; prey<=Nspecies; prey++){
	    for(pred=1; pred<=Nspecies; pred++){
		  dmatrix suittemp2 = suitability(area,pred).sub(prey*Nsizebins-(Nsizebins-1), prey*Nsizebins);
		  suittemp2.rowshift(1); //needed to match up array bounds
          M2(area,prey,t) += elem_div((elem_prod(intake(area,pred,yrct),N(area,pred,t)) * suittemp2) ,
                           (suitpreybio(area,pred,t) + 30000.0));    //Hall et al 2006 other prey too high
      }
    }
  }

//----------------------------------------------------------------------------------------
FUNCTION calc_fishing_mortality
//----------------------------------------------------------------------------------------

  //NOTE: Ftots are by area, species, and should be separated by fleet
  //not currently set up that way, assuming each fleet has same Ftot and they sum to F 
  //selectivities are not currently by area, assuming fleet selectivity same in each area
  for (area=1; area<=Nareas; area++){
      for(spp=1; spp<=Nspecies; spp++){
           for(fleet=1; fleet<=Nfleets; fleet++){
               for(int isizebin=1; isizebin<=Nsizebins; isizebin++) { //abeet added this loop to avoid compilation warnings
//  	    fishsel(area,spp,fleet) = 1/(1 + mfexp(-(fishsel_c(spp,fleet) +
//                                     (fishsel_d(spp,fleet)*lbinmidpt(spp)))));
            	  fishsel(area,spp,fleet,isizebin) = 1/(1 + mfexp(-(fishsel_c(spp,fleet) +
                                     (fishsel_d(spp,fleet)*lbinmidpt(spp,isizebin)))));
               }
               //        Ffl(area,spp,fleet,t-1) = fishsel(area,spp,fleet)*mfexp(Fyr(area,spp,fleet,yrct))/Nstepsyr; //AndyBeet
               // see also Fyr in initial_calcs
        Ffl(area,spp,fleet,t) = fishsel(area,spp,fleet)*(Fyr(area,spp,fleet,yrct))/Nstepsyr; //Andy Beet
        F(area,spp,t) += Ffl(area,spp,fleet,t);
      }
    }
  }

//----------------------------------------------------------------------------------------
FUNCTION calc_total_mortality
//----------------------------------------------------------------------------------------
 for (area=1; area<=Nareas; area++){
     for(spp=1; spp<=Nspecies; spp++){

       //mort components, with all fleets already in F
       // Split F up in calc_catch_etc
       Z(area,spp,t) = M1(area,spp) +  M2(area,spp,t) +  F(area,spp,t);
     }
 }

//----------------------------------------------------------------------------------------
FUNCTION calc_catch_etc
//----------------------------------------------------------------------------------------
 
  //calculate Catch numbers at size (C), total catch_biomass, and N/biomass eaten and dead of other causes

  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
      //temp vectors for holding proportions
      dvar_vector Fprop = elem_div(F(area,spp,t),Z(area,spp,t)); //prop of death due to fishing of each size class
      dvar_vector M2prop = elem_div(M2(area,spp,t),Z(area,spp,t)); //prop of death due to predation of each size class
      dvar_vector M1prop = elem_div(M1(area,spp),Z(area,spp,t)); // prop of death due to other mortality. M1 read in from Data file
      dvar_vector Ndeadtmp = elem_prod((1-exp(-Z(area,spp,t))),N(area,spp,t));// total number dead in each size class
      // note: Z = total mortality
      
      //these are numbers at size dying each timestep from fishing, predation, and M1 (other)
      C(area,spp,t) = elem_prod(Fprop, Ndeadtmp); 
      eatN(area,spp,t) = elem_prod(M2prop, Ndeadtmp); 
      otherDead(area,spp,t) = Ndeadtmp - C(area,spp,t) - eatN(area,spp,t);
     // C_tot(area,spp,yrct) += C(area,spp,t-1)
      //these are annual total biomass losses for comparison with production models
      eaten_biomass(area,spp,yrct) += sum(wtconv*elem_prod(eatN(area,spp,t),binavgwt(spp)));
      otherDead_biomass(area,spp,yrct) += sum(wtconv*elem_prod(otherDead(area,spp,t),binavgwt(spp)));
      
      //do fleet specific catch in numbers, biomass, sum for total catch
      for(fleet=1; fleet<=Nfleets; fleet++){
	  dvar_vector Fflprop = elem_div(Ffl(area,spp,fleet,t),F(area,spp,t));// proportion dead due to fleet in each size class. vec length= num classes
          Cfl(area,spp,fleet,t) = elem_prod(Fflprop, C(area,spp,t)); // numbers dying from fleet by sizeclass
          fleet_catch_biomass(area,spp,fleet,yrct) += sum(wtconv*elem_prod(Cfl(area,spp,fleet,t),binavgwt(spp)));
          catch_biomass(area,spp,yrct) += sum(wtconv*elem_prod(Cfl(area,spp,fleet,t),binavgwt(spp)));

          //dvar_vector ExploitationRate =  elem_prod(elem_div(Ffl(area,spp,fleet,t-1),Z(area,spp,t-1)),(1-exp(-Z(area,spp,t-1))));
          //dvar_vector roughExploitationRate = elem_div(Cfl(area,spp,fleet,t-1),N(area,spp,t-1));// catch(numbers) /population
          
          //add obs error for est fleet catch, sum for est total catch for simulations
          est_fleet_catch_biomass(area,spp,fleet,yrct) = fleet_catch_biomass(area,spp,fleet,yrct) * exp(catch_sigma(area,spp,fleet) 
                                         * catch_obsError(area,spp,fleet,yrct) 
                                         - 0.5 * catch_sigma(area,spp,fleet) * catch_sigma(area,spp,fleet)  ); //add obs error 
	      if (t % Nstepsyr == 0){//if we are just about to end the year //andybeet
                est_catch_biomass(area,spp,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
	      }//end if
      }//end fleet loop
    }//end species loop
  }//end area loop

 
 // aggregate catch to guild level at end of year Andy Beet
   if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){
   //  test<< "ttt = "<< t <<", yr =  "<<yrct <<endl;
   for(area=1; area<=Nareas; area++){
       for (int iguild=1; iguild<=Nguilds; iguild++) {
          for (spp=1; spp<=Nspecies; spp++) {
            for (fleet=1; fleet<=Nfleets; fleet++) {
                 if (guildMembers(spp) == iguild) {
                // cout<<iguild<<","<<spp<<","<<fleet<<endl;
                    est_fleet_catch_guild_biomass(area,iguild,fleet,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
                    est_catch_guild_biomass(area,iguild,yrct) += est_fleet_catch_biomass(area,spp,fleet,yrct);
                 }
            }     
          }
       }
   }
  } // end if t%

 
//----------------------------------------------------------------------------------------
FUNCTION calc_pop_dynamics
//----------------------------------------------------------------------------------------
 
  //for all older than recruits, 
  //pop is composed of survivors from previous size growing into current size and staying in area
  //plus survivors in current size not growing out of current size and staying in area
  //plus immigrants of current size from other areas
  //minus emigrants of current size to other areas
  //movement is not yet specified, so we leave out the immigration and emigration parts for now
   
  //N(area, spp,t,bin) +=sum(recruitment(area, spp,t, bin=1), 
  //                       N(area, spp,t,bin-1)*S(area,spp,t,bin-1)*growthprob_phi(area,spp,bin-1),
  //                       N(area, spp,t,bin)*S(area,spp,t,bin)*(1-growthprob_phi(area,spp,bin)))
  
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){


        // smallest size class. Survivors that stay in same size class
	N(area,spp,t,1) = N(area,spp,t-1,1)* exp(-Z(area,spp,t,1))*(1-growthprob_phi(area,spp,yrct,1));

        N_tot(area,spp,yrct,1) += N(area,spp,t,1); // running total for the year. Used for indices
      
      
      // For other size bins      
      //N = surviving and growing from smaller bin and surviving and staying in current bin                              
      for(int isize=2; isize<=Nsizebins; isize++){                              
	N(area,spp,t,isize) = N(area,spp,t-1,isize-1) * exp(-Z(area,spp,t,isize-1)) * growthprob_phi(area,spp,yrct,isize-1)
                            +  N(area,spp,t-1,isize)* exp(-Z(area,spp,t,isize))  * (1-growthprob_phi(area,spp,yrct,isize));
        N_tot(area,spp,yrct,isize) += N(area,spp,t,isize);// cumulate sum. averaged in indices

      }//end size loop
    }//end species loop
  }//end area loop


  
//----------------------------------------------------------------------------------------
FUNCTION calc_recruitment_SSB
//----------------------------------------------------------------------------------------
   
   
  //egg production(t) = sumover_length(fecundity*prop mature(t)*sexratio*N(t)) 
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
	  dvar_vector fecundmature = elem_prod(fecundity(area,spp), propmature(area, spp)(yrct));
          dvar_vector sexratioN = sexratio(area, spp) * N(area,spp)(t);  
          eggprod(area,spp)(yrct) += sum(elem_prod(fecundmature, sexratioN));  //accumulates eggs all year--appropriate?
      
          dvar_vector Nmature = elem_prod(propmature(area, spp)(yrct), N(area,spp)(t));
                //SSB(area,spp)(yrct) += sum(wtconv*elem_prod(Nmature,binavgwt(spp)));  //accumulates SSB all year; not appropriate
          SSB(area,spp)(yrct) = sum(wtconv*elem_prod(Nmature,binavgwt(spp)));  //SSB in this timestep, overwrites previous
          // Final SSB(year) = SSB in time step 5, 1- etc
    }   
  }
   
  //recruitment(t) =  recruitment_alpha  * pow (egg production(t-1),recruitment_shape) *
  //              exp(-recruitment_beta * egg production(t-1) + 
  //              sumover_?(recruitment_covwt * recruitment_cov(t)))
  if ((t % Nstepsyr == 0) && (yrct <= Nyrs)) {  // recruits enter at end of year 
    for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
		switch (rectype(spp)){
          case 1:	  				//egg production based recruitment, 3 par gamma (Ricker-ish) 
			eggprod(area,spp)(yrct) /= Nstepsyr; //average egg production for a single "spawning" timestep
			//eggprod(area,spp)(yrct) = recruitment_shape(area,spp)/recruitment_beta(area,spp);
			recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * pow(eggprod(area,spp)(yrct), recruitment_shape(area,spp)) *
                                          mfexp(-recruitment_beta(area,spp) * eggprod(area,spp)(yrct) + 
                                               recruitment_covwt(spp) * trans(recruitment_cov)(yrct));
		  break;
		  case 2:                   //SSB based recruitment, 3 par Deriso-Schnute; see Quinn & Deriso 1999 p 95
		    //SSB(area,spp)(yrct) /= Nstepsyr; //use? average spawning stock bio for a single "spawning" timestep, now SSB is at time t
    	    recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct) *
                                           pow((1-recruitment_beta(area,spp)*recruitment_shape(area,spp)*SSB(area,spp)(yrct)),
                                            (1/recruitment_shape(area,spp)));
                                     //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92       
            recruitment(area,spp)(yrct) *= mfexp(-recruitment_covwt(spp) * trans(recruitment_cov)(yrct));
		  break;
          case 3:	  				//SSB based recruitment, 3 par gamma (Ricker-ish) 
			//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
			recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * pow(SSB(area,spp)(yrct), recruitment_shape(area,spp)) *
                                          mfexp(-recruitment_beta(area,spp) * SSB(area,spp)(yrct) + 
                                               recruitment_covwt(spp) * trans(recruitment_cov)(yrct));
		  break;
          case 4:	  				//SSB based recruitment, 2 par Ricker 
			//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
			recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct) *
                                          mfexp(-recruitment_beta(area,spp) * SSB(area,spp)(yrct) + 
                                               recruitment_covwt(spp) * trans(recruitment_cov)(yrct));
		  break;
          case 5:	  				//SSB based recruitment, 2 par Beverton Holt 
			//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
			recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct) /
                                         (1 + (recruitment_beta(area,spp) * SSB(area,spp)(yrct))); 
                                     //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92       
                                      recruitment(area,spp)(yrct) *= mfexp(-recruitment_covwt(spp) * trans(recruitment_cov)(yrct));
		  break;


           case 6:
			//SSB(area,spp)(yrct) /= Nstepsyr; //average SSB for a single "spawning" timestep, now SSB is at time t
			recruitment(area,spp)(yrct) = recruitment_alpha(area,spp) * SSB(area,spp)(yrct) /
                                         (1 + (recruitment_beta(area,spp) * pow( SSB(area,spp)(yrct),recruitment_shape(area,spp)) )); 
                                     //"effective recruitment" with env covariates; see Quinn & Deriso 1999 p 92       
                                      recruitment(area,spp)(yrct) *= mfexp(-recruitment_covwt(spp) * trans(recruitment_cov)(yrct));

                 break;

           case 9:                   //Average recruitment plus devs--giving up on functional form
                       recruitment(area,spp)(yrct) = mfexp(avg_recruitment(area,spp)+recruitment_devs(area,spp,yrct));
		  break;
		  default:
            cout<<"undefined recruitment type, check .dat file"<<endl;
            exit(1);
		} //end switch

        if(stochrec(spp)){                //simulate devs around recruitment curve
           recruitment(area,spp)(yrct) *=  mfexp(recsigma(area,spp) * rec_procError(area,spp)(yrct)  
                                                  - 0.5 * recsigma(area,spp) * recsigma(area,spp));
        }  //end if stochastic

        // Now add recruitment to 1st size class
        N(area,spp,t,1) = N(area,spp,t,1) + recruitment(area,spp,yrct);
        

      }  //end spp                                            
    }  //end area
//    yrct++;  // start a new year
  }  //end if last timestep in year

            
//----------------------------------------------------------------------------------------
FUNCTION calc_movement
//----------------------------------------------------------------------------------------
  
  //not yet specified will probably have migration in array and add random too
  int probmovein = 0.0;        //will be an array for area to area movement
  int probmoveout = 0.0;       //as will this
  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
			  //N(area,spp,t) += N(area,spp,t) * probmovein(area,spp);
			 // N(area,spp,t) -= N(area,spp,t) * probmoveout(area,spp);
   // ******************************************************************************
   // weight/length is not linear. The following line this will underestimate B
      B(area,spp,t) = wtconv*elem_prod(N(area,spp,t),binavgwt(spp));  //do after movement
      for (int isize=1;isize<=Nsizebins;isize++) {
          // add up B over t for each year keeping size class structure. used in indices
          B_tot(area,spp,yrct,isize) += B(area,spp,t,isize);
      }
    }
  }
 
//----------------------------------------------------------------------------------------
FUNCTION calc_survey_abundance
//----------------------------------------------------------------------------------------
 //test<< t <<", yr =  "<<yrct <<endl;

  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
	   avByr(area,spp)(yrct) += sum(B(area,spp,t))/Nstepsyr; 
       est_survey_biomass(area,spp,yrct) =  avByr(area,spp,yrct)*survey_q(area,spp); //add surv q
       est_survey_biomass(area,spp,yrct) *= exp(surv_sigma(area,spp) 
                                         * surv_obsError(area,spp,yrct) 
                                         - 0.5 * surv_sigma(area,spp) * surv_sigma(area,spp)  ); //add obs error
              
              //    cout<<"sp "<<spp<<"  yrct= "<<yrct<<" t= "<<t<<" bio= "<< est_survey_biomass(area,spp,yrct)<<endl;
              
    }
  }
  
 // Added by Andy Beet
 // we need to sum up the biomass over each guild and check for excedences.
 // Do at end of year only. Used in assessment module
  if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){
  //  test<< "ttt = "<< t <<", yr =  "<<yrct <<endl;
   for(area = 1; area<=Nareas; area++){
       for (iguild=1; iguild<=Nguilds; iguild++) {
          for (spp=1; spp<=Nspecies; spp++) {
             if (guildMembers(spp) == iguild) {
               est_survey_guild_biomass(area,iguild,yrct) += est_survey_biomass(area,spp,yrct);
             }
          }
       }
   }
  } // end if t%


//----------------------------------------------------------------------------------------
FUNCTION calc_health_indices
//----------------------------------------------------------------------------------------
// Here we calculate several indices: measures of system health at the end of the year

 if ((t % Nstepsyr == 0) && (yrct <= Nyrs)){

// note: we could combine these indices into the same loop, but chose not to ease readability
// 1. Simpsons Diversity Index (Richness (number of species) and evenness(relative numbers))
//  we use mean N over Nstepsyr as out annual value of N
  test<< yrct<<endl;
   for(int iarea=1;iarea<=Nareas;iarea++){
      prob_species.initialize();
      dvariable N_total = 0;
      for (int isp=1; isp<=Nspecies;isp++) {
        prob_species(isp) = pow(sum(N_tot(iarea,isp,yrct))/Nstepsyr,2);
        N_total += sum(N_tot(iarea,isp,yrct))/Nstepsyr;
        test<<sum(N_tot(iarea,isp,yrct))/Nstepsyr<<endl;
      }
      simpsonsIndex(iarea,yrct) = sum(prob_species)/pow(N_total,2);
   }
// 2. Large Fish Index
//  i. LFI_Biomass = %biomass of largest sizeclass relative to total biomass for each species
// ii. LFI_Catch = same for catch data
//iii. LFI_Survey = same for survey data
//  ii. LFI_threshold = % biomass of fish > x cm relative to total biomass
   for (int iarea=1;iarea<=Nareas;iarea++) {
//       LF_Biomass = 0;
       for (int isp=1; isp<=Nspecies;isp++) {
      // cout<< B_tot(iarea,isp,yrct,Nsizebins)<<endl;
           LFI_Biomass(iarea,isp,yrct) = B_tot(iarea,isp,yrct,Nsizebins)/sum(B_tot(iarea,isp,yrct)); // large fish in top size category for each fish. Biomass
      //   LFI_Catch(iarea,isp,yrct) = C_tot(area,spp,yrct,Nsizebins)/sum(C_tot(area,spp,yrct));; // large fish in top size category for each fish. Catch
      //   LFI_Survey(iarea,isp,yrct) = S_tot(area,spp,yrct,Nsizebins)/sum(S_tot(area,spp,yrct));; // large fish in top size category for each fish. Catch
//           B_largestClass(isp) = 0.0;
//           B_total(isp) = sum(B_tot(iarea,isp,yrct));
//           for (int isize=1;isize<=Nsizebins;isize++){
//               S_total(isp) +=  est_survey_biomass()
//               if (LFI_size >= lbinmax(isp,isize)) { // small fish
//                  continue;
//               } else if (LFI_size < lbinmin(isp,isize) ) {
//                // entire interval is large fish
//                 LF_Biomass += B(iarea,isp,yrct,isize);
//               } else if ((LFI_size >= lbinmin(isp,isize)) && (LFI_size <= lbinmax(isp,isize))) {
//                 // part of interval is large fish. calculate the proportion of N then weight
//                 // this needs to be updated since length/weight is not linear
//                 dvariable propInterval = (LFI_size-lbinmin(isp,isize))/(lbinmax(isp,isize) - lbinmin(isp,isize));
//                 dvariable weight_g = propInterval*(wtbinmax(isp,isize) - wtbinmin(isp,isize))/2.0;
//                 LF_Biomass += wtconv*propInterval*(N_tot(iarea,isp,yrct,isize)/Nstepsyr)*weight_g;
//                
//               } else {
//                 test<<"Large Fish Index - issue check code"<<endl;
//                 exit(1);
//               }
//               if (isize == Nsizebins) {
//                  // in lagest size class
//                  B_largestClass(isp) +=  B(iarea,isp,t,isize);
//                  C_largestClass(isp) +=  B(iarea,isp,t,isize);
//                  S_largestClass(isp) +=   est_survey_biomass(iarea,isp,yrct);
//               }
//           }
//           LFI_Catch(iarea,isp,t) = C_largestClass(isp)/C_total(isp); // large fish in top size category for each fish
//           LFI_Survey(iarea,isp,t) = S_largestClass(isp)/S_total(isp); // large fish in top size category for each fish
      }
//       LFI_threshold(iarea,yrct) = LF_Biomass/sum(B_total); // large fish > threshold size
     }
   
   //exit(1);
  } // end of year if



//----------------------------------------------------------------------------------------
FUNCTION calc_assessment_strategy
//----------------------------------------------------------------------------------------

 if (AssessmentOn == 1) {

 // if end of year and enough years have passed to perform average.
 // every AssessmentPeriod we monitor stocks and adjust the effort for the future
 
  if ((t % Nstepsyr == 0) && (yrct <= (Nyrs-AssessmentPeriod))) {
    if( yrct % AssessmentPeriod == 0) {
       // We enter this loop every AssessmentPeriod years, at the last time period of the year.
      // now average the guild biomass values over AssessmentPeriod yrs and then we check to see if the levels exceed some threshhold
     for (area=1 ; area<=Nareas; area++) {
         for (iguild=1; iguild<=Nguilds; iguild++) {
             maxGuildThreshold(area,iguild) = Nthresholds; // set all to maximum worst case is that no change is made to effort
             for (iassess=1; iassess<=AssessmentPeriod;iassess++){
                  // calculate the mean biomass and catch over the Assessment period
                  est_survey_guild_biomass_assessment(area,iguild,yrct) += est_survey_guild_biomass(area,iguild,yrct-iassess+1)/AssessmentPeriod;
                  for (int ifleet=1;ifleet<=Nfleets;ifleet++) {// catch by fleet over last AssessmentPeriod Years
                      est_fleet_catch_guild_assessment(area,iguild,ifleet,yrct) += est_fleet_catch_guild_biomass(area,iguild,ifleet,yrct-iassess+1)/AssessmentPeriod;
                  }
             }
             // check to see if average < threshold (threshold_proportion * biomass at equilibrium)
             for (ithreshold=1; ithreshold<=Nthresholds; ithreshold++) {
             // test<<yrct<<","<<iguild<<","<<ithreshold<<endl;
                 if ((est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild)) <= threshold_proportion(ithreshold)) {
                   // test << est_survey_guild_biomass_assessment(area,iguild,yrct)/B0_guilds(area,iguild)<<endl;

                    maxGuildThreshold(area,iguild) = ithreshold;
                    // dont need to keep going for this guild since we've found the most severe case
                    break;
                  }

              }// threshold loop
         } // guild loop
     }  // area loop

     // we do the same thing but check for exceedances at the species level also
     // take the mean abundance over last AssessmentPeriod yrs for each species
     for (area=1; area<=Nareas; area++){
         for (spp=1; spp<=Nspecies; spp++){
             maxSpeciesThreshold(area,spp) = Nthresholds; // set all to safe level
             for (iassess=1; iassess<=AssessmentPeriod; iassess++){
                 // mean of last few years
                 est_survey_biomass_assessment(area,spp,yrct) +=  est_survey_biomass(area,spp,yrct-iassess+1)/AssessmentPeriod;
             }
             // now check for exceedances
             for (ithreshold=1; ithreshold<=Nthresholds; ithreshold++) {
                 if ((est_survey_biomass_assessment(area,spp,yrct)/B0(area,spp)) <= (threshold_proportion(ithreshold)+threshold_species(spp))) {
                    maxSpeciesThreshold(area,spp) = ithreshold;
                    // dont need to keep going for this species since we've found the most severe case
                    break;
                  }

             }// threshold loop
             test<<maxSpeciesThreshold(area,spp)<<", ";
         } // spp  loop
     } // area loop
     test<<endl;

     // now we have checked for exceedences we need to act on them.
     // calculate the new exploitation rate and then the new value of Effort.
     // note that if maxThreshold = Nthresholds we revert to max exploitation.
     // if any species <= threshold minimum exploitation rate.
     // First find the threshold each guild exeeded and the most severe of all combined
     for (area=1; area<=Nareas;area++) {
       int icount = 0;
        for (iguild=1;iguild<=Nguilds;iguild++) {
         //test<< maxGuildThreshold(area,iguild)<<endl;
            if (icount == 0) {
                 maxThreshold(area) =  maxGuildThreshold(area,iguild);
            } else {
                 maxThreshold(area) = min(maxThreshold(area),maxGuildThreshold(area,iguild));
            }
             icount++;
         }
         if (speciesDetection == 1) { // include species detection level in determining rate change
          for (spp=1; spp<=Nspecies; spp++){
             maxThreshold(area) = min(maxThreshold(area),maxSpeciesThreshold(area,spp));
          }
         }
         
     }
     for (area=1; area<=Nareas; area++){
         for (iguild=1; iguild<=Nguilds; iguild++) {
              test<<yrct<<","<<maxGuildThreshold(area,iguild) <<","<<maxThreshold(area)<<","<<threshold_proportion(maxThreshold(area))<<","<<exploitation_levels(maxThreshold(area))<<endl;
         }
     }

     // now we calculate the new effort for each fleet.//
     // Note that all fleets are impacted for any guild exceedance. This can and should change
     for (area=1 ; area<=Nareas ; area++) {
        for (int ifleet=1;ifleet<=Nfleets;ifleet++){
         if ( mean_fishery_q(area,ifleet) < 1e-29){ // a fleet doesn't fish a particluar guild. keep effort same
              // this will only happen at guild q not fleet q.
              effort_updated(area,ifleet) = obs_effort(area,ifleet,yrct);
         } else {
              effort_updated(area,ifleet) = exploitation_levels(maxThreshold(area))/mean_fishery_q(area,ifleet);
         }
//         test<<obs_effort(area,ifleet,yrct)<<","<<effort_updated(area,ifleet)<<","<<mean_fishery_q(area,ifleet)<<endl;
        }
     }
     // now effort is used just once in initial.calcs() to obtain the Fyr terms for the whole simulation  so
     // we need to use this new effort and create updated values for Fyr
     for (area=1; area<=Nareas ; area++) {
        for (spp=1; spp<=Nspecies; spp++) {
            for (int ifleet=1; ifleet<=Nfleets; ifleet++) {
               for (int iassess=1; iassess <= AssessmentPeriod; iassess++) { 
//                 obs_effort(area,ifleet,yrct+iassess) = value(effort_updated(area,ifleet));
                 
                Fyr(area,spp,ifleet,yrct+iassess) = fishery_q(area,spp,ifleet)*effort_updated(area,ifleet); //Andy Beet
                }
//                test<<Fyr(area,spp,ifleet,yrct+1)<<endl;
            }
        }
     }

      
   } //end if yrct
  }// end if t%

  } //end if assessment on



//----------------------------------------------------------------------------------------
FUNCTION write_simout
//----------------------------------------------------------------------------------------
  
  //send simulated biomass and catch data to csv for use in production model
      simout<<"rseed,"<<rseed<<endl;
      simout<<"BIOMASS"<<endl;
      for (area=1; area<=Nareas; area++){
   	    for(spp=1; spp<=Nspecies; spp++){
          simout<<"name_"<<spp;          
          for(yr=1; yr<=Nyrs; yr++){
             simout<<","<<est_survey_biomass(area,spp,yr);
          }
        simout<<endl;
        }
      }    
      simout<<"CATCH"<<endl;
      for (area=1; area<=Nareas; area++){
   	    for(spp=1; spp<=Nspecies; spp++){
          simout<<"name_"<<spp;
          for(yr=1; yr<=Nyrs; yr++){
              simout<<","<<est_catch_biomass(area,spp,yr);
          }
        simout<<endl;
        }
      }    

  
   
//----------------------------------------------------------------------------------------
FUNCTION evaluate_the_objective_function
//----------------------------------------------------------------------------------------
 
  //est and observed survey biomass and fishery catch are 3darrays(area,spp,yr)
  //fit matrices are area by spp
   
   resid_catch.initialize();
   resid_bio.initialize();
   totcatch_fit.initialize();
   totbio_fit.initialize();
   objfun_areaspp.initialize();

  for (area=1; area<=Nareas; area++){
  	for(spp=1; spp<=Nspecies; spp++){
		
       resid_catch(area,spp) = log(obs_catch_biomass(area,spp)+o)-log(est_catch_biomass(area,spp)+o);
       totcatch_fit(area,spp) = norm2(resid_catch(area,spp));
       
       resid_bio(area,spp) = log(obs_survey_biomass(area,spp)+o)-log(est_survey_biomass(area,spp)+o);
       totbio_fit(area,spp) = norm2(resid_bio(area,spp));
    }
  }
  //cout<<"resid_catch\n"<<resid_catch<<endl;
  //cout<<"totcatch_fit\n"<<totcatch_fit<<endl;
  //cout<<"totbio_fit\n"<<totbio_fit<<endl;
  
  objfun_areaspp = totcatch_fit + totbio_fit;
  //cout<<"objfun_areaspp\n"<<objfun_areaspp<<endl;  
  
  objfun = sum(objfun_areaspp);
       
//=======================================================================================
RUNTIME_SECTION
//=======================================================================================
  convergence_criteria 1.e-3 ,  1.e-4
  maximum_function_evaluations 1000

//=======================================================================================
TOP_OF_MAIN_SECTION
//=======================================================================================
  arrmblsize = 8000000;  //Increase amount of available dvar memory
  gradient_structure::set_CMPDIF_BUFFER_SIZE(6000000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);

//=======================================================================================
REPORT_SECTION
//=======================================================================================

  report << "EstNsize Estimated total numbers of fish " << endl;
  report << N << endl;
  report << "EstBsize Estimated total biomass of fish " << endl;
  report << B << endl;
  report << "EstRec Estimated recruitment " << endl;
  report << recruitment << endl; 
  report << "EstFsize Estimated fishing mortality " << endl;
  report << F << endl; 
  report << "EstM2size Estimated predation mortality " << endl;
  report << M2 << endl; 
  report << "EstSurvB Estimated survey biomass of fish " << endl;
  report << est_survey_biomass << endl;
  report << "ObsSurvB Observed survey biomass of fish " << endl;
  report << obs_survey_biomass << endl;
  report << "EstCatchB Estimated catch biomass of fish " << endl;
  report << est_catch_biomass << endl;
  report << "ObsCatchB Observed catch biomass of fish " << endl;
  report << obs_catch_biomass << endl;
 
