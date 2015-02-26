#ifndef KMC_PAR_INCLUDED
#define KMC_PAR_INCLUDED

const char   par_ltc[4]=		"BCC";

const int    par_nx=                       64;
const int    par_ny=                       64;
const int    par_nz=                       64;

const double par_compA=                  0.97; // composition of A atoms
const int    par_nV=                        1;

const double par_tend=                   1e13; 	// toal time (s)
const long long int par_tstep=           1e10; 	// toal timestep (give a minus step to ignore this quiterior to end the simulation)
const long long int par_confts=   par_tstep/5;	// timestep that output a conf file for restart later
const long long int par_step_write_his=	  1e4;

const bool par_isrestart=		    false;

const char   par_name_sol[20]=      "history.sol";
const char   par_name_vcc[20]=      "history.vcc";

// Ising model energy calculation parameters
const double par_temp=                       773.0; // 500c
const double par_beta= 1.0/par_temp/8.617332478e-5; // 1/kbT, units: eV, K

const double par_muA=			   6.1e+12; // units: s^-1 
const double par_muB=			   6.1e+12; // units: s^-1 

const double par_emA=				 0;
const double par_emB=				 0;

// bonding energy parameters
const double par_e1AA=                           0;
const double par_e1BB=                           0;
const double par_e1AB=                     0.04732;
const double par_e1AV=                     0.04732;
const double par_e1BV=                           0;
const double par_e1VV=                           0;

const double par_e2AA=                           0;
const double par_e2BB=                           0;
const double par_e2AB=                     0.01997;
const double par_e2AV=                     0.01997;
const double par_e2BV=                           0;
const double par_e2VV=                           0;

// trapping number: solute atom trapping and intersitial trapping							 
// 	const int par_trNsol= 3; // trapping number by solute atoms
// 	const int par_trNint= 3; // trapping number by interstitials
#endif // KMC_PAR_INCLUDED
