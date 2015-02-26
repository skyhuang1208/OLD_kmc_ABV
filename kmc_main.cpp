#include <iostream>
#include <fstream>
#include <ctime>
#include "kmc_par.h"
#include "kmc_system.h"
#include "kmc_events.h"
using namespace std;

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "########## Initializing System... ##########" << endl;
	
	int    ts_bg;
	double time_bg;
	
	int nA= 0;
	int nB= 0;
	int nV= 0;
	int nI= 0;
	
	class_system sys(par_nx, par_ny, par_nz, par_ltc, &nA, &nB, &nV, &nI);

	cout << "\n## Creating STATES array... ##" << endl;
	int states[par_nx][par_ny][par_nz];
	
	if(par_isrestart){
		cout << "RESTART FROM restart file..." << endl;
		if(nArg != 2) error(0, "when restart flag is true, nArg must be 2");
		sys.read_restart(Arg[1], &states[0][0][0], ts_bg, time_bg);
	}
	else{
		cout << "START FROM a random configuration..." << endl;
		ts_bg= 0; time_bg= 0;
		sys.init_states_array(par_nV, par_compA, &states[0][0][0]);
		sys.write_conf(ts_bg, time_bg, &states[0][0][0]);
		cout << "Output t0 conf files" << endl;
	}

	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events(par_step_write_his, par_nx, par_ny, par_nz, &nA, &nB, &nV, &nI, &states[0][0][0], 
			    sys.n1nbr, sys.v1nbr, sys.n2nbr, sys.v2nbr, par_name_sol, par_name_vcc, par_isrestart); 
	events.input_par(par_beta, par_muA, par_muB, par_emA, par_emB, 
			 par_e1AA, par_e1BB, par_e1AB, par_e1AV, par_e1BV, par_e1VV,
			 par_e2AA, par_e2BB, par_e2AB, par_e2AV, par_e2BV, par_e2VV);

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	long long int timestep=  ts_bg; 
	double totaltime= time_bg;
	int STEP_LOG= par_tstep/2000; cout << "TIMESTEP() TIME(s) ETOTAL(eV)";
	while((totaltime<= time_bg+par_tend) && (timestep != ts_bg+par_tstep)){
		// CALCULATION
		events.events_main(timestep, totaltime);

		// OUTPUT DATA
		if(0==timestep%STEP_LOG)
			cout << endl << timestep << " " << totaltime << " " << events.cal_energy(&states[0][0][0]) << " ";
		
		if(0==timestep%par_confts){
			sys.write_conf(timestep, totaltime, &states[0][0][0]);
			cout << "<Output conf files at: " << timestep << ">";
		}
	}

	// finalizing
	if(timestep%STEP_LOG != 0)   
		cout << endl << timestep << " " << totaltime << " " << events.cal_energy(&states[0][0][0]) << endl;
	if(timestep%par_confts != 0){ 
		sys.write_conf(timestep, totaltime, &states[0][0][0]); 
		cout << "Output conf files at: " << timestep;
	}

	int tfcpu= time(0);
	cout << "\n**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
