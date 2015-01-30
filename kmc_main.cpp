#include <iostream>
#include <fstream>
#include <ctime>
#include "kmc_par.h"
#include "kmc_system.h"
#include "kmc_events.h"
using namespace std;

#define STEP_LOG 10000

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
	class_events events(par_nx, par_ny, par_nz, &nA, &nB, &nV, &nI, &states[0][0][0], sys.n1nbr, sys.v1nbr, sys.n2nbr, sys.v2nbr, par_hisname, par_isrestart); 
	events.input_par(par_beta, par_muA, par_muB, par_emA, par_emB, par_consk1, par_consj1, par_consu1, par_consk2, par_consj2, par_consu2);

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	int    timestep=  ts_bg; 
	double totaltime= time_bg;
	cout << "TIMESTEP() TIME(s) ETOTAL(eV)" << endl;
	while((totaltime<= time_bg+par_tend) && (timestep != ts_bg+par_tstep)){
		// CALCULATION
		events.events_main(timestep, totaltime);

		// OUTPUT DATA
		if(0==timestep%STEP_LOG){
			cout << timestep << totaltime << sys.cal_energy(&states[0][0][0]) << endl;
		}
		if(0==timestep%par_confts){
			sys.write_conf(timestep, totaltime, &states[0][0][0]);
			cout << "Output conf files at: " << timestep << endl << endl;
		}
	}
	sys.write_conf(timestep, totaltime, &states[0][0][0]); cout << "Output conf files at: " << timestep << endl;

	int tfcpu= time(0);
	cout << "**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
