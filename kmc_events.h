#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class class_events{
	public:
		class_events(int nx_, int ny_, int nz_, int *nA_, int *nB_, int *nV_, int *nI_, int *states_, 
			     int n1nbr_, int v1nbr_[][3], int n2nbr_, int v2nbr_[][3], const char name_sol[], const char name_vcc[], bool isrestart):
		nx(nx_), ny(ny_), nz(nz_), nA(nA_), nB(nB_), nV(nV_), nI(nI_),
		states(states_), n1nbr(n1nbr_), v1nbr(v1nbr_), n2nbr(n2nbr_), v2nbr(v2nbr_)
		{
			init_list_vcc();

			if(isrestart){
				his_sol= fopen(name_sol, "a");
				his_vcc= fopen(name_vcc, "a");
			}
			else{
				his_sol= fopen(name_sol, "w");
				his_vcc= fopen(name_vcc, "w");
			}
			if(NULL==his_sol) error(2, "(class_events) the solute  history file was not opened!");
			if(NULL==his_vcc) error(2, "(class_events) the vacancy history file was not opened!");
	
			actions_sol[0].reserve(1000); 
			actions_sol[1].reserve(1000);
		}
		
		// functions
		void events_main(int& timestep, double& totaltime);
		void input_par(double beta_, double muA_, double muB_, double emA_, double emB_,
			       double cons_k1_, double cons_j1_, double cons_u1_, double cons_k2_, double cons_j2_, double cons_u2_);
		double cal_energy(int* const states_ce);
	private:
		int nx, ny, nz;		// system size 
		int *nA, *nB, *nV, *nI;	// # of atoms, vacancies, instls; pass by pointers

		FILE * his_sol;		// history file of solute atoms
		FILE * his_vcc;		// history file of vacancy and time: record every several steps
		int* const states;	// Don't change the address of this!
		vector <int> list_vcc;	// A list contains indexs of all vacancies
		vector <int> list_int;	// A list contains indexs of all interstitials
		
		vector <int> actions_sol[2]; // A list contains solute atom moves from [0] to [1]
		
		int n1nbr, n2nbr;	// number of neighbors
		int (*v1nbr)[3];	// indexes vectors of 1st neighbors
		int (*v2nbr)[3];	// indexes vectors of 2nd neighbors
		
		double beta; 	 	// energy parameters 
		double muA, muB;
		double emA, emB;
		double cons_k1, cons_j1, cons_u1; // Ising model
		double cons_k2, cons_j2, cons_u2; // Ising model
		bool is_e2nbr;

		//////functions for all events//////
		void write_hissol(int timestep, double totaltime, const vector<int> (&actions_sol)[2]);
		void write_hisvcc(int timestep, double totaltime);
		double cal_energy(int x1, int y1, int z1, int x2, int y2, int z2); 
		void init_list_vcc();

		//////functions for vacancies //////
		double vac_jump(vector <double> &v_rate, vector <int> &v_ivcc, vector <int> &v_inbr);
		void vac_recb(int vpos[3]);

		//////functions for interstitials////
		void int_motions();
		int  int_eval(int x_int, int y_int, int z_int);
		void int_jump(int x_begin, int y_begin, int z_begin);
		
		/////sorts of rules /////
		void rules_int_jump(int &typei, int &typef);	// initial to final
		void     rules_recb(int &typei, int &typev);	// interstitial and vacancy
};


#endif // KMC_EVENTS_INCLUDED
