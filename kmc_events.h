#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

class class_events{
	public:
		class_events(int nx_, int ny_, int nz_, int *nA_, int *nB_, int *nV_, int *nI_, int *states_, 
			     int n1nbr_, int v1nbr_[][3], int n2nbr_, int v2nbr_[][3], const char hisname[], bool isrestart):
		nx(nx_), ny(ny_), nz(nz_), 
		nA(nA_), nB(nB_), nV(nV_), nI(nI_),
		states(states_),
		n1nbr(n1nbr_), v1nbr(v1nbr_), n2nbr(n2nbr_), v2nbr(v2nbr_)
		{
			init_list_vcc();
			if(isrestart) 	out_his= fopen(hisname, "a");
			else		out_his= fopen(hisname, "w");
			if(NULL==out_his) error(2, "(class_events) the history file was not opened!");
		}
		
		// functions
		void events_main(int& timestep, double& totaltime);
		void input_par(double beta_, double muA_, double muB_, double emA_, double emB_,
			       double cons_k1_, double cons_j1_, double cons_u1_, double cons_k2_, double cons_j2_, double cons_u2_);
	private:
		int nx, ny, nz;		// system size 
		int *nA, *nB, *nV, *nI;	// # of atoms, vacancies, instls; pass by pointers

		FILE * out_his;		// OFstream of the history file
		int *states;		// Don't change the address of this!
		vector <int> list_vcc;	// A list contains indexs of all vacancies
		vector <int> list_int;	// A list contains indexs of all interstitials

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
		void write_his(int state, int i_ltcp);
		void write_his_time(double dt);
		double cal_energy(int *states_ce, int xlo, int xhi, int ylo, int yhi, int zlo, int zhi); // Ising model
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
