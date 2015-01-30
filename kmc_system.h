#ifndef KMC_SYSTEM_INCLUDED
#define KMC_SYSTEM_INCLUDED

#include <iostream>
#include <cstring>

#define MAX_NNBR 20
#define MAX_NSPE 20

using namespace std;

////////// Functions for the whole program //////////
void error(int nexit, string errinfo, int nnum=0, double num1=0, double num2=0);
double ran_generator();
int pbc(int x_, int nx_); // Periodic Boundary Condition
////////// Functions for the whole program //////////

class class_system{
	public:
		int n1nbr, n2nbr;	// number of neighbors
		int v1nbr[MAX_NNBR][3];	// indexes vectors of 1st neighbors
		int v2nbr[MAX_NNBR][3];	// indexes vectors of 2nd neighbors
		
		class_system(int nx_, int ny_, int nz_, const char type_ltc_[4], int *nA_, int *nB_, int *nV_, int *nI_): 
		nx(nx_), ny(ny_), nz(nz_), 
		nA(nA_), nB(nB_), nV(nV_), nI(nI_) // pointers
		{
			stpcpy(type_ltc, type_ltc_);

			ltc_constructor();	// execute constructor for vbra, n*nbr, and v*nbr
			
			// Printf the parameters
			cout << "Setting: " << endl << "nx= " << nx << ", ny= " << ny << ", nz= " << nz << endl;
			cout << "The crystal structure is: " << type_ltc << endl;
			for(int i=0; i<3; i++)
				cout << vbra[i][0] << " " << vbra[i][1] << " " << vbra[i][2] << endl;
			cout << "And the number of 1st neighbors: " << n1nbr << endl;
			cout << "    the number of 2nd neighbors: " << n2nbr << endl;
		}
		
		// Functions
		void init_states_array(int nVset, double compA, int* const states000);
		void read_restart(char name_restart[], int* const states000, int &ts_initial, double &time_initial);
		void write_conf(int timestep, double time, int *ptr_states);

	private:
		int       nx, ny, nz;
		char      type_ltc[4];	// crystal structure type
		double    vbra[3][3];	// coordinate vectors of bravice lattice

		int *nA;		// Pointer variables: all classes and functions use the same 
		int *nB;
		int *nV;	
		int *nI;	
		
		// functions
		void ltc_constructor();
};
#endif // KMC_SYSTEM_INCLUDED
