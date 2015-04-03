#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "kmc_system.h"

#define MAX_NNBR 20

using namespace std;

////////// Functions for the whole program //////////
void error(int nexit, string errinfo, int nnum, double num1, double num2){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo;
	switch(nnum){
		case 0:  cout << endl;				      break;
		case 1:  cout << ": " << num1 << endl;		      break;
		case 2:  cout << ": " << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
	exit(1); 
}

void error(int nexit, string errinfo, char c[]){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo << " " << c << endl;
	exit(1); 
}

double ran_generator(){
	static bool first= true;
	if(first){
		srand(time(NULL));
		first= false;
	}
	
	return rand()/((double) RAND_MAX+1.0);
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}

////////// Functions for the whole program //////////

void class_system::ltc_constructor(){
	double (*ptr_vbra)[3]; 
	int   (*ptr_v1nbr)[3];
	int   (*ptr_v2nbr)[3];
			
	// coordinate vectors of bravais lattices
	double vbra_bcc[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
		
	// BCC: index vectors of negibhor atoms: 1st nearest
	int v1nbr_bcc[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
 		              {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
 			      { 1,  1,  1}, {-1, -1, -1}	     };
	int v2nbr_bcc[6][3]= {{ 0,  1,  1}, { 1,  0,  1}, { 0,  1,  1},
 		              { 0, -1, -1}, {-1,  0, -1}, { 0, -1, -1}};

	// Choose ltc structure
	if     (strcmp(type_ltc, "SC ")==0){}
	else if(strcmp(type_ltc, "BCC")==0){ptr_vbra= vbra_bcc;	n1nbr=8; ptr_v1nbr= v1nbr_bcc; n2nbr=6; ptr_v2nbr= v2nbr_bcc;}
	else if(strcmp(type_ltc, "FCC")==0){}
	else if(strcmp(type_ltc, "HCP")==0){}
	else	error(1, "(ltc_constructor) coldn't match the lattice type", type_ltc);
			
	// assign array values
	for(int i=0; i<3; i ++){     // Bravice vectors
		for(int j=0; j<3; j ++){
			vbra[i][j]= (*(ptr_vbra+i))[j];
	}}
	for(int i=0; i<n1nbr; i ++){ // v1nbr
		for(int j=0; j<3; j ++){
			v1nbr[i][j]= (*(ptr_v1nbr+i))[j];
	}}
	for(int i=0; i<n2nbr; i ++){ // v2nbr
		for(int j=0; j<3; j ++){
			v2nbr[i][j]= (*(ptr_v2nbr+i))[j];
	}}
}

void class_system::init_states_array(int nVset, double compA, int* const states000){
	// STATE 0: vacancy, 1: A atom, -1: B atom

	int putV= (nx*ny*nz)/nVset;
	int Vcount= 0;

	for(int i=0; i<nx*ny*nz; i++){	
		if((i%putV==0) && (Vcount<nVset)){
			Vcount ++;
			*(states000 + i)= 0;
		}
		else{
			double ran= ran_generator();

			if(ran < compA) 	*(states000 + i)= +1;
			else			*(states000 + i)= -1;
		}
	}

	(*nV)= 0; (*nA)= 0; (*nB)= 0;
	////////// CHECK //////////
	for(int i=0; i<nx*ny*nz; i++){ 
		if(*(states000+i)==  0)		(*nV) ++;
		else if(*(states000+i)== +1)	(*nA) ++;
		else if(*(states000+i)== -1)	(*nB) ++;
		else				error(1, "(init_states_array) a state type is unrecognizable", 1, *(states000+i));
	}
	if(*nV != nVset) error(1, "(init_states_array) The number of vacancies is not nVset", 2, *nV, nVset);
#define TOL 0.01
	int nAtotal= *nA + *nB;
	if(abs((double) *nA/nAtotal-compA) > TOL) error(1, "(init_states_array) the composition of generated conf is inconsistent of compA", 1, *nA);
	////////// CHECK //////////
	
	cout << "The random solution configuration has been generated!" << endl;
	cout << "Vacancy: " << *nV << endl;
	cout << "Atype A: " << *nA << ", pct: " << 100* (double) *nA/(nAtotal) << "%" << endl;
	cout << "Atype B: " << *nB << ", pct: " << 100* (double) *nB/(nAtotal) << "%" << endl;
}

void class_system::read_restart(char name_restart[], int* const states000, int &ts_initial, double &time_initial){
	ifstream if_re(name_restart, ios::in);
	if(!if_re.is_open()) error(1, "(read_restart) the file is not opened!");
	
	int timestep;
	double time;

	int ntotal;
	if_re >> ntotal;
	if(ntotal != nx*ny*nz) error(1, "(read_restart) the input total ltc number isnt consistent", 2, ntotal, nx*ny*nz);
	
	char c_ltcp[5];
	if_re >> c_ltcp >> timestep >> time;
	if(strcmp(c_ltcp, "ltcp") !=0) error(1, "(read_restart) please input a ltcp file (ltcp at the second line)"); // check
	ts_initial= timestep;
	time_initial= time;

	*nV= 0; *nA= 0; *nB= 0; *nI= 0;
	for(int index=0; index<nx*ny*nz; index ++){	
		int type, i, j, k;
		if_re >> type >> i >> j >> k;
		if(index != i*ny*nz+j*nz+k) error(1, "(read_restart) the input index inconsistent");
		
		*(states000+index)= type;
		
		if     (type== 0) (*nV) ++;
		else if(type==+1) (*nA) ++;
		else if(type==-1) (*nB) ++;
		else		  error(1, "(read_restart) an unrecognizable state", 1, type);
	}
	if(*nV+*nA+*nB+*nI != nx*ny*nz) error(1, "(read_restart) the number inconsistent", 2, *nV+*nA+*nB+*nI, nx*ny*nz);

	cout << "The configuration has been generated from the restart file!" << endl;
	cout << "Vacancy: " << *nV << endl;
	cout << "Atype A: " << *nA << ", pct: " << 100* (double)*nA / ntotal << "%" << endl;
	cout << "Atype B: " << *nB << ", pct: " << 100* (double)*nB / ntotal << "%" << endl;
	cout << "Interstitials: " << *nI << endl;
	
	if_re.close();
}

void class_system::write_conf(long long int timestep, double time, int *ptr_states){
	ofstream of_xyz;
	ofstream of_ltcp;
	if(0==timestep){
		of_xyz.open("t0.xyz");
		of_ltcp.open("t0.ltcp");
	}
	else{
		char name_xyz[40], name_ltcp[40];
		sprintf(name_xyz, "%lld", timestep);  strcat(name_xyz, ".xyz");
		sprintf(name_ltcp, "%lld", timestep); strcat(name_ltcp, ".ltcp");
		
		of_xyz.open(name_xyz);
		of_ltcp.open(name_ltcp);
	}

	if(!of_xyz.is_open()) error(1, "(write_conf) xyz file is not opened!");
	if(!of_ltcp.is_open()) error(1, "(write_conf) ltcp file is not opened!");
	
	of_xyz << nx*ny*nz << "\n" << "xyz " << timestep << " " << time << "\n";
	of_ltcp << nx*ny*nz << "\n" << "ltcp " << timestep << " " << time << "\n";
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
				double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
				double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
				
				of_xyz  << *ptr_states << " " << x << " " << y << " " << z << "\n";
				of_ltcp << *ptr_states << " " << i << " " << j << " " << k << "\n";

				ptr_states ++;
	}}}
	
	of_xyz.close();
	of_ltcp.close();
}
