#include <cstdio>
#include <iostream>
#include <cmath>
#include "kmc_system.h"
#include "kmc_events.h"

using namespace std;

double class_events::cal_energy(int x1, int y1, int z1, int x2, int y2, int z2){ // calculate the energy of the input 2 ltcps 
	int dx= (abs(x1-x2)>= nx/2) ? nx-abs(x1-x2) :abs(x1-x2);
	int dy= (abs(y1-y2)>= ny/2) ? ny-abs(y1-y2) :abs(y1-y2);
	int dz= (abs(z1-z2)>= nz/2) ? nz-abs(z1-z2) :abs(z1-z2);
	int nn_check= dx+dy+dz;
	if(nn_check != 1 && nn_check != 3) error(2, "(cal_energy) ltcp 1 and 2 are not nn", 1, nn_check); //check; delete it
	
	int pos[2][3]= {{x1, y1, z1}, {x2, y2, z2}};

	int sum_k1=0, sum_j1=0, sum_u1=0;
	int sum_k2=0, sum_j2=0, sum_u2=0;
	for(int i=0; i<2; i ++){
		int xi= pbc(pos[i][0], nx);
		int yi= pbc(pos[i][1], ny);
		int zi= pbc(pos[i][2], nz);

		int state_i= *(states + xi*ny*nz + yi*nz + zi);
		if(0==state_i) continue;

		for(int a=0; a<n1nbr; a ++){ // 1st neighbors
			int xj= pbc(xi+(*(v1nbr+a))[0], nx);
			int yj= pbc(yi+(*(v1nbr+a))[1], ny);
			int zj= pbc(zi+(*(v1nbr+a))[2], nz);

			if((xj==x1) && (yj==y1) && (zj==z1)) continue; // avoid double count the (x1, y1, z1)-(x2, y2, z2) bond

			int state_j1= *(states + xj*ny*nz + yj*nz + zj);
			if(0==state_j1) continue;

			// calculating the energy for ABV system. If it's not ABV system, please modify it; also, e2nbr below 
			sum_k1 ++;			// i**2+j**2= 1
			sum_j1 += (state_i * state_j1); // i*j
			sum_u1 += (state_i + state_j1); // i**2*j+i*j**2= j+i
		}
				
		if(! is_e2nbr) continue;

		for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
			int xj= pbc(xi+(*(v2nbr+b))[0], nx);
			int yj= pbc(yi+(*(v2nbr+b))[1], ny);
			int zj= pbc(zi+(*(v2nbr+b))[2], nz);

			int state_j2= *(states + xj*ny*nz + yj*nz + zj);
			if(0==state_j2) continue;

			// e2nbr version of Ising model formulation
			sum_k2 ++;
			sum_j2 += (state_i * state_j2); 
			sum_u2 += (state_i + state_j2);
		}
	}

	double sum_energy = (cons_k1*sum_k1 + cons_j1*sum_j1 + cons_u1*sum_u1) + (cons_k2*sum_k2 + cons_j2*sum_j2 + cons_u2*sum_u2);
	return sum_energy;
}

double class_events::cal_energy(int* const states_ce){
	int sum_k1=0, sum_j1=0, sum_u1=0;
	int sum_k2=0, sum_j2=0, sum_u2=0;
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int xi= pbc(i, nx);
				int yi= pbc(j, ny);
				int zi= pbc(k, nz);

				int state_i= *(states_ce + xi*ny*nz + yi*nz + zi);
				if(0==state_i) continue;

				for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					int xj= pbc(xi+(*(v1nbr+a))[0], nx);
					int yj= pbc(yi+(*(v1nbr+a))[1], ny);
					int zj= pbc(zi+(*(v1nbr+a))[2], nz);

					int state_j1= *(states_ce + xj*ny*nz + yj*nz + zj);
					if(0==state_j1) continue;

					// calculating the energy for ABV system. If it's not ABV system, please modify it; also, e2nbr below 
					sum_k1 ++;			// i**2+j**2= 1
					sum_j1 += (state_i * state_j1); // i*j
					sum_u1 += (state_i + state_j1); // i**2*j+i*j**2= j+i
				}
				
				if(! is_e2nbr) continue;

				for(int b=0; b<n2nbr; b ++){ // 2nd neighbors
					int xj= pbc(xi+(*(v2nbr+b))[0], nx);
					int yj= pbc(yi+(*(v2nbr+b))[1], ny);
					int zj= pbc(zi+(*(v2nbr+b))[2], nz);

					int state_j2= *(states_ce + xj*ny*nz + yj*nz + zj);
					if(0==state_j2) continue;

					// e2nbr version of Ising model formulation
					sum_k2 ++;
					sum_j2 += (state_i * state_j2); 
					sum_u2 += (state_i + state_j2);
				}
	}}}

	double sum_energy = 0.5*((cons_k1*sum_k1 + cons_j1*sum_j1 + cons_u1*sum_u1) + (cons_k2*sum_k2 + cons_j2*sum_j2 + cons_u2*sum_u2));
	return sum_energy;
}

void class_events::rules_recb(int &typev, int &typei){
	if(typev !=0 || typei >=0)	error(2, "(rules_recb) wrong types-- ", 2, typev, typei);

//	switch(typei){
//		case -1:
//		case -2:
//		case -3:
//		case -4:
//	}
}

void class_events::rules_int_jump(int &typei, int &typea){
	if(typei >=0 || typea <=0)	error(2, "(rules_int_jump) wrong types-- ", 2, typei, typea);

//	switch(){
//		case -1:
//		case -2:
//		case -3:
//		case -4:
//	}
	
}
