#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_system.h"
#include "kmc_events.h"

using namespace std;

double class_events::vac_jump(vector <double> &v_rate, vector <int> &v_ivcc, vector <int> &v_inbr){
	int states_temp[nx][ny][nz];
	for(int i=0; i<nx; i ++){ // assign values to temp. states array
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				states_temp[i][j][k]= *(states + i*ny*nz + j*nz + k);
	}}}

	double sum_rate= 0;
	if(*nV != list_vcc.size()) error(2, "(vac_jump) vacancy number inconsistent");
	for(int ivcc=0; ivcc < *nV; ivcc ++){
		int i= (int) (list_vcc.at(ivcc)/nz)/ny;
		int j= (int) (list_vcc.at(ivcc)/nz)%ny;
		int k= (int)  list_vcc.at(ivcc)%nz;

		if(states_temp[i][j][k] !=0) error(2, "(vac_jump) there's an non-vacancy in the vacancy list");
		
		int xlo=i-2, xhi=i+2;
		int ylo=j-2, yhi=j+2;
		int zlo=k-2, zhi=k+2;
		double locale_before= cal_energy(&states_temp[0][0][0], xlo, xhi, ylo, yhi, zlo, zhi);

		for(int a=0; a<n1nbr; a ++){
			int x= pbc(i+(*(v1nbr+a))[0], nx);
			int y= pbc(j+(*(v1nbr+a))[1], ny);
			int z= pbc(k+(*(v1nbr+a))[2], nz);
			
			if(*(states+ x*ny*nz + y*nz + z) != 0){
				v_ivcc.push_back(ivcc);
				v_inbr.push_back(a);
				
				double em_xyz, mu_xyz;
				if(*(states+x*ny*nz+y*nz+z)==1){ em_xyz= emA; mu_xyz= muA;} 
				else			       { em_xyz= emB; mu_xyz= muB;}

				states_temp[i][j][k]= states_temp[x][y][z];
				states_temp[x][y][z]= 0;

				double ediff= cal_energy(&states_temp[0][0][0],xlo,xhi,ylo,yhi,zlo,zhi) - locale_before; // LOCALLY
				
				if(ediff >0) v_rate.push_back(mu_xyz * exp(-beta*(em_xyz+ediff)));
				else	     v_rate.push_back(mu_xyz * exp(-beta*(em_xyz)));
							
				sum_rate += v_rate.back();
				
				states_temp[x][y][z]= states_temp[i][j][k]; // restore states_temp array to the "original" states array
				states_temp[i][j][k]= 0;
			}
		}
	}

	return sum_rate;
}

void class_events::vac_recb(int vpos[3]){ // want to change to a n-fold MC? don't use this one, use the one below
	int vx= vpos[0];
	int vy= vpos[1];
	int vz= vpos[2];
	if(*(states+vx*ny*nz+vy*nz+vz) !=0) error(2, "(vac_recb) the input ltc point is not an vacancy: ", 1, *(states+vx*ny*nz+vy*nz+vz)); // check

	int int_nbr[n1nbr];
	int icount= 0;
	for(int i=0; i<n1nbr; i ++){
		int x= pbc(vx+(*(v1nbr+i))[0], nx); 
		int y= pbc(vy+(*(v1nbr+i))[1], ny); 
		int z= pbc(vz+(*(v1nbr+i))[2], nz); 

		if(*(states + x*ny*nz + y*nz + z)<0){
			icount ++;
			int_nbr[icount-1]= i;
		}
	}

	if(icount > 0){
		int ipicked= (int) (icount*ran_generator());
		
		int x= pbc(vx+(*(v1nbr+int_nbr[ipicked]))[0], nx); 
		int y= pbc(vy+(*(v1nbr+int_nbr[ipicked]))[1], ny); 
		int z= pbc(vz+(*(v1nbr+int_nbr[ipicked]))[2], nz); 
		if(*(states+x*ny*nz+y*nz+z) >=0) error(2, "(vac_recb) the selected ltc point is not an interstitial: ", 1, *(states+x*ny*nz+y*nz+z)); // check
		
		int type_vac= 0;
		int type_int= *(states + x*ny*nz + y*nz + z);

		cout << "  Rcombination: |0|(" << vx << " " << vy << " " << vz << ") ++ |" << type_int << "|(" << x << " " << y << " " << z << "), ";
		rules_recb(type_int, type_vac);		//!!!!! pass by reference, not value !!!!!
		*(states + vx*ny*nz + vy*nz + vz)= type_vac;
		*(states + x*ny*nz  + y*nz  + z) = type_int;
		cout << "|" << type_vac << "| and |" << type_int << "|" << endl;
				
//		(*nV) --; (*nI) --; (*nAtotal) +=2;

		error(2, "(vac_recb) are you sure you change the vcc_list to vector? if nV changes");
	}
}

