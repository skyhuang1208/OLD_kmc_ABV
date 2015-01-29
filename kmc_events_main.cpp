#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include "kmc_system.h"
#include "kmc_events.h"

using namespace std;

void class_events::events_main(int& timestep, double& totaltime){
	vector <double> v_rate; v_rate.reserve(*nV*8); // vacancy: rate(or probability)
	vector <int>    v_ivcc; v_ivcc.reserve(*nV*8); // vacancy: id in the vacancy list 
	vector <int>    v_inbr; v_inbr.reserve(*nV*8); // vacancy: id of the neighbor vector
	double vrate= vac_jump(v_rate, v_ivcc, v_inbr);

	double sum_rate= vrate; // sum_rate= v_rate + i_rate after we have interstitial jump events
	double ran= ran_generator();
	double acc_rate= 0; // accumulated rate

	for(int i=0; i<v_rate.size(); i ++){
		if( (ran > acc_rate) && (ran <= (acc_rate + v_rate.at(i)/sum_rate) ) ){
			int vx= (int) (list_vcc.at(v_ivcc.at(i))/nz)/ny;
			int vy= (int) (list_vcc.at(v_ivcc.at(i))/nz)%ny;
			int vz= (int)  list_vcc.at(v_ivcc.at(i))%nz;
			
			int x= pbc(vx+(*(v1nbr+v_inbr.at(i)))[0], nx);
			int y= pbc(vy+(*(v1nbr+v_inbr.at(i)))[1], ny);
			int z= pbc(vz+(*(v1nbr+v_inbr.at(i)))[2], nz);

			*(states + vx*ny*nz + vy*nz + vz)= *(states + x*ny*nz + y*nz + z); write_his(*(states+vx*ny*nz+vy*nz+vz), vx*ny*nz+vy*nz+vz);
			*(states +  x*ny*nz +  y*nz +  z)= 0;				   write_his(*(states+x*ny*nz+y*nz+z),     x*ny*nz+ y*nz+ z);

			list_vcc.at(v_ivcc.at(i))= x*ny*nz + y*nz + z;
		}
		
		acc_rate += v_rate.at(i)/sum_rate;
	}
	
	if(abs(acc_rate-1.0)>0.00000001) error(2, "(vac_jump_void) acc_rate isnt 100p at the end", 1, acc_rate); // if dont need this, put acc_rate += ... into an else
	if(*nA+*nB+*nV+*nI != nx*ny*nz)  error(2, "(main) numbers of ltc points arent consistent", 2, *nA+*nB+*nV+*nI, nx*ny*nz); // check

	double dt= 1.0/sum_rate;
	timestep ++;
	totaltime += dt;
	write_his_time(dt);
}

// functions in backupfun:
//	int vpos[3];
//	events.vac_jump_random(par_pr_vjump, vpos);
//	events.vac_recb(vpos);
//	events.int_motions();

void class_events::write_his(int state, int i_ltcp){
	fprintf(out_his, "%d %d\n", state, i_ltcp);
}

void class_events::write_his_time(double dt){
	fprintf(out_his, "dt: %e\n", dt);
}

void class_events::input_par(double beta_, double muA_, double muB_, double emA_, double emB_, 
			     double cons_k1_, double cons_j1_, double cons_u1_, double cons_k2_, double cons_j2_, double cons_u2_){
	beta= beta_;

	muA= muA_; muB= muB_;
	emA= emA_; emB= emB_;

	cons_k1= cons_k1_; cons_j1= cons_j1_; cons_u1= cons_u1_;
	cons_k2= cons_k2_; cons_j2= cons_j2_; cons_u2= cons_u2_;

	if((0==cons_k2) && (0==cons_j2) && (0==cons_u2)) is_e2nbr= false;
	else						 is_e2nbr= true;

	// print out the parameters to log file
	cout << "Parameters:" << endl; 
	
	cout << "beta= " << beta << endl;
	printf("mu= %f %f\n", muA, muB);
	printf("Em= %f %f\n", emA, emB);
	
	cout << "Ising formulation constants:" << endl;
	printf("Constant K: %f %f\n", cons_k1, cons_k2);
	printf("Constant J: %f %f\n", cons_j1, cons_j2);
	printf("Constant U: %f %f\n", cons_u1, cons_u2);
}

void class_events::init_list_vcc(){
	int n_check= 0;

	for(int i_ltcp=0; i_ltcp<nx*ny*nz; i_ltcp ++){
		if(0==*(states+i_ltcp)){
			n_check ++;
			list_vcc.push_back(i_ltcp);
		}
	}

	if(n_check != *nV) error(2, "(init_list_vcc) Vacancy number inconsistent", 2, n_check, *nV);
}
