#ifndef KMC_EVENTS_IMOTIONS_INCLUDED
#define KMC_EVENTS_IMOTIONS_INCLUDED
#include <iostream>
#include <cmath>
#include "kmc_system.h"
#include "kmc_events.h"

void class_events::int_motions(){
	bool IsMove;
	int i_ltcp= -1;
	do{
		IsMove= false;
		for(int i=0; i<nx; i ++){
			for(int j=0; j<ny; j ++){
				for(int k=0; k<nz; k ++){
					i_ltcp ++;
					if(*(states+i_ltcp)<0){
						if(0==int_eval(i, j, k)){
							IsMove= true;
							int_jump(i, j, k);
						}
					}
		}}}
	} while(IsMove);
}

int class_events::int_eval(int xi, int yi, int zi){
	if(*(states+xi*ny*nz+yi*nz+zi) >= 0) error(2, "(int_eval) input ltc point not a interstitial: ", 1, *(states+xi*ny*nz+yi*nz+zi)); // check

	// evaluate the status of a interstitial
	// 0: still moving(not recombined or trapped) 	
	// 1: recombined				(nI -=1)(nV -=1)
	// 2: trapped by solute neighbors	
	// 3: trapped by group of interstitials
	int vac_nbr[n1nbr];
	int count_v= 0;	// vacacny count
	int count_s= 0;	// solute atoms
	int count_i= 0; // interstitials
	for(int i=0; i<n1nbr; i ++){
		int x= pbc(xi+(*(v1nbr+i))[0], nx); 
		int y= pbc(yi+(*(v1nbr+i))[1], ny);
		int z= pbc(zi+(*(v1nbr+i))[2], nz);
	
		if     (*(states + x*ny*nz + y*nz + z)==0){ count_v ++;
			vac_nbr[count_v-1]= i;
		}
		else if(*(states + x*ny*nz + y*nz + z) >1)  count_s ++;
		else if(*(states + x*ny*nz + y*nz + z) <0)  count_i ++;
	}

	if(count_v >0){ // recombinations:    WEI-WAN-CHEN
		int ran= (int) (ran_generator()*count_v);
		int x= pbc(xi+(*(v1nbr+ran))[0], nx); 
		int y= pbc(yi+(*(v1nbr+ran))[1], ny);
		int z= pbc(zi+(*(v1nbr+ran))[2], nz);
		
		if(*(states+x*ny*nz+y*nz+z) != 0) error(2, "(int_eval) the selected ltc point is not a vacancy for recombination", 1, *(states+xi*ny*nz+yi*nz+zi)); // check
		int type_int= *(states + xi*ny*nz + yi*nz + zi);
		int type_vac= 0;

		cout << "  Rcombination: |" << type_int << "|(" << xi << " " << yi << " " << zi << ") ++ |0|(" << x << " " << y << " " << z << "), ";
		rules_recb(type_int, type_vac);		//!!!!! pass by reference, not value !!!!!
		*(states + xi*ny*nz + yi*nz + zi)= type_int;
		*(states +  x*ny*nz +  y*nz +  z)= type_vac;
		cout << "|" << type_int << "| and |" << type_vac << "|" << endl;
				
//		(*nV) --; (*nI) --; (*nAtotal) +=2;
		error(2, "(int_eval) are you sure you change the vcc_list to vector? if nV changes");
		
		return 1;
	}
//	else if(count_s >= trNsol) return 2;
//	else if(count_i >= trNint) return 3;
	else			   return 0;
}

void class_events::int_jump(int xi, int yi, int zi){
	if(*(states+xi*ny*nz+yi*nz+zi)>=0) error(2, "(int_jump) the input ltc point is not an interstitial: ", 1, *(states+xi*ny*nz+yi*nz+zi));

	int x_int= xi;
	int y_int= yi;
	int z_int= zi;

	int dir= (int) (ran_generator()*n1nbr);		// pick a number which is the index of picked direction vector
	int jump_line[3]= {*(v1nbr+dir)[0], *(v1nbr+dir)[1], *(v1nbr+dir)[2]};
	
	bool IsJump= true;
	while(IsJump){  // !!! PROBLEM: cannot jump to another interstitial!!
		int jump_sign= (int) pow(-1.0, (int) (ran_generator()*2.0));

		int x_jto= pbc(x_int + jump_sign*jump_line[0], nx);	// the position where the interstitial is about to jump into
		int y_jto= pbc(y_int + jump_sign*jump_line[1], ny);
		int z_jto= pbc(z_int + jump_sign*jump_line[2], nz);
		if(*(states+x_jto*ny*nz+y_jto*nz+z_jto)<=0) error(2, "(int_jump) the selected ltc point is not an atom: ", 1, *(states+x_jto*ny*nz+y_jto*nz+z_jto));

		int type_int= *(states + x_int*ny*nz + y_int*nz + z_int);
		int type_jto= *(states + x_jto*ny*nz + y_jto*nz + z_jto);

		cout << "  intrstl jump: |" << type_int << "|(" << x_int << " " << y_int << " " << z_int << ") -> |" << type_jto << "|(" << x_jto << " " << y_jto << " " << z_jto << "), ";
		rules_int_jump(type_int, type_jto);	//!!!!! pass by reference, not value !!!!!
		*(states + x_int*ny*nz + y_int*nz + z_int)= type_int;
		*(states + x_jto*ny*nz + y_jto*nz + z_jto)= type_jto;
		cout << "|" << type_int << "| and |" << type_jto << "|" << endl;

		x_int= x_jto;	// update the position of the interstitial
		y_int= y_jto;	
		z_int= z_jto;

		if(0 != int_eval(x_int, y_int, z_int)) IsJump= false;
	}
}
#endif // KMC_EVENTS_IMOTIONS_INCLUDED
