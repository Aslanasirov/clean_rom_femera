#if VERB_MAX > 10
#include <iostream>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctype.h>
#include <cstring>// std::memcpy
#include "femera.h"
//
int ElastOrtho3D::Setup( Elem* E ){
  IniRot();
  JacRot( E );
  JacT  ( E );
  const uint jacs_n = E->elip_jacs.size()/E->elem_n/ 10 ;
  const uint intp_n = E->gaus_n;
  this->tens_flop = uint(E->elem_n) * intp_n
    *( uint(E->elem_conn_n)* (36+18+3) + 2*54 + 27 );
  this->tens_band = uint(E->elem_n) *(
     sizeof(FLOAT_PHYS)*(3*uint(E->elem_conn_n)*3+ jacs_n*10)
    +sizeof(INT_MESH)*uint(E->elem_conn_n) );
  this->stif_flop = uint(E->elem_n)
    * 3*uint(E->elem_conn_n) *( 3*uint(E->elem_conn_n) );
  this->stif_band = uint(E->elem_n) *(
    sizeof(FLOAT_PHYS)* 3*uint(E->elem_conn_n) *( 3*uint(E->elem_conn_n) -1+2)
    +sizeof(INT_MESH) *uint(E->elem_conn_n) );
  return 0;
}
int ElastOrtho3D::ElemLinear( Elem* E, const INT_MESH e0, const INT_MESH ee, FLOAT_SOLV *part_f, const FLOAT_SOLV* part_u ){
	//printf("Test point for elas-ort3-bas2.cc\n");
	//FIXME Cleanup local variables.
	// Xiaoyu: Define number of RHS
	const int nRHS = 6;
	const int Dm = 3;//E->mesh_d;// Node (mesh) Dimension FIXME should be elem_d?
	// Xiaoyu: define number of grains
	const int dGrain = this->node_d / (nRHS*Dm);
	// Xiaoyu: Define Dnt as Dn total. The original Dn represent DOF/node in each eigen strain problem
	const int Dnt = nRHS*Dm;
	const int Dn = 3;//this->node_d;// this->node_d DOF/node
	const int Nj = 10;//Dm*Dm+1;// Jac inv & det
	const int Nc = E->elem_conn_n;// Number of nodes/element
	const int Ne = Dn*Nc;
	const int intp_n = int(E->gaus_n);
	#if VERB_MAX>11
	printf("Dim: %i, Elems:%i, IntPts:%i, Nodes/elem:%i\n", (int)mesh_d,(int)elem_n,(int)intp_n,(int)Nc);
	#endif
	//INT_MESH   conn[Nc];
	#if 0
	#ifdef FETCH_JAC
	FLOAT_MESH* jac[Nj];
	#endif
	#endif
	FLOAT_PHYS G[Ne], u[Ne],f[Ne];
	//FLOAT_PHYS det,
	FLOAT_PHYS H[9], S[9];//, B[9];
	//
	FLOAT_PHYS intp_shpg[intp_n*Ne];
	std::copy( &E->intp_shpg[0],// local copy
	&E->intp_shpg[intp_n*Ne], intp_shpg );
	FLOAT_PHYS wgt[intp_n];
	std::copy( &E->gaus_weig[0],
	&E->gaus_weig[intp_n], wgt );
	FLOAT_PHYS C[this->mtrl_matc.size()];
	std::copy( &this->mtrl_matc[0],
	&this->mtrl_matc[this->mtrl_matc.size()], C );
	FLOAT_PHYS D[]={
		C[0],C[3],C[5],0.0,0.0,0.0,
		C[3],C[1],C[4],0.0,0.0,0.0,
		C[5],C[4],C[2],0.0,0.0,0.0,
		0.0,0.0,0.0,C[6],0.0,0.0,
		0.0,0.0,0.0,0.0,C[7],0.0,
		0.0,0.0,0.0,0.0,0.0,C[8] };
	// printf("%15.13f, %15.13f, %15.13f, %15.13f, %15.13f, %15.13f\n", C[0], C[1], C[3], C[5], C[6], C[7]);
	const FLOAT_PHYS R[9] = {
	mtrl_rotc[0],mtrl_rotc[1],mtrl_rotc[2],
	mtrl_rotc[3],mtrl_rotc[4],mtrl_rotc[5],
	mtrl_rotc[6],mtrl_rotc[7],mtrl_rotc[8]};
	#if VERB_MAX>10
	printf( "Material [%u]:", (uint)mtrl_matc.size() );
	for(uint j=0;j<mtrl_matc.size();j++){
		//if(j%mesh_d==0){printf("\n");}
		printf("%+9.2e ",C[j]);
	}; printf("\n");
	#endif
	const   INT_MESH* RESTRICT Econn = &E->elem_conn[0];
	const FLOAT_MESH* RESTRICT Ejacs = &E->elip_jacs[0];
	const FLOAT_SOLV* RESTRICT sysu  = &part_u[0];
	FLOAT_SOLV* RESTRICT sysf  = &part_f[0];
	// Xiaoyu: rotated D matrix
	double R6x6[36], RD[36], Dr[36];
	// Rotation matrix for D
	int J1, J2, I1, I2;
	for(int JJ1 = 1; JJ1<=3; JJ1++){
		if (JJ1+1 <= 3) { J1 = JJ1+1;} else {J1 = JJ1+1-3;}
		if (JJ1+2 <= 3) { J2 = JJ1+2;} else {J2 = JJ1+2-3;}
		for(int II1 = 1; II1<=3; II1++){
			if (II1+1 <= 3) { I1 = II1+1;} else {I1 = II1+1-3;}
			if (II1+2 <= 3) { I2 = II1+2;} else {I2 = II1+2-3;}
			R6x6[6*(II1-1)+JJ1-1] = R[3*(II1-1)+JJ1-1] * R[3*(II1-1)+JJ1-1];
			R6x6[6*(II1-1)+JJ1+3-1] = 2.0 * R[3*(II1-1)+J1-1] * R[3*(II1-1)+J2-1];
			R6x6[6*(II1+3-1)+JJ1-1] = R[3*(I1-1)+JJ1-1] * R[3*(I2-1)+JJ1-1];
			R6x6[6*(II1+3-1)+JJ1+3-1] = R[3*(I1-1)+J1-1] * R[3*(I2-1)+J2-1] + R[3*(I1-1)+J2-1] * R[3*(I2-1)+J1-1];
		}
	}
	// printf( "R6x6 is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n",R6x6[0],R6x6[1],R6x6[2],R6x6[3],R6x6[4],R6x6[5],R6x6[6],R6x6[7],R6x6[8],R6x6[9],R6x6[10],R6x6[11]);
	// Rotate D matrix
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++){
			RD[6*i+j] = 0.0;
			for (int k = 0; k < 6; k++) {
				RD[6*i+j] += R6x6[6*i+k]*D[6*k+j];
	}}}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++){
			Dr[6*i+j] = 0.0;
			for (int k = 0; k < 6; k++) {
				Dr[6*i+j] += RD[6*i+k]*R6x6[6*j+k];
	}}}
	// ======================= ASLAN ====================
	// int Switch[6] = {0,1,2,5,4,3};// Ort3
	int Switch[6] = {0,1,2,5,3,4};
	double Drtemp[6*6];
	for(int i=0;i<nRHS*nRHS;i++){ Drtemp[i]=Dr[i]; };
	for (int irow=0; irow<nRHS; irow++){
		for (int icol=0; icol<nRHS; icol++){
			Drtemp[irow*nRHS+icol] = Dr[Switch[irow]*nRHS+Switch[icol]];
		}
	}
	for(int i=0;i<nRHS*nRHS;i++){ Dr[i]=Drtemp[i]; };
	// ======================= ASLAN ====================

	// printf( "Rotated D (Dr) is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n",Dr[0],Dr[1],Dr[2],Dr[3],Dr[4],Dr[5],Dr[6],Dr[7],Dr[8],Dr[9],Dr[10],Dr[11]);
	// Loop over dGrain and RHSs
	for (int iGrain = 0; iGrain < dGrain; iGrain++){
		// Xiaoyu: loop over each eigen strain problem
		for (int iRHS = 0; iRHS < nRHS; iRHS++){

			if(e0<ee){// Fetch first element data
				for (int i=0; i<Nc; i++){
					//std::memcpy( &u[Dn*i], &sysu[Econn[Nc*e0+i]*Dn], sizeof(FLOAT_SOLV)*Dn );
					// Xiaoyu: extract u for iRHS
					std::memcpy( & u[Dn*i],& sysu[Econn[Nc*e0+i]*Dnt*dGrain+iGrain*nRHS*Dn+iRHS*Dn], sizeof(FLOAT_SOLV)*Dn );
				}

			}// done fetching first element

			for(INT_MESH ie=e0;ie<ee;ie++){
				//jac = &Ejacs[Nj*ie];
				
				for (int i=0; i<Nc; i++){// Fetch the current output values
					//std::memcpy(& f[Dn*i],& sysf[Econn[Nc*ie+i]*Dn], sizeof(FLOAT_SOLV)*Dn );
					// Xiaoyu: extract f for iRHS
					std::memcpy( & f[Dn*i],& sysf[Econn[Nc*ie+i]*Dnt*dGrain+iGrain*nRHS*Dn+iRHS*Dn], sizeof(FLOAT_SOLV)*Dn );
				}
				
				//const uint d = 3;
                		
				for(int ip=0; ip<intp_n; ip++){
					
					// Xiaoyu: rotate C matrix
					//G = MatMul3x3xN( jac,shg );
					//A = MatMul3xNx3T( G,u );
					for(int i=0; i<(Dm*Dm) ; i++){ H[i]=0.0;}// H[i]=0.0; B[i]=0.0; };
					//for(int i=0; i<(Ne) ; i++){ G[i]=0.0; };
					//#pragma omp simd
					for(int i=0; i<Nc; i++){
						for(int k=0; k<Dm ; k++){
							G[Dm* i+k ]=0.0;
							for(int j=0; j<Dm ; j++){
								//G[Dm* k+i ] += jac[Dm* j+i ] * intp_shpg[ip*Ne+ Dm* k+j ];
								G[Dm* i+k ] += Ejacs[Nj*ie+ Dm* j+k ] * intp_shpg[ip*Ne+ Dm* i+j ];
							}
							for(int j=0; j<Dm ; j++){
								H[Dm* k+j ] += G[Dm* i+k ] * u[Dn* i+j ];
							}
						}
					}
					// printf("G is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n",G[0],G[1],G[2], G[3],G[4],G[5], G[6],G[7],G[8], G[9],G[10],G[11]);
					// printf("intp_shpg is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n", intp_shpg[ip*Ne+0],intp_shpg[ip*Ne+1],intp_shpg[ip*Ne+2] ,intp_shpg[ip*Ne+3],intp_shpg[ip*Ne+4],intp_shpg[ip*Ne+5], intp_shpg[ip*Ne+6],intp_shpg[ip*Ne+7],intp_shpg[ip*Ne+8]);
					// printf( "Ejacs[Nj*ie+ Dm* j+i ] is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n H is %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n",Ejacs[Nj*ie+0],Ejacs[Nj*ie+1],Ejacs[Nj*ie+2], Ejacs[Nj*ie+3],Ejacs[Nj*ie+4],Ejacs[Nj*ie+5], Ejacs[Nj*ie+6],Ejacs[Nj*ie+7],Ejacs[Nj*ie+8], Ejacs[Nj*ie+9], H[0], H[1], H[2], H[3], H[4], H[5], H[6], H[7], H[8]);
					if(ip==(intp_n-1)){
						if((ie+1)<ee){// Fetch stuff for the next iteration
							const INT_MESH* RESTRICT cnxt = &Econn[Nc*(ie+1)];
							for (int i=0; i<Nc; i++){
								#if 0
								std::memcpy( &jac, &Ejacs[Nj*(ie+1)], sizeof(FLOAT_MESH)*Nj);
								#endif
								//std::memcpy(&u[Dn*i],&sysu[cnxt[i]*Dn], sizeof(FLOAT_SOLV)*Dn);
								// Xiaoyu
								std::memcpy(&u[Dn*i],&sysu[cnxt[i]*Dnt*dGrain+iGrain*nRHS*Dn+iRHS*Dn], sizeof(FLOAT_SOLV)*Dn);
							}
						}// Done fetching next iter stuff
					}
					// Compute S
					const FLOAT_PHYS dw = Ejacs[Nj*ie+ 9 ] * wgt[ip];
					// double vH[] = {H[0],H[4],H[8],H[1]+H[3],H[5]+H[7],H[2]+H[6]};

					// S[0] = (Dr[0]*vH[0] + Dr[1] *vH[1] + Dr[2]*vH[2] + Dr[3]*vH[3] + Dr[4]*vH[4] + Dr[5]*vH[5])*dw;
					// S[4] = (Dr[6]*vH[0] + Dr[7] *vH[1] + Dr[8]*vH[2] + Dr[9]*vH[3]+ Dr[10]*vH[4]+ Dr[11]*vH[5])*dw;
					// S[8] = (Dr[12]*vH[0]+ Dr[13]*vH[1]+ Dr[14]*vH[2]+ Dr[15]*vH[3]+ Dr[16]*vH[4]+ Dr[17]*vH[5])*dw;
					// S[1] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					// S[5] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					// S[2] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					// S[3]=S[1]; S[7]=S[5]; S[6]=S[2];
//==================ASLAN======================================
					double vH[] = {H[0],H[4],H[8],H[1]+H[3],H[5]+H[7],H[2]+H[6]};
					// double vH[] = {H[0],H[4],H[8],H[5]+H[7],H[2]+H[6],H[1]+H[3]};
					// double vH[] = {H[0],H[4],H[8],H[2]+H[6],H[5]+H[7],H[1]+H[3]};
					// double vH[] = {H[0],H[4],H[8],H[1]+H[3],H[2]+H[6], H[5]+H[7]};
					// double vH[] = {H[0],H[4],H[8],H[2]+H[6],H[1]+H[3], H[5]+H[7]};

					S[0] = (Dr[0]*vH[0] + Dr[1] *vH[1] + Dr[2]*vH[2] + Dr[3]*vH[3] + Dr[4]*vH[4] + Dr[5]*vH[5])*dw;
					S[4] = (Dr[6]*vH[0] + Dr[7] *vH[1] + Dr[8]*vH[2] + Dr[9]*vH[3]+ Dr[10]*vH[4]+ Dr[11]*vH[5])*dw;
					S[8] = (Dr[12]*vH[0]+ Dr[13]*vH[1]+ Dr[14]*vH[2]+ Dr[15]*vH[3]+ Dr[16]*vH[4]+ Dr[17]*vH[5])*dw;
					
					S[1] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					S[5] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					S[2] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					S[3]=S[1]; S[7]=S[5]; S[6]=S[2];

					// S[5] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					// S[2] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					// S[1] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					// S[3]=S[1]; S[7]=S[5]; S[6]=S[2];

					// S[2] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					// S[5] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					// S[1] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					// S[3]=S[1]; S[7]=S[5]; S[6]=S[2];

					// S[1] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					// S[2] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					// S[5] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					// S[3]=S[1]; S[7]=S[5]; S[6]=S[2];

					// S[2] = (Dr[18]*vH[0]+ Dr[19]*vH[1]+ Dr[20]*vH[2]+ Dr[21]*vH[3]+ Dr[22]*vH[4]+ Dr[23]*vH[5])*dw;
					// S[1] = (Dr[24]*vH[0]+ Dr[25]*vH[1]+ Dr[26]*vH[2]+ Dr[27]*vH[3]+ Dr[28]*vH[4]+ Dr[29]*vH[5])*dw;
					// S[5] = (Dr[30]*vH[0]+ Dr[31]*vH[1]+ Dr[32]*vH[2]+ Dr[33]*vH[3]+ Dr[34]*vH[4]+ Dr[35]*vH[5])*dw;
					// S[3]=S[1]; S[7]=S[5]; S[6]=S[2];

					// Compute f
					for(int i=0; i<Nc; i++){
						for(int k=0; k<Dn; k++){
							for(int j=0; j<Dn; j++){
								f[Dn* i+k ] += G[Dm* i+j ] * S[Dn* j+k ];// 18*N FMA FLOP
								// f[Dn* i+k ] += G[Dm* i+j ] * S[Dn* k+j ];// 18*N FMA FLOP
					} } }
					// printf( "S is: %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n f is : %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e %16.14e, %16.14e, %16.14e, %16.14e, %16.14e, %16.14e \n",S[0],S[1],S[2], S[3],S[4],S[5], S[6],S[7],S[8], f[0],f[1],f[2], f[3],f[4],f[5], f[6],f[7],f[8], f[9],f[10],f[11]);
				}//end intp loop
				#if 1
				for (int i=0; i<Nc; i++){// Write output back to system vector
					//std::memcpy(& sysf[Econn[Nc*ie+i]*Dn],& f[Dn*i], sizeof(FLOAT_SOLV)*Dn );
					// Xiaoyu
					std::memcpy(& sysf[Econn[Nc*ie+i]*Dnt*dGrain+iGrain*nRHS*Dn+iRHS*Dn],& f[Dn*i], sizeof(FLOAT_SOLV)*Dn );
				}
				#else
				for (int i=0; i<Nc; i++){
					for(int j=0; j<3; j++){
						//part_f[3*Econn[Nc*ie+i]+j] += f[(3*i+j)];
						sysf[3*conn[i]+j] += f[3* i+j ];
				} }//---------------------------------------------------- N*3 =  3*N FLOP
				#endif
			}//end elem loop
		}// end RHS loop
	}
	return 0;
}
