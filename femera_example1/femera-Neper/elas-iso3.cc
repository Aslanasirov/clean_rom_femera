#if VERB_MAX > 10
#include <iostream>
#endif
#include <cstring>// std::memcpy
#include "femera.h"
//
int ElastIso3D::ElemLinear( Elem* ){ return 1; }//FIXME
int ElastIso3D::ElemJacobi( Elem* ){ return 1; }//FIXME
int ElastIso3D::BlocLinear( Elem* ,
  RESTRICT Phys::vals &, const RESTRICT Solv::vals & ){
  return 1;
  }
// Xiaoyu: print displacement
//int ElastIso3D::ElemStrainStress(std::ostream& , Elem* , FLOAT_SOLV* ){
//  return 1;
//}
int ElastIso3D::ElemStrainStress(std::ostream& of, Elem* E, FLOAT_SOLV* part_u){
	
	// Xiaoyu's modification
	printf("Test of ISO3 ElemStrainStress \n");

	const INT_MESH elem_n =E->elem_n;
	const int intp_n = int(E->gaus_n);
	// Xiaoyu: define number of RHS
	const int nRHS = 6;
	const int Dn = 3;// this->node_d DOF/node for one RHS
	// Xiaoyu: define number of grains
	const int dGrain = this->node_d / (nRHS*Dn);
	const int Nc = E->elem_conn_n;// Number of nodes/element
	const int Ne = Dn*Nc;
	FLOAT_PHYS u[Ne];
	//FLOAT_PHYS f[Ne];
	const   INT_MESH* RESTRICT Econn = &E->elem_conn[0];
	const FLOAT_SOLV* RESTRICT sysu  = &part_u[0];
	//const FLOAT_SOLV* RESTRICT sysf  = &part_f[0];
	INT_MESH e0=0, ee=elem_n;
	of << "ee is " << ee;
	of << '\n';
//	of << "Start to print nodal displacement!";
//	of << '\n';
	for(INT_MESH ie=e0;ie<ee;ie++){//=============================== Element Loop
		for(int ip=0; ip<intp_n; ip++){//=================== Integration Point Loop
			// Xiaoyu: print element id and integration point
			of << E->elem_glid[ie] <<","<< ip+1;
			of << '\n';
			// Xiaoyu: loop over grains
			for (int iGrain=0; iGrain<dGrain; iGrain++){
				// Xiaoyu: loop over each RHS
				for (int iRHS=0; iRHS<nRHS; iRHS++){
					// print iRHS
					//of << iRHS;
					//of << '\n';
					for (int i=0; i<Nc; i++){
						std::memcpy(&u[Dn*i],&sysu[Econn[Nc*ie+i]*Dn*nRHS*dGrain+iGrain*nRHS*Dn+iRHS*Dn],sizeof(FLOAT_SOLV)*Dn);
					}
					// Xiaoyu: print node displacement - Tet element, 4 nodes, 12 displacement
					of << u[0] << "," << u[1] << "," << u[2]  << "," << u[3] << ",";
					of << u[4] << "," << u[5] << "," << u[6]  << "," << u[7] << ",";
					of << u[8] << "," << u[9] << "," << u[10] << "," << u[11];
					of << '\n';
					
				}
			}
			//of << '\n';
			//fflush(stdout);
		}
	}
	
  return 1;
}
// Xiaoyu
int ElastIso3D::PijklCalculation(std::ostream& of, Elem* E, FLOAT_SOLV* part_u, int part_n){
//	printf("Test of ElastIso3D::PijklCalculation \n");
	// Define constants
	const int d = 3, d2=d*d, vn = 4;
	const int nRHS = 6;
	const int Dn = 3;
	const int dGrain = this->node_d / (nRHS*Dn);
	const INT_MESH elem_n =E->elem_n;
	const   INT_MESH* RESTRICT Econn = &E->elem_conn[0];
	const FLOAT_SOLV* RESTRICT sysu  = &part_u[0];
	FLOAT_PHYS B[12*6];
	FLOAT_PHYS u[Dn*vn];
	FLOAT_PHYS PhDisp[Dn*vn*nRHS*dGrain];
	FLOAT_PHYS PhPolar[6*6*dGrain];
	FLOAT_PHYS PhPolarOut[6*6*dGrain];
	FLOAT_PHYS phSegsPVec[6*6*dGrain*elem_n];
	FLOAT_PHYS Pijkl[6*6*part_n];// Portion of Pijkl for the current part
	int Switch[6] = {0,1,2,3,5,4};
	const FLOAT_MESH* RESTRICT Ejacs = &E->elip_jacs[0];
	for(int i=0;i<6*6*dGrain*int(elem_n);i++){ phSegsPVec[i]=0.0; };
	double PartVolume = 0.0;// Partition volume
	double weights[elem_n];// weight of phSegsPVec
//	of << "elem_n is " << elem_n << '\n';
	// Loop over elements
	INT_MESH e0=0, ee=elem_n;
	for(INT_MESH ie=e0;ie<ee;ie++){
		double ElVol = 0.0;
		// -------------- Step 1: Reconstruct B_matrix ---------------
		// compute vert
		FLOAT_MESH vert[vn*d];
		for(uint i=0; i<vn; i++){
			uint n=E->elem_conn[E->elem_conn_n*ie + i];
			for(uint j=0;j<d;j++){
				vert[d*i + j]=E->node_coor[d*n + j];
			}
		}
		// Compute jac
		FLOAT_MESH jac[d2];
		for(uint i=0;i< d;i++){
		      for(uint k=0;k< d;k++){
		    	  jac[d* i+k ]=0.0;
		    	  for(uint j=0;j<vn;j++){
		    		  jac[d* i+k ] += E->intp_shpg[d*j + i] * vert[d*j + k ];
		    	  }
		      }
		}
		// Compute corrected G
		FLOAT_PHYS Gc[12];
		const int Nj = 10;
		for(uint i=0;i<12;i++){ Gc[i]=0.0; };// initialize G for each integration point
		for(uint k=0;k<4;k++){
			for(uint i=0;i<3;i++){
				for(uint j=0;j<3;j++){
					Gc[4* i+k] += E->elip_jacs[Nj*ie+3* j+i] * E->intp_shpg[3* k+j];
		} } }
		// reconstruct B and D
		for(uint j=0; j<(12*6); j++){ B[j]=0.0; }
		for(uint j=0; j<4; j++){
			B[12*0 + 0+j*3] = Gc[4*0+j];// xx yy zz
			B[12*1 + 1+j*3] = Gc[4*1+j];
			B[12*2 + 2+j*3] = Gc[4*2+j];
			B[12*3 + 0+j*3] = Gc[4*1+j];// xy yx
			B[12*3 + 1+j*3] = Gc[4*0+j];
			B[12*4 + 1+j*3] = Gc[4*2+j];// yz zy
			B[12*4 + 2+j*3] = Gc[4*1+j];
			B[12*5 + 0+j*3] = Gc[4*2+j];// xz zx
			B[12*5 + 2+j*3] = Gc[4*0+j];
		}
		// -------------- Step 2: Extract PhDisp ---------------
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (int iRHS=0; iRHS<nRHS; iRHS++){
				for (int i=0; i<vn; i++){
					std::memcpy(&u[Dn*i],&sysu[Econn[vn*ie+i]*Dn*nRHS*dGrain+iGrain*nRHS*Dn+iRHS*Dn],sizeof(FLOAT_SOLV)*Dn);
				}
				for (int i=0; i<Dn*vn; i++){
					PhDisp[iGrain*Dn*vn*nRHS+i*nRHS+iRHS] = u[i];
				}
			}
		}
		// -------------- Step 3: Compute PhPolar ---------------
		for(int i=0;i<dGrain*6*6;i++){ PhPolar[i]=0.0; };
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (int irow=0; irow<6; irow++){
				for (int icol=0; icol<6; icol++){
					for (int i=0; i<Dn*vn; i++){
						PhPolar[iGrain*6*6+irow*6+icol] = PhPolar[iGrain*6*6+irow*6+icol] + B[irow*Dn*vn+i] * PhDisp[iGrain*Dn*vn*nRHS+i*6+icol];
					}
				}
			}
		}
		// -------------- Step 4: Switch PhPolar to obtain PhPolarOut ---------------
		for(int i=0;i<dGrain*6*6;i++){ PhPolarOut[i]=PhPolar[i]; };
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (int irow=0; irow<6; irow++){
				for (int icol=0; icol<6; icol++){
					PhPolarOut[iGrain*6*6+irow*6+icol] = PhPolar[iGrain*6*6+Switch[irow]*6+Switch[icol]];
				}
			}
		}
		// -------------- Step 5: Compute corresponding row of phSegsPVec ---------------
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (int irow=0; irow<6; irow++){
				for (int icol=0; icol<6; icol++){
					phSegsPVec[ie*6*6*dGrain+iGrain*6*6+irow*6+icol] = -PhPolarOut[iGrain*6*6+irow*6+icol];
				}
			}
		}
		// -------------- Step 6: Compute weight for the current element ---------------
		weights[ie] = 1.0/6.0 * Ejacs[10*ie+9];
		// -------------- Step 7: Compute element volume ---------------
		ElVol = (vert[9]-vert[0])* ( (vert[4]-vert[1])*(vert[8]-vert[2]) - (vert[5]-vert[2])*(vert[7]-vert[1]) ) +
			(vert[10]-vert[1])* ( (vert[5]-vert[2])*(vert[6]-vert[0]) - (vert[3]-vert[0])*(vert[8]-vert[2]) ) +
			(vert[11]-vert[2])* ( (vert[3]-vert[0])*(vert[7]-vert[1]) - (vert[4]-vert[1])*(vert[6]-vert[0]) );
		ElVol = ElVol / 6.0;
		PartVolume = PartVolume + ElVol;
		// ----------- Print element ID and B_matrix ------------
		/*
		of << E->elem_glid[ie];
		of << '\n';
		of << "PhDisp is ";
		of << '\n';
		of << PhDisp[0] << ", " << PhDisp[1] << ", " << PhDisp[2] << ", " << PhDisp[3] << ", " << PhDisp[4] << ", " << PhDisp[5] << '\n';
		of << PhDisp[6] << ", " << PhDisp[7] << ", " << PhDisp[8] << ", " << PhDisp[9] << ", " << PhDisp[10] << ", " << PhDisp[11] << '\n';
		of << PhDisp[12] << ", " << PhDisp[13] << ", " << PhDisp[14] << ", " << PhDisp[15] << ", " << PhDisp[16] << ", " << PhDisp[17] << '\n';
		of << PhDisp[18] << ", " << PhDisp[19] << ", " << PhDisp[20] << ", " << PhDisp[21] << ", " << PhDisp[22] << ", " << PhDisp[23] << '\n';
		of << PhDisp[24] << ", " << PhDisp[25] << ", " << PhDisp[26] << ", " << PhDisp[27] << ", " << PhDisp[28] << ", " << PhDisp[29] << '\n';
		of << PhDisp[30] << ", " << PhDisp[31] << ", " << PhDisp[32] << ", " << PhDisp[33] << ", " << PhDisp[34] << ", " << PhDisp[35] << '\n';
		of << "PhPolar is ";
		of << '\n';
		of << PhPolar[0] << ", " << PhPolar[1] << ", " << PhPolar[2] << ", " << PhPolar[3] << ", " << PhPolar[4] << ", " << PhPolar[5] << '\n';
		of << PhPolar[6] << ", " << PhPolar[7] << ", " << PhPolar[8] << ", " << PhPolar[9] << ", " << PhPolar[10] << ", " << PhPolar[11] << '\n';
		of << PhPolar[12] << ", " << PhPolar[13] << ", " << PhPolar[14] << ", " << PhPolar[15] << ", " << PhPolar[16] << ", " << PhPolar[17] << '\n';
		of << PhPolar[18] << ", " << PhPolar[19] << ", " << PhPolar[20] << ", " << PhPolar[21] << ", " << PhPolar[22] << ", " << PhPolar[23] << '\n';
		of << PhPolar[24] << ", " << PhPolar[25] << ", " << PhPolar[26] << ", " << PhPolar[27] << ", " << PhPolar[28] << ", " << PhPolar[29] << '\n';
		of << PhPolar[30] << ", " << PhPolar[31] << ", " << PhPolar[32] << ", " << PhPolar[33] << ", " << PhPolar[34] << ", " << PhPolar[35] << '\n';
		of << "PhPolarOut is ";
		of << '\n';
		of << PhPolarOut[0] << ", " << PhPolarOut[1] << ", " << PhPolarOut[2] << ", " << PhPolarOut[3] << ", " << PhPolarOut[4] << ", " << PhPolarOut[5] << '\n';
		of << PhPolarOut[6] << ", " << PhPolarOut[7] << ", " << PhPolarOut[8] << ", " << PhPolarOut[9] << ", " << PhPolarOut[10] << ", " << PhPolarOut[11] << '\n';
		of << PhPolarOut[12] << ", " << PhPolarOut[13] << ", " << PhPolarOut[14] << ", " << PhPolarOut[15] << ", " << PhPolarOut[16] << ", " << PhPolarOut[17] << '\n';
		of << PhPolarOut[18] << ", " << PhPolarOut[19] << ", " << PhPolarOut[20] << ", " << PhPolarOut[21] << ", " << PhPolarOut[22] << ", " << PhPolarOut[23] << '\n';
		of << PhPolarOut[24] << ", " << PhPolarOut[25] << ", " << PhPolarOut[26] << ", " << PhPolarOut[27] << ", " << PhPolarOut[28] << ", " << PhPolarOut[29] << '\n';
		of << PhPolarOut[30] << ", " << PhPolarOut[31] << ", " << PhPolarOut[32] << ", " << PhPolarOut[33] << ", " << PhPolarOut[34] << ", " << PhPolarOut[35] << '\n';
	//	of << PhDisp[36] << ", " << PhDisp[37] << ", " << PhDisp[38] << ", " << PhDisp[39] << ", " << PhDisp[40] << ", " << PhDisp[41] << '\n';
	//	of << PhDisp[42] << ", " << PhDisp[43] << ", " << PhDisp[44] << ", " << PhDisp[45] << ", " << PhDisp[46] << ", " << PhDisp[47] << '\n';
	//	of << PhDisp[48] << ", " << PhDisp[49] << ", " << PhDisp[50] << ", " << PhDisp[51] << ", " << PhDisp[52] << ", " << PhDisp[53] << '\n';
	//	of << PhDisp[54] << ", " << PhDisp[55] << ", " << PhDisp[56] << ", " << PhDisp[57] << ", " << PhDisp[58] << ", " << PhDisp[59] << '\n';
	//	of << PhDisp[60] << ", " << PhDisp[61] << ", " << PhDisp[62] << ", " << PhDisp[63] << ", " << PhDisp[64] << ", " << PhDisp[65] << '\n';
	//	of << PhDisp[66] << ", " << PhDisp[67] << ", " << PhDisp[68] << ", " << PhDisp[69] << ", " << PhDisp[70] << ", " << PhDisp[71] << '\n';
		of << "Element volume is " << ElVol;
		of << '\n';
		of << "weights[ie] is " << weights[ie] << '\n';
		*/
	}
	// -------------- Step 8: Compute Pijkl ---------------
	// weights * phSegsPVec
	FLOAT_PHYS WphSegsPVec[6*6*dGrain];
	for(int i=0;i<dGrain*6*6;i++){ WphSegsPVec[i]=0.0; };
	
	for (int iGrain=0; iGrain<dGrain; iGrain++){
		for (int irow=0; irow<6; irow++){
			for (int icol=0; icol<6; icol++){
				for(INT_MESH ie=e0;ie<ee;ie++){
					WphSegsPVec[iGrain*6*6+irow*6+icol] = WphSegsPVec[iGrain*6*6+irow*6+icol] + 
						weights[ie] * phSegsPVec[ie*6*6*dGrain+iGrain*6*6+irow*6+icol];
				}
			}
		}
	}
	// Compute Pijkl
	for (int iGrain=0; iGrain<part_n; iGrain++){
		for (int irow=0; irow<6; irow++){
			for (int icol=0; icol<6; icol++){
				Pijkl[iGrain*6*6+irow*6+icol] = 1.0/PartVolume * WphSegsPVec[iGrain*6*6+irow*6+icol];
			}
		}
	}
	/* Computation of Aijkl is performed in post-processing or reformating step */
//	of << "PartVolume is " << PartVolume << '\n';
	of.precision(17);
	of << PartVolume << '\n';
	/*
	of << "phSegsPVec, for ie=0, is " << '\n';
	for(INT_MESH ie=e0;ie<ee;ie++){
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (int irow=0; irow<nRHS; irow++){
				of << phSegsPVec[ie*dGrain*nRHS*nRHS+iGrain*nRHS*nRHS+irow*nRHS+0] << ", " << phSegsPVec[ie*dGrain*nRHS*nRHS+iGrain*nRHS*nRHS+irow*nRHS+1] << ", " << phSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+2] << ", "
				   << phSegsPVec[ie*dGrain*nRHS*nRHS+iGrain*nRHS*nRHS+irow*nRHS+3] << ", " << phSegsPVec[ie*dGrain*nRHS*nRHS+iGrain*nRHS*nRHS+irow*nRHS+4] << ", " << phSegsPVec[ie*dGrain*nRHS*nRHS+iGrain*nRHS*nRHS+irow*nRHS+5];
				of << '\n';
			}
		}
	}
	of << "weights is " << '\n';
	for(INT_MESH ie=e0;ie<ee;ie++){
		of << weights[ie] << ", ";
	}
	of << '\n';
	of << "WphSegsPVec is " << '\n';
	for (int iGrain=0; iGrain<dGrain; iGrain++){
		for (int irow=0; irow<nRHS; irow++){
			of << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+0] << ", " << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+1] << ", " << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+2] << ", "
			   << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+3] << ", " << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+4] << ", " << WphSegsPVec[iGrain*nRHS*nRHS+irow*nRHS+5];
			of << '\n';
		}
	}
	of << "Pijkl is " << '\n';
	*/
	for (int iGrain=0; iGrain<dGrain; iGrain++){
		for (int irow=0; irow<6; irow++){
			of << Pijkl[iGrain*6*6+irow*6+0] << ", " << Pijkl[iGrain*6*6+irow*6+1] << ", " << Pijkl[iGrain*6*6+irow*6+2] << ", "
			   << Pijkl[iGrain*6*6+irow*6+3] << ", " << Pijkl[iGrain*6*6+irow*6+4] << ", " << Pijkl[iGrain*6*6+irow*6+5];
			of << '\n';
		}
	}
	return 1;
}
int ElastIso3D::ElemJacobi( Elem*, FLOAT_SOLV*, const FLOAT_SOLV* ){
  return 1; }
int ElastIso3D::ElemJacNode(Elem* E, FLOAT_SOLV* part_d ){
#define DIAG_FROM_STIFF
  const uint Dm = 3;
  const uint Nj = 10;
  const uint Nc = E->elem_conn_n;
  const uint Ne = uint(Dm*Nc);// One row of stiffness matrix
  const uint elem_n = E->elem_n;
  const uint intp_n = E->gaus_n;
#ifdef DIAG_FROM_STIFF
  const uint Nk = Ne*Ne;// Elements of stiffness matrix
  FLOAT_PHYS elem_stiff[Nk];
#endif
  FLOAT_PHYS B[Ne*6];// 6 rows, Ne cols
  FLOAT_PHYS G[Ne];
  for(uint j=0; j<(Ne*6); j++){ B[j]=0.0; }
  const FLOAT_PHYS D[]={
    mtrl_matc[0],mtrl_matc[1],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[0],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[1],mtrl_matc[0],0.0,0.0,0.0,
    0.0,0.0,0.0,mtrl_matc[2],0.0,0.0,//FIXME Check these shear values.
    0.0,0.0,0.0,0.0,mtrl_matc[2],0.0,
    0.0,0.0,0.0,0.0,0.0,mtrl_matc[2] };
  for(uint ie=0;ie<elem_n;ie++){
#ifdef DIAG_FROM_STIFF
    for(uint i=0;i<Nk;i++){ elem_stiff[i]=0.0; }
#endif
    for(uint ip=0;ip<intp_n;ip++){
#if 1
      for(uint k=0; k<Nc; k++){
        for(uint i=0; i<Dm ; i++){ G[Nc* i+k ]=0.0;
          for(uint j=0; j<Dm ; j++){
            G[Nc* i+k ] += E->elip_jacs[Nj*ie+ Dm* j+i ]
                         * E->intp_shpg[ip*Ne+ Dm* k+j ];
      } } }
#else
      uint ig=ip*Ne;
      for(uint i=0;i<Ne;i++){ G[i]=0.0; }
      for(uint k=0;k<Nc;k++){
      for(uint i=0;i< 3;i++){
      for(uint j=0;j< 3;j++){
        G[Nc* i+k] += E->elip_jacs[Nj*ie+3* j+i] * E->intp_shpg[ig+3* k+j];
      } } }
#endif
      for(uint j=0; j<Nc; j++){
      // xx yy zz
        B[Ne*0 + 0+j*Dm] = G[Nc*0+j];
        B[Ne*1 + 1+j*Dm] = G[Nc*1+j];
        B[Ne*2 + 2+j*Dm] = G[Nc*2+j];
      // xy yx
        B[Ne*3 + 0+j*Dm] = G[Nc*1+j];
        B[Ne*3 + 1+j*Dm] = G[Nc*0+j];
      // yz zy
        B[Ne*4 + 1+j*Dm] = G[Nc*2+j];
        B[Ne*4 + 2+j*Dm] = G[Nc*1+j];
      // xz zx
        B[Ne*5 + 0+j*Dm] = G[Nc*2+j];
        B[Ne*5 + 2+j*Dm] = G[Nc*0+j];
      }
      FLOAT_PHYS dw = E->elip_jacs[Nj*ie+9] * E->gaus_weig[ip];
#ifdef DIAG_FROM_STIFF
      for(uint i=0; i<Ne; i++){
      for(uint l=0; l<Ne; l++){
      for(uint k=0; k<6 ; k++){
      for(uint j=0; j<6 ; j++){
        elem_stiff[Ne* i+l ] += B[Ne* j+i ] * D[6* k+j ] * B[Ne* k+l ] * dw;
      } } } }
#else
      for(uint n=0; n<Nc; n++){
      for(uint m=0; m<3 ; m++){ uint i=3*n+m;
      for(uint o=0; o<3 ; o++){ uint l=3*n+o;
      for(uint k=0; k<6 ; k++){
      for(uint j=0; j<6 ; j++){
        part_d[E->elem_conn[Nc*ie+n]*9+ 3* m+o ]
          += B[Ne* j+i ] * D[6* k+j ] * B[Ne* k+l ] * dw
#if 0
          *((m==o)?1.0:-1.0)
#endif
          ;
      } } } } }
#endif
    }//end intp loop
#ifdef DIAG_FROM_STIFF
    for(uint i=0; i<Nc; i++){
    for(uint k=0; k<3 ; k++){
    for(uint j=0; j<3 ; j++){
      part_d[E->elem_conn[Nc*ie+i]*9+ 3* k+j ]
        += elem_stiff[3*Ne*i +3*i   +Ne* k+j ];
    } } }
#endif
  }//end elem loop
  return 0;
}
//
int ElastIso3D::ElemNonlinear( Elem* E, const INT_MESH e0, const INT_MESH e1,
  FLOAT_SOLV* part_f, const FLOAT_SOLV* part_u, const FLOAT_SOLV*, bool ){
  return this->ElemLinear( E, e0,e1, part_f, part_u );
  }
int ElastIso3D::ElemStiff(Elem* E ){//FIXME should be ScatterStiff( E )
  const int Dm = 3;//E->mesh_d
  const int Dn = 3;//this->node_d;
  const int Nj = 10;
  const int Nc = int(E->elem_conn_n);
  const int Ne = Dm*Nc;
  const int Nr = Dn*Nc;// One row/col of stiffness matrix
#ifdef __INTEL_COMPILER
  const int Nk =int(Nr*(Nr + 1))/2;// Packed symmetric stiffness
#else
  const int Nk = Nr * Nr;// Elements of stiffness matrix
#endif
  const INT_MESH elem_n = int(E->elem_n);
  const int intp_n = int(E->gaus_n);
  //
  FLOAT_PHYS B[Ne*6];// 6 rows, Ne cols
  FLOAT_PHYS G[Ne];
  for(int j=0; j<(Ne*6); j++){ B[j]=0.0; }
  const FLOAT_PHYS D[]={
    mtrl_matc[0],mtrl_matc[1],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[0],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[1],mtrl_matc[0],0.0,0.0,0.0,
    0.0,0.0,0.0,mtrl_matc[2],0.0,0.0,//FIXME Check these shear values.
    0.0,0.0,0.0,0.0,mtrl_matc[2],0.0,
    0.0,0.0,0.0,0.0,0.0,mtrl_matc[2] };
  for(INT_MESH ie=0;ie<elem_n;ie++){
    for(int ip=0;ip<intp_n;ip++){
      int ig=ip*Ne;
      for(int i=0;i<Ne;i++){ G[i]=0.0; }
      for(int k=0;k<Nc;k++){
      for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        G[Nc* i+k] += E->elip_jacs[Nj*ie+3* j+i ] * E->intp_shpg[ig+3* k+j ];
      } } }
      for(int j=0; j<Nc; j++){
      // xx yy zz
        B[Ne*0 + 0+j*Dm] = G[Nc*0+j];
        B[Ne*1 + 1+j*Dm] = G[Nc*1+j];
        B[Ne*2 + 2+j*Dm] = G[Nc*2+j];
      // xy yx
        B[Ne*3 + 0+j*Dm] = G[Nc*1+j];
        B[Ne*3 + 1+j*Dm] = G[Nc*0+j];
      // yz zy
        B[Ne*4 + 1+j*Dm] = G[Nc*2+j];
        B[Ne*4 + 2+j*Dm] = G[Nc*1+j];
      // xz zx
        B[Ne*5 + 0+j*Dm] = G[Nc*2+j];
        B[Ne*5 + 2+j*Dm] = G[Nc*0+j];
      }
#if VERB_MAX>10
      printf( "[B]:");
      for(int j=0;j<B.size();j++){
        if(j%Ne==0){ printf("\n"); }
        printf("%+9.2e ",B[j]);
      } printf("\n");
#endif
      FLOAT_PHYS w = E->elip_jacs[Nj*ie+9] * E->gaus_weig[ip];
#ifdef __INTEL_COMPILER
      int ik=-1;
#endif
      for(int i=0; i<Nr; i++){
      for(int l=0; l<Nr; l++){
#ifdef __INTEL_COMPILER
        ik=ik+int(i<=l);
        if(i<=l){
#endif
      for(int k=0; k<6 ; k++){
      for(int j=0; j<6 ; j++){
#ifdef __INTEL_COMPILER
        // Use packed symmmetric matrix storage.
          elem_stiff[Nk*ie + ik ]+=B[Ne* j+i ] * D[6* k+j ] * B[Ne* k+l ] * w;
        }
#else
        elem_stiff[Nk*ie +Nr* i+l ]+=B[Ne* j+i ] * D[6* k+j ] * B[Ne* k+l ] * w;
#endif
      } } } }
    }// end intp loop
  }// End elem loop
  return 0;
}//============================================================== End ElemStiff
int ElastIso3D::ElemJacobi(Elem* E, FLOAT_SOLV* part_d ){
	// Xiaoyu: define # of RHS
	const uint nRHS = 6;
	const uint ndof   = 3;//this->node_d for single RHS
	// Xiaoyu: define number of UCPs with six RHSs
	const int dGrain = this->node_d / (nRHS*ndof);
	const uint  Nj = 10;
	const uint  Nc = E->elem_conn_n;
	const uint  Ne = uint(ndof*Nc);
	const uint elem_n = E->elem_n;
	const uint intp_n = E->gaus_n;
	//
	RESTRICT Phys::vals elem_diag(Ne);
	FLOAT_PHYS B[Ne*6];// 6 rows, Ne cols
	FLOAT_PHYS G[Ne];
	for(uint j=0; j<(Ne*6); j++){ B[j]=0.0; }
	const FLOAT_PHYS D[]={
		mtrl_matc[0],mtrl_matc[1],mtrl_matc[1],0.0,0.0,0.0,
		mtrl_matc[1],mtrl_matc[0],mtrl_matc[1],0.0,0.0,0.0,
		mtrl_matc[1],mtrl_matc[1],mtrl_matc[0],0.0,0.0,0.0,
		0.0,0.0,0.0,mtrl_matc[2],0.0,0.0,//FIXME Check these shear values.
		0.0,0.0,0.0,0.0,mtrl_matc[2],0.0,
		0.0,0.0,0.0,0.0,0.0,mtrl_matc[2] };
	for(uint ie=0;ie<elem_n;ie++){
		for(uint ip=0;ip<intp_n;ip++){
			for(uint i=0;i<Ne;i++){ G[i]=0.0; }
			for(uint k=0;k<Nc;k++){
				for(uint i=0;i< 3;i++){
					for(uint j=0;j< 3;j++){
						G[Nc* i+k] += E->elip_jacs[Nj*ie+3* j+i] * E->intp_shpg[ip*Ne+3* k+j];
			} } }
			for(uint j=0; j<Nc; j++){
				// xx yy zz
				B[Ne*0 + 0+j*ndof] = G[Nc*0+j];
				B[Ne*1 + 1+j*ndof] = G[Nc*1+j];
				B[Ne*2 + 2+j*ndof] = G[Nc*2+j];
				// xy yx
				B[Ne*3 + 0+j*ndof] = G[Nc*1+j];
				B[Ne*3 + 1+j*ndof] = G[Nc*0+j];
				// yz zy
				B[Ne*4 + 1+j*ndof] = G[Nc*2+j];
				B[Ne*4 + 2+j*ndof] = G[Nc*1+j];
				// xz zx
				B[Ne*5 + 0+j*ndof] = G[Nc*2+j];
				B[Ne*5 + 2+j*ndof] = G[Nc*0+j];
			}
			const FLOAT_PHYS dw = E->elip_jacs[Nj*ie+9] * E->gaus_weig[ip];
			for(uint i=0; i<Ne; i++){
				for(uint k=0; k<6 ; k++){
					for(uint j=0; j<6 ; j++){
						elem_diag[i]+= B[Ne*j + i] * D[6*j + k] * B[Ne*k + i] * dw;
			} } }
		}//end intp loop
		// Xiaoyu: assign to part_d
		for (int iGrain=0; iGrain<dGrain; iGrain++){
			for (uint iRHS=0; iRHS<nRHS; iRHS++){
				for (uint i=0; i<Nc; i++){
					for(uint j=0; j<ndof; j++){
						part_d[E->elem_conn[Nc*ie+i]*ndof*nRHS*dGrain+iGrain*nRHS*ndof+iRHS*ndof+j] += elem_diag[3*i+j];
				} }
			}
		}
	//	for (uint i=0; i<Nc; i++){
	//		for(uint j=0; j<3; j++){
	//			part_d[E->elem_conn[Nc*ie+i]*3+j] += elem_diag[3*i+j];
	//	} }
		elem_diag=0.0;
	}
	return 0;
}
int ElastIso3D::ElemRowSumAbs(Elem* E, FLOAT_SOLV* part_d ){
  const uint ndof   = 3;//this->node_d
  const uint elem_n = E->elem_n;
  const uint  Nc = E->elem_conn_n;
  const uint  Nj = 10;
  const uint  Ne = uint(ndof*Nc);
  const uint intp_n = E->gaus_n;
  //
  RESTRICT Phys::vals elem_sum(Ne);
  RESTRICT Phys::vals K(Ne*Ne);
  FLOAT_PHYS B[Ne*6];//6 rows, Ne cols
  FLOAT_PHYS G[Ne];
  for(uint j=0; j<(Ne*6); j++){ B[j]=0.0; };
  const FLOAT_PHYS D[]={
    mtrl_matc[0],mtrl_matc[1],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[0],mtrl_matc[1],0.0,0.0,0.0,
    mtrl_matc[1],mtrl_matc[1],mtrl_matc[0],0.0,0.0,0.0,
    0.0,0.0,0.0,mtrl_matc[2]*2.0,0.0,0.0,
    0.0,0.0,0.0,0.0,mtrl_matc[2]*2.0,0.0,
    0.0,0.0,0.0,0.0,0.0,mtrl_matc[2]*2.0};
  for(uint ie=0;ie<elem_n;ie++){
    for(uint ip=0;ip<intp_n;ip++){
      uint ig=ip*Ne;
      for(uint i=0;i<Ne;i++){ G[i]=0.0; };
      for(uint k=0;k<Nc;k++){
      for(uint i=0;i<3;i++){
      for(uint j=0;j<3;j++){
        G[3* i+k] += E->elip_jacs[Nj*ie+3* j+i] * E->intp_shpg[ig+3* k+j];
      } } }
      #if VERB_MAX>10
      printf( "Jacobian Inverse & Determinant:");
      for(uint j=0;j<d2;j++){
        if(j%3==0){printf("\n");}
        printf("%+9.2e",jac[j]);
      }; printf(" det:%+9.2e\n",det);
      #endif
      // xx yy zz
      for(uint j=0; j<Nc; j++){
        B[Ne*0 + 0+j*ndof] = G[Nc*0+j];
        B[Ne*1 + 1+j*ndof] = G[Nc*1+j];
        B[Ne*2 + 2+j*ndof] = G[Nc*2+j];
      // xy yx
        B[Ne*3 + 0+j*ndof] = G[Nc*1+j];
        B[Ne*3 + 1+j*ndof] = G[Nc*0+j];
      // yz zy
        B[Ne*4 + 1+j*ndof] = G[Nc*2+j];
        B[Ne*4 + 2+j*ndof] = G[Nc*1+j];
      // xz zx
        B[Ne*5 + 0+j*ndof] = G[Nc*2+j];
        B[Ne*5 + 2+j*ndof] = G[Nc*0+j];
      };
      #if VERB_MAX>10
      printf( "[B]:");
      for(uint j=0;j<B.size();j++){
        if(j%Ne==0){printf("\n");}
        printf("%+9.2e ",B[j]);
      }; printf("\n");
      #endif
      FLOAT_PHYS w = E->elip_jacs[Nj*ie+9] * E->gaus_weig[ip];
      for(uint i=0; i<Ne; i++){
      for(uint l=0; l<Ne; l++){
      for(uint j=0; j<6 ; j++){
      for(uint k=0; k<6 ; k++){
        K[Ne*i+l]+= B[Ne*j + i] * D[6*j + k] * B[Ne*k + l]*w;
      };};};};
    };//end intp loop
    for (uint i=0; i<Ne; i++){
      for(uint j=0; j<Ne; j++){
        elem_sum[i] += std::abs(K[Ne*i+j]);
      };};
    for (uint i=0; i<Nc; i++){
      for(uint j=0; j<3; j++){
        part_d[E->elem_conn[Nc*ie+i]*3+j] += elem_sum[3*i+j];
      };};
    K=0.0; elem_sum=0.0;
  };
  return 0;
};
int ElastIso3D::ElemStrain( Elem* E,FLOAT_SOLV* part_f ){
  //FIXME Clean up local variables.
  const uint ndof= 3;//this->node_d
  const uint  Nj =10;//,d2=9;//mesh_d*mesh_d;
  const INT_MESH elem_n = E->elem_n;
  const uint intp_n = uint(E->gaus_n);
  const uint     Nc = E->elem_conn_n;// Number of Nodes/Element
  const uint     Ne = ndof*Nc;
  //FLOAT_PHYS det;
  INT_MESH   conn[Nc];
  FLOAT_MESH jac[Nj];
  FLOAT_PHYS dw, G[Ne], f[Ne];
  FLOAT_PHYS H[9], S[9];
  //
  for(uint i=0; i< 9 ; i++){ H[i]=0.0; };
  H[0]=1.0; H[4]=1.0; H[8]=1.0;// unit pressure
  //
  FLOAT_PHYS intp_shpg[intp_n*Ne];
  std::copy( &E->intp_shpg[0],
             &E->intp_shpg[intp_n*Ne], intp_shpg );
  FLOAT_PHYS wgt[intp_n];
  std::copy( &E->gaus_weig[0],
             &E->gaus_weig[intp_n], wgt );
  FLOAT_PHYS C[this->mtrl_matc.size()];
  std::copy( &this->mtrl_matc[0],
             &this->mtrl_matc[this->mtrl_matc.size()], C );
  const auto Econn = &E->elem_conn[0];
  const auto Ejacs = &E->elip_jacs[0];
  //
  for(INT_MESH ie=0;ie<elem_n;ie++){
    std::memcpy( &conn, &Econn[Nc*ie], sizeof(  INT_MESH)*Nc);
    std::memcpy( &jac , &Ejacs[Nj*ie], sizeof(FLOAT_MESH)*Nj);
    //
    for(uint i=0;i<(Ne);i++){ f[i]=0.0; };
    for(uint ip=0; ip<intp_n; ip++){
      //G = MatMul3x3xN( jac,shg );
      //H = MatMul3xNx3T( G,u );// [H] Small deformation tensor
      for(uint k=0; k<Nc; k++){
        for(uint i=0; i<3 ; i++){ G[3* k+i ]=0.0;
          for(uint j=0; j<3 ; j++){
            G[(3* k+i) ] += jac[3* j+i ] * intp_shpg[ip*Ne+ 3* k+j ];
          };
        };
      };//------------------------------------------------- N*3*6*2 = 36*N FLOP
#if VERB_MAX>10
      printf( "Small Strains (Elem: %i):", ie );
      for(uint j=0;j<9;j++){
        if(j%mesh_d==0){printf("\n");}
        printf("%+9.2e ",H[j]);
      }; printf("\n");
#endif
      //det=jac[9 +Nj*l]; FLOAT_PHYS w = det * wgt[ip];
      dw = jac[9] * wgt[ip];
      //
      S[0]=(C[0]* H[0] + C[1]* H[4] + C[1]* H[8])*dw;//Sxx
      S[4]=(C[1]* H[0] + C[0]* H[4] + C[1]* H[8])*dw;//Syy
      S[8]=(C[1]* H[0] + C[1]* H[4] + C[0]* H[8])*dw;//Szz
      //
      S[1]=( H[1] + H[3] )*C[2]*dw;// S[3]= S[1];//Sxy Syx
      S[5]=( H[5] + H[7] )*C[2]*dw;// S[7]= S[5];//Syz Szy
      S[2]=( H[2] + H[6] )*C[2]*dw;// S[6]= S[2];//Sxz Szx
      S[3]=S[1]; S[7]=S[5]; S[6]=S[2];
      //------------------------------------------------------- 18+9 = 27 FLOP
      for(uint i=0; i<Nc; i++){
        for(uint k=0; k<3; k++){
          for(uint j=0; j<3; j++){
            f[(3* i+k) ] += G[(3* i+j) ] * S[(3* k+j) ];
      };};};//---------------------------------------------- N*3*6 = 18*N FLOP
#if VERB_MAX>10
      printf( "f:");
      for(uint j=0;j<Ne;j++){
        if(j%ndof==0){printf("\n");}
        printf("%+9.2e ",f[j]);
      }; printf("\n");
#endif
    };//end intp loop
    for (uint i=0; i<Nc; i++){
      for(uint j=0; j<3; j++){
        //part_f[3*conn[i]+j] += f[(3*i+j)];
        part_f[4*conn[i]+j] += std::abs( f[(3*i+j)] );
    }; };//--------------------------------------------------- N*3 =  3*N FLOP
  };//end elem loop
  return 0;
  };
#if 0
int ElastIso3D::ReadPartFMR( const char* fname, bool is_bin ){
  //FIXME This is not used. It's done in Mesh::ReadPartFMR...
  std::string s; if(is_bin){ s="binary";}else{s="ASCII";}
  if(is_bin){
    std::cout << "ERROR Could not open "<< fname << " for reading." <<'\n'
      << "ERROR Femera (fmr) "<< s <<" format not yet supported." <<'\n';
    return 1;
  }
  std::string fmrstring;
  std::ifstream fmrfile(fname);
  while( fmrfile >> fmrstring ){
    if(fmrstring=="$ElasticProperties"){//FIXME Deprecated
      int s=0; fmrfile >> s;
      mtrl_prop.resize(s);
      for(int i=0; i<s; i++){ fmrfile >> mtrl_prop[i]; }
      //this->MtrlProp2MatC();
      s=0; fmrfile >> s;
      if(s>0){
        mtrl_dirs.resize(s);
      for(int i=0; i<s; i++){ fmrfile >> mtrl_dirs[i]; mtrl_dirs[i]*=(PI/180.0) ;}
      }
    }
    if(fmrstring=="$Orientation"){// Material orientation (radians)
      int s=0; fmrfile >> s;
      if(s>0){
        mtrl_dirs.resize(s);
        for(int i=0; i<s; i++){ fmrfile >> mtrl_dirs[i]; mtrl_dirs[i]*=(PI/180.0) ;}
      }
    }
    //FIXME This parsing requires properties in a specific order
    auto tprop = mtrl_prop; auto tsz=tprop.size();
    if(fmrstring=="$Elastic"){// Elastic Constants
      int s=0; fmrfile >> s;
      mtrl_prop.resize(tsz+s);
      mtrl_prop[std::slice(tsz,tsz+s,1)] = tprop;
      for(int i=0; i<s; i++){ fmrfile >> mtrl_prop[i+tprop.size()]; }
    }
    if(fmrstring=="$ThermalExpansion"){// Thermal expansion
      int s=0; fmrfile >> s;
      mtrl_prop.resize(s + tprop.size());
      mtrl_prop[std::slice(tsz,tsz+s,1)] = tprop;
      for(int i=0; i<s; i++){ fmrfile >> mtrl_prop[i+tprop.size()]; }
    }
    if(fmrstring=="$ThermalConductivity"){// Thermal conductivity
      int s=0; fmrfile >> s;
      mtrl_prop.resize(s + tprop.size());
      mtrl_prop[std::slice(tsz,tsz+s,1)] = tprop;
      for(int i=0; i<s; i++){ fmrfile >> mtrl_prop[i+tprop.size()]; }
    }
  }
  return 0;
}
int ElastIso3D::SavePartFMR( const char* fname, bool is_bin ){
  std::string s; if(is_bin){ s="binary";}else{s="ASCII";};
  if(is_bin){
    std::cout << "ERROR Could not append "<< fname << "." <<'\n'
      << "ERROR Femera (fmr) "<< s <<" format not yet supported." <<'\n';
    return 1;
  };
  std::ofstream fmrfile;
  fmrfile.open(fname, std::ios_base::app);
  //
  fmrfile << "$ElasticProperties" <<'\n';
  fmrfile << mtrl_prop.size();
  for(uint i=0;i<mtrl_prop.size();i++){ fmrfile <<" "<< mtrl_prop[i]; };
  fmrfile << '\n';
  if(mtrl_dirs.size()>0){
    fmrfile << mtrl_dirs.size();
    for(uint i=0;i<mtrl_dirs.size();i++){ fmrfile <<" "<< mtrl_dirs[i]; };
  }; fmrfile <<'\n';
  return 0;
};
#endif
