# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>
# include <time.h>
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
extern void   StrainHS_Nikkhoo(float Stress[6], float Strain[6],   float X,        float Y,        float Z,        float P1[3],   float P2[3], float P3[3], float SS, float Ds, float Ts, float mu, float lambda);

void    RotateTensor_inKmatrix(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir);
void       GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3]);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
void     Build_K_Matrix(float *fK_SS, float *fK_SD, float *fK_SO, float *fK_DS, float *fK_DD, float *fK_DO, unsigned int *iTDg_SelfIndLoc, const float *fVDg_Epos,  const float *fVDg_Npos, const float *fVDg_Zpos,  const unsigned int *iTDg_V1, const unsigned int *iTDg_V2, const unsigned int *iTDg_V3, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, float fMDg_UnitSlip, float fMDg_ShearMod, float fMDg_Lambda, unsigned int *iTDl_StabType, const float *fTDl_Curr_DcVal, const float *fTDl_StatFric, const float *fTDl_DynFric, const float *fTDl_RefNormStrss, const int *iSTARTPOS, const int *iOFFSET, int iRANK, unsigned int iPatchNum, const float *fTDl_Area)
{	
    unsigned int    i,                  j,				iVectPos;
    unsigned int    globi;
    float  	        fX,                 fY,         	fZ; 
    float   		fP1s[3],           	fP2s[3],        fP3s[3];
    float   		fP1r[3],            fP2r[3],        fP3r[3];
    float   		fvNrm[3],           fvStk[3],       fvDip[3];
    float   		fStress[6],         fStressOut[6];
    float           fTemp1,				fTemp2;

	//FILE *fpIn;
	//if ((fpIn = fopen("KmatrixOutput.asc","w")) == NULL){   printf("Error -cant open *.flt file. in Main function -trying to open _Summary.flt here... file...\n");      exit(10);     }

    for (i = iOFFSET[iRANK]; i--;   )
    {	globi   = i + iSTARTPOS[iRANK];
		
    	fP1s[0] = fVDg_Epos[iTDg_V1[globi]];        fP1s[1] = fVDg_Npos[iTDg_V1[globi]];        fP1s[2] = fVDg_Zpos[iTDg_V1[globi]];      
        fP2s[0] = fVDg_Epos[iTDg_V2[globi]];        fP2s[1] = fVDg_Npos[iTDg_V2[globi]];        fP2s[2] = fVDg_Zpos[iTDg_V2[globi]]; 
        fP3s[0] = fVDg_Epos[iTDg_V3[globi]];        fP3s[1] = fVDg_Npos[iTDg_V3[globi]];        fP3s[2] = fVDg_Zpos[iTDg_V3[globi]]; 
          
        for (j = iPatchNum; j--;  ) /* going through the receivers */
        {   
   			fP1r[0]   = fVDg_Epos[iTDg_V1[j]];        fP1r[1] = fVDg_Npos[iTDg_V1[j]];        fP1r[2] = fVDg_Zpos[iTDg_V1[j]];      
        	fP2r[0]   = fVDg_Epos[iTDg_V2[j]];        fP2r[1] = fVDg_Npos[iTDg_V2[j]];        fP2r[2] = fVDg_Zpos[iTDg_V2[j]]; 
        	fP3r[0]   = fVDg_Epos[iTDg_V3[j]];        fP3r[1] = fVDg_Npos[iTDg_V3[j]];        fP3r[2] = fVDg_Zpos[iTDg_V3[j]]; 
        	fX        = fTDg_CentEpos[j];			  fY      = fTDg_CentNpos[j];			  fZ      = fTDg_CentZpos[j];
   			
   			iVectPos  = i*iPatchNum + j;

   			if (globi == j)     		{	iTDg_SelfIndLoc[globi] = iVectPos; 						} 
     
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);	

            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2s, fP3s, fMDg_UnitSlip, 0.0, 0.0, fMDg_ShearMod, fMDg_Lambda);  	   /* slip in stk */            
			//fprintf(fpIn,"i %d  j %d   %f   %f   %f   %f   %f   %f   ", i, j, fStress[0],fStress[3],fStress[5],fStress[1],fStress[2],fStress[4]);

            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);

            fK_SS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_SD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_SO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip; 
                    
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2s, fP3s, 0.0, fMDg_UnitSlip, 0.0, fMDg_ShearMod, fMDg_Lambda);  	  /* slip in dip */            
            //fprintf(fpIn,"%f   %f   %f   %f   %f   %f   \n", fStress[0],fStress[3],fStress[5],fStress[1],fStress[2],fStress[4]);

            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);

            fK_DS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_DD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_DO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;   
    	}	
    	/*-----------------------------------------------*/	          
    	/*determine friction difference, use together with reference normal stress to get stress drop; get corresponding slip amount, compare with Dc => assign type	*/
    	fTemp1 = (fTDl_DynFric[i] - fTDl_StatFric[i]) *fTDl_RefNormStrss[i];
        fTemp2 = fTemp1/(fK_SS[iTDg_SelfIndLoc[globi]] > fK_DD[iTDg_SelfIndLoc[globi]] ? fK_SS[iTDg_SelfIndLoc[globi]] : fK_DD[iTDg_SelfIndLoc[globi]]);//(0.5*(fK_SS[iTDg_SelfIndLoc[globi]] + fK_DD[iTDg_SelfIndLoc[globi]])));
            
        if (fTemp2 > fTDl_Curr_DcVal[i]) /* have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc*/
        {	iTDl_StabType[i] = 1; /* patch is unstable*/
        }
        else /*if (fTemp2 <= fTDl_Curr_DcVal[i]) -if stress drop i.e., corresponding slip is smaller than Dc - */
        {	if (fTemp1 < 0.0)  /*  is still weakening but just with slip that is lower than Dc => cond. stable*/
        	{	iTDl_StabType[i] = 2; /* patch is cond. stable*/
        	}
        	else /* no weakening but strengthening  => stable */
        	{	iTDl_StabType[i] = 3; /* patch is stable*/
        }	}
    	/*-----------------------------------------------*/	
    }
    //fclose(fpIn);
    return;
} 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
void RotateTensor_inKmatrix(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir)
{   
    float fTxx1,          fTxy1,         fTxz1,          fTyy1,         fTyz1,            fTzz1;
    float fTxx2,          fTxy2,         fTxz2,          fTyy2,         fTyz2,            fTzz2;
    float fA[9]; /* the "Rotation Matrix" */

    fTxx1       = fsig[0];             fTxy1       = fsig[1];            fTxz1       = fsig[2];
    fTyy1       = fsig[3];             fTyz1       = fsig[4];            fTzz1       = fsig[5];

    /* TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate system to x2y2z2 coordinate system. "A" is the transformation matrix, */
    /* whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The  coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose  of A (i.e., A') does the transformation from x2y2z2 into x1y1z1. */

    if (iRotDir == 1) /* from local to global */
    {   fA[0] = fvNrm[0];         fA[3] = fvStk[0];        fA[6] = fvDip[0];                
        fA[1] = fvNrm[1];         fA[4] = fvStk[1];        fA[7] = fvDip[1];
        fA[2] = fvNrm[2];         fA[5] = fvStk[2];        fA[8] = fvDip[2];
    }
    else             /* from global to local */
    {   fA[0] = fvNrm[0];         fA[3] = fvNrm[1];        fA[6] = fvNrm[2];
        fA[1] = fvStk[0];         fA[4] = fvStk[1];        fA[7] = fvStk[2];
        fA[2] = fvDip[0];         fA[5] = fvDip[1];        fA[8] = fvDip[2];
    }

    fTxx2 = fA[0]*fA[0]*fTxx1 +           2.0*fA[0]*fA[3] *fTxy1 +            2.0*fA[0]*fA[6] *fTxz1 +            2.0*fA[3]*fA[6] *fTyz1 + fA[3]*fA[3]*fTyy1 + fA[6]*fA[6]*fTzz1;
    fTyy2 = fA[1]*fA[1]*fTxx1 +           2.0*fA[1]*fA[4] *fTxy1 +            2.0*fA[1]*fA[7] *fTxz1 +            2.0*fA[4]*fA[7] *fTyz1 + fA[4]*fA[4]*fTyy1 + fA[7]*fA[7]*fTzz1;
    fTzz2 = fA[2]*fA[2]*fTxx1 +           2.0*fA[2]*fA[5] *fTxy1 +            2.0*fA[2]*fA[8] *fTxz1 +            2.0*fA[5]*fA[8] *fTyz1 + fA[5]*fA[5]*fTyy1 + fA[8]*fA[8]*fTzz1;  
    fTxy2 = fA[0]*fA[1]*fTxx1 + (fA[0]*fA[4]+ fA[1]*fA[3])*fTxy1 + (fA[0]*fA[7] + fA[1]*fA[6])*fTxz1 + (fA[7]*fA[3] + fA[6]*fA[4])*fTyz1 + fA[4]*fA[3]*fTyy1 + fA[6]*fA[7]*fTzz1;
    fTxz2 = fA[0]*fA[2]*fTxx1 + (fA[0]*fA[5]+ fA[2]*fA[3])*fTxy1 + (fA[0]*fA[8] + fA[2]*fA[6])*fTxz1 + (fA[8]*fA[3] + fA[6]*fA[5])*fTyz1 + fA[5]*fA[3]*fTyy1 + fA[6]*fA[8]*fTzz1;
    fTyz2 = fA[1]*fA[2]*fTxx1 + (fA[2]*fA[4]+ fA[1]*fA[5])*fTxy1 + (fA[2]*fA[7] + fA[1]*fA[8])*fTxz1 + (fA[7]*fA[5] + fA[8]*fA[4])*fTyz1 + fA[4]*fA[5]*fTyy1 + fA[7]*fA[8]*fTzz1;

    fsig_new[0] = fTxx2;        fsig_new[1] = fTxy2;      fsig_new[2] = fTxz2;
    fsig_new[3] = fTyy2;        fsig_new[4] = fTyz2;      fsig_new[5] = fTzz2;
 
    return;
}
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
void GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3])
{   /* this comes basically from Medhii's code... */
	float fTempfloat;
	float ftempVect1[3],       	ftempVect2[3];  
	float feY[3],			    feZ[3];

    feY[0] = 0.0;                feY[1] = 1.0;            feY[2] = 0.0;
    feZ[0] = 0.0;                feZ[1] = 0.0;            feZ[2] = 1.0;    

	ftempVect1[0] = fP2[0] -fP1[0];            ftempVect1[1] = fP2[1] -fP1[1];            ftempVect1[2] = fP2[2] -fP1[2];
    ftempVect2[0] = fP3[0] -fP1[0];            ftempVect2[1] = fP3[1] -fP1[1];            ftempVect2[2] = fP3[2] -fP1[2];
      
    fvNrm[0]      = ftempVect1[1]*ftempVect2[2] - ftempVect1[2]*ftempVect2[1];
    fvNrm[1]      = ftempVect1[2]*ftempVect2[0] - ftempVect1[0]*ftempVect2[2];
    fvNrm[2]      = ftempVect1[0]*ftempVect2[1] - ftempVect1[1]*ftempVect2[0];
    fTempfloat    = sqrtf(fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]);
    fvNrm[0]      = fvNrm[0]/fTempfloat;      fvNrm[1]     = fvNrm[1]/fTempfloat;      fvNrm[2]     = fvNrm[2]/fTempfloat;

	fvStk[0]      = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    fvStk[1]      = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    fvStk[2]      = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    /* For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) */
    fTempfloat    = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
	if (fTempfloat == 0.0)
    {   fvStk[0] = feY[0]*fvNrm[2];         fvStk[1] = feY[1]*fvNrm[2];             fvStk[2] = feY[2]*fvNrm[2];        
        /* For horizontal elements in case of half-space calculation!!! => Correct the strike vector of image dislocation only */
        if (fP1[2] > 0.0)
        {   fvStk[0] = -1.0*fvStk[0];       fvStk[1] = -1.0*fvStk[1];               fvStk[2] = -1.0*fvStk[2];
    }   }
	fTempfloat  = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    fvStk[0]    = fvStk[0]/fTempfloat;      fvStk[1] = fvStk[1]/fTempfloat;         fvStk[2] = fvStk[2]/fTempfloat;

    fvDip[0]    = fvNrm[1]*fvStk[2] - fvNrm[2]*fvStk[1];
    fvDip[1]    = fvNrm[2]*fvStk[0] - fvNrm[0]*fvStk[2];
    fvDip[2]    = fvNrm[0]*fvStk[1] - fvNrm[1]*fvStk[0];
    fTempfloat  = sqrtf(fvDip[0]*fvDip[0] +fvDip[1]*fvDip[1] +fvDip[2]*fvDip[2]);
    fvDip[0]    = fvDip[0]/fTempfloat;      fvDip[1] = fvDip[1]/fTempfloat;         fvDip[2] = fvDip[2]/fTempfloat;

    return;
}
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/