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
extern void   StrainHS_Nikkhoo(float Stress[6], float Strain[6], float X, float Y, float Z, float P1[3], float P2[3], float P3[3], float SS, float Ds, float Ts, const float mu, const float lambda);

void    RotateTensor_inKmatrix(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir);
void       GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3]);
void             GetPointOrder(float P1in[3], float P2in[3], float P3in[3], float fP2s[3], float fP3s[3]);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
void     Build_K_Matrix(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const  int iFltPtchNum, const  int iBndPtchNum, const float fMDg_UnitSlip, const float *fMDg_ShearMod, const float *fMDg_Lambda, const  int *iTDg_V1, const  int *iTDg_V2, const  int *iTDg_V3, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fTDl_Curr_DcVal, const float *fTDl_StatFric, const float *fTDl_DynFric, const float *fTDl_RefNormStrss, const float *fTDl_Area, float *fK_FF_SS, float *fK_FF_SD, float *fK_FF_SO, float *fK_FF_DS, float *fK_FF_DD, float *fK_FF_DO, float *fK_FF_OS, float *fK_FF_OD, float *fK_FF_OO, float *fK_FB_SS, float *fK_FB_SD, float *fK_FB_SO, float *fK_FB_DS, float *fK_FB_DD, float *fK_FB_DO, float *fK_BF_SS, float *fK_BF_SD, float *fK_BF_SO, float *fK_BF_DS, float *fK_BF_DD, float *fK_BF_DO, float *fK_BF_OS, float *fK_BF_OD, float *fK_BF_OO, float *fK_BB_SS, float *fK_BB_SD, float *fK_BB_SO, float *fK_BB_DS, float *fK_BB_DD, float *fK_BB_DO, float *fK_BB_OS, float *fK_BB_OD, float *fK_BB_OO, float *fKl_BB_SS, float *fKl_BB_DD, float *fKl_BB_OO,  int *iTDl_SelfLoc_F,  int *iTDl_SelfLoc_B,  int *iTDl_StabType)
{    
    int    i,                  j,              iVectPos;
    int    globi;
    float  fX,                 fY,             fZ; 
    float  fP1s[3],            fP2s[3],        fP3s[3];
    float  fP1r[3],            fP2r[3],        fP3r[3];
    float  fP2sOut[3],         fP3sOut[3];
    float  fvNrm[3],           fvStk[3],       fvDip[3];
    
    float  fStress[6],         fStressOut[6];
    float  fTemp1,             fTemp2;
 //   FILE   *fp1;
    /*----------------------------------------------------------------------------------*/
    /*----------------------------------------------------------------------------------*/     
 //   if ((fp1 = fopen("IntMatrix_m1.txt","w")) == NULL)         {       exit(10);     }
             
    for (j = 0; j < iFltPtchNum;  j++) /* going through the receivers */
    {  
        fP1r[0]   = fVDg_Epos[iTDg_V1[j]];        fP1r[1] = fVDg_Npos[iTDg_V1[j]];        fP1r[2] = fVDg_Zpos[iTDg_V1[j]];      
        fP2r[0]   = fVDg_Epos[iTDg_V2[j]];        fP2r[1] = fVDg_Npos[iTDg_V2[j]];        fP2r[2] = fVDg_Zpos[iTDg_V2[j]]; 
        fP3r[0]   = fVDg_Epos[iTDg_V3[j]];        fP3r[1] = fVDg_Npos[iTDg_V3[j]];        fP3r[2] = fVDg_Zpos[iTDg_V3[j]]; 
       fX        = fTDg_CentEpos[j];             fY      = fTDg_CentNpos[j];             fZ      = fTDg_CentZpos[j];
        /*-----------------------------------------------*/           
        GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);   
     //   fprintf(stdout,"%f     %f    %f        %f     %f    %f  \n",fvNrm[0], fvNrm[1], fvNrm[2], fvStk[0], fvStk[1], fvStk[2]);
        
        /*-----------------------------------------------*/           
  //    fX        = fTDg_CentEpos[j]-fvNrm[0]*fMDg_UnitSlip;             fY      = fTDg_CentNpos[j]-fvNrm[1]*fMDg_UnitSlip;             fZ      = fTDg_CentZpos[j]-fvNrm[2]*fMDg_UnitSlip;
            
        for (i = 0; i < iOFFSET_F[iRANK]; i++)
        {    
            globi     = i + iSTARTPOS_F[iRANK];
            iVectPos  = i*iFltPtchNum + j;
            
            if (globi == j)             {    iTDl_SelfLoc_F[i] = iVectPos;                         } 
        
            fP1s[0] = fVDg_Epos[iTDg_V1[globi]];        fP1s[1] = fVDg_Npos[iTDg_V1[globi]];        fP1s[2] = fVDg_Zpos[iTDg_V1[globi]];      
            fP2s[0] = fVDg_Epos[iTDg_V2[globi]];        fP2s[1] = fVDg_Npos[iTDg_V2[globi]];        fP2s[2] = fVDg_Zpos[iTDg_V2[globi]]; 
            fP3s[0] = fVDg_Epos[iTDg_V3[globi]];        fP3s[1] = fVDg_Npos[iTDg_V3[globi]];        fP3s[2] = fVDg_Zpos[iTDg_V3[globi]]; 
            /*-----------------------------------------------*/             
            GetPointOrder(fP1s, fP2s, fP3s, fP2sOut, fP3sOut); //to make sure that fault normal will point upward
            /*-----------------------------------------------*/       
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, fMDg_UnitSlip, 0.0, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in stk */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_FF_SS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_FF_SD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_FF_SO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip; 
     //       fprintf(fp1,"%f    %f    %f       %f    %f    %f    ", fX, fY, fZ, fK_FF_SS[iVectPos], fK_FF_SD[iVectPos] , fK_FF_SO[iVectPos] );
            
            /*-----------------------------------------------*/               
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, fMDg_UnitSlip, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);        /* slip in dip */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_FF_DS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_FF_DD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_FF_DO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
    //         fprintf(fp1,"%f    %f    %f    ", fK_FF_DS[iVectPos], fK_FF_DD[iVectPos] , fK_FF_DO[iVectPos] );
                
            /*-----------------------------------------------*/       
             StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, 0.0, fMDg_UnitSlip, fMDg_ShearMod[0], fMDg_Lambda[0]);        /* slip in dip */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_FF_OS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_FF_OD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_FF_OO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;   
  //           fprintf(fp1,"%f    %f    %f    ", fK_FF_OS[iVectPos], fK_FF_OD[iVectPos] , fK_FF_OO[iVectPos] );
             
            /*-----------------------------------------------*/ 
 //            fprintf(fp1,"\n");                       
        }   
       
        /*-----------------------------------------------*/        
        /*-----------------------------------------------*/        
        for (i = 0; i < iOFFSET_B[iRANK]; i++)
        {    fprintf(stdout,"Iam here even though I don't have boundary faults right now");
            globi     = i + iSTARTPOS_B[iRANK] + iFltPtchNum; //the plus iFltPtchNum is here b/c the boundary values are following the fault values in the respective list
            iVectPos  = i*iFltPtchNum + j;
            
            fP1s[0] = fVDg_Epos[iTDg_V1[globi]];        fP1s[1] = fVDg_Npos[iTDg_V1[globi]];        fP1s[2] = fVDg_Zpos[iTDg_V1[globi]];      
            fP2s[0] = fVDg_Epos[iTDg_V2[globi]];        fP2s[1] = fVDg_Npos[iTDg_V2[globi]];        fP2s[2] = fVDg_Zpos[iTDg_V2[globi]]; 
            fP3s[0] = fVDg_Epos[iTDg_V3[globi]];        fP3s[1] = fVDg_Npos[iTDg_V3[globi]];        fP3s[2] = fVDg_Zpos[iTDg_V3[globi]]; 
            /*-----------------------------------------------*/             
            GetPointOrder(fP1s, fP2s, fP3s, fP2sOut, fP3sOut);
            /*-----------------------------------------------*/       
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, fMDg_UnitSlip, 0.0, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in stk */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BF_SS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_BF_SD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BF_SO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/               
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, fMDg_UnitSlip, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);        /* slip in dip */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BF_DS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_BF_DD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BF_DO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/                 
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, 0.0, fMDg_UnitSlip, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in normal */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BF_OS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_BF_OD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BF_OO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/             
    }   }
//    fclose(fp1);
            for (i = 0; i < iOFFSET_F[iRANK]; i++)
        {
            
            if ((fK_FF_SS[iTDl_SelfLoc_F[i]] > 0.0) || ( fK_FF_DD[iTDl_SelfLoc_F[i]] > 0.0 )  || (fK_FF_OO[iTDl_SelfLoc_F[i]] > 0.0))
            {
                fprintf(stdout,"%f    %f     %f   \n",fK_FF_SS[iTDl_SelfLoc_F[i]],  fK_FF_DD[iTDl_SelfLoc_F[i]], fK_FF_OO[iTDl_SelfLoc_F[i]]);
            }
        }
    
    /*----------------------------------------------------------------------------------*/    
    /*----------------------------------------------------------------------------------*/    
    for (j = 0; j < iBndPtchNum;  j++) /* going through the receivers */
    {   
        fP1r[0]   = fVDg_Epos[iTDg_V1[j+iFltPtchNum]];        fP1r[1] = fVDg_Npos[iTDg_V1[j+iFltPtchNum]];        fP1r[2] = fVDg_Zpos[iTDg_V1[j+iFltPtchNum]];      
        fP2r[0]   = fVDg_Epos[iTDg_V2[j+iFltPtchNum]];        fP2r[1] = fVDg_Npos[iTDg_V2[j+iFltPtchNum]];        fP2r[2] = fVDg_Zpos[iTDg_V2[j+iFltPtchNum]]; 
        fP3r[0]   = fVDg_Epos[iTDg_V3[j+iFltPtchNum]];        fP3r[1] = fVDg_Npos[iTDg_V3[j+iFltPtchNum]];        fP3r[2] = fVDg_Zpos[iTDg_V3[j+iFltPtchNum]]; 
        fX        = fTDg_CentEpos[j+iFltPtchNum];             fY      = fTDg_CentNpos[j+iFltPtchNum];             fZ      = fTDg_CentZpos[j+iFltPtchNum];
        /*-----------------------------------------------*/           
        GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);    
        /*-----------------------------------------------*/
        /*-----------------------------------------------*/                   
        for (i = 0; i < iOFFSET_B[iRANK]; i++)
        {    
            globi     = i + iSTARTPOS_B[iRANK] +iFltPtchNum;
            iVectPos  = i*iBndPtchNum + j;
            
            if (globi == (j+iFltPtchNum))          {    iTDl_SelfLoc_B[i] = iVectPos;                         } 
    
            fP1s[0] = fVDg_Epos[iTDg_V1[globi]];        fP1s[1] = fVDg_Npos[iTDg_V1[globi]];        fP1s[2] = fVDg_Zpos[iTDg_V1[globi]];      
            fP2s[0] = fVDg_Epos[iTDg_V2[globi]];        fP2s[1] = fVDg_Npos[iTDg_V2[globi]];        fP2s[2] = fVDg_Zpos[iTDg_V2[globi]]; 
            fP3s[0] = fVDg_Epos[iTDg_V3[globi]];        fP3s[1] = fVDg_Npos[iTDg_V3[globi]];        fP3s[2] = fVDg_Zpos[iTDg_V3[globi]]; 
            /*-----------------------------------------------*/             
            GetPointOrder(fP1s, fP2s, fP3s, fP2sOut, fP3sOut);
            /*-----------------------------------------------*/          
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, fMDg_UnitSlip, 0.0, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in stk */            
           /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BB_SS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_BB_SD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BB_SO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/               
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, fMDg_UnitSlip, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);        /* slip in dip */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BB_DS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_BB_DD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BB_DO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/                 
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, 0.0, fMDg_UnitSlip, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in normal */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_BB_OS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_BB_OD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_BB_OO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;
            /*-----------------------------------------------*/   
            if (globi == (j+iFltPtchNum))       {       fKl_BB_SS[i] = fK_BB_SS[iVectPos];       fKl_BB_DD[i] = fK_BB_DD[iVectPos];      fKl_BB_OO[i] = fK_BB_OO[iVectPos];     }
        }
        /*-----------------------------------------------*/        
        /*-----------------------------------------------*/            
        for (i = 0; i < iOFFSET_F[iRANK]; i++)
        {    
            globi     = i + iSTARTPOS_F[iRANK];
            iVectPos  = i*iBndPtchNum + j;
             
            fP1s[0] = fVDg_Epos[iTDg_V1[globi]];        fP1s[1] = fVDg_Npos[iTDg_V1[globi]];        fP1s[2] = fVDg_Zpos[iTDg_V1[globi]];      
            fP2s[0] = fVDg_Epos[iTDg_V2[globi]];        fP2s[1] = fVDg_Npos[iTDg_V2[globi]];        fP2s[2] = fVDg_Zpos[iTDg_V2[globi]]; 
            fP3s[0] = fVDg_Epos[iTDg_V3[globi]];        fP3s[1] = fVDg_Npos[iTDg_V3[globi]];        fP3s[2] = fVDg_Zpos[iTDg_V3[globi]]; 
            /*-----------------------------------------------*/             
            GetPointOrder(fP1s, fP2s, fP3s, fP2sOut, fP3sOut);
            /*-----------------------------------------------*/   
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, fMDg_UnitSlip, 0.0, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);         /* slip in stk */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_FB_SS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip; 
            fK_FB_SD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_FB_SO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip; 
            /*-----------------------------------------------*/               
            StrainHS_Nikkhoo(fStress, fStressOut, fX, fY, fZ, fP1s, fP2sOut, fP3sOut, 0.0, fMDg_UnitSlip, 0.0, fMDg_ShearMod[0], fMDg_Lambda[0]);        /* slip in dip */            
            /*-----------------------------------------------*/ 
            RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
            /*-----------------------------------------------*/ 
            fK_FB_DS[iVectPos]  = fStressOut[1]/fMDg_UnitSlip;
            fK_FB_DD[iVectPos]  = fStressOut[2]/fMDg_UnitSlip;
            fK_FB_DO[iVectPos]  = fStressOut[0]/fMDg_UnitSlip;   
            /*-----------------------------------------------*/                              
    }   }
    /*----------------------------------------------------------------------------------*/    
    /*----------------------------------------------------------------------------------*/    
    for (i = 0; i < iOFFSET_F[iRANK]; i++)
    {   /*determine friction difference, use together with reference normal stress to get stress drop; get corresponding slip amount, compare with Dc => assign type    */        
        fTemp1 = (fTDl_DynFric[i] - fTDl_StatFric[i]) *fTDl_RefNormStrss[i];
        fTemp2 = fTemp1/(fK_FF_SS[iTDl_SelfLoc_F[i]] < fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]]);
        
        if (fTemp2 > fTDl_Curr_DcVal[i]) /* have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc*/
        {         iTDl_StabType[i] = 1; /* patch is unstable*/
        }
        else /*if (fTemp2 <= fTDl_Curr_DcVal[i]) -if stress drop i.e., corresponding slip is smaller than Dc - */
        {    if (fTemp1 < 0.0)  /*  is still weakening but just with slip that is lower than Dc => cond. stable*/
             {    iTDl_StabType[i] = 2; /* patch is cond. stable*/
             }
             else /* no weakening but strengthening  => stable */
             {    iTDl_StabType[i] = 3; /* patch is stable*/ 
   }    }   
 
   }  
    /*-----------------------------------------------*/    
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
    /*-----------------------------------------------*/  
    fTxx1       = fsig[0];             fTxy1       = fsig[1];            fTxz1       = fsig[2];
    fTyy1       = fsig[3];             fTyz1       = fsig[4];            fTzz1       = fsig[5];
    /*-----------------------------------------------*/  
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
    /*-----------------------------------------------*/  
    fTxx2 = fA[0]*fA[0]*fTxx1 +           2.0*fA[0]*fA[3] *fTxy1 +            2.0*fA[0]*fA[6] *fTxz1 +            2.0*fA[3]*fA[6] *fTyz1 + fA[3]*fA[3]*fTyy1 + fA[6]*fA[6]*fTzz1;
    fTyy2 = fA[1]*fA[1]*fTxx1 +           2.0*fA[1]*fA[4] *fTxy1 +            2.0*fA[1]*fA[7] *fTxz1 +            2.0*fA[4]*fA[7] *fTyz1 + fA[4]*fA[4]*fTyy1 + fA[7]*fA[7]*fTzz1;
    fTzz2 = fA[2]*fA[2]*fTxx1 +           2.0*fA[2]*fA[5] *fTxy1 +            2.0*fA[2]*fA[8] *fTxz1 +            2.0*fA[5]*fA[8] *fTyz1 + fA[5]*fA[5]*fTyy1 + fA[8]*fA[8]*fTzz1;  
    fTxy2 = fA[0]*fA[1]*fTxx1 + (fA[0]*fA[4]+ fA[1]*fA[3])*fTxy1 + (fA[0]*fA[7] + fA[1]*fA[6])*fTxz1 + (fA[7]*fA[3] + fA[6]*fA[4])*fTyz1 + fA[4]*fA[3]*fTyy1 + fA[6]*fA[7]*fTzz1;
    fTxz2 = fA[0]*fA[2]*fTxx1 + (fA[0]*fA[5]+ fA[2]*fA[3])*fTxy1 + (fA[0]*fA[8] + fA[2]*fA[6])*fTxz1 + (fA[8]*fA[3] + fA[6]*fA[5])*fTyz1 + fA[5]*fA[3]*fTyy1 + fA[6]*fA[8]*fTzz1;
    fTyz2 = fA[1]*fA[2]*fTxx1 + (fA[2]*fA[4]+ fA[1]*fA[5])*fTxy1 + (fA[2]*fA[7] + fA[1]*fA[8])*fTxz1 + (fA[7]*fA[5] + fA[8]*fA[4])*fTyz1 + fA[4]*fA[5]*fTyy1 + fA[7]*fA[8]*fTzz1;
    /*-----------------------------------------------*/  
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
    float ftempVect1[3],         ftempVect2[3];  
    float feY[3],                feZ[3];
    /*-----------------------------------------------*/  
    feY[0] = 0.0;                              feY[1] = 1.0;                              feY[2] = 0.0;
    feZ[0] = 0.0;                              feZ[1] = 0.0;                              feZ[2] = 1.0;        
    ftempVect1[0] = fP2[0] -fP1[0];            ftempVect1[1] = fP2[1] -fP1[1];            ftempVect1[2] = fP2[2] -fP1[2];
    ftempVect2[0] = fP3[0] -fP1[0];            ftempVect2[1] = fP3[1] -fP1[1];            ftempVect2[2] = fP3[2] -fP1[2];
    /*-----------------------------------------------*/  
    fvNrm[0]      = ftempVect1[1]*ftempVect2[2] - ftempVect1[2]*ftempVect2[1];
    fvNrm[1]      = ftempVect1[2]*ftempVect2[0] - ftempVect1[0]*ftempVect2[2];
    fvNrm[2]      = ftempVect1[0]*ftempVect2[1] - ftempVect1[1]*ftempVect2[0];
    fTempfloat    = sqrtf(fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]);
    fvNrm[0]      = fvNrm[0]/fTempfloat;      fvNrm[1]     = fvNrm[1]/fTempfloat;      fvNrm[2]     = fvNrm[2]/fTempfloat;
    
    if (fvNrm[2] < 0.0)     {   fvNrm[0] = -fvNrm[0];            fvNrm[1] = -fvNrm[1];            fvNrm[2] = -fvNrm[2];     }
    /*-----------------------------------------------*/  
    fvStk[0]      = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    fvStk[1]      = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    fvStk[2]      = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    /* For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) */
    fTempfloat    = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    if (fTempfloat < 0.0)
    {   fvStk[0] = 0.0;                 fvStk[1] = feY[1]*fvNrm[2];                 fvStk[2] = 0.0;        
        /* For horizontal elements in case of half-space calculation!!! => Correct the strike vector of image dislocation only */
        if (fP1[2] > 0.0)
        {   fvStk[0] = -fvStk[0];       fvStk[1] = -fvStk[1];                       fvStk[2] = -fvStk[2];
    }   }
    fTempfloat  = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    fvStk[0]    = fvStk[0]/fTempfloat;      fvStk[1] = fvStk[1]/fTempfloat;         fvStk[2] = fvStk[2]/fTempfloat;
    /*-----------------------------------------------*/  
    fvDip[0]    = fvNrm[1]*fvStk[2] - fvNrm[2]*fvStk[1];
    fvDip[1]    = fvNrm[2]*fvStk[0] - fvNrm[0]*fvStk[2];
    fvDip[2]    = fvNrm[0]*fvStk[1] - fvNrm[1]*fvStk[0];
    fTempfloat  = sqrtf(fvDip[0]*fvDip[0] +fvDip[1]*fvDip[1] +fvDip[2]*fvDip[2]);
    fvDip[0]    = fvDip[0]/fTempfloat;      fvDip[1] = fvDip[1]/fTempfloat;         fvDip[2] = fvDip[2]/fTempfloat;
    /*-----------------------------------------------*/  
    return;
}
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
void GetPointOrder(float fP1in[3], float fP2in[3], float fP3in[3], float fP2s[3], float fP3s[3])
{
    float tempVect1[3],                     tempVect2[3],                           Vnorm[3];
    
    tempVect1[0] = fP2in[0] -fP1in[0];      tempVect1[1] = fP2in[1] -fP1in[1];      tempVect1[2] = fP2in[2] -fP1in[2];
    tempVect2[0] = fP3in[0] -fP1in[0];      tempVect2[1] = fP3in[1] -fP1in[1];      tempVect2[2] = fP3in[2] -fP1in[2];
      
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
           
    if (Vnorm[2] < 0.0)
    {   
        fP2s[0]      = fP3in[0];            fP2s[1]      = fP3in[1];                fP2s[2]      = fP3in[2];
        fP3s[0]      = fP2in[0];            fP3s[1]      = fP2in[1];                fP3s[2]      = fP2in[2];
    }
    else
    {   fP2s[0]      = fP2in[0];            fP2s[1]      = fP2in[1];                fP2s[2]      = fP2in[2];
        fP3s[0]      = fP3in[0];            fP3s[1]      = fP3in[1];                fP3s[2]      = fP3in[2];
    }
    return;
}
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/