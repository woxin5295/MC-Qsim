#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
struct MDstruct
{   int         iRANK,              		iSIZE;
    int       	iF_BASEelem,        		iF_ADDelem,				iB_BASEelem,       			iB_ADDelem;  
    int       	*ivF_OFFSET,        		*ivF_START,				*ivB_OFFSET,        		*ivB_START;
	int       	*ivF_OFFSET2,        		*ivF_START2;

    int       	iRunNum,            		iSeedStart,				iUseProp,           		iUsePropTest,           iUsePSeis;
    float       fViscRelTime,     		    fAftrSlipTime,			fISeisStep,         		fTStpInSec,				fRecLgth;		
	float		fHealFact;

    int       	iFPNum,             		iBPNum,					iFVNum,             		iBVNum; 
    int       	iFSegNum,          		 	iBSegNum,				iFUsdGrd;
    
    int       	iMaxIterat,      			iEQcntr,      			iMaxMRFlgth;
    int       	iMaxSTFlgth,				iWritePos,          	iReadPos;
    int      	iChgBtwEQs,					iGlobTTmax;
    float     	fLegLgth,          			fUnitSlip,				fDeltT,             		fCutSlip;
	float 		fFltLegs,					fBndLegs,               fMinMag4Prop;

    float       fg,							fMedDense,     			fAddNrmStrs,				fVp;
    float       fPoisson,           		fLambda,        		fShearMod,					fVs;
    float       fISeisTStp,         		fTimeYears,             fVpVsRatio;   

    char        cInputName[512];
};
//------------------------------------------------------------------
struct SGstruct
{	float     	*fvSlipRate,				*fvSlipRake,			*fvStrsRate;
	int         *ivFricLawUSED;
    float       *fvRefStaFr,				*fvRefDynFr,			*fvRefCritDist;
	float       *fvRefStaFr_vari,			*fvRefDynFr_vari,		*fvRefCritDist_vari;       
};
//------------------------------------------------------------------
struct TRstruct
{	int     	*ivFL_StabT,				*ivFG_Activated,		*ivFG_Ptch_t0;
    int    		*ivFL_FricLaw,				*ivFL_SelfLoc,			*ivBL_SelfLoc;			

    float    	*fvFL_RefNrmStrs,			*fvFL_Area;

    float    	*fvFL_RefStaFric,           *fvFL_RefDynFric,    	*fvFL_RefDcVal;
    float    	*fvFL_RefStaFric_vari,      *fvFL_RefDynFric_vari,	*fvFL_RefDcVal_vari,		*fvFL_AccumSlp;
    float    	*fvFL_StaFric,              *fvFL_DynFric,			*fvFL_CurFric,              *fvFL_B4_Fric,			*fvFL_TempRefFric,		*fvFL_CurDcVal;

    float    	*fvFL_PSeis_T0_F,          	*fvFL_PSeis_T0_N,     	*fvFL_PSeis_Step,			*fvBL_CurFric;
    float  		*fvBL_PSeis_T0_F,           *fvBL_PSeis_T0_N,    	*fvBL_PSeis_Step;

	int      	*ivFL_SegID_temp,			*ivFG_FltID_temp,		*ivFG_Flagged_temp;
	float    	*fvFL_SlipRate_temp,      	*fvFL_SlipRake_temp; 
    int    		*ivFG_V1_temp,              *ivFG_V2_temp,        	*ivFG_V3_temp;
    int    		*ivBG_V1_temp,              *ivBG_V2_temp,       	*ivBG_V3_temp;
    float    	*fvFG_CentE_temp,          	*fvFG_CentN_temp,    	*fvFG_CentZ_temp;
    float    	*fvBG_CentE_temp,           *fvBG_CentN_temp,    	*fvBG_CentZ_temp;

	gsl_matrix_int	    *imFGL_TTP,         *imFGL_TTS;
	gsl_matrix_int    	*imFGL_NextP,		*imFGL_NextS; 
    gsl_matrix_float    *fmFGL_SrcRcvH,    	*fmFGL_SrcRcvV,      	*fmFGL_SrcRcvN;    

	gsl_vector_float	*fvFL_StrsRateStk, 	*fvFL_StrsRateDip;
	gsl_vector_float	*fvFL_CurStrsH,		*fvFL_CurStrsV,			*fvFL_CurStrsN;
	gsl_vector_float	*fvBL_CurStrsH,		*fvBL_CurStrsV,			*fvBL_CurStrsN;

	gsl_vector_float    *fvFL_B4_StrsH,		*fvFL_B4_StrsV,		    *fvFL_B4_StrsN;    
};
//------------------------------------------------------------------
struct VTstruct
{	float  		*fvFG_PosE_temp,          	*fvFG_PosN_temp,      	*fvFG_PosZ_temp;
	float   	*fvBG_PosE_temp,           	*fvBG_PosN_temp,      	*fvBG_PosZ_temp;
};
//------------------------------------------------------------------
struct Kstruct
{	gsl_matrix_float 	*FF_SS,				*FF_SD,						*FF_SO; //first index is "receiver"; second is "source"
	gsl_matrix_float 	*FF_DS,				*FF_DD,						*FF_DO; //this means: FB_DS => Fault is receiver, Boundary the source; DipStress is induced by StrikeSlip 
	gsl_matrix_float 	*FF_OS,				*FF_OD,						*FF_OO;

	gsl_matrix_float 	*FB_SS,				*FB_SD,						*FB_SO;
	gsl_matrix_float 	*FB_DS,				*FB_DD,						*FB_DO;
	gsl_matrix_float 	*FB_OS,				*FB_OD,						*FB_OO;

	gsl_matrix_float 	*BF_SS,				*BF_SD,						*BF_SO;
	gsl_matrix_float 	*BF_DS,				*BF_DD,						*BF_DO;
	gsl_matrix_float 	*BF_OS,				*BF_OD,						*BF_OO;

	gsl_matrix_float 	*BB_SS,				*BB_SD,						*BB_SO;
	gsl_matrix_float 	*BB_DS,				*BB_DD,						*BB_DO;
	gsl_matrix_float 	*BB_OS,				*BB_OD,						*BB_OO;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
extern void   StrainHS_Nikkhoo(float Stress[6], float Strain[6], float X, float Y, float Z, float P1[3], float P2[3], float P3[3], float SS, float Ds, float Ts, const float mu, const float lambda);

void        GetPtOrderAndPlane(float fP1in[3], float fP2in[3], float fP3in[3], float fP2out[3], float fP3out[3],float fNrm[3], float fd[1], float fArea[1]);
void       GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3]);
float            GetDist_2Pnts(float fT[3], float fL[3]);
void    RotateTensor_inKmatrix(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir);
int               GetSamePtNum(int   iSrcTrigs[3], int iRcvTrigs[3]);
int               PtLiesInTrig(float fP1[3], float fP2[3], float fP3[3], float fNrm[3], float fd[1], float fArea[1], float fT[3]); 
float                  GetArea(float fP1[3], float fP2[3], float fP3[3]);
float     GetLotPunkt_PktEbene(float fNrm[3], float fd[1], float fT[3], float fL[3]);
float           GetMinFraction(float fP1[3], float fP2[3], float fP3[3], float fT[3]);
void     GetLotPunkt_PktGerade(float fP1[3], float fP2[3], float fT[3], float fL[3]);
void              GetNewPoints(float CurFrac, float fP1in[3], float fP2in[3], float fP3in[3], float fP1out[3], float fP2out[3], float fP3out[3]);
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void     Build_K_Matrix( struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct  Kstruct *K)
{    
    int     i,                  j,              k,                  iGlobPos,           iTrigNum = 10;         
    float   fTemp0,             fTemp1,         fTemp2,             fTemp3,             fTemp4,         fTemp5;
    float   fTemp6,             fTemp7,         fTemp8,             fTestDist,          fTestFrac1,     fTestFrac2;
    float   fPartlMom,          fPartlSlip,     fCurFrac,           fTstArea,           fSrcLgth,       fRcvLgth;
    int     iPntCntr,           iFlagPair,      iSrcFltID,          iRcvFltID,          iSrcTrigs[3],   iRcvTrigs[3];
    float   fP1s[3],            fP2s[3],        fP3s[3];
    float   fP1r[3],            fP2r[3],        fP3r[3];
    float   fP1n[3],            fP2n[3],        fP3n[3];
    float   fP2sOut[3],         fP3sOut[3],     fP2rOut[3],         fP3rOut[3];
    float   fvNrm[3],           fvStk[3],       fvDip[3];
    float   fSrcPt[3],          fRcvPt[3],      fTestPt[3];
    float   fStress[6],         fStressOut[6];
    float   fnSrc[3],           fnRcv[3],       fdSrc[1],           fdRcv[1],           fAreaSrc[1],    fAreaRcv[1];

    int     iOverRuled = 0,     iCntFlags_F = 0;
    float   fDistFrac   = 1.0; //if distance between observation point and source triangle (it's legs) is less than this fraction of the average length of this triangle's legs... => then...
    float   fMaxLgth    = (MD->fFltLegs > MD->fBndLegs) ? MD->fFltLegs : MD->fBndLegs;
    //----------------------------------------------------------------------------------          
    for (i = 0; i < MD->iFPNum;  i++)   //going through the sources
    {   iSrcTrigs[0] = TR->ivFG_V1_temp[i];                    iSrcTrigs[1] = TR->ivFG_V2_temp[i];                 iSrcTrigs[2] = TR->ivFG_V3_temp[i];
        iSrcFltID    = TR->ivFG_FltID_temp[i];
        fP1s[0]      = VT->fvFG_PosE_temp[iSrcTrigs[0]];       fP1s[1]  = VT->fvFG_PosN_temp[iSrcTrigs[0]];        fP1s[2]  = VT->fvFG_PosZ_temp[iSrcTrigs[0]];      
        fP2s[0]      = VT->fvFG_PosE_temp[iSrcTrigs[1]];       fP2s[1]  = VT->fvFG_PosN_temp[iSrcTrigs[1]];        fP2s[2]  = VT->fvFG_PosZ_temp[iSrcTrigs[1]]; 
        fP3s[0]      = VT->fvFG_PosE_temp[iSrcTrigs[2]];       fP3s[1]  = VT->fvFG_PosN_temp[iSrcTrigs[2]];        fP3s[2]  = VT->fvFG_PosZ_temp[iSrcTrigs[2]];                
        fSrcPt[0]    = TR->fvFG_CentE_temp[i];                 fSrcPt[1]= TR->fvFG_CentN_temp[i];                  fSrcPt[2]= TR->fvFG_CentZ_temp[i];
        fSrcLgth     = GetDist_2Pnts(fP1s, fP2s);              fSrcLgth+= GetDist_2Pnts(fP2s, fP3s);               fSrcLgth+= GetDist_2Pnts(fP3s, fP1s);          
        fSrcLgth    /= 3.0;
        //-----------------------------------------------   
        GetPtOrderAndPlane(fP1s, fP2s, fP3s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc); //to make sure that fault normal will point upward
        //----------------------------------------------- 
        fPartlMom = fAreaSrc[0] *MD->fUnitSlip /(float)(iTrigNum);
        //-----------------------------------------------
        for (j = 0; j < MD->ivF_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivF_START[MD->iRANK];     
            iRcvTrigs[0] = TR->ivFG_V1_temp[iGlobPos];              iRcvTrigs[1] = TR->ivFG_V2_temp[iGlobPos];          iRcvTrigs[2] = TR->ivFG_V3_temp[iGlobPos];               
            iRcvFltID    = TR->ivFG_FltID_temp[iGlobPos];
            fP1r[0]      = VT->fvFG_PosE_temp[iRcvTrigs[0]];        fP1r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[0]];        fP1r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[0]];      
            fP2r[0]      = VT->fvFG_PosE_temp[iRcvTrigs[1]];        fP2r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[1]];        fP2r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[1]]; 
            fP3r[0]      = VT->fvFG_PosE_temp[iRcvTrigs[2]];        fP3r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[2]];        fP3r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[2]]; 
            fRcvPt[0]    = TR->fvFG_CentE_temp[iGlobPos];           fRcvPt[1]= TR->fvFG_CentN_temp[iGlobPos];           fRcvPt[2]= TR->fvFG_CentZ_temp[iGlobPos];    
            fRcvLgth     = GetDist_2Pnts(fP1r, fP2r);               fRcvLgth+= GetDist_2Pnts(fP2r, fP3r);               fRcvLgth+= GetDist_2Pnts(fP3r, fP1r);          
            fRcvLgth    /= 3.0;
            //-----------------------------------------------   
            GetPtOrderAndPlane(fP1r, fP2r, fP3r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv); //to make sure that fault normal will point upward
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2rOut, fP3rOut);   //to rotate into the coordinate system of the receiver!
            //-----------------------------------------------       
            fTemp0   = 0.0;           fTemp1 = 0.0;           fTemp2 = 0.0; 
            fTemp3   = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;
            fTemp6   = 0.0;           fTemp7 = 0.0;           fTemp8 = 0.0;
            //-------------------------------  
            fTestDist= GetDist_2Pnts(fSrcPt, fRcvPt);            

            if ((iGlobPos == i) || (fTestDist > 2.0*fMaxLgth))     //is "self" or if distance is far enough away to not need the extra distance check 
            {   if (iGlobPos == i)      {      TR->ivFL_SelfLoc[j] = iGlobPos;          }               
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, MD->fUnitSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp1 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp2 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, MD->fUnitSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp4 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp5 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, 0.0, MD->fUnitSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp7 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp8 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
            }          
            else
            {   //------------------------------------------------------------------------------------------ 
                iFlagPair = 0;
                iPntCntr  = GetSamePtNum(iSrcTrigs, iRcvTrigs);
                if (iPntCntr == 0)              
                {   //----------------------------------------------
                    fTestDist = GetLotPunkt_PktEbene(fnRcv, fdRcv, fSrcPt, fTestPt); //TestPt is here the "lotpunkt"  i.e., the projection of SrcPt onto the triangle plane
                    iPntCntr  = PtLiesInTrig(fP1r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv, fTestPt); 
                    if ((iPntCntr > 0)     && (fTestDist <= 1.0*fDistFrac*fRcvLgth))     {    iFlagPair = 1;      }
                    //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestDist = GetLotPunkt_PktEbene(fnSrc, fdSrc, fRcvPt, fTestPt);
                        iPntCntr  = PtLiesInTrig(fP1s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc, fTestPt); //now the "TestPt is projection point of RcvPt onto triangle plane of source
                        if ((iPntCntr > 0) && (fTestDist <= 1.0*fDistFrac*fSrcLgth))     {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestFrac1 = GetMinFraction(fP1s, fP2sOut, fP3sOut, fRcvPt);
                        fTestFrac2 = GetMinFraction(fP1r, fP2rOut, fP3rOut, fSrcPt);
                        if ( (fTestFrac1 < fDistFrac) || (fTestFrac2 < fDistFrac) )     {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                }
                if (iFlagPair == 1)
                {   if (iSrcFltID == iRcvFltID)         {       iOverRuled++;           TR->ivFG_Flagged_temp[i] = 2;           }
                    else                                {       iCntFlags_F++;          TR->ivFG_Flagged_temp[i] = 1;           }
                }
                //------------------------------------------------------------------------------------------                   
                if ((iFlagPair == 0) || (iSrcFltID == iRcvFltID))
                {   //-----------------------------------------------
                    for (k = 0; k < iTrigNum; k++)
                    {   
                        fCurFrac   = (float)(k+1)/(float)iTrigNum;
                        GetNewPoints(fCurFrac, fP1s, fP2sOut, fP3sOut, fP1n, fP2n, fP3n);
                        fTstArea   = GetArea(fP1n,fP2n,fP3n);
                        fPartlSlip = fPartlMom/fTstArea;
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, fPartlSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp0 += fStressOut[0];        fTemp1 += fStressOut[1];           fTemp2 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, fPartlSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp3 += fStressOut[0];        fTemp4 += fStressOut[1];           fTemp5 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, 0.0, fPartlSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp6 += fStressOut[0];        fTemp7 += fStressOut[1];           fTemp8 += fStressOut[2];                         
                    }
                    fTemp0 = fTemp0/MD->fUnitSlip *1.0e-6;          fTemp1 = fTemp1/MD->fUnitSlip *1.0e-6;          fTemp2 = fTemp2/MD->fUnitSlip *1.0e-6;
                    fTemp3 = fTemp3/MD->fUnitSlip *1.0e-6;          fTemp4 = fTemp4/MD->fUnitSlip *1.0e-6;          fTemp5 = fTemp5/MD->fUnitSlip *1.0e-6;
                    fTemp6 = fTemp6/MD->fUnitSlip *1.0e-6;          fTemp7 = fTemp7/MD->fUnitSlip *1.0e-6;          fTemp8 = fTemp8/MD->fUnitSlip *1.0e-6;                              
            }   }   
            //------------------------------------------------------------------------------------------        
            gsl_matrix_float_set(K->FF_SO, i, j, fTemp0);   gsl_matrix_float_set(K->FF_SS, i, j, fTemp1);   gsl_matrix_float_set(K->FF_SD, i, j, fTemp2);  
            gsl_matrix_float_set(K->FF_DO, i, j, fTemp3);   gsl_matrix_float_set(K->FF_DS, i, j, fTemp4);   gsl_matrix_float_set(K->FF_DD, i, j, fTemp5);  
            gsl_matrix_float_set(K->FF_OO, i, j, fTemp6);   gsl_matrix_float_set(K->FF_OS, i, j, fTemp7);   gsl_matrix_float_set(K->FF_OD, i, j, fTemp8);       
            //------------------------------------------------------------------------------------------        
        }    
        //------------------------------------------------------------------------------------------
        for (j = 0; j < MD->ivB_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivB_START[MD->iRANK];     
            iRcvTrigs[0] = TR->ivBG_V1_temp[iGlobPos];          iRcvTrigs[1] = TR->ivBG_V2_temp[iGlobPos];          iRcvTrigs[2] = TR->ivBG_V3_temp[iGlobPos];               
            fP1r[0]  = VT->fvBG_PosE_temp[iRcvTrigs[0]];        fP1r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[0]];        fP1r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[0]];      
            fP2r[0]  = VT->fvBG_PosE_temp[iRcvTrigs[1]];        fP2r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[1]];        fP2r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[1]]; 
            fP3r[0]  = VT->fvBG_PosE_temp[iRcvTrigs[2]];        fP3r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[2]];        fP3r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[2]]; 
            fRcvPt[0]= TR->fvBG_CentE_temp[iGlobPos];           fRcvPt[1]= TR->fvBG_CentN_temp[iGlobPos];           fRcvPt[2]= TR->fvBG_CentZ_temp[iGlobPos];    
            fRcvLgth     = GetDist_2Pnts(fP1r, fP2r);               fRcvLgth+= GetDist_2Pnts(fP2r, fP3r);               fRcvLgth+= GetDist_2Pnts(fP3r, fP1r);          
            fRcvLgth    /= 3.0;
            //-----------------------------------------------   
            GetPtOrderAndPlane(fP1r, fP2r, fP3r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv); //to make sure that fault normal will point upward
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2rOut, fP3rOut);   //to rotate into the coordinate system of the receiver!
            //-----------------------------------------------             
            fTemp0   = 0.0;           fTemp1 = 0.0;           fTemp2 = 0.0; 
            fTemp3   = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;
            fTemp6   = 0.0;           fTemp7 = 0.0;           fTemp8 = 0.0;
            //-------------------------------  
            fTestDist= GetDist_2Pnts(fSrcPt, fRcvPt);            

            if (fTestDist > 2.0*fMaxLgth)    //is "self" or if distance is far enough away to not need the extra distance check 
            {               
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, MD->fUnitSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp1 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp2 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, MD->fUnitSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp4 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp5 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, 0.0, MD->fUnitSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp7 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp8 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
            }          
            else
            {   //------------------------------------------------------------------------------------------ 
                iFlagPair  = 0;
                fTestDist = GetLotPunkt_PktEbene(fnRcv, fdRcv, fSrcPt, fTestPt); //TestPt is here the "lotpunkt"  i.e., the projection of SrcPt onto the triangle plane
                iPntCntr  = PtLiesInTrig(fP1r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv, fTestPt); 
                if ((iPntCntr > 0)     && (fTestDist <= 1.0*fDistFrac*fRcvLgth))     {    iFlagPair = 1;      }
                //----------------------------------------------
                if (iFlagPair == 0)
                {   fTestDist = GetLotPunkt_PktEbene(fnSrc, fdSrc, fRcvPt, fTestPt);
                    iPntCntr  = PtLiesInTrig(fP1s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc, fTestPt); //now the "TestPt is projection point of RcvPt onto triangle plane of source
                    if ((iPntCntr > 0) && (fTestDist <= 1.0*fDistFrac*fSrcLgth))     {    iFlagPair = 1;      }
                }
                //----------------------------------------------
                if (iFlagPair == 0)
                {   fTestFrac1 = GetMinFraction(fP1s, fP2sOut, fP3sOut, fRcvPt);
                    fTestFrac2 = GetMinFraction(fP1r, fP2rOut, fP3rOut, fSrcPt);
                    if ( (fTestFrac1 < fDistFrac) || (fTestFrac2 < fDistFrac) )     {    iFlagPair = 1;      }
                }
                //------------------------------------------------------------------------------------------                   
                if (iFlagPair == 0) 
                {   //-----------------------------------------------
                    for (k = 0; k < iTrigNum; k++)
                    {   
                        fCurFrac   = (float)(k+1)/(float)iTrigNum;
                        GetNewPoints(fCurFrac, fP1s, fP2sOut, fP3sOut, fP1n, fP2n, fP3n);
                        fTstArea   = GetArea(fP1n,fP2n,fP3n);
                        fPartlSlip = fPartlMom/fTstArea;
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, fPartlSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp0 += fStressOut[0];        fTemp1 += fStressOut[1];           fTemp2 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, fPartlSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp3 += fStressOut[0];        fTemp4 += fStressOut[1];           fTemp5 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, 0.0, fPartlSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp6 += fStressOut[0];        fTemp7 += fStressOut[1];           fTemp8 += fStressOut[2];                         
                    }
                    fTemp0 = fTemp0/MD->fUnitSlip *1.0e-6;          fTemp1 = fTemp1/MD->fUnitSlip *1.0e-6;          fTemp2 = fTemp2/MD->fUnitSlip *1.0e-6;
                    fTemp3 = fTemp3/MD->fUnitSlip *1.0e-6;          fTemp4 = fTemp4/MD->fUnitSlip *1.0e-6;          fTemp5 = fTemp5/MD->fUnitSlip *1.0e-6;
                    fTemp6 = fTemp6/MD->fUnitSlip *1.0e-6;          fTemp7 = fTemp7/MD->fUnitSlip *1.0e-6;          fTemp8 = fTemp8/MD->fUnitSlip *1.0e-6;                                    
            }   } 
            //------------------------------------------------------------------------------------------        
            gsl_matrix_float_set(K->FB_SO, i, j, fTemp0);   gsl_matrix_float_set(K->FB_SS, i, j, fTemp1);   gsl_matrix_float_set(K->FB_SD, i, j, fTemp2);  
            gsl_matrix_float_set(K->FB_DO, i, j, fTemp3);   gsl_matrix_float_set(K->FB_DS, i, j, fTemp4);   gsl_matrix_float_set(K->FB_DD, i, j, fTemp5);  
            gsl_matrix_float_set(K->FB_OO, i, j, fTemp6);   gsl_matrix_float_set(K->FB_OS, i, j, fTemp7);   gsl_matrix_float_set(K->FB_OD, i, j, fTemp8);       
            //------------------------------------------------------------------------------------------        
        } 
    }
    //----------------------------------------------------------------------------------  
    //----------------------------------------------------------------------------------          
    for (i = 0; i < MD->iBPNum;  i++)   //going through the sources
    {   iSrcTrigs[0] = TR->ivBG_V1_temp[i];                 iSrcTrigs[1] = TR->ivBG_V2_temp[i];                 iSrcTrigs[2] = TR->ivBG_V3_temp[i];

        fP1s[0]      = VT->fvBG_PosE_temp[iSrcTrigs[0]];       fP1s[1]  = VT->fvBG_PosN_temp[iSrcTrigs[0]];        fP1s[2]  = VT->fvBG_PosZ_temp[iSrcTrigs[0]];      
        fP2s[0]      = VT->fvBG_PosE_temp[iSrcTrigs[1]];       fP2s[1]  = VT->fvBG_PosN_temp[iSrcTrigs[1]];        fP2s[2]  = VT->fvBG_PosZ_temp[iSrcTrigs[1]]; 
        fP3s[0]      = VT->fvBG_PosE_temp[iSrcTrigs[2]];       fP3s[1]  = VT->fvBG_PosN_temp[iSrcTrigs[2]];        fP3s[2]  = VT->fvBG_PosZ_temp[iSrcTrigs[2]];                
        fSrcPt[0]    = TR->fvBG_CentE_temp[i];                 fSrcPt[1]= TR->fvBG_CentN_temp[i];                  fSrcPt[2]= TR->fvBG_CentZ_temp[i];
        fSrcLgth     = GetDist_2Pnts(fP1s, fP2s);              fSrcLgth+= GetDist_2Pnts(fP2s, fP3s);               fSrcLgth+= GetDist_2Pnts(fP3s, fP1s);          
        fSrcLgth    /= 3.0;
        //-----------------------------------------------   
        GetPtOrderAndPlane(fP1s, fP2s, fP3s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc); //to make sure that fault normal will point upward
        //----------------------------------------------- 
        fPartlMom = fAreaSrc[0] *MD->fUnitSlip /(float)(iTrigNum);
        //-----------------------------------------------
        for (j = 0; j < MD->ivB_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivB_START[MD->iRANK];     
            iRcvTrigs[0] = TR->ivBG_V1_temp[iGlobPos];          iRcvTrigs[1] = TR->ivBG_V2_temp[iGlobPos];          iRcvTrigs[2] = TR->ivBG_V3_temp[iGlobPos];               
            
            fP1r[0]      = VT->fvBG_PosE_temp[iRcvTrigs[0]];        fP1r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[0]];        fP1r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[0]];      
            fP2r[0]      = VT->fvBG_PosE_temp[iRcvTrigs[1]];        fP2r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[1]];        fP2r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[1]]; 
            fP3r[0]      = VT->fvBG_PosE_temp[iRcvTrigs[2]];        fP3r[1]  = VT->fvBG_PosN_temp[iRcvTrigs[2]];        fP3r[2]  = VT->fvBG_PosZ_temp[iRcvTrigs[2]]; 
            fRcvPt[0]    = TR->fvBG_CentE_temp[iGlobPos];           fRcvPt[1]= TR->fvBG_CentN_temp[iGlobPos];           fRcvPt[2]= TR->fvBG_CentZ_temp[iGlobPos];    
            fRcvLgth     = GetDist_2Pnts(fP1r, fP2r);               fRcvLgth+= GetDist_2Pnts(fP2r, fP3r);               fRcvLgth+= GetDist_2Pnts(fP3r, fP1r);          
            fRcvLgth    /= 3.0;
            //-----------------------------------------------   
            GetPtOrderAndPlane(fP1r, fP2r, fP3r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv); //to make sure that fault normal will point upward
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2rOut, fP3rOut);   //to rotate into the coordinate system of the receiver!
            //-----------------------------------------------             
            fTemp0   = 0.0;           fTemp1 = 0.0;           fTemp2 = 0.0; 
            fTemp3   = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;
            fTemp6   = 0.0;           fTemp7 = 0.0;           fTemp8 = 0.0;
            //-------------------------------  
            fTestDist= GetDist_2Pnts(fSrcPt, fRcvPt);            

            if ((iGlobPos == i) || (fTestDist > 2.0*fMaxLgth))     //is "self" or if distance is far enough away to not need the extra distance check 
            {   if (iGlobPos == i)      {      TR->ivBL_SelfLoc[j] = iGlobPos;          }               
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, MD->fUnitSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp1 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp2 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, MD->fUnitSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp4 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp5 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, 0.0, MD->fUnitSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp7 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp8 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
            }          
            else
            {   //------------------------------------------------------------------------------------------ 
                iFlagPair = 0;
                iPntCntr  = GetSamePtNum(iSrcTrigs, iRcvTrigs);
                if (iPntCntr == 0)              
                {   fTestDist = GetLotPunkt_PktEbene(fnRcv, fdRcv, fSrcPt, fTestPt); //TestPt is here the "lotpunkt"  i.e., the projection of SrcPt onto the triangle plane
                    iPntCntr  = PtLiesInTrig(fP1r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv, fTestPt); 
                    if ((iPntCntr > 0)     && (fTestDist <= 1.0*fDistFrac*fRcvLgth))        {    iFlagPair = 1;      }
                    //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestDist = GetLotPunkt_PktEbene(fnSrc, fdSrc, fRcvPt, fTestPt);
                        iPntCntr  = PtLiesInTrig(fP1s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc, fTestPt); //now the "TestPt is projection point of RcvPt onto triangle plane of source
                        if ((iPntCntr > 0) && (fTestDist <= 1.0*fDistFrac*fSrcLgth))        {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestFrac1 = GetMinFraction(fP1s, fP2sOut, fP3sOut, fRcvPt);
                        fTestFrac2 = GetMinFraction(fP1r, fP2rOut, fP3rOut, fSrcPt);
                        if ( (fTestFrac1 < fDistFrac) || (fTestFrac2 < fDistFrac) )         {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                }
                //------------------------------------------------------------------------------------------                   
                if (iFlagPair == 0) 
                {   //-----------------------------------------------
                    for (k = 0; k < iTrigNum; k++)
                    {   
                        fCurFrac   = (float)(k+1)/(float)iTrigNum;
                        GetNewPoints(fCurFrac, fP1s, fP2sOut, fP3sOut, fP1n, fP2n, fP3n);
                        fTstArea   = GetArea(fP1n,fP2n,fP3n);
                        fPartlSlip = fPartlMom/fTstArea;
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, fPartlSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp0 += fStressOut[0];        fTemp1 += fStressOut[1];           fTemp2 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, fPartlSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp3 += fStressOut[0];        fTemp4 += fStressOut[1];           fTemp5 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, 0.0, fPartlSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp6 += fStressOut[0];        fTemp7 += fStressOut[1];           fTemp8 += fStressOut[2];                         
                    }
                    fTemp0 = fTemp0/MD->fUnitSlip *1.0e-6;          fTemp1 = fTemp1/MD->fUnitSlip *1.0e-6;          fTemp2 = fTemp2/MD->fUnitSlip *1.0e-6;
                    fTemp3 = fTemp3/MD->fUnitSlip *1.0e-6;          fTemp4 = fTemp4/MD->fUnitSlip *1.0e-6;          fTemp5 = fTemp5/MD->fUnitSlip *1.0e-6;
                    fTemp6 = fTemp6/MD->fUnitSlip *1.0e-6;          fTemp7 = fTemp7/MD->fUnitSlip *1.0e-6;          fTemp8 = fTemp8/MD->fUnitSlip *1.0e-6;                              
            }   }   
            //------------------------------------------------------------------------------------------        
            gsl_matrix_float_set(K->BB_SO, i, j, fTemp0);   gsl_matrix_float_set(K->BB_SS, i, j, fTemp1);   gsl_matrix_float_set(K->BB_SD, i, j, fTemp2);  
            gsl_matrix_float_set(K->BB_DO, i, j, fTemp3);   gsl_matrix_float_set(K->BB_DS, i, j, fTemp4);   gsl_matrix_float_set(K->BB_DD, i, j, fTemp5);  
            gsl_matrix_float_set(K->BB_OO, i, j, fTemp6);   gsl_matrix_float_set(K->BB_OS, i, j, fTemp7);   gsl_matrix_float_set(K->BB_OD, i, j, fTemp8);       
            //------------------------------------------------------------------------------------------        
        }    
        //------------------------------------------------------------------------------------------
        for (j = 0; j < MD->ivF_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivF_START[MD->iRANK];     
            iRcvTrigs[0] = TR->ivFG_V1_temp[iGlobPos];          iRcvTrigs[1] = TR->ivFG_V2_temp[iGlobPos];          iRcvTrigs[2] = TR->ivFG_V3_temp[iGlobPos];               
            fP1r[0]  = VT->fvFG_PosE_temp[iRcvTrigs[0]];        fP1r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[0]];        fP1r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[0]];      
            fP2r[0]  = VT->fvFG_PosE_temp[iRcvTrigs[1]];        fP2r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[1]];        fP2r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[1]]; 
            fP3r[0]  = VT->fvFG_PosE_temp[iRcvTrigs[2]];        fP3r[1]  = VT->fvFG_PosN_temp[iRcvTrigs[2]];        fP3r[2]  = VT->fvFG_PosZ_temp[iRcvTrigs[2]]; 
            fRcvPt[0]= TR->fvFG_CentE_temp[iGlobPos];           fRcvPt[1]= TR->fvFG_CentN_temp[iGlobPos];           fRcvPt[2]= TR->fvFG_CentZ_temp[iGlobPos];    
            fRcvLgth     = GetDist_2Pnts(fP1r, fP2r);               fRcvLgth+= GetDist_2Pnts(fP2r, fP3r);               fRcvLgth+= GetDist_2Pnts(fP3r, fP1r);          
            fRcvLgth    /= 3.0;
            //-----------------------------------------------   
            GetPtOrderAndPlane(fP1r, fP2r, fP3r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv); //to make sure that fault normal will point upward
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2rOut, fP3rOut);   //to rotate into the coordinate system of the receiver!

            //-----------------------------------------------             
            fTemp0   = 0.0;           fTemp1 = 0.0;           fTemp2 = 0.0; 
            fTemp3   = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;
            fTemp6   = 0.0;           fTemp7 = 0.0;           fTemp8 = 0.0;
            //-------------------------------  
            fTestDist= GetDist_2Pnts(fSrcPt, fRcvPt);            

            if (fTestDist > 2.0*fMaxLgth)    //is "self" or if distance is far enough away to not need the extra distance check 
            {   //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, MD->fUnitSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp1 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp2 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, MD->fUnitSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp4 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp5 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
                //------------------------------------------------------------------------------------------   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2sOut, fP3sOut, 0.0, 0.0, MD->fUnitSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = fStressOut[0]/MD->fUnitSlip *1.0e-6;           fTemp7 = fStressOut[1]/MD->fUnitSlip *1.0e-6;           fTemp8 = fStressOut[2]/MD->fUnitSlip *1.0e-6; 
            }          
            else
            {   //------------------------------------------------------------------------------------------ 
                iFlagPair  = 0;
                iPntCntr = GetSamePtNum(iSrcTrigs, iRcvTrigs);
                if (iPntCntr == 0)              
                {   //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestDist = GetLotPunkt_PktEbene(fnRcv, fdRcv, fSrcPt, fTestPt); //TestPt is here the "lotpunkt"  i.e., the projection of SrcPt onto the triangle plane
                        iPntCntr  = PtLiesInTrig(fP1r, fP2rOut, fP3rOut, fnRcv, fdRcv, fAreaRcv, fTestPt); 
                        if ((iPntCntr > 0) && (fTestDist <= 1.0*fDistFrac*fRcvLgth))             {    iFlagPair = 1;      }
                        //----------------------
                        fTestDist = GetLotPunkt_PktEbene(fnSrc, fdSrc, fRcvPt, fTestPt);
                        iPntCntr  = PtLiesInTrig(fP1s, fP2sOut, fP3sOut, fnSrc, fdSrc, fAreaSrc, fTestPt); //now the "TestPt is projection point of RcvPt onto triangle plane of source
                        if ((iPntCntr > 0) && (fTestDist <= 1.0*fDistFrac*fSrcLgth))             {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                    if (iFlagPair == 0)
                    {   fTestFrac1 = GetMinFraction(fP1s, fP2sOut, fP3sOut, fRcvPt);
                        fTestFrac2 = GetMinFraction(fP1r, fP2rOut, fP3rOut, fSrcPt);
                        if ( (fTestFrac1 < fDistFrac) || (fTestFrac2 < fDistFrac) )     {    iFlagPair = 1;      }
                    }
                    //----------------------------------------------
                }
                //------------------------------------------------------------------------------------------                   
                if (iFlagPair == 0) 
                {   //-----------------------------------------------
                    for (k = 0; k < iTrigNum; k++)
                    {   
                        fCurFrac   = (float)(k+1)/(float)iTrigNum;
                        GetNewPoints(fCurFrac, fP1s, fP2sOut, fP3sOut, fP1n, fP2n, fP3n);
                        fTstArea   = GetArea(fP1n,fP2n,fP3n);
                        fPartlSlip = fPartlMom/fTstArea;
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, fPartlSlip, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp0 += fStressOut[0];        fTemp1 += fStressOut[1];           fTemp2 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, fPartlSlip, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp3 += fStressOut[0];        fTemp4 += fStressOut[1];           fTemp5 += fStressOut[2]; 
                        //-----------------------------------------------  
                        StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1n, fP2n, fP3n, 0.0, 0.0, fPartlSlip, MD->fShearMod, MD->fLambda);         // slip in stk         
                        //-----------------------------------------------  
                        RotateTensor_inKmatrix(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                        //----------------------------------------------- 
                        fTemp6 += fStressOut[0];        fTemp7 += fStressOut[1];           fTemp8 += fStressOut[2];                         
                    }
                    fTemp0 = fTemp0/MD->fUnitSlip *1.0e-6;          fTemp1 = fTemp1/MD->fUnitSlip *1.0e-6;          fTemp2 = fTemp2/MD->fUnitSlip *1.0e-6;
                    fTemp3 = fTemp3/MD->fUnitSlip *1.0e-6;          fTemp4 = fTemp4/MD->fUnitSlip *1.0e-6;          fTemp5 = fTemp5/MD->fUnitSlip *1.0e-6;
                    fTemp6 = fTemp6/MD->fUnitSlip *1.0e-6;          fTemp7 = fTemp7/MD->fUnitSlip *1.0e-6;          fTemp8 = fTemp8/MD->fUnitSlip *1.0e-6;                              
            }   }   
            //------------------------------------------------------------------------------------------        
            gsl_matrix_float_set(K->BF_SO, i, j, fTemp0);   gsl_matrix_float_set(K->BF_SS, i, j, fTemp1);   gsl_matrix_float_set(K->BF_SD, i, j, fTemp2);  
            gsl_matrix_float_set(K->BF_DO, i, j, fTemp3);   gsl_matrix_float_set(K->BF_DS, i, j, fTemp4);   gsl_matrix_float_set(K->BF_DD, i, j, fTemp5);  
            gsl_matrix_float_set(K->BF_OO, i, j, fTemp6);   gsl_matrix_float_set(K->BF_OS, i, j, fTemp7);   gsl_matrix_float_set(K->BF_OD, i, j, fTemp8);       
            //------------------------------------------------------------------------------------------        
        } 
    }
    //------------------------------------------------------------------   
    //------------------------------------------------------------------  
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++)
    {   // determine friction difference, use together with reference normal stress to get stress drop; get corresponding slip amount, compare with Dc => assign type         
        fTemp0 = gsl_matrix_float_get(K->FF_SS, TR->ivFL_SelfLoc[i], i);
		fTemp1 = gsl_matrix_float_get(K->FF_DD, TR->ivFL_SelfLoc[i], i);
        fTemp2 = (TR->fvFL_DynFric[i] - TR->fvFL_StaFric[i]) *-1.0*TR->fvFL_RefNrmStrs[i]; //multiply with -1.0 to get normal stress now to compressive positive
        fTemp3 = fTemp2/((fTemp0 > fTemp1) ? fTemp0 : fTemp1);

        if (fTemp3 > TR->fvFL_CurDcVal[i]) // have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc
        {        TR->ivFL_StabT[i] = 1; // patch is unstable
        }
        else
        {   if (fTemp2 < 0.0)  //  is still weakening but just with slip that is lower than Dc => cond. stable
            {    TR->ivFL_StabT[i] = 2; // patch is cond. stable
            }
            else // no weakening but strengthening  => stable/
            {    TR->ivFL_StabT[i] = 3; // patch is stable 
    }   }   }
    fprintf(stdout,"iRANK: %d flagged interactions: %d (fault)   and overruled %d \n",MD->iRANK, iCntFlags_F, iOverRuled);
    //------------------------------------------------------------------ 
    return;
} 
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void  GetPtOrderAndPlane(float fP1in[3], float fP2in[3], float fP3in[3], float fP2out[3], float fP3out[3], float fNrm[3], float fd[1], float fArea[1])
{   float fTemp;
    float fP1toP2[3],       fP1toP3[3];

    fP1toP2[0] = fP2in[0] - fP1in[0];           fP1toP2[1] = fP2in[1] - fP1in[1];               fP1toP2[2] = fP2in[2] - fP1in[2];
    fP1toP3[0] = fP3in[0] - fP1in[0];           fP1toP3[1] = fP3in[1] - fP1in[1];               fP1toP3[2] = fP3in[2] - fP1in[2];
 
    fNrm[0]    = fP1toP2[1]*fP1toP3[2] - fP1toP2[2]*fP1toP3[1];
    fNrm[1]    = fP1toP2[2]*fP1toP3[0] - fP1toP2[0]*fP1toP3[2];
    fNrm[2]    = fP1toP2[0]*fP1toP3[1] - fP1toP2[1]*fP1toP3[0];
    //----------------------
    fArea[0]   = 0.5*sqrtf( fNrm[0]*fNrm[0] +fNrm[1]*fNrm[1] +fNrm[2]*fNrm[2]); //get area of the triangle
    //----------------------
    fTemp      = sqrtf(fNrm[0]*fNrm[0] +fNrm[1]*fNrm[1] +fNrm[2]*fNrm[2]); //get length of normal vector => to normalize the normal vector so that it has a length of 1...
    fNrm[0]   /= fTemp;                         fNrm[1]   /= fTemp;                             fNrm[2]   /= fTemp;
    //----------------------
    if (fNrm[2] < 0.0) //flip the points and the normal vector (accordingly) if the normal vector is pointing downwards
    {   fP2out[0]      = fP3in[0];            fP2out[1]      = fP3in[1];             fP2out[2]      = fP3in[2];
        fP3out[0]      = fP2in[0];            fP3out[1]      = fP2in[1];             fP3out[2]      = fP2in[2];
        fNrm[0]       *= -1.0;                fNrm[1]       *= -1.0;                 fNrm[2]       *= -1.0; 
    }
    else //else, don't flip...
    {   fP2out[0]      = fP2in[0];            fP2out[1]      = fP2in[1];                fP2out[2]      = fP2in[2];
        fP3out[0]      = fP3in[0];            fP3out[1]      = fP3in[1];                fP3out[2]      = fP3in[2];
    }
    //----------------------
    fd[0]      = fNrm[0]*fP1in[0] + fNrm[1]*fP1in[1] +fNrm[2]*fP1in[2]; //part of the plane equation
    //----------------------
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3])
{   // this comes basically from Medhii's code... 
    float fTempfloat;
    float ftempVect1[3],         ftempVect2[3],             feZ[3];
    //-----------------------------------------------  
    feZ[0]        = 0.0;                       feZ[1]        = 0.0;                       feZ[2]        = 1.0;        
    ftempVect1[0] = fP2[0] -fP1[0];            ftempVect1[1] = fP2[1] -fP1[1];            ftempVect1[2] = fP2[2] -fP1[2];
    ftempVect2[0] = fP3[0] -fP1[0];            ftempVect2[1] = fP3[1] -fP1[1];            ftempVect2[2] = fP3[2] -fP1[2];
    //-----------------------------------------------  
    fvNrm[0]      = ftempVect1[1]*ftempVect2[2] - ftempVect1[2]*ftempVect2[1];
    fvNrm[1]      = ftempVect1[2]*ftempVect2[0] - ftempVect1[0]*ftempVect2[2];
    fvNrm[2]      = ftempVect1[0]*ftempVect2[1] - ftempVect1[1]*ftempVect2[0];
    fTempfloat    = sqrtf(fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]);
    fvNrm[0]      = fvNrm[0]/fTempfloat;      fvNrm[1]     = fvNrm[1]/fTempfloat;      fvNrm[2]     = fvNrm[2]/fTempfloat;
    
    if (fvNrm[2] < 0.0)     {   fvNrm[0]       *= -1.0;                fvNrm[1]       *= -1.0;                 fvNrm[2]       *= -1.0;      }
    //-----------------------------------------------  
    fvStk[0]      = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    fvStk[1]      = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    fvStk[2]      = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) 
    fTempfloat    = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    if (fabs(fTempfloat) < FLT_EPSILON)
    {   fvStk[0]    = 0.0;                 fvStk[1] = 1.0;                 fvStk[2] = 0.0;
        fTempfloat  = 1.0;        
    }
    fvStk[0]    = fvStk[0]/fTempfloat;      fvStk[1] = fvStk[1]/fTempfloat;         fvStk[2] = fvStk[2]/fTempfloat;
    //-----------------------------------------------  
    fvDip[0]    = fvNrm[1]*fvStk[2] - fvNrm[2]*fvStk[1];
    fvDip[1]    = fvNrm[2]*fvStk[0] - fvNrm[0]*fvStk[2];
    fvDip[2]    = fvNrm[0]*fvStk[1] - fvNrm[1]*fvStk[0];
    fTempfloat  = sqrtf(fvDip[0]*fvDip[0] +fvDip[1]*fvDip[1] +fvDip[2]*fvDip[2]);
    fvDip[0]    = fvDip[0]/fTempfloat;      fvDip[1] = fvDip[1]/fTempfloat;         fvDip[2] = fvDip[2]/fTempfloat;
    //-----------------------------------------------  
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetDist_2Pnts(float fT[3], float fL[3])
{   float TheDist;
    TheDist   = sqrtf( (fT[0]-fL[0])*(fT[0]-fL[0]) + (fT[1]-fL[1])*(fT[1]-fL[1]) + (fT[2]-fL[2])*(fT[2]-fL[2]) ); 
    return TheDist;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void RotateTensor_inKmatrix(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir)
{   
    float fTxx1,          fTxy1,         fTxz1,          fTyy1,         fTyz1,            fTzz1;
    float fTxx2,          fTxy2,         fTxz2,          fTyy2,         fTyz2,            fTzz2;
    float fA[9]; // the "Rotation Matrix" 
    //-----------------------------------------------  
    fTxx1       = fsig[0];             fTxy1       = fsig[1];            fTxz1       = fsig[2];
    fTyy1       = fsig[3];             fTyz1       = fsig[4];            fTzz1       = fsig[5];
    //-----------------------------------------------  
    // TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate system to x2y2z2 coordinate system. "A" is the transformation matrix, 
    // whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The  coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose  of A (i.e., A') does the transformation from x2y2z2 into x1y1z1. 
    if (iRotDir == 1) // from local to global 
    {   fA[0] = fvNrm[0];         fA[3] = fvStk[0];        fA[6] = fvDip[0];                
        fA[1] = fvNrm[1];         fA[4] = fvStk[1];        fA[7] = fvDip[1];
        fA[2] = fvNrm[2];         fA[5] = fvStk[2];        fA[8] = fvDip[2];
    }
    else             // from global to local 
    {   fA[0] = fvNrm[0];         fA[3] = fvNrm[1];        fA[6] = fvNrm[2];
        fA[1] = fvStk[0];         fA[4] = fvStk[1];        fA[7] = fvStk[2];
        fA[2] = fvDip[0];         fA[5] = fvDip[1];        fA[8] = fvDip[2];
    }
    //-----------------------------------------------  
    fTxx2 = fA[0]*fA[0]*fTxx1 +           2.0*fA[0]*fA[3] *fTxy1 +            2.0*fA[0]*fA[6] *fTxz1 +            2.0*fA[3]*fA[6] *fTyz1 + fA[3]*fA[3]*fTyy1 + fA[6]*fA[6]*fTzz1;
    fTyy2 = fA[1]*fA[1]*fTxx1 +           2.0*fA[1]*fA[4] *fTxy1 +            2.0*fA[1]*fA[7] *fTxz1 +            2.0*fA[4]*fA[7] *fTyz1 + fA[4]*fA[4]*fTyy1 + fA[7]*fA[7]*fTzz1;
    fTzz2 = fA[2]*fA[2]*fTxx1 +           2.0*fA[2]*fA[5] *fTxy1 +            2.0*fA[2]*fA[8] *fTxz1 +            2.0*fA[5]*fA[8] *fTyz1 + fA[5]*fA[5]*fTyy1 + fA[8]*fA[8]*fTzz1;  
    fTxy2 = fA[0]*fA[1]*fTxx1 + (fA[0]*fA[4]+ fA[1]*fA[3])*fTxy1 + (fA[0]*fA[7] + fA[1]*fA[6])*fTxz1 + (fA[7]*fA[3] + fA[6]*fA[4])*fTyz1 + fA[4]*fA[3]*fTyy1 + fA[6]*fA[7]*fTzz1;
    fTxz2 = fA[0]*fA[2]*fTxx1 + (fA[0]*fA[5]+ fA[2]*fA[3])*fTxy1 + (fA[0]*fA[8] + fA[2]*fA[6])*fTxz1 + (fA[8]*fA[3] + fA[6]*fA[5])*fTyz1 + fA[5]*fA[3]*fTyy1 + fA[6]*fA[8]*fTzz1;
    fTyz2 = fA[1]*fA[2]*fTxx1 + (fA[2]*fA[4]+ fA[1]*fA[5])*fTxy1 + (fA[2]*fA[7] + fA[1]*fA[8])*fTxz1 + (fA[7]*fA[5] + fA[8]*fA[4])*fTyz1 + fA[4]*fA[5]*fTyy1 + fA[7]*fA[8]*fTzz1;
    //-----------------------------------------------  
    fsig_new[0] = fTxx2;        fsig_new[1] = fTxy2;      fsig_new[2] = fTxz2;
    fsig_new[3] = fTyy2;        fsig_new[4] = fTyz2;      fsig_new[5] = fTzz2;
 
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int GetSamePtNum(int iSrcTrigs[3], int iRcvTrigs[3])
{   int iSamePtNum = 0;

    iSamePtNum = (iSrcTrigs[0] == iRcvTrigs[0]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[0] == iRcvTrigs[1]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[0] == iRcvTrigs[2]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[1] == iRcvTrigs[0]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[1] == iRcvTrigs[1]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[1] == iRcvTrigs[2]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[2] == iRcvTrigs[0]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[2] == iRcvTrigs[1]) ? iSamePtNum+1 : iSamePtNum;
    iSamePtNum = (iSrcTrigs[2] == iRcvTrigs[2]) ? iSamePtNum+1 : iSamePtNum;
    return iSamePtNum;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int  PtLiesInTrig(float fP1[3], float fP2[3], float fP3[3], float fNrm[3], float fd[1], float fArea[1], float fT[3])
{   int IsInPlane = 0;
    float fArea1,       fArea2,         fArea3,           fTestVal;

    fTestVal = fNrm[0]*fT[0] + fNrm[1]*fT[1] +fNrm[2]*fT[2];

    if (fTestVal == fd[0]) //test point is in the same plane that the triangle forms
    {   fArea1  = GetArea(fP1, fP2, fT);        fArea2  = GetArea(fP2, fP3, fT);                  fArea3  = GetArea(fP3, fP1, fT);
        
        if ((fArea1+fArea2+fArea3) <= fArea[0])
        {   IsInPlane = 2;
            if ((fArea1 == 0.0) || (fArea2 == 0.0) || (fArea3 == 0.0))      
            {   IsInPlane = 1;
    }   }   }

    return IsInPlane;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetArea(float fP1[3], float fP2[3], float fP3[3])
{   float fP1toP2[3],       fP1toP3[3],      fNrm[3];
    float fArea0;
 
    fP1toP2[0] = fP2[0] - fP1[0];              fP1toP2[1] = fP2[1] - fP1[1];               fP1toP2[2] = fP2[2] - fP1[2];
    fP1toP3[0] = fP3[0] - fP1[0];              fP1toP3[1] = fP3[1] - fP1[1];               fP1toP3[2] = fP3[2] - fP1[2];
    fNrm[0]    = fP1toP2[1]*fP1toP3[2] - fP1toP2[2]*fP1toP3[1];
    fNrm[1]    = fP1toP2[2]*fP1toP3[0] - fP1toP2[0]*fP1toP3[2];
    fNrm[2]    = fP1toP2[0]*fP1toP3[1] - fP1toP2[1]*fP1toP3[0];
	fArea0     = 0.5*sqrtf( fNrm[0]*fNrm[0] +fNrm[1]*fNrm[1] +fNrm[2]*fNrm[2]);
    return fArea0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetLotPunkt_PktEbene(float fNrm[3], float fd[1], float fT[3], float fL[3])
{
    float s,   dist;
    s          = (fd[0] -fNrm[0]*fT[0] -fNrm[1]*fT[1] -fNrm[2]*fT[2]) / (fNrm[0]*fNrm[0] +fNrm[1]*fNrm[1] +fNrm[2]*fNrm[2]); 
    fL[0]      = fT[0] +s*fNrm[0];
    fL[1]      = fT[1] +s*fNrm[1];
    fL[2]      = fT[2] +s*fNrm[2];
    dist       = sqrtf((fT[0]-fL[0])*(fT[0]-fL[0]) + (fT[1]-fL[1])*(fT[1]-fL[1]) + (fT[2]-fL[2])*(fT[2]-fL[2]) );
    return dist;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
 float GetMinFraction(float fP1s[3], float fP2s[3], float fP3s[3], float fT[3])
 {  //https://www.mathematik-oberstufe.de/vektoren/a/abstand-punkt-gerade-lot.html
    float  fDist0 = 0.0,    fDist1,           fDist2,         fDist3,            fDist4;
    float  fThisFrac = INFINITY; //just a starting value
    float  fL[3];
    //-------------------------------------------------------------
    fDist0   += GetDist_2Pnts(fP1s, fP2s);          fDist0   += GetDist_2Pnts(fP2s, fP3s);
    fDist0   += GetDist_2Pnts(fP3s, fP1s);          fDist0   /= 3.0;
    //first, test the distance of receiver point to all the vertices (P1, P2, P3); 
    fDist4    = GetDist_2Pnts(fP1s, fT);
    fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac;    
    fDist4    = GetDist_2Pnts(fP2s, fT);
    fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac;  
    fDist4    = GetDist_2Pnts(fP3s, fT);
    fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac;
    //-------------------------------------------------------------
    //for P2-P1 vector i.e., triangle leg
    fDist1    = GetDist_2Pnts(fP1s, fP2s);
    GetLotPunkt_PktGerade(fP1s, fP2s, fT, fL);
    fDist2    = GetDist_2Pnts(fL, fP1s);
    fDist3    = GetDist_2Pnts(fL, fP2s);
    fDist4    = GetDist_2Pnts(fL, fT);
    if ((fDist2 <= fDist1) && (fDist3 <= fDist1)) /*if both Schnittpunkt-Point lengths are equal or smaller than Point to Point length, then the Schnittpunkt lies between them and I have to check the distance*/
    {   fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac;  
    }
    //for P2-P3 vector i.e., triangle leg
    fDist1    = GetDist_2Pnts(fP2s, fP3s);
    GetLotPunkt_PktGerade(fP2s, fP3s, fT, fL);
    fDist2    = GetDist_2Pnts(fL, fP2s);
    fDist3    = GetDist_2Pnts(fL, fP3s);
    fDist4    = GetDist_2Pnts(fL, fT);
    if ((fDist2 <= fDist1) && (fDist3 <= fDist1)) /*if both Schnittpunkt-Point lengths are equal or smaller than Point to Point length, then the Schnittpunkt lies between them and I have to check the distance*/
    {   fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac; 
    }
    //for P3-P1 vector i.e., triangle leg
    fDist1    = GetDist_2Pnts(fP3s, fP1s);
    GetLotPunkt_PktGerade(fP3s, fP1s, fT, fL);
    fDist2    = GetDist_2Pnts(fL, fP3s);
    fDist3    = GetDist_2Pnts(fL, fP1s);
    fDist4    = GetDist_2Pnts(fL, fT);
    if ((fDist2 <= fDist1) && (fDist3 <= fDist1)) /*if both Schnittpunkt-Point lengths are equal or smaller than Point to Point length, then the Schnittpunkt lies between them and I have to check the distance*/
    {   fThisFrac = (fDist4/fDist0 < fThisFrac) ? fDist4/fDist0 : fThisFrac;   
    }
    //-------------------------------------------------------------
    return fThisFrac;
 }
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetLotPunkt_PktGerade(float fP1[3], float fP2[3], float fT[3], float fL[3])
{   float fP1toP2[3];
    float s;
    
    fP1toP2[0] = fP2[0] - fP1[0];              fP1toP2[1] = fP2[1] - fP1[1];               fP1toP2[2] = fP2[2] - fP1[2];
    s          = ( -(fP1[0]-fT[0])*fP1toP2[0] -(fP1[1]-fT[1])*fP1toP2[1] -(fP1[2]-fT[2])*fP1toP2[2]) / (fP1toP2[0]*fP1toP2[0] +fP1toP2[1]*fP1toP2[1] +fP1toP2[2]*fP1toP2[2]);
    fL[0]      = fP1[0] +s*fP1toP2[0];
    fL[1]      = fP1[1] +s*fP1toP2[1];
    fL[2]      = fP1[2] +s*fP1toP2[2];

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetNewPoints(float CurFrac, float fP1s[3], float fP2s[3], float fP3s[3], float fP1n[3], float fP2n[3], float fP3n[3])
{   float Cent[3],  dCent2P1[3], dCent2P2[3], dCent2P3[3];

    Cent[0]  = (fP1s[0]+fP2s[0]+fP3s[0])/3.0;   Cent[1]  = (fP1s[1]+fP2s[1]+fP3s[1])/3.0;       Cent[2]  = (fP1s[2]+fP2s[2]+fP3s[2])/3.0;

    dCent2P1[0]= fP1s[0] - Cent[0];             dCent2P1[1]= fP1s[1] - Cent[1];                 dCent2P1[2]= fP1s[2] - Cent[2];
    dCent2P2[0]= fP2s[0] - Cent[0];             dCent2P2[1]= fP2s[1] - Cent[1];                 dCent2P2[2]= fP2s[2] - Cent[2];
    dCent2P3[0]= fP3s[0] - Cent[0];             dCent2P3[1]= fP3s[1] - Cent[1];                 dCent2P3[2]= fP3s[2] - Cent[2];

    fP1n[0] = Cent[0]+CurFrac*dCent2P1[0];      fP1n[1] = Cent[1]+CurFrac*dCent2P1[1];          fP1n[2] = Cent[2]+CurFrac*dCent2P1[2];
    fP2n[0] = Cent[0]+CurFrac*dCent2P2[0];      fP2n[1] = Cent[1]+CurFrac*dCent2P2[1];          fP2n[2] = Cent[2]+CurFrac*dCent2P2[2]; 
    fP3n[0] = Cent[0]+CurFrac*dCent2P3[0];      fP3n[1] = Cent[1]+CurFrac*dCent2P3[1];          fP3n[2] = Cent[2]+CurFrac*dCent2P3[2];

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//