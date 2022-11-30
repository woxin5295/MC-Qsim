#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
struct MDstruct
{   int         iRANK,              		iSIZE;
    int       	iF_BASEelem,        		iF_ADDelem,				iB_BASEelem,       			iB_ADDelem;  
    int       	*ivF_OFFSET,        		*ivF_START,				*ivB_OFFSET,        		*ivB_START;
	int         *ivF_ModOFFs,				*ivB_ModOFFs;

    int       	iRunNum,            		iSeedStart,				iUseProp,           		iUsePSeis;
    int       	iFPNum,             		iBPNum,					iFVNum,             		iBVNum; 
    int       	iFSegNum,          		 	iBSegNum,				iAlsoBoundForPostSeis;
    int       	iMaxIterat,      			iEQcntr,      			iMaxMRFlgth,				iMaxSTFlgth;				
    int      	iChgBtwEQs,					iGlobTTmax,				iWritePos,					iUseVpVs;
   
    float     	fLegLgth,          			fUnitSlip,				fDeltT;
	float 		fFltLegs,					fBndLegs,               fMinMag4Prop,				fCutStFrac;

	float       fAftrSlipTime,				fISeisStep,         	fRecLgth,					fHealFact;
	float       fDeepRelaxTime,				fPSeis_Step;

    float       fMedDense,     				fAddNrmStrs,			fVp,						fVs;
    float       fPoisson,           		fLambda,        		fShearMod,					fMeanStiffness;					
    float       fISeisTStp,         		fVpVsRatio,				fDcfrac4Threshold;	
    long long   lTimeYears,					lRecLgth;
    
    char        cInputName[512];
};
//------------------------------------------------------------------
struct TRstruct
{	int     	*ivFL_StabT,				*ivFL_Activated,		*ivFL_Ptch_t0,				*ivFL_FricLaw;	

    float    	*fvFL_RefNrmStrs,			*fvFL_Area;
	float		*fvFL_SelfStiffStk,			*fvFL_SelfStiffDip,		*fvFL_MeanSelfStiff;
	float		*fvBL_SelfStiffStk,			*fvBL_SelfStiffDip,		*fvBL_SelfStiffOpn;

    float    	*fvFL_RefStaFric,           *fvFL_RefDynFric,    	*fvFL_RefDcVal;
    float    	*fvFL_RefStaFric_vari,      *fvFL_RefDynFric_vari,	*fvFL_RefDcVal_vari;
    float    	*fvFL_StaFric,              *fvFL_DynFric,			*fvFL_CurFric;              
	float       *fvFL_B4_Fric,			    *fvFL_TempRefFric,		*fvFL_CurDcVal;

    float    	*fvFL_PSeis_T0_F,          	*fvFL_AccumSlp,			*fvFL_CutStress;
	float       *fvBL_PSeis_T0_S,			*fvBL_PSeis_T0_N;
	float    	*fvFL_SlipRate_temp,      	*fvFL_SlipRake_temp;
	 
	int      	*ivFG_SegID_temp,			*ivFG_FltID_temp,		*ivFG_Flagged_temp;
    int    		*ivFG_V1_temp,              *ivFG_V2_temp,        	*ivFG_V3_temp;
    float       *fvFG_StressRatetemp,		*fvFG_SlipRatetemp,     *fvFG_Raketemp;
    float       *fvFG_MaxTransient;
    	
    int    		*ivBG_V1_temp,              *ivBG_V2_temp,       	*ivBG_V3_temp,				 *ivBG_SegID_temp;
    float    	*fvFG_CentE_temp,          	*fvFG_CentN_temp,    	*fvFG_CentZ_temp;
    float    	*fvBG_CentE_temp,           *fvBG_CentN_temp,    	*fvBG_CentZ_temp;    
    float       *fvFL_StaFricMod_temp, 		*fvFL_DynFricMod_temp,	*fvFL_NrmStrsMod_temp,	     *fvFL_DcMod_temp;

	gsl_matrix_int	    *imFGL_TTP,         *imFGL_NextP;
	gsl_matrix_int    	*imFGL_TTS,			*imFGL_NextS; 
    gsl_matrix_float    *fmFGL_SrcRcvH,    	*fmFGL_SrcRcvV,      	*fmFGL_SrcRcvN;    

	gsl_vector_float	*fvFL_StrsRateStk, 	*fvFL_StrsRateDip;
	gsl_vector_float	*fvFL_CurStrsH,		*fvFL_CurStrsV,			*fvFL_CurStrsN;
	gsl_vector_float    *fvFL_B4_StrsH,		*fvFL_B4_StrsV,		    *fvFL_B4_StrsN;    
	gsl_vector_float	*fvBL_CurStrsH,		*fvBL_CurStrsV,			*fvBL_CurStrsN;
};
//------------------------------------------------------------------
struct VTstruct
{	float   	*fvFG_VlX_temp,             *fvFG_VlY_temp;
	float  		*fvFG_PosE_temp,          	*fvFG_PosN_temp,      	*fvFG_PosZ_temp,			*fvFG_Hght_temp;
	float   	*fvBG_PosE_temp,           	*fvBG_PosN_temp,      	*fvBG_PosZ_temp;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float   fMinOf2(float val1, float val2);   
float   fMaxOf2(float val1, float val2);   
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void LoadInput(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT)
{	
    //------------------------------------------------------------------
	gsl_rng *fRandN; // this is pretty much straight from the GSL reference, the default RNG has good performance, so no need to change
    const gsl_rng_type *RandT;
	gsl_rng_env_setup(); 
	RandT   = gsl_rng_default;
    fRandN  = gsl_rng_alloc(RandT);
	unsigned long RSeed =(unsigned long)(MD->iRANK + MD->iSeedStart);   
    gsl_rng_set(fRandN, RSeed);
	//-----------------------------------------------------------------
	int		i,            			iGlobPos,				iTemp0;
    float   fTemp0,					fTemp1,					fTemp2,				fTemp3;
    char    ctempVals[512],         cAppend[512];
	char    cFileName1[512], 		cFileName2[512],		cFileName3[512],	cFileName4[512];
    FILE    *fp1,               	*fp2,           		*fp3; 
	//-------------------------------------
	//reading/loading the remaining information from txt and dat files.. this step is done by each rank individually; also, some values are read locally (every rank only reads the information relevant for it) while other values are read globally
    strcpy(cFileName1,MD->cInputName);      strcat(cFileName1,"_Summary_RoughStrength.txt"); 
	strcpy(cFileName2,MD->cInputName);      strcat(cFileName2,"_");               	sprintf(cAppend, "%d",MD->iRunNum); 	strcat(cFileName2,cAppend);     	strcat(cFileName2,"_Roughn.dat");
	strcpy(cFileName3,MD->cInputName);      strcat(cFileName3,"_");             	sprintf(cAppend, "%d",MD->iRunNum); 	strcat(cFileName3,cAppend);     	strcat(cFileName3,"_Strgth.dat"); 
	strcpy(cFileName4,MD->cInputName);      strcat(cFileName4,"_BNDtrig.dat");
    //-----------------------------------------------------------------
    if ((fp1 = fopen(cFileName1,"r")) == NULL)          { 	printf("Error -cant open *.flt file. LoadInputParameter function...\n");      exit(10);     }
    
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fMedDense);                                          }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fAddNrmStrs);    MD->fAddNrmStrs *= -1.0;            }// keep at MPa but make compression negative ==> in MATLAB and the input files, normal stress is still positive for compression => is switched here to compression == negative convention
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fShearMod);      MD->fShearMod   *=  1.0E+9;         }// convert from GPa to Pa (for actual computation of K-matrix and also for magnitude calculation, I need Pa and not MPa)
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fPoisson);                                           }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %d",&MD->iChgBtwEQs);                                         }
     
    MD->fLambda = (2.0*MD->fShearMod*MD->fPoisson)/(1.0-2.0*MD->fPoisson); //elastic parameters => for K-matrix calculation, remember that the code uses a linear elastic half-space
    fTemp0      = (MD->fMedDense > 0.0) ? MD->fMedDense : 2700.0; // if zero density is used for depth gradient (e.g., to have a depth-independent normal stress) then I still!! need to define a density for the wave propagation speed.... => is done here...
    MD->fVp     = sqrtf((MD->fLambda +2.0*MD->fShearMod)/fTemp0); // in m/s
    MD->fVs     = sqrtf(MD->fShearMod/fTemp0); // in m/s 

	//-----------------------------------------------------------------
	if ((fp2 = fopen(cFileName2,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName2);      exit(10);     }
    if ((fp3 = fopen(cFileName3,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName3);      exit(10);     }    
	//-----------------------------------------------------------------
	if (fread(&iTemp0, sizeof(int),  1,fp2) != 1)  		{	exit(10);	}// have those already 
    if (fread(&iTemp0, sizeof(int),  1,fp2) != 1) 		{	exit(10);	}// have those already 
    if (fread(&fTemp0, sizeof(float),1,fp2) != 1)  		{	exit(10);	}// have those already 
    if (fread(&fTemp0, sizeof(float),1,fp2) != 1) 		{	exit(10);	}// have those already 
 	//-----------------------------------------------------
	if (fread(TR->ivFG_V1_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	} //reading the geometric information about the faults
    if (fread(TR->ivFG_V2_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->ivFG_V3_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
	if (fread(TR->ivFG_SegID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->ivFG_FltID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)	    {	exit(10);	}    
    
    if (fread(TR->fvFG_StressRatetemp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)	{	exit(10);	}    
    if (fread(TR->fvFG_SlipRatetemp,   sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)	{	exit(10);	}    
    if (fread(TR->fvFG_Raketemp,       sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)	{	exit(10);	}    
 
	//-------------------------------------
	if (fread(VT->fvFG_VlX_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)	    {	exit(10);	}
	if (fread(VT->fvFG_VlY_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)	    {	exit(10);	}
    //-------------------------------------
    if (fread(VT->fvFG_PosE_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    if (fread(VT->fvFG_PosN_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    if (fread(VT->fvFG_PosZ_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    if (fread(VT->fvFG_Hght_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)		{	exit(10);	}
    //----------------------------------------
    if (fread(TR->fvFG_CentE_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
    if (fread(TR->fvFG_CentN_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}   
    if (fread(TR->fvFG_CentZ_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)		{	exit(10);	}
	//------------------------------------- 	
    fclose(fp2);
    //-----------------------------------------------------------------
    fseek(fp3, (1L*sizeof(int)),       SEEK_CUR); // contains the patch number again -> skipped
    
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefStaFric,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDynFric,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefNrmStrs,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefDcVal,        sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	//------------------------------------- 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefStaFric_vari, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDynFric_vari, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_RefDcVal_vari,   sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	//------------------------------------- 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_StaFricMod_temp, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_DynFricMod_temp, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])	{	exit(10);	}
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_NrmStrsMod_temp, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	
	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->fvFL_DcMod_temp,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	

	fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
	if (fread(TR->ivFL_FricLaw,         sizeof( int),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])			{	exit(10);	}
	fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
	//-------------------------------------
    fclose(fp3);
    //-----------------------------------------------------------------
	if (MD->iBPNum > 0)
    {   if ((fp1 = fopen(cFileName4,"rb")) == NULL)  
     	{   printf("Warning -cant open  %s LoadInputParameter function...\n Continue with boundary surface",cFileName4);         
     		MD->iBPNum = 0;
     		MD->iBVNum = 0;
    	}
		else
		{	if (fread(&iTemp0, sizeof(int),  1,fp1) != 1)  		{	exit(10);	}// have those already 
    		if (fread(&iTemp0, sizeof(int),  1,fp1) != 1) 		{	exit(10);	}// have those already 
    		if (fread(&fTemp0, sizeof(float),1,fp1) != 1)  		{	exit(10);	}// have those already 
    		if (fread(&fTemp0, sizeof(float),1,fp1) != 1) 		{	exit(10);	}// have those already 
        	//-------------------------------------
			if (fread(TR->ivBG_V1_temp,    sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
    		if (fread(TR->ivBG_V2_temp,    sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	} 
    		if (fread(TR->ivBG_V3_temp,    sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	} 
    		if (fread(TR->ivBG_SegID_temp, sizeof(int), MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	} 
        	//-------------------------------------
        	fseek(fp1, 2L*sizeof(float)*(long)MD->iBVNum,   SEEK_CUR); // this would be local coordintates (from gridding) of the vertices..; not needed    
        	//-------------------------------------
	    	if (fread(VT->fvBG_PosE_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    		if (fread(VT->fvBG_PosN_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    		if (fread(VT->fvBG_PosZ_temp,  sizeof(float),MD->iBVNum,fp1) != MD->iBVNum)	{	exit(10);	}
    		//----------------------------------------
    		if (fread(TR->fvBG_CentE_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
    		if (fread(TR->fvBG_CentN_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}    
    		if (fread(TR->fvBG_CentZ_temp, sizeof(float),MD->iBPNum,fp1) != MD->iBPNum)	{	exit(10);	}
			//-------------------------------------  
        	fclose(fp1);   
    }	}
	//-----------------------------------------------------------------
    for (i = 0; i < (MD->iFVNum); i++)        //converting all km to meters, also shifting the indizes => Matlab starting with "1" while C starting with "0" => hence the -1 used here....    
    {   VT->fvFG_PosE_temp[i]     *= 1000.0;            VT->fvFG_PosN_temp[i]     *= 1000.0;         VT->fvFG_PosZ_temp[i]     *= 1000.0; 				VT->fvFG_Hght_temp[i]     *= 1000.0; 
    }
	for (i = 0; i < (MD->iBVNum); i++)            
    {   VT->fvBG_PosE_temp[i]     *= 1000.0;            VT->fvBG_PosN_temp[i]     *= 1000.0;         VT->fvBG_PosZ_temp[i]     *= 1000.0; 
    }
	for (i = 0; i < (MD->iFPNum); i++)            
    {	TR->ivFG_V1_temp[i]       -= 1;             	TR->ivFG_V2_temp[i]       -= 1;              TR->ivFG_V3_temp[i]       -= 1;                  
        TR->fvFG_CentE_temp[i]    *= 1000.0;         	TR->fvFG_CentN_temp[i]    *= 1000.0;         TR->fvFG_CentZ_temp[i]    *= 1000.0; 
        TR->ivFG_SegID_temp[i]    -= 1; 
	}
	for (i = 0; i < (MD->iBPNum); i++)            
    {	TR->ivBG_V1_temp[i]        -= 1;                TR->ivBG_V2_temp[i]       -= 1;              TR->ivBG_V3_temp[i]       -= 1;                  
        TR->fvBG_CentE_temp[i]  *= 1000.0;         		TR->fvBG_CentN_temp[i]    *= 1000.0;         TR->fvBG_CentZ_temp[i]    *= 1000.0; 
	}
	//-----------------------------------------------------------------

    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) //here the current static/dynamic friction coefficients are determines, also the current Dc value; further, the stressing rates are assigned
    {   
		iGlobPos                    = i + MD->ivF_START[MD->iRANK];
   		fTemp0                      = cosf(TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0) *TR->fvFG_StressRatetemp[iGlobPos];
   		fTemp1                      = sinf(TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0) *TR->fvFG_StressRatetemp[iGlobPos];
        gsl_vector_float_set(TR->fvFL_StrsRateStk, i, fTemp0); 
        gsl_vector_float_set(TR->fvFL_StrsRateDip, i, fTemp1);  
    
		TR->fvFL_SlipRate_temp[i]   = TR->fvFG_SlipRatetemp[iGlobPos]/1000; //was defined in mm/yr, butneed in m/y => /1000;
        TR->fvFL_SlipRake_temp[i]   = TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0; //rake angle in radian
  		//------------------------//------------------------
  		TR->fvFL_RefDcVal_vari[i]   = (TR->fvFL_RefDcVal[i]  *TR->fvFL_RefDcVal_vari[i]  /100.0) / (TR->fvFL_RefDcVal[i]  +TR->fvFL_DcMod_temp[i]);
  		TR->fvFL_RefStaFric_vari[i] = (TR->fvFL_RefStaFric[i]*TR->fvFL_RefStaFric_vari[i]/100.0) / (TR->fvFL_RefStaFric[i]+TR->fvFL_StaFricMod_temp[i]);
  		TR->fvFL_RefDynFric_vari[i] = (TR->fvFL_RefDynFric[i]*TR->fvFL_RefDynFric_vari[i]/100.0) / (TR->fvFL_RefDynFric[i]+TR->fvFL_DynFricMod_temp[i]);
  
  		TR->fvFL_RefStaFric[i]      = TR->fvFL_RefStaFric[i] + TR->fvFL_StaFricMod_temp[i];
  		TR->fvFL_RefDynFric[i]      = TR->fvFL_RefDynFric[i] + TR->fvFL_DynFricMod_temp[i];
  		TR->fvFL_RefNrmStrs[i]      = TR->fvFL_RefNrmStrs[i] + TR->fvFL_NrmStrsMod_temp[i];
  		TR->fvFL_RefDcVal[i]        = TR->fvFL_RefDcVal[i]   + TR->fvFL_DcMod_temp[i];
    	//------------------------//------------------------	
		if (TR->ivFL_FricLaw[i] == 1)		{	TR->fvFL_RefDcVal[i] = 0.0;		TR->fvFL_RefDcVal_vari[i]  = 0.0;			}
		//------------------------
        fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
		TR->fvFL_CurDcVal[i]       = TR->fvFL_RefDcVal[i]  *(1.0 + TR->fvFL_RefDcVal_vari[i]*fTemp2);
        //------------------------
		fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        TR->fvFL_StaFric[i]        = TR->fvFL_RefStaFric[i]*(1.0 + TR->fvFL_RefStaFric_vari[i]*fTemp2);     
        TR->fvFL_CurFric[i]        = TR->fvFL_StaFric[i];
		//------------------------
		fTemp1                     = (TR->fvFL_RefStaFric[i] - TR->fvFL_RefDynFric[i])/TR->fvFL_RefStaFric[i]; //reference friction change as fraction of static coefficient
		fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        
		fTemp3                     = fTemp1*(1.0 + TR->fvFL_RefDynFric_vari[i]*fTemp2); //calculation of dyn friction coefficient looks a bit more complicated, this is because this value is defined relative to the static coefficent i.e., as a fraction change of that coefficient
        TR->fvFL_DynFric[i]        = TR->fvFL_StaFric[i]*(1.0 - fTemp3);
    }
    //-----------------------------------------------------------------

	if (MD->iUseVpVs == 1)           {   MD->fVpVsRatio = MD->fVp/MD->fVs;                 }			else                            {   MD->fVpVsRatio = 1.0;                			}              
	if (MD->iUseProp == 1)           {   MD->fDeltT     = 1.0*MD->fLegLgth/MD->fVp;        }			else                            {   MD->fDeltT     = FLT_MAX;      	     			}

    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void DefineMoreParas(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, int iPlot2Screen)
{	int 	i,				j,				iGlobPos;
	int    	iTemp0,			iTemp1;
	float	fTemp0,			fTemp1,			fTemp2,			fTemp3,			fTemp4,			fTemp5,			fTemp6;
	float   fTemp7,			fTemp8;
	float	fP1[3],			fP2[3],			fP3[3];
	float	fP1P2[3],		fP1P3[3],		fSrcRcvVect[3];
	float	fvNrm[3],		fvStk[3],		fvDip[3];	
	float   feZ[3];
	//--------------------------------------------------------------------
	for (i = 0; i < MD->ivF_OFFSET[MD->iRANK];  i++)
    {   iGlobPos = i + MD->ivF_START[MD->iRANK];	
		fP1[0]   = VT->fvFG_PosE_temp[TR->ivFG_V1_temp[iGlobPos]];			fP1[1]   = VT->fvFG_PosN_temp[TR->ivFG_V1_temp[iGlobPos]];			fP1[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V1_temp[iGlobPos]];
		fP2[0]   = VT->fvFG_PosE_temp[TR->ivFG_V2_temp[iGlobPos]];			fP2[1]   = VT->fvFG_PosN_temp[TR->ivFG_V2_temp[iGlobPos]];			fP2[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V2_temp[iGlobPos]];
		fP3[0]   = VT->fvFG_PosE_temp[TR->ivFG_V3_temp[iGlobPos]];			fP3[1]   = VT->fvFG_PosN_temp[TR->ivFG_V3_temp[iGlobPos]];			fP3[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V3_temp[iGlobPos]];
		fP1P2[0] = fP2[0] - fP1[0];                fP1P2[1] = fP2[1] - fP1[1];                fP1P2[2] = fP2[2] - fP1[2];
        fP1P3[0] = fP3[0] - fP1[0];                fP1P3[1] = fP3[1] - fP1[1];                fP1P3[2] = fP3[2] - fP1[2];
        fvNrm[0] = fP1P2[1]*fP1P3[2] - fP1P2[2]*fP1P3[1];
        fvNrm[1] = fP1P2[2]*fP1P3[0] - fP1P2[0]*fP1P3[2];
        fvNrm[2] = fP1P2[0]*fP1P3[1] - fP1P2[1]*fP1P3[0];
		//--------------------------------------------------------------------
		TR->fvFL_Area[i]       = 0.5*sqrtf( fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]); //is needed for moment/magnitue calculation
		//--------------------------------------------------------------------
		//determine local coordinate system for currently selected (local) fault patch
    	feZ[0] = 0.0;                       feZ[1] = 0.0;                               feZ[2] = 1.0;        
    	//-----------------------------------------  
    	fTemp0      = sqrtf(fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]);
    	fvNrm[0]    = fvNrm[0]/fTemp0; 	fvNrm[1] = fvNrm[1]/fTemp0;      		fvNrm[2] = fvNrm[2]/fTemp0;
    	//----------------------------------------- 
    	fvStk[0]    = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    	fvStk[1]    = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    	fvStk[2]    = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    	// For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) 
    	fTemp0  = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    	if (fabs(fTemp0) < FLT_EPSILON)
		{   fvStk[0]= 0.0;                 		fvStk[1] = 1.0;        				fvStk[2] = 0.0;             }
		else
		{	fvStk[0]    = fvStk[0]/fTemp0;      fvStk[1] = fvStk[1]/fTemp0;         fvStk[2] = fvStk[2]/fTemp0;	}
    	//----------------------------------------- 
    	fvDip[0]    = fvNrm[1]*fvStk[2] - fvNrm[2]*fvStk[1];
    	fvDip[1]    = fvNrm[2]*fvStk[0] - fvNrm[0]*fvStk[2];
    	fvDip[2]    = fvNrm[0]*fvStk[1] - fvNrm[1]*fvStk[0];
    	fTemp0      = sqrtf(fvDip[0]*fvDip[0] +fvDip[1]*fvDip[1] +fvDip[2]*fvDip[2]);
    	fvDip[0]    = fvDip[0]/fTemp0;      fvDip[1] = fvDip[1]/fTemp0;         fvDip[2] = fvDip[2]/fTemp0;
		//-------------------------------------------------------------------- 
		for (j = 0; j < MD->iFPNum; j++)
        {   
			fTemp7   = sqrtf( ( (TR->fvFG_CentE_temp[j]-TR->fvFG_CentE_temp[iGlobPos])*(TR->fvFG_CentE_temp[j]-TR->fvFG_CentE_temp[iGlobPos]) ) + ( (TR->fvFG_CentN_temp[j]-TR->fvFG_CentN_temp[iGlobPos])*(TR->fvFG_CentN_temp[j]-TR->fvFG_CentN_temp[iGlobPos]) ) + ( (TR->fvFG_CentZ_temp[j]-TR->fvFG_CentZ_temp[iGlobPos])*(TR->fvFG_CentZ_temp[j]-TR->fvFG_CentZ_temp[iGlobPos])) )/MD->fVp; // the distance between both 
            fTemp8   = fTemp7*MD->fVpVsRatio; //these are the travel times from current "source" to current "receiver"
 			//----------------------------------------- 
			iTemp0   = (int)(fTemp7/MD->fDeltT); //means it's rounding downm, cutting off whatever floating point contribution the travel times had
			iTemp1   = (int)(fTemp8/MD->fDeltT); 
            gsl_matrix_int_set(TR->imFGL_TTP, i, j, iTemp0);
			gsl_matrix_int_set(TR->imFGL_TTS, i, j, iTemp1);
            //----------------------------------------- 
            MD->iGlobTTmax   = (MD->iGlobTTmax > iTemp1) ? MD->iGlobTTmax : iTemp1;     
			//----------------------------------------- 
            fSrcRcvVect[0]   = TR->fvFG_CentE_temp[j] - TR->fvFG_CentE_temp[iGlobPos]; // this is vector from current patch (source, local) to receiver (global); this one points therefore from source to the receiver
            fSrcRcvVect[1]   = TR->fvFG_CentN_temp[j] - TR->fvFG_CentN_temp[iGlobPos]; //will be needed for vector projection => when using rupture propagation => which component of slip vector is in direct line towards receiver patch (mode II) and which component is perpendicular to that (mode III)
            fSrcRcvVect[2]   = TR->fvFG_CentZ_temp[j] - TR->fvFG_CentZ_temp[iGlobPos]; //this is therefore for rupture propagation and the transient signals that arrise when using different velocities for mode II and mode III
            
            if (i == iGlobPos)
            {  	gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, 0.0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, 0.0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, 0.0);
			}
            else
            {   fTemp0          = sqrtf(fSrcRcvVect[0]*fSrcRcvVect[0] +fSrcRcvVect[1]*fSrcRcvVect[1] +fSrcRcvVect[2]*fSrcRcvVect[2]);
                fSrcRcvVect[0] /= fTemp0;                fSrcRcvVect[1] /= fTemp0;                fSrcRcvVect[2] /= fTemp0; 
				fTemp0	        = fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2]; //this rotates the vector from current source to receiver= fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2];
				fTemp1          = fvStk[0]*fSrcRcvVect[0] + fvStk[1]*fSrcRcvVect[1] + fvStk[2]*fSrcRcvVect[2]; //into local coordinate system of the source patch(that is the idea...)
				fTemp2		    = fvDip[0]*fSrcRcvVect[0] + fvDip[1]*fSrcRcvVect[1] + fvDip[2]*fSrcRcvVect[2];         

                gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, fTemp0);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, fTemp1);
				gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, fTemp2);       
				//----------------------------------------- 
				fTemp3 = sqrtf(fTemp1*fTemp1 +fTemp2*fTemp2);
				fTemp4 = fMinOf2(fabs(fTemp1),fabs(fTemp2)) / fMaxOf2(fabs(fTemp1),fabs(fTemp2));
				fTemp5 = fMinOf2(fabs(fTemp0),fabs(fTemp3)) / fMaxOf2(fabs(fTemp0),fabs(fTemp3));
				fTemp6 = sqrtf((fTemp4*fTemp4 + fTemp5*fTemp5))*(fTemp8-fTemp7);
				TR->fvFG_MaxTransient[j] = (fTemp6 > TR->fvFG_MaxTransient[j]) ? fTemp6 : TR->fvFG_MaxTransient[j];			
				//----------------------------------------- 	 
    }   }	}

	if ((iPlot2Screen == 1) && (MD->iRANK == 0))
    {	fprintf(stdout,"Medium density:              %5.2f\n",MD->fMedDense);
		fprintf(stdout,"Shear modulus:               %5.2e\n",MD->fShearMod);
 		fprintf(stdout,"Poisson ratio:               %5.2f\n",MD->fPoisson);
 		fprintf(stdout,"Lambda:                      %5.2e\n",MD->fLambda);
		fprintf(stdout,"P-wave velocity:             %5.2f\n",MD->fVp);
		fprintf(stdout,"S-wave velocity:             %5.2f\n",MD->fVs);
		fprintf(stdout,"Added normal stress (MPa):   %5.2f\n",MD->fAddNrmStrs);
		fprintf(stdout,"Change friction between EQs: %d\n",MD->iChgBtwEQs);
	}

	return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float   fMinOf2(float val1, float val2)
{	float MinVal;
	MinVal = (val1 < val2) ? val1 : val2;
	return MinVal;
}
float  fMaxOf2(float val1, float val2) 
{	float MaxVal;
	MaxVal = (val1 > val2) ? val1 : val2;
	return MaxVal;
}