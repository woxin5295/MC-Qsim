#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
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
    float     	fLegLgth,          			fUnitSlip,				fDeltT,             		fCutStrss;
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
//------------------------------------------------------------------
struct EQstruct
{	int 	iStillOn,			iEndCntr,			iActFPNum;
	int 	iCmbFPNum,			iMRFLgth,			iTotlRuptT;
	
	float	fSeisPot,			fMaxSlip,			fMaxDTau;
	
	int 	*ivR_WrtStrtPos,	*ivL_ActPtchID,		*ivL_t0ofPtch,		*ivL_StabType;
	float	*fvL_PtchSlpH,		*fvL_PtchSlpV,		*fvL_PtchDTau,		*fvM_MRFvals;

	gsl_vector_float			*fvG_EQslipH,		*fvG_EQslipV;
	gsl_matrix_float			*fmG_STF_H,			*fmG_STF_V;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void           LoadInput(struct MDstruct *MD, struct SGstruct *SG, struct TRstruct *TR, struct VTstruct *VT);
void     DefineMoreParas(struct MDstruct *MD, struct SGstruct *SG, struct TRstruct *TR, struct VTstruct *VT);
void      Build_K_Matrix(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct  Kstruct *K);
float            TheSign(float TestVal);
float GetUpdatedFriction(float B4Fric, float RefFric, float CurFric, float DynFric, int FricLaw, int StabType, float CurD_c, float AccumSlip, float PrevSlip, float HealingFact);
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int main(int argc, char **argv)
{   if (argc != 2 )
    {   fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -np 4 ./MCQsim RunParaFile.txt");
    }
    //------------------------------------------------------------------
    struct 	MDstruct MD; //initializing the structures
    struct 	SGstruct SG;
    struct 	TRstruct TR;
    struct 	VTstruct VT;
    struct 	Kstruct  K;
	struct 	EQstruct EQ;
    
    int    	iPlot2Screen   = 1; //for testing and to ensure that code runs ok
	int		iMinPtch4Cat   = 5; //to not include those small events from saved catalog
    int     iUseVpAndVs    = 1;
	float   fMinSlipForCat = 0.01; //at least 1mm slip on patch to be put into EQcatalog
    //------------------------------------------------------------------
    MPI_Init(&argc, &argv); //initializing MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &MD.iRANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MD.iSIZE);
    MPI_Status  STATUS;
    MPI_Offset  OFFSETall; //this offset is for writing to file
    MPI_File    fp_MPIOUT;
    //------------------------------------------------------------------
    int     i,        	j,            	k,				iGlobPos;
    int    	iTemp0,     iTemp1,			iTemp2,			STFcnt;				
     
    float   fUsedMemory = 0.0; //is just a proxy, idea is to get a feeling for how much memory this run actually uses (per rank)         
    float  	fTemp0,   	fTemp1,      	fTemp2,			fTemp3,  		fTemp4;
	float   fTemp5,		fTemp6,		 	fTemp7,			fTemp8,			fTemp9;
	float   fTemp10,	fTemp11,		fTemp12;
	
	gsl_vector_float 	*fvFG_Temp0,	*fvFG_Temp1,	*fvFG_Temp2;
	gsl_vector_float 	*fvFL_Temp0,	*fvFL_Temp1,	*fvFL_Temp2;
	gsl_vector_float 	*fvBG_Temp0,	*fvBG_Temp1,	*fvBG_Temp2;
	gsl_vector_float 	*fvBG_Temp3,	*fvBG_Temp4,	*fvBG_Temp5;
	gsl_vector_float 	*fvBL_Temp0,	*fvBL_Temp1,	*fvBL_Temp2;

    char 	ctempVals[512],     cFile1_In[512],   cFile2_In[512],   	cFile1_Out[512],    cFile2_Out[512],		cAppend[512];
    FILE 	*fpIn,				*fpPre;
   	//------------------------------------------------------------------
    clock_t timer;
	timer = clock();
    srand(time(0)); 
    //------------------------------------------------------------------
    MD.ivF_START     = (int *) calloc(MD.iSIZE, sizeof(int)); //these starts and offsets relate to how the code accesses "local" and "global" vectors => use the start and offset to 
	MD.ivF_START2    = (int *) calloc(MD.iSIZE, sizeof(int)); //locate where each rank will put/get data from when accessing a "global" vector
    MD.ivB_START     = (int *) calloc(MD.iSIZE, sizeof(int));
    MD.ivF_OFFSET    = (int *) calloc(MD.iSIZE, sizeof(int));
	MD.ivF_OFFSET2   = (int *) calloc(MD.iSIZE, sizeof(int));
    MD.ivB_OFFSET    = (int *) calloc(MD.iSIZE, sizeof(int));
    //------------------------------------------------------------------
    strcpy(cFile1_In,  argv[1]);// opening and reading the "run file" that contains specifics about this model run
    if ((fpIn = fopen(cFile1_In,"r")) == NULL)      {   fprintf(stdout,"Error -cant open %s file. in Main function \n",cFile1_In);      exit(10);     }

    if (fgets(ctempVals, 512, fpIn) != NULL)        {                       	                            }                            
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %s", MD.cInputName);          }             
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iRunNum);            }
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iFUsdGrd);           }                
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iUsePSeis);          } 
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iUseProp);           }           
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fMinMag4Prop);       }
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fISeisStep);         }           
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fViscRelTime);       }
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fAftrSlipTime);      }            
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fRecLgth);           }
	if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD.fHealFact);          }
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iSeedStart);         }
    if (MD.iSeedStart <= 0)				            {   MD.iSeedStart = rand();                             }

	MD.fHealFact = (MD.fHealFact < 0.0)	? 0.0 : MD.fHealFact;
	MD.fHealFact = (MD.fHealFact > 1.0)	? 1.0 : MD.fHealFact;
    //------------------------------------------------------------------
    if ((iPlot2Screen == 1) && (MD.iRANK == 0)) //making sure that the data were imported from file correctly
    {	fprintf(stdout,"Number of RANKS: %d\n",MD.iSIZE);					
		fprintf(stdout,"System info: Byte Size for FLOAT: %lu     INT: %lu    \n\n", sizeof(float), sizeof(int));
		fprintf(stdout,"FileName:           %s\n",MD.cInputName);				    
		fprintf(stdout,"RunNumber:          %d\n",MD.iRunNum);	
    	fprintf(stdout,"UsedGrid:           %d\n",MD.iFUsdGrd);					
		fprintf(stdout,"UsePostSeis:        %d\n",MD.iUsePSeis);
        fprintf(stdout,"UseRuptProp:        %d\n",MD.iUseProp);	                
		fprintf(stdout,"MinMag2UseRuptProp: %f\n",MD.fMinMag4Prop);
    	fprintf(stdout,"IntSeisTStep:       %3.1f\n",MD.fISeisStep);
    	fprintf(stdout,"ViscDeepRlx:        %3.1f\n",MD.fViscRelTime);			
		fprintf(stdout,"ViscAftSlip:        %3.1f\n",MD.fAftrSlipTime);
    	fprintf(stdout,"RecLength:          %5.1f\n",MD.fRecLgth);				
		fprintf(stdout,"CoSeisHealFraction: %5.1f\n",MD.fHealFact);	
		fprintf(stdout,"SeedLocation:       %5d\n\n",MD.iSeedStart);
    }
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    MD.fCutStrss    = 0.05; //stress on patch has to exceed current strength by this amount (in MPa) in order to slip/release stress in the next iteration step
    MD.iGlobTTmax   = 0; //max. global travel time => is used when "use rupture propagation is used
    MD.iEQcntr      = 0; //counting the number of events that were written to file
	MD.iMaxIterat   = 100; //this is for case where I want to use boundary faults for loading (using them to get stressing rate on faults)
	MD.iMaxMRFlgth  = 10000; //some "randomly" high number to ensure that the entire MRF will fit into this vector..
   	MD.fg           = 9.81; //graviational acceleration
	MD.fISeisTStp   = MD.fISeisStep/365.25; //fraction of a year of the interseismic time step
	MD.fTStpInSec   = MD.fISeisTStp*31536000.0; //that fraction of a year in seconds
    MD.iUsePropTest = 0; //this is for speed up, before actually using rupture propagation (which is slow) go through the event without => only if larger than certain magnitude is propagation actually used
    EQ.iEndCntr     = 0; //this counter is used to check if an event is really over (at the end of EQ while loop)
	EQ.fMaxDTau     = 0.0; //max change in shear stress during the event
	EQ.fMaxSlip     = 0.0; //max in-plane slip of the event
    //------------------------------------------------------------------
	//------------------------------------------------------------------
	gsl_rng *fRandN; // this is pretty much straight from the GSL reference, the default RNG has good performance, so no need to change
    const gsl_rng_type *RandT;
	gsl_rng_env_setup(); 
	RandT   = gsl_rng_default;
    fRandN  = gsl_rng_alloc(RandT);
	unsigned long RSeed =(unsigned long)(MD.iRANK + MD.iSeedStart);   //so, every rank has its own random number sequence, and that sequences differs from the sequences of the other ranks
    gsl_rng_set(fRandN, RSeed);
    //------------------------------------------------------------------
    strcpy(cFile1_In,  MD.cInputName);                  strcat(cFile1_In,"_FLT.txt");
    strcpy(cFile2_In,  MD.cInputName);                  strcat(cFile2_In,"_BND.txt");
    strcpy(cFile1_Out, MD.cInputName);                  strcat(cFile1_Out,"_");    	    sprintf(cAppend, "%d",MD.iRunNum); strcat(cFile1_Out,cAppend);   strcat(cFile1_Out,"_");       sprintf(cAppend, "%d",MD.iFUsdGrd);  strcat(cFile1_Out,cAppend);     strcat(cFile1_Out,"_Catalog.dat");
    strcpy(cFile2_Out, MD.cInputName);                  strcat(cFile2_Out,"_");    	    sprintf(cAppend, "%d",MD.iRunNum); strcat(cFile2_Out,cAppend);   strcat(cFile2_Out,"_");       sprintf(cAppend, "%d",MD.iFUsdGrd);  strcat(cFile2_Out,cAppend);     strcat(cFile2_Out,"_PreRunData.dat");
	//------------------------------------------------------------------
    iTemp0 = 0;
    if ((fpIn = fopen(cFile1_In,"r")) == NULL)      {   fprintf(stdout,"Error -cant open %s file. in Main function \n",cFile1_In);      exit(10);       }
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD.iFSegNum);                                                       }//opening the txt files to read the number of fault segments (and other geometry/hierarchy related info)
    if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d", &iTemp0);                                                            }//number of grids that were defined in MATLAB GUI
    //---------------------------
    SG.fvSlipRate = (float *) calloc(MD.iFSegNum, sizeof(float)); //reading segment-related information => here the long-term slip i.e., stressing rate and the corresponding rake
    SG.fvSlipRake = (float *) calloc(MD.iFSegNum, sizeof(float)); 
    SG.fvStrsRate = (float *) calloc(MD.iFSegNum, sizeof(float)); 
    //---------------------------
    MD.iFUsdGrd   = (MD.iFUsdGrd  < iTemp0) ? MD.iFUsdGrd : iTemp0; //if i picked a grid-to-use that is larger than number of grids that were generated, then use "last grid that is generated"
    //---------------------------
    for (i = 0; i < iTemp0; i++)
    {   if ((i+1) == MD.iFUsdGrd)
        {   if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %d %d",&MD.iFPNum,  &MD.iFVNum);  }
            if (fgets(ctempVals, 512, fpIn) != NULL)        {   sscanf(ctempVals,"%*s %e %e",&MD.fFltLegs,&fTemp1);     }
        }
        else
        {   if (fgets(ctempVals, 512, fpIn) != NULL)        {                   }
            if (fgets(ctempVals, 512, fpIn) != NULL)        {                   }
    }   }
    MD.fFltLegs *= 1.0E+3; // now it is in meters 
    //---------------------------
    for (i = 0; i < MD.iFSegNum; i++)
    {   for (j = 0; j < 10; j++)                         
        {  if (fgets(ctempVals, 512, fpIn) != NULL)         {                   }                                
        }
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %e", &SG.fvSlipRate[i]);        	   SG.fvSlipRate[i] /= 1.0E+3;         }// now in m/yr
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %e", &SG.fvSlipRake[i]);         	   SG.fvSlipRake[i] *= M_PI/180.0;     }// now in radiants
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %e", &SG.fvStrsRate[i]);         	   SG.fvStrsRate[i] *= 1.0;            }// now in MPa/yr
        for (j = 0; j < 4; j++)                          
        {   if (fgets(ctempVals, 512, fpIn) != NULL)        {                   }   
        }
        for (j = 0; j < iTemp0; j++)                     
        {   if (fgets(ctempVals, 512, fpIn) != NULL)        {                   }    
            if (fgets(ctempVals, 512, fpIn) != NULL)        {             	    }
    }   }
    fclose(fpIn);
    //------------------------------------------------------------------
    if ((fpIn = fopen(cFile2_In,"r")) == NULL) //if I dont have a boundary fault file to load => then go here
    {   MD.iBSegNum = 0;                                    iTemp1      = 0;
        MD.iBPNum   = 0;                                    MD.iBVNum   = 0;										MD.fBndLegs  = MD.fFltLegs*10.0; //the last thing here is only to have a non-zero value for the BndLeg length that is also larger than the FltLeg length
    }
    else
    {   iTemp1 = 0;
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %d",  &MD.iBSegNum);                  }
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %d",  &iTemp1);                       }//this is grid number (always 1)
        if (iTemp1 < 1)                             	    {   fprintf(stdout,"Error -boundary fault file was opened but BC faults are not meshed/gridded. \n");      exit(10);     }
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %d  %d", &MD.iBPNum,  &MD.iBVNum);    }
        if (fgets(ctempVals, 512, fpIn) != NULL)            {   sscanf(ctempVals,"%*s %e  %e", &MD.fBndLegs, &fTemp1);      }        
        fclose(fpIn);
    }
    MD.fBndLegs *= 1.0E+3; // now it is in meters 
    //------------------------------------------------------------------
    MD.fLegLgth  = (MD.fFltLegs <= MD.fBndLegs) ? MD.fFltLegs : MD.fBndLegs;
 	MD.fUnitSlip = 1.0E-4*MD.fLegLgth; //for a 1000m fault patch (leg length) the 1e-4 means a test slip of 0.1m
	
    //------------------------------------------------------------------
    MD.iF_BASEelem = (int)(MD.iFPNum/MD.iSIZE); 
    MD.iF_ADDelem  = (int)(MD.iFPNum%MD.iSIZE);
    MD.iB_BASEelem = (int)(MD.iBPNum/MD.iSIZE); 
    MD.iB_ADDelem  = (int)(MD.iBPNum%MD.iSIZE);
    //---------------------------
    for (i = 0; i < MD.iSIZE;     i++)      {   MD.ivF_OFFSET[i]     = MD.iF_BASEelem;                                   }
    for (i = 0; i < MD.iF_ADDelem;i++)      {   MD.ivF_OFFSET[i]    += 1;                                                }
    for (i = 1; i < MD.iSIZE;     i++)      {   MD.ivF_START[i]      = MD.ivF_START[i-1] + MD.ivF_OFFSET[i-1];           }
    
    for (i = 0; i < MD.iSIZE;     i++)      {   MD.ivB_OFFSET[i]     = MD.iB_BASEelem;                                   }
    for (i = 0; i < MD.iB_ADDelem;i++)      {   MD.ivB_OFFSET[i]    += 1;                                                }
    for (i = 1; i < MD.iSIZE;     i++)      {   MD.ivB_START[i]      = MD.ivB_START[i-1] + MD.ivB_OFFSET[i-1];           }
	
	for (i = 0; i < MD.iSIZE;     i++)	    {   MD.ivF_OFFSET2[i]    = MD.ivF_OFFSET[i];
												MD.ivF_START2[i]     = MD.ivF_START[i];									 }
	//------------------------------------------------------------------
    if ((iPlot2Screen == 1) && (MD.iRANK == 0))
    {	fprintf(stdout,"FaultPatchNumber %d     FaultVertexNumber %d\n", MD.iFPNum, MD.iFVNum);
    	fprintf(stdout,"BoundPatchNumber %d     BoundVertexNumber %d\n", MD.iBPNum, MD.iBVNum);
    	fprintf(stdout,"MeanLegLenth     %f\n",MD.fLegLgth);   
    }
    //------------------------------------------------------------------
    fvFG_Temp0    = gsl_vector_float_calloc(MD.iFPNum);											fvFG_Temp1    = gsl_vector_float_calloc(MD.iFPNum);											fvFG_Temp2    = gsl_vector_float_calloc(MD.iFPNum);		
    fvBG_Temp0    = gsl_vector_float_calloc(MD.iBPNum);											fvBG_Temp1    = gsl_vector_float_calloc(MD.iBPNum);											fvBG_Temp2    = gsl_vector_float_calloc(MD.iBPNum);		
	fvBG_Temp3    = gsl_vector_float_calloc(MD.iBPNum);											fvBG_Temp4    = gsl_vector_float_calloc(MD.iBPNum);											fvBG_Temp5    = gsl_vector_float_calloc(MD.iBPNum);
    fvFL_Temp0    = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);							fvFL_Temp1    = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);							fvFL_Temp2    = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);
    fvBL_Temp0    = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);							fvBL_Temp1    = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);							fvBL_Temp2    = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);
	//------------------------------------------------------------------
	SG.ivFricLawUSED        = ( int *)  calloc(MD.iFSegNum,   sizeof( int));
    SG.fvRefStaFr           = (float *) calloc(MD.iFSegNum,   sizeof(float));					SG.fvRefStaFr_vari      = (float *) calloc(MD.iFSegNum,   sizeof(float));
    SG.fvRefDynFr           = (float *) calloc(MD.iFSegNum,   sizeof(float));					SG.fvRefDynFr_vari      = (float *) calloc(MD.iFSegNum,   sizeof(float));
    SG.fvRefCritDist        = (float *) calloc(MD.iFSegNum,   sizeof(float));					SG.fvRefCritDist_vari   = (float *) calloc(MD.iFSegNum,   sizeof(float));
    //------------------------------------------------------------------  
	TR.ivFG_Activated       = ( int *)  calloc(MD.iFPNum,     sizeof( int));					TR.ivFG_Ptch_t0         = ( int *)  calloc(MD.iFPNum,     sizeof( int));

    TR.ivFL_FricLaw         = ( int *)  calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof( int));		TR.ivFL_StabT           = ( int *)  calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof( int));	
	TR.ivFL_SelfLoc         = ( int *)  calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof( int));		TR.ivBL_SelfLoc         = ( int *)  calloc(MD.ivB_OFFSET[MD.iRANK],   sizeof( int));

	TR.fvFL_Area           	= (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_RefNrmStrs      = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_RefStaFric     	= (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_RefDynFric      = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_RefDcVal        = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_RefStaFric_vari = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_RefDynFric_vari = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_RefDcVal_vari   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	
	TR.fvFL_StaFric         = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_DynFric         = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_CurFric         = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_CurDcVal        = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_B4_Fric         = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_TempRefFric     = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));	
	TR.fvFL_AccumSlp        = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	//----------------------------------	
	TR.fvFL_StrsRateStk     = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);					TR.fvFL_StrsRateDip     = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);
	TR.fvFL_B4_StrsH        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);					TR.fvFL_B4_StrsV        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);                 TR.fvFL_B4_StrsN        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);
	TR.fvFL_CurStrsH        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);					TR.fvFL_CurStrsV        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);					TR.fvFL_CurStrsN        = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);

	TR.fvFL_PSeis_T0_F     	= (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_PSeis_T0_N      = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvFL_PSeis_Step      = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		
	//----------------------------------
	TR.fvBL_CurStrsH        = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);					TR.fvBL_CurStrsV        = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);					TR.fvBL_CurStrsN        = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);

	TR.fvBL_PSeis_T0_F      = (float *) calloc(MD.ivB_OFFSET[MD.iRANK],   sizeof(float));		TR.fvBL_PSeis_T0_N      = (float *) calloc(MD.ivB_OFFSET[MD.iRANK],   sizeof(float));
	TR.fvBL_PSeis_Step      = (float *) calloc(MD.ivB_OFFSET[MD.iRANK],   sizeof(float));		TR.fvBL_CurFric         = (float *) calloc(MD.ivB_OFFSET[MD.iRANK],   sizeof(float));		
	//--------------------------------------
	TR.ivFL_SegID_temp      = ( int *)  calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof( int));	
	TR.fvFL_SlipRate_temp   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));		TR.fvFL_SlipRake_temp   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK],   sizeof(float));
	TR.ivFG_FltID_temp      = ( int *)  calloc(MD.iFPNum,   sizeof( int));						TR.ivFG_Flagged_temp    = ( int *)  calloc(MD.iFPNum,   sizeof( int));
	TR.ivFG_V1_temp         = ( int *)  calloc(MD.iFPNum,   sizeof( int));						TR.ivFG_V2_temp         = ( int *)  calloc(MD.iFPNum,   sizeof( int));						TR.ivFG_V3_temp         = ( int *)  calloc(MD.iFPNum,   sizeof( int));
	TR.ivBG_V1_temp         = ( int *)  calloc(MD.iBPNum,   sizeof( int));						TR.ivBG_V2_temp         = ( int *)  calloc(MD.iBPNum,   sizeof( int));						TR.ivBG_V3_temp         = ( int *)  calloc(MD.iBPNum,   sizeof( int));

	TR.fvFG_CentE_temp  	= (float *) calloc(MD.iFPNum,   sizeof(float));						TR.fvFG_CentN_temp      = (float *) calloc(MD.iFPNum,   sizeof(float));						TR.fvFG_CentZ_temp      = (float *) calloc(MD.iFPNum,   sizeof(float));	
	TR.fvBG_CentE_temp   	= (float *) calloc(MD.iBPNum,   sizeof(float));						TR.fvBG_CentN_temp      = (float *) calloc(MD.iBPNum,   sizeof(float));						TR.fvBG_CentZ_temp      = (float *) calloc(MD.iBPNum,   sizeof(float));
	//--------------------------------------
	TR.imFGL_TTP            = gsl_matrix_int_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);  		TR.imFGL_TTS            = gsl_matrix_int_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);  
	TR.imFGL_NextP          = gsl_matrix_int_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]); 		TR.imFGL_NextS          = gsl_matrix_int_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]); 
    //these two tell me the index for source-receiver pair => how much of STF at source has been "seen" at receiver already... takes more memory but is worth it b/c it allows to make the STF shorter and ensures that it will not "blow up" when events are too long...
	TR.fmFGL_SrcRcvH        = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);		TR.fmFGL_SrcRcvV        = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);		TR.fmFGL_SrcRcvN        = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);
	//------------------------------------------------------------------
    VT.fvFG_PosE_temp       = (float *) calloc(MD.iFVNum,   sizeof(float));						VT.fvFG_PosN_temp       = (float *) calloc(MD.iFVNum,   sizeof(float));						VT.fvFG_PosZ_temp       = (float *) calloc(MD.iFVNum,   sizeof(float));	
	VT.fvBG_PosE_temp       = (float *) calloc(MD.iBVNum,   sizeof(float));						VT.fvBG_PosN_temp       = (float *) calloc(MD.iBVNum,   sizeof(float));						VT.fvBG_PosZ_temp       = (float *) calloc(MD.iBVNum,   sizeof(float));
	//-------------------------------------------------------------------------------------
	//first dimension is number of rows (slow direction) => is the SOURCES; second dimension is number of columns (fast direction) => is the RECEIVERS
	K.FF_SS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]); 					K.FF_SD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);						K.FF_SO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);
    K.FF_DS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);						K.FF_DD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);						K.FF_DO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);
    K.FF_OS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);						K.FF_OD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);						K.FF_OO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivF_OFFSET[MD.iRANK]);
   	
   	K.FB_SS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]); 					K.FB_SD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);						K.FB_SO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);
    K.FB_DS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);						K.FB_DD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);						K.FB_DO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);
    K.FB_OS  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);						K.FB_OD  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);						K.FB_OO  = gsl_matrix_float_calloc(MD.iFPNum, MD.ivB_OFFSET[MD.iRANK]);
   
 	K.BF_SS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_SD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_SO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);
    K.BF_DS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_DD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_DO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);
    K.BF_OS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_OD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);						K.BF_OO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivF_OFFSET[MD.iRANK]);
   	
 	K.BB_SS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_SD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_SO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);
    K.BB_DS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_DD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_DO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);
    K.BB_OS  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_OD  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);						K.BB_OO  = gsl_matrix_float_calloc(MD.iBPNum, MD.ivB_OFFSET[MD.iRANK]);
 	//-------------------------------------------------------------------------------------
	EQ.ivR_WrtStrtPos = (  int *) calloc(MD.iSIZE,                sizeof( int));				EQ.ivL_ActPtchID  = (  int *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof( int));
	EQ.ivL_t0ofPtch   = (  int *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof( int));				EQ.ivL_StabType   = (  int *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof( int));
	EQ.fvL_PtchSlpH   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof(float));  				EQ.fvL_PtchSlpV   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof(float));  
  	EQ.fvL_PtchDTau   = (float *) calloc(MD.ivF_OFFSET[MD.iRANK], sizeof(float));  				EQ.fvM_MRFvals    = (float *) calloc(MD.iMaxMRFlgth,          sizeof(float));
	EQ.fvG_EQslipH    = gsl_vector_float_calloc(MD.iFPNum);										EQ.fvG_EQslipV    = gsl_vector_float_calloc(MD.iFPNum);
	
	fUsedMemory += (float)(sizeof(float)*(MD.iFPNum+MD.iBPNum)*(MD.iFPNum+MD.iBPNum)*9.0 * 1.0e-9); //this is this looking at the largest parts i.e., the stiffness matrix and later adding the STFpart and the few temp vectors...
	fUsedMemory /= (float)MD.iSIZE; //know how much each RANK gets (approximately)
	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------
	
	LoadInput(&MD, &SG, &TR, &VT);

	//----------------------------------
    if (iUseVpAndVs  == 1)          {   MD.fVpVsRatio = MD.fVp/MD.fVs;                  }
    else                            {   MD.fVpVsRatio = 1.0;                			}              
	if (MD.iUseProp == 1)           {   MD.fDeltT     = 1.0*MD.fLegLgth/MD.fVp;         }
    else                            {   MD.fDeltT     = FLT_MAX;      	     			}
	//----------------------------------

	DefineMoreParas(&MD, &SG, &TR, &VT);

	//----------------------------------
	MPI_Allreduce(MPI_IN_PLACE, &MD.iGlobTTmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    MD.iMaxSTFlgth = MD.iGlobTTmax +1;//use the +1 here to ensure I didn't mess this up and to make sure that the overwriting of the STF does not happen too fast (b/c I'm looping over the STF to save memory)
	//----------------------------------
	for (i = 0; i < MD.iSIZE;     i++)		{	MD.ivF_OFFSET2[i]   *= MD.iMaxSTFlgth;			MD.ivF_START2[i]   *= MD.iMaxSTFlgth;			}
	//----------------------------------
	EQ.fmG_STF_H = gsl_matrix_float_calloc(MD.iFPNum, MD.iMaxSTFlgth);							EQ.fmG_STF_V = gsl_matrix_float_calloc(MD.iFPNum, MD.iMaxSTFlgth);
	fUsedMemory += (float)(sizeof(float)*MD.iFPNum*MD.iMaxSTFlgth*2.0*1.0e-9); //this is the size of the STF matrix
	fUsedMemory += (float)(sizeof(float)*(MD.iFPNum+MD.iBPNum)*4.0*1.0e-9); //rough approx. for the size of the temporary vectors that are used
	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------    
	if ((iPlot2Screen == 1) && (MD.iRANK == 0))
    {	fprintf(stdout,"Medium density:  %5.2f\n",MD.fMedDense);
		fprintf(stdout,"Shear modulus:   %e\n",MD.fShearMod);
 		fprintf(stdout,"Poisson ratio:   %5.2f\n",MD.fPoisson);
 		fprintf(stdout,"Lambda:          %e\n",MD.fLambda);
		fprintf(stdout,"P-wave velocity: %5.2f\n",MD.fVp);
		fprintf(stdout,"S-wave velocity: %5.2f\n",MD.fVs);
		fprintf(stdout,"Added normal stress (MPa):   %5.2f\n",MD.fAddNrmStrs*1.0e-6);
		fprintf(stdout,"Change friction between EQs: %d\n",MD.iChgBtwEQs);
		fprintf(stdout,"\nValues per Segment:\n");
		for (i = 0; i < MD.iFSegNum; i++)
		{	fprintf(stdout,"Segment# %d   SegmentSlipRate (m/yr):       %5.5f\n",i,   SG.fvSlipRate[i]);
			fprintf(stdout,"Segment# %d   SegmentSlipRake:              %5.2f\n",i,   SG.fvSlipRake[i]);
			fprintf(stdout,"Segment# %d   SegmentStressRate (MPa):      %5.2f\n",i, SG.fvStrsRate[i]);
			fprintf(stdout,"Segment# %d   Friction law used:            %d\n", i, SG.ivFricLawUSED[i]);
			fprintf(stdout,"Segment# %d   Reference static friction:    %4.4f\n", i, SG.fvRefStaFr[i]);
			fprintf(stdout,"Segment# %d   Static friction variation:    %4.4f\n", i, SG.fvRefStaFr_vari[i]);
			fprintf(stdout,"Segment# %d   Reference dynamic friction:   %4.4f\n", i, SG.fvRefDynFr[i]);
			fprintf(stdout,"Segment# %d   Dynamic friction variation:   %4.4f\n", i, SG.fvRefDynFr_vari[i]);
			fprintf(stdout,"Segment# %d   Critical slip distance (m):   %4.6f\n", i, SG.fvRefCritDist[i]*MD.fLegLgth);
			fprintf(stdout,"Segment# %d   Critical slip variation:      %4.4f\n\n", i, SG.fvRefCritDist_vari[i]);
 		}	
		fprintf(stdout,"Other stuff after loading input and defining more parameters\n");
		fprintf(stdout,"MaxSTFlgth: %d      and  %d\n",MD.iMaxSTFlgth, MD.iGlobTTmax);
		fprintf(stdout,"Delta time: %4.4f   and Unit slip: %4.4f \n\n",MD.fDeltT, MD.fUnitSlip);
		fprintf(stdout,"APPROXIMATE MEMORY USAGE PER RANK: %2.5f Gb\n",fUsedMemory);
	}
 	//-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------	
   	if ((iPlot2Screen == 1) && (MD.iRANK == 0))	{	fprintf(stdout,"Building K-matrix....\n");	}
    
	Build_K_Matrix(&MD, &TR, &VT, &K); //later on, if helpful, this could also be stored to file and the loaded here...
	MPI_Barrier( MPI_COMM_WORLD );

	if ((iPlot2Screen == 1) && (MD.iRANK == 0))	{	fprintf(stdout,"K-matrix done\n");	}
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------	
	for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)     {			if (fabs(TR.fvFL_SlipRate_temp[i]) > 0.0)    {   iTemp0 = 1;    }      			}
	MPI_Allreduce(MPI_IN_PLACE, &iTemp0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, TR.ivFG_Flagged_temp, MD.iFPNum, MPI_INT, MPI_MAX, MPI_COMM_WORLD); //for output to write out what fault patches have been flagged
	//-----------------------------------------------
 	if ((iPlot2Screen == 1) && (MD.iRANK == 0))            
    {   fTemp0 = MD.fCutStrss / -gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[0], 0); 
		fprintf(stdout,"Threshold stress amount (per increment) to leave EQ iteration (MPa): %4.4f. This is approx. %4.4f m of slip (looking only at one example patch).\n",MD.fCutStrss, fTemp0);             
        fprintf(stdout,"Use slip boundary condition? %d\n",iTemp0);
    }	
	//------------------------------------------------------------------
    if (iTemp0 == 1) //have a slip boundary condition on at least one patch
    {   fTemp0 = 0.0;
		//--------------------------------------------------
        for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
        {   iGlobPos = i + MD.ivF_START[MD.iRANK];
            fTemp0  += fabs(TR.fvFL_SlipRate_temp[i]); 
			fTemp1   =  -1.0*cosf(TR.fvFL_SlipRake_temp[i])*TR.fvFL_SlipRate_temp[i];
            fTemp2   =  -1.0*sinf(TR.fvFL_SlipRake_temp[i])*TR.fvFL_SlipRate_temp[i];
				
			gsl_vector_float_set(fvFG_Temp0, iGlobPos,fTemp1);  //the strike slip component            
			gsl_vector_float_set(fvFG_Temp1, iGlobPos,fTemp2);  //the dip slip component    
		}
		MPI_Allreduce(MPI_IN_PLACE, &fTemp0, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); 
		fTemp0 /= (float)MD.iFPNum;
		if ((iPlot2Screen == 1) &&(MD.iRANK == 0))		{			fprintf(stdout,"Prescribed average slip-rate on faults (m/yr): %5.5f \n",fTemp0);			}
		//--------------------------------------------------
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);			
		//--------------------------------------------------
		if (MD.iBPNum > 0) //I'm using boundary box faults => use the slip boundary condition on the EQ faults to load the boundary faults; at the same time, also "load/release" stress on the EQfaults
        {   //--------------------------------------------------
			gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SS, fvFG_Temp0, 0.0, fvBL_Temp0); 				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SD, fvFG_Temp0, 0.0, fvBL_Temp1); 				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SO, fvFG_Temp0, 0.0, fvBL_Temp2); 
			gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);									gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
			gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DS, fvFG_Temp1, 0.0, fvBL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DD, fvFG_Temp1, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DO, fvFG_Temp1, 0.0, fvBL_Temp2); 
			gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);									gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
			//--------------------------------------------------
            //actually only need the top part to load the boundaries, BUT also do the regular back-slipping (following) to use that value to scale the boundary fault loading, ensuring the rates are consistent
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SS, fvFG_Temp0, 0.0, fvFL_Temp0); 				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SD, fvFG_Temp0, 0.0, fvFL_Temp1);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);									gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DS, fvFG_Temp1, 0.0, fvFL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DD, fvFG_Temp1, 0.0, fvFL_Temp1);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);									gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			//----------------------------
			fTemp0  = gsl_blas_sasum(TR.fvFL_CurStrsH);
			fTemp1  = gsl_blas_sasum(TR.fvFL_CurStrsV);
			fTemp6  = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1);
			MPI_Allreduce(MPI_IN_PLACE, &fTemp6, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //fTemp6 contains the total/combined stressing rate from back-slip 
			fTemp6 /= (float)MD.iFPNum;  //again, this value is later used in a scaling step (when boundary fault loading is used...)
			//--------------------------------------------------
			gsl_vector_float_set_zero(TR.fvFL_CurStrsH);
			gsl_vector_float_set_zero(TR.fvFL_CurStrsV);
			//--------------------------------------------------
			for (k = 0; k < MD.iMaxIterat; k++)		
			{	//-------------------------------------	
				for (i = 0; i < MD.ivB_OFFSET[MD.iRANK]; i++)  //going through the boundary faults to iteratively release the stress that was put onto them by the EQfaults => the corresponding slip to achieve that is stored and defines the boundary fault loading 
            	{	iGlobPos = i + MD.ivB_START[MD.iRANK];
					fTemp0   = -1.0*gsl_vector_float_get(TR.fvBL_CurStrsH,i) / gsl_matrix_float_get(K.BB_SS, TR.ivBL_SelfLoc[i], i);
					fTemp1   = -1.0*gsl_vector_float_get(TR.fvBL_CurStrsV,i) / gsl_matrix_float_get(K.BB_DD, TR.ivBL_SelfLoc[i], i);
					fTemp2   = -1.0*gsl_vector_float_get(TR.fvBL_CurStrsN,i) / gsl_matrix_float_get(K.BB_OO, TR.ivBL_SelfLoc[i], i);
						
					gsl_vector_float_set(fvBG_Temp0, iGlobPos, fTemp0); 
					gsl_vector_float_set(fvBG_Temp1, iGlobPos, fTemp1);
					gsl_vector_float_set(fvBG_Temp2, iGlobPos, fTemp2);
				}	
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp0->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp1->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp2->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
				//-------------------------------------
				gsl_vector_float_add(fvBG_Temp3, fvBG_Temp0); //collecting the total slip in strike/dip/normal => used later for loading of fault patches
				gsl_vector_float_add(fvBG_Temp4, fvBG_Temp1);
				gsl_vector_float_add(fvBG_Temp5, fvBG_Temp2);
				//----------------------
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SS, fvBG_Temp0, 0.0, fvBL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SD, fvBG_Temp0, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SO, fvBG_Temp0, 0.0, fvBL_Temp2); 
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);									gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DS, fvBG_Temp1, 0.0, fvBL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DD, fvBG_Temp1, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DO, fvBG_Temp1, 0.0, fvBL_Temp2); 
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);									gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OS, fvBG_Temp2, 0.0, fvBL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OD, fvBG_Temp2, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OO, fvBG_Temp2, 0.0, fvBL_Temp2); 
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);									gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				//--------------------------------------------------
			}
			//--------------------------------------------------
			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_SS, fvBG_Temp3, 0.0, fvFL_Temp0); 					gsl_blas_sgemv(CblasTrans, 1.0, K.BF_SD, fvBG_Temp3, 0.0, fvFL_Temp1); 
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);										gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_DS, fvBG_Temp4, 0.0, fvFL_Temp0);					gsl_blas_sgemv(CblasTrans, 1.0, K.BF_DD, fvBG_Temp4, 0.0, fvFL_Temp1); 
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);										gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_OS, fvBG_Temp5, 0.0, fvFL_Temp0);					gsl_blas_sgemv(CblasTrans, 1.0, K.BF_OD, fvBG_Temp5, 0.0, fvFL_Temp1);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);										gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			//----------------------------
			fTemp0  = gsl_blas_sasum(TR.fvFL_CurStrsH); //here I now determine how much stressing (combined) was actually put onto the EQfaults by letting the boundary faults do the loading
			fTemp1  = gsl_blas_sasum(TR.fvFL_CurStrsV); //ideally, both should be the same?! -the classic back-slip stressing rate and the boundary box stressing rate, but they aren't => hence the scaling
			fTemp4  = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1);//so, here is the second value for the scaling ratio determined
			MPI_Allreduce(MPI_IN_PLACE, &fTemp4, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //fTemp4 contains the total/combined stressing rate from boundary fault slip
			fTemp4 /= (float)MD.iFPNum;  
			fTemp5  = fTemp6/fTemp4; //this is the actual scaling factor to be applied...

			//--------------------------------------------------
			if ((iPlot2Screen == 1) &&(MD.iRANK == 0))		{			fprintf(stdout,"Stressing Rate (from Boundary): %5.5f    Stressing Rate (from Back-Slip): %5.5f    Scale Factor: %5.5f\n",fTemp4, fTemp6, fTemp5);			}
			//--------------------------------------------------
			gsl_vector_float_scale(TR.fvFL_CurStrsH, fTemp5); //these are just temporary place holders (not currrent stress but stressing rate put in here) at the moment...
			gsl_vector_float_scale(TR.fvFL_CurStrsV, fTemp5); //the values are actually stressing rates, will be assigned within the next few lines...
		}
		else
		{   //this here is the normal/classic back-slip method (if no boundary faults are used but a slip-boundary condition is defined)
            gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SS, fvFG_Temp0, 0.0, fvFL_Temp0);					gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SD, fvFG_Temp0, 0.0, fvFL_Temp1);
            gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);										gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DS, fvFG_Temp1, 0.0, fvFL_Temp0);					gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DD, fvFG_Temp1, 0.0, fvFL_Temp1);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);										gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);
		}
        
        //now, the values from the temporary file are copied/added to Stressing Rate; the "adding" part allows to have mixed boundary conditions -not that this is a good idea... but it is possible to define one fault with stressing rate and another with slip-rate and then get combined stressing rate for each
		//gsl_vector_float_add(TR.fvFL_StrsRateStk, TR.fvFL_CurStrsH);			gsl_vector_float_add(TR.fvFL_StrsRateDip, TR.fvFL_CurStrsV);
        //TEST
        for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)  //going through the boundary faults to iteratively release the stress that was put onto them by the EQfaults => the corresponding slip to achieve that is stored and defines the boundary fault loading 
        {   fTemp0 = sqrtf( gsl_vector_float_get(TR.fvFL_CurStrsH,i)*gsl_vector_float_get(TR.fvFL_CurStrsH,i) + gsl_vector_float_get(TR.fvFL_CurStrsV,i)*gsl_vector_float_get(TR.fvFL_CurStrsV,i) );
            fTemp1 = fTemp0*cosf(TR.fvFL_SlipRake_temp[i]);
            fTemp2 = fTemp0*sinf(TR.fvFL_SlipRake_temp[i]);
            gsl_vector_float_set(TR.fvFL_StrsRateStk,i, fTemp1);
            gsl_vector_float_set(TR.fvFL_StrsRateDip,i, fTemp2);
        }
        //TEST

	}	
	//--------------------------------------------------
	fTemp5 = 0.0;
	for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
	{	fTemp5 += sqrtf(gsl_vector_float_get(TR.fvFL_StrsRateStk, i)*gsl_vector_float_get(TR.fvFL_StrsRateStk, i) + gsl_vector_float_get(TR.fvFL_StrsRateDip, i)*gsl_vector_float_get(TR.fvFL_StrsRateDip, i)); //the amount of applied stressing rate (per interseismic time step)
    }
	//------------------------------------------------------------------------------------------
	MPI_Allreduce(MPI_IN_PLACE, &fTemp5, 1 , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //total applied stressing rate

    fTemp5 /= (float)MD.iFPNum;			
    if ((iPlot2Screen == 1) && (MD.iRANK == 0))	
	{	    fprintf(stdout,"Resulting average stressing-rate on faults (MPa/yr): %5.5f and per time step: %5.5f\n",fTemp5, (fTemp5*MD.fISeisTStp));  
	}
	gsl_vector_float_scale(TR.fvFL_StrsRateStk, MD.fISeisTStp ); //reference stressing rate is now including the time step => only need to add that amount
	gsl_vector_float_scale(TR.fvFL_StrsRateDip, MD.fISeisTStp ); //when applicable (when stepping forward by one interseis. time step
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	if (MD.iRANK == 0)           
    {	if ((fpPre = fopen(cFile2_Out,"wb")) == NULL)      {   		exit(10);       }
    	fwrite(&MD.iFPNum,           sizeof(int),        1, fpPre);  
		fwrite(TR.ivFG_FltID_temp,   sizeof(int),MD.iFPNum, fpPre); //the fault id
		fwrite(TR.ivFG_Flagged_temp, sizeof(int),MD.iFPNum, fpPre); //the "flagged" status for interaction => if removed/overridden etc...
	}
	for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
    {   iGlobPos = i + MD.ivF_START[MD.iRANK];
		gsl_vector_float_set(fvFG_Temp0, iGlobPos, gsl_vector_float_get(TR.fvFL_StrsRateStk, i));  //the strike slip component            
		gsl_vector_float_set(fvFG_Temp1, iGlobPos, gsl_vector_float_get(TR.fvFL_StrsRateDip, i));  //the dip slip component    
	}
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);	
    if (MD.iRANK == 0)           
    {	fwrite(fvFG_Temp0->data, sizeof(float), MD.iFPNum, fpPre); //the stressing rate in strike direction
		fwrite(fvFG_Temp1->data, sizeof(float), MD.iFPNum, fpPre); //the stressing rate in dip direction
	}
	for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
    {   iGlobPos = i + MD.ivF_START[MD.iRANK];
		fTemp0   = TR.fvFL_StaFric[i]                        *-1.0*TR.fvFL_RefNrmStrs[i]; //this is the static strength of the patch
		fTemp1   = (TR.fvFL_StaFric[i] - TR.fvFL_DynFric[i]) *-1.0*TR.fvFL_RefNrmStrs[i]; //this is the stress drop of the patch
		fTemp2   = TR.fvFL_CurDcVal[i];

		gsl_vector_float_set(fvFG_Temp0, iGlobPos, fTemp0);          
		gsl_vector_float_set(fvFG_Temp1, iGlobPos, fTemp1);   
		gsl_vector_float_set(fvFG_Temp2, iGlobPos, fTemp2);
	}	
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);	
	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp2->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    if (MD.iRANK == 0)           
    {	fwrite(fvFG_Temp0->data, sizeof(float), MD.iFPNum, fpPre); //the patch strength
		fwrite(fvFG_Temp1->data, sizeof(float), MD.iFPNum, fpPre); //the patch stress drop when having full drop
		fwrite(fvFG_Temp2->data, sizeof(float), MD.iFPNum, fpPre); //the Dc value
		fclose(fpPre);
	}
	//------------------------------------------------------------------------------------------	
	//------------------------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
	//------------------------------------------------------------------------------------------	
	//------------------------------------------------------------------------------------------
	//open the file here once and then keep open (if not kept open, the continuous opening and closing slows things down and also might cause code to crash
	MPI_File_delete(cFile1_Out, MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, cFile1_Out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIOUT);
    if (MD.iRANK == 0)           
    {   MPI_File_write(fp_MPIOUT, &MD.iEQcntr,   1, MPI_INT,      &STATUS);
        MPI_File_write(fp_MPIOUT, &MD.iFUsdGrd,  1, MPI_INT,      &STATUS);            
        MPI_File_write(fp_MPIOUT, &MD.fShearMod, 1, MPI_FLOAT,    &STATUS);    
        MPI_File_write(fp_MPIOUT, &MD.fDeltT,    1, MPI_FLOAT,    &STATUS);           
    }
    MPI_Barrier( MPI_COMM_WORLD );
    OFFSETall = 2*(sizeof(int) + sizeof(float));
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------		
	free(SG.fvSlipRate);			free(SG.fvSlipRake);			free(SG.fvStrsRate);			free(SG.fvRefStaFr);			free(SG.fvRefDynFr);			free(SG.fvRefCritDist);
	free(SG.fvRefStaFr_vari);		free(SG.fvRefDynFr_vari);		free(SG.fvRefCritDist_vari);
	free(TR.fvFL_SlipRate_temp);	free(TR.fvFL_SlipRake_temp);	free(TR.ivFG_V1_temp);			free(TR.ivFG_V2_temp);			free(TR.ivFG_V3_temp);
	free(TR.ivBG_V1_temp);			free(TR.ivBG_V2_temp);			free(TR.ivBG_V3_temp);			free(TR.fvFG_CentE_temp);		free(TR.fvFG_CentN_temp);		free(TR.fvFG_CentZ_temp);
	free(TR.fvBG_CentE_temp);		free(TR.fvBG_CentN_temp);		free(TR.fvBG_CentZ_temp);		free(TR.ivFL_SegID_temp);	
	free(VT.fvFG_PosE_temp);		free(VT.fvFG_PosN_temp);		free(VT.fvFG_PosZ_temp);		free(VT.fvBG_PosE_temp);		free(VT.fvBG_PosN_temp);		free(VT.fvBG_PosZ_temp);	
	free(TR.ivFG_FltID_temp);	
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//load all fault patches to be fully loaded -then subtracting some fraction to make it not too artificial in the beginning of the catalog
    for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
    {	fTemp0            = gsl_vector_float_get(TR.fvFL_StrsRateStk, i);   
		fTemp1            = gsl_vector_float_get(TR.fvFL_StrsRateDip, i);   
		fTemp2            = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1); //this is combined loading per time step
        fTemp3            = (TR.fvFL_CurFric[i]*-1.0*TR.fvFL_RefNrmStrs[i]);//this is strength of element; the "-1" is here to make normal stress compression positive again (to get positive strength value)
        fTemp4            = fTemp3/fTemp2; //this is number of loading steps I'd need to reach the strength level
      	fTemp5            = (float)(gsl_rng_uniform(fRandN) *0.09 +0.90); //is a random number between 0.90 and 0.99 => so many loading steps to use (fraction of temp2)
		gsl_vector_float_set(TR.fvFL_CurStrsH, i, fTemp0*fTemp4*fTemp5); //stressing-rate * number of stressing/loading steps * random number between 0.9 and 0.99
		gsl_vector_float_set(TR.fvFL_CurStrsV, i, fTemp1*fTemp4*fTemp5);
		gsl_vector_float_set(TR.fvFL_CurStrsN, i, TR.fvFL_RefNrmStrs[i]);
    }   
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
    MD.fTimeYears = 0.0;
    while (MD.fTimeYears <= MD.fRecLgth)
    {   
		MD.fTimeYears      += MD.fISeisTStp;
		//--------------------------------------------------------------
		//first, the regular loading step, is first applied to ALL fault patches, using the respective reference stressing rate
		//this is done to all! (non-boundary) patches => if a patch is aseismic, then this stress will be further redistributed in a following step
     	gsl_vector_float_add(TR.fvFL_CurStrsH, TR.fvFL_StrsRateStk);   
		gsl_vector_float_add(TR.fvFL_CurStrsV, TR.fvFL_StrsRateDip);  	
		//--------------------------------------------------------------
		if (MD.iUsePSeis == 1) //  https://en.wikipedia.org/wiki/Stress_relaxation    http://web.mit.edu/course/3/3.11/www/modules/visco.pdf  => page 9ff; using Maxwell spring-dashpot model    
      	{	iTemp0 = 1;	
			//-------------------------------------------------------------
			for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
			{	//-----------------------------
                //first doing the normal stress component, this will be released regardless if above or below reference value => always want's to go back to reference value...
				iGlobPos = i + MD.ivF_START[MD.iRANK];    
				TR.fvFL_PSeis_Step[i] += 1.0;                   	
                fTemp0   = expf(-1.0*(MD.fTStpInSec*TR.fvFL_PSeis_Step[i])/MD.fAftrSlipTime); //factor/fraction from decay function
             
			    fTemp1   = fabs(fTemp0 *(TR.fvFL_PSeis_T0_N[i]       - TR.fvFL_RefNrmStrs[i])); // this is "allowed"  excess normal stress value (relative to reference normal stress and time in postseismic cycle)
                fTemp2   = gsl_vector_float_get(TR.fvFL_CurStrsN, i) - TR.fvFL_RefNrmStrs[i]; // this is the current excess normal stress (i.e. deviation, positive or negative relative to reference stress) that is present					  
                if (fTemp1 < fabs(fTemp2)) //if the allowed excess (in absolute terms) is smaller than the currently present excess (in absolute terms) => then I have excess to release
				{	fTemp3   = fTemp2 - fTemp1*TheSign(fTemp2); //this is the amount of normal stress that needs to be released via post-seismic relaxation in the current time step (this part being normal stresses..)
                    fTemp4   = -1.0*fTemp3 / gsl_matrix_float_get(K.FF_OO, TR.ivFL_SelfLoc[i], i); // fTemp4 is the corresponding slip amount => opening slip required to release that excess/lack-of normal stress
 
                    fTemp9   = gsl_vector_float_get(TR.fvFL_CurStrsN, i) - fTemp3; //this is the current/next normal stress value on the patch that is permitted by  => need for strength calculation
                    gsl_vector_float_set(fvFG_Temp2, iGlobPos, fTemp4); //set that offset amount into temp global fault vector		
				}
				else
				{	fTemp9   = gsl_vector_float_get(TR.fvFL_CurStrsN, i); // in that case the normal stress is not exceeding (in absolute terms) the allowed excess as defined by post-seismic relax. => use the curent normal stress for current-strength calculation
					gsl_vector_float_set(fvFG_Temp2, iGlobPos, 0.0);    
				}
				//-----------------------------
				if (TR.ivFL_StabT[i] != 1) //if fault element is NOT unstable, then not only the nomal component but also the shear component can be adjusted postseismically	
				{	
					fTemp1   = fabs(fTemp0 *(TR.fvFL_PSeis_T0_F[i]  - TR.fvFL_StaFric[i])); // this is "allowed"  excess friction coeff. value (relative to static friction and time in postseismic cycle)
                	fTemp2   = TR.fvFL_CurFric[i]                   - TR.fvFL_StaFric[i]; // this is the current excess friction coefficient (i.e. deviation, positive or negative relative to reference/static friction) that is present					
                	if (fTemp1 < fabs(fTemp2)) //if the allowed excess (in absolute terms) is smaller than the currently present excess (in absolute terms) => then I have to modify the friction coefficient (either "healing" or "relaxation")
					{	//this part here only occurs/happens/is possible for stable patches, b/c for unstable and cond. stable patches, the CurFric is SET to StaFric after slip occurrence => the difference between both is zero...
				    	TR.fvFL_CurFric[i] =  TR.fvFL_CurFric[i] - fTemp2 + fTemp1*TheSign(fTemp2); //this is the updated/new friction value; combined with the new normal stress, I compute the current/new strength and therefore determine what excess shear stress I have and must release
					}
					//----------------------------
					fTemp7 = gsl_vector_float_get(TR.fvFL_CurStrsH, i); //the currently applied shear stress in strike-direction
                    fTemp8 = gsl_vector_float_get(TR.fvFL_CurStrsV, i); //same for dip-direction
                    fTemp3 = sqrtf(fTemp7*fTemp7 + fTemp8*fTemp8); //this is currently applied shear stress, combined
					fTemp2 = fTemp3 - (TR.fvFL_CurFric[i]*-1.0*fTemp9); //this is currently applied shear stress minus current strength; ==> current excess shear stress above  current/new strength
					
					if (fTemp2 > 0.0)  //if the current stress (fTemp3) exceeds the current strength then ftemp2 is positive and the excess amount is released
					{	fTemp5   = fTemp7/fTemp3 *fTemp2; //this is the amount of in-strike shear that is release
						fTemp6   = fTemp8/fTemp3 *fTemp2; //this is the amount of in-dip shear that is relased
						fTemp5   = -1.0*fTemp5 / gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[i], i); //now converted to the corresponding slip
						fTemp6   = -1.0*fTemp6 / gsl_matrix_float_get(K.FF_DD, TR.ivFL_SelfLoc[i], i); 

						gsl_vector_float_set(fvFG_Temp0, iGlobPos, fTemp5); 
						gsl_vector_float_set(fvFG_Temp1, iGlobPos, fTemp6);
					}
					else
					{	gsl_vector_float_set(fvFG_Temp0, iGlobPos, 0.0); 
						gsl_vector_float_set(fvFG_Temp1, iGlobPos, 0.0); 
				}   }
				else
				{	gsl_vector_float_set(fvFG_Temp0, iGlobPos, 0.0); 
					gsl_vector_float_set(fvFG_Temp1, iGlobPos, 0.0); 
			}   }
			//-------------------------------------------------------------
            for (i = 0; i < MD.ivB_OFFSET[MD.iRANK]; i++) 
            {	iGlobPos = i + MD.ivB_START[MD.iRANK];
				TR.fvBL_PSeis_Step[i] += 1.0;   
				fTemp0   = expf(-1.0*(MD.fTStpInSec *TR.fvBL_PSeis_Step[i])/MD.fViscRelTime);   // this is the expected value (fraction) for strength relative to initial value after last earthquake

				fTemp1   = fabs(fTemp0 *(TR.fvBL_PSeis_T0_N[i]       - 0.0)); // this is "allowed"  excess normal stress value (relative to reference normal stress and time in postseismic cycle)
                fTemp2   = gsl_vector_float_get(TR.fvBL_CurStrsN, i) - 0.0; // this is the current excess normal stress (i.e. deviation, positive or negative relative to reference stress) that is present					  
                if (fTemp1 < fabs(fTemp2)) //if the allowed excess (in absolute terms) is smaller than the currently present excess (in absolute terms) => then I have excess to release
				{	fTemp3   = fTemp2 - fTemp1*TheSign(fTemp2); //this is the amount of normal stress that needs to be released via post-seismic relaxation in the current time step (this part being normal stresses..)
                    fTemp4   = -1.0*fTemp3 / gsl_matrix_float_get(K.BB_OO, TR.ivBL_SelfLoc[i], i); // fTemp4 is the corresponding slip amount => opening slip required to release that excess/lack-of normal stress
 
                    fTemp9   = gsl_vector_float_get(TR.fvBL_CurStrsN, i) - fTemp3; //this is the current/next normal stress value on the patch that is permitted by  => need for strength calculation
                    gsl_vector_float_set(fvBG_Temp2, iGlobPos, fTemp4); //set that offset amount into temp global fault vector		
				}
				else
				{	fTemp9   = gsl_vector_float_get(TR.fvBL_CurStrsN, i); // in that case the normal stress is not exceeding (in absolute terms) the allowed excess as defined by post-seismic relax. => use the curent normal stress for current-strength calculation
					gsl_vector_float_set(fvBG_Temp2, iGlobPos, 0.0);    
				}
                //-----------------------------
				fTemp1   = fabs(fTemp0 *(TR.fvBL_PSeis_T0_F[i]  - 0.0)); // this is "allowed"  excess friction coeff. value (relative to static friction and time in postseismic cycle)
                fTemp2   = TR.fvBL_CurFric[i]                   - 0.0; // this is the current excess friction coefficient (i.e. deviation, positive or negative relative to reference/static friction) that is present					
                if (fTemp1 < fabs(fTemp2)) //if the allowed excess (in absolute terms) is smaller than the currently present excess (in absolute terms) => then I have to modify the friction coefficient (either "healing" or "relaxation")
				{	
				    TR.fvBL_CurFric[i] =  TR.fvBL_CurFric[i] - fTemp2 + fTemp1*TheSign(fTemp2); //this is the updated/new friction value; combined with the new normal stress, I compute the current/new strength and therefore determine what excess shear stress I have and must release
				}
				//----------------------------
				fTemp7 = gsl_vector_float_get(TR.fvBL_CurStrsH, i); //the currently applied shear stress in strike-direction
                fTemp8 = gsl_vector_float_get(TR.fvBL_CurStrsV, i); //same for dip-direction
                fTemp3 = sqrtf(fTemp7*fTemp7 + fTemp8*fTemp8); //this is currently applied shear stress, combined
				fTemp2 = fTemp3 - (TR.fvBL_CurFric[i]*-1.0*fTemp9); //this is currently applied shear stress minus current strength; ==> current excess shear stress above  current/new strength
					
				if (fTemp2 > 0.0)  //if the current stress (fTemp3) exceeds the current strength then ftemp2 is positive and the excess amount is released
				{	fTemp5   = fTemp7/fTemp3 *fTemp2; //this is the amount of in-strike shear that is release
					fTemp6   = fTemp8/fTemp3 *fTemp2; //this is the amount of in-dip shear that is relased
					fTemp5   = -1.0*fTemp5 / gsl_matrix_float_get(K.BB_SS, TR.ivBL_SelfLoc[i], i); //now converted to the corresponding slip
					fTemp6   = -1.0*fTemp6 / gsl_matrix_float_get(K.BB_DD, TR.ivBL_SelfLoc[i], i); 

					gsl_vector_float_set(fvBG_Temp0, iGlobPos, fTemp5); 
					gsl_vector_float_set(fvBG_Temp1, iGlobPos, fTemp6);
				}
				else
				{	gsl_vector_float_set(fvBG_Temp0, iGlobPos, 0.0); 
					gsl_vector_float_set(fvBG_Temp1, iGlobPos, 0.0); 	
		}	}	} 	
		else
		{	//---------------------------
			iTemp0 = 0;
			for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
        	{   
				iGlobPos = i + MD.ivF_START[MD.iRANK];
				if (TR.ivFL_StabT[i] != 1) 
                {	fTemp0   = gsl_vector_float_get(TR.fvFL_CurStrsH, i);
					fTemp1   = gsl_vector_float_get(TR.fvFL_CurStrsV, i);
					fTemp2   = sqrtf(fTemp0*fTemp0 + fTemp1*fTemp1);
					fTemp3   = fTemp2 - TR.fvFL_CurFric[i]*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //this is the excess stress (is excess if value > 0) I currently have with respect to curr strength

            		if (fTemp3 >  MD.fCutStrss) //if I have excess shear stress above current strength and I am not looking at an unstable patch
            		{   iTemp0   = 1;
						fTemp4   = -1.0*(fTemp3/fTemp2*fTemp0) / gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[i], i); // slip amount to release excess horizontal shear stress
                    	fTemp5   = -1.0*(fTemp3/fTemp2*fTemp1) / gsl_matrix_float_get(K.FF_DD, TR.ivFL_SelfLoc[i], i); // slip amount to release excess vertical shear stress 
						gsl_vector_float_set(fvFG_Temp0, iGlobPos,fTemp4);  //the strike slip component            
						gsl_vector_float_set(fvFG_Temp1, iGlobPos,fTemp5);  //the dip slip component  
						gsl_vector_float_set(fvFG_Temp2, iGlobPos,   0.0);
                	}
                	else
                	{   gsl_vector_float_set(fvFG_Temp0, iGlobPos,   0.0);  //the strike slip component            
						gsl_vector_float_set(fvFG_Temp1, iGlobPos,   0.0);  //the dip slip component  
						gsl_vector_float_set(fvFG_Temp2, iGlobPos,   0.0);
				}	}
				else
				{ 	gsl_vector_float_set(fvFG_Temp0, iGlobPos,   0.0);  //the strike slip component            
					gsl_vector_float_set(fvFG_Temp1, iGlobPos,   0.0);  //the dip slip component  
					gsl_vector_float_set(fvFG_Temp2, iGlobPos,   0.0);
		}	}   }	    
		MPI_Allreduce(MPI_IN_PLACE, &iTemp0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		//--------------------------------------------------------------
		if (iTemp0 == 1)
		{	//----------------------------------------------------------	
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);	
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp2->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
			//---------------------------
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SS, fvFG_Temp0, 0.0, fvFL_Temp0); 				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SD, fvFG_Temp0, 0.0, fvFL_Temp1); 				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SO, fvFG_Temp0, 0.0, fvFL_Temp2); 
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);									gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);									gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DS, fvFG_Temp1, 0.0, fvFL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DD, fvFG_Temp1, 0.0, fvFL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DO, fvFG_Temp1, 0.0, fvFL_Temp2);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);									gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);									gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_OS, fvFG_Temp2, 0.0, fvFL_Temp0);				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_OD, fvFG_Temp2, 0.0, fvFL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_OO, fvFG_Temp2, 0.0, fvFL_Temp2);
			gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);									gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);									gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
			//----------------------------------------------------------
			if ((MD.iBPNum > 0) && (MD.iUsePSeis == 1))
        	{	MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp0->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp1->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);	
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvBG_Temp2->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);	
				//---------------------------
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SS, fvFG_Temp0, 0.0, fvBL_Temp0); 			gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SD, fvFG_Temp0, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SO, fvFG_Temp0, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DS, fvFG_Temp1, 0.0, fvBL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DD, fvFG_Temp1, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DO, fvFG_Temp1, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_OS, fvFG_Temp2, 0.0, fvBL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.FB_OD, fvFG_Temp2, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_OO, fvFG_Temp2, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				//---------------------------
				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_SS, fvBG_Temp0, 0.0, fvFL_Temp0); 			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_SD, fvBG_Temp0, 0.0, fvFL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_SO, fvBG_Temp0, 0.0, fvFL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_DS, fvBG_Temp1, 0.0, fvFL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_DD, fvBG_Temp1, 0.0, fvFL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_DO, fvBG_Temp1, 0.0, fvFL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_OS, fvBG_Temp2, 0.0, fvFL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.BF_OD, fvBG_Temp2, 0.0, fvFL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BF_OO, fvBG_Temp2, 0.0, fvFL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				//---------------------------
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SS, fvBG_Temp0, 0.0, fvBL_Temp0); 			gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SD, fvBG_Temp0, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_SO, fvBG_Temp0, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DS, fvBG_Temp1, 0.0, fvBL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DD, fvBG_Temp1, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_DO, fvBG_Temp1, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OS, fvBG_Temp2, 0.0, fvBL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OD, fvBG_Temp2, 0.0, fvBL_Temp1);				gsl_blas_sgemv(CblasTrans, 1.0, K.BB_OO, fvBG_Temp2, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);								gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);									gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				//---------------------------
		}	}
		//--------------------------------------------------------------
		EQ.iStillOn = 0;
        for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
        { 	if (TR.ivFL_StabT[i] == 1 )
           	{   fTemp0   = gsl_vector_float_get(TR.fvFL_CurStrsH, i);
				fTemp1   = gsl_vector_float_get(TR.fvFL_CurStrsV, i);
				fTemp2   = sqrtf(fTemp0*fTemp0 + fTemp1*fTemp1);
				fTemp3   = fTemp2 - TR.fvFL_CurFric[i]*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //this is the excess stress (is excess if value > 0) I currently have with respect to curr strength
            	
				if (fTemp3  >  MD.fCutStrss) //I have excess stress => can start an EQ
				{	EQ.iStillOn  = 1;      		TR.ivFG_Activated[i + MD.ivF_START[MD.iRANK]]  = 1;					TR.ivFG_Ptch_t0[i + MD.ivF_START[MD.iRANK]]  = 0;  	    //Earthquake begins; I set this patch to "activated" and write its start time (==0) 
        }   }	}   
        MPI_Allreduce(MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);	 
		//------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------
		if (EQ.iStillOn == 1)
        {	//-----------------------------------------
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Activated, MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD); //collect and synchronize which patches have been activated already
			MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Ptch_t0,   MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD); //also collect and synchronize the patch start times (time of patch "activation")

			gsl_vector_float_memcpy(TR.fvFL_B4_StrsH, TR.fvFL_CurStrsH); //store the pre-EQ stress on fault for stress drop;
			gsl_vector_float_memcpy(TR.fvFL_B4_StrsV, TR.fvFL_CurStrsV); //put it into "B4_Strs" vector
            gsl_vector_float_memcpy(TR.fvFL_B4_StrsN, TR.fvFL_CurStrsN);
            gsl_vector_float_set_zero(fvFG_Temp0);       	gsl_vector_float_set_zero(fvFG_Temp1);             	gsl_vector_float_set_zero(fvFG_Temp2); 
			gsl_vector_float_set_zero(EQ.fvG_EQslipH);		gsl_vector_float_set_zero(EQ.fvG_EQslipV);	

			for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
        	{   EQ.ivL_ActPtchID[i] = 0;					EQ.ivL_t0ofPtch[i]     = 0;							EQ.ivL_StabType[i]  = 0;			
				EQ.fvL_PtchSlpH[i]  = 0.0;					EQ.fvL_PtchSlpV[i]     = 0.0;						EQ.fvL_PtchDTau[i]  = 0.0;
                TR.fvFL_B4_Fric[i] = TR.fvFL_CurFric[i];	TR.fvFL_TempRefFric[i] = TR.fvFL_CurFric[i];		TR.fvFL_AccumSlp[i] = 0.0;
            }     
			for (i = 0; i < MD.iSIZE; i++)          {      	EQ.ivR_WrtStrtPos[i] = 0;   						}            
			for (i = 0; i < MD.iMaxMRFlgth; i++)	{		EQ.fvM_MRFvals[i]    = 0.0;							}		
			EQ.iEndCntr   = 0;      EQ.iActFPNum = 0;		EQ.iCmbFPNum = 0;           EQ.iMRFLgth = -1;		
			EQ.iTotlRuptT =-1;		EQ.fSeisPot  = 0.0;		EQ.fMaxSlip  = 0.0;		    EQ.fMaxDTau = 0.0;
        	//------------------------------------------------------------------------------------------
			//------------------------------------------------------------------------------------------  
    		// EARTHQUAKE ITERATION LOOP STARTS
        	while (EQ.iStillOn == 1)
        	{   
            	EQ.iTotlRuptT++;    //continuously count time since initiation, set "ongoing" to FALSE => only if more slip on patches is added, it gets to be reset            
                EQ.iStillOn = 0; 
				//----------------------------------------------------
            	for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
            	{  	
					iGlobPos = i + MD.ivF_START[MD.iRANK];    

					if (TR.ivFG_Activated[iGlobPos] == 1)
                    {	fTemp2             = sqrtf(gsl_vector_float_get(fvFG_Temp0, iGlobPos) *gsl_vector_float_get(fvFG_Temp0, iGlobPos)  +  gsl_vector_float_get(fvFG_Temp1, iGlobPos) *gsl_vector_float_get(fvFG_Temp1, iGlobPos));
						//--------------------------------------
						TR.fvFL_CurFric[i] = GetUpdatedFriction(TR.fvFL_B4_Fric[i], TR.fvFL_TempRefFric[i], TR.fvFL_CurFric[i], TR.fvFL_DynFric[i], TR.ivFL_FricLaw[i], TR.ivFL_StabT[i], TR.fvFL_CurDcVal[i], TR.fvFL_AccumSlp[i], fTemp2, MD.fHealFact);

					//	TR.fvFL_CurFric[i] = TR.fvFL_DynFric[i];
						//--------------------------------------
					}
					//-----------------------------------------------------------------------
					// determine amount of excess stress (if any available) and the corresponding slip amount
                	fTemp0 = gsl_vector_float_get(TR.fvFL_CurStrsH, i);
					fTemp1 = gsl_vector_float_get(TR.fvFL_CurStrsV, i);
					fTemp2 = sqrtf(fTemp0*fTemp0 + fTemp1*fTemp1);
					fTemp3 = fTemp2 - TR.fvFL_CurFric[i]*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //this is the excess stress (is excess if value > 0) I currently have with respect to curr/updated strength
					//--------------------------
                	if (   (TR.ivFG_Activated[iGlobPos] == 1) || ( (TR.ivFG_Activated[iGlobPos] == 0) && (fTemp3 > MD.fCutStrss) )   )
                	{  	
						if (TR.ivFG_Activated[iGlobPos] == 0)                  {    TR.ivFG_Activated[iGlobPos] = 1;           TR.ivFG_Ptch_t0[iGlobPos] = EQ.iTotlRuptT;           }
		
                    	if (fTemp3 > MD.fCutStrss)                       
                    	{  	EQ.iStillOn = 1;                       
                        	EQ.iMRFLgth = EQ.iTotlRuptT;   						
							fTemp4      = -1.0*(fTemp3/fTemp2*fTemp0) / gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[i], i); // slip amount to release excess horizontal shear stress
                			fTemp5      = -1.0*(fTemp3/fTemp2*fTemp1) / gsl_matrix_float_get(K.FF_DD, TR.ivFL_SelfLoc[i], i); // slip amount to release excess vertical shear stress 
							fTemp6      = sqrtf(fTemp4*fTemp4 +fTemp5*fTemp5); //the slip amount needed to release the excess...  
							
							if ((TR.fvFL_AccumSlp[i] == 0.0) && (TR.ivFL_StabT[i] !=3))		{			TR.fvFL_TempRefFric[i] = TR.fvFL_CurFric[i];					}
							
							TR.fvFL_AccumSlp[i] += fTemp6;

							gsl_vector_float_set(fvFG_Temp0,   iGlobPos, fTemp4); //also put that slip into "total slip at source patch" list => for event slip
							gsl_vector_float_set(fvFG_Temp1,   iGlobPos, fTemp5); //so, this will go into EVENTslip...

							if (EQ.iMRFLgth < MD.iMaxMRFlgth)	//write out the moment-rate-function, but only to a maximum of maxMRFlength => can still do the earthquake but cutoff the MRF...
                        	{   EQ.fvM_MRFvals[EQ.iTotlRuptT]    += fTemp6 *TR.fvFL_Area[i] *MD.fShearMod; //this is done locally here, will be synchronized once the whole EQ is over!
                    	}   }
                    	else //this is necessary here b/c I can get into that if-statement b/c I got activated but might not have any amount of slip left => set value to zero...
                    	{ 	TR.fvFL_AccumSlp[i] = 0.0;
							gsl_vector_float_set(fvFG_Temp0,   iGlobPos, 0.0);
                            gsl_vector_float_set(fvFG_Temp1,   iGlobPos, 0.0);
				    }   }
                    else
                    {  	TR.fvFL_AccumSlp[i] = 0.0;
						gsl_vector_float_set(fvFG_Temp0,   iGlobPos, 0.0); //all this is necessary b/c I cannot reset these TempVectors to zero in a single statement; reason is that it would
				        gsl_vector_float_set(fvFG_Temp1,   iGlobPos, 0.0); //be needed/contains info for case that velocity weakening is used...; have to make sure that both tempvectors are fully defined/rewritten "by hand"
				}	}
				//-----------------------------------------------------------------------
				MPI_Allreduce( MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Activated, MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD);	
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Ptch_t0,   MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD); 	
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data,  MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT,    MPI_COMM_WORLD); //temporary storage for event slip (horizontal)
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data,  MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT,    MPI_COMM_WORLD); //and vertical	
				gsl_vector_float_add(EQ.fvG_EQslipH, fvFG_Temp0);
				gsl_vector_float_add(EQ.fvG_EQslipV, fvFG_Temp1);
            	//-----------------------------------------------------------------------
				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SS, fvFG_Temp0, 0.0, fvFL_Temp0); 			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SD, fvFG_Temp0, 0.0, fvFL_Temp1); 			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_SO, fvFG_Temp0, 0.0, fvFL_Temp2); 
				gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);								gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);								gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DS, fvFG_Temp1, 0.0, fvFL_Temp0);			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DD, fvFG_Temp1, 0.0, fvFL_Temp1);			gsl_blas_sgemv(CblasTrans, 1.0, K.FF_DO, fvFG_Temp1, 0.0, fvFL_Temp2);
				gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);								gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp1);								gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
				//-----------------------------------------------------------------------
			}			
			//-----------------------------------------------------------------------
			fTemp0       = 0.0; //this is temp/test magnitude (i.e., combined seismic potential)
			//----------------------
    		for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
        	{  	iGlobPos = i + MD.ivF_START[MD.iRANK]; 
				if (TR.ivFG_Activated[iGlobPos] == 1)     
            	{ 	fTemp1  = gsl_vector_float_get(EQ.fvG_EQslipH, iGlobPos); //horizontal slip of patch
					fTemp2  = gsl_vector_float_get(EQ.fvG_EQslipV, iGlobPos); //vertical slip of patch
	
					if ( sqrtf(fTemp1*fTemp1 +fTemp2*fTemp2) >= fMinSlipForCat) //minimum slip amount to be put into the list
					{	fTemp0  += sqrtf(fTemp1*fTemp1 +fTemp2*fTemp2)*TR.fvFL_Area[i]; //this is seis. potential of that patch (added to other patches from same ranke)
        	}   } 	}
			//-----------------------------------------------------------------------
			MPI_Allreduce(MPI_IN_PLACE,		&fTemp0, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //combined seismic potential that is released in this event
			//------------------------------------	
            fTemp1 = (log10f(fTemp0*MD.fShearMod)-9.1)/1.5; 	
           	//-----------------------------------------------------------------------
			//-----------------------------------------------------------------------
			if (  (MD.iUseProp == 1) && (fTemp1 >= MD.fMinMag4Prop) )
			{	//------------------------------------		 
                if ((iPlot2Screen == 1) && (MD.iRANK == 0))				{	fprintf(stdout,"Test Mag:  %5.2f   ",fTemp1);					}
				//------------------------------------
                gsl_vector_float_memcpy(TR.fvFL_CurStrsH, TR.fvFL_B4_StrsH); //reset stress to pre-EQ kind...and also all the other things to zero i.e., pre-rupture conditions.
			    gsl_vector_float_memcpy(TR.fvFL_CurStrsV, TR.fvFL_B4_StrsV);
                gsl_vector_float_memcpy(TR.fvFL_CurStrsN, TR.fvFL_B4_StrsN);
                //------------------------------------		
        	    gsl_matrix_float_set_zero(EQ.fmG_STF_H);			gsl_matrix_float_set_zero(EQ.fmG_STF_V);
                gsl_matrix_int_set_zero(TR.imFGL_NextP);			gsl_matrix_int_set_zero(TR.imFGL_NextS);
				gsl_vector_float_set_zero(EQ.fvG_EQslipH);			gsl_vector_float_set_zero(EQ.fvG_EQslipV);
				gsl_vector_float_set_zero(fvFG_Temp0);              gsl_vector_float_set_zero(fvFG_Temp1);                  gsl_vector_float_set_zero(fvFG_Temp2);  

                //------------------------------------		
                for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 	{   TR.fvFL_CurFric[i]  = TR.fvFL_B4_Fric[i];				TR.fvFL_TempRefFric[i] = TR.fvFL_B4_Fric[i];			TR.fvFL_AccumSlp[i] = 0.0;				}
                //--------------------------------------------------------------    
                for (i = 0; i < MD.iFPNum; i++) 		{			TR.ivFG_Ptch_t0[i]   = 0;			TR.ivFG_Activated[i] = 0;				}
			    for (i = 0; i < MD.iMaxMRFlgth; i++)	{			EQ.fvM_MRFvals[i]    = 0.0;													}		
   			    //--------------------------------------------------------------  
                EQ.iMRFLgth = -1;									EQ.iTotlRuptT =-1;		
			    //--------------------------------------------------------------    
				EQ.iStillOn = 0;
				for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
				{ 	if (TR.ivFL_StabT[i] == 1 )
					{   fTemp0   = gsl_vector_float_get(TR.fvFL_CurStrsH, i);
						fTemp1   = gsl_vector_float_get(TR.fvFL_CurStrsV, i);
						fTemp2   = sqrtf(fTemp0*fTemp0 + fTemp1*fTemp1);
						fTemp3   = fTemp2 - TR.fvFL_CurFric[i]*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //this is the excess stress (is excess if value > 0) I currently have with respect to curr strength
						
						if (fTemp3  >  MD.fCutStrss) //I have excess stress => how much slip is that?
						{	EQ.iStillOn  = 1;      		TR.ivFG_Activated[i + MD.ivF_START[MD.iRANK]]  = 1;					TR.ivFG_Ptch_t0[i + MD.ivF_START[MD.iRANK]]  = 0;  	    
				}   }	}   
        		MPI_Allreduce(MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);	 
				
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Activated, MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD);
				MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Ptch_t0,   MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD);      

 				//--------------------------------------------------------------    
				while (EQ.iStillOn == 1)
        		{	
					EQ.iTotlRuptT++;                
                	EQ.iStillOn = 0; //continuously count time since initiation, set "ongoing" to FALSE => only if more slip on patches is added, it gets to be reset
					//----------------------------------------------------
					for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
            		{  	
						iGlobPos = i + MD.ivF_START[MD.iRANK];    

						if (TR.ivFG_Activated[iGlobPos] == 1)
                   	 	{	fTemp2   = sqrtf(gsl_vector_float_get(fvFG_Temp0, iGlobPos) *gsl_vector_float_get(fvFG_Temp0, iGlobPos)  +  gsl_vector_float_get(fvFG_Temp1, iGlobPos) *gsl_vector_float_get(fvFG_Temp1, iGlobPos));
							//--------------------------------------
							TR.fvFL_CurFric[i] = GetUpdatedFriction(TR.fvFL_B4_Fric[i], TR.fvFL_TempRefFric[i], TR.fvFL_CurFric[i], TR.fvFL_DynFric[i], TR.ivFL_FricLaw[i], TR.ivFL_StabT[i], TR.fvFL_CurDcVal[i], TR.fvFL_AccumSlp[i], fTemp2, MD.fHealFact);
							//--------------------------------------
						//	TR.fvFL_CurFric[i] = TR.fvFL_DynFric[i];
						}
						//-----------------------------------------------------------------------
						// determine amount of excess stress (if any available) and the corresponding slip amount
                		fTemp0 = gsl_vector_float_get(TR.fvFL_CurStrsH, i);
						fTemp1 = gsl_vector_float_get(TR.fvFL_CurStrsV, i);
						fTemp2 = sqrtf(fTemp0*fTemp0 + fTemp1*fTemp1);
						fTemp3 = fTemp2 - TR.fvFL_CurFric[i]*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //this is the excess stress (is excess if value > 0) I currently have with respect to curr/updated strength
						//--------------------------
                		if (   (TR.ivFG_Activated[iGlobPos] == 1) || ( (TR.ivFG_Activated[iGlobPos] == 0) && (fTemp3 > MD.fCutStrss) )   )
                		{  
							STFcnt = EQ.iTotlRuptT -TR.ivFG_Ptch_t0[iGlobPos];
							if (TR.ivFG_Activated[iGlobPos] == 0)                  {    STFcnt = 0;     	TR.ivFG_Activated[iGlobPos] = 1;           TR.ivFG_Ptch_t0[iGlobPos] = EQ.iTotlRuptT;           }                                                        

                        	iTemp0 = STFcnt%MD.iMaxSTFlgth;   //gives me remainder of wPos/MaxSTFlength => allows me to loop/overwrite the STF vectors 
		
                    		if (fTemp3 > MD.fCutStrss)                       
                    		{  	EQ.iStillOn = 1;                       
                        		EQ.iMRFLgth = EQ.iTotlRuptT;   						
								fTemp4      = -1.0*(fTemp3/fTemp2*fTemp0) / gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[i], i); // slip amount to release excess horizontal shear stress
                				fTemp5      = -1.0*(fTemp3/fTemp2*fTemp1) / gsl_matrix_float_get(K.FF_DD, TR.ivFL_SelfLoc[i], i); // slip amount to release excess vertical shear stress 
								fTemp6      = sqrtf(fTemp4*fTemp4 +fTemp5*fTemp5); //the slip amount needed to release the excess...  
							
								if ((TR.fvFL_AccumSlp[i] == 0.0) && (TR.ivFL_StabT[i] !=3))		{			TR.fvFL_TempRefFric[i] = TR.fvFL_CurFric[i];					}
							
								TR.fvFL_AccumSlp[i] += fTemp6;
								gsl_matrix_float_set(EQ.fmG_STF_H, iGlobPos, iTemp0, fTemp4); //put into the STFmatrix, which has dimenstion [FPNUM, MaxSTFlgth]
								gsl_matrix_float_set(EQ.fmG_STF_V, iGlobPos, iTemp0, fTemp5);
								gsl_vector_float_set(fvFG_Temp0,   iGlobPos, fTemp4); //also put that slip into "total slip at source patch" list => for event slip
								gsl_vector_float_set(fvFG_Temp1,   iGlobPos, fTemp5); //so, this will go into EVENTslip...

								if (EQ.iMRFLgth < MD.iMaxMRFlgth)	//write out the moment-rate-function, but only to a maximum of maxMRFlength => can still do the earthquake but cutoff the MRF...
                        		{   EQ.fvM_MRFvals[EQ.iTotlRuptT]    += fTemp6 *TR.fvFL_Area[i] *MD.fShearMod; //this is done locally here, will be synchronized once the whole EQ is over!
                    		}   }
                    		else //this is necessary here b/c I can get into that if-statement b/c I got activated but might not have any amount of slip left => set value to zero...
                    		{ 	TR.fvFL_AccumSlp[i] = 0.0;
							    gsl_matrix_float_set(EQ.fmG_STF_H, iGlobPos, iTemp0, 0.0);
						        gsl_matrix_float_set(EQ.fmG_STF_V, iGlobPos, iTemp0, 0.0);// it is necessary b/c that position would remain its previous value, and I am looping over this stuff => would reuse previous slip amount   
								gsl_vector_float_set(fvFG_Temp0,   iGlobPos, 0.0);
                            	gsl_vector_float_set(fvFG_Temp1,   iGlobPos, 0.0);
				   	 	}   }
                    	else
                    	{  	TR.fvFL_AccumSlp[i] = 0.0;
						    gsl_vector_float_set(fvFG_Temp0,   iGlobPos, 0.0); //all this is necessary b/c I cannot reset these TempVectors to zero in a single statement; reason is that it would
				        	gsl_vector_float_set(fvFG_Temp1,   iGlobPos, 0.0); //be needed/contains info for case that velocity weakening is used...; have to make sure that both tempvectors are fully defined/rewritten "by hand"
                	}   }
					//-----------------------------------------------------------------------
					MPI_Allreduce( MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Activated,  MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD);	
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, TR.ivFG_Ptch_t0,    MD.ivF_OFFSET, MD.ivF_START, MPI_INT, MPI_COMM_WORLD); 	
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, EQ.fmG_STF_H->data, MD.ivF_OFFSET2, MD.ivF_START2, MPI_FLOAT, MPI_COMM_WORLD);
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, EQ.fmG_STF_V->data, MD.ivF_OFFSET2, MD.ivF_START2, MPI_FLOAT, MPI_COMM_WORLD);
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp0->data,   MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT,    MPI_COMM_WORLD); //temporary storage for event slip (horizontal)
					MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, fvFG_Temp1->data,   MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT,    MPI_COMM_WORLD); //and vertical	
					gsl_vector_float_add(EQ.fvG_EQslipH, fvFG_Temp0);
					gsl_vector_float_add(EQ.fvG_EQslipV, fvFG_Temp1);
					//-----------------------------------------------------------------------
					for (i = 0; i < MD.iFPNum; i++) //goes through all the SOURCES!
					{	if (TR.ivFG_Activated[i] == 1) //if that source did actually get already activated in the current event... otherwise I can skip over it...
						{	//-----------------------------------------------------------------------	
							STFcnt  = EQ.iTotlRuptT - TR.ivFG_Ptch_t0[i]; //the "time" since the rupture started a source
							for (j = 0; j < MD.ivF_OFFSET[MD.iRANK]; j++) //this is my local receiver...
							{	//-----------------------------------------------------------------------
								//first the "P-wave" signal (PROjection of slip vector onto SrcRcvVect)
								iTemp0   = gsl_matrix_int_get(TR.imFGL_NextP, i, j); 		//"pointer" location for the STFmatrix => to know what part of the source STF a given receiver has already seen/used
								iTemp1   = gsl_matrix_int_get(TR.imFGL_TTP,   i, j); 		//for stress change calculation; ACTUALLY, better way is not pointer but overall how much of the STF has already been processed... (this value, NextP DOES NOT LOOP AROUND)
								for (k = iTemp0; k <= (STFcnt - iTemp1); k++)				//2nd part is negative if signal has not arrived, then the loop is skipped
								{	iTemp2 = k%MD.iMaxSTFlgth; 								//this is the actual "read position" => where to look in the STFmatrix
									fTemp0 = gsl_matrix_float_get(EQ.fmG_STF_H, i, iTemp2); //this is the next horizontal slip
									fTemp1 = gsl_matrix_float_get(EQ.fmG_STF_V, i, iTemp2); //this is the next vertical slip
							
									fTemp2 = gsl_vector_float_get(TR.fvFL_CurStrsH, j); 	//I'm doing it this way b/c I don't know how to do a " += " with a gsl_vector...
									fTemp3 = gsl_vector_float_get(TR.fvFL_CurStrsV, j);
									fTemp4 = gsl_vector_float_get(TR.fvFL_CurStrsN, j); 	//remember, there is no actual (static) motion in normal direction, but a transient slip signal when vP/vS is used...
                          	
									if ((j + MD.ivF_START[MD.iRANK]) == i)
                            		{   fTemp10 = fTemp0*gsl_matrix_float_get(K.FF_SS, i, j) + fTemp1*gsl_matrix_float_get(K.FF_DS, i, j); //this is on "self", => use actual slip from STF to determine new stress state i.e., the stress change due to slip
										fTemp11 = fTemp0*gsl_matrix_float_get(K.FF_SD, i, j) + fTemp1*gsl_matrix_float_get(K.FF_DD, i, j);
										fTemp12 = 0.0;
									}
									else    
						    		{	//http://sites.science.oregonstate.edu/math/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/dotprod/dotprod.html
										//https://en.wikipedia.org/wiki/Vector_projection
										fTemp9 = fTemp0*gsl_matrix_float_get(TR.fmFGL_SrcRcvH, i, j) + fTemp1*gsl_matrix_float_get(TR.fmFGL_SrcRcvV, i, j); //"normally", this should also include the normal component => but that one is always zero! => don't need it herethis is basically a vector projection, using scalar projection using dot product
										//fTemp9 is the length/projection of the slip vector (fTemp0 and fTemp1 being the horizontal and vertical component) onto the SrcRcv vector (which sits at source and points to receiver, is normalized; ALSO: that vector is in local coordinates of the source patch! -done in "AddMoreParas")	
										fTemp5 = fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvH, i, j); //this part here is actually a really cool step; fTemp9 is length of slip amount in direction of source-receiver vector (can be negative)
										fTemp6 = fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvV, i, j); //this "length" is multiplied with local orientation of source-receiver vector to get the transient slip values to be used...
										fTemp7 = fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvN, i, j); //these are currently slip values... => now need to convert to stress change

										fTemp10= fTemp5*gsl_matrix_float_get(K.FF_SS, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DS, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OS, i, j);
										fTemp11= fTemp5*gsl_matrix_float_get(K.FF_SD, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DD, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OD, i, j);
										fTemp12= fTemp5*gsl_matrix_float_get(K.FF_SO, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DO, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OO, i, j);
									}	
									gsl_vector_float_set(TR.fvFL_CurStrsH, j, (fTemp2+fTemp10));
									gsl_vector_float_set(TR.fvFL_CurStrsV, j, (fTemp3+fTemp11));
									gsl_vector_float_set(TR.fvFL_CurStrsN, j, (fTemp4+fTemp12));
								}
								iTemp0 = (STFcnt - iTemp1) >= 0 ? (STFcnt - iTemp1 +1) : 0; 
								gsl_matrix_int_set(TR.imFGL_NextP, i, j, iTemp0);
								//-----------------------------------------------------------------------
								//second the "S-wave" signal (REjection of slip vector onto SrcRcvVect)
								iTemp0   = gsl_matrix_int_get(TR.imFGL_NextS, i, j); 
								iTemp1   = gsl_matrix_int_get(TR.imFGL_TTS,   i, j);
								for (k = iTemp0; k <= (STFcnt - iTemp1); k++)
								{	iTemp2 = k%MD.iMaxSTFlgth; 								//this is the actual "read position" => where to look in the STFmatrix
    								fTemp0 = gsl_matrix_float_get(EQ.fmG_STF_H, i, iTemp2); //this is the next horizontal slip
									fTemp1 = gsl_matrix_float_get(EQ.fmG_STF_V, i, iTemp2); //this is the next vertical slip
							
									fTemp2 = gsl_vector_float_get(TR.fvFL_CurStrsH, j); 	//I'm doing it this way b/c I don't know how to do a " += " with a gsl_vector...
									fTemp3 = gsl_vector_float_get(TR.fvFL_CurStrsV, j);
									fTemp4 = gsl_vector_float_get(TR.fvFL_CurStrsN, j); 	//remember, there is no actual (static) motion in normal direction, but a transient slip signal when vP/vS is used...
                          	
									if (j + MD.ivF_START[MD.iRANK] == i)
                            		{	fTemp10= 0.0;
										fTemp11= 0.0;
										fTemp12= 0.0;	
									}
									else
									{	fTemp9 = fTemp0*gsl_matrix_float_get(TR.fmFGL_SrcRcvH, i, j) + fTemp1*gsl_matrix_float_get(TR.fmFGL_SrcRcvV, i, j);
								
										fTemp5 = fTemp0 - fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvH, i, j);
										fTemp6 = fTemp1 - fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvV, i, j);
										fTemp7 = 0.0    - fTemp9*gsl_matrix_float_get(TR.fmFGL_SrcRcvN, i, j);

										fTemp10= fTemp5*gsl_matrix_float_get(K.FF_SS, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DS, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OS, i, j);
										fTemp11= fTemp5*gsl_matrix_float_get(K.FF_SD, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DD, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OD, i, j);
										fTemp12= fTemp5*gsl_matrix_float_get(K.FF_SO, i, j) + fTemp6*gsl_matrix_float_get(K.FF_DO, i, j) + fTemp7*gsl_matrix_float_get(K.FF_OO, i, j);
									}
									gsl_vector_float_set(TR.fvFL_CurStrsH, j, (fTemp2+fTemp10));
									gsl_vector_float_set(TR.fvFL_CurStrsV, j, (fTemp3+fTemp11));
									gsl_vector_float_set(TR.fvFL_CurStrsN, j, (fTemp4+fTemp12));
								}
                     			iTemp0 = (STFcnt - iTemp1) >= 0 ? (STFcnt - iTemp1 +1) : 0;   
								gsl_matrix_int_set(TR.imFGL_NextS, i, j, iTemp0);
								//-----------------------------------------------------------------------
					}	}	}	
					if (EQ.iStillOn == 0) // make sure that all the signal is out of the system => even if no slip occurred in last step on any patch; there may still be stress in the system that has not reached a receiver => wait until they all got their share 
            		{ 	EQ.iEndCntr += 1;
                	   	if (EQ.iEndCntr < MD.iMaxSTFlgth)    	{  		EQ.iStillOn     = 1;         } 		
               		}     
            		else 
					{  	EQ.iEndCntr  = 0;              
			}	}	}
			//-----------------------------------------------------------------------
			//-----------------------------------------------------------------------
			fTemp0       = 0.0;	
			EQ.iMRFLgth  = (EQ.iMRFLgth < MD.iMaxMRFlgth) ? EQ.iMRFLgth : MD.iMaxMRFlgth;
			//----------------------
    		for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
        	{  	iGlobPos = i + MD.ivF_START[MD.iRANK]; 
				if (TR.ivFG_Activated[iGlobPos] == 1)     
            	{ 	
					fTemp1  = gsl_vector_float_get(EQ.fvG_EQslipH, iGlobPos); //horizontal slip of patch
					fTemp2  = gsl_vector_float_get(EQ.fvG_EQslipV, iGlobPos); //vertical slip of patch
					fTemp3  = gsl_vector_float_get(TR.fvFL_B4_StrsH, i);				fTemp4  = gsl_vector_float_get(TR.fvFL_B4_StrsV, i);
					fTemp5  = gsl_vector_float_get(TR.fvFL_CurStrsH, i);				fTemp6  = gsl_vector_float_get(TR.fvFL_CurStrsV, i);

					if ( sqrtf(fTemp1*fTemp1 +fTemp2*fTemp2) >= fMinSlipForCat) //minimum slip amount to be put into the list
					{	fTemp0                        += sqrtf(fTemp1*fTemp1 +fTemp2*fTemp2)*TR.fvFL_Area[i]; //this is seis. potential of that patch (added to other patches from same ranke)
						EQ.ivL_ActPtchID[EQ.iActFPNum] = iGlobPos;
						EQ.ivL_t0ofPtch[EQ.iActFPNum]  = TR.ivFG_Ptch_t0[iGlobPos];
						EQ.fvL_PtchSlpH[EQ.iActFPNum]  = fTemp1;
      					EQ.fvL_PtchSlpV[EQ.iActFPNum]  = fTemp2;
               			EQ.fvL_PtchDTau[EQ.iActFPNum]  = sqrtf(fTemp3*fTemp3 +fTemp4*fTemp4) - sqrtf(fTemp5*fTemp5 +fTemp6*fTemp6); //stress DROP is therefore POSITIVE
                   	 	EQ.ivL_StabType[EQ.iActFPNum]  = TR.ivFL_StabT[i];
						EQ.fMaxDTau                    =      (EQ.fvL_PtchDTau[EQ.iActFPNum] > EQ.fMaxDTau) ?      EQ.fvL_PtchDTau[EQ.iActFPNum] : EQ.fMaxDTau;
						EQ.fMaxSlip                    = (sqrt(fTemp1*fTemp1 +fTemp2*fTemp2) > EQ.fMaxSlip) ? sqrt(fTemp1*fTemp1 +fTemp2*fTemp2) : EQ.fMaxSlip;  
						EQ.iActFPNum++;  
        	}   } 	}
			//-----------------------------------------------------------------------
			MPI_Allreduce(&EQ.iActFPNum, 	&EQ.iCmbFPNum,    1      , MPI_INT,   MPI_SUM, MPI_COMM_WORLD);
        	MPI_Allreduce(MPI_IN_PLACE,		&fTemp0,          1      , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //combined seismic potential that is released in this event
        	MPI_Allreduce(MPI_IN_PLACE,		&EQ.iMRFLgth,     1      , MPI_INT,   MPI_MAX, MPI_COMM_WORLD); //the length of the MRF
        	MPI_Allreduce(MPI_IN_PLACE,		&EQ.fMaxSlip,     1      , MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); //this is maximum slip
        	MPI_Allreduce(MPI_IN_PLACE,		&EQ.fMaxDTau,     1      , MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); //this is maximum stress drop
        	MPI_Allreduce(MPI_IN_PLACE, EQ.fvM_MRFvals,MD.iMaxMRFlgth, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
         	//------------------------------------			
            //------------------------------------		
            fTemp1 = (log10f(fTemp0*MD.fShearMod)-9.1)/1.5; 
            //------------------------------------		
            //------------------------------------		
			if ( (MD.iBPNum > 0) && (MD.iUsePSeis == 1) )
			{	gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SS, EQ.fvG_EQslipH, 0.0, fvBL_Temp0); 
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SD, EQ.fvG_EQslipH, 0.0, fvBL_Temp1);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_SO, EQ.fvG_EQslipH, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);
				gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);
				gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DS, EQ.fvG_EQslipV, 0.0, fvBL_Temp0);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DD, EQ.fvG_EQslipV, 0.0, fvBL_Temp1);
				gsl_blas_sgemv(CblasTrans, 1.0, K.FB_DO, EQ.fvG_EQslipV, 0.0, fvBL_Temp2);
				gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);			
				gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp1);
				gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp2);
			}
			//-----------------------------------------------------------------------					    
            if (EQ.iCmbFPNum >= iMinPtch4Cat)
			{	
                MD.iEQcntr++;  

                MPI_Allgather(&EQ.iActFPNum, 1, MPI_INT, EQ.ivR_WrtStrtPos, 1, MPI_INT,  MPI_COMM_WORLD);   
				for (i = 1;     i < MD.iSIZE; i++)           {    EQ.ivR_WrtStrtPos[i] += EQ.ivR_WrtStrtPos[i-1];             }    
        		for (i = (MD.iSIZE-1); i > 0; i--)           {    EQ.ivR_WrtStrtPos[i]  = EQ.ivR_WrtStrtPos[i-1];             }    
        		EQ.ivR_WrtStrtPos[0]  = 0;
				//-----------------------------------------------------------------------
				if ((MD.iRANK == 0)&&(iPlot2Screen == 1) && (fTemp1 > 4.5))
				{   iTemp0 = (int)fTemp1 -4; //just for nicer plotting...
					fprintf(stdout,"%6d  Earthquake time %5.2f    act patches: %5d   MRF length: %4d   MaxSlip: %4.2f     MaxStressDrop: %4.2f    ",MD.iEQcntr, MD.fTimeYears,  EQ.iCmbFPNum, EQ.iMRFLgth, EQ.fMaxSlip, EQ.fMaxDTau);
					for (i = 0; i < iTemp0; i++) 	{		fprintf(stdout,"    ");			}
					fprintf(stdout,"Magn: %3.2f\n",fTemp1);              
				}
      			//-----------------------------------------------------------------------
        		if (MD.iRANK == 0)
        		{	MPI_File_write_at(fp_MPIOUT,             0,                               &MD.iEQcntr,        1,       MPI_INT,   &STATUS);
				    MPI_File_write_at(fp_MPIOUT,   OFFSETall,                                 &MD.fTimeYears,     1,       MPI_FLOAT, &STATUS);  //Earthquake time
            	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +  sizeof(float),                &fTemp1,            1,       MPI_FLOAT, &STATUS);  //Earthquake magnitude
            	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +2*sizeof(float),                &EQ.iCmbFPNum,      1,       MPI_INT,   &STATUS);  //#of fault patches participating in EQ
            	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +2*sizeof(float) +1*sizeof(int), &EQ.iMRFLgth,       1,       MPI_INT,   &STATUS);  //length of moment rate function "time steps"
            	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +2*sizeof(float) +2*sizeof(int), EQ.fvM_MRFvals, EQ.iMRFLgth, MPI_FLOAT, &STATUS);  //moment rate function values..
				}
				OFFSETall += 2*sizeof(int) +(2 +EQ.iMRFLgth)*sizeof(float);
				//--------------------------------------------
				MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(int),  EQ.ivL_ActPtchID,EQ.iActFPNum, MPI_INT,   &STATUS);
        	    OFFSETall += EQ.iCmbFPNum*sizeof(int);
        	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(int),  EQ.ivL_t0ofPtch, EQ.iActFPNum, MPI_INT,   &STATUS);
        	    OFFSETall += EQ.iCmbFPNum*sizeof(int);
        	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(float),EQ.fvL_PtchDTau, EQ.iActFPNum, MPI_FLOAT, &STATUS);
        	    OFFSETall += EQ.iCmbFPNum*sizeof(float);
        	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(float),EQ.fvL_PtchSlpH, EQ.iActFPNum, MPI_FLOAT, &STATUS);
        	    OFFSETall += EQ.iCmbFPNum*sizeof(float);
        	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(float),EQ.fvL_PtchSlpV, EQ.iActFPNum, MPI_FLOAT, &STATUS);
        	    OFFSETall += EQ.iCmbFPNum*sizeof(float);
        	    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ.ivR_WrtStrtPos[MD.iRANK]*sizeof(int),  EQ.ivL_StabType, EQ.iActFPNum, MPI_INT,   &STATUS);
			    OFFSETall += EQ.iCmbFPNum*sizeof(int);

			    MPI_Barrier( MPI_COMM_WORLD );
		    }
            //-----------------------------------------------------------------------
		    //----------------------------------------------------------------------- 
		    for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
            {  	iGlobPos = i + MD.ivF_START[MD.iRANK]; 
			    //-----------------------------------------
			    if (TR.ivFG_Activated[iGlobPos] == 1)     
           	    {	if (MD.iChgBtwEQs == 0)
					{	TR.fvFL_CurFric[i]  = TR.fvFL_StaFric[i];				}
					else
				    {	//this is copied from LoadData where it is done first
					    fTemp0              = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
					    TR.fvFL_CurDcVal[i] = TR.fvFL_RefDcVal[i]  *(1.0 + TR.fvFL_RefDcVal_vari[i]*fTemp0);
        			    //------------------------
					    fTemp0              = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        			    TR.fvFL_StaFric[i]  = TR.fvFL_RefStaFric[i]*(1.0 + TR.fvFL_RefStaFric_vari[i]*fTemp0); 
					    //------------------------
					    fTemp0              = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
					    fTemp1              = (TR.fvFL_RefStaFric[i] - TR.fvFL_RefDynFric[i])/TR.fvFL_RefStaFric[i]; //reference friction change as fraction of static coefficient
				
					    fTemp2              = fTemp1*(1.0 + TR.fvFL_RefDynFric_vari[i]*fTemp0);
        			    TR.fvFL_DynFric[i]  = TR.fvFL_StaFric[i]*(1.0 - fTemp2);
					    //--------------------------------------------------------------
					    fTemp0 = gsl_matrix_float_get(K.FF_SS, TR.ivFL_SelfLoc[i], i);
					    fTemp1 = gsl_matrix_float_get(K.FF_DD, TR.ivFL_SelfLoc[i], i);
					    fTemp2 = (TR.fvFL_DynFric[i] - TR.fvFL_StaFric[i]) *-1.0*TR.fvFL_RefNrmStrs[i]; //gives a negative shear stress b/c of (dynF - staF)
        			    fTemp3 = fTemp2/((fTemp0 > fTemp1) ? fTemp0 : fTemp1); // b/c "self" has negative sign, the two negatives give a positive slip amount... (for weakening case, when dyn < stat fric) 
        		
        			    if (fTemp3 > TR.fvFL_CurDcVal[i]) 	{ 	TR.ivFL_StabT[i] = 1; 	TR.fvFL_CurFric[i]  = TR.fvFL_StaFric[i];							}
        			    else	
						{   if (fTemp2 < 0.0)  				{ 	TR.ivFL_StabT[i] = 2; 	TR.fvFL_CurFric[i]  = TR.fvFL_StaFric[i];							}
           					else 							{   TR.ivFL_StabT[i] = 3; 																		}      
        	    }	}	}		
			    //-----------------------------------------
			    if (MD.iUsePSeis == 1)
				{	TR.fvFL_PSeis_Step[i] = 0.0;
					TR.fvFL_PSeis_T0_N[i] = gsl_vector_float_get(TR.fvFL_CurStrsN,i);
					TR.fvFL_PSeis_T0_F[i] = TR.fvFL_CurFric[i];
				}
				else
				{  	TR.fvFL_CurFric[i] = TR.fvFL_StaFric[i]; //this happens for all unstable and cond. stable patches anyways; but stable ones could evolve posteismically => if not used, set value to static one...
					gsl_vector_float_set(TR.fvFL_CurStrsN, i, TR.fvFL_RefNrmStrs[i]);	
				}
			    //-----------------------------------------				
		    }
		    for (i = 0; i < MD.ivB_OFFSET[MD.iRANK]; i++) 
        	{	if (MD.iUsePSeis == 1)
				{	TR.fvBL_PSeis_Step[i] = 0.0;
					TR.fvBL_PSeis_T0_N[i] = fabs(gsl_vector_float_get(TR.fvBL_CurStrsN,i));	
					TR.fvBL_PSeis_T0_F[i] = sqrtf(gsl_vector_float_get(TR.fvBL_CurStrsH,i)*gsl_vector_float_get(TR.fvBL_CurStrsH,i)  + gsl_vector_float_get(TR.fvBL_CurStrsV,i)*gsl_vector_float_get(TR.fvBL_CurStrsV,i))  /  TR.fvBL_PSeis_T0_N[i] ;   
					TR.fvBL_CurFric[i]    = TR.fvBL_PSeis_T0_F[i];          
		    }	}
			//-----------------------------------------
			for (i = 0; i < MD.iFPNum; i++) 		{			TR.ivFG_Ptch_t0[i]   = 0;			TR.ivFG_Activated[i] = 0;				}
		    //--------------------------------------------------------------    
        	MPI_Barrier( MPI_COMM_WORLD );  
        }
	}
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
	if (MD.iRANK == 0)		{			MPI_File_write_at(fp_MPIOUT, 0, &MD.iEQcntr,   1, MPI_INT,      &STATUS);					}
 	MPI_File_close(&fp_MPIOUT);
    
	//-------------------------------------------------------------------------------------
    if (MD.iRANK == 0)             
    {   timer = clock() - timer;     
        double time_taken;
        time_taken  = ((double)timer)/CLOCKS_PER_SEC;
        time_taken /= 60.0;
        fprintf(stdout,"Total RunTime in minutes: %6.2f\n",time_taken);
		fprintf(stdout,"Times iSize =>  total of %6.2f CPU hours\n",(time_taken*(float)MD.iSIZE)/60.0);
    }
   	//-------------------------------------------------------------------------------------
	MPI_Finalize();
    return 0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float      TheSign(float TestVal)
{   float Sign;
    Sign = (TestVal < 0.0) ? -1.0 : 1.0;
    return Sign;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetUpdatedFriction(float B4Fric, float RefFric, float CurFric, float DynFric, int FricLaw, int StabType, float CurD_c, float AccumSlip, float PrevSlip, float HealingFact)
{	float NextFric = CurFric,		fTemp;
	//-----------------------------------
 	if (PrevSlip > 0.0)	
	{
		if      (FricLaw == 1)		{				NextFric = DynFric;																																}

		else if (FricLaw == 2)		{				fTemp    = AccumSlip/CurD_c < 1.0 ? AccumSlip/CurD_c : 1.0;				NextFric = DynFric + (1.0 - fTemp)*(RefFric - DynFric);					}

		else if (FricLaw == 3)		{				fTemp    = PrevSlip/CurD_c < 1.0 ? PrevSlip/CurD_c : 1.0;				NextFric = DynFric + (1.0 - fTemp)*(RefFric - DynFric);					}
    }
	else	{		if (StabType != 3)		{		NextFric = CurFric + HealingFact*(B4Fric -CurFric); 	}				}					

	return NextFric;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//