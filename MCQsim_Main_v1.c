#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>

#define PLOT2SCREEN             1
#define MINPATCHNUM4CAT         1
#define USEBOUNDARY4POSTSEIS    1

#define USEVPVSVELOCITY         0
#define WRITESTFOFEACHELEMENT   0
#define WRITEWITHANDWITHOUTPROP 0

#define CUTDISTANCE             1.0//for stiffness matrix -> when to flag values
#define OVERSHOOTFRAC           1.20//dynamic overshoot factor;
#define FRAC2STARTRUPT          1.02//VALUE MUST BE >= 1.0; determines by how much much stress must exceed static strength to initate rupture; value of 1.1 means 10%; 1.05 means 5%
#define LOADINGSTEPINTSEIS      0.05//this is fraction of a day => 0.02 == 1/50'th of a day ~30min
#define LOADSTEPS_POW2          12
#define MAXITERATION4BOUNDARY   200
#define MAXMOMRATEFUNCLENGTH    5000

#define MAX(A,B) ((A) > (B)  ? (A) : (B))
#define MIN(A,B) ((A) < (B)  ? (A) : (B))
#define SIGN(A)  ((A) >(0.0) ? (1.0) : (-1.0))

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
struct MDstruct
{   char    cInputName[512];

    int         iRANK,                      iSIZE;
    int         iRunNum,                    iUsePSeis,              iUseProp,               iSeedStart;
    int         iF_BASEelem,                iF_ADDelem,             iB_BASEelem,            iB_ADDelem;
    int         iFPNum,                     iBPNum,                 iFVNum,                 iBVNum;
    int         iFSegNum,                   iBSegNum,               iChgBtwEQs,             iMaxSTFlgth;
    int         iGlobTTmax,                 iWritePos,              iLoadSteps,             iStepNum,           iEQcntr;
    int         iTimeYears,                 iRecLgth,               iSTFcntr;

    float       fAftrSlipTime,              fDeepRelaxTime,         fRecLgth,               fHealFact;
    float       fFltLegs,                   fBndLegs,               fMinMag4Prop,           fPSeis_Step;
    float       fMedDense,                  fAddNrmStrs,            fVp,                    fVs;
    float       fPoisson,                   fLambda,                fShearMod,              fMeanStiffness,     fNextStep;
    float       fIntSeisDeltT_InYrs,        fVpVsRatio,             fUnitSlipF,             fUnitSlipB,         fDeltT,         fRealDeltT;

    int         *ivF_OFFSET,                *ivF_START,             *ivB_OFFSET,            *ivB_START;
    int         *iv_OFFST2,                 *iv_STRT2;
    int         *iv_OFFST3,                 *iv_STRT3;
};
//------------------------------------------------------------------
struct TRstruct
{   int                 *ivFL_StabT,                *ivFL_Ptch_t0,          *ivFL_Ptch_tDur;
    int                 *ivFL_Activated,            *ivFL_FricLaw;
    int                 *ivFG_SegID_temp,           *ivFG_FltID_temp,       *ivFG_Flagged_temp;
    int                 *ivFG_V1_temp,              *ivFG_V2_temp,          *ivFG_V3_temp;
    float               *fvFG_CentE_temp,           *fvFG_CentN_temp,       *fvFG_CentZ_temp;
    float               *fvFG_StressRatetemp,       *fvFG_SlipRatetemp,     *fvFG_Raketemp;
    float               *fvFL_SlipRate_temp,        *fvFL_SlipRake_temp;

    int                 *ivBG_V1_temp,              *ivBG_V2_temp,          *ivBG_V3_temp,              *ivBG_SegID_temp;
    float               *fvBG_CentE_temp,           *fvBG_CentN_temp,       *fvBG_CentZ_temp;

    gsl_vector_float    *fvFL_RefNrmStrs,           *fvFL_Area;
    gsl_vector_float    *fvFL_SelfStiffStk,         *fvFL_SelfStiffDip,     *fvFL_MeanSelfStiff;
    gsl_vector_float    *fvBL_SelfStiffStk,         *fvBL_SelfStiffDip,     *fvBL_SelfStiffOpn;

    gsl_vector_float    *fvFL_RefStaFric,           *fvFL_RefDynFric,       *fvFL_RefDcVal;
    gsl_vector_float    *fvFL_RefStaFric_vari,      *fvFL_RefDynFric_vari,  *fvFL_RefDcVal_vari;
    gsl_vector_float    *fvFL_StaFric,              *fvFL_DynFric,          *fvFL_ArrFric,              *fvFL_CurFric;
    gsl_vector_float    *fvFL_B4_Fric,              *fvFL_TempRefFric,      *fvFL_CurDcVal;
    gsl_vector_float    *fvFL_StaFricMod_temp,      *fvFL_DynFricMod_temp,  *fvFL_NrmStrsMod_temp,      *fvFL_DcMod_temp;

    gsl_vector_float    *fvFL_PSeis_T0_F,           *fvFL_CurSlpH,          *fvFL_CurSlpV,              *fvFL_AccumSlp;
    gsl_vector_float    *fvBL_PSeis_T0_S,           *fvBL_PSeis_T0_N,       *fvFL_OverShotStress;

    gsl_vector_float    *fvFL_StrsRateStk,          *fvFL_StrsRateDip,      *fvFL_MaxSlip;
    gsl_vector_float    *fvFL_CurStrsH,             *fvFL_CurStrsV,         *fvFL_CurStrsN;
    gsl_vector_float    *fvFL_B4_StrsH,             *fvFL_B4_StrsV,         *fvFL_B4_StrsN;
    gsl_vector_float    *fvBL_B4_StrsH,             *fvBL_B4_StrsV,         *fvBL_B4_StrsN;
    gsl_vector_float    *fvBL_CurStrsH,             *fvBL_CurStrsV,         *fvBL_CurStrsN;

    gsl_matrix_int      *imFGL_TTP,                 *imFGL_TTS;
    gsl_matrix_float    *fmFGL_SrcRcvH,             *fmFGL_SrcRcvV,         *fmFGL_SrcRcvN;
};
//------------------------------------------------------------------
struct VTstruct
{   float       *fvFG_VlX_temp,             *fvFG_VlY_temp,         *fvFG_Hght_temp;
    float       *fvFG_PosE_temp,            *fvFG_PosN_temp,        *fvFG_PosZ_temp;
    float       *fvBG_PosE_temp,            *fvBG_PosN_temp,        *fvBG_PosZ_temp;
};
//------------------------------------------------------------------
struct Kstruct
{   gsl_matrix_float    *FF_SS,         *FF_SD,     *FF_SO; //first index is "receiver"; second is "source"; first is long; second is short
    gsl_matrix_float    *FF_DS,         *FF_DD,     *FF_DO; //this means: FB_DS => Fault is receiver, Boundary the source; Dipslip causing strike_stress

    gsl_matrix_float    *FFp_SS,        *FFp_SD,    *FFp_SO; //first index is "receiver"; second is "source"; first is long; second is short
    gsl_matrix_float    *FFp_DS,        *FFp_DD,    *FFp_DO; //this means: FB_DS => Fault is receiver, Boundary the source; Dipslip causing strike_stress

    gsl_matrix_float    *FFs_SS,        *FFs_SD,    *FFs_SO; //first index is "receiver"; second is "source"; first is long; second is short
    gsl_matrix_float    *FFs_DS,        *FFs_DD,    *FFs_DO; //this means: FB_DS => Fault is receiver, Boundary the source; Dipslip causing strike_stress

    gsl_matrix_float    *FB_SS,         *FB_SD,     *FB_SO;
    gsl_matrix_float    *FB_DS,         *FB_DD,     *FB_DO;

    gsl_matrix_float    *BF_SS,         *BF_SD;
    gsl_matrix_float    *BF_DS,         *BF_DD;
    gsl_matrix_float    *BF_OS,         *BF_OD;

    gsl_matrix_float    *BB_SS,         *BB_SD,     *BB_SO;
    gsl_matrix_float    *BB_DS,         *BB_DD,     *BB_DO;
    gsl_matrix_float    *BB_OS,         *BB_OD,     *BB_OO;
};
//------------------------------------------------------------------
struct EQstruct
{   int     iStillOn,           iEndCntr,           iActFPNum;
    int     iCmbFPNum,          iTotlRuptT,         iMRFLgth,           iCombSTFoffs;

    float   fMaxSlip,           fMaxDTau;

    int     *ivR_WrtStrtPos,    *ivR_StfStrtPos,    *ivL_ActPtchID,     *ivL_t0ofPtch,      *ivL_tDurPtch,      *ivL_StabType;
    float   *fvL_PtchSlpH,      *fvL_PtchSlpV,      *fvL_PtchDTau;

    gsl_vector_float            *fvL_EQslipH,       *fvL_EQslipV,       *fvM_MRFvals;
    gsl_matrix_float            *fmFSL_H,           *fmFSL_V,           *fmFSL_N;
    gsl_matrix_float            *fmSTF_slip,        *fmSTF_strength;
};
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
extern void          StrainHS_Nikkhoo(float Stress[6], float Strain[6], float X, float Y, float Z, float P1[3], float P2[3], float P3[3], float SS, float Ds, float Ts, const float mu, const float lambda);

void              InitializeVariables(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct EQstruct *EQ,  struct Kstruct *K, char **argv);
int                         LoadInput(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct EQstruct *EQ);
void                   Build_K_Matrix(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct  Kstruct *K);
void GetSlipLoadingAndWritePreRunData(struct MDstruct *MD, struct VTstruct *VT, struct TRstruct *TR, struct Kstruct *K,char *cFile2_Out,int HaveBoundarySlip);
MPI_Offset         StartEQcatalogFile(MPI_File fp_MPIOUT,  struct MDstruct *MD, struct VTstruct *VT, struct TRstruct *TR);
void                   FreeSomeMemory(struct TRstruct *TR, struct VTstruct *VT);
void                  PreloadTheFault(struct MDstruct *MD, struct TRstruct *TR, float LoadFactor);
float              GetUpdatedFriction(int StabType, float B4Fric, float RefFric, float CurFric, float ArrFric, float CurD_c, float AccumSlip, float PrevSlip, float HealingFact);
MPI_Offset           WriteCatalogFile(MPI_File fp_MPIOUT, MPI_Offset OffsetAll, struct MDstruct *MD, struct TRstruct *TR, struct EQstruct *EQ, int iPrevEQtime, float Magn);
MPI_Offset           WriteEQ_STFsFile(MPI_File fp_STFOUT, MPI_Offset OffsetAll, struct MDstruct *MD, struct TRstruct *TR, struct EQstruct *EQ, float Magn);
void              ResetFaultParameter(struct MDstruct *MD, struct TRstruct *TR, gsl_rng *fRandN);

int                IntPow(int Base, int Exp);
void          GetVertices(struct VTstruct *VT,struct TRstruct *TR, int iGlobPos, int ForBoundary, float fP1[3], float fP2[3], float fP3[3]);
void            GetVector(float fP1[3], float fP2[3], float fP1P2[3]);
void      NormalizeVector(float fVect[3]);
void         CrossProduct(float fA[3],  float fB[3],  float fvNrm[3]);
float         GetDistance(float fT[3], float fL[3]);
float    GetMinDist_8Pnts(float fP1[3], float fP1a[3],float P1b[3], float P1c[3], float fP2[3], float fP2a[3],float P2b[3], float P2c[3]);
void  GetLocKOS_inKmatrix(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3]);
void         RotateTensor(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir);
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int main(int argc, char **argv)
{   if ( (argc  > 3 ) || (argc < 2) )           
    {   fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -np 4 ./MCQsim RunParaFile.txt -optionally adding more file name here");       }
    //------------------------------------------------------------------
    struct  MDstruct MD;//initializing the structures
    struct  TRstruct TR;
    struct  VTstruct VT;
    struct  Kstruct  K;
    struct  EQstruct EQ;
    //------------------------------------------------------------------
    MPI_Init(&argc, &argv); //initializing MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &MD.iRANK);
    MPI_Comm_size(MPI_COMM_WORLD, &MD.iSIZE);
    MPI_Status  STATUS;
    MPI_Offset  OFFSETall; //this offset is for writing to file, which is done by all ranks
    MPI_Offset  OFFSETstf; //this offset is for writing to file, which is done by all ranks
    MPI_File    fp_MPIOUT; //file pointer for writing the catalog
    MPI_File    fp_STFOUT; //file pointer for writing the catalog
    //------------------------------------------------------------------
    // initialize variables -those that are not part of structures
    int         i,          j,          iTemp0,     iGlobPos,       iActElemNum,    iTpos0,         iTpos1,     iPrevEQtime;
    float       fTemp0,     fTemp1,     fTemp2,     fTemp3,         fTemp4,         fTemp5,         fTemp7,     fTemp8; 
    float       fTestMag;
    char        cFile1_Out[512],        cFile2_Out[512],            cFile3_Out[512],                cAppend[512];
    clock_t     timer;
    //------------------------------------------------------------------------------------
    srand(time(0)); 
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------  
    InitializeVariables(&MD, &TR, &VT, &EQ, &K, argv);//some pre-processing steps such as initializing variables/vectors etc...
    //------------------------------------------------------------------    
    //------------------------------------------------------------------
    gsl_vector_int   *ivFL_Temp0   = gsl_vector_int_calloc(  MD.ivF_OFFSET[MD.iRANK]);      
    gsl_vector_float *fvFL_Temp0   = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);          gsl_vector_float *fvFL_Temp1   = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);
    gsl_vector_float *fvFL_Temp2   = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);          gsl_vector_float *fvFL_Temp3   = gsl_vector_float_calloc(MD.ivF_OFFSET[MD.iRANK]);

    gsl_vector_float *fvBL_Temp0   = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);          gsl_vector_float *fvBL_Temp1   = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);
    gsl_vector_float *fvBL_Temp2   = gsl_vector_float_calloc(MD.ivB_OFFSET[MD.iRANK]);
    
    gsl_vector_int   *ivFG_Temp0   = gsl_vector_int_calloc(  MD.iFPNum);
    gsl_vector_float *fvFG_Temp0   = gsl_vector_float_calloc(MD.iFPNum);                        gsl_vector_float *fvFG_Temp1   = gsl_vector_float_calloc(MD.iFPNum);                
    
    gsl_vector_float *fvBG_Temp0   = gsl_vector_float_calloc(MD.iBPNum);                        gsl_vector_float *fvBG_Temp1   = gsl_vector_float_calloc(MD.iBPNum);    
    gsl_vector_float *fvBG_Temp2   = gsl_vector_float_calloc(MD.iBPNum);                                                            
    //------------------------------------------------------------------
    strcpy(cFile1_Out, MD.cInputName);                  strcat(cFile1_Out,"_");         sprintf(cAppend, "%d",MD.iRunNum);  strcat(cFile1_Out,cAppend);    strcat(cFile1_Out,"_Catalog.dat");
    strcpy(cFile2_Out, MD.cInputName);                  strcat(cFile2_Out,"_");         sprintf(cAppend, "%d",MD.iRunNum);  strcat(cFile2_Out,cAppend);    strcat(cFile2_Out,"_PreRunData.dat");
    strcpy(cFile3_Out, MD.cInputName);                  strcat(cFile3_Out,"_");         sprintf(cAppend, "%d",MD.iRunNum);  strcat(cFile3_Out,cAppend);    strcat(cFile3_Out,"_STFCat.dat");
    //------------------------------------------------------------------  
    gsl_rng *fRandN; // this is pretty much straight from the GSL reference, the default RNG has good performance, so no need to change
    const gsl_rng_type *RandT;
    gsl_rng_env_setup(); 
    RandT   = gsl_rng_default;
    fRandN  = gsl_rng_alloc(RandT);
    unsigned long RSeed =(unsigned long)(MD.iRANK + MD.iSeedStart);   //so, every rank has its own random number sequence ==> that sequence differs from the sequences of the other ranks
    gsl_rng_set(fRandN, RSeed);
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    iTemp0 = LoadInput(&MD, &TR, &VT, &EQ);
    MPI_Allreduce(MPI_IN_PLACE, &iTemp0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);  
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------  
    if (MD.iRANK == 0)              {       fprintf(stdout,"start Kmatrix");            }
    Build_K_Matrix(&MD, &TR, &VT, &K); //later on, if helpful, this could also be stored to file and the loaded here...
    if (MD.iRANK == 0)              {       fprintf(stdout,"finish Kmatrix\n");         }
    //------------------------------------------------------------------------------------  
    //------------------------------------------------------------------------------------      
    GetSlipLoadingAndWritePreRunData(&MD, &VT, &TR, &K, cFile2_Out, iTemp0);    
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //open the EQ catalog file here once and then keep open (if not kept open, the continuous opening and closing slows things down and also might cause code to crash
    MPI_File_delete(cFile1_Out, MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, cFile1_Out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIOUT);
    OFFSETall = StartEQcatalogFile(fp_MPIOUT, &MD, &VT, &TR);
    
    if (WRITESTFOFEACHELEMENT == 1)
    {   MPI_File_delete(cFile3_Out, MPI_INFO_NULL);
        MPI_File_open(MPI_COMM_WORLD, cFile3_Out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_STFOUT);
        OFFSETstf = StartEQcatalogFile(fp_STFOUT, &MD, &VT, &TR);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    FreeSomeMemory(&TR, &VT);
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    PreloadTheFault(&MD, &TR, 0.93);
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    timer         = clock();
    MD.iRecLgth   = (int)(MD.fRecLgth/MD.fIntSeisDeltT_InYrs);
    MD.iTimeYears = 0;
    iPrevEQtime   = 0;
    if (MD.iRANK == 0)               {  fprintf(stdout,"Start earthquake catalog. Time steps: %d \n",MD.iRecLgth);          }
    
    while (MD.iTimeYears <= MD.iRecLgth)
    {   //-------------------------------------------------
        gsl_vector_float_memcpy(TR.fvFL_B4_StrsH, TR.fvFL_CurStrsH);        gsl_vector_float_memcpy(TR.fvFL_B4_StrsV, TR.fvFL_CurStrsV);            
        gsl_vector_float_memcpy(TR.fvFL_B4_Fric,  TR.fvFL_CurFric);
        gsl_vector_float_memcpy(TR.fvBL_B4_StrsH, TR.fvBL_CurStrsH);        gsl_vector_float_memcpy(TR.fvBL_B4_StrsV, TR.fvBL_CurStrsV);  
        gsl_vector_float_memcpy(TR.fvBL_B4_StrsN, TR.fvBL_CurStrsN);  
        //-------------------------------------------------
        MD.fNextStep = (float)MD.iStepNum; //e.g., 1024
        fTemp7       = MD.fNextStep;
        for (j = 0; j < MD.iLoadSteps; j++) //e.g., 10
        {   //------------------------------------------------- 
            gsl_vector_float_memcpy(TR.fvFL_CurStrsH, TR.fvFL_B4_StrsH);    gsl_vector_float_memcpy(TR.fvFL_CurStrsV, TR.fvFL_B4_StrsV);
            gsl_vector_float_memcpy(TR.fvFL_CurFric,  TR.fvFL_B4_Fric);  
            gsl_vector_float_memcpy(TR.fvBL_CurStrsH, TR.fvBL_B4_StrsH);    gsl_vector_float_memcpy(TR.fvBL_CurStrsV, TR.fvBL_B4_StrsV);  
            gsl_vector_float_memcpy(TR.fvBL_CurStrsN, TR.fvBL_B4_StrsN);        
            //-------------------------------------------------
            gsl_vector_float_memcpy(fvFL_Temp0, TR.fvFL_StrsRateStk);       gsl_vector_float_scale(fvFL_Temp0, MD.fNextStep);           gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);   
            gsl_vector_float_memcpy(fvFL_Temp0, TR.fvFL_StrsRateDip);       gsl_vector_float_scale(fvFL_Temp0, MD.fNextStep);           gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0);   
            //----------------------------------------------------------------------------
            if (MD.iUsePSeis == 1) //  https://en.wikipedia.org/wiki/Stress_relaxation    http://web.mit.edu/course/3/3.11/www/modules/visco.pdf  => page 9ff; using Maxwell spring-dashpot model    
            {   //------------------------------------------------------------- 
                fTemp0          = expf(-1.0*(MD.fIntSeisDeltT_InYrs*(MD.fPSeis_Step+MD.fNextStep))/MD.fAftrSlipTime); //factor/fraction from decay function
                for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)           
                {   fTemp1 = gsl_vector_float_get(TR.fvFL_StaFric, i) +fTemp0*gsl_vector_float_get(TR.fvFL_PSeis_T0_F, i);      
                    gsl_vector_float_set(TR.fvFL_CurFric, i, fTemp1);           
                }   
                //-------------------------------------------------------------
                fTemp0          = expf(-1.0*(MD.fIntSeisDeltT_InYrs*(MD.fPSeis_Step+MD.fNextStep))/MD.fDeepRelaxTime); //factor/fraction from decay function
                //---------------------------
                if ( (USEBOUNDARY4POSTSEIS == 1) &&  (MD.iBPNum > 0) )
                {   iTemp0 = 0;
                    for (i = 0; i < MD.ivB_OFFSET[MD.iRANK]; i++)   
                    {   //------------------------------------------
                        fTemp1 = fTemp0*gsl_vector_float_get(TR.fvBL_PSeis_T0_S, i);    
                        fTemp3 = sqrtf(gsl_vector_float_get(TR.fvBL_CurStrsH, i)*gsl_vector_float_get(TR.fvBL_CurStrsH, i) + gsl_vector_float_get(TR.fvBL_CurStrsV, i)*gsl_vector_float_get(TR.fvBL_CurStrsV, i) );
                        fTemp2 = fTemp0*gsl_vector_float_get(TR.fvBL_PSeis_T0_N, i); 
                        fTemp4 = fabs(gsl_vector_float_get(TR.fvBL_CurStrsN, i));
                        //------------------------------------------
                        if ((fTemp3 - fTemp1) <= 0.01)
                        {   gsl_vector_float_set(fvBL_Temp0, i, 0.0);       
                            gsl_vector_float_set(fvBL_Temp1, i, 0.0);                     
                        }
                        else
                        {   fTemp5 = -1.0*((fTemp3-fTemp1)/fTemp3 *gsl_vector_float_get(TR.fvBL_CurStrsH, i)) / gsl_vector_float_get(TR.fvBL_SelfStiffStk, i);              
                            gsl_vector_float_set(fvBL_Temp0, i, fTemp5); 
                            fTemp5 = -1.0*((fTemp3-fTemp1)/fTemp3 *gsl_vector_float_get(TR.fvBL_CurStrsV, i)) / gsl_vector_float_get(TR.fvBL_SelfStiffDip, i);          
                            gsl_vector_float_set(fvBL_Temp1, i, fTemp5); 
                            iTemp0 = 1;
                        }                     
                        //------------------------------------------
                        if ((fTemp4 - fTemp2) <= 0.01)
                        {   gsl_vector_float_set(fvBL_Temp2, i, 0.0);                                                                   
                        }
                        else  
                        {   fTemp5 = -1.0*((fTemp4-fTemp2)/fTemp4 *gsl_vector_float_get(TR.fvBL_CurStrsN, i)) /  gsl_vector_float_get(TR.fvBL_SelfStiffOpn, i);             
                            gsl_vector_float_set(fvBL_Temp2, i, fTemp5); 
                            iTemp0 = 1; 
                    }   }
                    //--------------------------------------------------------------
                    MPI_Allreduce(MPI_IN_PLACE, &iTemp0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    //--------------------------------------------------------------
                    if (iTemp0 == 1)
                    {   MPI_Allgatherv(fvBL_Temp0->data, MD.ivB_OFFSET[MD.iRANK], MPI_FLOAT, fvBG_Temp0->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
                        MPI_Allgatherv(fvBL_Temp1->data, MD.ivB_OFFSET[MD.iRANK], MPI_FLOAT, fvBG_Temp1->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
                        MPI_Allgatherv(fvBL_Temp2->data, MD.ivB_OFFSET[MD.iRANK], MPI_FLOAT, fvBG_Temp2->data, MD.ivB_OFFSET, MD.ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
                    
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_SS->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);                         
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_SD->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);
                        gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0);             
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_SS->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);      
                        gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);                     
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_SD->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp0);                                     
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_SO->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp0);                             
                    
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_DS->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);      
                        gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);                                                     
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_DD->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);                                          
                        gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0); 
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_DS->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);                             
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_DD->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp0);                                             
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_DO->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp0);                                 

                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_OS->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp2->data, 1, 0.0, fvFL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);                                                         
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.BF_OD->data, MD.ivF_OFFSET[MD.iRANK], fvBG_Temp2->data, 1, 0.0, fvFL_Temp0->data, 1);                                          
                        gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0); 
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_OS->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);                                                         
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_OD->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp0);                                     
                        cblas_sgemv(CblasRowMajor,CblasTrans, MD.iBPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.BB_OO->data, MD.ivB_OFFSET[MD.iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);  
                        gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp0);                                 
            }   }   }
            //---------------------------       
            gsl_vector_float_set_zero(fvFL_Temp0);                              
            gsl_vector_float_set_zero(fvFL_Temp1);  
            iTemp0 = 0; 
            for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
            {   if (TR.ivFL_StabT[i] != 1) 
                {   fTemp3   = sqrtf(gsl_vector_float_get(TR.fvFL_CurStrsH, i)*gsl_vector_float_get(TR.fvFL_CurStrsH, i) + gsl_vector_float_get(TR.fvFL_CurStrsV, i)*gsl_vector_float_get(TR.fvFL_CurStrsV, i));
                    fTemp4   = fTemp3 - (gsl_vector_float_get(TR.fvFL_CurFric, i)*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i)); //this is the excess stress (is excess if value > 0) I currently have with respect to curr strength
                    if (fTemp4 > 0.01) //if I have excess shear stress above current strength and I am not looking at an unstable patch
                    {   
                        fTemp5   = (-1.0*((fTemp4/fTemp3)*gsl_vector_float_get(TR.fvFL_CurStrsH, i)) )/ gsl_vector_float_get(TR.fvFL_SelfStiffStk, i); // slip amount to release excess horizontal shear stress
                        gsl_vector_float_set(fvFL_Temp0, i, fTemp5);  //the strike slip component    
                        fTemp5   = (-1.0*((fTemp4/fTemp3)*gsl_vector_float_get(TR.fvFL_CurStrsV, i)) )/ gsl_vector_float_get(TR.fvFL_SelfStiffDip, i); // slip amount to release excess vertical shear stress        
                        gsl_vector_float_set(fvFL_Temp1, i, fTemp5);  //the dip slip component  
                        iTemp0   = 1;
            }   }   }      
            MPI_Allreduce(MPI_IN_PLACE, &iTemp0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            //--------------------------------------------------------------
            if (iTemp0 == 1)
            {   MPI_Allgatherv(fvFL_Temp0->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Allgatherv(fvFL_Temp1->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);         
                
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_SS->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);      
                gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);                                                     
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_SD->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);                      
                gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0); 

                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_DS->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);  
                gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp0);                                                         
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_DD->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);                      
                gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp0);                                  
            }   
            //--------------------------------------------------------------
            EQ.iStillOn = 0;
            for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)
            {   if (TR.ivFL_StabT[i] == 1 )
                {   fTemp3   = sqrtf(gsl_vector_float_get(TR.fvFL_CurStrsH, i)*gsl_vector_float_get(TR.fvFL_CurStrsH, i) + gsl_vector_float_get(TR.fvFL_CurStrsV, i)*gsl_vector_float_get(TR.fvFL_CurStrsV, i)); //currently applied stress
                    fTemp4   = gsl_vector_float_get(TR.fvFL_CurFric, i)*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i); //current strength

                if (fTemp3 >= fTemp4*FRAC2STARTRUPT)    {           EQ.iStillOn  = 1;                                   }   
            }   }     
            MPI_Allreduce(MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);       
            //----------------------------------------------------------------------------
            fTemp7 /= 2.0;
            if (EQ.iStillOn == 0)
            {   if      (j == 0)                        {   break;                                            }
                else if (j  < (MD.iLoadSteps-1))        {   MD.fNextStep = MD.fNextStep + fTemp7;             }
            }
            else    
            {   if      (j  < (MD.iLoadSteps-2))        {   MD.fNextStep = MD.fNextStep - fTemp7;             }       
        }   }
        MD.fPSeis_Step +=      MD.fNextStep;
        MD.iTimeYears  += (int)MD.fNextStep; 
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        if (EQ.iStillOn == 1)
        {   //-----------------------------------------
            gsl_vector_float_memcpy(TR.fvFL_B4_StrsH, TR.fvFL_CurStrsH);        gsl_vector_float_memcpy(TR.fvFL_B4_StrsV, TR.fvFL_CurStrsV);            
            gsl_vector_float_memcpy(TR.fvFL_B4_StrsN, TR.fvFL_CurStrsN);        gsl_vector_float_memcpy(TR.fvFL_B4_Fric,  TR.fvFL_CurFric);
            
            gsl_vector_float_set_zero(TR.fvFL_AccumSlp);                        gsl_vector_float_set_zero(EQ.fvM_MRFvals);                          
            gsl_vector_float_set_zero(EQ.fvL_EQslipH);                          gsl_vector_float_set_zero(EQ.fvL_EQslipV);                          
            gsl_matrix_float_set_zero(EQ.fmSTF_slip);                           gsl_matrix_float_set_zero(EQ.fmSTF_strength);
            gsl_vector_float_set_zero(fvFL_Temp0);                              gsl_vector_float_set_zero(fvFL_Temp1);  


            for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)               
            {   fTemp0   = gsl_vector_float_get(TR.fvFL_DynFric, i)*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i)  - 1.0*gsl_vector_float_get(TR.fvFL_MeanSelfStiff,i) *gsl_vector_float_get(TR.fvFL_CurDcVal,i);
                fTemp1   = gsl_vector_float_get(TR.fvFL_StaFric, i)*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i);
                fTemp0   = MAX(fTemp0, fTemp1) / (-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i));
                gsl_vector_float_set(TR.fvFL_CurFric, i,     fTemp0);
                gsl_vector_float_set(TR.fvFL_TempRefFric, i, fTemp0);  
                
                fTemp1   = fabs(gsl_vector_float_get(TR.fvFL_CurFric, i) - gsl_vector_float_get(TR.fvFL_DynFric, i))*-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i);
                
                fTemp5   = MD.fVs*((fTemp1*1.0e+6)/MD.fShearMod);
                fTemp5   = (fTemp5 > 0.01) ? fTemp5 : 0.01; //the 0.01 is the minimum slip "velocity" (1cm/s)
                fTemp5   = fTemp5*MD.fRealDeltT;
                gsl_vector_float_set(TR.fvFL_MaxSlip, i, fTemp5); 
            }


            EQ.iActFPNum  =  0;               
            EQ.iMRFLgth   = -1;                 
            EQ.iTotlRuptT = -1;     
            //------------------------------------------------------------------------------------------
            //------------------------------------------------------------------------------------------  
            // EARTHQUAKE ITERATION LOOP STARTS
            while (EQ.iStillOn == 1)
            {   EQ.iTotlRuptT++;    //continuously count time since initiation, set "ongoing" to FALSE => only if more slip on patches is added in the coming iteration, it gets to be reset            
                EQ.iStillOn = 0;  
                //----------------------------------------------------
                for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
                {   
                    fTemp0   = sqrtf(gsl_vector_float_get(TR.fvFL_CurStrsH, i)*gsl_vector_float_get(TR.fvFL_CurStrsH, i) + gsl_vector_float_get(TR.fvFL_CurStrsV, i)*gsl_vector_float_get(TR.fvFL_CurStrsV, i));
                    fTemp7   = (gsl_vector_float_get(TR.fvFL_CurStrsN, i) < 0.0 ) ? gsl_vector_float_get(TR.fvFL_CurStrsN, i) : 0.0 ;                           
                    gsl_vector_float_set(TR.fvFL_CurStrsN, i, fTemp7);
                    //------------------------------------------------
                    if (TR.ivFL_Activated[i] == 1)
                    {   fTemp7 = sqrtf(gsl_vector_float_get(fvFL_Temp0, i) *gsl_vector_float_get(fvFL_Temp0, i)  +  gsl_vector_float_get(fvFL_Temp1, i) *gsl_vector_float_get(fvFL_Temp1, i));
                        fTemp7 = GetUpdatedFriction(TR.ivFL_StabT[i] , gsl_vector_float_get(TR.fvFL_B4_Fric,i), gsl_vector_float_get(TR.fvFL_TempRefFric,i), gsl_vector_float_get(TR.fvFL_CurFric,i), gsl_vector_float_get(TR.fvFL_ArrFric,i), gsl_vector_float_get(TR.fvFL_CurDcVal,i), gsl_vector_float_get(TR.fvFL_AccumSlp,i), fTemp7, MD.fHealFact);
                        gsl_vector_float_set(TR.fvFL_CurFric, i, fTemp7);  
                        if (WRITESTFOFEACHELEMENT == 1)                                                                        
                        {   gsl_matrix_float_set(EQ.fmSTF_strength, i, EQ.iTotlRuptT, gsl_vector_float_get(TR.fvFL_CurFric, i));    
                    }   }
                    //------------------------------------------------
                    fTemp7 = (gsl_vector_float_get(TR.fvFL_CurFric, i) *-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i));
                    fTemp8 = (gsl_vector_float_get(TR.fvFL_ArrFric, i) *-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i));
                    //------------------------------------------------
                    if ((TR.ivFL_Activated[i] == 0) && (fTemp0 >= fTemp7*FRAC2STARTRUPT)) 
                    {   TR.ivFL_Activated[i]      = 1;  
                        TR.ivFL_Ptch_t0[i]        = EQ.iTotlRuptT; 
                        EQ.iActFPNum++;                      
                        if (WRITESTFOFEACHELEMENT == 1)                                                                        
                        {   gsl_matrix_float_set(EQ.fmSTF_strength, i, EQ.iTotlRuptT, gsl_vector_float_get(TR.fvFL_CurFric, i)); 
                    }   }             
                    //------------------------------------------------
                    fTemp1 = (fTemp0 - fTemp7); //amount of excess stress above current friction level
                    fTemp2 = (fTemp0 - fTemp8); //amount of excess stress above arrest friction level                    
                    //------------------------------------------------
                    if ((TR.ivFL_Activated[i] == 1) && (fTemp1 > 0.0) && (fTemp2 >= gsl_vector_float_get(TR.fvFL_OverShotStress, i)))                     
                    {   EQ.iStillOn          = 1;                       
                        EQ.iMRFLgth          = EQ.iTotlRuptT;       
                        TR.ivFL_Ptch_tDur[i] = EQ.iTotlRuptT - TR.ivFL_Ptch_t0[i] +1;   

                        fTemp5      = gsl_vector_float_get(TR.fvFL_MaxSlip, i);

                        fTemp2      = -1.0*(fTemp1/fTemp0 *gsl_vector_float_get(TR.fvFL_CurStrsH, i)) / gsl_vector_float_get(TR.fvFL_SelfStiffStk, i); // slip amount to release excess horizontal shear stress
                        fTemp2      = fTemp5 < fabs(fTemp2) ?   SIGN(fTemp2)*fTemp5 : fTemp2;
                        gsl_vector_float_set(fvFL_Temp0, i, fTemp2); 
                                        
                        fTemp3      = -1.0*(fTemp1/fTemp0 *gsl_vector_float_get(TR.fvFL_CurStrsV, i)) / gsl_vector_float_get(TR.fvFL_SelfStiffDip, i); // slip amount to release excess vertical shear stress 
                        fTemp3      = fTemp5  < fabs(fTemp3) ?  SIGN(fTemp3)*fTemp5 : fTemp3;
                        gsl_vector_float_set(fvFL_Temp1, i, fTemp3);    
                                                                    
                        fTemp4      = sqrtf(fTemp2*fTemp2 +fTemp3*fTemp3);      
                        if (WRITESTFOFEACHELEMENT == 1)                                                                        
                        {   gsl_matrix_float_set(EQ.fmSTF_slip, i, EQ.iTotlRuptT, fTemp4);          }                    

                        if (gsl_vector_float_get(TR.fvFL_AccumSlp,i) == 0.0)                {   gsl_vector_float_set(TR.fvFL_TempRefFric, i, gsl_vector_float_get(TR.fvFL_CurFric, i));                 }
                        fTemp5      = gsl_vector_float_get(TR.fvFL_AccumSlp, i) +fTemp4;
                        gsl_vector_float_set(TR.fvFL_AccumSlp, i, fTemp5);
                        fTemp5      = gsl_vector_float_get(EQ.fvM_MRFvals, EQ.iMRFLgth) + fTemp4 *gsl_vector_float_get(TR.fvFL_Area,i) *MD.fShearMod;
                        gsl_vector_float_set(EQ.fvM_MRFvals, EQ.iMRFLgth, fTemp5);
                    }
                    else //this is necessary here b/c I can get into that if-statement b/c I got activated but might not have any amount of slip left => set value to zero...
                    {   if (TR.ivFL_StabT[i] != 3)                  {   gsl_vector_float_set(TR.fvFL_AccumSlp, i, 0.0);         }
                        gsl_vector_float_set(fvFL_Temp0, i, 0.0);
                        gsl_vector_float_set(fvFL_Temp1, i, 0.0);  
                }   }
                if (EQ.iStillOn == 1)           
                {   gsl_vector_float_add(EQ.fvL_EQslipH, fvFL_Temp0);           
                    gsl_vector_float_add(EQ.fvL_EQslipV, fvFL_Temp1);                           
                }
                //-----------------------------------------------------------------------
                MPI_Allreduce( MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                
                if (EQ.iStillOn == 1)
                {   
                    MPI_Allgatherv(fvFL_Temp0->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD); 
                    MPI_Allgatherv(fvFL_Temp1->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD); 
                    
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_SS->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp2->data, 1);
                    gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp2);                             
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_SD->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp2->data, 1);
                    gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp2);                     
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_SO->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp2->data, 1);      
                    gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
                                
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_DS->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp2->data, 1);  
                    gsl_vector_float_add(TR.fvFL_CurStrsH, fvFL_Temp2);                                                     
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_DD->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp2->data, 1);
                    gsl_vector_float_add(TR.fvFL_CurStrsV, fvFL_Temp2);                                                 
                    cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivF_OFFSET[MD.iRANK], 1.0, K.FF_DO->data, MD.ivF_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp2->data, 1);      
                    gsl_vector_float_add(TR.fvFL_CurStrsN, fvFL_Temp2);
                }
                if (EQ.iTotlRuptT >= MAXMOMRATEFUNCLENGTH-1)        {       EQ.iStillOn     = 0;         }          //a hard stop   
            }
            //----------------------------------------------------------------------------
            //----------------------------------------------------------------------------
            //----------------------------------------------------------------------------
            fTemp0       = 0.0; //this is temp/test magnitude (i.e., combined seismic potential)
            for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
            {   if (TR.ivFL_Activated[i] == 1)     
                {   fTemp3  = sqrtf(gsl_vector_float_get(EQ.fvL_EQslipH, i)*gsl_vector_float_get(EQ.fvL_EQslipH, i) + gsl_vector_float_get(EQ.fvL_EQslipV, i)*gsl_vector_float_get(EQ.fvL_EQslipV, i));
                    fTemp0 += fTemp3*gsl_vector_float_get(TR.fvFL_Area, i); //this is seis. potential of that patch (added to other patches from same ranke)
            }   } 
            //----------------------------------------------------------------------------
            MPI_Allreduce(MPI_IN_PLACE,     &fTemp0, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //combined seismic potential that is released in this event
            fTestMag = (log10f(fTemp0*MD.fShearMod)-9.1)/1.5; 
            //----------------------------------------------------------------------------      
            EQ.iCmbFPNum = 0; 
            MPI_Allreduce(&EQ.iActFPNum,    &EQ.iCmbFPNum,    1      , MPI_INT,   MPI_SUM, MPI_COMM_WORLD);
        
            if (((WRITEWITHANDWITHOUTPROP == 1) && (EQ.iCmbFPNum >= MINPATCHNUM4CAT))||((MD.iUseProp == 0) && (EQ.iCmbFPNum >= MINPATCHNUM4CAT)) || ((MD.iUseProp == 1) && (EQ.iCmbFPNum >= MINPATCHNUM4CAT) && (fTestMag < MD.fMinMag4Prop)))
            {   OFFSETall = WriteCatalogFile(fp_MPIOUT, OFFSETall, &MD, &TR, &EQ, iPrevEQtime, fTestMag);   
                
                if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
                {   OFFSETstf = WriteEQ_STFsFile(fp_STFOUT, OFFSETstf, &MD, &TR, &EQ, fTestMag);
            }   }  
            //----------------------------------------------------------------------------
            //----------------------------------------------------------------------------
            //----------------------------------------------------------------------------
            if (  (MD.iUseProp == 1) && (fTestMag >= MD.fMinMag4Prop) )
            {   gsl_vector_float_memcpy(TR.fvFL_CurStrsH, TR.fvFL_B4_StrsH); //reset stress to pre-EQ kind...and also all the other things to zero i.e., pre-rupture conditions.
                gsl_vector_float_memcpy(TR.fvFL_CurStrsV, TR.fvFL_B4_StrsV);
                gsl_vector_float_memcpy(TR.fvFL_CurStrsN, TR.fvFL_B4_StrsN);        
                //------------------------------------      
                gsl_matrix_float_set_zero(EQ.fmFSL_H);              gsl_matrix_float_set_zero(EQ.fmFSL_V);              gsl_matrix_float_set_zero(EQ.fmFSL_N);
                gsl_matrix_float_set_zero(EQ.fmSTF_slip);           gsl_matrix_float_set_zero(EQ.fmSTF_strength);
                gsl_vector_float_set_zero(EQ.fvL_EQslipH);          gsl_vector_float_set_zero(EQ.fvL_EQslipV);          
                gsl_vector_float_set_zero(TR.fvFL_AccumSlp);        gsl_vector_float_set_zero(EQ.fvM_MRFvals);  
                gsl_vector_float_set_zero(fvFL_Temp2);              gsl_vector_float_set_zero(fvFL_Temp3);      
                //------------------------------------                    
                for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++)   {   TR.ivFL_Ptch_t0[i]   = 0;       TR.ivFL_Ptch_tDur[i]   = 0;     TR.ivFL_Activated[i] = 0;           }           
                EQ.iStillOn = 1;            EQ.iActFPNum = 0;       EQ.iMRFLgth = -1;               EQ.iTotlRuptT =-1;              EQ.iEndCntr   = 0;  
                //------------------------------------
                while (EQ.iStillOn == 1)
                {   EQ.iTotlRuptT++;           
                    EQ.iStillOn = 0; 
                    iActElemNum = 0;    
                    iTpos0      = (EQ.iTotlRuptT+MD.iMaxSTFlgth-1) % MD.iMaxSTFlgth;    //adding MD.iMaxSTFlgth is necessary to ensure that I can loop back and grab last value, even when look position would be negative...                   

                    for (i = 0; i < MD.iSIZE; i++)              {    MD.iv_STRT2[i]  = 0;           MD.iv_OFFST2[i] = 0;    } 
                    //------------------------------------
                    for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
                    {   //------------------------------------------------
                        iGlobPos = i + MD.ivF_START[MD.iRANK];  
                        fTemp0   = gsl_vector_float_get(TR.fvFL_CurStrsH, i) + gsl_matrix_float_get(EQ.fmFSL_H, iTpos0, i);
                        fTemp1   = gsl_vector_float_get(TR.fvFL_CurStrsV, i) + gsl_matrix_float_get(EQ.fmFSL_V, iTpos0, i);
                        fTemp2   = gsl_vector_float_get(TR.fvFL_CurStrsN, i) + gsl_matrix_float_get(EQ.fmFSL_N, iTpos0, i);
                        gsl_vector_float_set(TR.fvFL_CurStrsH, i, fTemp0);
                        gsl_vector_float_set(TR.fvFL_CurStrsV, i, fTemp1);
                        gsl_vector_float_set(TR.fvFL_CurStrsN, i, fTemp2);
                        gsl_matrix_float_set(EQ.fmFSL_H, iTpos0, i, 0.0);
                        gsl_matrix_float_set(EQ.fmFSL_V, iTpos0, i, 0.0);
                        gsl_matrix_float_set(EQ.fmFSL_N, iTpos0, i, 0.0);
                        //------------------------------------------------
                        fTemp0   = sqrtf(gsl_vector_float_get(TR.fvFL_CurStrsH, i)*gsl_vector_float_get(TR.fvFL_CurStrsH, i) + gsl_vector_float_get(TR.fvFL_CurStrsV, i)*gsl_vector_float_get(TR.fvFL_CurStrsV, i));
                        fTemp7   = (gsl_vector_float_get(TR.fvFL_CurStrsN, i) < 0.0 ) ? gsl_vector_float_get(TR.fvFL_CurStrsN, i) : 0.0 ;                           
                        gsl_vector_float_set(TR.fvFL_CurStrsN, i, fTemp7);
                        //------------------------------------------------
                        if (TR.ivFL_Activated[i] == 1)
                        {   fTemp7 = sqrtf(gsl_vector_float_get(fvFL_Temp2, i) *gsl_vector_float_get(fvFL_Temp2, i)  +  gsl_vector_float_get(fvFL_Temp3, i) *gsl_vector_float_get(fvFL_Temp3, i));
                            fTemp7 = GetUpdatedFriction(TR.ivFL_StabT[i] , gsl_vector_float_get(TR.fvFL_B4_Fric,i), gsl_vector_float_get(TR.fvFL_TempRefFric,i), gsl_vector_float_get(TR.fvFL_CurFric,i), gsl_vector_float_get(TR.fvFL_ArrFric,i), gsl_vector_float_get(TR.fvFL_CurDcVal,i), gsl_vector_float_get(TR.fvFL_AccumSlp,i), fTemp7, MD.fHealFact);
                            gsl_vector_float_set(TR.fvFL_CurFric, i, fTemp7);
                            if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
                            {   gsl_matrix_float_set(EQ.fmSTF_strength, i, EQ.iTotlRuptT, gsl_vector_float_get(TR.fvFL_CurFric, i));
                        }   }
                        //------------------------------------------------
                        fTemp7 = (gsl_vector_float_get(TR.fvFL_CurFric, i) *-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i));
                        fTemp8 = (gsl_vector_float_get(TR.fvFL_ArrFric, i) *-1.0*gsl_vector_float_get(TR.fvFL_CurStrsN, i));
                        //------------------------------------------------
                        if ((TR.ivFL_Activated[i] == 0) && (fTemp0 >= fTemp7*FRAC2STARTRUPT)) 
                        {   TR.ivFL_Activated[i]      = 1;            
                            TR.ivFL_Ptch_t0[i]        = EQ.iTotlRuptT; 
                            EQ.iActFPNum++;  
                            if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
                            {   gsl_matrix_float_set(EQ.fmSTF_strength, i, EQ.iTotlRuptT, gsl_vector_float_get(TR.fvFL_CurFric, i));                        
                        }   }
                        fTemp1 = (fTemp0 - fTemp7); //amount of excess stress
                        fTemp2 = (fTemp0 - fTemp8);//amount of excess stress above dynamic friction level   
                        //------------------------------------------------
                        if ((TR.ivFL_Activated[i] == 1) && (fTemp1 > 0.0) && (fTemp2 >= gsl_vector_float_get(TR.fvFL_OverShotStress, i)))       
                        {   EQ.iStillOn          = 1;                       
                            EQ.iMRFLgth          = EQ.iTotlRuptT;  
                            TR.ivFL_Ptch_tDur[i] = EQ.iTotlRuptT - TR.ivFL_Ptch_t0[i] +1;     
                            
                            fTemp5      = gsl_vector_float_get(TR.fvFL_MaxSlip, i);
                                            
                            gsl_vector_int_set(ivFL_Temp0, MD.iv_OFFST2[MD.iRANK], iGlobPos);   
                            //------------------------------------------------  
                            fTemp2      = -1.0*(fTemp1/fTemp0 *gsl_vector_float_get(TR.fvFL_CurStrsH, i)) / gsl_vector_float_get(TR.fvFL_SelfStiffStk, i); // slip amount to release excess horizontal shear stress
                            fTemp2      = fTemp5 < fabs(fTemp2) ?   SIGN(fTemp2)*fTemp5 : fTemp2;
                            gsl_vector_float_set(fvFL_Temp2, i, fTemp2); 
                            gsl_vector_float_set(fvFL_Temp0, MD.iv_OFFST2[MD.iRANK], fTemp2);      
                            
                            fTemp3      = -1.0*(fTemp1/fTemp0 *gsl_vector_float_get(TR.fvFL_CurStrsV, i)) / gsl_vector_float_get(TR.fvFL_SelfStiffDip, i); // slip amount to release excess vertical shear stress 
                            fTemp3      = fTemp5 < fabs(fTemp3) ?   SIGN(fTemp3)*fTemp5 : fTemp3; 
                            gsl_vector_float_set(fvFL_Temp3, i, fTemp3);
                            gsl_vector_float_set(fvFL_Temp1, MD.iv_OFFST2[MD.iRANK], fTemp3);
                            //------------------------------------------------
                            fTemp4      = sqrtf(fTemp2*fTemp2 +fTemp3*fTemp3);
                            if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
                            {   gsl_matrix_float_set(EQ.fmSTF_slip, i, EQ.iTotlRuptT, fTemp4);          	}
                            
                            if (gsl_vector_float_get(TR.fvFL_AccumSlp,i) == 0.0)                {   gsl_vector_float_set(TR.fvFL_TempRefFric, i, gsl_vector_float_get(TR.fvFL_CurFric, i));                 }
                            fTemp5      = gsl_vector_float_get(TR.fvFL_AccumSlp, i) +fTemp4;
                            gsl_vector_float_set(TR.fvFL_AccumSlp, i, fTemp5);
                            fTemp5      = gsl_vector_float_get(EQ.fvM_MRFvals, EQ.iMRFLgth) + fTemp4 *gsl_vector_float_get(TR.fvFL_Area,i) *MD.fShearMod;
                            gsl_vector_float_set(EQ.fvM_MRFvals, EQ.iMRFLgth, fTemp5);  
                            //------------------------------------------------      
                            MD.iv_OFFST2[MD.iRANK] = MD.iv_OFFST2[MD.iRANK]+1;
                        }   
                        else //this is necessary here b/c I can get into that if-statement b/c I got activated but might not have any amount of slip left => set value to zero...
                        {   if (TR.ivFL_StabT[i] != 3)                  {   gsl_vector_float_set(TR.fvFL_AccumSlp, i, 0.0);         }
                            gsl_vector_float_set(fvFL_Temp2, i, 0.0);
                            gsl_vector_float_set(fvFL_Temp3, i, 0.0);
                        }
                    }
                    //--------------------------------------------------------------------
                    MPI_Allreduce(MPI_IN_PLACE, MD.iv_OFFST2, MD.iSIZE, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
                    //---------------------------------------
                    for (i = 1; i < MD.iSIZE; i++)      {   MD.iv_STRT2[i]   = MD.iv_STRT2[i-1] + MD.iv_OFFST2[i-1];    } 
                    iActElemNum = MD.iv_STRT2[MD.iSIZE-1] +MD.iv_OFFST2[MD.iSIZE-1];
                                        
                    MPI_Allgatherv(ivFL_Temp0->data, MD.iv_OFFST2[MD.iRANK], MPI_INT,   ivFG_Temp0->data, MD.iv_OFFST2, MD.iv_STRT2, MPI_INT,   MPI_COMM_WORLD);
                    MPI_Allgatherv(fvFL_Temp0->data, MD.iv_OFFST2[MD.iRANK], MPI_FLOAT, fvFG_Temp0->data, MD.iv_OFFST2, MD.iv_STRT2, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Allgatherv(fvFL_Temp1->data, MD.iv_OFFST2[MD.iRANK], MPI_FLOAT, fvFG_Temp1->data, MD.iv_OFFST2, MD.iv_STRT2, MPI_FLOAT, MPI_COMM_WORLD);
                    //--------------------------------------------------------------------
                    for (j = 0; j < iActElemNum; j++) //the sources of current iteration step
                    {   iTemp0 = gsl_vector_int_get(ivFG_Temp0, j);
                        for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) //the receivers
                        {   iGlobPos = i + MD.ivF_START[MD.iRANK];
                            if (iTemp0 != iGlobPos)
                            {   iTpos0 = (EQ.iTotlRuptT + gsl_matrix_int_get(TR.imFGL_TTP, iTemp0, i)) % MD.iMaxSTFlgth;
                                iTpos1 = (EQ.iTotlRuptT + gsl_matrix_int_get(TR.imFGL_TTS, iTemp0, i)) % MD.iMaxSTFlgth;
                                //http://sites.science.oregonstate.edu/math/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/dotprod/dotprod.html
                                //https://en.wikipedia.org/wiki/Vector_projection
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFp_SS, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFp_DS, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_H, iTpos0, i);          
                                gsl_matrix_float_set(EQ.fmFSL_H, iTpos0, i, fTemp0);
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFp_SD, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFp_DD, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_V, iTpos0, i);          
                                gsl_matrix_float_set(EQ.fmFSL_V, iTpos0, i, fTemp0);
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFp_SO, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFp_DO, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_N, iTpos0, i);          
                                gsl_matrix_float_set(EQ.fmFSL_N, iTpos0, i, fTemp0);
                                //----------------------------------------------------
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFs_SS, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFs_DS, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_H, iTpos1, i);          
                                gsl_matrix_float_set(EQ.fmFSL_H, iTpos1, i, fTemp0);
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFs_SD, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFs_DD, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_V, iTpos1, i);          
                                gsl_matrix_float_set(EQ.fmFSL_V, iTpos1, i, fTemp0);
                                fTemp0  = gsl_vector_float_get(fvFG_Temp0, j)* gsl_matrix_float_get(K.FFs_SO, iTemp0, i) + gsl_vector_float_get(fvFG_Temp1, j)* gsl_matrix_float_get(K.FFs_DO, iTemp0, i);
                                fTemp0 += gsl_matrix_float_get(EQ.fmFSL_N, iTpos1, i);          
                                gsl_matrix_float_set(EQ.fmFSL_N, iTpos1, i, fTemp0);
                                //----------------------------------------------------      
                            }
                            else
                            {   iTpos0  = EQ.iTotlRuptT % MD.iMaxSTFlgth; //travel time to "self" is of course zero; hence, use the current iteration step to write the stress change
                                fTemp4  = gsl_vector_float_get(fvFG_Temp0, j)*gsl_vector_float_get(TR.fvFL_SelfStiffStk, i);
                                fTemp4 += gsl_matrix_float_get(EQ.fmFSL_H, iTpos0, i);          
                                gsl_matrix_float_set(EQ.fmFSL_H, iTpos0, i, fTemp4);
                                
                                fTemp4  = gsl_vector_float_get(fvFG_Temp1, j)*gsl_vector_float_get(TR.fvFL_SelfStiffDip, i);
                                fTemp4 += gsl_matrix_float_get(EQ.fmFSL_V, iTpos0, i);          
                                gsl_matrix_float_set(EQ.fmFSL_V, iTpos0, i, fTemp4);
                    }   }   }   
                    //--------------------------------------------------------------------
                    MPI_Allreduce( MPI_IN_PLACE, &EQ.iStillOn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    if (EQ.iStillOn == 1)           
                    {   EQ.iEndCntr = 0;        
                        gsl_vector_float_add(EQ.fvL_EQslipH, fvFL_Temp2);           
                        gsl_vector_float_add(EQ.fvL_EQslipV, fvFL_Temp3);               
                    }
                    else
                    {   EQ.iEndCntr++;
                        if (EQ.iEndCntr < MD.iMaxSTFlgth)               {       EQ.iStillOn     = 1;         }      
                    }           
                    if (EQ.iTotlRuptT >= MAXMOMRATEFUNCLENGTH-1)        {       EQ.iStillOn     = 0;         }          //a hard stop
                    //--------------------------------------------------------------------
                }   
                //------------------------------------------------------------------------
                //------------------------------------------------------------------------
                //------------------------------------------------------------------------              
                fTemp0       = 0.0; //this is temp/test magnitude (i.e., combined seismic potential)
                for (i = 0; i < MD.ivF_OFFSET[MD.iRANK]; i++) 
                {   if (TR.ivFL_Activated[i] == 1)     
                    {   fTemp3  = sqrtf(gsl_vector_float_get(EQ.fvL_EQslipH, i)*gsl_vector_float_get(EQ.fvL_EQslipH, i) + gsl_vector_float_get(EQ.fvL_EQslipV, i)*gsl_vector_float_get(EQ.fvL_EQslipV, i));
                        fTemp0 += fTemp3*gsl_vector_float_get(TR.fvFL_Area, i); //this is seis. potential of that patch (added to other patches from same ranke)
                }    } 
                //----------------------------------------------------------------------------
                MPI_Allreduce(MPI_IN_PLACE,     &fTemp0, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //combined seismic potential that is released in this event
                fTestMag = (log10f(fTemp0*MD.fShearMod)-9.1)/1.5; 
        
                EQ.iCmbFPNum = 0; 
                MPI_Allreduce(&EQ.iActFPNum,    &EQ.iCmbFPNum,    1      , MPI_INT,   MPI_SUM, MPI_COMM_WORLD);
                //----------------------------------------------------------------------------      
                if (EQ.iCmbFPNum >= MINPATCHNUM4CAT)
                {   OFFSETall = WriteCatalogFile(fp_MPIOUT, OFFSETall, &MD, &TR, &EQ, iPrevEQtime, fTestMag);   
                    if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
                    {   OFFSETstf = WriteEQ_STFsFile(fp_STFOUT, OFFSETstf, &MD, &TR, &EQ, fTestMag); 
            }   }   }
            //-----------------------------------------------------------------------------------------------------------------------
            //-----------------------------------------------------------------------------------------------------------------------
            //-----------------------------------------------------------------------------------------------------------------------
            //-----------------------------------------------------------------------------------------------------------------------
            if (EQ.iCmbFPNum >= MINPATCHNUM4CAT)        {       iPrevEQtime = MD.iTimeYears;        }
            //-----------------------------------------
            ResetFaultParameter(&MD, &TR, fRandN);
            //-----------------------------------------         
            MD.fPSeis_Step = 0.0;
            //--------------------------------------------------------------    
            if ((USEBOUNDARY4POSTSEIS == 1) && (MD.iBPNum > 0) && (MD.iUsePSeis == 1))
            {   
                MPI_Allgatherv(EQ.fvL_EQslipH->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp0->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Allgatherv(EQ.fvL_EQslipV->data, MD.ivF_OFFSET[MD.iRANK], MPI_FLOAT, fvFG_Temp1->data, MD.ivF_OFFSET, MD.ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
            
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_SS->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);  
                gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);                                     
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_SD->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);  
                gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp0);                             
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_SO->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);      
                gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp0);
                
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_DS->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);  
                gsl_vector_float_add(TR.fvBL_CurStrsH, fvBL_Temp0);                                     
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_DD->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);  
                gsl_vector_float_add(TR.fvBL_CurStrsV, fvBL_Temp0);                             
                cblas_sgemv(CblasRowMajor,CblasTrans, MD.iFPNum, MD.ivB_OFFSET[MD.iRANK], 1.0, K.FB_DO->data, MD.ivB_OFFSET[MD.iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);      
                gsl_vector_float_add(TR.fvBL_CurStrsN, fvBL_Temp0);
                
                for (i = 0; i < MD.ivB_OFFSET[MD.iRANK]; i++) 
                {   fTemp1 = sqrtf(gsl_vector_float_get(TR.fvBL_CurStrsH, i)*gsl_vector_float_get(TR.fvBL_CurStrsH, i) + gsl_vector_float_get(TR.fvBL_CurStrsV, i)*gsl_vector_float_get(TR.fvBL_CurStrsV, i) );
                    gsl_vector_float_set(TR.fvBL_PSeis_T0_S, i, fTemp1); 
                    gsl_vector_float_set(TR.fvBL_PSeis_T0_N, i, fabs(gsl_vector_float_get(TR.fvBL_CurStrsN, i)));
            }   }
            //--------------------------------------------------------------    
            MPI_Barrier( MPI_COMM_WORLD );  
        }
    }
    //-------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------
    if (MD.iRANK == 0)      {           MPI_File_write_at(fp_MPIOUT, 0, &MD.iEQcntr,   1, MPI_INT,      &STATUS);                   }
    
    MPI_Barrier( MPI_COMM_WORLD );  
    MPI_File_close(&fp_MPIOUT);
    if ((WRITESTFOFEACHELEMENT == 1) && (fTestMag >= MD.fMinMag4Prop))
    {    MPI_File_close(&fp_STFOUT);            } 
    MPI_Finalize();
    //-------------------------------------------------------------------------------------
    if (MD.iRANK == 0)             
    {   timer = clock() - timer;     
        double time_taken;
        time_taken  = ((double)timer)/CLOCKS_PER_SEC;
        time_taken /= 60.0;
        fprintf(stdout,"Total RunTime in minutes: %6.2f\n",time_taken);
        fprintf(stdout,"Times iSize =>  total of %6.2f CPU hours\n\n",(time_taken*(float)MD.iSIZE)/60.0);
        fprintf(stdout,"STFcounter: %d\n",MD.iSTFcntr);
    }
    //-------------------------------------------------------------------------------------
    return 0;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void InitializeVariables(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct EQstruct *EQ, struct Kstruct *K, char **argv)
{   int     i;
    float   fTemp0;
    char    ctempVals[512],     cAppend[512];
    char    cFile1_In[512],     cFile2_In[512],         cFile3_In[512]; 
    FILE    *fp1,               *fp2,                   *fp3;           
    //-----------------------------------------------------------------------------------
    strcpy(cFile1_In,  argv[1]);// opening and reading the "run file" that contains specifics about this model run
    if ((fp1 = fopen(cFile1_In,"r")) == NULL)      {   fprintf(stdout,"Error -cant open %s file in initializeVariables function \n",cFile1_In);      exit(10);     }

    if (fgets(ctempVals, 512, fp1) != NULL)        {                                                            }                            
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %s",  MD->cInputName);             }             
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD->iRunNum);                }
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD->iUsePSeis);              } 
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD->iUseProp);               }           
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD->fMinMag4Prop);           }
    MD->fIntSeisDeltT_InYrs = LOADINGSTEPINTSEIS/365.25;
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD->fAftrSlipTime);          }  
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD->fDeepRelaxTime);         }            
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD->fRecLgth);               }
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %f", &MD->fHealFact);              }
    if (fgets(ctempVals, 512, fp1) != NULL)        {   sscanf(ctempVals,"%*s %d", &MD->iSeedStart);             }
    if (MD->iSeedStart <= 0)                       {   MD->iSeedStart = rand();                                 }
    
    fclose(fp1);
    //------------------------------------------------------------------------------------
    MD->fHealFact = MAX(MD->fHealFact, 0.0); //this is to ensure that the healing factor is between 0 and 1
    MD->fHealFact = MIN(MD->fHealFact, 1.0);
    
    //------------------------------------------------------------------------------------
    strcpy(cFile2_In,MD->cInputName);      strcat(cFile2_In,"_");                   sprintf(cAppend, "%d",MD->iRunNum);     strcat(cFile2_In,cAppend);      strcat(cFile2_In,"_Roughn.dat"); 
    strcpy(cFile3_In,MD->cInputName);      strcat(cFile3_In,"_BNDtrig.dat");  
    if ((fp2 = fopen(cFile2_In,"rb"))     == NULL)     {    printf("Error -cant open *_Roughn.dat file.  in Initialize variables...   %s\n", cFile1_In);      exit(10);     }
    
    if (fread(&MD->iFPNum,   sizeof(int),   1, fp2) != 1)       {   exit(10);   }// this is currently read PatchNum 
    if (fread(&MD->iFVNum,   sizeof(int),   1, fp2) != 1)       {   exit(10);   }// this is currently read VertexNum
    if (fread(&MD->fFltLegs, sizeof(float), 1, fp2) != 1)       {   exit(10);   }// this is currently read mean leg length
    if (fread(&fTemp0,       sizeof(float), 1, fp2) != 1)       {   exit(10);   }// this is currently read mean patch area
    
    fclose(fp2);
    //------------------------------------------------------------------------------------
    if ((fp3 = fopen(cFile3_In,"rb"))     == NULL)     
    {   printf("Warning -cant open *_BNDtrig.dat file. in Initialize variables...\n");      
        MD->iBPNum   = 0;
        MD->iBVNum   = 0;
        MD->fBndLegs = 1.0e+9;
    }
    else
    {   if (fread(&MD->iBPNum,   sizeof(int),   1, fp3) != 1)       {   exit(10);   }// this is currently read PatchNum 
        if (fread(&MD->iBVNum,   sizeof(int),   1, fp3) != 1)       {   exit(10);   }// this is currently read VertexNum
        if (fread(&MD->fBndLegs, sizeof(float), 1, fp3) != 1)       {   exit(10);   }// this is currently mean leg length
        if (fread(&fTemp0,       sizeof(float), 1, fp3) != 1)       {   exit(10);   }// this is currently mean patch area 
        fclose(fp3);
    }
    //------------------------------------------------------------------
    MD->fFltLegs    *= 1.0E+3; // now it is in meters 
    MD->fBndLegs    *= 1.0E+3; // now it is in meters 
    //------------------------------------------------------------------
    MD->fUnitSlipF   = MD->fFltLegs *1.0E-4; //take the smaller of the two leg lengths and then a fraction of it (e.g., 1/10000) as value for unit slip => if leg length is 1000m, then the unit slip is 0.1m => idea is to ensure "infinitesimal" slip/strain conditions
    MD->fUnitSlipB   = MD->fBndLegs *1.0E-4;
    MD->iGlobTTmax   = 0; //max. global travel time => is used when "use rupture propagation" is used
    MD->fPSeis_Step  = 0.0;
    MD->iEQcntr      = 0; //counting the number of events that were written to file
    MD->iSTFcntr     = 0;
    MD->iLoadSteps   = LOADSTEPS_POW2;
    MD->iStepNum     = IntPow(2, MD->iLoadSteps);

    EQ->iEndCntr     = 0; //this counter is used to check if an event is really over (at the end of EQ while loop)
    EQ->fMaxDTau     = 0.0; //max change in shear stress during the event
    EQ->fMaxSlip     = 0.0; //max in-plane slip of the event
    //------------------------------------------------------------------------------------
    MD->ivF_START    = (int *) calloc(MD->iSIZE, sizeof(int)); //these starts and offsets relate to how the code accesses "local" and "global" vectors => use the start and offset to 
    MD->ivB_START    = (int *) calloc(MD->iSIZE, sizeof(int));
    MD->ivF_OFFSET   = (int *) calloc(MD->iSIZE, sizeof(int));//locate where each rank will put/get data from when accessing a "global" vector
    MD->ivB_OFFSET   = (int *) calloc(MD->iSIZE, sizeof(int));

    MD->iv_STRT2     = (int *) calloc(MD->iSIZE, sizeof(int)); 
    MD->iv_OFFST2    = (int *) calloc(MD->iSIZE, sizeof(int));
    MD->iv_STRT3     = (int *) calloc(MD->iSIZE, sizeof(int)); 
    MD->iv_OFFST3    = (int *) calloc(MD->iSIZE, sizeof(int)); 

    MD->iF_BASEelem  = (int)(MD->iFPNum/MD->iSIZE); 
    MD->iF_ADDelem   = (int)(MD->iFPNum%MD->iSIZE);
    MD->iB_BASEelem  = (int)(MD->iBPNum/MD->iSIZE); 
    MD->iB_ADDelem   = (int)(MD->iBPNum%MD->iSIZE);
    //---------------------------
    for (i = 0; i < MD->iSIZE;     i++)      {   MD->ivF_OFFSET[i]     = MD->iF_BASEelem;                                   }
    for (i = 0; i < MD->iF_ADDelem;i++)      {   MD->ivF_OFFSET[i]    += 1;                                                 }
    for (i = 1; i < MD->iSIZE;     i++)      {   MD->ivF_START[i]      = MD->ivF_START[i-1] + MD->ivF_OFFSET[i-1];          }
    
    for (i = 0; i < MD->iSIZE;     i++)      {   MD->ivB_OFFSET[i]     = MD->iB_BASEelem;                                   }
    for (i = 0; i < MD->iB_ADDelem;i++)      {   MD->ivB_OFFSET[i]    += 1;                                                 }
    for (i = 1; i < MD->iSIZE;     i++)      {   MD->ivB_START[i]      = MD->ivB_START[i-1] + MD->ivB_OFFSET[i-1];          }
    //------------------------------------------------------------------------------------  
    TR->ivFL_StabT           = ( int *)  calloc(MD->ivF_OFFSET[MD->iRANK],   sizeof( int));     
    TR->ivFL_Activated       = ( int *)  calloc(MD->ivF_OFFSET[MD->iRANK],   sizeof( int));     
    TR->ivFL_Ptch_t0         = ( int *)  calloc(MD->ivF_OFFSET[MD->iRANK],   sizeof( int));
    TR->ivFL_Ptch_tDur       = ( int *)  calloc(MD->ivF_OFFSET[MD->iRANK],   sizeof( int));
    TR->ivFL_FricLaw         = ( int *)  calloc(MD->ivF_OFFSET[MD->iRANK],   sizeof( int));     
    
    TR->ivFG_SegID_temp      = ( int *)  calloc(MD->iFPNum,   sizeof( int));    
    TR->ivFG_FltID_temp      = ( int *)  calloc(MD->iFPNum,   sizeof( int));
    TR->ivFG_Flagged_temp    = ( int *)  calloc(MD->iFPNum,   sizeof( int));
    TR->ivFG_V1_temp         = ( int *)  calloc(MD->iFPNum,   sizeof( int));                        
    TR->ivFG_V2_temp         = ( int *)  calloc(MD->iFPNum,   sizeof( int));                        
    TR->ivFG_V3_temp         = ( int *)  calloc(MD->iFPNum,   sizeof( int));
    TR->fvFG_CentE_temp      = (float *) calloc(MD->iFPNum,   sizeof(float));                       
    TR->fvFG_CentN_temp      = (float *) calloc(MD->iFPNum,   sizeof(float));                       
    TR->fvFG_CentZ_temp      = (float *) calloc(MD->iFPNum,   sizeof(float));   
    TR->fvFG_StressRatetemp  = (float *) calloc(MD->iFPNum,   sizeof(float));                       
    TR->fvFG_SlipRatetemp    = (float *) calloc(MD->iFPNum,   sizeof(float));
    TR->fvFG_Raketemp        = (float *) calloc(MD->iFPNum,   sizeof(float));  
    TR->fvFL_SlipRate_temp   = (float *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof(float));  
    TR->fvFL_SlipRake_temp   = (float *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof(float));  
    
    TR->ivBG_SegID_temp      = ( int *)  calloc(MD->iBPNum,   sizeof( int));
    TR->ivBG_V1_temp         = ( int *)  calloc(MD->iBPNum,   sizeof( int));                        
    TR->ivBG_V2_temp         = ( int *)  calloc(MD->iBPNum,   sizeof( int));                        
    TR->ivBG_V3_temp         = ( int *)  calloc(MD->iBPNum,   sizeof( int));    
    TR->fvBG_CentE_temp      = (float *) calloc(MD->iBPNum,   sizeof(float));                       
    TR->fvBG_CentN_temp      = (float *) calloc(MD->iBPNum,   sizeof(float));                       
    TR->fvBG_CentZ_temp      = (float *) calloc(MD->iBPNum,   sizeof(float));
    //--------------------------------------    
    TR->fvFL_RefNrmStrs      = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_Area            = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_SelfStiffStk    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);          
    TR->fvFL_SelfStiffDip    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_MeanSelfStiff   = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvBL_SelfStiffStk    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);      
    TR->fvBL_SelfStiffDip    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);      
    TR->fvBL_SelfStiffOpn    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);  

    TR->fvFL_RefStaFric      = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_RefDynFric      = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_RefDcVal        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);      
    TR->fvFL_RefStaFric_vari = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_RefDynFric_vari = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_RefDcVal_vari   = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    
    TR->fvFL_StaFric         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_DynFric         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_ArrFric         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_CurFric         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_B4_Fric         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);      
    TR->fvFL_TempRefFric     = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_CurDcVal        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]); 
    TR->fvFL_MaxSlip         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]); 
    
    TR->fvFL_StaFricMod_temp = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);          
    TR->fvFL_DynFricMod_temp = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_NrmStrsMod_temp = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);              
    TR->fvFL_DcMod_temp      = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    
    TR->fvFL_PSeis_T0_F      = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_CurSlpH         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_CurSlpV         = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_AccumSlp        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);  
    TR->fvFL_OverShotStress  = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvBL_PSeis_T0_S      = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);  
    TR->fvBL_PSeis_T0_N      = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);  
    
    TR->fvFL_StrsRateStk     = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                  
    TR->fvFL_StrsRateDip     = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_B4_StrsH        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                  
    TR->fvFL_B4_StrsV        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                  
    TR->fvFL_B4_StrsN        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvFL_CurStrsH        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                  
    TR->fvFL_CurStrsV        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                  
    TR->fvFL_CurStrsN        = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    TR->fvBL_B4_StrsH        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                  
    TR->fvBL_B4_StrsV        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                  
    TR->fvBL_B4_StrsN        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);
    TR->fvBL_CurStrsH        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                  
    TR->fvBL_CurStrsV        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                  
    TR->fvBL_CurStrsN        = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);
    //--------------------------------------
    TR->imFGL_TTP            = gsl_matrix_int_calloc(  MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);          
    TR->imFGL_TTS            = gsl_matrix_int_calloc(  MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);  
    TR->fmFGL_SrcRcvH        = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);      
    TR->fmFGL_SrcRcvV        = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);      
    TR->fmFGL_SrcRcvN        = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);
    //------------------------------------------------------------------
    VT->fvFG_VlX_temp        = (float *) calloc(MD->iFVNum,   sizeof(float));                       
    VT->fvFG_VlY_temp        = (float *) calloc(MD->iFVNum,   sizeof(float));                       
    VT->fvFG_Hght_temp       = (float *) calloc(MD->iFVNum,   sizeof(float));   
    VT->fvFG_PosE_temp       = (float *) calloc(MD->iFVNum,   sizeof(float));                       
    VT->fvFG_PosN_temp       = (float *) calloc(MD->iFVNum,   sizeof(float));                       
    VT->fvFG_PosZ_temp       = (float *) calloc(MD->iFVNum,   sizeof(float));   
    VT->fvBG_PosE_temp       = (float *) calloc(MD->iBVNum,   sizeof(float));                       
    VT->fvBG_PosN_temp       = (float *) calloc(MD->iBVNum,   sizeof(float));                       
    VT->fvBG_PosZ_temp       = (float *) calloc(MD->iBVNum,   sizeof(float));
    //------------------------------------------------------------------------------------
    //first dimension is number of rows (slow direction) => is the RECEIVERS; second dimension is number of columns (fast direction) => is the SOURCES
    K->FF_SS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                 // immer kurz -lang 
    K->FF_SD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                 // SS => source (kurz) = strike; receiver (lang) = strike
    K->FF_SO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                 // SD => source (kurz) = strike; receiver (lang) = dip
    K->FF_DS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                 // FF => first is source, second is receiver    
    K->FF_DD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    
    K->FF_DO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);
    if (MD->iUseProp == 1)
    {   K->FFp_SS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                // immer kurz -lang 
        K->FFp_SD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // SS => source (kurz) = strike; receiver (lang) = strike
        K->FFp_SO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // SD => source (kurz) = strike; receiver (lang) = dip
        K->FFp_DS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // FF => first is source, second is receiver    
        K->FFp_DD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                       
        K->FFp_DO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);
   
        K->FFs_SS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                // immer kurz -lang 
        K->FFs_SD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // SS => source (kurz) = strike; receiver (lang) = strike
        K->FFs_SO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // SD => source (kurz) = strike; receiver (lang) = dip
        K->FFs_DS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                    // FF => first is source, second is receiver    
        K->FFs_DD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);                       
        K->FFs_DO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivF_OFFSET[MD->iRANK]);
    }
    K->FB_SS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->FB_SD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->FB_SO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);
    K->FB_DS  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->FB_DD  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->FB_DO  = gsl_matrix_float_calloc(MD->iFPNum, MD->ivB_OFFSET[MD->iRANK]);
 
    K->BF_SS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    K->BF_SD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    K->BF_DS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    K->BF_DD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    K->BF_OS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    K->BF_OD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivF_OFFSET[MD->iRANK]);                     
    
    K->BB_SS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_SD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_SO  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);
    K->BB_DS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_DD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_DO  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);
    K->BB_OS  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_OD  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);                     
    K->BB_OO  = gsl_matrix_float_calloc(MD->iBPNum, MD->ivB_OFFSET[MD->iRANK]);
    //------------------------------------------------------------------------------------
    EQ->ivR_WrtStrtPos = (  int *) calloc(MD->iSIZE,                 sizeof( int)); 
    EQ->ivR_StfStrtPos = (  int *) calloc(MD->iSIZE,                 sizeof( int));                 
    EQ->ivL_ActPtchID  = (  int *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof( int));
    EQ->ivL_t0ofPtch   = (  int *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof( int)); 
    EQ->ivL_tDurPtch   = (  int *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof( int));                 
    EQ->ivL_StabType   = (  int *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof( int));
    EQ->fvL_PtchSlpH   = (float *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof(float));                
    EQ->fvL_PtchSlpV   = (float *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof(float));  
    EQ->fvL_PtchDTau   = (float *) calloc(MD->ivF_OFFSET[MD->iRANK], sizeof(float));                
    
    EQ->fvL_EQslipH    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);                        
    EQ->fvL_EQslipV    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);
    
    EQ->fvM_MRFvals    = gsl_vector_float_calloc(MAXMOMRATEFUNCLENGTH);
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    if (MD->iRANK == 0) //making sure that the data were imported from file correctly
    {   fprintf(stdout,"Number of RANKS: %d\n",MD->iSIZE);                  
        fprintf(stdout,"System info: Byte Size for FLOAT: %lu     INT: %lu    \n\n", sizeof(float), sizeof(int));
        fprintf(stdout,"FileName:           %s\n",MD->cInputName);                  fprintf(stdout,"RunNumber:          %d\n",MD->iRunNum); 
        fprintf(stdout,"UsePostSeis:        %d\n",MD->iUsePSeis);                   fprintf(stdout,"UseRuptProp:        %d\n",MD->iUseProp);                    
        fprintf(stdout,"MinMag2UseRuptProp: %f\n",MD->fMinMag4Prop);                fprintf(stdout,"IntSeisTStep:       %2.3f\n",MD->fIntSeisDeltT_InYrs*365.25);                   
        fprintf(stdout,"ViscAftSlip:        %3.1f\n",MD->fAftrSlipTime);            fprintf(stdout,"ViscDeepRelax:      %3.1f\n",MD->fDeepRelaxTime);
        fprintf(stdout,"RecLength:          %5.1f\n",MD->fRecLgth);                 fprintf(stdout,"CoSeisHealFraction: %5.1f\n",MD->fHealFact);                
        fprintf(stdout,"SeedLocation:       %5d\n",MD->iSeedStart);
        fTemp0  = 100 - expf(-1.0/MD->fAftrSlipTime)*100.0; //factor/fraction from decay function
        fprintf(stdout,"Fractional post-seismic change during first year: %2.4f percent released \n\n",fTemp0);

        fprintf(stdout,"FaultPatchNumber %d     FaultVertexNumber %d\n", MD->iFPNum, MD->iFVNum);
        fprintf(stdout,"BoundPatchNumber %d     BoundVertexNumber %d\n", MD->iBPNum, MD->iBVNum);
        if (MD->iBPNum > 1)
        {   fprintf(stdout,"MeanLegLenth      %f   %f\n", MD->fFltLegs, MD->fBndLegs);   }
        else
        {   fprintf(stdout,"MeanLegLenth      %f     \n", MD->fFltLegs);    }
    }   
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int LoadInput(struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct EQstruct *EQ)
{   
    gsl_rng *fRandN; // this is pretty much straight from the GSL reference, the default RNG has good performance, so no need to change
    const gsl_rng_type *RandT;
    gsl_rng_env_setup(); 
    RandT   = gsl_rng_default;
    fRandN  = gsl_rng_alloc(RandT);
    unsigned long RSeed =(unsigned long)(MD->iRANK + MD->iSeedStart);   
    gsl_rng_set(fRandN, RSeed);
    //------------------------------------------------------------------------------------
    int     iHasSlip = 0;
    int     i,                  j,                  iGlobPos;
    int     iTemp0,             iTemp1;
    float   fCombMemory;
    float   fTemp0,             fTemp1,             fTemp2,             fTemp3;
    float   fTemp4,             fTemp6,             fTemp7,             fTemp8;
    float   fP1[3],             fP2[3],             fP3[3],             fP1P2[3],       fP1P3[3];       
    float   fvNrm[3],           fvStk[3],           fvDip[3],           feZ[3],         fSrcRcvVect[3];
    
    char    ctempVals[512],     cAppend[512];
    char    cFileName1[512],    cFileName2[512],    cFileName3[512],    cFileName4[512];
    
    FILE    *fp1,                   *fp2,                   *fp3,               *fp4; 
    //------------------------------------------------------------------------------------
    //reading/loading the remaining information from txt and dat files.. this step is done by each rank individually; also, some values are read locally (every rank only reads the information relevant for it) while other values are read globally
    strcpy(cFileName1,MD->cInputName);      strcat(cFileName1,"_Summary_RoughStrength.txt"); 
    strcpy(cFileName2,MD->cInputName);      strcat(cFileName2,"_");                 sprintf(cAppend, "%d",MD->iRunNum);     strcat(cFileName2,cAppend);         strcat(cFileName2,"_Roughn.dat");
    strcpy(cFileName3,MD->cInputName);      strcat(cFileName3,"_");                 sprintf(cAppend, "%d",MD->iRunNum);     strcat(cFileName3,cAppend);         strcat(cFileName3,"_Strgth.dat"); 
    strcpy(cFileName4,MD->cInputName);      strcat(cFileName4,"_BNDtrig.dat");
    //------------------------------------------------------------------------------------
    if ((fp1 = fopen(cFileName1,"r")) == NULL)          {   printf("Error -cant open *.flt file. LoadInputParameter function...\n");      exit(10);     }
    
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fMedDense);                                          }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fAddNrmStrs);    MD->fAddNrmStrs *= -1.0;            }// keep at MPa but make compression negative ==> in MATLAB and the input files, normal stress is still positive for compression => is switched here to compression == negative convention
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fShearMod);      MD->fShearMod   *=  1.0E+9;         }// convert from GPa to Pa (for actual computation of K-matrix and also for magnitude calculation, I need Pa and not MPa)
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %e",&MD->fPoisson);                                           }
    if (fgets(ctempVals, 512, fp1) != NULL)             {   sscanf(ctempVals,"%*s %d",&MD->iChgBtwEQs);                                         }
    
    fclose(fp1);
    //------------------------------------------------------------------------------------     
    MD->fLambda = (2.0*MD->fShearMod*MD->fPoisson)/(1.0-2.0*MD->fPoisson); //elastic parameters => for K-matrix calculation, remember that the code uses a linear elastic half-space
    fTemp0      = (MD->fMedDense > 0.0) ? MD->fMedDense : 2700.0; // if zero density is used for depth gradient (e.g., to have a depth-independent normal stress) then I still!! need to define a density for the wave propagation speed.... => is done here...
    MD->fVp     = sqrtf((MD->fLambda +2.0*MD->fShearMod)/fTemp0); // in m/s
    MD->fVs     = sqrtf(MD->fShearMod/fTemp0); // in m/s 
    //-----------------------------------------------------------------
    if ((fp2 = fopen(cFileName2,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName2);      exit(10);     }
   
    if (fread(&iTemp0, sizeof(int),  1,fp2) != 1)       {   exit(10);   }// have those already 
    if (fread(&iTemp0, sizeof(int),  1,fp2) != 1)       {   exit(10);   }// have those already 
    if (fread(&fTemp0, sizeof(float),1,fp2) != 1)       {   exit(10);   }// have those already 
    if (fread(&fTemp0, sizeof(float),1,fp2) != 1)       {   exit(10);   }// have those already 
    //-----------------------------------------------------
    if (fread(TR->ivFG_V1_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)       {   exit(10);   } //reading the geometric information about the faults
    if (fread(TR->ivFG_V2_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)       {   exit(10);   }
    if (fread(TR->ivFG_V3_temp,    sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)       {   exit(10);   }
    if (fread(TR->ivFG_SegID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)       {   exit(10);   }
    if (fread(TR->ivFG_FltID_temp, sizeof(int),MD->iFPNum,fp2) != MD->iFPNum)       {   exit(10);   }    
    //-----------------------------------------------------
    if (fread(TR->fvFG_StressRatetemp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum) {   exit(10);   }    
    if (fread(TR->fvFG_SlipRatetemp,   sizeof(float),MD->iFPNum,fp2) != MD->iFPNum) {   exit(10);   }    
    if (fread(TR->fvFG_Raketemp,       sizeof(float),MD->iFPNum,fp2) != MD->iFPNum) {   exit(10);   }     
    //-----------------------------------------------------
    if (fread(VT->fvFG_VlX_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)       {   exit(10);   }
    if (fread(VT->fvFG_VlY_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)       {   exit(10);   }
    //-----------------------------------------------------
    if (fread(VT->fvFG_PosE_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)      {   exit(10);   }
    if (fread(VT->fvFG_PosN_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)      {   exit(10);   }
    if (fread(VT->fvFG_PosZ_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)      {   exit(10);   }
    if (fread(VT->fvFG_Hght_temp, sizeof(float),MD->iFVNum,fp2) != MD->iFVNum)      {   exit(10);   }
    //-----------------------------------------------------
    if (fread(TR->fvFG_CentE_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)     {   exit(10);   }
    if (fread(TR->fvFG_CentN_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)     {   exit(10);   }   
    if (fread(TR->fvFG_CentZ_temp, sizeof(float),MD->iFPNum,fp2) != MD->iFPNum)     {   exit(10);   }
    //-------------------------------------     
    fclose(fp2);
    //------------------------------------------------------------------------------------        
    if ((fp3 = fopen(cFileName3,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName3);      exit(10);     }    
    
    fseek(fp3, (1L*sizeof(int)),       SEEK_CUR); // contains the patch number again -> skipped
    
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefStaFric->data,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDynFric->data,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefNrmStrs->data,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])            {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
    
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDcVal->data,        sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])            {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
    //------------------------------------- 
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefStaFric_vari->data, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDynFric_vari->data, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_RefDcVal_vari->data,   sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])            {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
    //------------------------------------- 
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_StaFricMod_temp->data, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
 
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_DynFricMod_temp->data, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])    {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum

    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_NrmStrsMod_temp->data, sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])            {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
    
    fseek(fp3, (1L*sizeof(float)*(long)(MD->ivF_START[MD->iRANK])), SEEK_CUR); //skip the part that is in front of segId for specific RANK
    if (fread(TR->fvFL_DcMod_temp->data,      sizeof(float),MD->ivF_OFFSET[MD->iRANK],fp3) != MD->ivF_OFFSET[MD->iRANK])            {   exit(10);   }
    fseek(fp3, (1L*sizeof(float)*(long)(MD->iFPNum - MD->ivF_START[MD->iRANK] - MD->ivF_OFFSET[MD->iRANK])), SEEK_CUR); //skip over to end of SegIDs => in total i will have moved by MD->iFPNum
    //-------------------------------------
    fclose(fp3);
    //------------------------------------------------------------------------------------     
    if (MD->iBPNum > 0)
    {   if ((fp4 = fopen(cFileName4,"rb")) == NULL)  
        {   printf("Warning -cant open  %s LoadInputParameter function...\n Continue with boundary surface",cFileName4);         
            MD->iBPNum = 0;
            MD->iBVNum = 0;
        }
        else
        {   if (fread(&iTemp0, sizeof(int),  1,fp4) != 1)       {   exit(10);   }// have those already 
            if (fread(&iTemp0, sizeof(int),  1,fp4) != 1)       {   exit(10);   }// have those already 
            if (fread(&fTemp0, sizeof(float),1,fp4) != 1)       {   exit(10);   }// have those already 
            if (fread(&fTemp0, sizeof(float),1,fp4) != 1)       {   exit(10);   }// have those already 
            //-------------------------------------
            if (fread(TR->ivBG_V1_temp,    sizeof(int), MD->iBPNum,fp4) != MD->iBPNum)  {   exit(10);   }
            if (fread(TR->ivBG_V2_temp,    sizeof(int), MD->iBPNum,fp4) != MD->iBPNum)  {   exit(10);   } 
            if (fread(TR->ivBG_V3_temp,    sizeof(int), MD->iBPNum,fp4) != MD->iBPNum)  {   exit(10);   } 
            if (fread(TR->ivBG_SegID_temp, sizeof(int), MD->iBPNum,fp4) != MD->iBPNum)  {   exit(10);   } 
            //-------------------------------------
            fseek(fp1, 2L*sizeof(float)*(long)MD->iBVNum,   SEEK_CUR); // this would be local coordintates (from gridding) of the vertices..; not needed    
            //-------------------------------------
            if (fread(VT->fvBG_PosE_temp,  sizeof(float),MD->iBVNum,fp4) != MD->iBVNum) {   exit(10);   }
            if (fread(VT->fvBG_PosN_temp,  sizeof(float),MD->iBVNum,fp4) != MD->iBVNum) {   exit(10);   }
            if (fread(VT->fvBG_PosZ_temp,  sizeof(float),MD->iBVNum,fp4) != MD->iBVNum) {   exit(10);   }
            //----------------------------------------
            if (fread(TR->fvBG_CentE_temp, sizeof(float),MD->iBPNum,fp4) != MD->iBPNum) {   exit(10);   }
            if (fread(TR->fvBG_CentN_temp, sizeof(float),MD->iBPNum,fp4) != MD->iBPNum) {   exit(10);   }    
            if (fread(TR->fvBG_CentZ_temp, sizeof(float),MD->iBPNum,fp4) != MD->iBPNum) {   exit(10);   }
            //-------------------------------------  
            fclose(fp4);   
    }   }
    //------------------------------------------------------------------------------------     
    for (i = 0; i < (MD->iFVNum); i++)        //converting all km to meters, also shifting the indizes => Matlab starting with "1" while C starting with "0" => hence the -1 used here....    
    {   VT->fvFG_PosE_temp[i]     *= 1000.0;            VT->fvFG_PosN_temp[i]     *= 1000.0;         VT->fvFG_PosZ_temp[i]     *= 1000.0;               VT->fvFG_Hght_temp[i]     *= 1000.0; 
    }
    for (i = 0; i < (MD->iBVNum); i++)            
    {   VT->fvBG_PosE_temp[i]     *= 1000.0;            VT->fvBG_PosN_temp[i]     *= 1000.0;         VT->fvBG_PosZ_temp[i]     *= 1000.0; 
    }
    for (i = 0; i < (MD->iFPNum); i++)            
    {   TR->ivFG_V1_temp[i]       -= 1;                 TR->ivFG_V2_temp[i]       -= 1;              TR->ivFG_V3_temp[i]       -= 1;                  
        TR->fvFG_CentE_temp[i]    *= 1000.0;            TR->fvFG_CentN_temp[i]    *= 1000.0;         TR->fvFG_CentZ_temp[i]    *= 1000.0; 
        TR->ivFG_SegID_temp[i]    -= 1; 
    }
    for (i = 0; i < (MD->iBPNum); i++)            
    {   TR->ivBG_V1_temp[i]        -= 1;                TR->ivBG_V2_temp[i]       -= 1;              TR->ivBG_V3_temp[i]       -= 1;                  
        TR->fvBG_CentE_temp[i]  *= 1000.0;              TR->fvBG_CentN_temp[i]    *= 1000.0;         TR->fvBG_CentZ_temp[i]    *= 1000.0; 
    }
    //------------------------------------------------------------------------------------               
    gsl_vector_float_scale(TR->fvFL_RefDcVal_vari,   0.01);             gsl_vector_float_mul(TR->fvFL_RefDcVal_vari, TR->fvFL_RefDcVal);                        
              
    gsl_vector_float_scale(TR->fvFL_RefStaFric_vari, 0.01);             gsl_vector_float_mul(TR->fvFL_RefStaFric_vari, TR->fvFL_RefStaFric);                        
        
    gsl_vector_float_scale(TR->fvFL_RefDynFric_vari, 0.01);             gsl_vector_float_mul(TR->fvFL_RefDynFric_vari, TR->fvFL_RefDynFric);  
       
    gsl_vector_float_scale(TR->fvFL_NrmStrsMod_temp, -1.0);

    gsl_vector_float_add(TR->fvFL_RefStaFric, TR->fvFL_StaFricMod_temp);
    gsl_vector_float_add(TR->fvFL_RefDynFric, TR->fvFL_DynFricMod_temp);
    gsl_vector_float_add(TR->fvFL_RefNrmStrs, TR->fvFL_NrmStrsMod_temp);
    gsl_vector_float_add(TR->fvFL_RefDcVal,   TR->fvFL_DcMod_temp); 
    //------------------------------------------------------------------------------------     
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) //here the current static/dynamic friction coefficients are determines, also the current Dc value; further, the stressing rates are assigned
    {   
        iGlobPos                    = i + MD->ivF_START[MD->iRANK];
        fTemp0                      = cosf(TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0) *TR->fvFG_StressRatetemp[iGlobPos];
        fTemp1                      = sinf(TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0) *TR->fvFG_StressRatetemp[iGlobPos];
        gsl_vector_float_set(TR->fvFL_StrsRateStk, i, fTemp0); 
        gsl_vector_float_set(TR->fvFL_StrsRateDip, i, fTemp1);  
    
        TR->fvFL_SlipRate_temp[i]   = TR->fvFG_SlipRatetemp[iGlobPos]/1000.0; //was defined in mm/yr, butneed in m/y => /1000;
        TR->fvFL_SlipRake_temp[i]   = TR->fvFG_Raketemp[iGlobPos]*M_PI/180.0; //rake angle in radian
                
        if (fabs(TR->fvFL_SlipRate_temp[i]) > 0.0)          {       iHasSlip = 1;       } 
        //------------------------//------------------------    
        fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        fTemp1                     = gsl_vector_float_get(TR->fvFL_RefDcVal, i)    + gsl_vector_float_get(TR->fvFL_RefDcVal_vari,   i)*fTemp2;
        gsl_vector_float_set(TR->fvFL_CurDcVal,i, fTemp1);
        //------------------------
        fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        fTemp1                     = gsl_vector_float_get(TR->fvFL_RefStaFric, i)  + gsl_vector_float_get(TR->fvFL_RefStaFric_vari, i)*fTemp2; 
        gsl_vector_float_set(TR->fvFL_StaFric, i, fTemp1);     
        gsl_vector_float_set(TR->fvFL_CurFric, i, fTemp1);
        //------------------------
        fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
        fTemp3                     = gsl_vector_float_get(TR->fvFL_RefDynFric, i)  + gsl_vector_float_get(TR->fvFL_RefDynFric_vari, i)*fTemp2; 
        gsl_vector_float_set(TR->fvFL_DynFric, i, fTemp3);
        fTemp4                     = fTemp3 - (OVERSHOOTFRAC-1.0) *fabs(fTemp1 -fTemp3); //this friction is the "arrest friction" i.e., the dynamic overshoot friction value (the "lowest posssible")     
        gsl_vector_float_set(TR->fvFL_ArrFric, i, fTemp4);
        
        fTemp6                     = (fTemp3 - fTemp4)*-1.0*gsl_vector_float_get(TR->fvFL_RefNrmStrs,i); //this is "dynamic fric" minus "arrest fric" multiplied with normal stress => gives amount of excess above arrest fric to continue sliding => corresponds to "dyn fric" level
        gsl_vector_float_set(TR->fvFL_OverShotStress, i, fTemp6);
        //-------------------------------------------------------------------------------
        GetVertices(VT,TR,iGlobPos, 0, fP1, fP2, fP3); 
        GetVector(fP1,fP2, fP1P2);
        GetVector(fP1,fP3, fP1P3);
        CrossProduct(fP1P2,fP1P3, fvNrm);
        //-----------------------------------------          
        fTemp0 = 0.5*sqrtf( fvNrm[0]*fvNrm[0] +fvNrm[1]*fvNrm[1] +fvNrm[2]*fvNrm[2]); //is needed for moment/magnitue calculation
        gsl_vector_float_set(TR->fvFL_Area, i, fTemp0);    
    }
    //------------------------------------------------------------------------------------     
    if (USEVPVSVELOCITY == 1)        {   MD->fVpVsRatio = MD->fVp/MD->fVs;       }          else                            {   MD->fVpVsRatio = 1.0;                           }              
    MD->fRealDeltT = 0.7*MD->fFltLegs/MD->fVp;
    if (MD->iUseProp == 1)           {   MD->fDeltT     = MD->fRealDeltT;        }          else                            {   MD->fDeltT     = FLT_MAX;                       }
     
    //------------------------------------------------------------------------------------     
    //------------------------------------------------------------------------------------
    for (i = 0; i < MD->iFPNum;  i++) //sources
    {   //----------------------------------------- 
        GetVertices(VT,TR,i, 0, fP1, fP2, fP3); 
        GetVector(fP1,fP2, fP1P2);
        GetVector(fP1,fP3, fP1P3);
        CrossProduct(fP1P2,fP1P3, fvNrm);   
        NormalizeVector(fvNrm);
        //----------------------------------------- 
        //determine local coordinate system for currently selected (local) fault patch
        feZ[0] = 0.0;                       feZ[1] = 0.0;                               feZ[2] = 1.0;        
        //-----------------------------------------  
        CrossProduct(feZ,fvNrm,fvStk);
        // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) 
        fTemp0  = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
        if (fabs(fTemp0) < FLT_EPSILON)
        {   fvStk[0]= 0.0;                      fvStk[1] = 1.0;                     fvStk[2] = 0.0;             }
        else
        {   fvStk[0]    = fvStk[0]/fTemp0;      fvStk[1] = fvStk[1]/fTemp0;         fvStk[2] = fvStk[2]/fTemp0; }
        //----------------------------------------- 
        CrossProduct(fvNrm, fvStk, fvDip);
        NormalizeVector(fvDip);
        //-------------------------------------------------------------------- 
        for (j = 0; j < MD->ivF_OFFSET[MD->iRANK]; j++) //receiver
        {   
            iGlobPos         = j + MD->ivF_START[MD->iRANK];    
            fSrcRcvVect[0]   = TR->fvFG_CentE_temp[iGlobPos] - TR->fvFG_CentE_temp[i]; // this is vector from current patch (source, global) to receiver (local); this one points therefore from source to the receiver
            fSrcRcvVect[1]   = TR->fvFG_CentN_temp[iGlobPos] - TR->fvFG_CentN_temp[i]; //will be needed for vector projection => when using rupture propagation => which component of slip vector is in direct line towards receiver patch (mode II) and which component is perpendicular to that (mode III)
            fSrcRcvVect[2]   = TR->fvFG_CentZ_temp[iGlobPos] - TR->fvFG_CentZ_temp[i]; //this is therefore for rupture propagation and the transient signals that arrise when using different velocities for mode II and mode III

            fTemp7   =(1.0/MD->fVp)*sqrtf( (fSrcRcvVect[0]*fSrcRcvVect[0]) + (fSrcRcvVect[1]*fSrcRcvVect[1]) + (fSrcRcvVect[2]*fSrcRcvVect[2]) ); // Vp travel time between the pair 
            fTemp8   = fTemp7*MD->fVpVsRatio; //this is the Vs travel tim
            //----------------------------------------- 
            iTemp0   = (int)(fTemp7/MD->fDeltT); //means it's rounding downm, cutting off whatever floating point contribution the travel times had
            iTemp1   = (int)(fTemp8/MD->fDeltT); 
            gsl_matrix_int_set(TR->imFGL_TTP, i, j, iTemp0);
            gsl_matrix_int_set(TR->imFGL_TTS, i, j, iTemp1);
            //----------------------------------------- 
            MD->iGlobTTmax   = MAX(MD->iGlobTTmax, iTemp1);   
            //----------------------------------------- 
         
            if (iGlobPos == i)
            {   gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, 0.0);
                gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, 0.0);
                gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, 0.0);
            }
            else
            {   NormalizeVector(fSrcRcvVect);
                fTemp0          = fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2]; //this rotates the current source-to-receiver vector (in global coordinates)= fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2];
                fTemp1          = fvStk[0]*fSrcRcvVect[0] + fvStk[1]*fSrcRcvVect[1] + fvStk[2]*fSrcRcvVect[2]; //into the local coordinate system of the source patch(that is the idea...)
                fTemp2          = fvDip[0]*fSrcRcvVect[0] + fvDip[1]*fSrcRcvVect[1] + fvDip[2]*fSrcRcvVect[2];         

                gsl_matrix_float_set(TR->fmFGL_SrcRcvN, i, j, fTemp0);
                gsl_matrix_float_set(TR->fmFGL_SrcRcvH, i, j, fTemp1);
                gsl_matrix_float_set(TR->fmFGL_SrcRcvV, i, j, fTemp2);           
    }   }   }
    //------------------------------------------------------------------------------------     
    //------------------------------------------------------------------------------------
    if (MD->iRANK == 0)
    {   fprintf(stdout,"Medium density:              %5.2f\n",MD->fMedDense);
        fprintf(stdout,"Shear modulus:               %5.2e\n",MD->fShearMod);
        fprintf(stdout,"Poisson ratio:               %5.2f\n",MD->fPoisson);
        fprintf(stdout,"Lambda:                      %5.2e\n",MD->fLambda);
        fprintf(stdout,"P-wave velocity:             %5.2f\n",MD->fVp);
        fprintf(stdout,"S-wave velocity:             %5.2f\n",MD->fVs);
        fprintf(stdout,"Added normal stress (MPa):   %5.2f\n",MD->fAddNrmStrs);
        fprintf(stdout,"Change friction between EQs: %d\n",MD->iChgBtwEQs);
    }
    //------------------------------------------------------------------------------------
    MPI_Allreduce(MPI_IN_PLACE, &MD->iGlobTTmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    MD->iMaxSTFlgth = MD->iGlobTTmax +1;//use the +1 here to ensure I didn't mess this up and to make sure that the overwriting of the STF does not happen too fast (b/c I'm looping over the STF to save memory)
//--------------------------------
//--------------------------------
//--------------------------------
//--------------------------------
if (WRITESTFOFEACHELEMENT == 1) {   MD->iMaxSTFlgth =  MAXMOMRATEFUNCLENGTH;    }
//--------------------------------
//--------------------------------
//--------------------------------
//--------------------------------
    if (MD->iRANK == 0) {           fprintf(stdout,"Other stuff after loading input and defining more parameters\nMaxSTFlgth: %d      and  %d\nDelta time: %4.4f   and Unit slip: %4.4f \n\n", MD->iMaxSTFlgth, MD->iGlobTTmax, MD->fRealDeltT, MD->fUnitSlipF);                }
    //----------------------------------
    EQ->fmFSL_H = gsl_matrix_float_calloc(MD->iMaxSTFlgth, MD->ivF_OFFSET[MD->iRANK]);
    EQ->fmFSL_V = gsl_matrix_float_calloc(MD->iMaxSTFlgth, MD->ivF_OFFSET[MD->iRANK]);
    EQ->fmFSL_N = gsl_matrix_float_calloc(MD->iMaxSTFlgth, MD->ivF_OFFSET[MD->iRANK]);      
    
    EQ->fmSTF_slip     = gsl_matrix_float_calloc(MD->ivF_OFFSET[MD->iRANK], MD->iMaxSTFlgth);       
    EQ->fmSTF_strength = gsl_matrix_float_calloc(MD->ivF_OFFSET[MD->iRANK], MD->iMaxSTFlgth);   
    //------------------------------------------------------------------------------------     
    //------------------------------------------------------------------------------------
    fCombMemory  = 0.0; 
    fCombMemory += (float)(sizeof(  int)*( 3*MD->ivF_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*(37*MD->ivF_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*( 9*MD->ivB_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(  int)*( 6*MD->iFPNum));
    fCombMemory += (float)(sizeof(  int)*( 4*MD->iBPNum));
    fCombMemory += (float)(sizeof(float)*( 6*MD->iFPNum));
    fCombMemory += (float)(sizeof(float)*( 3*MD->iBPNum));
    fCombMemory += (float)(sizeof(  int)*( 2*MD->iFPNum *MD->ivF_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*(23*MD->iFPNum *MD->ivF_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*( 6*MD->iFPNum *MD->ivB_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*( 9*MD->iBPNum *MD->ivB_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)*( 6*MD->iBPNum *MD->ivF_OFFSET[MD->iRANK]));
    fCombMemory += (float)(sizeof(float)* (5*MD->iMaxSTFlgth* MD->ivF_OFFSET[MD->iRANK]));
    
    fCombMemory /= 1.0E+9;
    
    if (MD->iRANK == 0) {           fprintf(stdout,"Approx. memory per rank: %3.4fGb => %3.4fGb in total\n",fCombMemory, fCombMemory*(float)MD->iSIZE);     }
    return iHasSlip;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void Build_K_Matrix( struct MDstruct *MD, struct TRstruct *TR, struct VTstruct *VT, struct  Kstruct *K)
{   
    //------------------------------------------------------------------------------------  
    int     i,                  j,              iSrcFltID,          iRcvFltID,          iGlobPos,       iCntFlags_F = 0;         
    float   fTemp0,             fTemp1,         fTemp2,             fTemp3,             fTemp4,         fTemp5;
    float   fTemp6,             fTemp7,         fTemp8,             fTempB,             fTempC;
    float   fTempD,             fTempE,         fTempF,             fTempG,             fMinDist,       fSrcLgth;
    float   fP1s[3],            fP2s[3],        fP3s[3],            fP1r[3],            fP2r[3],        fP3r[3];
    float   fvNrm[3],           fvStk[3],       fvDip[3],           fSrcPt[3],          fRcvPt[3],      fStress[6],         fStressOut[6];
    //------------------------------------------------------------------------------------          
    for (i = 0; i < MD->iFPNum;  i++)   //going through the sources
    {   GetVertices(VT, TR, i, 0, fP1s, fP2s, fP3s); 
        iSrcFltID    = TR->ivFG_FltID_temp[i];     
        fSrcPt[0]    = TR->fvFG_CentE_temp[i];              fSrcPt[1] = TR->fvFG_CentN_temp[i];             fSrcPt[2] = TR->fvFG_CentZ_temp[i];
        fSrcLgth     = GetDistance(fP1s, fP2s);             fSrcLgth += GetDistance(fP2s, fP3s);            fSrcLgth += GetDistance(fP3s, fP1s);                    fSrcLgth /= 3.0;    
        //-----------------------------------------------
        for (j = 0; j < MD->ivF_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivF_START[MD->iRANK];  
            GetVertices(VT, TR, iGlobPos, 0, fP1r, fP2r, fP3r);  
            iRcvFltID    = TR->ivFG_FltID_temp[iGlobPos];        
            fRcvPt[0]    = TR->fvFG_CentE_temp[iGlobPos];   fRcvPt[1]= TR->fvFG_CentN_temp[iGlobPos];       fRcvPt[2]= TR->fvFG_CentZ_temp[iGlobPos];    
            //----------------------------------------------- 
            fMinDist     = GetMinDist_8Pnts(fSrcPt, fP1s, fP2s, fP3s, fRcvPt,  fP1r, fP2r, fP3r);     //the distance between the 6 (2 times three vertices) vertices of the triangles and their center locations; I want/need to ensure that patch that are too close to each other are flagged i.e., not used in that interaction to avoid numerical instability (Medhi's code is not entirely artefact free after all...)     
            //-----------------------------------------------   
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);   //to rotate into the coordinate system of the receiver!   //to rotate into the coordinate system of the receiver!
            //-------------------------------  
            fTemp0 = 0.0;       fTemp1 = 0.0;       fTemp2 = 0.0;       fTemp3 = 0.0;       fTemp4 = 0.0;       fTemp5 = 0.0;           fTemp6 = 0.0;       fTemp7 = 0.0;       fTemp8 = 0.0;       

            if ( (iSrcFltID == iRcvFltID) || (fMinDist > CUTDISTANCE*fSrcLgth) ) //if "self" or distance is larger than 1.0 UsedLgLength's then I can proceed as normal
            {   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, MD->fUnitSlipF, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = (fStressOut[0]/MD->fUnitSlipF)*1.0e-6;           fTemp1 = (fStressOut[1]/MD->fUnitSlipF)*1.0e-6;           fTemp2 = (fStressOut[2]/MD->fUnitSlipF)*1.0e-6; 
                //-----------------------------------------------  
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, MD->fUnitSlipF, 0.0, MD->fShearMod, MD->fLambda);         // slip in dip      
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = (fStressOut[0]/MD->fUnitSlipF)*1.0e-6;           fTemp4 = (fStressOut[1]/MD->fUnitSlipF)*1.0e-6;           fTemp5 = (fStressOut[2]/MD->fUnitSlipF)*1.0e-6; 
                //-----------------------------------------------  
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, 0.0, MD->fUnitSlipF, MD->fShearMod, MD->fLambda);         // slip in open         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = (fStressOut[0]/MD->fUnitSlipF)*1.0e-6;           fTemp7 = (fStressOut[1]/MD->fUnitSlipF)*1.0e-6;           fTemp8 = (fStressOut[2]/MD->fUnitSlipF)*1.0e-6; 
                //-----------------------------------------------  
                if (iGlobPos == i)     {      gsl_vector_float_set(TR->fvFL_SelfStiffStk, j, fTemp1);               gsl_vector_float_set(TR->fvFL_SelfStiffDip, j, fTemp5);        }       
            } 
            else //if too close, then go here
            {   TR->ivFG_Flagged_temp[i] = 1;           iCntFlags_F++;
            }  
            gsl_matrix_float_set(K->FF_SO, i, j,  fTemp0);   gsl_matrix_float_set(K->FF_SS, i, j,  fTemp1);    gsl_matrix_float_set(K->FF_SD, i, j,  fTemp2);  
            gsl_matrix_float_set(K->FF_DO, i, j,  fTemp3);   gsl_matrix_float_set(K->FF_DS, i, j,  fTemp4);    gsl_matrix_float_set(K->FF_DD, i, j,  fTemp5);     
           
            if (MD->iUseProp == 1)
            {   //http://sites.science.oregonstate.edu/math/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/dotprod/dotprod.html
                //https://en.wikipedia.org/wiki/Vector_projection
                //----------------------------------------
                fTempB = gsl_matrix_float_get(TR->fmFGL_SrcRcvH, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvH, i, j); //these are the components of slip vector in the  
                fTempC = gsl_matrix_float_get(TR->fmFGL_SrcRcvH, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvV, i, j); //source-receiver direction (in local source element coordinate system)
                fTempD = gsl_matrix_float_get(TR->fmFGL_SrcRcvH, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvN, i, j); //in other words, this is the H/V/N P-wave component
                fTempE = fTempB*fTemp0 + fTempC*fTemp3 + fTempD*fTemp6; //slip_stk*K_so +slip_dip*K_do +slip_opn*K_oo
                fTempF = fTempB*fTemp1 + fTempC*fTemp4 + fTempD*fTemp7; //slip_stk*K_ss +slip_dip*K_ds +slip_opn*K_os
                fTempG = fTempB*fTemp2 + fTempC*fTemp5 + fTempD*fTemp8; //slip_stk*K_sd +slip_dip*K_dd +slip_opn*K_od               
                gsl_matrix_float_set(K->FFp_SO, i, j,  fTempE);   gsl_matrix_float_set(K->FFp_SS, i, j,  fTempF);    gsl_matrix_float_set(K->FFp_SD, i, j,  fTempG);  
                //--------------------------------------------------------------------
                fTempB = 1.0 - fTempB; //this is the S-component; being perpendicular to P component
                fTempC = 0.0 - fTempC; //hence, S-comp = OrigSlipVect - P-comp
                fTempD = 0.0 - fTempD;
                    
                fTempE = fTempB*fTemp0 + fTempC*fTemp3 + fTempD*fTemp6; 
                fTempF = fTempB*fTemp1 + fTempC*fTemp4 + fTempD*fTemp7;
                fTempG = fTempB*fTemp2 + fTempC*fTemp5 + fTempD*fTemp8;                     
                gsl_matrix_float_set(K->FFs_SO, i, j,  fTempE);   gsl_matrix_float_set(K->FFs_SS, i, j,  fTempF);    gsl_matrix_float_set(K->FFs_SD, i, j,  fTempG);  
                //--------------------------------------------------------------------  
                //--------------------------------------------------------------------  
                fTempB = gsl_matrix_float_get(TR->fmFGL_SrcRcvV, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvH, i, j); //these are the components of slip vector in the 
                fTempC = gsl_matrix_float_get(TR->fmFGL_SrcRcvV, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvV, i, j); //source-receiver direction (in local source element coordinate system)
                fTempD = gsl_matrix_float_get(TR->fmFGL_SrcRcvV, i, j)*gsl_matrix_float_get(TR->fmFGL_SrcRcvN, i, j); //in other words, this is the H/V/N P-wave component
                fTempE = fTempB*fTemp0 + fTempC*fTemp3 + fTempD*fTemp6; //slip_stk*K_so +slip_dip*K_do +slip_opn*K_oo
                fTempF = fTempB*fTemp1 + fTempC*fTemp4 + fTempD*fTemp7; //slip_stk*K_ss +slip_dip*K_ds +slip_opn*K_os
                fTempG = fTempB*fTemp2 + fTempC*fTemp5 + fTempD*fTemp8; //slip_stk*K_sd +slip_dip*K_dd +slip_opn*K_od               
                gsl_matrix_float_set(K->FFp_DO, i, j,  fTempE);   gsl_matrix_float_set(K->FFp_DS, i, j,  fTempF);    gsl_matrix_float_set(K->FFp_DD, i, j,  fTempG);  
                //--------------------------------------------------------------------
                fTempB =  0.0 - fTempB; //this is the S-component; being perpendicular to P component
                fTempC =  1.0 - fTempC; //hence, S-comp = OrigSlipVect - P-comp
                fTempD =  0.0 - fTempD;
                    
                fTempE = fTempB*fTemp0 + fTempC*fTemp3 + fTempD*fTemp6; 
                fTempF = fTempB*fTemp1 + fTempC*fTemp4 + fTempD*fTemp7;
                fTempG = fTempB*fTemp2 + fTempC*fTemp5 + fTempD*fTemp8;                     
                gsl_matrix_float_set(K->FFs_DO, i, j,  fTempE);   gsl_matrix_float_set(K->FFs_DS, i, j,  fTempF);    gsl_matrix_float_set(K->FFs_DD, i, j,  fTempG);  
                //--------------------------------------------------------------------  
        }   }
        //--------------------------------------------------------------------------------
        for (j = 0; j <  MD->ivB_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivB_START[MD->iRANK];     
            GetVertices(VT, TR, iGlobPos, 1, fP1r, fP2r, fP3r);  
            fRcvPt[0]    = TR->fvBG_CentE_temp[iGlobPos];   fRcvPt[1] = TR->fvBG_CentN_temp[iGlobPos];      fRcvPt[2] = TR->fvBG_CentZ_temp[iGlobPos];    
            //----------------------------------------------- 
            fMinDist     = GetMinDist_8Pnts(fSrcPt, fP1s, fP2s, fP3s, fRcvPt,  fP1r, fP2r, fP3r);     //the distance between the 6 (2 times three vertices) vertices of the triangles and their center locations; I want/need to ensure that patch that are too close to each other are flagged i.e., not used in that interaction to avoid numerical instability (Medhi's code is not entirely artefact free after all...)     
            //-----------------------------------------------   
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);   //to rotate into the coordinate system of the receiver!
            //-------------------------------  
            
            fTemp0 = 0.0;       fTemp1 = 0.0;       fTemp2 = 0.0;       fTemp3 = 0.0;       fTemp4 = 0.0;       fTemp5 = 0.0;       
            
            if (fMinDist > CUTDISTANCE*fSrcLgth) //if "self" or distance is larger than 1.0 UsedLgLength's then I can proceed as normal
            {   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, MD->fUnitSlipF, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = (fStressOut[0]/MD->fUnitSlipF)*1.0e-6;           fTemp1 = (fStressOut[1]/MD->fUnitSlipF)*1.0e-6;           fTemp2 = (fStressOut[2]/MD->fUnitSlipF)*1.0e-6; 
                //-----------------------------------------------  
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, MD->fUnitSlipF, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = (fStressOut[0]/MD->fUnitSlipF)*1.0e-6;           fTemp4 = (fStressOut[1]/MD->fUnitSlipF)*1.0e-6;           fTemp5 = (fStressOut[2]/MD->fUnitSlipF)*1.0e-6; 
                //-----------------------------------------------  
            }  
            else
            {   TR->ivFG_Flagged_temp[i] = 1;           iCntFlags_F++;
            }    
            gsl_matrix_float_set(K->FB_SO, i, j,  fTemp0);   gsl_matrix_float_set(K->FB_SS, i, j,  fTemp1);   gsl_matrix_float_set(K->FB_SD, i, j,  fTemp2);  
            gsl_matrix_float_set(K->FB_DO, i, j,  fTemp3);   gsl_matrix_float_set(K->FB_DS, i, j,  fTemp4);   gsl_matrix_float_set(K->FB_DD, i, j,  fTemp5);            
        }   
    }
    //------------------------------------------------------------------------------------  
    //------------------------------------------------------------------------------------          
    for (i = 0; i <MD->iBPNum; i++)   //going through the sources
    {   
        GetVertices(VT, TR, i, 1, fP1s, fP2s, fP3s); 
        iSrcFltID    = TR->ivBG_SegID_temp[i];     
        fSrcPt[0]    = TR->fvBG_CentE_temp[i];              fSrcPt[1]= TR->fvBG_CentN_temp[i];              fSrcPt[2]= TR->fvBG_CentZ_temp[i];
        fSrcLgth     = GetDistance(fP1s, fP2s);             fSrcLgth+= GetDistance(fP2s, fP3s);             fSrcLgth+= GetDistance(fP3s, fP1s);                     fSrcLgth    /= 3.0;    
        //----------------------------------------------- 
        for (j = 0; j < MD->ivB_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivB_START[MD->iRANK];       
            GetVertices(VT, TR, iGlobPos, 1, fP1r, fP2r, fP3r);  
            iRcvFltID    = TR->ivBG_SegID_temp[iGlobPos];        
            fRcvPt[0]    = TR->fvBG_CentE_temp[iGlobPos];   fRcvPt[1]= TR->fvBG_CentN_temp[iGlobPos];       fRcvPt[2]= TR->fvBG_CentZ_temp[iGlobPos];        
            //-----------------------------------------------
            fMinDist     = GetMinDist_8Pnts(fSrcPt, fP1s, fP2s, fP3s, fRcvPt,  fP1r, fP2r, fP3r);            
            //-----------------------------------------------
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);   //to rotate into the coordinate system of the receiver!
           //----------------------------------------------- 
              
            fTemp0 = 0.0;       fTemp1 = 0.0;       fTemp2 = 0.0;           fTemp3 = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;       fTemp6 = 0.0;       fTemp7 = 0.0;           fTemp8 = 0.0; 

            if ( (iSrcFltID == iRcvFltID) || (fMinDist >= CUTDISTANCE*fSrcLgth) )    //is "self" or if distance is far enough away to not need the extra distance check 
            {   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, MD->fUnitSlipB, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp1 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp2 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //----------------------------------------------- 
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, MD->fUnitSlipB, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp4 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp5 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //-----------------------------------------------
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, 0.0, MD->fUnitSlipB, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp7 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp8 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //----------------------------------------------- 
                if (iGlobPos == i)     
                {       gsl_vector_float_set(TR->fvBL_SelfStiffOpn, j, fTemp6);           gsl_vector_float_set(TR->fvBL_SelfStiffStk, j, fTemp1);               gsl_vector_float_set(TR->fvBL_SelfStiffDip, j, fTemp5);       
                }
            }              
            gsl_matrix_float_set(K->BB_SO, i, j,  fTemp0);   gsl_matrix_float_set(K->BB_SS, i, j,  fTemp1);   gsl_matrix_float_set(K->BB_SD, i, j,  fTemp2);  
            gsl_matrix_float_set(K->BB_DO, i, j,  fTemp3);   gsl_matrix_float_set(K->BB_DS, i, j,  fTemp4);   gsl_matrix_float_set(K->BB_DD, i, j,  fTemp5);  
            gsl_matrix_float_set(K->BB_OO, i, j,  fTemp6);   gsl_matrix_float_set(K->BB_OS, i, j,  fTemp7);   gsl_matrix_float_set(K->BB_OD, i, j,  fTemp8);           
        }    
        //--------------------------------------------------------------------------------
        for (j = 0; j < MD->ivF_OFFSET[MD->iRANK]; j++)// going through the receivers
        {   iGlobPos     = j + MD->ivF_START[MD->iRANK];      
            GetVertices(VT, TR, iGlobPos, 0, fP1r, fP2r, fP3r);        
            fRcvPt[0]    = TR->fvFG_CentE_temp[iGlobPos];   fRcvPt[1]= TR->fvFG_CentN_temp[iGlobPos];       fRcvPt[2]= TR->fvFG_CentZ_temp[iGlobPos];        
            //-----------------------------------------------
            fMinDist     = GetMinDist_8Pnts(fSrcPt, fP1s, fP2s, fP3s, fRcvPt,  fP1r, fP2r, fP3r);            
            //-----------------------------------------------
            GetLocKOS_inKmatrix(fvNrm, fvStk, fvDip, fP1r, fP2r, fP3r);   //to rotate into the coordinate system of the receiver!
           //----------------------------------------------- 
            
            fTemp0 = 0.0;       fTemp1 = 0.0;       fTemp2 = 0.0;           fTemp3 = 0.0;           fTemp4 = 0.0;           fTemp5 = 0.0;       fTemp6 = 0.0;       fTemp7 = 0.0;           fTemp8 = 0.0; 

            if (fMinDist > CUTDISTANCE*fSrcLgth) //if "self" or distance is larger than 1.0 UsedLgLength's then I can proceed as normal
            {   
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, MD->fUnitSlipB, 0.0, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp0 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp1 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp2 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //-----------------------------------------------
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, MD->fUnitSlipB, 0.0, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp3 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp4 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp5 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //-----------------------------------------------  
                StrainHS_Nikkhoo(fStress, fStressOut, fRcvPt[0], fRcvPt[1], fRcvPt[2], fP1s, fP2s, fP3s, 0.0, 0.0, MD->fUnitSlipB, MD->fShearMod, MD->fLambda);         // slip in stk         
                //-----------------------------------------------  
                RotateTensor(fStressOut, fStress, fvNrm, fvStk, fvDip, 0);
                //-----------------------------------------------  
                fTemp6 = (fStressOut[0]/MD->fUnitSlipB)*1.0e-6;           fTemp7 = (fStressOut[1]/MD->fUnitSlipB)*1.0e-6;           fTemp8 = (fStressOut[2]/MD->fUnitSlipB)*1.0e-6; 
                //----------------------------------------------- 
            }
            else
            {   TR->ivFG_Flagged_temp[iGlobPos] = 1;           iCntFlags_F++;
            }
            gsl_matrix_float_set(K->BF_SS, i, j,  fTemp1);  gsl_matrix_float_set(K->BF_SD, i, j,  fTemp2);   
            gsl_matrix_float_set(K->BF_DS, i, j,  fTemp4);  gsl_matrix_float_set(K->BF_DD, i, j,  fTemp5);     
            gsl_matrix_float_set(K->BF_OS, i, j,  fTemp7);  gsl_matrix_float_set(K->BF_OD, i, j,  fTemp8);         
    }   }
    //------------------------------------------------------------------------------------    
    //------------------------------------------------------------------------------------   
    gsl_vector_float_memcpy(TR->fvFL_MeanSelfStiff, TR->fvFL_SelfStiffStk);
    gsl_vector_float_add(TR->fvFL_MeanSelfStiff,    TR->fvFL_SelfStiffDip);
    gsl_vector_float_scale(TR->fvFL_MeanSelfStiff, 0.5);

    MD->fMeanStiffness = cblas_sasum(MD->ivF_OFFSET[MD->iRANK], TR->fvFL_MeanSelfStiff->data, 1);
    MD->fMeanStiffness = MD->fMeanStiffness/(float)MD->ivF_OFFSET[MD->iRANK];

    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++)
    {   // determine friction difference, use together with reference normal stress to get stress drop; get corresponding slip amount, compare with Dc => assign stability type                 
        //------------------------  
        fTemp2 = (gsl_vector_float_get(TR->fvFL_DynFric, i) - gsl_vector_float_get(TR->fvFL_StaFric, i)) *-1.0*gsl_vector_float_get(TR->fvFL_RefNrmStrs, i); //gives a negative shear stress b/c of (dynF - staF)
        fTemp3 = fTemp2/gsl_vector_float_get(TR->fvFL_MeanSelfStiff, i); // b/c "self" has negative sign, the two negatives give a positive slip amount... (for weakening case, when dyn < stat fric) 
        if (fTemp3 > gsl_vector_float_get(TR->fvFL_CurDcVal, i))    {   TR->ivFL_StabT[i] = 1;      }
        else    
        {   if (fTemp2 < 0.0)                                       {   TR->ivFL_StabT[i] = 2;      }
            else                                                    {   TR->ivFL_StabT[i] = 3;      }      
    }   }
    //------------------------------------------------------------------------------------  
    //------------------------------------------------------------------------------------  
    MPI_Allreduce(MPI_IN_PLACE, TR->ivFG_Flagged_temp, MD->iFPNum, MPI_INT, MPI_MAX, MPI_COMM_WORLD); //for output, to write out which fault patches have been flagged
    MPI_Allreduce(MPI_IN_PLACE, &MD->fMeanStiffness, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MD->fMeanStiffness /= (float)MD->iSIZE;
    return;
} 
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetSlipLoadingAndWritePreRunData(struct MDstruct *MD, struct VTstruct *VT, struct TRstruct *TR, struct Kstruct *K, char *cFile2_Out, int HaveBoundarySlip)
{
    int   i, k;
    float fTemp0, fTemp1;
    float fTempA, fTempB, fTempC;
    FILE *fpPre;
    //------------------------------------------------------------------------------------  
    gsl_vector_float    *fvFL_Temp0    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);    
    gsl_vector_float    *fvFL_Temp1    = gsl_vector_float_calloc(MD->ivF_OFFSET[MD->iRANK]);    
    
    gsl_vector_int      *ivFG_Temp0    = gsl_vector_int_calloc(MD->iFPNum);
    gsl_vector_float    *fvFG_Temp0    = gsl_vector_float_calloc(MD->iFPNum);                               
    gsl_vector_float    *fvFG_Temp1    = gsl_vector_float_calloc(MD->iFPNum);                                       
    gsl_vector_float    *fvFG_Temp2    = gsl_vector_float_calloc(MD->iFPNum);                                           
                                
    gsl_vector_float    *fvBL_Temp0    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                    
    gsl_vector_float    *fvBL_Temp1    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                            
    gsl_vector_float    *fvBL_Temp2    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);
    gsl_vector_float    *fvBL_Temp3    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                    
    gsl_vector_float    *fvBL_Temp4    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);                            
    gsl_vector_float    *fvBL_Temp5    = gsl_vector_float_calloc(MD->ivB_OFFSET[MD->iRANK]);
    gsl_vector_float    *fvBG_Temp0    = gsl_vector_float_calloc(MD->iBPNum);                               
    gsl_vector_float    *fvBG_Temp1    = gsl_vector_float_calloc(MD->iBPNum);                                       
    gsl_vector_float    *fvBG_Temp2    = gsl_vector_float_calloc(MD->iBPNum);       
    gsl_vector_float    *fvBG_Temp3    = gsl_vector_float_calloc(MD->iBPNum);                               
    gsl_vector_float    *fvBG_Temp4    = gsl_vector_float_calloc(MD->iBPNum);                                       
    gsl_vector_float    *fvBG_Temp5    = gsl_vector_float_calloc(MD->iBPNum);           
    //------------------------------------------------------------------------------------ 
    if (HaveBoundarySlip == 1)
    {  for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) 
        {   fTemp0   =  -1.0*cosf(TR->fvFL_SlipRake_temp[i])*TR->fvFL_SlipRate_temp[i]; //the -1.0 is here b/c this is back-slip...
            fTemp1   =  -1.0*sinf(TR->fvFL_SlipRake_temp[i])*TR->fvFL_SlipRate_temp[i];
            gsl_vector_float_set(fvFL_Temp0, i, fTemp0);  //the strike slip component            
            gsl_vector_float_set(fvFL_Temp1, i, fTemp1);  //the dip slip component    
        }
        //--------------------------------------------------------------------------------  
        //this here is the normal/classic back-slip method (if no boundary faults are used but a slip-boundary condition is defined)
        
        MPI_Allgatherv(fvFL_Temp0->data, MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp0->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fvFL_Temp1->data, MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp1->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
        //-------------------------------------
        cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->FF_SS->data, MD->ivF_OFFSET[MD->iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);    
        gsl_vector_float_add(TR->fvFL_CurStrsH, fvFL_Temp0);                                                    
        cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->FF_SD->data, MD->ivF_OFFSET[MD->iRANK], fvFG_Temp0->data, 1, 0.0, fvFL_Temp0->data, 1);                        
        gsl_vector_float_add(TR->fvFL_CurStrsV, fvFL_Temp0);        
        
        cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->FF_DS->data, MD->ivF_OFFSET[MD->iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);                                
        gsl_vector_float_add(TR->fvFL_CurStrsH, fvFL_Temp0);
        cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->FF_DD->data, MD->ivF_OFFSET[MD->iRANK], fvFG_Temp1->data, 1, 0.0, fvFL_Temp0->data, 1);                        
        gsl_vector_float_add(TR->fvFL_CurStrsV, fvFL_Temp0);                                                                                                    
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        if (MD->iBPNum > 0)
        {   //----------------------------------------------------------------------------
            fTempA = 0.0;
            for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++)     
            {   if (TR->ivFL_StabT[i] < 3)      //add all induced stresses from patches that are NOT stable (use cond. stable and unstable)
                {   fTempA += sqrtf(gsl_vector_float_get(TR->fvFL_CurStrsH, i)*gsl_vector_float_get(TR->fvFL_CurStrsH, i) + gsl_vector_float_get(TR->fvFL_CurStrsV, i)*gsl_vector_float_get(TR->fvFL_CurStrsV, i) );
            }   }
            MPI_Allreduce(MPI_IN_PLACE, &fTempA, 1 , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //total applied stressing rate
            //----------------------------------------------------------------------------
            gsl_vector_float_scale(fvFG_Temp0, -1.0);                                       gsl_vector_float_scale(fvFG_Temp1, -1.0);    //it appears that I need to use the actual slip direction here...;instead of "backslip" i mean
            //----------------------------------------------------------------------------
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_SS->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);
            gsl_vector_float_add(TR->fvBL_CurStrsH, fvBL_Temp0);                                
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_SD->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);    
            gsl_vector_float_add(TR->fvBL_CurStrsV, fvBL_Temp0);                    
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_SO->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);                            
            gsl_vector_float_add(TR->fvBL_CurStrsN, fvBL_Temp0);    
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_DS->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);                                
            gsl_vector_float_add(TR->fvBL_CurStrsH, fvBL_Temp0);
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_DD->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);    
            gsl_vector_float_add(TR->fvBL_CurStrsV, fvBL_Temp0);                                            
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iFPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->FB_DO->data, MD->ivB_OFFSET[MD->iRANK], fvFG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);                            
            gsl_vector_float_add(TR->fvBL_CurStrsN, fvBL_Temp0);    
                        
            gsl_vector_float_set_zero(TR->fvFL_CurStrsH);                                   gsl_vector_float_set_zero(TR->fvFL_CurStrsV);
                        
            for (k = 0; k < MAXITERATION4BOUNDARY; k++) // release the stress on boundary elements iteratively to determine corresponding slip-rate on boundary fault elements
            {   //--------------------------------------------------
                gsl_vector_float_memcpy(fvBL_Temp0, TR->fvBL_CurStrsH);                     gsl_vector_float_div(fvBL_Temp0, TR->fvBL_SelfStiffStk);                    gsl_vector_float_scale(fvBL_Temp0, -1.0);
                gsl_vector_float_memcpy(fvBL_Temp1, TR->fvBL_CurStrsV);                     gsl_vector_float_div(fvBL_Temp1, TR->fvBL_SelfStiffDip);                    gsl_vector_float_scale(fvBL_Temp1, -1.0);
                gsl_vector_float_memcpy(fvBL_Temp2, TR->fvBL_CurStrsN);                     gsl_vector_float_div(fvBL_Temp2, TR->fvBL_SelfStiffOpn);                    gsl_vector_float_scale(fvBL_Temp2, -1.0);
                gsl_vector_float_add(fvBL_Temp3, fvBL_Temp0);                               gsl_vector_float_add(fvBL_Temp4, fvBL_Temp1);                               gsl_vector_float_add(fvBL_Temp5, fvBL_Temp2);   
                
                MPI_Allgatherv(fvBL_Temp0->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp0->data, MD->ivB_OFFSET, MD->ivB_START ,MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Allgatherv(fvBL_Temp1->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp1->data, MD->ivB_OFFSET, MD->ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Allgatherv(fvBL_Temp2->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp2->data, MD->ivB_OFFSET, MD->ivB_START ,MPI_FLOAT, MPI_COMM_WORLD);
                
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_SS->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);    
                gsl_vector_float_add(TR->fvBL_CurStrsH, fvBL_Temp0);                            
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_SD->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);    
                gsl_vector_float_add(TR->fvBL_CurStrsV, fvBL_Temp0);                    
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_SO->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp0->data, 1, 0.0, fvBL_Temp0->data, 1);                            
                gsl_vector_float_add(TR->fvBL_CurStrsN, fvBL_Temp0);
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_DS->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);    
                gsl_vector_float_add(TR->fvBL_CurStrsH, fvBL_Temp0);                                                    
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_DD->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);    
                gsl_vector_float_add(TR->fvBL_CurStrsV, fvBL_Temp0);                                            
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_DO->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp1->data, 1, 0.0, fvBL_Temp0->data, 1);                            
                gsl_vector_float_add(TR->fvBL_CurStrsN, fvBL_Temp0);
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_OS->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);    
                gsl_vector_float_add(TR->fvBL_CurStrsH, fvBL_Temp0);                                                    
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_OD->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);
                gsl_vector_float_add(TR->fvBL_CurStrsV, fvBL_Temp0);                                                
                cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivB_OFFSET[MD->iRANK], 1.0, K->BB_OO->data, MD->ivB_OFFSET[MD->iRANK], fvBG_Temp2->data, 1, 0.0, fvBL_Temp0->data, 1);                            
                gsl_vector_float_add(TR->fvBL_CurStrsN, fvBL_Temp0);
            }
            //--------------------------------------------------
            
            MPI_Allgatherv(fvBL_Temp3->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp3->data, MD->ivB_OFFSET, MD->ivB_START ,MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Allgatherv(fvBL_Temp4->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp4->data, MD->ivB_OFFSET, MD->ivB_START, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Allgatherv(fvBL_Temp5->data, MD->ivB_OFFSET[MD->iRANK], MPI_FLOAT, fvBG_Temp5->data, MD->ivB_OFFSET, MD->ivB_START ,MPI_FLOAT, MPI_COMM_WORLD);
            //--------------------------------------------------
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_SS->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp3->data, 1, 0.0, fvFL_Temp0->data, 1);    
            gsl_vector_float_add(TR->fvFL_CurStrsH, fvFL_Temp0);                            
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_SD->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp3->data, 1, 0.0, fvFL_Temp0->data, 1);                                            
            gsl_vector_float_add(TR->fvFL_CurStrsV, fvFL_Temp0);            
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_DS->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp4->data, 1, 0.0, fvFL_Temp0->data, 1);
            gsl_vector_float_add(TR->fvFL_CurStrsH, fvFL_Temp0);                                                            
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_DD->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp4->data, 1, 0.0, fvFL_Temp0->data, 1);                                            
            gsl_vector_float_add(TR->fvFL_CurStrsV, fvFL_Temp0);            
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_OS->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp5->data, 1, 0.0, fvFL_Temp0->data, 1);        
            gsl_vector_float_add(TR->fvFL_CurStrsH, fvFL_Temp0);                                                    
            cblas_sgemv(CblasRowMajor,CblasTrans, MD->iBPNum, MD->ivF_OFFSET[MD->iRANK], 1.0, K->BF_OD->data, MD->ivF_OFFSET[MD->iRANK], fvBG_Temp5->data, 1, 0.0, fvFL_Temp0->data, 1);                                            
            gsl_vector_float_add(TR->fvFL_CurStrsV, fvFL_Temp0);
            //----------------------------------------------------------------------------
            fTempB = 0;
            for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++)     
            {   if (TR->ivFL_StabT[i] < 3)      //add all induced stresses from patches that are NOT stable (use cond. stable and unstable)
                {   fTempB += sqrtf(gsl_vector_float_get(TR->fvFL_CurStrsH, i)*gsl_vector_float_get(TR->fvFL_CurStrsH, i) + gsl_vector_float_get(TR->fvFL_CurStrsV, i)*gsl_vector_float_get(TR->fvFL_CurStrsV, i) );
            }   }
            MPI_Allreduce(MPI_IN_PLACE, &fTempB, 1 , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //total applied stressing rate
            //----------------------------------------------------------------------------
            fTempC = fTempA/fTempB; //this is the actual scaling factor to be applied...
            gsl_vector_float_scale(TR->fvFL_CurStrsH, fTempC); //these are just temporary place holders (not currrent stress but stressing rate put in here) at the moment...
            gsl_vector_float_scale(TR->fvFL_CurStrsV, fTempC); //the values are actually stressing rates, will be assigned within the next few lines...
        }
        gsl_vector_float_add(TR->fvFL_StrsRateStk, TR->fvFL_CurStrsH);          
        gsl_vector_float_add(TR->fvFL_StrsRateDip, TR->fvFL_CurStrsV);          
    }
    //------------------------------------------------------------------------------------
    fTempA = 0.0;
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++)
    {   fTempA += sqrtf(gsl_vector_float_get(TR->fvFL_StrsRateStk, i)*gsl_vector_float_get(TR->fvFL_StrsRateStk, i) +gsl_vector_float_get(TR->fvFL_StrsRateDip, i)*gsl_vector_float_get(TR->fvFL_StrsRateDip, i)); //the amount of applied stressing rate (per interseismic time step)
    }
    MPI_Allreduce(MPI_IN_PLACE, &fTempA, 1 , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); //total applied stressing rate
    fTempA /= (float)MD->iFPNum;            
    if (MD->iRANK == 0)             {       fprintf(stdout,"Resulting average stressing-rate on faults (MPa/yr): %5.5f and per time step: %5.6f\n",fTempA, (fTempA*MD->fIntSeisDeltT_InYrs));           }
    //------------------------------------------------------------------------------------
    gsl_vector_float_scale(TR->fvFL_StrsRateStk, MD->fIntSeisDeltT_InYrs ); //reference stressing rate is now including the time step => only need to add that amount
    gsl_vector_float_scale(TR->fvFL_StrsRateDip, MD->fIntSeisDeltT_InYrs ); //when applicable (when stepping forward by one interseis. time step)
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------  
    MPI_Barrier( MPI_COMM_WORLD );
    //following, the pre-run output is written; this is really helpful for debugging etc... and also to check for example that the stressing-rate distribution on the fault is making sense (is possible that grid resolution of boundary fault is low and those patches are too close to EQfaults, negatively affecting the loading function)
    if (MD->iRANK == 0)           
    {   if ((fpPre = fopen(cFile2_Out,"wb")) == NULL)      {        exit(10);       }

        fwrite( &MD->iFPNum,        sizeof(int),            1, fpPre);          fwrite( &MD->iFVNum,        sizeof(int),           1, fpPre); 
        fwrite( &MD->iBPNum,        sizeof(int),            1, fpPre);          fwrite( &MD->iBVNum,        sizeof(int),           1, fpPre); 
              
        fwrite( TR->ivFG_V1_temp,   sizeof(int),   MD->iFPNum, fpPre);          fwrite( TR->ivFG_V2_temp,   sizeof(int),   MD->iFPNum, fpPre);          fwrite( TR->ivFG_V3_temp,   sizeof(int),   MD->iFPNum, fpPre);
        fwrite( TR->ivFG_SegID_temp,sizeof(int),   MD->iFPNum, fpPre);          fwrite( TR->ivFG_FltID_temp,sizeof(int),   MD->iFPNum, fpPre);

        fwrite( VT->fvFG_PosE_temp, sizeof(float), MD->iFVNum, fpPre);          fwrite( VT->fvFG_PosN_temp, sizeof(float), MD->iFVNum, fpPre);
        fwrite( VT->fvFG_PosZ_temp, sizeof(float), MD->iFVNum, fpPre);          fwrite( VT->fvFG_Hght_temp, sizeof(float), MD->iFVNum, fpPre);
        
        fwrite( TR->fvFG_CentE_temp,sizeof(int),   MD->iFPNum, fpPre);          fwrite( TR->fvFG_CentN_temp,sizeof(int),   MD->iFPNum, fpPre);          fwrite( TR->fvFG_CentZ_temp,sizeof(int),   MD->iFPNum, fpPre);
        
        fwrite( TR->ivBG_V1_temp,   sizeof(int),   MD->iBPNum, fpPre);          fwrite( TR->ivBG_V2_temp,   sizeof(int),   MD->iBPNum, fpPre);          fwrite( TR->ivBG_V3_temp,   sizeof(int),   MD->iBPNum, fpPre);

        fwrite( VT->fvBG_PosE_temp, sizeof(float), MD->iBVNum, fpPre);          fwrite( VT->fvBG_PosN_temp, sizeof(float), MD->iBVNum, fpPre);          fwrite( VT->fvBG_PosZ_temp, sizeof(float), MD->iBVNum, fpPre);  
    }
    //-----------------------------------------------------------
    MPI_Allgatherv(TR->fvFL_StrsRateStk->data, MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp0->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgatherv(TR->fvFL_StrsRateDip->data, MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp1->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    if (MD->iRANK == 0)           
    {   fwrite(TR->ivFG_Flagged_temp, sizeof(int),  MD->iFPNum, fpPre); //the "flagged" status for interaction => if removed/overridden etc...
        fwrite(fvFG_Temp0->data,     sizeof(float), MD->iFPNum, fpPre); //the stressing rate in strike direction
        fwrite(fvFG_Temp1->data,     sizeof(float), MD->iFPNum, fpPre); //the stressing rate in dip direction
    }
    //-----------------------------------------------------------
    gsl_vector_float_memcpy(fvFL_Temp0, TR->fvFL_StaFric);                                                                  gsl_vector_float_mul(fvFL_Temp0, TR->fvFL_RefNrmStrs);      gsl_vector_float_scale(fvFL_Temp0, -1.0);
    gsl_vector_float_memcpy(fvFL_Temp1, TR->fvFL_StaFric);      gsl_vector_float_sub(fvFL_Temp1, TR->fvFL_DynFric);         gsl_vector_float_mul(fvFL_Temp1, TR->fvFL_RefNrmStrs);      gsl_vector_float_scale(fvFL_Temp1, -1.0);
    
    MPI_Allgatherv(TR->ivFL_StabT,          MD->ivF_OFFSET[MD->iRANK], MPI_INT,   ivFG_Temp0->data, MD->ivF_OFFSET, MD->ivF_START, MPI_INT,   MPI_COMM_WORLD);
    MPI_Allgatherv(fvFL_Temp0->data,        MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp0->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgatherv(fvFL_Temp1->data,        MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp1->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgatherv(TR->fvFL_CurDcVal->data, MD->ivF_OFFSET[MD->iRANK], MPI_FLOAT, fvFG_Temp2->data, MD->ivF_OFFSET, MD->ivF_START, MPI_FLOAT, MPI_COMM_WORLD);
    if (MD->iRANK == 0)           
    {   fwrite(fvFG_Temp0->data,     sizeof(float), MD->iFPNum, fpPre); //the patch strength
        fwrite(fvFG_Temp1->data,     sizeof(float), MD->iFPNum, fpPre); //the patch stress drop when having full drop
        fwrite(fvFG_Temp2->data,     sizeof(float), MD->iFPNum, fpPre); //the Dc value
        fwrite(ivFG_Temp0->data,     sizeof(int),   MD->iFPNum, fpPre); //the patch stability type (1 = unstable, 2 = cond. stable, 3 = stable)
    }
    //-----------------------------------------------------------
    if (MD->iRANK == 0)           
    {   fwrite(fvBG_Temp3->data, sizeof(float), MD->iBPNum, fpPre); //the slip/stressing rate in strike direction
        fwrite(fvBG_Temp4->data, sizeof(float), MD->iBPNum, fpPre); //the slip/stressing rate in dip direction
        fwrite(fvBG_Temp5->data, sizeof(float), MD->iBPNum, fpPre); //the slip/stressing rate in opening direction
        fclose(fpPre);
    }
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    return; 
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
MPI_Offset StartEQcatalogFile(MPI_File fp_MPIOUT, struct MDstruct *MD, struct VTstruct *VT, struct TRstruct *TR)
{   MPI_Status STATUS;  
        if (MD->iRANK == 0)           
    {   MPI_File_write(fp_MPIOUT, &MD->iEQcntr,                  1, MPI_INT,    &STATUS);       MPI_File_write(fp_MPIOUT, &MD->iFPNum,                   1, MPI_INT,    &STATUS);
        MPI_File_write(fp_MPIOUT, &MD->iFVNum,                   1, MPI_INT,    &STATUS);       MPI_File_write(fp_MPIOUT, &MD->fShearMod,                1, MPI_FLOAT,  &STATUS);    
        MPI_File_write(fp_MPIOUT, &MD->fRealDeltT,                   1, MPI_FLOAT,  &STATUS);      

        MPI_File_write(fp_MPIOUT, TR->ivFG_V1_temp,     MD->iFPNum, MPI_INT,    &STATUS);       MPI_File_write(fp_MPIOUT, TR->ivFG_V2_temp,     MD->iFPNum, MPI_INT,    &STATUS);   
        MPI_File_write(fp_MPIOUT, TR->ivFG_V3_temp,     MD->iFPNum, MPI_INT,    &STATUS);       MPI_File_write(fp_MPIOUT, TR->ivFG_SegID_temp,  MD->iFPNum, MPI_INT,    &STATUS);   
        MPI_File_write(fp_MPIOUT, TR->ivFG_FltID_temp,  MD->iFPNum, MPI_INT,    &STATUS);

        MPI_File_write(fp_MPIOUT, VT->fvFG_VlX_temp,    MD->iFVNum, MPI_FLOAT,  &STATUS);       MPI_File_write(fp_MPIOUT, VT->fvFG_VlY_temp,    MD->iFVNum, MPI_FLOAT,  &STATUS);

        MPI_File_write(fp_MPIOUT, VT->fvFG_PosE_temp,   MD->iFVNum, MPI_FLOAT,  &STATUS);       MPI_File_write(fp_MPIOUT, VT->fvFG_PosN_temp,   MD->iFVNum, MPI_FLOAT,  &STATUS);
        MPI_File_write(fp_MPIOUT, VT->fvFG_PosZ_temp,   MD->iFVNum, MPI_FLOAT,  &STATUS);       MPI_File_write(fp_MPIOUT, VT->fvFG_Hght_temp,   MD->iFVNum, MPI_FLOAT,  &STATUS);

        MPI_File_write(fp_MPIOUT, TR->fvFG_CentE_temp,  MD->iFPNum, MPI_FLOAT,  &STATUS);       MPI_File_write(fp_MPIOUT, TR->fvFG_CentN_temp,  MD->iFPNum, MPI_FLOAT,  &STATUS);
        MPI_File_write(fp_MPIOUT, TR->fvFG_CentZ_temp,  MD->iFPNum, MPI_FLOAT,  &STATUS);
    }        
    MPI_Offset OFFSETAll = 3*sizeof(int) + 2*sizeof(float) +5*MD->iFPNum*sizeof(int) +6*MD->iFVNum*sizeof(float) +3*MD->iFPNum*sizeof(float);
    return OFFSETAll;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void FreeSomeMemory(struct TRstruct *TR, struct VTstruct *VT)
{   //------------------------------------------------------------------------------------      
    free(TR->fvFL_SlipRate_temp);   free(TR->fvFL_SlipRake_temp);   free(TR->ivFG_Flagged_temp);    free(VT->fvFG_VlX_temp);        free(VT->fvFG_VlY_temp);        free(VT->fvFG_Hght_temp);
    free(TR->ivFG_V1_temp);         free(TR->ivFG_V2_temp);         free(TR->ivFG_V3_temp);         free(TR->ivFG_FltID_temp);      free(TR->ivFG_SegID_temp);      free(TR->ivBG_SegID_temp);  
    free(TR->fvFG_StressRatetemp);  free(TR->fvFG_SlipRatetemp);     free(TR->fvFG_Raketemp);
    free(TR->ivBG_V1_temp);         free(TR->ivBG_V2_temp);         free(TR->ivBG_V3_temp);         
    free(TR->fvFG_CentE_temp);      free(TR->fvFG_CentN_temp);      free(TR->fvFG_CentZ_temp);      
    free(TR->fvBG_CentE_temp);      free(TR->fvBG_CentN_temp);      free(TR->fvBG_CentZ_temp);      
    free(TR->fvFL_StaFricMod_temp); free(TR->fvFL_DynFricMod_temp); free(TR->fvFL_NrmStrsMod_temp); free(TR->fvFL_DcMod_temp);
    free(VT->fvFG_PosE_temp);       free(VT->fvFG_PosN_temp);       free(VT->fvFG_PosZ_temp);       free(VT->fvBG_PosE_temp);       free(VT->fvBG_PosN_temp);       free(VT->fvBG_PosZ_temp);   
    //------------------------------------------------------------------------------------  
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void PreloadTheFault(struct MDstruct *MD, struct TRstruct *TR, float LoadFactor)
{   int i;
    float fTemp0, fTemp1,   fTemp2;
    //-------------------------------
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) 
    {   fTemp0 = sqrtf( gsl_vector_float_get(TR->fvFL_StrsRateStk,i)*gsl_vector_float_get(TR->fvFL_StrsRateStk,i) + gsl_vector_float_get(TR->fvFL_StrsRateDip,i)*gsl_vector_float_get(TR->fvFL_StrsRateDip,i) );
        fTemp1 = gsl_vector_float_get(TR->fvFL_CurFric, i)*-1.0*gsl_vector_float_get(TR->fvFL_RefNrmStrs, i);
        
        fTemp2 = fTemp1/fTemp0 *LoadFactor *gsl_vector_float_get(TR->fvFL_StrsRateStk,i); 
        gsl_vector_float_set(TR->fvFL_CurStrsH, i, fTemp2);
        fTemp2 = fTemp1/fTemp0 *LoadFactor *gsl_vector_float_get(TR->fvFL_StrsRateDip,i); 
        gsl_vector_float_set(TR->fvFL_CurStrsV, i, fTemp2);
  
        gsl_vector_float_set(TR->fvFL_CurStrsN, i, gsl_vector_float_get(TR->fvFL_RefNrmStrs, i));   
    }
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetUpdatedFriction(int StabType, float B4Fric, float RefFric, float CurFric, float ArrFric, float CurD_c, float AccumSlip, float PrevSlip, float HealingFact)
{   float NextFric = CurFric,       fTemp;
    //-----------------------------------
    if (PrevSlip > 0.0) 
    {   fTemp    = AccumSlip/CurD_c < 1.0 ? AccumSlip/CurD_c : 1.0;                     
        NextFric = ArrFric + (1.0 - fTemp)*(RefFric - ArrFric);             
    }
    else
    {   if (StabType == 3)          {               NextFric = CurFric;                                             }
        else                        {               NextFric = CurFric + HealingFact*(B4Fric -CurFric);             }
    }

    return NextFric;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
MPI_Offset WriteCatalogFile(MPI_File fp_MPIOUT, MPI_Offset OFFSETall, struct MDstruct *MD, struct TRstruct *TR, struct EQstruct *EQ, int iPrevEQtime, float fMagn)
{   
    int    i, iTemp0;
    int iAddToDuration = 0;
    float  fTemp1, fTemp2,  fTemp5, fTemp6, fTemp3, fTemp4;
    double dEQtime, dEQtimeDiff;
    
    MPI_Status STATUS;
    //----------------------------------------------------------------
    EQ->iActFPNum    = 0;                           
    EQ->fMaxSlip     = 0.0;                 
    EQ->fMaxDTau     = 0.0;     
    //----------------------------------------------------------------
    MD->iEQcntr++;          
    dEQtime     = (double)(MD->iTimeYears)              *(double)MD->fIntSeisDeltT_InYrs;
    dEQtimeDiff = (double)(MD->iTimeYears - iPrevEQtime)*(double)MD->fIntSeisDeltT_InYrs;
    //----------------------------------------------------------------      
    for (i = 0; i < MD->iSIZE; i++)                 {       EQ->ivR_WrtStrtPos[i] = 0;                  EQ->ivR_StfStrtPos[i] = 0;          }     
    //----------------------
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) 
    {   
        if (TR->ivFL_Activated[i] == 1)     
        {     
            fTemp1  = gsl_vector_float_get(EQ->fvL_EQslipH,   i); //horizontal slip of patch
            fTemp2  = gsl_vector_float_get(EQ->fvL_EQslipV,   i); //vertical slip of patch
            fTemp3  = gsl_vector_float_get(TR->fvFL_B4_StrsH, i);               fTemp4  = gsl_vector_float_get(TR->fvFL_B4_StrsV, i);
            fTemp5  = gsl_vector_float_get(TR->fvFL_CurStrsH, i);               fTemp6  = gsl_vector_float_get(TR->fvFL_CurStrsV, i);

            EQ->ivL_ActPtchID[EQ->iActFPNum] = i + MD->ivF_START[MD->iRANK]; 
            EQ->ivL_t0ofPtch[EQ->iActFPNum]  = TR->ivFL_Ptch_t0[i];
            EQ->ivL_tDurPtch[EQ->iActFPNum]  = TR->ivFL_Ptch_tDur[i]+iAddToDuration;
            EQ->fvL_PtchDTau[EQ->iActFPNum]  = sqrtf(fTemp3*fTemp3 +fTemp4*fTemp4) - sqrtf(fTemp5*fTemp5 +fTemp6*fTemp6); //stress DROP is therefore POSITIVE
            //EQ->fvL_PtchDTau[EQ->iActFPNum]  = sqrtf(fTemp5*fTemp5 +fTemp6*fTemp6); //stress DROP is therefore POSITIVE
            EQ->fvL_PtchSlpH[EQ->iActFPNum]  = fTemp1;
            EQ->fvL_PtchSlpV[EQ->iActFPNum]  = fTemp2;          
            EQ->ivL_StabType[EQ->iActFPNum]  = TR->ivFL_StabT[i];           
            EQ->ivR_WrtStrtPos[MD->iRANK]   += 1;
            EQ->ivR_StfStrtPos[MD->iRANK]   += 2*(TR->ivFL_Ptch_tDur[i]+iAddToDuration)*sizeof(float);
                
            EQ->fMaxDTau                     =    (EQ->fvL_PtchDTau[EQ->iActFPNum] > EQ->fMaxDTau) ?    EQ->fvL_PtchDTau[EQ->iActFPNum] : EQ->fMaxDTau;
            EQ->fMaxSlip                     = (sqrt(fTemp1*fTemp1 +fTemp2*fTemp2) > EQ->fMaxSlip) ? sqrt(fTemp1*fTemp1 +fTemp2*fTemp2) : EQ->fMaxSlip;  
            EQ->iActFPNum++;
        }   
    }   
    //-----------------------------------------------------------------------
    MPI_Allreduce(&EQ->iActFPNum, &EQ->iCmbFPNum,      1       , MPI_INT,   MPI_SUM, MPI_COMM_WORLD); //this is maximum stress drop 
    MPI_Allreduce(MPI_IN_PLACE, &EQ->iMRFLgth,         1       , MPI_INT,   MPI_MAX, MPI_COMM_WORLD); //the length of the MRF
    MPI_Allreduce(MPI_IN_PLACE, EQ->ivR_WrtStrtPos, MD->iSIZE  , MPI_INT,   MPI_SUM, MPI_COMM_WORLD); //the length of the MRF
    MPI_Allreduce(MPI_IN_PLACE, EQ->ivR_StfStrtPos, MD->iSIZE  , MPI_INT,   MPI_SUM, MPI_COMM_WORLD); //the length of the MRF   
    MPI_Allreduce(MPI_IN_PLACE, &EQ->fMaxSlip,         1       , MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); //this is maximum slip
    MPI_Allreduce(MPI_IN_PLACE, &EQ->fMaxDTau,         1       , MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); //this is maximum stress drop
    MPI_Allreduce(MPI_IN_PLACE, EQ->fvM_MRFvals->data,MAXMOMRATEFUNCLENGTH, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    //-----------------------------------------------------------------------     
    for (i = 1; i < MD->iSIZE; i++)           {    EQ->ivR_WrtStrtPos[i] += EQ->ivR_WrtStrtPos[i-1];          EQ->ivR_StfStrtPos[i] += EQ->ivR_StfStrtPos[i-1];   }    
    EQ->iCombSTFoffs =  EQ->ivR_StfStrtPos[MD->iSIZE-1];
    for (i = (MD->iSIZE-1); i > 0; i--)       {    EQ->ivR_WrtStrtPos[i]  = EQ->ivR_WrtStrtPos[i-1];          EQ->ivR_StfStrtPos[i]  = EQ->ivR_StfStrtPos[i-1];   }    
    EQ->ivR_WrtStrtPos[0]  = 0;
    EQ->ivR_StfStrtPos[0]  = 0;                 
    //-----------------------------------------------------------------------
    if ((MD->iRANK == 0) && (PLOT2SCREEN == 1))
    {   iTemp0 = (int)fMagn - 4; //just for nicer plotting...
        fprintf(stdout,"# %6d   time: %5.3f   D(time): %2.4f   MRF: %4d   max(u): %4.2f   max(tau): %4.2f   Element #: %5d   ",MD->iEQcntr, dEQtime, dEQtimeDiff, EQ->iMRFLgth, EQ->fMaxSlip, EQ->fMaxDTau,EQ->iCmbFPNum);
        for (i = 0; i < iTemp0; i++)    {       fprintf(stdout,"    ");         }
        fprintf(stdout,"M:  %3.2f \n", fMagn);              
    }
    //-----------------------------------------------------------------------
    if (MD->iRANK == 0)
    {   MPI_File_write_at(fp_MPIOUT,        0,                                                   &MD->iEQcntr,             1,       MPI_INT,    &STATUS);
        MPI_File_write_at(fp_MPIOUT,   OFFSETall,                                                &dEQtime,                 1,       MPI_DOUBLE, &STATUS);  //Earthquake time
        MPI_File_write_at(fp_MPIOUT,   OFFSETall + sizeof(double),                               &fMagn,                   1,       MPI_FLOAT,  &STATUS);  //Earthquake magnitude
        MPI_File_write_at(fp_MPIOUT,   OFFSETall + sizeof(double) +sizeof(float),                &EQ->iCmbFPNum,           1,       MPI_INT,    &STATUS);  //#of fault patches participating in EQ
        MPI_File_write_at(fp_MPIOUT,   OFFSETall + sizeof(double) +sizeof(float) +1*sizeof(int), &EQ->iMRFLgth,            1,       MPI_INT,    &STATUS);  //length of moment rate function "time steps"
        MPI_File_write_at(fp_MPIOUT,   OFFSETall + sizeof(double) +sizeof(float) +2*sizeof(int), EQ->fvM_MRFvals->data, EQ->iMRFLgth, MPI_FLOAT,&STATUS);  //moment rate function values..
    }
    OFFSETall += 2*sizeof(int) +(1 +EQ->iMRFLgth)*sizeof(float) +sizeof(double);
    //--------------------------------------------
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_ActPtchID,EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(int);
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_t0ofPtch, EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(int);
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(float),EQ->fvL_PtchDTau, EQ->iActFPNum, MPI_FLOAT, &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(float);
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(float),EQ->fvL_PtchSlpH, EQ->iActFPNum, MPI_FLOAT, &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(float);
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(float),EQ->fvL_PtchSlpV, EQ->iActFPNum, MPI_FLOAT, &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(float);
    MPI_File_write_at(fp_MPIOUT,   OFFSETall +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_StabType, EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETall += EQ->iCmbFPNum*sizeof(int);

    MPI_Barrier( MPI_COMM_WORLD );
            
    return OFFSETall;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
MPI_Offset WriteEQ_STFsFile(MPI_File fp_STFOUT, MPI_Offset OFFSETstf, struct MDstruct *MD, struct TRstruct *TR, struct EQstruct *EQ, float fMagn)
{   //--------------------------------------------
    int i,  j;
    int iStartPos, TempOffs = 0;
    int iLPatchPos;
    
    float  fNextVal;
    double dEQtime = (double)(MD->iTimeYears)*(double)MD->fIntSeisDeltT_InYrs;
    
    MPI_Status STATUS;

    gsl_vector_float *STF1   = gsl_vector_float_calloc(MD->iMaxSTFlgth);
    gsl_vector_float *STF2   = gsl_vector_float_calloc(MD->iMaxSTFlgth);
        
    MD->iSTFcntr++;  
    //--------------------------------------------
    if (MD->iRANK == 0)
    {   MPI_File_write_at(fp_STFOUT,        0,                                                   &MD->iSTFcntr,            1,       MPI_INT,    &STATUS);
        MPI_File_write_at(fp_STFOUT,   OFFSETstf,                                                &dEQtime,                 1,       MPI_DOUBLE, &STATUS);  //Earthquake time
        MPI_File_write_at(fp_STFOUT,   OFFSETstf + sizeof(double),                               &fMagn,                   1,       MPI_FLOAT,  &STATUS);  //Earthquake magnitude
        MPI_File_write_at(fp_STFOUT,   OFFSETstf + sizeof(double) + sizeof(float),               &EQ->iCmbFPNum,           1,       MPI_INT,    &STATUS);  //#of fault patches participating in EQ 
        //at header of each event, I need to write the length of each STF...
    }
    OFFSETstf += sizeof(int) + sizeof(float) + sizeof(double);
    //--------------------------------------------
    MPI_File_write_at(fp_STFOUT,   OFFSETstf +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_ActPtchID,EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETstf += EQ->iCmbFPNum*sizeof(int);
    MPI_File_write_at(fp_STFOUT,   OFFSETstf +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_t0ofPtch, EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETstf += EQ->iCmbFPNum*sizeof(int);
    MPI_File_write_at(fp_STFOUT,   OFFSETstf +EQ->ivR_WrtStrtPos[MD->iRANK]*sizeof(int),  EQ->ivL_tDurPtch, EQ->iActFPNum, MPI_INT,   &STATUS);
    OFFSETstf += EQ->iCmbFPNum*sizeof(int);
    //--------------------------------------------
    iStartPos  = OFFSETstf + EQ->ivR_StfStrtPos[MD->iRANK];
    for (i = 0; i < EQ->iActFPNum; i++)
    {   iLPatchPos = (EQ->ivL_ActPtchID[i] - MD->ivF_START[MD->iRANK]);
        for (j = 0; j < EQ->ivL_tDurPtch[i]; j++)
        {   fNextVal = gsl_matrix_float_get(EQ->fmSTF_slip,     iLPatchPos, EQ->ivL_t0ofPtch[i]+j);             gsl_vector_float_set(STF1, j, fNextVal);
            fNextVal = gsl_matrix_float_get(EQ->fmSTF_strength, iLPatchPos, EQ->ivL_t0ofPtch[i]+j);             gsl_vector_float_set(STF2, j, fNextVal);
        }
        MPI_File_write_at(fp_STFOUT,   iStartPos + TempOffs , STF1->data,  EQ->ivL_tDurPtch[i], MPI_FLOAT,   &STATUS);
        TempOffs  += EQ->ivL_tDurPtch[i]*sizeof(float);
        MPI_File_write_at(fp_STFOUT,   iStartPos + TempOffs , STF2->data,  EQ->ivL_tDurPtch[i], MPI_FLOAT,   &STATUS);
        TempOffs  += EQ->ivL_tDurPtch[i]*sizeof(float);
    }   
    OFFSETstf += EQ->iCombSTFoffs;
    //--------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    return OFFSETstf;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void ResetFaultParameter(struct MDstruct *MD, struct TRstruct *TR, gsl_rng *fRandN)
{   int   i;
    float fTemp0, fTemp1, fTemp2, fTemp3, fTemp4, fTemp6;
    //------------------------
    for (i = 0; i < MD->ivF_OFFSET[MD->iRANK]; i++) 
    {   //-----------------------------------------
        if (TR->ivFL_Activated[i] == 1)     
        {   if (MD->iChgBtwEQs == 1)
            {   fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
                fTemp1                     = gsl_vector_float_get(TR->fvFL_RefDcVal, i)    + gsl_vector_float_get(TR->fvFL_RefDcVal_vari,   i)*fTemp2;
                gsl_vector_float_set(TR->fvFL_CurDcVal,i, fTemp1);
                //------------------------
                fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
                fTemp1                     = gsl_vector_float_get(TR->fvFL_RefStaFric, i)  + gsl_vector_float_get(TR->fvFL_RefStaFric_vari, i)*fTemp2; 
                gsl_vector_float_set(TR->fvFL_StaFric, i, fTemp1);     
                gsl_vector_float_set(TR->fvFL_CurFric, i, fTemp1);
                //------------------------
                fTemp2                     = (float)(gsl_rng_uniform(fRandN) *2.0 -1.0); //is a random number between -1 and 1
                fTemp3                     = gsl_vector_float_get(TR->fvFL_RefDynFric, i)  + gsl_vector_float_get(TR->fvFL_RefDynFric_vari, i)*fTemp2; 
                gsl_vector_float_set(TR->fvFL_DynFric, i, fTemp3);
                fTemp4                     = fTemp3 - (OVERSHOOTFRAC-1.0) *fabs(fTemp1 -fTemp3); //this friction is the "arrest friction" i.e., the dynamic overshoot friction value (the "lowest posssible")     
                gsl_vector_float_set(TR->fvFL_ArrFric, i, fTemp4);
        
                fTemp6                     = (fTemp3 - fTemp4)*-1.0*gsl_vector_float_get(TR->fvFL_RefNrmStrs,i); //this is "dynamic fric" minus "arrest fric" multiplied with normal stress => gives amount of excess above arrest fric to continue sliding => corresponds to "dyn fric" level
                gsl_vector_float_set(TR->fvFL_OverShotStress, i, fTemp6);
                //------------------------ 
                fTemp2 = (gsl_vector_float_get(TR->fvFL_DynFric,i) - gsl_vector_float_get(TR->fvFL_StaFric,i)) *-1.0*gsl_vector_float_get(TR->fvFL_RefNrmStrs,i); //gives a negative shear stress b/c of (dynF - staF)
                fTemp3 = fTemp2/gsl_vector_float_get(TR->fvFL_MeanSelfStiff,i); // b/c "self" has negative sign, the two negatives give a positive slip amount... (for weakening case, when dyn < stat fric) 
                if (fTemp3 > gsl_vector_float_get(TR->fvFL_CurDcVal,i))     {   TR->ivFL_StabT[i] = 1;      }
                else    
                {   if (fTemp2 < 0.0)                                       {   TR->ivFL_StabT[i] = 2;      }
                    else                                                    {   TR->ivFL_StabT[i] = 3;      }      
                }
            }       
        }
        TR->ivFL_Ptch_t0[i]        = 0;     
        TR->ivFL_Ptch_tDur[i]      = 0;     
        TR->ivFL_Activated[i]      = 0; 
        //-----------------------------------------
        gsl_vector_float_set(TR->fvFL_CurFric,  i, gsl_vector_float_get(TR->fvFL_StaFric,    i));   
        gsl_vector_float_set(TR->fvFL_CurStrsN, i, gsl_vector_float_get(TR->fvFL_RefNrmStrs, i));   
        //-----------------------------------------
        fTemp0 = sqrtf(gsl_vector_float_get(TR->fvFL_CurStrsH,i)*gsl_vector_float_get(TR->fvFL_CurStrsH,i) + gsl_vector_float_get(TR->fvFL_CurStrsV,i)*gsl_vector_float_get(TR->fvFL_CurStrsV,i) );
        fTemp1 = gsl_vector_float_get(TR->fvFL_CurFric,i)*-1.0*gsl_vector_float_get(TR->fvFL_CurStrsN,i);
        //-----------------------------------------
        if (MD->iUsePSeis == 1)
        {   if (fTemp0 > fTemp1)
            {   fTemp1 = fTemp0/(-1.0*gsl_vector_float_get(TR->fvFL_CurStrsN,i));
                gsl_vector_float_set(TR->fvFL_CurFric, i, fTemp1);
                        
                fTemp2 = fTemp1 - gsl_vector_float_get(TR->fvFL_StaFric, i);
                gsl_vector_float_set(TR->fvFL_PSeis_T0_F,i, fTemp2);
        }   }   
    }
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
int IntPow(int Base, int Exp)
{   int i, p;
    p = 1;
    for (i = 1; i<= Exp; i++)   {   p = p * Base;           }
    return p;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetVertices(struct VTstruct *VT,struct TRstruct *TR, int iGlobPos, int ForBoundary, float fP1[3], float fP2[3], float fP3[3])
{   
    if (ForBoundary == 0)
    {   fP1[0]   = VT->fvFG_PosE_temp[TR->ivFG_V1_temp[iGlobPos]];          fP1[1]   = VT->fvFG_PosN_temp[TR->ivFG_V1_temp[iGlobPos]];          fP1[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V1_temp[iGlobPos]];
        fP2[0]   = VT->fvFG_PosE_temp[TR->ivFG_V2_temp[iGlobPos]];          fP2[1]   = VT->fvFG_PosN_temp[TR->ivFG_V2_temp[iGlobPos]];          fP2[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V2_temp[iGlobPos]];
        fP3[0]   = VT->fvFG_PosE_temp[TR->ivFG_V3_temp[iGlobPos]];          fP3[1]   = VT->fvFG_PosN_temp[TR->ivFG_V3_temp[iGlobPos]];          fP3[2]   = VT->fvFG_PosZ_temp[TR->ivFG_V3_temp[iGlobPos]];
    }
    else
    {   fP1[0]   = VT->fvBG_PosE_temp[TR->ivBG_V1_temp[iGlobPos]];          fP1[1]   = VT->fvBG_PosN_temp[TR->ivBG_V1_temp[iGlobPos]];          fP1[2]   = VT->fvBG_PosZ_temp[TR->ivBG_V1_temp[iGlobPos]];
        fP2[0]   = VT->fvBG_PosE_temp[TR->ivBG_V2_temp[iGlobPos]];          fP2[1]   = VT->fvBG_PosN_temp[TR->ivBG_V2_temp[iGlobPos]];          fP2[2]   = VT->fvBG_PosZ_temp[TR->ivBG_V2_temp[iGlobPos]];
        fP3[0]   = VT->fvBG_PosE_temp[TR->ivBG_V3_temp[iGlobPos]];          fP3[1]   = VT->fvBG_PosN_temp[TR->ivBG_V3_temp[iGlobPos]];          fP3[2]   = VT->fvBG_PosZ_temp[TR->ivBG_V3_temp[iGlobPos]];  
    }
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void GetVector(float fP1[3], float fP2[3], float fP1P2[3])
{   fP1P2[0] = fP2[0] - fP1[0];                                         fP1P2[1] = fP2[1] - fP1[1];                                         fP1P2[2] = fP2[2] - fP1[2];  
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void NormalizeVector(float fVect[3])
{   float fVectLength,      fV[3];
    fV[0] = fVect[0];
    fV[1] = fVect[1];
    fV[2] = fVect[2];
    fVectLength = sqrtf(fV[0]*fV[0] +fV[1]*fV[1] +fV[2]*fV[2]);
    fVect[0]    = fV[0]/fVectLength;                                    fVect[1] = fV[1]/fVectLength;                                   fVect[2] = fV[2]/fVectLength;
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
void CrossProduct(float fA[3],  float fB[3],  float fvNrm[3])
{   fvNrm[0] = fA[1]*fB[2] - fA[2]*fB[1];                               fvNrm[1] = fA[2]*fB[0] - fA[0]*fB[2];                               fvNrm[2] = fA[0]*fB[1] - fA[1]*fB[0];
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetDistance(float fT[3], float fL[3])
{   float TheDist;
    TheDist   = sqrtf( (fT[0]-fL[0])*(fT[0]-fL[0]) + (fT[1]-fL[1])*(fT[1]-fL[1]) + (fT[2]-fL[2])*(fT[2]-fL[2]) ); 
    return TheDist;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
float GetMinDist_8Pnts(float fP1[3], float fP1a[3], float fP1b[3], float fP1c[3], float fP2[3], float fP2a[3], float fP2b[3], float fP2c[3])
{   float MinDist,      TstDist;

    MinDist   = sqrtf( (fP1[0]- fP2[0])*(fP1[0]- fP2[0]) + (fP1[1] -fP2[1])*(fP1[1] -fP2[1]) + (fP1[2] -fP2[2])*(fP1[2] -fP2[2]) ); 
    TstDist   = sqrtf( (fP1[0]-fP2a[0])*(fP1[0]-fP2a[0]) + (fP1[1]-fP2a[1])*(fP1[1]-fP2a[1]) + (fP1[2]-fP2a[2])*(fP1[2]-fP2a[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1[0]-fP2b[0])*(fP1[0]-fP2b[0]) + (fP1[1]-fP2b[1])*(fP1[1]-fP2b[1]) + (fP1[2]-fP2b[2])*(fP1[2]-fP2b[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1[0]-fP2c[0])*(fP1[0]-fP2c[0]) + (fP1[1]-fP2c[1])*(fP1[1]-fP2c[1]) + (fP1[2]-fP2c[2])*(fP1[2]-fP2c[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;

    TstDist   = sqrtf( (fP1a[0]- fP2[0])*(fP1a[0]- fP2[0]) + (fP1a[1] -fP2[1])*(fP1a[1] -fP2[1]) + (fP1a[2] -fP2[2])*(fP1a[2] -fP2[2]) ); 
    TstDist   = sqrtf( (fP1a[0]-fP2a[0])*(fP1a[0]-fP2a[0]) + (fP1a[1]-fP2a[1])*(fP1a[1]-fP2a[1]) + (fP1a[2]-fP2a[2])*(fP1a[2]-fP2a[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1a[0]-fP2b[0])*(fP1a[0]-fP2b[0]) + (fP1a[1]-fP2b[1])*(fP1a[1]-fP2b[1]) + (fP1a[2]-fP2b[2])*(fP1a[2]-fP2b[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1a[0]-fP2c[0])*(fP1a[0]-fP2c[0]) + (fP1a[1]-fP2c[1])*(fP1a[1]-fP2c[1]) + (fP1a[2]-fP2c[2])*(fP1a[2]-fP2c[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;

    TstDist   = sqrtf( (fP1b[0]- fP2[0])*(fP1b[0]- fP2[0]) + (fP1b[1] -fP2[1])*(fP1b[1] -fP2[1]) + (fP1b[2] -fP2[2])*(fP1b[2] -fP2[2]) ); 
    TstDist   = sqrtf( (fP1b[0]-fP2a[0])*(fP1b[0]-fP2a[0]) + (fP1b[1]-fP2a[1])*(fP1b[1]-fP2a[1]) + (fP1b[2]-fP2a[2])*(fP1b[2]-fP2a[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1b[0]-fP2b[0])*(fP1b[0]-fP2b[0]) + (fP1b[1]-fP2b[1])*(fP1b[1]-fP2b[1]) + (fP1b[2]-fP2b[2])*(fP1b[2]-fP2b[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1b[0]-fP2c[0])*(fP1b[0]-fP2c[0]) + (fP1b[1]-fP2c[1])*(fP1b[1]-fP2c[1]) + (fP1b[2]-fP2c[2])*(fP1b[2]-fP2c[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;

    TstDist   = sqrtf( (fP1c[0]- fP2[0])*(fP1c[0]- fP2[0]) + (fP1c[1] -fP2[1])*(fP1c[1] -fP2[1]) + (fP1c[2] -fP2[2])*(fP1c[2] -fP2[2]) ); 
    TstDist   = sqrtf( (fP1c[0]-fP2a[0])*(fP1c[0]-fP2a[0]) + (fP1c[1]-fP2a[1])*(fP1c[1]-fP2a[1]) + (fP1c[2]-fP2a[2])*(fP1c[2]-fP2a[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1c[0]-fP2b[0])*(fP1c[0]-fP2b[0]) + (fP1c[1]-fP2b[1])*(fP1c[1]-fP2b[1]) + (fP1c[2]-fP2b[2])*(fP1c[2]-fP2b[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;
    TstDist   = sqrtf( (fP1c[0]-fP2c[0])*(fP1c[0]-fP2c[0]) + (fP1c[1]-fP2c[1])*(fP1c[1]-fP2c[1]) + (fP1c[2]-fP2c[2])*(fP1c[2]-fP2c[2]) );     MinDist = (MinDist < TstDist) ? MinDist : TstDist;

    return MinDist;
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
    
 //   if (fvNrm[2] < 0.0)     {   fvNrm[0]       *= -1.0;                fvNrm[1]       *= -1.0;                 fvNrm[2]       *= -1.0;      }
    //-----------------------------------------------  
    fvStk[0]      = feZ[1]*fvNrm[2] - feZ[2]*fvNrm[1];
    fvStk[1]      = feZ[2]*fvNrm[0] - feZ[0]*fvNrm[2];
    fvStk[2]      = feZ[0]*fvNrm[1] - feZ[1]*fvNrm[0];
    // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction) 
    fTempfloat    = sqrtf(fvStk[0]*fvStk[0] +fvStk[1]*fvStk[1] +fvStk[2]*fvStk[2]);
    if (fabs(fTempfloat) < FLT_EPSILON)
    {   fvStk[0]    = 0.0;                      fvStk[1] = 1.0;                         fvStk[2] = 0.0;                     } 
    else
    {   fvStk[0]    = fvStk[0]/fTempfloat;      fvStk[1] = fvStk[1]/fTempfloat;         fvStk[2] = fvStk[2]/fTempfloat;     }
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
void RotateTensor(float fsig_new[6], const float fsig[6], const float fvNrm[3], const float fvStk[3], const float fvDip[3], int iRotDir)
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
    fsig_new[0] = fTxx2;        fsig_new[1] = fTxy2;      fsig_new[2] = fTxz2; //normal, strike, dip direction
    fsig_new[3] = fTyy2;        fsig_new[4] = fTyz2;      fsig_new[5] = fTzz2;
 
    return;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//