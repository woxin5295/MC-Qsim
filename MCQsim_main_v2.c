/* inialize the libraries*/
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>
# include <time.h>
# include <mpi.h>
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
float              ran0_inmain(long *idum);
extern void LoadInputParameter(char **argv, const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int iRealizNum,const int iUsdGrid, const int iFltSegmNum, const int iFltPtchNum, const int iBndPtchNum, const int iFltVertNum, const int iBndVertNum, const float *fSD_StrssRate, const float *fSD_SlipRate, const float *fSD_SlipRake,  int *iMD_ChgFricBtwEQs, float *fMD_AddNrmStrss, float *fMD_Vp, float *fMD_Vs, float *fMD_Poisson, float *fMD_Lambda, float *fMD_ShearMod, float *fMD_MedDense,  int *iSD_FricLawUSED,  float *fSD_RefStatFric, float *fSD_RefStatFr_vari, float *fSD_RefDynFric, float *fSD_RefDynFr_vari, float *fSD_CritSlipDist, float *fSD_CritSlipD_vari, int *iTDg_V1,  int *iTDg_V2,  int *iTDg_V3,   int *iTDl_SegID, float *fVDg_Epos, float *fVDg_Npos, float *fVDg_Zpos,float *fTDg_CentEpos, float *fTDg_CentNpos, float *fTDg_CentZpos,float *fTDl_StatFric, float *fTDl_DynFric,int   *iTDl_StabType, float *fTDl_RefStrssRateStk, float *fTDl_RefStrssRateDip, float *fTDl_SlipRate, float *fTDl_SlipRake, float *fTDl_CurrFric);
extern void    DefineMoreParas(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const  int iFltPtchNum, const  int iBndPtchNum, const  int *iTDg_V1, const  int *iTDg_V2, const  int *iTDg_V3, const  int *iTDl_SegID, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fSD_CritSlipDist, const float *fSD_CritSlipD_vari, const float *fMD_Vp, const float *fMD_VpVsRatio, const float *fMD_MedDense, const float *fMD_AddNrmStrss, const float fMD_g, const float fdeltTincr, const float *fRandVector, const  int iRandNumber,  int *iRandPos, int *iTDlg_TravTimesP,  int *iTDlg_TravTimesS,  int *iMD_GlobTTmax, float *fTDl_Area, float *fTDl_RefNormStrss, float *fTDl_Curr_DcVal, float *fTDlg_LocSrcRcv_H, float *fTDlg_LocSrcRcv_V, float *fTDlg_LocSrcRcv_N);

extern void     Build_K_Matrix(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const  int iFltPtchNum, const  int iBndPtchNum, const float fMDg_UnitSlip, const float *fMDg_ShearMod, const float *fMDg_Lambda, const  int *iTDg_V1, const  int *iTDg_V2, const  int *iTDg_V3, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fTDl_Curr_DcVal, const float *fTDl_StatFric, const float *fTDl_DynFric, const float *fTDl_RefNormStrss, const float *fTDl_Area, float *fK_FF_SS, float *fK_FF_SD, float *fK_FF_SO, float *fK_FF_DS, float *fK_FF_DD, float *fK_FF_DO, float *fK_FF_OS, float *fK_FF_OD, float *fK_FF_OO, float *fK_FB_SS, float *fK_FB_SD, float *fK_FB_SO, float *fK_FB_DS, float *fK_FB_DD, float *fK_FB_DO, float *fK_BF_SS, float *fK_BF_SD, float *fK_BF_SO, float *fK_BF_DS, float *fK_BF_DD, float *fK_BF_DO, float *fK_BF_OS, float *fK_BF_OD, float *fK_BF_OO, float *fK_BB_SS, float *fK_BB_SD, float *fK_BB_SO, float *fK_BB_DS, float *fK_BB_DD, float *fK_BB_DO, float *fK_BB_OS, float *fK_BB_OD, float *fK_BB_OO, float *fKl_BB_SS, float *fKl_BB_DD, float *fKl_BB_OO,  int *iTDl_SelfLoc_F,  int *iTDl_SelfLoc_B,  int *iTDl_StabType);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int main(int argc, char **argv)
{   if (argc != 7)  {   fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -n 4 ./Run_MCQsim InputName RealizationNumber GridNumber WhichDeltT2Use UsePostSeis EQrecordLength\n");        exit(10);           }
    /*-------------------------------------------------------------------*/
    int       i,              j,                  k,               stfpos;
    int       iRANK,          iSIZE,              iPrintAll;
    /*-------------------------------------------------------------------*/
    iPrintAll = 1;
    //printf("The minimum value of INT = %d\n", INT_MIN);
    //printf("The maximum value of INT = %d\n", INT_MAX);
    MPI_Init(&argc, &argv);               /* start MPI processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &iRANK);/* current proc. id    */
    MPI_Comm_size(MPI_COMM_WORLD, &iSIZE);/* # of processes      */ 

    MPI_Status status;
    MPI_Offset offset1; /*this stuff here is for file output*/
    MPI_Offset offset2; /*this stuff here is for file output*/
    MPI_File fp_MPIout1;
    /*-------------------------------------------------------------------*/
    int    iOFFSET_F[iSIZE];       for (i = 0; i < iSIZE; i++)                     {       iOFFSET_F[i]    = 0;                              }
    int    iSTARTPOS_F[iSIZE];     for (i = 0; i < iSIZE; i++)                     {       iSTARTPOS_F[i]  = 0;                              }
    int    iOFFSET_B[iSIZE];       for (i = 0; i < iSIZE; i++)                     {       iOFFSET_B[i]    = 0;                              }
    int    iSTARTPOS_B[iSIZE];     for (i = 0; i < iSIZE; i++)                     {       iSTARTPOS_B[i]  = 0;                              }
    /*-------------------------------------------------------------------*/
    int    iUsdGrid,           iDeltT2Use,             iUsePostSeismic,        iRealizNum;
    float  fRecLgth;
    char   *retch;
    /*----------------------------*/
    sscanf(argv[2],"%i",  &iRealizNum);
    sscanf(argv[3],"%i",  &iUsdGrid);                           sscanf(argv[4],"%i",  &iDeltT2Use);
    sscanf(argv[5],"%i",  &iUsePostSeismic);                    sscanf(argv[6],"%f",  &fRecLgth);
    /*
      argv[1] is the file name -but only the "core" of the filename => like "AqabaFaults" for all models like "AqabaFaults_1_Roughn.dat" etc...
      argv[2] is the realization number (if I created 100 random models... -> which one to use..)
      !if [2] is zero, then I'll use the planar model!
      
      argv[3] grid resolution to use (if more than one was used and stored in the binary files => which one to use
      argv[4] switch to decide whether to store large-earthquake kinematic models
      argv[5] this defines the magnitude for those "large-earthquakes" => EQs larger than this will be stored...
      argv[6] whether to use deltT based on v_p (== 1) or to use inf deltT (instantaneous signal distribution,  == 0)
      argv[7] this is length of record (in years)
    */  
    /*-------------------------------------------------------------------*/
    int    iFltSegmNum,        iFltGridNum,        iFltPtchNum,                iFltVertNum; 
    int    iBndSegmNum,        iBndGridNum,        iBndPtchNum,                iBndVertNum;
    float  fMeanLegLgth,       fMeanBndLegLgth,    fDummy;
    float  *fSD_SlipRate,      *fSD_SlipRake,      *fSD_StrssRate;

    char   ctempVals[512],     cFile1_In[512],     cFile2_In[512],             cFile1_Out[512];   /* out1 is for the overall catalog; out2 is the stf for each patch that had slip => both are binaries!!!      */
    FILE   *fpIn;  
    /*----------------------------*/
    strcpy(cFile1_In,  argv[1]);                            strcat(cFile1_In,"_FLT.txt");
    strcpy(cFile2_In,  argv[1]);                            strcat(cFile2_In,"_BND.txt");
    strcpy(cFile1_Out, argv[1]);                            strcat(cFile1_Out,"_");             strcat(cFile1_Out,argv[2]);         strcat(cFile1_Out,"_Catalog.dat");
    /*----------------------------------------------------------------------------------*/  
    if ((fpIn = fopen(cFile1_In,"r")) == NULL)          {   fprintf(stdout,"Error -cant open %s file. in Main function \n",cFile1_In);      exit(10);     }
    retch =fgets(ctempVals, 512, fpIn);                     sscanf(ctempVals,"%*s %d",  &iFltSegmNum); 
    retch =fgets(ctempVals, 512, fpIn);                     sscanf(ctempVals,"%*s %d",  &iFltGridNum);    
    /*----------------------------*/
    fSD_SlipRate  = (float *) calloc(iFltSegmNum, sizeof(float));
    fSD_SlipRake  = (float *) calloc(iFltSegmNum, sizeof(float));
    fSD_StrssRate = (float *) calloc(iFltSegmNum, sizeof(float));
    /*----------------------------*/
    iUsdGrid = (iUsdGrid  < iFltGridNum) ? iUsdGrid : iFltGridNum; /*if i picked a grid-to-use that is larger than number of grids that were generated, then use "last grid that is generated"*/
    /*----------------------------*/
    for (i = 0; i < iFltGridNum; i++)
    {   if ((i+1) == iUsdGrid)
        {   retch =fgets(ctempVals, 512, fpIn);                 sscanf(ctempVals,"%*s %d  %d", &iFltPtchNum,  &iFltVertNum);
            retch =fgets(ctempVals, 512, fpIn);                 sscanf(ctempVals,"%*s %e  %e", &fMeanLegLgth, &fDummy);
        }
        else
        {   retch =fgets(ctempVals, 512, fpIn);
            retch =fgets(ctempVals, 512, fpIn);
    }   }
    fMeanLegLgth *= 1.0E+3; /* now it is in meters */
    /*----------------------------*/
    for (i = 0; i < iFltSegmNum; i++)
    {   for (j = 0; j < 10; j++)                         {   retch =fgets(ctempVals, 512, fpIn);                                                                        }
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %e", &fSD_SlipRate[i]);              fSD_SlipRate[i]  /=1.0E+3; /* now in m/yr*/
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %e", &fSD_SlipRake[i]);              fSD_SlipRake[i]  *= M_PI/180.0; /* now in radiants*/
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %e", &fSD_StrssRate[i]);             fSD_StrssRate[i] *=1.0E+6; /* now in Pa/yr*/
        for (j = 0; j < 4; j++)                          {   retch =fgets(ctempVals, 512, fpIn);                                                                        }
        for (j = 0; j < iFltGridNum; j++)                {   retch =fgets(ctempVals, 512, fpIn);                        retch = fgets(ctempVals, 512, fpIn);            }
        if ((iPrintAll > 0) && (iRANK == 0))
        {   fprintf(stdout,"Seg# %d   SlipRate: %f     SlipRake: %f      StrssRate: %f\n",i, fSD_SlipRate[i],    fSD_SlipRake[i],      fSD_StrssRate[i]); 
        }
    }
    fclose(fpIn);
    /*----------------------------------------------------------------------------------*/  
    if ((fpIn = fopen(cFile2_In,"r")) == NULL)
    {   fprintf(stdout,"_BND.txt file was not found i.e., opened; continue without BC faults.... \n");
        iBndSegmNum = 0;                                     iBndGridNum = 0;
        iBndPtchNum = 0;                                     iBndVertNum = 0;
        fMeanBndLegLgth = 0.0;
    }
    else
    {   retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %d",  &iBndSegmNum);
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %d",  &iBndGridNum);
        if (iBndGridNum < 1)                             {   fprintf(stdout,"Error -boundary fault file was opened but BC faults are not meshed/gridded. \n");      exit(10);     }
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %d  %d", &iBndPtchNum,     &iBndVertNum);
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %e  %e", &fMeanBndLegLgth, &fDummy);
        
        fclose(fpIn);
    }
    fMeanBndLegLgth *= 1.0E+3; /* now it is in meters */
    
    /*----------------------------------------------------------------------------------*/  
    fMeanBndLegLgth = (fMeanBndLegLgth <=      0.0    ) ? fMeanLegLgth : fMeanBndLegLgth; /*set BndLegLength to MeanLegLgth if BndLegLength is equal or smaller than 0.0*/
    fMeanLegLgth    = (fMeanBndLegLgth >= fMeanLegLgth) ? fMeanLegLgth : fMeanBndLegLgth; /*set MeanLegLength to be min value of LegLength and BndLegLength*/
    /*----------------------------------------------------------------------------------*/  
    int    iBASEelem_F;          iBASEelem_F = ( int)(iFltPtchNum/iSIZE); 
    int    iADDelem_F;           iADDelem_F  = ( int)(iFltPtchNum%iSIZE);
    int    iBASEelem_B;          iBASEelem_B = ( int)(iBndPtchNum/iSIZE); 
    int    iADDelem_B;           iADDelem_B  = ( int)(iBndPtchNum%iSIZE);
    /*----------------------------*/
    for (i = 0; i < iSIZE;      i++)      {   iOFFSET_F[i]     = iBASEelem_F;                                         }
    for (i = 0; i < iADDelem_F; i++)      {   iOFFSET_F[i]    += 1;                                                   }
    for (i = 1; i < iSIZE;      i++)      {   iSTARTPOS_F[i]   = iSTARTPOS_F[i-1] + iOFFSET_F[i-1];                   }
    for (i = 0; i < iSIZE;      i++)      {   iOFFSET_B[i]     = iBASEelem_B;                                         }
    for (i = 0; i < iADDelem_B; i++)      {   iOFFSET_B[i]    += 1;                                                   }
    for (i = 1; i < iSIZE;      i++)      {   iSTARTPOS_B[i]   = iSTARTPOS_B[i-1] + iOFFSET_B[i-1];                   }
    /*----------------------------------------------------------------------------------*/  
    if (iPrintAll > 0)
    {   fprintf(stdout,"My Rank: %d  TriNum %d    BASEelem  %d   ADDelem  %d  STARTPOS  %d OFFSET   %d; start: %d     end: %d\n",iRANK, iFltPtchNum, iBASEelem_F, iADDelem_F, iSTARTPOS_F[iRANK], iOFFSET_F[iRANK], iSTARTPOS_F[iRANK],  (iSTARTPOS_F[iRANK]+iOFFSET_F[iRANK])  );  
    }

    int    iEQcounter      = 0,                 iMaxMRFlength  = 10000; 
    float  fMD_UnitSlip   = 1.0E-4*fMeanLegLgth; /* slip in meter => if meanleglength is 1000m, then a slip of .1m is used */
    float  fMD_g           = 9.81,              fViscosity     = 2.0E+15; /* in Pa seconds; value doesn't seem right, but need something that corresponds to patch stiffness => to give proper char. time scale for post-seismic decay*/ 
    float  fIntSeisTimeStp = 2.0/365.25,        fTimeStpInSecs = fIntSeisTimeStp*31536000.0; /* now the time is in seconds */
    /* stress in Pa ==> earth tides change stress by about 5000Pa (= 0.005MPa); atmospheric pressure is ~100kPa => 2000Pa == 2% of atmospheric pressure change */   
    int    iVectPos,                            iMaxSTFLength,                              STFcnt;
    int    iwritePos,                           ireadPos;
    float  fdeltTincr,                          fTimeYears,                                 fMD_CutStrss; 
    float fTempSlipStk,                         fTempSlipDip;
    float fTemp,                                fTemp1,                                     fTemp2;   
    float fTemp3,                               fStkSlip,                                   fDipSlip;            
    float fPostSeisIncrStk,                     fPostSeisIncrDip,                           fPostSeisIncrNrm;           
    float fUsedStkSlip,                         fUsedDipSlip,                               fUsedNrmSlip;        
    /*-------------------------------------------------------------------*/ 
    int    iRandNumber     = 1000000; 
    float  *fRandVector;               fRandVector              = (float *)  calloc(iRandNumber,  sizeof(float));
    int    *iRandPos;                  iRandPos                 = (int   *)  calloc(       1   ,  sizeof( int));  
    long   *iSeedVal;                  iSeedVal                 = (long *)   calloc(       1   ,  sizeof(long));  
    time_t  t;
    srand((unsigned) time(&t)+iRANK*10);  /* initialize the random number generator; will generate 1E+6 random numbers => should be enough*/
    iSeedVal[0] = rand();
    for (i = 0; i < iRandNumber; i++ )       {           fRandVector[i] = ran0_inmain(iSeedVal);       } /* using a generator as described in "Numerical Recipes, value between 0  and 1; page 279*/
    /*-------------------------------------------------------------------*/
    float  *fMD_AddNrmStrss;           fMD_AddNrmStrss          = (float *)  calloc(1, sizeof(float));
    float  *fMD_VpVsRatio;             fMD_VpVsRatio            = (float *)  calloc(1, sizeof(float));
    float  *fMD_ShearMod;              fMD_ShearMod             = (float *)  calloc(1, sizeof(float));
    float  *fMD_Vp;                    fMD_Vp                   = (float *)  calloc(1, sizeof(float));
    float  *fMD_Vs;                    fMD_Vs                   = (float *)  calloc(1, sizeof(float));
    float  *fMD_Poisson;               fMD_Poisson              = (float *)  calloc(1, sizeof(float));
    float  *fMD_Lambda;                fMD_Lambda               = (float *)  calloc(1, sizeof(float));
    float  *fMD_MedDense;              fMD_MedDense             = (float *)  calloc(1, sizeof(float));
    //float  *fMD_CritSlipVelo;          fMD_CritSlipVelo         = (float *)  calloc(1, sizeof(float));    
    int    *iMD_ChgFricBtwEQs;         iMD_ChgFricBtwEQs        = (int   *)  calloc(1, sizeof( int));  
    int    *iMD_GlobTTmax;             iMD_GlobTTmax            = (int   *)  calloc(1, sizeof( int));    
    /*-------------------------------------------------------------------*/ 
    int    *iSD_FricLawUSED;           iSD_FricLawUSED          = (int   *)  calloc(iFltSegmNum,   sizeof( int));
    float  *fSD_RefStatFric;           fSD_RefStatFric          = (float *)  calloc(iFltSegmNum,   sizeof(float));
    float  *fSD_RefStatFr_vari;        fSD_RefStatFr_vari       = (float *)  calloc(iFltSegmNum,   sizeof(float));
    float  *fSD_RefDynFric;            fSD_RefDynFric           = (float *)  calloc(iFltSegmNum,   sizeof(float));
    float  *fSD_RefDynFr_vari;         fSD_RefDynFr_vari        = (float *)  calloc(iFltSegmNum,   sizeof(float));
    float  *fSD_CritSlipDist;          fSD_CritSlipDist         = (float *)  calloc(iFltSegmNum,   sizeof(float));
    float  *fSD_CritSlipD_vari;        fSD_CritSlipD_vari       = (float *)  calloc(iFltSegmNum,   sizeof(float));
    /*-------------------------------------------------------------------*/ 
    int    *iTDl_StabType;             iTDl_StabType            = ( int *)   calloc(iOFFSET_F[iRANK],   sizeof( int));
    int    *iTDl_StabTypeChgd;         iTDl_StabTypeChgd        = ( int *)   calloc(iOFFSET_F[iRANK],   sizeof( int));
    int    *iTDl_SegID;                iTDl_SegID               = ( int *)   calloc(iOFFSET_F[iRANK],   sizeof( int));
    float  *fTDl_SlipRate;             fTDl_SlipRate            = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_SlipRake;             fTDl_SlipRake            = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_RefStrssRateStk;      fTDl_RefStrssRateStk     = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));/*Ref and CurStrss rate can be different b/c of stability type changes/effects...*/
    float  *fTDl_RefStrssRateDip;      fTDl_RefStrssRateDip     = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_CurStrssRateStk;      fTDl_CurStrssRateStk     = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_CurStrssRateDip;      fTDl_CurStrssRateDip     = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_StatFric;             fTDl_StatFric            = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_CurrFric;             fTDl_CurrFric            = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_DynFric;              fTDl_DynFric             = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssAt_t0_S_F;  fTDl_PSeisStrssAt_t0_S_F = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssAt_t0_N_F;  fTDl_PSeisStrssAt_t0_N_F = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssTStep_F;    fTDl_PSeisStrssTStep_F   = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssAt_t0_S_B;  fTDl_PSeisStrssAt_t0_S_B = (float *)  calloc(iOFFSET_B[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssAt_t0_N_B;  fTDl_PSeisStrssAt_t0_N_B = (float *)  calloc(iOFFSET_B[iRANK],   sizeof(float));
    float  *fTDl_PSeisStrssTStep_B;    fTDl_PSeisStrssTStep_B   = (float *)  calloc(iOFFSET_B[iRANK],   sizeof(float));
    float  *fTDl_RefNormStrss;         fTDl_RefNormStrss        = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_Area;                 fTDl_Area                = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    float  *fTDl_Curr_DcVal;           fTDl_Curr_DcVal          = (float *)  calloc(iOFFSET_F[iRANK],   sizeof(float));
    /*-------------------------------------------------------------------*/ 
    int    *iTempValue;                iTempValue               = ( int *)   calloc(       1,         sizeof( int));  
    /*-------------------------------------------------------------------*/ 
    float  *fTDg_CentEpos;             fTDg_CentEpos            = (float *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof(float));
    float  *fTDg_CentNpos;             fTDg_CentNpos            = (float *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof(float));
    float  *fTDg_CentZpos;             fTDg_CentZpos            = (float *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof(float));

    int    *iTDg_V1;                   iTDg_V1                  = (int   *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof( int)); /*those are "T" b/c they will be freed later on -> only used during model setup e.g., calculation of K-matrix*/   
    int    *iTDg_V2;                   iTDg_V2                  = (int   *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof( int));
    int    *iTDg_V3;                   iTDg_V3                  = (int   *)  calloc((iFltPtchNum + iBndPtchNum),  sizeof( int));
    /*-------------------------------------------------------------------*/ 
    float  *fVDg_Epos;                 fVDg_Epos                = (float *)  calloc((iFltVertNum + iBndVertNum),  sizeof(float));
    float  *fVDg_Npos;                 fVDg_Npos                = (float *)  calloc((iFltVertNum + iBndVertNum),  sizeof(float));
    float  *fVDg_Zpos;                 fVDg_Zpos                = (float *)  calloc((iFltVertNum + iBndVertNum),  sizeof(float));
    /*-------------------------------------------------------------------*/ 
    int    *iTDlg_TravTimesP;          iTDlg_TravTimesP         = (int   *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof( int));
    int    *iTDlg_TravTimesS;          iTDlg_TravTimesS         = (int   *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof( int));
    float  *fTDlg_LocSrcRcv_H;         fTDlg_LocSrcRcv_H        = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float));
    float  *fTDlg_LocSrcRcv_V;         fTDlg_LocSrcRcv_V        = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float)); 
    float  *fTDlg_LocSrcRcv_N;         fTDlg_LocSrcRcv_N        = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float));  
    /*-------------------------------------------------------------------*/ 
    int    *iTDl_SelfLoc_F;           iTDl_SelfLoc_F            = ( int *)   calloc(iOFFSET_F[iRANK], sizeof( int));     
    int    *iTDl_SelfLoc_B;           iTDl_SelfLoc_B            = ( int *)   calloc(iOFFSET_B[iRANK], sizeof( int));     
    /*-------------------------------------------------------------------*/ 
    float  *fK_FF_SS;                 fK_FF_SS                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));            
    float  *fK_FF_SD;                 fK_FF_SD                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));                    
    float  *fK_FF_SO;                 fK_FF_SO                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));            
    float  *fK_FF_DS;                 fK_FF_DS                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));       
    float  *fK_FF_DD;                 fK_FF_DD                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));                           
    float  *fK_FF_DO;                 fK_FF_DO                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));    
    float  *fK_FF_OS;                 fK_FF_OS                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));       
    float  *fK_FF_OD;                 fK_FF_OD                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));                           
    float  *fK_FF_OO;                 fK_FF_OO                  = (float *)  calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));    
    /*use full interaction for FF; reason is the transient signals that can cause "opening" => needs to be computed to change normal stresses*/
    /*-------------------------------------------------------------------*/ 
    float  *fK_FB_SS;                 fK_FB_SS                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));            
    float  *fK_FB_SD;                 fK_FB_SD                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));                    
    float  *fK_FB_SO;                 fK_FB_SO                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));            
    float  *fK_FB_DS;                 fK_FB_DS                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));       
    float  *fK_FB_DD;                 fK_FB_DD                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));                           
    float  *fK_FB_DO;                 fK_FB_DO                  = (float *)  calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));  
    /*this uses not the transient signal but the actual/final slip on fault patches => which is always just in-plane motion*/  
    /*-------------------------------------------------------------------*/ 
    float  *fK_BF_SS;                 fK_BF_SS                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));            
    float  *fK_BF_SD;                 fK_BF_SD                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));                    
    float  *fK_BF_SO;                 fK_BF_SO                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));   
    float  *fK_BF_DS;                 fK_BF_DS                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));       
    float  *fK_BF_DD;                 fK_BF_DD                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));                           
    float  *fK_BF_DO;                 fK_BF_DO                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));                        
    float  *fK_BF_OS;                 fK_BF_OS                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));   
    float  *fK_BF_OD;                 fK_BF_OD                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));                           
    float  *fK_BF_OO;                 fK_BF_OO                  = (float *)  calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));    
    /*this is also full b/c I allow the postseismic part to change/relax normal stress changes that were produced by previous earthquakes*/                                 
    /*-------------------------------------------------------------------*/     
    float  *fK_BB_SS;                 fK_BB_SS                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));            
    float  *fK_BB_SD;                 fK_BB_SD                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));    
    float  *fK_BB_SO;                 fK_BB_SO                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));                    
    float  *fK_BB_DS;                 fK_BB_DS                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));       
    float  *fK_BB_DD;                 fK_BB_DD                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));    
    float  *fK_BB_DO;                 fK_BB_DO                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));
    float  *fK_BB_OS;                 fK_BB_OS                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));       
    float  *fK_BB_OD;                 fK_BB_OD                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));    
    float  *fK_BB_OO;                 fK_BB_OO                  = (float *)  calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));        
    /*need full interaction here for computing the stressing rate in case that slip rates on fault along with boundary faults are defined.*/
    float  *fKl_BB_SS;                fKl_BB_SS                 = (float *)  calloc(iOFFSET_B[iRANK],            sizeof(float)); 
    float  *fKl_BB_DD;                fKl_BB_DD                 = (float *)  calloc(iOFFSET_B[iRANK],            sizeof(float)); 
    float  *fKl_BB_OO;                fKl_BB_OO                 = (float *)  calloc(iOFFSET_B[iRANK],            sizeof(float));                         
    /*the BB matrix can be deleted once the stressing rates are found; however, I still need the "self-stiffness" of the boundary patches to know how much slip to use for given stress amount in postseismic relaxation*/   
    int    *iTDlg_Pdone;              iTDlg_Pdone               = (int *)    calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(int));    
    int    *iTDlg_Sdone;              iTDlg_Sdone               = (int *)    calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(int));    
    /*these two tell me the index for source-receiver pair => how much of STF at source has been "seen" at receiver already... takes more memory but is worth it b/c it allows to make the STF shorter and ensures that it will not "blow up" when events are too long...*/
    float  *fTDg_StrssRteChgStk;      fTDg_StrssRteChgStk       = (float *)  calloc(iFltPtchNum,      sizeof(float)); 
    float  *fTDg_StrssRteChgDip;      fTDg_StrssRteChgDip       = (float *)  calloc(iFltPtchNum,      sizeof(float)); 
    float  *fTDl_CharRlxTime_F;       fTDl_CharRlxTime_F        = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float ));         
    float  *fTDl_CharRlxTime_B;       fTDl_CharRlxTime_B        = (float *)  calloc(iOFFSET_B[iRANK], sizeof(float )); 
    /*---------------------------------------*/                       
    int   *iEQ_EQongoing;              iEQ_EQongoing            = ( int *)   calloc(1 ,sizeof( int)); 
    int   *iEQ_ChangedStab;            iEQ_ChangedStab          = ( int *)   calloc(1, sizeof( int));       
    int   *iEQ_EQEndCntr;              iEQ_EQEndCntr            = ( int *)   calloc(1 ,sizeof( int));    
    int   *iEQ_ActPtchNum;             iEQ_ActPtchNum           = ( int *)   calloc(1, sizeof( int));  
    int   *iEQ_CmbPtchNum;             iEQ_CmbPtchNum           = ( int *)   calloc(1, sizeof( int));  
    int   *iEQ_MRFlength;              iEQ_MRFlength            = ( int *)   calloc(1, sizeof( int));  
    int   *iEQ_TotalRuptT;             iEQ_TotalRuptT           = ( int *)   calloc(1, sizeof( int));    
    float *fEQ_SeisPot;                fEQ_SeisPot              = (float *)  calloc(1, sizeof(float));
    int   *iEQ_WrtStartPos;            iEQ_WrtStartPos          = ( int *)   calloc(iSIZE, sizeof( int));
    
    int   *iEQl_ActPtchID;             iEQl_ActPtchID           = ( int *)   calloc(iOFFSET_F[iRANK],  sizeof( int));                         
    int   *iEQl_t0ofPtch;              iEQl_t0ofPtch            = ( int *)   calloc(iOFFSET_F[iRANK],  sizeof( int));    
    int   *iEQl_StabType;              iEQl_StabType            = ( int *)   calloc(iOFFSET_F[iRANK],  sizeof( int));    

    float *fEQl_SlipHofPtch;           fEQl_SlipHofPtch         = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));            
    float *fEQl_SlipVofPtch;           fEQl_SlipVofPtch         = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));
    float *fEQl_DeltTofPtch;           fEQl_DeltTofPtch         = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));
    float *fEQ_MRFvals;                fEQ_MRFvals              = (float *)  calloc(iMaxMRFlength,    sizeof(float));
    
    int   *iTDl_WasActivated;          iTDl_WasActivated        = (int *)    calloc(iOFFSET_F[iRANK], sizeof( int)); 
    int   *iTDl_t0;                    iTDl_t0                  = (int *)    calloc(iOFFSET_F[iRANK], sizeof( int));                       
    float *fTDl_StrssB4_H;             fTDl_StrssB4_H           = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));                 
    float *fTDl_StrssB4_V;             fTDl_StrssB4_V           = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));     
    float *fTDl_StrssB4_N;             fTDl_StrssB4_N           = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float)); 
    float *fTDl_CurStrss_H;            fTDl_CurStrss_H          = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));                 
    float *fTDl_CurStrss_V;            fTDl_CurStrss_V          = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));     
    float *fTDl_CurStrss_N;            fTDl_CurStrss_N          = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));     
    float *fTDg_DeltStrssH;            fTDg_DeltStrssH          = (float *)  calloc(iFltPtchNum,      sizeof(float));                 
    float *fTDg_DeltStrssV;            fTDg_DeltStrssV          = (float *)  calloc(iFltPtchNum,      sizeof(float));     
    float *fTDg_DeltStrssN;            fTDg_DeltStrssN          = (float *)  calloc(iFltPtchNum,      sizeof(float));               
    float *fTDg_TempBnd_H;             fTDg_TempBnd_H           = (float *)  calloc(iBndPtchNum,      sizeof(float)); 
    float *fTDg_TempBnd_V;             fTDg_TempBnd_V           = (float *)  calloc(iBndPtchNum,      sizeof(float)); 
    float *fTDg_TempBnd_N;             fTDg_TempBnd_N           = (float *)  calloc(iBndPtchNum,      sizeof(float)); 
    
    float *fTDl_StrssOnBnd_H;          fTDl_StrssOnBnd_H        = (float *)  calloc(iOFFSET_B[iRANK], sizeof(float));
    float *fTDl_StrssOnBnd_V;          fTDl_StrssOnBnd_V        = (float *)  calloc(iOFFSET_B[iRANK], sizeof(float));
    float *fTDl_StrssOnBnd_N;          fTDl_StrssOnBnd_N        = (float *)  calloc(iOFFSET_B[iRANK], sizeof(float));
    
    float *fTDl_EventSlipH;            fTDl_EventSlipH          = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float));     
    float *fTDl_EventSlipV;            fTDl_EventSlipV          = (float *)  calloc(iOFFSET_F[iRANK], sizeof(float)); 
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/  
    LoadInputParameter(argv, iRANK, iSTARTPOS_F, iOFFSET_F, iRealizNum, iUsdGrid, iFltSegmNum, iFltPtchNum, iBndPtchNum, iFltVertNum, iBndVertNum, fSD_StrssRate, fSD_SlipRate, fSD_SlipRake, iMD_ChgFricBtwEQs, fMD_AddNrmStrss, fMD_Vp, fMD_Vs, fMD_Poisson, fMD_Lambda, fMD_ShearMod, fMD_MedDense, iSD_FricLawUSED, fSD_RefStatFric, fSD_RefStatFr_vari, fSD_RefDynFric, fSD_RefDynFr_vari, fSD_CritSlipDist, fSD_CritSlipD_vari, iTDg_V1, iTDg_V2, iTDg_V3, iTDl_SegID, fVDg_Epos, fVDg_Npos, fVDg_Zpos, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fTDl_StatFric, fTDl_DynFric, iTDl_StabType, fTDl_RefStrssRateStk, fTDl_RefStrssRateDip, fTDl_SlipRate, fTDl_SlipRake, fTDl_CurrFric);
 
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/  
    fMD_VpVsRatio[0] = fMD_Vp[0]/fMD_Vs[0];
    /*-------------------------------------------------------------------*/    
    if (iDeltT2Use == 1)    {       fdeltTincr = 2.0*fMeanLegLgth/fMD_Vp[0];                 fdeltTincr = floorf(fdeltTincr*100.0)/100.0;           }
    else                    {       fdeltTincr = 1.0E+09;                                                                                           }
    if ((iPrintAll > 0) && (iRANK == 0))        
    {   fprintf(stdout,"delta T (coseismic) %f   fMeanLegLgth %f  RecLength %f\n",fdeltTincr,fMeanLegLgth,fRecLgth);  
        fTemp  = (0.5*fMD_ShearMod[0])/(fMeanLegLgth/3.0); //this is in conceptually based on stuff presented in Dieterich, 1992; but the coefficients are picked to fit data (for triangles) better
        fDummy = fViscosity/fTemp /31536000.0;
        fprintf(stdout,"approx. char. relaxation time in years (to get shear to 1/e (~35percent) of initial shear) %f \n",fDummy);   
    }
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/  
    DefineMoreParas(iRANK, iSTARTPOS_F, iOFFSET_F, iSTARTPOS_B, iOFFSET_B, iFltPtchNum, iBndPtchNum, iTDg_V1, iTDg_V2, iTDg_V3, iTDl_SegID, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fVDg_Epos, fVDg_Npos, fVDg_Zpos, fSD_CritSlipDist, fSD_CritSlipD_vari, fMD_Vp, fMD_VpVsRatio, fMD_MedDense, fMD_AddNrmStrss, fMD_g, fdeltTincr, fRandVector, iRandNumber, iRandPos, iTDlg_TravTimesP, iTDlg_TravTimesS, iMD_GlobTTmax, fTDl_Area, fTDl_RefNormStrss, fTDl_Curr_DcVal, fTDlg_LocSrcRcv_H, fTDlg_LocSrcRcv_V, fTDlg_LocSrcRcv_N);
    
//    if (iPrintAll > 0)         
//    {   for (i = 0; i < iOFFSET_F[iRANK]; i++)           
//        {   for (j = 0; j < iFltPtchNum;  j++) 
//            {   iVectPos  = i*iFltPtchNum + j;       
//                fprintf(stdout,"local coords: %d   %d    %f    %f    %f \n",i, j, fTDlg_LocSrcRcv_N[iVectPos] , fTDlg_LocSrcRcv_H[iVectPos] , fTDlg_LocSrcRcv_V[iVectPos] );
//    }   }   }
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/     
    MPI_Allreduce(MPI_IN_PLACE, iMD_GlobTTmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    
    iMaxSTFLength  = iMD_GlobTTmax[0] +1;
    if (iRANK == 0)            {         fprintf(stdout,"MaxSTFlength =  %d\n",  iMaxSTFLength);    }
    

    float *fTDlg_STF_H;        fTDlg_STF_H       = (float *)  calloc(iOFFSET_F[iRANK]*iMaxSTFLength, sizeof(float)); 
    float *fTDlg_STF_V;        fTDlg_STF_V       = (float *)  calloc(iOFFSET_F[iRANK]*iMaxSTFLength, sizeof(float));     
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/ 
    if (iRANK == 0)            {         fprintf(stdout,"Building K-matrix....\n");                 }
    
    Build_K_Matrix(iRANK, iSTARTPOS_F, iOFFSET_F, iSTARTPOS_B, iOFFSET_B, iFltPtchNum, iBndPtchNum, fMD_UnitSlip, fMD_ShearMod, fMD_Lambda, iTDg_V1, iTDg_V2, iTDg_V3, fVDg_Epos, fVDg_Npos, fVDg_Zpos, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fTDl_Curr_DcVal, fTDl_StatFric, fTDl_DynFric, fTDl_RefNormStrss, fTDl_Area, fK_FF_SS, fK_FF_SD, fK_FF_SO, fK_FF_DS, fK_FF_DD, fK_FF_DO, fK_FF_OS, fK_FF_OD, fK_FF_OO, fK_FB_SS, fK_FB_SD, fK_FB_SO, fK_FB_DS, fK_FB_DD, fK_FB_DO, fK_BF_SS, fK_BF_SD, fK_BF_SO, fK_BF_DS, fK_BF_DD, fK_BF_DO, fK_BF_OS, fK_BF_OD, fK_BF_OO, fK_BB_SS, fK_BB_SD, fK_BB_SO, fK_BB_DS, fK_BB_DD, fK_BB_DO, fK_BB_OS, fK_BB_OD, fK_BB_OO, fKl_BB_SS, fKl_BB_DD, fKl_BB_OO, iTDl_SelfLoc_F,  iTDl_SelfLoc_B, iTDl_StabType);

    fMD_CutStrss    = round(fabs(1.0E-6*fMeanLegLgth*fK_FF_SS[0]));
   
    if (iRANK == 0)            {         fprintf(stdout,"K-matrix done \n");                        }    
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/ 
    if ((iPrintAll > 0) && (iRANK == 0))        
    {   fprintf(stdout,"fMD_CutStrss: %f   , representing ~ %f m of slip\n",fMD_CutStrss, 1.0E-6*fMeanLegLgth);
    }
  
    for (i = 0; i < iOFFSET_F[iRANK]; i++)           {       if (fabs(fTDl_SlipRate[i]) > 0.0)       {   iTempValue[0] = 1;    }      }

    MPI_Allreduce(MPI_IN_PLACE, iTempValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   
    if (iTempValue[0] == 1)
    {    
      //this will be done beginning of next week... have look over everything and it seems all fine; clarified a few points that were mysterious to me;...
      //I had this part before (in MATLAB); so, bringing it here should hopefully not take too long
       // void  GetTheStressingRates(
       //the idea is that I now know the RefStrssRateStk and Dip; 
       //use slip on faults and interaction matrix (the "_FB_", "_BB_", and "_BF_") to get first the induced stress at boundary faults, then release via slip and then use that slip
       //for loading on interaction faults; this is added to the existing stressing rate
    //i have this in some matlab function... should be able to more or less copy it out... could also "limit" the number of iterations; although this part should be reasonably fast anyways...

    }
    /*-------------------------------------------------------------------*/   
    for (i = 0; i < iOFFSET_F[iRANK]; i++)
    {   if (iTDl_StabType[i] != 1 )/* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
        {    
            fTempSlipStk =  fTDl_RefStrssRateStk[i]/fK_FF_SS[iTDl_SelfLoc_F[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
            fTempSlipDip =  fTDl_RefStrssRateDip[i]/fK_FF_DD[iTDl_SelfLoc_F[i]];
            for (j = 0; j < iFltPtchNum;  j++)
            {    if (iTDl_StabType[j] == 1) /* apply the creep slip to load the  unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */    
                {   iVectPos         = i*iFltPtchNum + j;                  
                    fTDg_StrssRteChgStk[j] += fTempSlipStk*fK_FF_SS[iVectPos] + fTempSlipDip*fK_FF_DS[iVectPos];
                    fTDg_StrssRteChgDip[j] += fTempSlipStk*fK_FF_SD[iVectPos] + fTempSlipDip*fK_FF_DD[iVectPos];     
        }   }   }
        fTDl_CharRlxTime_F[i] = fabs(fViscosity/(fK_FF_SS[iTDl_SelfLoc_F[i]] < fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]]));  
    }
    for (i = 0; i < iOFFSET_B[iRANK]; i++) 
    {   fTDl_CharRlxTime_B[i] = fabs(fViscosity/(fK_BB_SS[iTDl_SelfLoc_B[i]] < fK_BB_DD[iTDl_SelfLoc_B[i]] ? fK_BB_SS[iTDl_SelfLoc_B[i]] : fK_BB_DD[iTDl_SelfLoc_B[i]])); 
    }
    /*-------------------------------------------------------------------*/
    MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgStk, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
    MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgDip, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);    
    /*-------------------------------------------------------------------*/
    for (i = 0; i < iOFFSET_F[iRANK]; i++) 
    {   
        if (iTDl_StabType[i] == 1) /* if it is an unstable atch/cell... =>apply the creep slip to load */
        {   fTDl_CurStrssRateStk[i]  = fTDl_RefStrssRateStk[i] +fTDg_StrssRteChgStk[i+iSTARTPOS_F[iRANK]];
            fTDl_CurStrssRateDip[i]  = fTDl_RefStrssRateDip[i] +fTDg_StrssRteChgDip[i+iSTARTPOS_F[iRANK]];  
        }   
        else
        {   fTDl_CurStrssRateStk[i]  = 0.0;
            fTDl_CurStrssRateDip[i]  = 0.0;
    }   }
    /*-------------------------------------------------------------------*/  
    /*-------------------------------------------------------------------*/
    free(fMD_Vp);           free(fMD_Vs);           free(fMD_Poisson);         
    free(fMD_MedDense);     free(fMD_Lambda);       free(fMD_AddNrmStrss);
    free(fSD_StrssRate);    free(fSD_SlipRate);     free(fSD_SlipRake);   
    free(iTDg_V1);          free(iTDg_V2);          free(iTDg_V3); 
    free(fVDg_Epos);        free(fVDg_Npos);        free(fVDg_Zpos); 
    free(fTDg_CentEpos);    free(fTDg_CentNpos);    free(fTDg_CentZpos);
    
    free(fK_BB_SD);         free(fK_BB_SO);         free(fK_BB_DS);
    free(fK_BB_DO);         free(fK_BB_OS);         free(fK_BB_OD);  
    free(fK_BB_SS);         free(fK_BB_DD);         free(fK_BB_OO);        
    free(iMD_GlobTTmax);  

    /*-------------------------------------------------------------------*/
    MPI_File_delete(cFile1_Out, MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, cFile1_Out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
    if (iRANK == 0)           
    {   
        MPI_File_write(fp_MPIout1, &iEQcounter,   1, MPI_INT,      &status);
        MPI_File_write(fp_MPIout1, &iUsdGrid,     1, MPI_INT,      &status);            
        MPI_File_write(fp_MPIout1, fMD_ShearMod,  1, MPI_FLOAT,    &status);    
        MPI_File_write(fp_MPIout1, &fdeltTincr,   1, MPI_FLOAT,    &status);           
    }
    MPI_File_close(&fp_MPIout1);
    MPI_Barrier( MPI_COMM_WORLD );
    
    /*load all fault patches to be fully loaded -subtracting 0 to 5MPa from that.... */
    for (i = 0; i < iOFFSET_F[iRANK]; i++) 
    {   fTemp             = sqrtf(fTDl_RefStrssRateStk[i]*fTDl_RefStrssRateStk[i] +fTDl_RefStrssRateDip[i]*fTDl_RefStrssRateDip[i]);
        fTemp1            = (fTDl_RefStrssRateStk[i] < 0.0) ?  -1.0 : 1.0; /*need to get the sign =>b/c loading could be in negative direction I and I want to pre-stress but not above the strength level; hence need to add/subtract depending on direction of stressing rate*/
        fTDl_StrssB4_H[i] = fTDl_RefStrssRateStk[i]/fTemp *(fTDl_CurrFric[i]*fTDl_RefNormStrss[i]);
        fTDl_StrssB4_H[i]-= fTDl_StrssB4_H[i]*fTemp1*(fRandVector[iRandPos[0]]*0.01+0.05);                    iRandPos[0]++;    if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     }       
        fTemp1            = (fTDl_RefStrssRateDip[i] < 0.0) ?  -1.0 : 1.0;
        fTDl_StrssB4_V[i] = fTDl_RefStrssRateDip[i]/fTemp *(fTDl_CurrFric[i]*fTDl_RefNormStrss[i]);
        fTDl_StrssB4_V[i]-= fTDl_StrssB4_V[i]*fTemp1*(fRandVector[iRandPos[0]]*0.01 +0.05);                   iRandPos[0]++;    if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     }   
        fTDl_StrssB4_N[i] = 0.0;
    }   
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
    fTimeYears = 0.0;
    while (fTimeYears <= fRecLgth)
    {   fTimeYears      += fIntSeisTimeStp;
        iTempValue[0] = 0;
        /*------------------------------------------------------------------------------*/
        if (iUsePostSeismic == 1) /*  https://en.wikipedia.org/wiki/Stress_relaxation    http://web.mit.edu/course/3/3.11/www/modules/visco.pdf  => page 9ff; using Maxwell spring-dashpot model    */
        {   //fprintf(stdout,"-");
            /*-------------------------------------------------------------------*/
            for (j = 0; j < iFltPtchNum;  j++)         {    fTDg_DeltStrssH[j]   = 0.0;                fTDg_DeltStrssV[j]   = 0.0;             fTDg_DeltStrssN[j]   = 0.0;            } 
            /*-------------------------------------------------------------------*/
            for (i = 0; i < iOFFSET_F[iRANK]; i++) 
            {      
                fTemp            = expf(-(fTimeStpInSecs*fTDl_PSeisStrssTStep_F[i])/fTDl_CharRlxTime_F[i]);
                fPostSeisIncrNrm = fTDl_StrssB4_N[i] -fTDl_PSeisStrssAt_t0_N_F[i]*fTemp;/*is current deviation from RefNrmStrss minus the expected value after given relaxation time; the difference needs to be released*/    
                fTemp1           = fTDl_PSeisStrssAt_t0_S_F[i] - (fTDl_RefNormStrss[i] +fTDl_PSeisStrssAt_t0_N_F[i]*fTemp)*fTDl_CurrFric[i]; /*this is amount of shear stress when event stopped minus the updated frictional strength; it determines the amount of excess stress, assuming the updated frictional strength was present at the time t0*/
                fTemp2           = sqrtf(fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);
                fTemp3           = fTemp2 - fTemp1*fTemp;/*this is currently applied amount of shear stress minus the expected value -which comes from the assumed excess and the time decay; defines how much stress can be removed*/
                /*remember that fTemp is the decay fraction (value between 0 and 1); multiplied with fTemp1 determines the expected value of shear stress at given time; this is subtracted from the actual amount fTemp2 to get difference that can be released*/
                if ((fTemp3 > 0.0) && (iTDl_StabType[i] != 1))  /*if there is shear stress above the current strength value and if the patch is not-stick-slip....*/ 
                {   fprintf(stdout," AM I EVEN GETTING TO THIS?");
                    fPostSeisIncrStk   = fTemp3/fTemp2 *fTDl_StrssB4_H[i]; /*gives me component to release in current time step...*/
                    fPostSeisIncrDip   = fTemp3/fTemp2 *fTDl_StrssB4_V[i];
                }
                else
                {   fPostSeisIncrStk   = 0.0;
                    fPostSeisIncrDip   = 0.0;
                }   /*so, if I use postseismic, then I can release excess normal stress on all patches, and shear stress on non-stick-slip patches*/
                fTDl_PSeisStrssTStep_F[i] += 1.0;     
                    
                fTemp1           = -1.0*fPostSeisIncrStk/fK_FF_SS[iTDl_SelfLoc_F[i]]; /* slip in horizontal */
                fTemp2           = -1.0*fPostSeisIncrDip/fK_FF_DD[iTDl_SelfLoc_F[i]]; /* slip in vertical  */
                fTemp3           = -1.0*fPostSeisIncrNrm/fK_FF_OO[iTDl_SelfLoc_F[i]]; /* slip in normal  */  
                
                if (sqrt(fTemp1*fTemp1 +fTemp2*fTemp2 +fTemp3*fTemp3) > 0.01*fMD_UnitSlip) /*if slip is less then 1% of the unitslip (which in turn is 0.1% of mean leg length (per given increment) then skip...*/
                {   for (j = 0; j < iFltPtchNum;  j++) /* here is then the relaxation i.e., stress redistribution according to the amount of loading that will be used this time step */
                    {   /*do I need to test there that only the patches with proper stab. type are getting this?*/
                        iVectPos            = i*iFltPtchNum + j;        
                        fTDg_DeltStrssH[j] += fTemp1*fK_FF_SS[iVectPos] +fTemp2*fK_FF_DS[iVectPos] +fTemp3*fK_FF_OS[iVectPos];
                        fTDg_DeltStrssV[j] += fTemp1*fK_FF_SD[iVectPos] +fTemp2*fK_FF_DD[iVectPos] +fTemp3*fK_FF_OD[iVectPos];
                        fTDg_DeltStrssN[j] += fTemp1*fK_FF_SO[iVectPos] +fTemp2*fK_FF_DO[iVectPos] +fTemp3*fK_FF_OO[iVectPos];
            }   }   }   

            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssN, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            /*-------------------------------------------------------------------*/  
            for (i = 0; i < iOFFSET_F[iRANK]; i++) 
            {    fTDl_StrssB4_H[i]   += fTDg_DeltStrssH[i+iSTARTPOS_F[iRANK]];                
                 fTDl_StrssB4_V[i]   += fTDg_DeltStrssV[i+iSTARTPOS_F[iRANK]];    
                 fTDl_StrssB4_N[i]   += fTDg_DeltStrssN[i+iSTARTPOS_F[iRANK]];    
            }      
            /*-------------------------------------------------------------------*/
            for (j = 0; j < iFltPtchNum;  j++)         {    fTDg_DeltStrssH[j]   = 0.0;                fTDg_DeltStrssV[j]   = 0.0;                 fTDg_DeltStrssN[j]   = 0.0; } 
            /*-------------------------------------------------------------------*/
            for (i = 0; i < iOFFSET_B[iRANK]; i++)  //then the same thing for the boundary patches...
            {   
                fTemp            = expf(-(fTimeStpInSecs*fTDl_PSeisStrssTStep_B[i])/fTDl_CharRlxTime_B[i]);
                fPostSeisIncrNrm = fTDl_StrssOnBnd_N[i] -fTDl_PSeisStrssAt_t0_N_B[i]*fTemp;/*is current deviation from RefNrmStrss minus the expected value after given relaxation time; the difference needs to be released*/
                fTemp1           = fTDl_PSeisStrssAt_t0_S_B[i];
                fTemp2           = sqrtf(fTDl_StrssOnBnd_H[i]*fTDl_StrssOnBnd_H[i] +fTDl_StrssOnBnd_V[i]*fTDl_StrssOnBnd_V[i]);
                fTemp3           = fTemp2 - fTemp1*fTemp;
                
                if (fTemp3 > 0.0) 
                {   fPostSeisIncrStk = fTemp3/fTemp2 *fTDl_StrssOnBnd_H[i];
                    fPostSeisIncrDip = fTemp3/fTemp2 *fTDl_StrssOnBnd_V[i];
                }
                else
                {   fPostSeisIncrStk = 0.0;
                    fPostSeisIncrDip = 0.0;
                }
                fTDl_PSeisStrssTStep_B[i] += 1.0;  
                 
                fTemp1           = -1.0*fPostSeisIncrStk/fKl_BB_SS[i]; /* slip in horizontal */
                fTemp2           = -1.0*fPostSeisIncrDip/fKl_BB_DD[i]; /* slip in vertical  */
                fTemp3           = -1.0*fPostSeisIncrNrm/fKl_BB_OO[i]; /* slip in vertical  */
                    
                if (sqrt(fTemp1*fTemp1 +fTemp2*fTemp2 +fTemp3*fTemp3) > 0.01*fMD_UnitSlip)
                {   fTDl_StrssOnBnd_H[i] -= fPostSeisIncrStk; /*subtract stress on boundary patch; */
                    fTDl_StrssOnBnd_V[i] -= fPostSeisIncrDip;/*is done this way b/c the boundary patches are not in the fltptchnum list*/
                    fTDl_StrssOnBnd_N[i] -= fPostSeisIncrNrm;
                
                    for (j = 0; j < iFltPtchNum;  j++) /* here is then the relaxation i.e., stress redistribution according to the amount of loading that will be used this time step */
                    {   /*do I need to test there that only the patches with proper stab. type are getting this?*/
                        iVectPos            = i*iFltPtchNum + j;
                         fTDg_DeltStrssH[j] += fTemp1*fK_BF_SS[iVectPos] +fTemp2*fK_BF_DS[iVectPos] +fTemp3*fK_BF_OS[iVectPos];
                        fTDg_DeltStrssV[j] += fTemp1*fK_BF_SD[iVectPos] +fTemp2*fK_BF_DD[iVectPos] +fTemp3*fK_BF_OD[iVectPos];
                        fTDg_DeltStrssN[j] += fTemp1*fK_BF_SO[iVectPos] +fTemp2*fK_BF_DO[iVectPos] +fTemp3*fK_BF_OO[iVectPos]; 
            }   }   }

            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssN, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);   
            /*-------------------------------------------------------------------*/
            for (i = 0; i < iOFFSET_F[iRANK]; i++) 
            {    fTDl_StrssB4_H[i]   += fTDg_DeltStrssH[i+iSTARTPOS_F[iRANK]];                
                 fTDl_StrssB4_V[i]   += fTDg_DeltStrssV[i+iSTARTPOS_F[iRANK]];      
                 fTDl_StrssB4_N[i]   += fTDg_DeltStrssN[i+iSTARTPOS_F[iRANK]];       
            }   
        }     
        /*------------------------------------------------------------------------------*/
        for (i = 0; i < iOFFSET_F[iRANK]; i++) 
        {   if (iTDl_StabType[i] == 1 )
            {   fTDl_StrssB4_H[i] += fTDl_CurStrssRateStk[i]*fIntSeisTimeStp;
                fTDl_StrssB4_V[i] += fTDl_CurStrssRateDip[i]*fIntSeisTimeStp;
                fTemp              = sqrtf (fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);
                fTemp1             = fTemp - fTDl_CurrFric[i]*(fTDl_RefNormStrss[i] +fTDl_StrssB4_N[i]); /*this is the amount of excess shear stress above current frictional strength*/
                /*now, typically I should also check if the slip that can be released is > DcVal... but that is actually not necessary b/c I only get into this part of code with stick-slip patches; they have by definition more slip than Dc (otherwise they would be cond. stable or stable patches) */
                if (fTemp1 > fMD_CutStrss) /* if my excess stress is at least larger than CutStrss value; */
                {    iTempValue[0]     = 1;        iTDl_WasActivated[i]  = 1;          iTDl_t0[i]            = 0;
        }   }   }
        MPI_Allreduce(iTempValue, iEQ_EQongoing, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        /*-------------------------------------------------------------------*/
        if (iEQ_EQongoing[0] == 1)    /* EARTHQUAKE STARTS */    
        {   iEQ_TotalRuptT[0] = -1;                              iEQ_MRFlength[0]      = -1;                          iEQcounter++;    
            
            for (i = 0; i < iOFFSET_F[iRANK]; i++) 
            {   fTDl_CurStrss_H[i] = fTDl_StrssB4_H[i];          fTDl_CurStrss_V[i]    = fTDl_StrssB4_V[i];           fTDl_CurStrss_N[i]    = fTDl_RefNormStrss[i] +fTDl_StrssB4_N[i];                          
            }
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
            /* EARTHQUAKE ITERATION LOOP STARTS */
            while (iEQ_EQongoing[0] == 1)  /*main variable to tell if I can leave the loop*/
            {   
                iEQ_TotalRuptT[0]++;                iEQ_EQongoing[0] = 0;/*continuously count time since initiation, set "ongoing" to FALSE => only if more slip on patches is added, it gets to be reset*/
                
                for (i = 0; i < iFltPtchNum;  i++)          {        fTDg_DeltStrssH[i] = 0.0;            fTDg_DeltStrssV[i] = 0.0;            fTDg_DeltStrssN[i] = 0.0;       }
                
                for (i = 0; i < iOFFSET_F[iRANK]; i++) 
                {    /* determine new friction coefficient */
                    
                    if (iTDl_WasActivated[i] == 1)
                    {   STFcnt    = iEQ_TotalRuptT[0] -iTDl_t0[i]; 
                        if      (iSD_FricLawUSED[ iTDl_SegID[ i ] ] == 1) /* classic static/dynamic friction */
                        {    
                            fTDl_CurrFric[i] = fTDl_DynFric[i]; 
                        }
                        else if (iSD_FricLawUSED[ iTDl_SegID[ i ] ] == 2) /* slip weakening -using accumulated slip of current event (on that patch) */
                        {   
                            fTemp            = sqrtf(fTDl_EventSlipH[i]*fTDl_EventSlipH[i] +fTDl_EventSlipV[i]*fTDl_EventSlipV[i]);
                            fTDl_CurrFric[i] = fTDl_DynFric[i] + (fTDl_StatFric[i]- fTDl_DynFric[i])* (1.0 - fTemp/fTDl_Curr_DcVal[i]); /* the first part is the current friction then the second part defines delta_mue => change in friction coefficient; and the second one is the  */
                            if      ((iTDl_StabType[i] != 3) && (fTDl_CurrFric[i] < fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
                            else if ((iTDl_StabType[i] == 3) && (fTDl_CurrFric[i] > fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }   

                        }
                        else if (iSD_FricLawUSED[ iTDl_SegID[ i ] ] == 3) /* velocity weakening -using slip increment from next iterative slip (with deltT that equates to velocity) */
                        {   
                            ireadPos         = (STFcnt == 0) ? 0 : (STFcnt%iMaxSTFLength -1);
                            
                            fTemp            = sqrtf(fTDlg_STF_H[i*iMaxSTFLength+ireadPos]*fTDlg_STF_H[i*iMaxSTFLength+ireadPos] + fTDlg_STF_V[i*iMaxSTFLength+ireadPos]*fTDlg_STF_V[i*iMaxSTFLength+ireadPos]);
                        //    if (iDeltT2Use == 1)    {       fTemp           /= fdeltTincr;              } /*this is a slip velocity now*/
                        //    else                    {       fTemp            = fTemp;                   } /*assumes that Dc is then defined not as critical slip, but critical velocity => Dc is here a velocity!!!; if DeltT is not used then I compare slip from last step with Dc slip; if DeltT is used, divide that slip by time incremeent and compare with Dc_velocity*/
                            fTDl_CurrFric[i] = fTDl_DynFric[i] + (fTDl_StatFric[i]- fTDl_DynFric[i])* (1.0 - fTemp/fTDl_Curr_DcVal[i]); 
                            if      ((iTDl_StabType[i] != 3) && (fTDl_CurrFric[i] < fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
                            else if ((iTDl_StabType[i] == 3) && (fTDl_CurrFric[i] > fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }    
                    }   }
                    else
                    {   STFcnt    = 0;   
                    }
                    /*-------------------------------------------------------*/
                    /* determine amount of excess stress (if any available) and the corresponding slip amount */
                    fTemp  = sqrtf(fTDl_CurStrss_H[i]*fTDl_CurStrss_H[i] + fTDl_CurStrss_V[i]*fTDl_CurStrss_V[i]);  /* this is the applied shear stress*/
                    fTemp1 = (fTDl_CurStrss_N[i] > 0.0) ? fTDl_CurStrss_N[i] : 0.0;/*nice and simple idea here: want to make sure that the normal stress on fault is not becoming smaller than zero => use the smaller of the two values; this allows for proper "evolution" of normal stress due to different signal speed of P and S*/    
                    fTemp2 = fTemp - fTDl_CurrFric[i]*fTemp1; /* this is the amount of stress above current friction level in combination with normal stress at the location*/                    
                    fTemp3 = (fTDl_DynFric[i] *fTemp1 - fTemp)/(fK_FF_SS[iTDl_SelfLoc_F[i]] < fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]]); /* this is releasable stress divided by self-stiffness => gives slip amount that would happen if patch fails "alone".... */                    

                    if ((iTDl_WasActivated[i] == 1)  || ((iTDl_WasActivated[i] == 0)  && (fTemp2 > fMD_CutStrss) && (fTemp3 > fTDl_Curr_DcVal[i])))
                    {    
                        if (iTDl_WasActivated[i] == 0)                  {   iTDl_WasActivated[i] = 1;           iTDl_t0[i] = iEQ_TotalRuptT[0];        STFcnt    = 0;   }
                        
                        fTemp2    = (fTemp2 > fMD_CutStrss) ? fTemp2 : 0.0;
                        
                        if (fTemp2 > fMD_CutStrss)                       
                        {  
                            iEQ_EQongoing[0] = 1;                       
                            iEQ_MRFlength[0] = iEQ_TotalRuptT[0];              
                            iwritePos        = STFcnt%iMaxSTFLength; /*gives me remainder of wPos/MaxSTFlength => allows me to loop/overwrite the STF vectors */

                            fTDlg_STF_H[i*iMaxSTFLength+iwritePos] = -1.0*(fTemp2/fTemp *fTDl_CurStrss_H[i]) /fK_FF_SS[iTDl_SelfLoc_F[i]]; /* horizontal slip at patch in current iteration (turn "n") */
                            fTDlg_STF_V[i*iMaxSTFLength+iwritePos] = -1.0*(fTemp2/fTemp *fTDl_CurStrss_V[i]) /fK_FF_DD[iTDl_SelfLoc_F[i]]; /* the "-1" is here because the stiffness matrix stuff gives the amount of slip nessessary to MAKE the observed stress; but I want to RELEASE it => opposite direction */
                            fTDl_EventSlipH[i]                    += fTDlg_STF_H[i*iMaxSTFLength+iwritePos];                
                            fTDl_EventSlipV[i]                    += fTDlg_STF_V[i*iMaxSTFLength+iwritePos];   
                        
                            if (iEQ_MRFlength[0] < iMaxMRFlength) /*write out the moment-rate-function, but only to a maximum of maxMRFlength => can still do the earthquake but cutoff the MRF...*/
                            {   fEQ_MRFvals[iEQ_TotalRuptT[0]]    += sqrtf(fTDlg_STF_H[i*iMaxSTFLength+iwritePos]*fTDlg_STF_H[i*iMaxSTFLength+iwritePos] +fTDlg_STF_V[i*iMaxSTFLength+iwritePos]*fTDlg_STF_V[i*iMaxSTFLength+iwritePos])*fTDl_Area[i] *fMD_ShearMod[0];
                        }   }
                        else /*seems unnecessary to add here if fTemp2 is zero -but is! necessary -even if it is to set that position to be zero!*/
                        {   iwritePos  = STFcnt%iMaxSTFLength; /*gives me remainder of wPos/MaxSTFlength => allows me to loop/overwrite the STF vectors */
                            fTDlg_STF_H[i*iMaxSTFLength+iwritePos] = 0.0; /* it is necessary b/c that position would remain its previous value, and I am looping over this stuff => would reuse previous slip amount */
                            fTDlg_STF_V[i*iMaxSTFLength+iwritePos] = 0.0;
                    }   }
                    /*---------------------------------------------------------*/
                    /* determine modified stress due to slip on fault patches */
                    if (iTDl_WasActivated[i] == 1) 
                    {   /*---------------------------------------------------------*/                         
                        for (j = 0; j < iFltPtchNum;  j++) /*taking the receivers..*/
                        {
                            iVectPos  = i*iFltPtchNum + j;    
                            /*---------------------------------------------------------*/    
                            for (stfpos = iTDlg_Pdone[iVectPos]; stfpos <= (STFcnt - iTDlg_TravTimesP[iVectPos]); stfpos++) /*2nd part is negative if signal has not arrived, then the loop is skipped*/
                            {  
                                ireadPos = stfpos%iMaxSTFLength;
                                fStkSlip = fTDlg_STF_H[i*iMaxSTFLength +ireadPos];
                                fDipSlip = fTDlg_STF_V[i*iMaxSTFLength +ireadPos];
                                fTemp    = sqrt(fStkSlip*fStkSlip + fDipSlip*fDipSlip);   /*might be that this part of the STF is without slip, then it is not necessary to go trough the next steps*/
                                if (fTemp > FLT_EPSILON)
                                { 
                                   if (i+iSTARTPOS_F[iRANK] == j)
                                    {   fUsedStkSlip        = fStkSlip; /*this part here is actually a really cool step; fTemp is length of slip amount in direction of source-receiver vector (can be negative)*/
                                        fUsedDipSlip        = fDipSlip;
                                        fTDg_DeltStrssH[j] += fUsedStkSlip*fK_FF_SS[iVectPos] +fUsedDipSlip*fK_FF_DS[iVectPos];
                                        fTDg_DeltStrssV[j] += fUsedStkSlip*fK_FF_SD[iVectPos] +fUsedDipSlip*fK_FF_DD[iVectPos];
                                    }
                                    else
                                    {   fTemp               = fStkSlip*fTDlg_LocSrcRcv_H[iVectPos] + fDipSlip*fTDlg_LocSrcRcv_V[iVectPos]; 
                                        fUsedStkSlip        = fTDlg_LocSrcRcv_H[iVectPos]*fTemp; /*this part here is actually a really cool step; fTemp is length of slip amount in direction of source-receiver vector (can be negative)*/
                                        fUsedDipSlip        = fTDlg_LocSrcRcv_V[iVectPos]*fTemp; /*this "length" is multiplied with local orientation of source-receiver vector to get the transient slip values to be used...*/
                                        fUsedNrmSlip        = fTDlg_LocSrcRcv_N[iVectPos]*fTemp;
                                        fTDg_DeltStrssH[j] += fUsedStkSlip*fK_FF_SS[iVectPos] +fUsedDipSlip*fK_FF_DS[iVectPos] +fUsedNrmSlip*fK_FF_OS[iVectPos];
                                        fTDg_DeltStrssV[j] += fUsedStkSlip*fK_FF_SD[iVectPos] +fUsedDipSlip*fK_FF_DD[iVectPos] +fUsedNrmSlip*fK_FF_OD[iVectPos];
                                        fTDg_DeltStrssN[j] += fUsedStkSlip*fK_FF_SO[iVectPos] +fUsedDipSlip*fK_FF_DO[iVectPos] +fUsedNrmSlip*fK_FF_OO[iVectPos];   
                            }   }   }                      
                            iTDlg_Pdone[iVectPos] = (STFcnt - iTDlg_TravTimesP[iVectPos]) >= 0 ? (STFcnt - iTDlg_TravTimesP[iVectPos]+1) : 0;        
                            /*---------------------------------------------------------*/                
                            for (stfpos = iTDlg_Sdone[iVectPos]; stfpos <= (STFcnt - iTDlg_TravTimesS[iVectPos]); stfpos++) /*2nd part is negative if signal has not arrived, then the loop is skipped*/
                            {   
                                ireadPos = stfpos%iMaxSTFLength;
                                fStkSlip = fTDlg_STF_H[i*iMaxSTFLength +ireadPos];
                                fDipSlip = fTDlg_STF_V[i*iMaxSTFLength +ireadPos];   
                                fTemp    = sqrt(fStkSlip*fStkSlip + fDipSlip*fDipSlip);                                
                                if (fTemp > FLT_EPSILON)
                                {   if ((i + iSTARTPOS_F[iRANK]) != j)    /*this should be at "self" location -> then I cannot use the src-rcv-vector b/c it is zero...*/         
                                    {   fTemp               = fStkSlip*fTDlg_LocSrcRcv_H[iVectPos] + fDipSlip*fTDlg_LocSrcRcv_V[iVectPos]; 
                                        fUsedStkSlip        = fStkSlip - fTDlg_LocSrcRcv_H[iVectPos]*fTemp; /*subtract the P-direction (as above) from the actual slip vector => this gives the S-velocity component*/
                                        fUsedDipSlip        = fDipSlip - fTDlg_LocSrcRcv_V[iVectPos]*fTemp; 
                                        fUsedNrmSlip        =    0.0   - fTDlg_LocSrcRcv_N[iVectPos]*fTemp; 
                                                                           
                                        fTDg_DeltStrssH[j] += fUsedStkSlip*fK_FF_SS[iVectPos] +fUsedDipSlip*fK_FF_DS[iVectPos] +fUsedNrmSlip*fK_FF_OS[iVectPos];
                                        fTDg_DeltStrssV[j] += fUsedStkSlip*fK_FF_SD[iVectPos] +fUsedDipSlip*fK_FF_DD[iVectPos] +fUsedNrmSlip*fK_FF_OD[iVectPos];
                                        fTDg_DeltStrssN[j] += fUsedStkSlip*fK_FF_SO[iVectPos] +fUsedDipSlip*fK_FF_DO[iVectPos] +fUsedNrmSlip*fK_FF_OO[iVectPos];
                            }   }   }    
                            iTDlg_Sdone[iVectPos] = (STFcnt - iTDlg_TravTimesS[iVectPos]) >= 0 ? (STFcnt- iTDlg_TravTimesS[iVectPos]+1) : 0; 
                }   }   }   
                /*---------------------------------------------------------------*/
                /* last was done "locally" now transfer information to all patches*/
                MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssN, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, iEQ_EQongoing,      1        , MPI_INT,MPI_MAX, MPI_COMM_WORLD);
                /* combine "stress before EQ" and "EQ generated" changes to get "Currently applied stress" */    
                for (i = 0; i < iOFFSET_F[iRANK]; i++) 
                {    
                    fTDl_CurStrss_H[i] += fTDg_DeltStrssH[i+iSTARTPOS_F[iRANK]];
                    fTDl_CurStrss_V[i] += fTDg_DeltStrssV[i+iSTARTPOS_F[iRANK]];
                    fTDl_CurStrss_N[i] += fTDg_DeltStrssN[i+iSTARTPOS_F[iRANK]];
                }   
                /*---------------------------------------------------------------*/
                if (iEQ_EQongoing[0] == 0) /* make sure that all the signal is out of the system => even if no slip occurred in last step on any patch; there may still be stress in the system that has not reached a receiver => wait until they all got their share */
                {   
                    iEQ_EQongoing[0]   = 1;                       iEQ_EQEndCntr[0]  = iEQ_EQEndCntr[0] +1;
                    if (iEQ_EQEndCntr[0] >= iMaxSTFLength)    {   iEQ_EQongoing[0]  = 0;                            }   
                }
                else                                          {    iEQ_EQEndCntr[0] = 0;                            }
                MPI_Allreduce(MPI_IN_PLACE, iEQ_EQongoing,      1        , MPI_INT,MPI_MAX, MPI_COMM_WORLD);                        
                /*---------------------------------------------------------------*/
            }  /* EARTHQUAKE ITERATION LOOP EXITED */            
            
            for (i = 0; i < iOFFSET_F[iRANK]; i++)    /*here I load the boundary faults, if they are used...; this will then be released via post-seis. deformation (if turned on); otherwise it will not be used...*/
            {   if (iTDl_WasActivated[i] == 1)
                {    for (j = 0; j < iBndPtchNum;  j++)
                    {   iVectPos = i*iBndPtchNum +j;
                        fTDg_TempBnd_H[j] += fTDl_EventSlipH[i]*fK_FB_SS[iVectPos] +fTDl_EventSlipV[i]*fK_FB_DS[iVectPos];
                        fTDg_TempBnd_V[j] += fTDl_EventSlipH[i]*fK_FB_SD[iVectPos] +fTDl_EventSlipV[i]*fK_FB_DD[iVectPos];
                        fTDg_TempBnd_N[j] += fTDl_EventSlipH[i]*fK_FB_SO[iVectPos] +fTDl_EventSlipV[i]*fK_FB_DO[iVectPos];
            }   }   }
            MPI_Allreduce(MPI_IN_PLACE, fTDg_TempBnd_H, iBndPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_TempBnd_V, iBndPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, fTDg_TempBnd_N, iBndPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            for (i = 0; i < iOFFSET_B[iRANK]; i++) 
            {    
                fTDl_StrssOnBnd_H[i] += fTDg_TempBnd_H[i+iSTARTPOS_B[iRANK]];
                fTDl_StrssOnBnd_V[i] += fTDg_TempBnd_V[i+iSTARTPOS_B[iRANK]];
                fTDl_StrssOnBnd_N[i] += fTDg_TempBnd_N[i+iSTARTPOS_B[iRANK]];
            }          
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
            iEQ_ActPtchNum[0]  = 0;
            iEQ_MRFlength[0]   = (iEQ_MRFlength[0] < iMaxMRFlength) ? iEQ_MRFlength[0] : iMaxMRFlength;
            
            for (i = 0; i < iOFFSET_F[iRANK]; i++)    /* now I define all the metrics of the event (its size, its MRF, etc...) */
            {   if (iTDl_WasActivated[i] == 1)     
                {   fEQ_SeisPot[0]                     += sqrtf(fTDl_EventSlipH[i]*fTDl_EventSlipH[i] + fTDl_EventSlipV[i]*fTDl_EventSlipV[i] ) *fTDl_Area[i];    
                    iEQl_ActPtchID[iEQ_ActPtchNum[0]]   = i+iSTARTPOS_F[iRANK];
                    iEQl_t0ofPtch[iEQ_ActPtchNum[0]]    = iTDl_t0[i];       
                    fEQl_DeltTofPtch[iEQ_ActPtchNum[0]] = sqrtf((fTDl_CurStrss_H[i]-fTDl_StrssB4_H[i])*(fTDl_CurStrss_H[i]-fTDl_StrssB4_H[i]) + (fTDl_CurStrss_V[i]-fTDl_StrssB4_V[i])*(fTDl_CurStrss_V[i]-fTDl_StrssB4_V[i]));
                    fEQl_SlipHofPtch[iEQ_ActPtchNum[0]] = fTDl_EventSlipH[i];    
                    fEQl_SlipVofPtch[iEQ_ActPtchNum[0]] = fTDl_EventSlipV[i];    
                    iEQl_StabType[iEQ_ActPtchNum[0]]    = iTDl_StabType[i];        
                    iEQ_ActPtchNum[0]++;       
            }   } 
            /*------------------------------------------------------------------ */  
            MPI_Allreduce(iEQ_ActPtchNum,  iEQ_CmbPtchNum,       1          ,    MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,    iEQ_MRFlength,        1          ,    MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,    fEQ_SeisPot,          1          ,    MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,    fEQ_MRFvals,     iMaxMRFlength   ,    MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
            
            MPI_Allgather(iEQ_ActPtchNum,      1,      MPI_INT,  iEQ_WrtStartPos, 1, MPI_INT,  MPI_COMM_WORLD);
            
            for (i = 1;     i < iSIZE; i++)           {    iEQ_WrtStartPos[i] += iEQ_WrtStartPos[i-1];             }    
            for (i = (iSIZE-1); i > 0; i--)           {    iEQ_WrtStartPos[i]  = iEQ_WrtStartPos[i-1];             }    
            iEQ_WrtStartPos[0]  = 0;
            
            /* outoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutout */
            /*------------------------------------------------------------------ */
            MPI_File_open(MPI_COMM_WORLD, cFile1_Out,MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
            MPI_File_get_size(fp_MPIout1, &offset1);
            MPI_Barrier( MPI_COMM_WORLD );
    
            if (iRANK == 0)
            {
                fTemp2 = (log10f(fEQ_SeisPot[0]*fMD_ShearMod[0])-9.1)/1.5; 
                MPI_File_write_at(fp_MPIout1,   offset1,                                       &fTimeYears,          1,            MPI_FLOAT, &status); 
                MPI_File_write_at(fp_MPIout1,   offset1 +  sizeof(float),                      &fTemp2,              1,            MPI_FLOAT, &status);
                MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float),                      iEQ_CmbPtchNum,       1,            MPI_INT,   &status);
                MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float) +1*sizeof(int),       iEQ_MRFlength,        1,            MPI_INT,   &status);
                MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float) +2*sizeof(int),       fEQ_MRFvals,    iEQ_MRFlength[0],   MPI_FLOAT, &status);

                fprintf(stdout,"Earthquake time %f    Magn: %f act patches: %d   MRF length: %d   \n",fTimeYears,  fTemp2, iEQ_CmbPtchNum[0],iEQ_MRFlength[0]);                
            }

            offset2 = offset1 +2*sizeof(int) +(2 +iEQ_MRFlength[0])*sizeof(float);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(int),  iEQl_ActPtchID,    iEQ_ActPtchNum[0], MPI_INT,   &status);
            offset2 = offset2 +iEQ_CmbPtchNum[0]*sizeof(int);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(int),  iEQl_t0ofPtch,     iEQ_ActPtchNum[0], MPI_INT,   &status);
            offset2 = offset2 +iEQ_CmbPtchNum[0]*sizeof(int);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(float), fEQl_DeltTofPtch, iEQ_ActPtchNum[0], MPI_FLOAT, &status);
            offset2 = offset2 +iEQ_CmbPtchNum[0]*sizeof(float);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(float), fEQl_SlipHofPtch, iEQ_ActPtchNum[0], MPI_FLOAT, &status);
            offset2 = offset2 +iEQ_CmbPtchNum[0]*sizeof(float);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(float), fEQl_SlipVofPtch, iEQ_ActPtchNum[0], MPI_FLOAT, &status);
            offset2 = offset2 +iEQ_CmbPtchNum[0]*sizeof(float);
            MPI_File_write_at(fp_MPIout1,   offset2 +iEQ_WrtStartPos[iRANK]*sizeof(int),   iEQl_StabType,    iEQ_ActPtchNum[0], MPI_INT,   &status);
            
            MPI_File_close(&fp_MPIout1);
            MPI_Barrier( MPI_COMM_WORLD );
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */    
            iEQ_ChangedStab[0] = 0; /* this will record if stability type was changed*/
            for (i = 0; i < iOFFSET_F[iRANK]; i++)  /* RESET ALL VALUES FOR PATCHES THAT MOVED */
            {
                if (iTDl_WasActivated[i] == 1)     
                {   if (iMD_ChgFricBtwEQs[0] == 0) //IMPORTANT -FOR TESTING I CHANGED THIS (FLIPPED)
                    {    fTDl_CurrFric[i]     = fTDl_StatFric[i];
                    }
                    else   
                    {   fTDl_StatFric[i]   = fSD_RefStatFric[iTDl_SegID[i]]  *(1.0 +(fRandVector[iRandPos[0]]*2.0 -1.0) *fSD_RefStatFr_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     }       
                        fTDl_CurrFric[i]   = fTDl_StatFric[i];
                        fTDl_DynFric[i]    = fSD_RefDynFric[iTDl_SegID[i]]   *(1.0 +(fRandVector[iRandPos[0]]*2.0 -1.0) *fSD_RefDynFr_vari[ iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     }  
                        fTDl_Curr_DcVal[i] = fSD_CritSlipDist[iTDl_SegID[i]] *(1.0 +(fRandVector[iRandPos[0]]*2.0 -1.0) *fSD_CritSlipD_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     }      
                        /* test now again for stability type, while doing that check if any type was changed, if not then I don't have to recompute, if change is from cond stable to stable then I only change label, only if it goes from unstable to another state or vice versa, is it that I need to update loading*/

                        fTemp1 = (fTDl_DynFric[i] - fTDl_StatFric[i]) *fTDl_RefNormStrss[i]; /* this is the releasable stress amount*/
                        fTemp2 = fTemp1/(fK_FF_SS[iTDl_SelfLoc_F[i]] < fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]]);/* this is corresponding slip, sign is correct -> assuming that K-matrix at self-induced location is ALWAYS negative*/
        
                        if (fTemp2 > fTDl_Curr_DcVal[i]) /* have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc*/
                        {   if (iTDl_StabType[i] != 1)              {   iEQ_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
                            iTDl_StabType[i] = 1; /* patch is unstable*/
                        }
                        else /*if (fTemp2 <= fTDl_Curr_DcVal[i]) -if stress drop i.e., corresponding slip is smaller than Dc - */
                        {   if (fTemp1 < 0.0)  /*  is still weakening but just with slip that is lower than Dc => cond. stable*/
                            {   if (iTDl_StabType[i] != 2)          {   iEQ_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
                                iTDl_StabType[i] = 2; /* patch is cond. stable*/
                            }
                            else /* no weakening but strengthening  => stable */
                            {   if (iTDl_StabType[i] != 3)          {   iEQ_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
                                iTDl_StabType[i] = 3; /* patch is stable*/ 
                }   }   }   }
                for (k = 0; k < iMaxSTFLength; k++)  {  iVectPos = i*iMaxSTFLength + k;       fTDlg_STF_H[iVectPos] = 0.0;            fTDlg_STF_V[iVectPos] = 0.0;        }     
                for (k = 0; k < iFltPtchNum;   k++)  {  iVectPos = i*iFltPtchNum + k;         iTDlg_Pdone[iVectPos] = 0;              iTDlg_Sdone[iVectPos] = 0;          }
                 
                iEQl_ActPtchID[i]    = 0;               iEQl_t0ofPtch[i]     = 0;             fEQl_DeltTofPtch[i]   = 0.0;            fEQl_SlipHofPtch[i]   = 0.0;
                fEQl_SlipVofPtch[i]  = 0.0;             iTDl_WasActivated[i] = 0;             iTDl_t0[i]            = 0;             
                fTDl_EventSlipH[i]   = 0.0;             fTDl_EventSlipV[i]   = 0.0;                        
            }
            iEQ_ActPtchNum[0]        = 0;               iEQ_CmbPtchNum[0]    = 0;             iEQ_MRFlength[0]      = 0;                                        
            iEQ_EQongoing[0]         = 0;               fEQ_SeisPot[0]       = 0.0;           iEQ_EQEndCntr[0]      = 0;
            
            for (i = 0; i < iSIZE; i++)          {      iEQ_WrtStartPos[i]   = 0;                                                                                         }
            for (i = 0; i < iMaxMRFlength; i++)  {      fEQ_MRFvals[i]       = 0.0;           fEQ_MRFvals[i]       = 0.0;                                                 }
            for (j = 0; j < iBndPtchNum; j++)    {      fTDg_TempBnd_H[j]    = 0.0;           fTDg_TempBnd_V[j]    = 0.0;             fTDg_TempBnd_N[j]     = 0.0;        }
            /*------------------------------------------------------------------------- */      
            MPI_Allreduce(MPI_IN_PLACE, iEQ_ChangedStab,       1          ,    MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                      
            if (iEQ_ChangedStab[0] == 1) /* next, update the current stressing rate*/
            {  fprintf(stdout," I DID MAKE IT HERE");       exit(9);
                /*-------------------------------------------------------------------*/
                for (j = 0; j < iFltPtchNum; j++)     {       fTDg_StrssRteChgStk[j]   = 0.0;         fTDg_StrssRteChgDip[j]   = 0.0;         }
                for (i = 0; i < iOFFSET_F[iRANK]; i++) 
                {        
                    if (iTDl_StabType[i] != 1)/* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
                    {   fTempSlipStk = fTDl_RefStrssRateStk[i]/fK_FF_SS[iTDl_SelfLoc_F[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
                        fTempSlipDip = fTDl_RefStrssRateDip[i]/fK_FF_DD[iTDl_SelfLoc_F[i]];
    
                        for (j = 0; j < iFltPtchNum;  j++)
                        {    if (iTDl_StabType[i] == 1) /* apply the creep slip to load the unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */    
                            {   iVectPos        = i*iFltPtchNum + j;                  
                                fTDg_StrssRteChgStk[j] += fTempSlipStk*fK_FF_SS[iVectPos] + fTempSlipDip*fK_FF_DS[iVectPos];
                                fTDg_StrssRteChgDip[j] += fTempSlipStk*fK_FF_SD[iVectPos] + fTempSlipDip*fK_FF_DD[iVectPos];        
                }   }   }   }
                /*-------------------------------------------------------------------*/
                MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgStk, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
                MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgDip, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);   
                /*-------------------------------------------------------------------*/
                for (i = 0; i < iOFFSET_F[iRANK]; i++) 
                {   if (iTDl_StabType[i] == 1) /* if it is a stick-slip patch/cell... =>apply the creep slip to load the  locked/strong patches    */
                    {   fTDl_CurStrssRateStk[i]  = fTDl_RefStrssRateStk[i] +fTDg_StrssRteChgStk[i+iSTARTPOS_F[iRANK]];            
                        fTDl_CurStrssRateDip[i]  = fTDl_RefStrssRateDip[i] +fTDg_StrssRteChgDip[i+iSTARTPOS_F[iRANK]];          
                    }   
                    else
                    {   fTDl_CurStrssRateStk[i]  = 0.0;
                        fTDl_CurStrssRateDip[i]  = 0.0;             
            }   }   }
            /*------------------------------------------------------------------------- */        
            for (i = 0; i < iOFFSET_F[iRANK]; i++) 
            {   if (iUsePostSeismic == 1)
                {   fTDl_StrssB4_H[i]           = fTDl_CurStrss_H[i];
                    fTDl_StrssB4_V[i]           = fTDl_CurStrss_V[i];
                    fTDl_StrssB4_N[i]           = fTDl_CurStrss_V[i] - fTDl_RefNormStrss[i]; 
                    fTDl_PSeisStrssAt_t0_S_F[i] = sqrtf(fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);  
                    fTDl_PSeisStrssAt_t0_N_F[i] = fTDl_StrssB4_N[i];   
                    fTDl_PSeisStrssTStep_F[i]   = 0.0;
                }
                else /*no postseismic, everything above static strength is removed....uses the orientation of the currently applied stress*/
                {   fTemp             = sqrtf(fTDl_CurStrss_H[i]*fTDl_CurStrss_H[i] +fTDl_CurStrss_V[i]*fTDl_CurStrss_V[i]);  /*applied shear stress*/
                    fTemp1            = fTDl_CurStrss_H[i]/fTemp *(fTDl_StatFric[i]*fTDl_RefNormStrss[i]);/*stressValue in H corresponding to strength and given aspect ratio of */
                    fTemp2            = fTDl_CurStrss_V[i]/fTemp *(fTDl_StatFric[i]*fTDl_RefNormStrss[i]);
                    
                    fTDl_StrssB4_H[i] = (fabs(fTDl_CurStrss_H[i]) < fabs(fTemp1)) ? fTDl_CurStrss_H[i] : fTemp1;
                    fTDl_StrssB4_V[i] = (fabs(fTDl_CurStrss_V[i]) < fabs(fTemp2)) ? fTDl_CurStrss_V[i] : fTemp2;
                    fTDl_StrssB4_N[i] = 0.0;
                    /*means that I enforce that all patches (regardless of stabtype) have shear stress that is at or below static strength*/
                    /*for stick-slip patches, that should just be the applied stress; for non-stick-slip I remove the excess stress above static strength...*/        
            }   }
            /*------------------------------------------------------------------ */
            for (i = 0; i < iOFFSET_B[iRANK]; i++) 
            {
                if (iUsePostSeismic == 1)
                {   
                    fTDl_PSeisStrssAt_t0_S_B[i] = sqrtf(fTDl_StrssOnBnd_H[i]*fTDl_StrssOnBnd_H[i] +fTDl_StrssOnBnd_V[i]*fTDl_StrssOnBnd_V[i]);  
                    fTDl_PSeisStrssAt_t0_N_B[i] = fTDl_StrssOnBnd_N[i];
                    fTDl_PSeisStrssTStep_B[i]   = 0.0;
            }   }
            /*------------------------------------------------------------------ */
        }
        MPI_Barrier( MPI_COMM_WORLD );
    }
    
    MPI_File_open(MPI_COMM_WORLD, cFile1_Out, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
    if (iRANK == 0)            {  MPI_File_write_at(fp_MPIout1, 0, &iEQcounter, 1, MPI_INT,  &status);            }
    MPI_File_close(&fp_MPIout1);
    MPI_Barrier( MPI_COMM_WORLD );
       
    MPI_Finalize();
    return 0;
}
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0_inmain(long *idum) /* returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to initialize or re-initialize the sequence.*/
{    
    long k;
    float ans;
    
    *idum ^= MASK;
    k      = (*idum)/IQ;
    *idum  = IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0)  *idum += IM;     
    ans    = AM*(*idum);
    *idum ^= MASK;
    
    if ((ans < 0.0) || (ans > 1.0))
    {   fprintf(stdout, "The random number generator is not behaving properly.... gives values below zero i.e., above one \n"); exit(10);
    }    
    
    return ans;
}
