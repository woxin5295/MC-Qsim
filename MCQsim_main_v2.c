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
extern void LoadInputParameter(const char **argv, const unsigned int iRealizNum, const unsigned int iUsdGrid, const unsigned int iFltSegmNum, const unsigned int iFltPtchNum, const unsigned int iBndPtchNum, const unsigned int iFltVertNum, const unsigned int iBndVertNum, const float *fSD_StrssRate, const float *fSD_SlipRate, const float *fSD_SlipRake, unsigned int *iMD_ChgFricBtwEQs, float *fMD_AddNrmStrss, float *fMD_Vp, float *fMD_Vs, float *fMD_Poisson, float *fMD_Lambda, float *fMD_ShearMod, float *fMD_MedDense, float *fMD_CritSlipVelo, unsigned int *iSD_FricLawUSED,  float *fSD_RefStatFric, float *fSD_RefStatFr_vari, float *fSD_RefDynFric, float *fSD_RefDynFr_vari, float *fSD_CritSlipDist, float *fSD_CritSlipD_vari, float *fTDg_CentEpos, float *fTDg_CentNpos, float *fTDg_CentZpos, unsigned int *iTDg_V1, unsigned int *iTDg_V2, unsigned int *iTDg_V3, unsigned int *iTDg_StabType, unsigned int *iTDg_SegID, float *fTDg_SlipRate, float *fTDg_SlipRake, float *fTDg_StrsRteStk, float *fTDg_StrsRteDip, float *fTDg_StatFric, float *fTDg_DynFric, float *fVDg_Epos, float *fVDg_Npos, float *fVDg_Zpos);						
extern void  DefineMoreParas(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const unsigned ind iFltPtchNum, const unsigned ind iBndPtchNum, const unsigned int *iTDg_V1, const unsigned int *iTDg_V2, const unsigned int *iTDg_V3, const unsigned int *iTDl_SegID, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fSD_CritSlipDist, const float *fSD_CritSlipD_vari, const float fMD_Vp, const float fMD_VpVsRatio, const float fMD_MedDense, const float fMD_AddNrmStrss, const float fMD_g, const float fdeltTincr, const float *fRandVector, const unsigned int iRandNumber, unsigned int *iRandPos, unsigned int *iTDlg_TravTimes, unsigned int *iTDg_GlobTTmax,  float *fTDl_Area, float *fTDl_RefNormStrss, float *fTDl_Curr_DcVal, float *fTDlg_LocSrcRcv_H, float *fTDlg_LocSrcRcv_V, float *fTDlg_LocSrcRcv_N);
extern void   Build_K_Matrix(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const unsigned int iFltPtchNum, const unsigned int iBndPtchNum, const float fMD_UnitSlip, const float fMD_ShearMod, const float fMD_Lambda, const unsigned int *iTDg_V1, const unsigned int *iTDg_V2, const unsigned int *iTDg_V3, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fTDl_Curr_DcVal, const float *fTDl_StatFric, const float *fTDl_DynFric, const float *fTDl_RefNormStrss, const float *fTDl_Area, const float *fTDl_RotMatF, const float *fTDl_RotMatB, float *fK_FF_SS, float *fK_FF_SD, float *fK_FF_SO, float *fK_FF_DS, float *fK_FF_DD, float *fK_FF_DO, float *fK_FF_OS, float *fK_FF_OD, float *fK_FF_OO, float *fK_FB_SS, float *fK_FB_SD, float *fK_FB_SO, float *fK_FB_DS, float *fK_FB_DD, float *fK_FB_DO, float *fK_BF_SS, float *fK_BF_SD, float *fK_BF_DS, float *fK_BF_DD, float *fK_BF_OS, float *fK_BF_OD, float *fK_BB_SS, float *fK_BB_SD, float *fK_BB_SO, float *fK_BB_DS, float *fK_BB_DD, float *fK_BB_DO, float *fK_BB_OS, float *fK_BB_OD, float *fK_BB_OO, unsigned int *iTDl_SelfLoc_F, unsigned int *iTDl_SelfLoc_B, unsigned int *iTDl_StabType);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int main(int argc, char **argv)
{   if (argc != 7)	{	fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -n 4 ./Run_EQsim InputName RealizationNumber GridNumber(Faults) WhichDeltT2Use UsePostSeis EQrecordLength\n");		exit(10);			}
	/*-------------------------------------------------------------------*/
    unsigned int	i,				j,              k;
    int	    		iRANK,         	iSIZE;
    /*-------------------------------------------------------------------*/
    MPI_Init(&argc, &argv);               /* start MPI processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &iRANK);/* current proc. id    */
    MPI_Comm_size(MPI_COMM_WORLD, &iSIZE);/* # of processes      */ 

    MPI_Status status;
    MPI_Offset offset1; /*this stuff here is for file output*/
    MPI_Offset offset2; /*this stuff here is for file output*/
    /*-------------------------------------------------------------------*/
    int 	iOFFSET_F[iSIZE];	for (i = 0; i < iSIZE; i++)             		{       iOFFSET_F[i]    = 0;                              }
    int 	iSTARTPOS_F[iSIZE];     for (i = 0; i < iSIZE; i++)         			{       iSTARTPOS_F[i]  = 0;                              }
    int 	iOFFSET_B[iSIZE];		for (i = 0; i < iSIZE; i++)             	{       iOFFSET_B[i]    = 0;                              }
    int 	iSTARTPOS_B[iSIZE];     for (i = 0; i < iSIZE; i++)         			{       iSTARTPOS_B[i]  = 0;                              }
    /*-------------------------------------------------------------------*/
    unsigned int iUsdGrid,	     iWrteKinMod,     iDeltT2Use,       iUsePostSeismic,        iRealizNum;
    float        fKinModMinMag,      fRecLgth;
    char         *retch;
    /*----------------------------*/
    sscanf(argv[2],"%u",  &iRealizNum);
  	sscanf(argv[3],"%u",  &iUsdGrid);    						sscanf(argv[4],"%u",  &iDeltT2Use);	                        
  	sscanf(argv[5],"%u",  &iUsePostSeismic);			        sscanf(argv[6],"%f",  &fRecLgth); 			
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
    unsigned int    iFltSegmNum,	iFltGridNum,	    iFltPtchNum,                iFltVertNum; 
    unsigned int    iBndSegmNum,        iBndGridNum,        iBndPtchNum,                iBndVertNum;
    unsigned int    iCmbPtchNum,        iCmbVertNum;
    float  	    fMeanLegLgth,	fMeanBndLegLgth,    fDummy;
    float           *fSD_SlipRate,      *fSD_SlipRake,      *fSD_StrssRate;
  	
  	char			ctempVals[512],		cFile1_In[512],		cFile2_In[512],             cFile1_Out[512];   /* out1 is for the overall catalog; out2 is the stf for each patch that had slip => both are binaries!!!  	*/
	FILE    *fpIn;  
	/*----------------------------*/
    strcpy(cFile1_In,  argv[1]);              			   		strcat(cFile1_In,"_FLT.txt");
    strcpy(cFile2_In,  argv[1]);              			   		strcat(cFile2_In,"_BND.txt");
    strcpy(cFile1_Out, argv[1]);              					strcat(cFile1_Out,"_");				strcat(cFile1_Out,argv[2]);		strcat(cFile1_Out,"_Catalog.dat");
    /*----------------------------------------------------------------------------------*/  
    if ((fpIn = fopen(cFile1_In,"r")) == NULL)       		{   fprintf(stdout,"Error -cant open %s file. in Main function \n",cFile1_In);      exit(10);     }
	retch =fgets(ctempVals, 512, fpIn);							sscanf(ctempVals,"%*s %d",  &iFltSegmNum); 
    retch =fgets(ctempVals, 512, fpIn);          				sscanf(ctempVals,"%*s %d",  &iFltGridNum);    
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
        for (j = 0; j < iFltGridNum; j++)			   	 {   retch =fgets(ctempVals, 512, fpIn);                          retch =fgets(ctempVals, 512, fpIn);           } 
    }
    fclose(fpIn);
    /*----------------------------------------------------------------------------------*/  
    if ((fpIn = fopen(cFile2_In,"r")) == NULL)       		
    {   fprintf(stdout,"_BND.txt file was not found i.e., opened; continue without BC faults.... \n");  
        iBndSegmNum = 0;                                         iBndGridNum = 0;
        iBndPtchNum = 0;                                         iBndVertNum = 0;
        fMeanBndLegLgth = 0.0;
    }
	else
	{   retch =fgets(ctempVals, 512, fpIn);					 sscanf(ctempVals,"%*s %d",  &iBndSegmNum); 
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %d",  &iBndGridNum); 
        if (iBndGridNum < 1)                            {   fprintf(stdout,"Error -boundary fault file was opened but BC faults are not meshed/gridded. \n");      exit(10);     }
	    retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %d  %d", &iBndPtchNum,     &iBndVertNum); 
        retch =fgets(ctempVals, 512, fpIn);                  sscanf(ctempVals,"%*s %e  %e", &fMeanBndLegLgth, &fDummy); 
	}
	fMeanBndLegLgth *= 1.0E+3; /* now it is in meters */
	fclose(fpIn);
	/*----------------------------------------------------------------------------------*/  
	fMeanBndLegLgth = (fMeanBndLegLgth <=      0.0    ) ? fMeanLegLgth : fMeanBndLegLgth; /*set BndLegLength to MeanLegLgth if BndLegLength is equal or smaller than 0.0*/
	fMeanLegLgth    = (fMeanBndLegLgth >= fMeanLegLgth) ? fMeanLegLgth : fMeanBndLegLgth; /*set MeanLegLength to be min value of LegLength and BndLegLength*/
    /*----------------------------------------------------------------------------------*/  
    iCmbPtchNum     = iFltPtchNum + iBndPtchNum;
    iCmbVertNum     = iFltVertNum + iBndVertNum;
    /*----------------------------------------------------------------------------------*/  
    unsigned int    iBASEelem_F;		iBASEelem_F = (unsigned int)(iFltPtchNum/iSIZE); 
    unsigned int    iADDelem_F;			iADDelem_F  = (unsigned int)(iFltPtchNum%iSIZE);
    unsigned int    iBASEelem_B;	    iBASEelem_B = (unsigned int)(iBndPtchNum/iSIZE); 
    unsigned int    iADDelem_B;			iADDelem_B  = (unsigned int)(iBndPtchNum%iSIZE);
    /*----------------------------*/
    for (i = 0; i < iSIZE;      i++)          				{ 	iOFFSET_F[i]     = iBASEelem_F;                                         }
    for (i = 0; i < iADDelem_F; i++)       					{   iOFFSET_F[i]    += 1;                                                   }   
    for (i = 1; i < iSIZE;      i++)          				{   iSTARTPOS_F[i]   = iSTARTPOS_F[i-1] + iOFFSET_F[i-1];                   }
    for (i = 0; i < iSIZE;      i++)          				{ 	iOFFSET_B[i]     = iBASEelem_B;                                         }
    for (i = 0; i < iADDelem_B; i++)       					{   iOFFSET_B[i]    += 1;                                                   }   
    for (i = 1; i < iSIZE;      i++)          				{   iSTARTPOS_B[i]   = iSTARTPOS_B[i-1] + iOFFSET_B[i-1];                   }    

   // fprintf(stdout,"My Rank: %d  TriNum %d    BASEelem  %d   ADDelem  %d  STARTPOS  %d OFFSET   %d; start: %d     end: %d\n",iRANK, iCmbPtchNum, iBASEelem_A, iADDelem_A, iSTARTPOS_A[iRANK], iOFFSET_A[iRANK], iSTARTPOS_A[iRANK],  (iSTARTPOS_A[iRANK]+iOFFSET_A[iRANK])  );
    /*-------------------------------------------------------------------*/ 
	float           *fMD_AddNrmStrss;       fMD_AddNrmStrss    = (float *)        calloc(1, sizeof(float));    
    float           *fMD_VpVsRatio;         fMD_VpVsRatio      = (float *)        calloc(1, sizeof(float));	
    float	        *fMD_ShearMod;          fMD_ShearMod       = (float *)        calloc(1, sizeof(float));		
    float           *fMD_Vp;                fMD_Vp             = (float *)        calloc(1, sizeof(float));	   
    float           *fMD_Vs;                fMD_Vs             = (float *)        calloc(1, sizeof(float));	
    float	        *fMD_Poisson;           fMD_Poisson        = (float *)        calloc(1, sizeof(float));						
    float	        *fMD_Lambda;            fMD_Lambda         = (float *)        calloc(1, sizeof(float));						             	
    float		    *fMD_MedDense;          fMD_MedDense       = (float *)        calloc(1, sizeof(float));	
    float           *fMD_CritSlipVelo;      fMD_CritSlipVelo   = (float *)        calloc(1, sizeof(float));	  
	unsigned int	*iMD_ChgFricBtwEQs;     iMD_ChgFricBtwEQs  = (unsigned int *) calloc(1, sizeof(unsigned int));  
	
    unsigned int 	*iSD_FricLawUSED;       iSD_FricLawUSED    = (unsigned int *) calloc(iFltSegmNum, sizeof(unsigned int));	
    float           *fSD_RefStatFric;       fSD_RefStatFric    = (float *)        calloc(iFltSegmNum, sizeof(float));                       
    float           *fSD_RefStatFr_vari;    fSD_RefStatFr_vari = (float *)        calloc(iFltSegmNum, sizeof(float));		
    float 	        *fSD_RefDynFric;        fSD_RefDynFric     = (float *)        calloc(iFltSegmNum, sizeof(float));                       
    float 	        *fSD_RefDynFr_vari;     fSD_RefDynFr_vari  = (float *)        calloc(iFltSegmNum, sizeof(float));		
	float 	        *fSD_CritSlipDist;      fSD_CritSlipDist   = (float *)        calloc(iFltSegmNum, sizeof(float));                       
	float 	        *fSD_CritSlipD_vari;    fSD_CritSlipD_vari = (float *)        calloc(iFltSegmNum, sizeof(float));		    				   

	unsigned int	*iTDg_SegID;             iTDg_SegID          = (unsigned int *) calloc(iFltPtchNum,   sizeof(unsigned int));
	float           *fTDg_SlipRate;          fTDg_SlipRate       = (float *)        calloc(iFltPtchNum,   sizeof(float));                  
	float           *fTDg_SlipRake;          fTDg_SlipRake       = (float *)        calloc(iFltPtchNum,   sizeof(float));
	float           *fTDg_StrsRteStk;        fTDg_StrsRteStk     = (float *)        calloc(iFltPtchNum,   sizeof(float));                  
	float           *fTDg_StrsRteDip;        fTDg_StrsRteDip     = (float *)        calloc(iFltPtchNum,   sizeof(float));			
    float 		    *fTDg_StatFric;          fTDg_StatFric       = (float *)        calloc(iFltPtchNum,   sizeof(float));                       
	float 	        *fTDg_DynFric;           fTDg_DynFric        = (float *)        calloc(iFltPtchNum,   sizeof(float));	
    unsigned int	*iTDg_StabType;          iTDg_StabType       = (unsigned int *) calloc(iFltPtchNum,   sizeof(unsigned int));    
	
	float           *fTDg_CentEpos;          fTDg_CentEpos       = (float *)        calloc((iFltPtchNum + iBndPtchNum),   sizeof(float));	                    
	float 	        *fTDg_CentNpos;          fTDg_CentNpos       = (float *)        calloc((iFltPtchNum + iBndPtchNum),   sizeof(float));					
	float 	        *fTDg_CentZpos;          fTDg_CentZpos       = (float *)        calloc((iFltPtchNum + iBndPtchNum),   sizeof(float));
	      
    unsigned int	*iTDg_V1;                iTDg_V1	         = (unsigned int *) calloc((iFltPtchNum + iBndPtchNum),   sizeof(unsigned int));      /*those are "T" b/c they will be freed later on -> only used during model setup e.g., calculation of K-matrix*/   
    unsigned int	*iTDg_V2;                iTDg_V2             = (unsigned int *) calloc((iFltPtchNum + iBndPtchNum),   sizeof(unsigned int));    	   
    unsigned int    *iTDg_V3;                iTDg_V3 	 	     = (unsigned int *) calloc((iFltPtchNum + iBndPtchNum),   sizeof(unsigned int));  
						
	float           *fVDg_Epos;              fVDg_Epos           = (float *)        calloc((iFltVertNum + iBndVertNum),  sizeof(float));	                    
	float           *fVDg_Npos;              fVDg_Npos           = (float *)        calloc((iFltVertNum + iBndVertNum),  sizeof(float));
	float 	        *fVDg_Zpos;              fVDg_Zpos           = (float *)        calloc((iFltVertNum + iBndVertNum),  sizeof(float));	
    /*-----------------------------------*/
	if (iRANK == 0) 
    {	
    
        LoadInputParameter(argv, iRealizNum, iUsdGrid, iFltSegmNum, iFltPtchNum, iBndPtchNum, iFltVertNum, iBndVertNum, fSD_StrssRate, fSD_SlipRate, fSD_SlipRake, iMD_ChgFricBtwEQs, fMD_AddNrmStrss, fMD_Vp, fMD_Vs, fMD_Poisson, fMD_Lambda, fMD_ShearMod, fMD_MedDense, fMD_CritSlipVelo, iSD_FricLawUSED, fSD_RefStatFric, fSD_RefStatFr_vari, fSD_RefDynFric, fSD_RefDynFr_vari, fSD_CritSlipDist, fSD_CritSlipD_vari, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, iTDg_V1, iTDg_V2, iTDg_V3, iTDg_StabType, iTDg_SegID, fTDg_SlipRate, fTDg_SlipRake, fTDg_StrsRteStk, fTDg_StrsRteDip, fTDg_StatFric, fTDg_DynFric, fVDg_Epos, fVDg_Npos, fVDg_Zpos);						

	}
    /*-----------------------------------*/
    float           fIntSeisTimeStp = 2.0/365.25; /* == 1 is one-day steps; is fraction of a year*/
    float           fTimeStpInSecs  = fIntSeisTimeStp*31536000.0; /* now the time is in seconds */
    unsigned int    *iTDl_StabType;             iTDl_StabType        = (unsigned int *) calloc(iOFFSET_F[iRANK],   sizeof(unsigned int));
    unsigned int    *iTDl_StabTypeChgd;         iTDl_StabTypeChgd    = (unsigned int *) calloc(iOFFSET_F[iRANK],   sizeof(unsigned int)); 
    unsigned int    *iTDl_SegID;                iTDl_SegID           = (unsigned int *) calloc(iOFFSET_F[iRANK],   sizeof(unsigned int));
    float           *fTDl_SlipRate;             fTDl_SlipRate        = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_SlipRake;             fTDl_SlipRake        = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));   
    float           *fTDl_RefStrssRateStk;      fTDl_RefStrssRateStk = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float)); /*Ref and CurStrss rate can be different b/c of stability type changes/effects...*/
    float           *fTDl_RefStrssRateDip;      fTDl_RefStrssRateDip = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_CurStrssRateStk;      fTDl_CurStrssRateStk = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_CurStrssRateDip;      fTDl_CurStrssRateDip = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));         
    float           *fTDl_StatFric;             fTDl_StatFric        = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_CurrFric;             fTDl_CurrFric        = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_DynFric;              fTDl_DynFric         = (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));
    float           *fTDl_PostStrssat_t0_F;     fTDl_PostStrssat_t0_F= (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float));    
    float           *fTDl_PostStrssTStep_F;     fTDl_PostStrssTStep_F= (float *)        calloc(iOFFSET_F[iRANK],   sizeof(float)); 
    float           *fTDl_PostStrssat_t0_B;     fTDl_PostStrssat_t0_B= (float *)        calloc(iOFFSET_B[iRANK],   sizeof(float));    
    float           *fTDl_PostStrssTStep_B;     fTDl_PostStrssTStep_B= (float *)        calloc(iOFFSET_B[iRANK],   sizeof(float)); 


    fMD_VpVsRatio[0] = fMD_Vp[0]/fMD_Vs[0];
    /*----------------------------*/
    MPI_Bcast(iMD_ChgFricBtwEQs,           1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(fMD_AddNrmStrss,             1, MPI_FLOAT,    0, MPI_COMM_WORLD);	
	MPI_Bcast(iSD_FricLawUSED,   iFltSegmNum, MPI_UNSIGNED, 0, MPI_COMM_WORLD);			

	MPI_Bcast(fMD_Vp,                         1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(fMD_Vs,                         1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(fMD_VpVsRatio,                  1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(fMD_Poisson,                    1, MPI_FLOAT, 0, MPI_COMM_WORLD);	
    MPI_Bcast(fMD_Lambda,                     1, MPI_FLOAT, 0, MPI_COMM_WORLD);                   
    MPI_Bcast(fMD_ShearMod,                   1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fMD_MedDense,                   1, MPI_FLOAT, 0, MPI_COMM_WORLD);                    		
	
    MPI_Bcast(fSD_RefStatFric,      iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	                
    MPI_Bcast(fSD_RefStatFr_vari,   iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
    MPI_Bcast(fSD_RefDynFric,       iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);                   
    MPI_Bcast(fSD_RefDynFr_vari,    iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fSD_CritSlipDist,     iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);			        
    MPI_Bcast(fSD_CritSlipD_vari,   iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
    MPI_Bcast(fSD_CritSlipVelo,     iFltSegmNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
	
	MPI_Bcast(iTDg_V1,    (iFltPtchNum + iBndPtchNum), MPI_UNSIGNED, 0, MPI_COMM_WORLD);			
	MPI_Bcast(iTDg_V2,    (iFltPtchNum + iBndPtchNum), MPI_UNSIGNED, 0, MPI_COMM_WORLD);			
    MPI_Bcast(iTDg_V3,    (iFltPtchNum + iBndPtchNum), MPI_UNSIGNED, 0, MPI_COMM_WORLD);	   
	MPI_Bcast(fTDg_CentEpos, (iFltPtchNum + iBndPtchNum), MPI_FLOAT, 0, MPI_COMM_WORLD);			        
	MPI_Bcast(fTDg_CentNpos, (iFltPtchNum + iBndPtchNum), MPI_FLOAT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(fTDg_CentZpos, (iFltPtchNum + iBndPtchNum), MPI_FLOAT, 0, MPI_COMM_WORLD);	
	
    MPI_Bcast(fVDg_Epos,     (iFltVertNum + iBndVertNum), MPI_FLOAT, 0, MPI_COMM_WORLD);	                
	MPI_Bcast(fVDg_Npos,     (iFltVertNum + iBndVertNum), MPI_FLOAT, 0, MPI_COMM_WORLD);	        
	MPI_Bcast(fVDg_Zpos,     (iFltVertNum + iBndVertNum), MPI_FLOAT, 0, MPI_COMM_WORLD);			
	
	MPI_Scatterv(iTDg_StabType,     iOFFSET_F, iSTARTPOS_F, MPI_UNSIGNED, iTDl_StabType,     iOFFSET_F[iRANK], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Scatterv(iTDg_SegID,        iOFFSET_F, iSTARTPOS_F, MPI_UNSIGNED, iTDl_SegID,        iOFFSET_F[iRANK], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
    MPI_Scatterv(fTDg_SlipRate,     iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_SlipRate,        iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	MPI_Scatterv(fTDg_SlipRake,     iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_SlipRake,        iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	MPI_Scatterv(fTDg_StrsRteStk,   iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_RefStrssRateStk, iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	MPI_Scatterv(fTDg_StrsRteDip,   iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_RefStrssRateDip, iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);

	MPI_Scatterv(fTDg_StatFric,     iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_StatFric,        iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD); 
	MPI_Scatterv(fTDg_StatFric,     iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_CurrFric,        iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD); 
	MPI_Scatterv(fTDg_DynFric,      iOFFSET_F, iSTARTPOS_F, MPI_FLOAT, fTDl_DynFric,         iOFFSET_F[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	/*-------------------------------------------------------------------*/    
	float fdeltTincr; /* this one is for time steps during rupture => depends on distance of patches... but this needs calibration */
	  
	if (iDeltT2Use == 1)    {	fdeltTincr = 2.0*fMeanLegLgth/fMD_Vp[0];           fdeltTincr = floor(fdeltTincr*10.0)/10.0;            }
	else                    {	fdeltTincr = 1.0E+09;                                                                                   }
	
	if (iRANK == 0)		{		fprintf(stdout,"delta T (coseismic) %f   fMeanLegLgth %f  RecLength %f\n",fdeltTincr,fMeanLegLgth,fRecLgth);			}
    /*-------------------------------------------------------------------*/
	unsigned int    iRandNumber = 1000000; 
	float           fMD_g      = 9.81;
	float           fViscosity  = 1.0E+13; /* in Pa seconds; value doesn't seem right, but need something that corresponds to patch stiffness => to give proper char. time scale for post-seismic decay*/
	long            *iSeedVal;          iSeedVal          = (long *)         calloc(       1   ,                  sizeof(long)); 
	
	unsigned int    *iTDlg_TravTimes;       iTDlg_TravTimes     = (unsigned int *) calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(unsigned int));
    float           *fTDlg_LocSrcRcv_H;     fTDlg_LocSrcRcv_H   = (float *)        calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float));
    float           *fTDlg_LocSrcRcv_V;     fTDlg_LocSrcRcv_V   = (float *)        calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float)); 
    float           *fTDlg_LocSrcRcv_N;     fTDlg_LocSrcRcv_N   = (float *)        calloc(iOFFSET_F[iRANK]*iFltPtchNum,   sizeof(float));  
    unsigned int    *iTDg_GlobTTmax;        iTDg_GlobTTmax      = (unsigned int *) calloc(       1   ,                    sizeof(unsigned int));
	float           *fTDl_RefNormStrss;     fTDl_RefNormStrss   = (float *)        calloc(iOFFSET_F[iRANK],               sizeof(float));
    float           *fTDl_Area;             fTDl_Area           = (float *)        calloc(iOFFSET_F[iRANK],               sizeof(float));
    float           *fTDl_Curr_DcVal;       fTDl_Curr_DcVal     = (float *)        calloc(iOFFSET_F[iRANK],               sizeof(float));

	float           *fRandVector;           fRandVector         = (float *)        calloc(iRandNumber,                    sizeof(float));
	unsigned int    *iRandPos;              iRandPos            = (unsigned int *) calloc(       1   ,                    sizeof(unsigned int)); 
	  
	time_t  t;
    srand((unsigned) time(&t)+iRANK*10);  /* initialize the random number generator; will generate 1E+6 random numbers => should be enough*/
	iSeedVal[0] = rand();
	for (i = 0; i < iRandNumber; i++ )       {           fRandVector[i] = ran0_inmain(iSeedVal);       } /* using a generator as described in "Numerical Recipes, page 279*/
	/*-------------------------------------------------------------------*/
	
    DefineMoreParas(iRANK, iSTARTPOS_F, iOFFSET_F, iSTARTPOS_B, iOFFSET_B, iFltPtchNum, iBndPtchNum, iTDg_V1, iTDg_V2, iTDg_V3, iTDl_SegID, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fVDg_Epos, fVDg_Npos, fVDg_Zpos, fSD_CritSlipDist, fSD_CritSlipD_vari, fMD_Vp, fMD_VpVsRatio, fMD_MedDense, fMD_AddNrmStrss, fMD_g, fdeltTincr, fRandVector, iRandNumber, iRandPos, iTDlg_TravTimes, iTDg_GlobTTmax, fTDl_Area, fTDl_RefNormStrss, fTDl_Curr_DcVal, fTDlg_LocSrcRcv_H, fTDlg_LocSrcRcv_V, fTDlg_LocSrcRcv_N);

    MPI_Allreduce(MPI_IN_PLACE, iTDg_GlobTTmax, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD); 
    
	/*-------------------------------------------------------------------*/ 
	/* once I do include postseismic deformation or after slip, I'll have to add the additional patches here => iPatchNum +iBoundPatchNum */
	int          iVectPos;
	unsigned int iMaxSTFLength = iTDg_GlobTTmax[0] +1;
	float        fMD_UnitSlip = 1.0E-4*fMeanLegLgth; /* slip in meter => if meanleglength is 1000m, then a slip of .1m is used */
	float        *fK_FF_SS;                fK_FF_SS           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	        
	float        *fK_FF_SD;                fK_FF_SD           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));					
    float        *fK_FF_SO;                fK_FF_SO           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	        
    float        *fK_FF_DS;                fK_FF_DS           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	   
    float        *fK_FF_DD;                fK_FF_DD           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));						   
    float        *fK_FF_DO;                fK_FF_DO           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	
    float        *fK_FF_OS;                fK_FF_OS           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	   
    float        *fK_FF_OD;                fK_FF_OD           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));						   
    float        *fK_FF_OO;                fK_FF_OO           = (float *)          calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	
    /*for fault-fault i allow change in tensile, but no slip in that direction, change in tensile(normal comp.) will clamp/unclamp the fault...
      that change in normal stress will only last for the event and then instantly revert back to reference value; better would be to 
      have it revet back to reference value in a non-instant way -like diffusing it away... -will do that maybe later; for now it will jump directly back*/
    float        *fK_FB_SS;                fK_FB_SS           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));	        
	float        *fK_FB_SD;                fK_FB_SD           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));					
    float        *fK_FB_SO;                fK_FB_SO           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));	        
    float        *fK_FB_DS;                fK_FB_DS           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));	   
    float        *fK_FB_DD;                fK_FB_DD           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));						   
    float        *fK_FB_DO;                fK_FB_DO           = (float *)          calloc(iOFFSET_F[iRANK]*iBndPtchNum,sizeof(float));	
    /*as the comment above, no opening/closing motion on the faults...*/
    float        *fK_BF_SS;                fK_BF_SS           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));	        
	float        *fK_BF_SD;                fK_BF_SD           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));					
    float        *fK_BF_DS;                fK_BF_DS           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));	   
    float        *fK_BF_DD;                fK_BF_DD           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));						   
    float        *fK_BF_OS;                fK_BF_OS           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));						
    float        *fK_BF_OD;                fK_BF_OD           = (float *)          calloc(iOFFSET_B[iRANK]*iFltPtchNum,sizeof(float));	        
    /*when boundary faults are slipping to release postseismic signal, they are not allowed to modify the normal stress on patch...
      that is mainly b/c I have no method right now to bring it back to its reference value... meaning it would stay at higher/lower value
      and doing this multiple times might totally clamp or open a fault, not realistic; if a time-diffusion of normal stress variation would
      be included, then this might be handeled... but will not happen now... so, i allow slip in all directions on B but only look at induced values in strike and dip direction*/      
    float        *fK_BB_SS;                fK_BF_SS           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	        
	float        *fK_BB_SD;                fK_BF_SD           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	
	float        *fK_BB_SO;                fK_BF_SO           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));					
    float        *fK_BB_DS;                fK_BF_DS           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	   
    float        *fK_BB_DD;                fK_BF_DD           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	
    float        *fK_BB_DO;                fK_BF_DO           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));
    float        *fK_BB_OS;                fK_BF_OS           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	   
    float        *fK_BB_OD;                fK_BF_OD           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));	
    float        *fK_BB_OO;                fK_BF_OO           = (float *)          calloc(iOFFSET_B[iRANK]*iBndPtchNum,sizeof(float));							   
   /*for interaction of boundary faults (only needed to get stressing rates when fault-slip-rates have been defined) I want to determine the slip
     distribution that is associated with the prescribed stressing on boundary (which comes from slip on faults), since the boundaries are
     frictionless there is no clamping or unclamping, the opening-mode stress would have to be removed by opening mode motion, so, all components
     are used here */
   
    unsigned int *iTDlg_Pdone;         iTDlg_Pdone      = (unsigned int *) calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	
    unsigned int *iTDlg_Sdone;         iTDlg_Sdone      = (unsigned int *) calloc(iOFFSET_F[iRANK]*iFltPtchNum,sizeof(float));	
    /*these two tell me the index for source-receiver pair => how much of STF at source has been "seen" at receiver already...*/
    unsigned int *iTempValue;          iTempValue       = (unsigned int *) calloc(1,                sizeof(unsigned int)); 
    unsigned int *iTDl_SelfLoc_F;      iTDl_SelfLoc_F   = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int)); 	
    unsigned int *iTDl_SelfLoc_B;      iTDl_SelfLoc_B   = (unsigned int *) calloc(iOFFSET_B[iRANK], sizeof(unsigned int)); 	
    		
    
    if (iRANK == 0)			{ 		fprintf(stdout,"Building K-matrix....\n");			}
	/*-------------------------------------------------------------------*/
	
    Build_K_Matrix(iRANK, iSTARTPOS_F, OFFSET_F, iSTARTPOS_B, iOFFSET_B, iFltPtchNum, iBndPtchNum, fMD_UnitSlip, fMD_ShearMod, fMD_Lambda, iTDg_V1, iTDg_V2, iTDg_V3, fVDg_Epos, fVDg_Npos, fVDg_Zpos, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fTDl_Curr_DcVal, fTDl_StatFric, fTDl_DynFric, fTDl_RefNormStrss, fTDl_Area, fTDl_RotMatF, fTDl_RotMatB, fK_FF_SS, fK_FF_SD, fK_FF_SO, fK_FF_DS, fK_FF_DD, fK_FF_DO, fK_FF_OS, fK_FF_OD, fK_FF_OO, fK_FB_SS, fK_FB_SD, fK_FB_SO, fK_FB_DS, fK_FB_DD, fK_FB_DO, fK_BF_SS, fK_BF_SD, fK_BF_DS, fK_BF_DD, K_BF_OS, fK_BF_OD, fK_BB_SS, fK_BB_SD, fK_BB_SO, fK_BB_DS, fK_BB_DD, fK_BB_DO, fK_BB_OS, fK_BB_OD, fK_BB_OO, iTDl_SelfLoc_F, iTDl_SelfLoc_B, iTDl_StabType);

    if (iRANK == 0)			{ 		fprintf(stdout,"K-matrix done \n");			        }	
	/*-------------------------------------------------------------------*/  
	float *fTDg_StrssRteChgStk;     fTDg_StrssRteChgStk   = (float *) calloc(iFltPtchNum,      sizeof(float)); 
	float *fTDg_StrssRteChgDip;     fTDg_StrssRteChgDip   = (float *) calloc(iFltPtchNum,      sizeof(float)); 
	float *fTDl_CharRlxTime_F;      fTDl_CharRlxTime_F    = (float *) calloc(iOFFSET_F[iRANK], sizeof(float )); 		
	float *fTDl_CharRlxTime_B;      fTDl_CharRlxTime_B    = (float *) calloc(iOFFSET_B[iRANK], sizeof(float )); 
	float fTempSlipStk,             fTempSlipDip;
	/*-------------------------------------------------------------------*/   
    for (i = iOFFSET_F[iRANK]; i--;  )           {       if (fabs(fTDl_SlipRate[i]) > 0.0)       {   iTempValue[0] = 1;    }                    }
    MPI_Allreduce(MPI_IN_PLACE iTempValue, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
   
    if (iTempValue[0] == 1)
    {    
      
      
      
      
      
      
       // void  GetTheStressingRates(
       //the idea is that I now know the RefStrssRateStk and Dip; 
       //use slip on faults and interaction matrix (the "_FB_", "_BB_", and "_BF_") to get first the induced stress at boundary faults, then release via slip and then use that slip
       //for loading on interaction faults; this is added to the existing stressing rate
 
 
    //i have this in some matlab function... should be able to more or less copy it out... could also "limit" the number of iterations; although this part should be reasonably fast anyways...
    
 
 
 
 
 
 
    }
    /*-------------------------------------------------------------------*/   
    for (i = iOFFSET_F[iRANK]; i--;   )
    {   if (iTDl_StabType[i] > 1 )/* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
    	{	
    	    fTempSlipStk =  fTDl_RefStrssRateStk[i]/fK_FF_SS[iTDl_SelfLoc_F[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
    		fTempSlipDip =  fTDl_RefStrssRateDip[i]/fK_FF_DD[iTDl_SelfLoc_F[i]];
    
    		for (j = iFltPtchNum; j--;    )
    		{	if (iTDl_StabType[j] == 1) /* apply the creep slip to load the  unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */	
    			{   iVectPos         = i*iFltPtchNum + j;  				
    				fTDg_StrssRteChgStk[j] += fTempSlipStk*fK_FF_SS[iVectPos] + fTempSlipDip*fK_FF_DS[iVectPos];
    				fTDg_StrssRteChgDip[j] += fTempSlipStk*fK_FF_SD[iVectPos] + fTempSlipDip*fK_FF_DD[iVectPos];	
	        }   }   
	    } 
	    fTDl_CharRlxTime_F[i] = fabs(fViscosity/(fK_FF_SS[iTDl_SelfLoc_F[i]] > fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]])); 
	}
	for (i = iOFFSET_B[iRANK]; i--;   )
    {   fTDl_CharRlxTime_B[i] = fabs(fViscosity/(fK_BB_SS[iTDl_SelfLoc_B[i]] > fK_BB_DD[iTDl_SelfLoc_B[i]] ? fK_BB_SS[iTDl_SelfLoc_B[i]] : fK_BB_DD[iTDl_SelfLoc_B[i]])); 
    }
    /*-------------------------------------------------------------------*/
    MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgStk, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
    MPI_Allreduce(MPI_IN_PLACE, fTDg_StrssRteChgDip, iFltPtchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);    
    /*-------------------------------------------------------------------*/
	for (i = iOFFSET_F[iRANK]; i--;    )
    {	
        if (iTDl_StabType[i] == 1) /* if it is an unstable atch/cell... =>apply the creep slip to load */
    	{	fTDl_CurStrssRateStk[i]  = fTDl_RefStrssRateStk[i] +fTDg_StrssRteChgStk[i+iSTARTPOS_F[iRANK]];			
    		fTDl_CurStrssRateDip[i]  = fTDl_RefStrssRateDip[i] +fTDg_StrssRteChgDip[i+iSTARTPOS_F[iRANK]];	      
	    }   
	    else
	    {   fTDl_CurStrssRateStk[i]  = 0.0;
	        fTDl_CurStrssRateDip[i]  = 0.0;
	}   }	
    /*-------------------------------------------------------------------*/
	/* now the main preparation work is done, define remaining variables before going into the EQ iteration loops*/    
	unsigned int    iwritePos,         ireadPos; 
	float           fTimeYears;
	unsigned int    iEQcounter     = 0;
	unsigned int    iMaxMRFlength  = 10000;
    float           fMD_CutStrss   = 2000.0; /* stress in Pa ==> earth tides change stress by about 5000Pa (= 0.005MPa); atmospheric pressure is ~100kPa => 2000Pa == 2% of atmospheric pressure change */
    //this "CutStrss is also a way to stop rupture i.e., make re-rupture more difficult... "healing" => making it higher could be good to ensure that the iterations don't become extreme long and time consuming
    //another way could be to link this again to Dc value
    //leave the value as is for now... maybe add another variable to do the "healing"        
    float         fTemp,              fTemp1,                 fTemp2,             fTemp3;
    float         fStkSlip,           fDipSlip;
    float         fPostSeisIncrStk,   fPostSeisIncrDip,       fPostSeisIncrNrm;               
    /*---------------------------------------*/ 				
	unsigned int *iMD_EQongoing;      iMD_EQongoing    = (unsigned int *) calloc(1 ,sizeof(unsigned int));	
	unsigned int *iEndOfEventCntr;    iEndOfEventCntr  = (unsigned int *) calloc(1 ,sizeof(unsigned int));	
	unsigned int *iEQg_ActPtchNum;	  iEQg_ActPtchNum  = (unsigned int *) calloc(1, sizeof(unsigned int));	
	unsigned int *iEQg_MRFlength;     iEQg_MRFlength   = (unsigned int *) calloc(1, sizeof(unsigned int));  
    unsigned int *iMD_TotalRuptT;     iMD_TotalRuptT  = (unsigned int *) calloc(1, sizeof(unsigned int));	
    unsigned int *iMD_ChangedStab;    iMD_ChangedStab  = (unsigned int *) calloc(1, sizeof(unsigned int));
    	
	float        *fEQg_SeisPot;       fEQg_SeisPot     = (float *)        calloc(1, sizeof(float));
    
    unsigned int *iEQg_WrtStartPos;   iEQg_WrtStartPos = (unsigned int *) calloc(iSIZE, sizeof(unsigned int));
     
    unsigned int *iTDl_WasActivated;  iTDl_WasActivated= (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int)); 
	unsigned int *iTDl_t0;            iTDl_t0          = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int)); 		  
	unsigned int *iTDl_STFcnt;        iTDl_STFcnt      = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int)); 			
    unsigned int *iTDl_PosOfLastSlp;  iTDl_PosOfLastSlp= (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int));  
	unsigned int *iEQl_ActPtchID;     iEQl_ActPtchID   = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int));					
	unsigned int *iEQl_t0ofPtch;      iEQl_t0ofPtch	   = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int));	
	unsigned int *iEQl_StabType;      iEQl_StabType    = (unsigned int *) calloc(iOFFSET_F[iRANK], sizeof(unsigned int));

	float *fTDl_StrssB4_H;           fTDl_StrssB4_H    = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 				
	float *fTDl_StrssB4_V;           fTDl_StrssB4_V    = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));     
	float *fTDl_StrssB4_N;           fTDl_StrssB4_N    = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 
	float *fTDl_CurStrss_H;          fTDl_CurStrss_H   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 				
	float *fTDl_CurStrss_V;          fTDl_CurStrss_V   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));     
	float *fTDl_CurStrss_N;          fTDl_CurStrss_N   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 	
    float *fTDl_DeltStrssH;          fTDl_DeltStrssH   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 				
	float *fTDl_DeltStrssV;          fTDl_DeltStrssV   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));     
	float *fTDl_DeltStrssN;          fTDl_DeltStrssN   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 
	float *fTDg_DeltStrssH;          fTDg_DeltStrssH   = (float *)         calloc(iFltPtchNum,      sizeof(float)); 				
	float *fTDg_DeltStrssV;          fTDg_DeltStrssV   = (float *)         calloc(iFltPtchNum,      sizeof(float));     
	float *fTDg_DeltStrssN;          fTDg_DeltStrssN   = (float *)         calloc(iFltPtchNum,      sizeof(float)); 		
	float *fg_TempVal1_F;            fg_TempVal1_F     = (float *)         calloc(iFltPtchNum,      sizeof(float)); 
	float *fg_TempVal2_F;            fg_TempVal2_F     = (float *)         calloc(iFltPtchNum,      sizeof(float)); 
	float *fg_TempVal3_F;            fg_TempVal3_F     = (float *)         calloc(iFltPtchNum,      sizeof(float)); 
	
	float *fg_TempVal1_B;            fg_TempVal1_B     = (float *)        calloc(iBndPtchNum,      sizeof(float)); 
	float *fg_TempVal2_B;            fg_TempVal2_B     = (float *)        calloc(iBndPtchNum,      sizeof(float)); 
	
	float *fTDl_StrssOnBnd_H;        fTDl_StrssOnBnd_H = (float *)         calloc(iOFFSET_B[iRANK], sizeof(float));
	float *fTDl_StrssOnBnd_V;        fTDl_StrssOnBnd_V = (float *)         calloc(iOFFSET_B[iRANK], sizeof(float));
	float *fTDl_StrssOnBnd_N;        fTDl_StrssOnBnd_N = (float *)         calloc(iOFFSET_B[iRANK], sizeof(float));
	
	float *fTDl_EventSlipH;          fTDl_EventSlipH   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));     
	float *fTDl_EventSlipV;          fTDl_EventSlipV   = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float)); 
	float *fEQl_SlipHofPtch;         fEQl_SlipHofPtch  = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));    		
	float *fEQl_SlipVofPtch;         fEQl_SlipVofPtch  = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));
	float *fEQl_DeltTofPtch;         fEQl_DeltTofPtch  = (float *)         calloc(iOFFSET_F[iRANK], sizeof(float));
	
	float *fTDlg_STF_H;              fTDlg_STF_H       = (float *)         calloc(iOFFSET_F[iRANK]*iMaxSTFLength, sizeof(float)); 
	float *fTDlg_STF_V;              fTDlg_STF_V       = (float *)         calloc(iOFFSET_F[iRANK]*iMaxSTFLength, sizeof(float)); 
		
    float *fEQl_MRFvals;             fEQl_MRFvals      = (float *)         calloc(iMaxMRFlength, sizeof(float));
    float *fEQg_MRFvals;             fEQg_MRFvals      = (float *)         calloc(iMaxMRFlength, sizeof(float));			
		
    /*-------------------------------------------------------------------*/  
    free(fMD_Vp);           free(fMD_Vs);           free(fMD_Poisson);         
    free(fMD_MedDense);     free(fMD_Lambda);       free(fMD_AddNrmStrss);
    free(fSD_StrssRate);    free(fSD_SlipRate);     free(fSD_SlipRake);
    free(iTDg_SegID);       free(fTDg_SlipRate);    free(fTDg_SlipRake); 
    free(fTDg_StrsRteStk);  free(fTDg_StrsRteDip);  
    free(fTDg_StatFric);    free(fTDg_DynFric);     free(iTDg_StabType);
    
    free(iTDg_V1);          free(iTDg_V2);          free(iTDg_V3); 
    free(fVDg_Epos);        free(fVDg_Npos);        free(fVDg_Zpos); 
    
    free(fK_BB_SS);         free(fK_BB_SD);         free(fK_BB_SO);
    free(fK_BB_DS);         free(fK_BB_DD);         free(fK_BB_DO);
    free(fK_BB_OS);         free(fK_BB_OD);         free(fK_BB_OO);
    free(iTempValue);
    //maybe have to i.e., can free more variabiles...   
    /*-------------------------------------------------------------------*/
    MPI_File fp_MPIout1;
	/*-------------------------------------------------------------------*/
    MPI_File_delete(cFile1_Out, MPI_INFO_NULL);
	MPI_File_open(MPI_COMM_WORLD, cFile1_Out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
	if (iRANK == 0)			
	{	MPI_File_write(fp_MPIout1, &iEQcounter,   1, MPI_UNSIGNED, &status);
		MPI_File_write(fp_MPIout1, &iUsdGrid,     1, MPI_UNSIGNED, &status);			
		MPI_File_write(fp_MPIout1, fMD_ShearMod,  1, MPI_FLOAT,    &status);	
		MPI_File_write(fp_MPIout1, &fdeltTincr,   1, MPI_FLOAT,    &status);	
	}
	MPI_File_close(&fp_MPIout1);	
	
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

    fTimeYears = 0.0;
    while (fTimeYears < fRecLgth)
    {	
		fTimeYears       += fIntSeisTimeStp;
		iMD_EQongoing[0] = 0;
		/*------------------------------------------------------------------------------*/
        if (iUsePostSeismic == 1) /*  https://en.wikipedia.org/wiki/Stress_relaxation    http://web.mit.edu/course/3/3.11/www/modules/visco.pdf      */
		{   /*-------------------------------------------------------------------*/
		    for (j = iFltPtchNum; j--;    )         {    fTDg_DeltStrssH[j]   = 0.0;                fTDg_DeltStrssV[j]   = 0.0;            fTDg_DeltStrssN[j]   = 0.0;             } 
	        /*-------------------------------------------------------------------*/
		    for (i = iOFFSET_F[iRANK]; i--;    ) 
            {	if ((iTDl_StabType[i] > 1) && (fTDl_PostStrssat_t0_F[i] > FLT_EPSILON)) /* don't use 0.0 but epsilon of float!*/
		        {   
		            fTemp1                 = sqrtf(fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);/*that's the currently applied shear stress*/		            
		            fTemp                  = fTemp1 - fTDl_CurrFric[i]*(fTDl_RefNormStrss[i] +fTDl_StrssB4_N[i]); /*difference between applied shear stress and product of friction and normal stress*/
		            fTemp2                 = fTDl_PostStrssat_t0_F[i]*expf(-(fTimeStpInSecs*fTDl_PostStrssTStep_F[i])/fTDl_CharRlxTime_F[i]);
		            fTDl_PostStrssTStep_F[i] += 1.0;
		              
		            if (fTemp > fTemp2) /*if a patch is not locked and if the applied shear stress exceeds the patches strength, then release that amount via slip*/
		            {   fPostSeisIncrStk   = (fTemp-fTemp2)/fTemp *fTDl_StrssB4_H[i]; /*take of a fraction of B4 shear stress, same amount on _H and _V => orientation doesn't change*/
		                fPostSeisIncrDip   = (fTemp-fTemp2)/fTemp *fTDl_StrssB4_V[i];
		                fTemp1             = -1.0*fPostSeisIncrStk/fK_FF_SS[iTDl_SelfLoc_F[i]]; /* slip in horizontal */
		                fTemp2             = -1.0*fPostSeisIncrDip/fK_FF_DD[iTDl_SelfLoc_F[i]]; /* slip in vertical  */
		                
		                for (j = iFltPtchNum; j--;   ) /* here is then the relaxation i.e., stress redistribution according to the amount of loading that will be used this time step */
				        {   /*do I need to test there that only the patches with proper stab. type are getting this?*/
				            iVectPos            = i*iFltPtchNum + j;		
        	                fTDg_DeltStrssH[j] += fTemp1*fK_FF_SS[iVectPos] +fTemp2*fK_FF_DS[iVectPos];
		                    fTDg_DeltStrssV[j] += fTemp1*fK_FF_SD[iVectPos] +fTemp2*fK_FF_DD[iVectPos];
            }   }   }   }
            
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    /*-------------------------------------------------------------------*/  
    	    for (i = iOFFSET_F[iRANK]; i--;    )
            {	fTDl_StrssB4_H[i]   += fTDg_DeltStrssH[i+iSTARTPOS_F[iRANK]];                
                fTDl_StrssB4_V[i]   += fTDg_DeltStrssV[i+iSTARTPOS_F[iRANK]];	
            }      
            /*-------------------------------------------------------------------*/
		    for (j = iFltPtchNum; j--;    )         {    fTDg_DeltStrssH[j]   = 0.0;                fTDg_DeltStrssV[j]   = 0.0;            fTDg_DeltStrssN[j]   = 0.0;             } 
            /*-------------------------------------------------------------------*/
            for (i = iOFFSET_B[iRANK]; i--;    ) //then the same thing for the boundary patches...
            {   
                fTemp                  = sqrtf(fTDl_StrssOnBnd_H[i]*fTDl_StrssOnBnd_H[i] +fTDl_StrssOnBnd_V[i]*fTDl_StrssOnBnd_V[i] +fTDl_StrssOnBnd_N[i]*fTDl_StrssOnBnd_N[i]);
                fTemp2                 = fTDl_PostStrssat_t0_B[i]*expf(-(fTimeStpInSecs*fTDl_PostStrssTStep_B[i])/fTDl_CharRlxTime_B[i]);
                fTDl_PostStrssTStep_B[i] += 1.0;
                
                if (fTemp > fTemp2) /*if a patch is not locked and if the applied shear stress exceeds the patches strength, then release that amount via slip*/
		        {   fPostSeisIncrStk   = (fTemp-fTemp2)/fTemp *fTDl_StrssOnBnd_H[i]; /*take of a fraction of B4 shear stress, same amount on _H and _V => orientation doesn't change*/
		            fPostSeisIncrDip   = (fTemp-fTemp2)/fTemp *fTDl_StrssOnBnd_V[i];
		            fPostSeisIncrNrm   = (fTemp-fTemp2)/fTemp *fTDl_StrssOnBnd_N[i];
		            fTemp1             = -1.0*fPostSeisIncrStk/fK_BB_SS[iTDl_SelfLoc_B[i]]; /* slip in horizontal */
		            fTemp2             = -1.0*fPostSeisIncrDip/fK_BB_DD[iTDl_SelfLoc_B[i]]; /* slip in vertical  */
		            fTemp3             = -1.0*fPostSeisIncrNrm/fK_BB_OO[iTDl_SelfLoc_B[i]]; /* slip in vertical  */
		            
		            fTDl_StrssOnBnd_H[i] -= fPostSeisIncrStk; /*subtract stress on boundary patch; */
		            fTDl_StrssOnBnd_V[i] -= fPostSeisIncrDip;
		            fTDl_StrssOnBnd_N[i] -= fPostSeisIncrNrm;
		                
                    for (j = iFltPtchNum; j--;   ) /* here is then the relaxation i.e., stress redistribution according to the amount of loading that will be used this time step */
				    {   
				        /*do I need to test there that only the patches with proper stab. type are getting this?*/
				        iVectPos      = i*iFltPtchNum + j;		
        	            fTDg_DeltStrssH[j] += fTemp1*fK_BF_SS[iVectPos] +fTemp2*fK_BF_DS[iVectPos] +fTemp3*fK_BF_OS[iVectPos];
		                fTDg_DeltStrssV[j] += fTemp1*fK_BF_SD[iVectPos] +fTemp2*fK_BF_DD[iVectPos] +fTemp3*fK_BF_OD[iVectPos];                
            }   }   }         
            MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    /*-------------------------------------------------------------------*/
    	    for (i = iOFFSET_F[iRANK]; i--;    )
            {	fTDl_StrssB4_H[i]   += fTDg_DeltStrssH[i+iSTARTPOS_F[iRANK]];                
                fTDl_StrssB4_V[i]   += fTDg_DeltStrssV[i+iSTARTPOS_F[iRANK]];	   
        }   }     
		/*------------------------------------------------------------------------------*/
		for (i = iOFFSET_F[iRANK]; i--;    ) 
        {	if (iTDl_StabType[i] == 1 )
        	{   fTDl_StrssB4_H[i] += fTDl_CurStrssRateStk[i]*fIntSeisTimeStp;
    			fTDl_StrssB4_V[i] += fTDl_CurStrssRateDip[i]*fIntSeisTimeStp;
    			fTemp1             = sqrtf (fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);
    			fTemp              = fTemp1 - fTDl_CurrFric[i]*(fTDl_RefNormStrss[i] +fTDl_StrssB4_N[i]); 
    			
       			if (fTemp >= fMD_CutStrss) /*if my excess stress is at least larger than CutStrss value; */
    			{	iMD_EQongoing[0]     = 1;		iTDl_WasActivated[i]  = 1;          iTDl_t0[i]            = 0;	  
        }   }   }
        MPI_Allreduce(MPI_IN_PLACE  iMD_EQongoing, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    	/*-------------------------------------------------------------------*/
    	if (iMD_EQongoing[0] == 1)	/* EARTHQUAKE STARTS */	
    	{   iMD_TotalRuptT[0] = -1;                iEQcounter++;	
    		for (i = iOFFSET_F[iRANK]; i--;    ) 
            {   
                fTDl_CurStrss_H[i] = fTDl_StrssB4_H[i];         fTDl_CurStrss_V[i]   = fTDl_StrssB4_V[i];           fTDl_CurStrss_N[i]   = (fTDl_RefNormStrss[i] +fTDl_StrssB4_N[i]);	    		    
    		    for (j = iFltPtchNum; j--;    )
    		    {   iVectPos = i*iFltPtchNum + j;	            iTDlg_Pdone[iVectPos] = 0;                           iTDlg_Sdone[iVectPos] = 0;	    
    		}   }
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
            /* EARTHQUAKE ITERATION LOOP STARTS */
    		while (iMD_EQongoing[0] == 1)  /*main variable to tell if I can leave the loop*/
    		{   	
				iMD_TotalRuptT[0]++;            iMD_EQongoing[0] = 0;           /*continuously count time since initiation, set "ongoing" to FALSE => only if more slip on patches is added, it gets to be reset*/
    			
    			for (i = iFltPtchNum;  i--;  )          {	    fTDg_DeltStrssH[i] = 0.0;	        fTDg_DeltStrssV[i] = 0.0;	        fTDg_DeltStrssN[i] = 0.0;       }
    			
    			for (i = iOFFSET_F[iRANK]; i--;    )
                {	/* determine new friction coefficient */
                    if (iTDl_WasActivated[i] == 1)
    				{	
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
    					    fTemp            = sqrtf(fTDl_STF_H[i*iTDg_GlobTTmax[0]+iTDl_STFcnt[i]]*fTDl_STF_H[i*iTDg_GlobTTmax[0]+iTDl_STFcnt[i]] + fTDl_STF_V[i*iTDg_GlobTTmax[0]+iTDl_STFcnt[i]]*fTDl_STF_V[i*iTDg_GlobTTmax[0]+iTDl_STFcnt[i]]);
							fTDl_CurrFric[i] = fTDl_DynFric[i] + (fTDl_StatFric[i]- fTDl_DynFric[i])* (1.0 - fTemp/fTDl_Curr_DcVal[i]); 
    						if      ((iTDl_StabType[i] != 3) && (fTDl_CurrFric[i] < fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
    						else if ((iTDl_StabType[i] == 3) && (fTDl_CurrFric[i] > fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }	
    				}   }
    				/*-------------------------------------------------------*/
    				/* determine amount of excess stress (if any available) and the corresponding slip amount */
					fTemp  = sqrtf(fTDl_CurrStrss_H[i]*fTDl_CurrStrss_H[i] + fTDl_CurrStrss_V[i]*fTDl_CurrStrss_V[i]);  /* this is the applied shear stress*/
					fTemp3 = (fTDl_CurrStrss_N[i] > 0.0) ? fTDl_CurrStrss_N[i] : 0.0;
					fTemp1 = fTemp - fTDl_CurrFric[i]*fTemp3; /* this is the amount of stress above current friction level in combination with normal stress at the location*/					
					fTemp2 = (fTDl_DynFric[i] *fTemp3 - fTemp)/(fK_FF_SS[iTDl_SelfLoc_F[i]] > fK_FF_DD[iTDl_SelfLoc_F[i]] ? fK_FF_SS[iTDl_SelfLoc_F[i]] : fK_FF_DD[iTDl_SelfLoc_F[i]]); /* this is releasable stress divided by self-stiffness => gives slip amount that would happen if patch fails "alone".... */					



					if ((iTDl_WasActivated[i] == 1)  || ((iTDl_WasActivated[i] == 0)  && (fTemp1 > fMD_CutStrss) && (fTemp2 >= fTDl_Curr_DcVal[i])))
					{	
						iwritePos = iTDl_STFcnt[i]%iMaxSTFLength; /*gives me remainder of wPos/MaxSTFlength => allows me to loop/overwrite the STF vectors */
						fTemp1    = (fTemp1 >fMD_CutStrss) ? fTemp1 : 0.0;
						if (iTDl_WasActivated[i] == 0)		    {	iTDl_WasActivated[i] = 1;	        iTDl_t0[i]         = iMD_TotalRuptT[0];		        }	
						if (fTemp1 > 0.0)                       {   iMD_EQongoing[0]     = 1;           iEndOfEventCntr[0] = 0;                             }
						
						fTDlg_STF_H[i*iMaxSTFLength+iwritePos] = -1.0*(fTemp1/fTemp *fTDl_CurrStrss_H[i]) /fK_FF_SS[iTDl_SelfLoc_F[i]]; /* horizontal slip at patch in current iteration (turn "n") */
						fTDlg_STF_V[i*iMaxSTFLength+iwritePos] = -1.0*(fTemp1/fTemp *fTDl_CurrStrss_V[i]) /fK_FF_DD[iTDl_SelfLoc_F[i]]; /* the "-1" is here because the stiffness matrix stuff gives the amount of slip nessessary to MAKE the observed stress; but I want to RELEASE it => opposite direction */
						fTDl_EventSlipH[i]                    += fTDlg_STF_H[i*iMaxSTFLength+iwritePos];				
						fTDl_EventSlipV[i]                    += fTDlg_STF_V[i*iMaxSTFLength+iwritePos];	
					}
					/*---------------------------------------------------------*/
					/* determine modified stress due to slip on fault patches */
					if (iTDl_WasActivated[i] == 1) 
    				{   /*---------------------------------------------------------*/ 						
						for (j = iFltPtchNum; j--;   )
					    {
						    iVectPos = i*iFltPtchNum + j;	
						    /*---------------------------------------------------------*/ 		
						    for (k = iTDlg_Pdone[iVectPos]; k <= (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]); k++) /*2nd part is negative if signal has not arrived, then the loop is skipped*/
						    {   
						        ireadPos = k%iMaxSTFLength;
						        fStkSlip = fTDlg_STF_H[i*iMaxSTFLength +ireadPos];
						        fDipSlip = fTDlg_STF_V[i*iMaxSTFLength +ireadPos];   
                                fTemp    = sqrt(fStkSlip*fStkSlip + fDipSlip*fDipSlip);						        
    				            if (fTemp > 0.0)
    				            {
   		                            if ((i + iSTARTPOS_F[iRANK]) == j)    /*this should be at "self" location -> then I cannot use the src-rcv-vector b/c it is zero...*/ 		
   		                            {   fTDg_DeltStrssH[j] += fStkSlip*fK_FF_SS[iVectPos] +fDipSlip*fK_FF_DS[iVectPos];
								        fTDg_DeltStrssV[j] += fStkSlip*fK_FF_SD[iVectPos] +fDipSlip*fK_FF_DD[iVectPos];  		                        						
   		                            } 
						            else
						            {   fTemp               = fStkSlip*fTDlg_LocSrcRcv_H[iVectPos] + fDipSlip*fTDlg_LocSrcRcv_V[iVectPos]; /*b/c SrcRcv vector is normalized, this length is already the length of the component in "P-direction" i.e., in SrcRcv direciton*/		        
						                fUsedStkSlip        = fTDlg_LocSrcRcv_H[iVectPos]*fTemp; 
						                fUsedDipSlip        = fTDlg_LocSrcRcv_V[iVectPos]*fTemp; 
						                fUsedNrmSlip        = fTDlg_LocSrcRcv_N[iVectPos]*fTemp; 
						        
						                fTDg_DeltStrssH[j] += fUsedStkSlip*fK_FF_SS[iVectPos] +fUsedDipSlip*fK_FF_DS[iVectPos] +fUsedNrmSlip*fK_FF_OS[iVectPos];
								        fTDg_DeltStrssV[j] += fUsedStkSlip*fK_FF_SD[iVectPos] +fUsedDipSlip*fK_FF_DD[iVectPos] +fUsedNrmSlip*fK_FF_OD[iVectPos];
								        fTDg_DeltStrssN[j] += fUsedStkSlip*fK_FF_SO[iVectPos] +fUsedDipSlip*fK_FF_DO[iVectPos] +fUsedNrmSlip*fK_FF_OO[iVectPos];
						    }   }   }  
						    iTDlg_Pdone[iVectPos] = (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]) >= 0 ? (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]+1) : 0; 
							/*---------------------------------------------------------*/			    
						    for (k = iTDlg_Sdone[iVectPos]; k <= (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]); k++) /*2nd part is negative if signal has not arrived, then the loop is skipped*/
						    {   
						        ireadPos = k%iMaxSTFLength;
						        fStkSlip = fTDlg_STF_H[i*iMaxSTFLength +ireadPos];
						        fDipSlip = fTDlg_STF_V[i*iMaxSTFLength +ireadPos];   
                                fTemp    = sqrt(fStkSlip*fStkSlip + fDipSlip*fDipSlip);						        
    				            if (fTemp > 0.0)
    				            {
   		                            if ((i + iSTARTPOS_F[iRANK]) == j)    /*this should be at "self" location -> then I cannot use the src-rcv-vector b/c it is zero...*/ 		
   		                            {   fTDg_DeltStrssH[j] += fStkSlip*fK_FF_SS[iVectPos] +fDipSlip*fK_FF_DS[iVectPos];
								        fTDg_DeltStrssV[j] += fStkSlip*fK_FF_SD[iVectPos] +fDipSlip*fK_FF_DD[iVectPos];  		                        						
   		                            } 
						            else
						            {   fTemp               = fStkSlip*fTDlg_LocSrcRcv_H[iVectPos] + fDipSlip*fTDlg_LocSrcRcv_V[iVectPos]; /*b/c SrcRcv vector is normalized, this length is already the length of the component in "P-direction" i.e., in SrcRcv direciton*/		        						                
						                fUsedStkSlip        = fStkSlip - fTDlg_LocSrcRcv_H[iVectPos]*fTemp; 
						                fUsedDipSlip        = fDipSlip - fTDlg_LocSrcRcv_V[iVectPos]*fTemp; 
						                fUsedNrmSlip        =    0.0   - fTDlg_LocSrcRcv_N[iVectPos]*fTemp; 
						        
						                fTDg_DeltStrssH[j] += fUsedStkSlip*fK_FF_SS[iVectPos] +fUsedDipSlip*fK_FF_DS[iVectPos] +fUsedNrmSlip*fK_FF_OS[iVectPos];
								        fTDg_DeltStrssV[j] += fUsedStkSlip*fK_FF_SD[iVectPos] +fUsedDipSlip*fK_FF_DD[iVectPos] +fUsedNrmSlip*fK_FF_OD[iVectPos];
								        fTDg_DeltStrssN[j] += fUsedStkSlip*fK_FF_SO[iVectPos] +fUsedDipSlip*fK_FF_DO[iVectPos] +fUsedNrmSlip*fK_FF_OO[iVectPos];
						    }   }   }  
						    iTDlg_Sdone[iVectPos] = (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]) >= 0 ? (iTDl_STFcnt[i] - iTDlg_TravTimes[iVectPos]+1) : 0; 
						    /*---------------------------------------------------------*/ 		
						}
						/*---------------------------------------------------------*/   	
						iTDl_STFcnt[i]++;
    			}	}
    			/*---------------------------------------------------------------*/
    			/* last was done "locally" now transfer information to all patches*/
    			MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssH, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssV, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(MPI_IN_PLACE, fTDg_DeltStrssN, iFltPtchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(MPI_IN_PLACE, iMD_EQongoing,      1        , MPI_UNSIGNED,MPI_MAX, MPI_COMM_WORLD);
    			/* combine "stress before EQ" and "EQ generated" changes to get "Currently applied stress" */	
    			for (i = iOFFSET[iRANK]; i--;    )
                {	
					fTDl_CurStrss_H[i] += fTDg_DeltStrssH[i+iSTARTPOS[iRANK]];
        			fTDl_CurStrss_V[i] += fTDg_DeltStrssV[i+iSTARTPOS[iRANK]];
        			fTDl_CurStrss_N[i] += fTDg_DeltStrssN[i+iSTARTPOS[iRANK]];
                }
        		/*---------------------------------------------------------------*/
        		if (iMD_EQongoing[0] == 0) /* make sure that all the signal is out of the system => even if no slip occurred in last step on any patch; there may still be stress in the system that has not reached a receiver => wait until they all got their share */
        		{   
        		    iMD_EQongoing[0]   = 1;                         iEndOfEventCntr[0] = iEndOfEventCntr[0] +1;
        		    if (iEndOfEventCntr[0] >= iMaxSTFLength)    {   iMD_EQongoing[0]   = 0;                             }   
        		}	
        		MPI_Allreduce(MPI_IN_PLACE, iMD_EQongoing,      1        , MPI_UNSIGNED,MPI_MAX, MPI_COMM_WORLD);				    
				}
    		} 
    		/* EARTHQUAKE ITERATION LOOP EXITED */
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
    		
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */   		
    		iEQl_ActPtchNum[0] = 0;
    		for (i = iOFFSET[iRANK]; i--;   )   /* now I define all the metrics of the event (its size, its MRF, etc...) */
            {   if (iTDl_WasActivated[i] == 1) 	
    			{	
    			    for (k = iTDl_PosOfLastSlp[i]; k--;    )       
    			    {	fEQl_MRFvals[iTDl_t0[i]+k] += sqrtf( fTDl_STF_H[i*iTDg_GlobTTmax[0]+k]*fTDl_STF_H[i*iTDg_GlobTTmax[0]+k] + fTDl_STF_V[i*iTDg_GlobTTmax[0]+k]*fTDl_STF_V[i*iTDg_GlobTTmax[0]+k]) *fTDl_Area[i];			}
    				
    				iTemp1            = iTDl_t0[i] + iTDl_PosOfLastSlp[i] +1; /* the plus 1 is here because it is actual length and not position in vector */
    				iEQl_MRFlength[0] = (iEQl_MRFlength[0] > iTemp1) ? iEQl_MRFlength[0] : iTemp1;
    			
    				fEQl_SeisPot[0]                     += sqrtf(fTDl_EventSlipH[i]*fTDl_EventSlipH[i] + fTDl_EventSlipV[i]*fTDl_EventSlipV[i] ) *fTDl_Area[i];	
					iEQl_ActPtchID[iEQl_ActPtchNum[0]]   = i+iSTARTPOS[iRANK];
		    		iEQl_t0ofPtch[iEQl_ActPtchNum[0]]    = iTDl_t0[i];
		    		iEQl_StabType[iEQl_ActPtchNum[0]]    = iTDl_StabType[i];		
					fEQl_DeltTofPtch[iEQl_ActPtchNum[0]] = sqrtf((fTDl_LoclStrssH[i]-fTDl_StrssB4_H[i])*(fTDl_LoclStrssH[i]-fTDl_StrssB4_H[i]) + (fTDl_LoclStrssV[i]-fTDl_StrssB4_V[i])*(fTDl_LoclStrssV[i]-fTDl_StrssB4_V[i]));
					fEQl_SlipHofPtch[iEQl_ActPtchNum[0]] = fTDl_EventSlipH[i];	
					fEQl_SlipVofPtch[iEQl_ActPtchNum[0]] = fTDl_EventSlipV[i];			
					iEQl_ActPtchNum[0]++;		
    		}   }	
    		/*------------------------------------------------------------------ */
    		MPI_Allreduce(iEQl_ActPtchNum, iEQg_ActPtchNum,       1          ,    MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    		MPI_Allreduce(iEQl_MRFlength,  iEQg_MRFlength,        1          ,    MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    		MPI_Allreduce(fEQl_SeisPot,    fEQg_SeisPot,          1          ,    MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
    		MPI_Allreduce(fEQl_MRFvals,    fEQg_MRFvals,     iMaxMRFlength   ,    MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
    		MPI_Allgather(iEQl_ActPtchNum,  1,  MPI_UNSIGNED, iEQg_WrtStartPos, 1, MPI_UNSIGNED,  MPI_COMM_WORLD);
	
    		for (i = 1; i < iSIZE; i++)		{	iEQg_WrtStartPos[i] += iEQg_WrtStartPos[i-1]; 			}	
    		for (i = iSIZE; i > 0; i--)		{	iEQg_WrtStartPos[i]  = iEQg_WrtStartPos[i-1]; 			}	
    		iEQg_WrtStartPos[0]  = 0;	
    		/* outoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutoutout */
    		/*------------------------------------------------------------------ */
    		MPI_File_open(MPI_COMM_WORLD, cFile1_Out,MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
    		MPI_File_get_size(fp_MPIout1, &offset1);
    		MPI_Barrier( MPI_COMM_WORLD );
    			
			if (iRANK == 0)
			{	
			    fTemp2 = (log10f(fEQg_SeisPot[0]*fMD_ShearMod[0])-9.1)/1.5; 
				MPI_File_write_at(fp_MPIout1,   offset1,                                       &fTimeYears,           1,              MPI_FLOAT, &status); 
				MPI_File_write_at(fp_MPIout1,   offset1 +  sizeof(float),                      &fTemp2,               1,              MPI_FLOAT, &status);
				MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float),                      iEQg_ActPtchNum,       1,              MPI_INT,   &status);
				MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float) +1*sizeof(int),       iEQg_MRFlength,        1,              MPI_INT,   &status);
				MPI_File_write_at(fp_MPIout1,   offset1 +2*sizeof(float) +2*sizeof(int),       fEQg_MRFvals,      iEQg_MRFlength[0],  MPI_FLOAT, &status);
			
				fprintf(stdout,"Earthquake time %f    Magn: %f act patches: %d   MFR length: %d  \n",fTimeYears,  fTemp2, iEQg_ActPtchNum[0],iEQg_MRFlength[0]);				
			}
			offset2 = offset1 +2*sizeof(int) +(2 +iEQg_MRFlength[0])*sizeof(float);
			MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(int),  iEQl_ActPtchID,    iEQl_ActPtchNum[0], MPI_INT,   &status);
			offset2 = offset2 +iEQg_ActPtchNum[0]*sizeof(int);
			MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(int),  iEQl_t0ofPtch,     iEQl_ActPtchNum[0], MPI_INT,   &status);
			offset2 = offset2 +iEQg_ActPtchNum[0]*sizeof(int);
			MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(float), fEQl_DeltTofPtch,  iEQl_ActPtchNum[0],MPI_FLOAT, &status);
			offset2 = offset2 +iEQg_ActPtchNum[0]*sizeof(float);
			MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(float), fEQl_SlipHofPtch, iEQl_ActPtchNum[0], MPI_FLOAT, &status);
			offset2 = offset2 +iEQg_ActPtchNum[0]*sizeof(float);
			MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(float), fEQl_SlipVofPtch, iEQl_ActPtchNum[0], MPI_FLOAT, &status);
		    offset2 = offset2 +iEQg_ActPtchNum[0]*sizeof(float);
		    MPI_File_write_at(fp_MPIout1,   offset2 +iEQg_WrtStartPos[iRANK]*sizeof(int),   iEQl_StabType,    iEQl_ActPtchNum[0], MPI_INT,   &status);
		
		
			MPI_File_close(&fp_MPIout1);
			/*---------------------------------------- */
			if ((fTemp2 > fKindModMinMag) &&(iWrteKinMod == 1))
			{ /* here will come the output of srf-file*/
			
			}
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */	
            iMD_ChangedStab[0] = 0; /* this will record if stability type was changed*/
            for (i = iOFFSET[iRANK]; i--;   ) /* RESET ALL VALUES FOR PATCHES THAT MOVED */
            {
                /*------------------------------------------------------------- */
                if (iTDl_WasActivated[i] == 1) 	
    			{   if (iMD_ChgFricBtwEQs[0] == 0)
		 		    {	fTDl_CurrFric[i]     = fTDl_StatFric[i];
		 		    }
		 		    else   
		 		    {	
		 		        fTDl_StatFric[i]   = fSD_RefStatFric[iTDl_SegID[i]]  *(1.0 +fRandVector[iRandPos[0]] *fSD_RefStatFr_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }       
                        fTDl_CurrFric[i]   = fTDl_StatFric[i];
                        fTDl_DynFric[i]    = fSD_RefDynFric[iTDl_SegID[i]]   *(1.0 +fRandVector[iRandPos[0]] *fSD_RefDynFr_vari[ iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }  
		 			    fTDl_Curr_DcVal[i] = fSD_CritSlipDist[iTDl_SegID[i]] *(1.0 +fRandVector[iRandPos[0]] *fSD_CritSlipD_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }      
		 		        /* test now again for stability type, while doing that check if any type was changed, if not then I don't have to recompute, if change is from cond stable to stable then I only change label, only if it goes from unstable to another state or vice versa, is it that I need to update loading*/

		 		        fTemp1 = (fTDl_DynFric[i] - fTDl_StatFric[i]) *fTDl_RefNormStrss[i]; /* this is the releasable stress amount*/
                        fTemp2 = fTemp1/(fK_SS[iTDl_SelfLoc_F[i]] > fK_DD[iTDl_SelfLoc_F[i]] ? fK_SS[iTDl_SelfLoc_F[i]] : fK_DD[iTDl_SelfLoc_F[i]]);//(0.5*(fK_SS[iTDl_SelfLoc_F[i]] + fK_DD[iTDl_SelfLoc_F[i]])); /* this is corresponding slip, sign is correct -> assuming that K-matrix at self-induced location is ALWAYS negative*/
        
                        if (fTemp2 > fTDl_Curr_DcVal[i]) /* have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc*/
                        {	if (iTDl_StabType[i] != 1)              {   iMD_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
                            iTDl_StabType[i] = 1; /* patch is unstable*/
                        }
                        else /*if (fTemp2 <= fTDl_Curr_DcVal[i]) -if stress drop i.e., corresponding slip is smaller than Dc - */
                        {	if (fTemp1 < 0.0)  /*  is still weakening but just with slip that is lower than Dc => cond. stable*/
        	                {	if (iTDl_StabType[i] == 1)          {   iMD_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
        	                    iTDl_StabType[i] = 2; /* patch is cond. stable*/
        	                }
        	                else /* no weakening but strengthening  => stable */
        	                {	if (iTDl_StabType[i] == 1)          {   iMD_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
        	                    iTDl_StabType[i] = 3; /* patch is stable*/
                    }   }   }   
                }
                /*------------------------------------------------------------- */
		 		for (k = iTDg_GlobTTmax[0]; k--;   )     {	fTDl_STF_H[i*iTDg_GlobTTmax[0]+k] = 0.0;	    fTDl_STF_V[i*iTDg_GlobTTmax[0]+k] = 0.0;		} 	
		 		
		 		iEQl_ActPtchID[i]    = 0;               iEQl_t0ofPtch[i]     = 0;		        iTDl_WasActivated[i] = 0;				
		 		iTDl_t0[i]           = 0;               iTDl_STFcnt[i]         = 0;               iTDl_PosOfLastSlp[i] = 0;
		 		fEQl_DeltTofPtch[i]  = 0.0;             fEQl_SlipHofPtch[i]  = 0.0;			    fEQl_SlipVofPtch[i]  = 0.0;
    			fTDl_EventSlipH[i]   = 0.0;			    fTDl_EventSlipV[i]   = 0.0;
    		}
    		iEQl_ActPtchNum[0] 	= 0;		            iEQg_ActPtchNum[0] 	 = 0;               iEQl_MRFlength[0]    = 0;						
    		iEQg_MRFlength[0]   = 0;                    
    		iMD_EQongoing[0]   = 0;                    iMD_EQongoing[0]    = 0;               
    		  fEQl_SeisPot[0]      = 0.0;		        fEQg_SeisPot[0]      = 0.0;
    		
    		for (i = iSIZE; i--;    )				{	iEQg_WrtStartPos[i] = 0;											    }
    		for (i = iMaxMRFlength; i--;    )	    {	fEQl_MRFvals[i]     = 0.0; 			    fEQg_MRFvals[i]      = 0.0;		}
            /*------------------------------------------------------------------------- */
            MPI_Allreduce(iMD_ChangedStab, iMD_ChangedStab,       1          ,    MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    		      	
		 	if (iMD_ChangedStab[0] == 1) /* next, update the current stressing rate*/
		 	{   /*-------------------------------------------------------------------*/
	            for (j = iPatchNum; j--;    )
                {   fTempVal1_F[j]   = 0.0;            fTempVal2_F[j]   = 0.0;
	            }
	            for (i = iOFFSET[iRANK]; i--;   )
                {		
    	            if ((iTDl_StabType[i] != 1) && ( iTDl_StabTypeChgd[i] == 1))/* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
    	            {	fTempSlipStk =  fTDl_RefStrssRateStk[i]/fK_SS[iTDl_SelfLoc_F[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
    		            fTempSlipDip =  fTDl_RefStrssRateDip[i]/fK_DD[iTDl_SelfLoc_F[i]];
    
    		            for (j = iPatchNum; j--;    )
    		            {	if (iTDl_StabType[i] == 1) /* apply the creep slip to load the  unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */	
    			            {   iVectPos                = i*iPatchNum + j;  				
    				            fTempVal1_F[j] += fTempSlipStk*fK_SS[iVectPos] + fTempSlipDip*fK_DS[iVectPos];
    				            fTempVal2_F[j] += fTempSlipStk*fK_SD[iVectPos] + fTempSlipDip*fK_DD[iVectPos];	    
	            }   }   }   }
                /*-------------------------------------------------------------------*/
                MPI_Allreduce(fTempVal1_F, fTDg_StrssRteChgStk, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
                MPI_Allreduce(fTempVal2_F, fTDg_StrssRteChgDip, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);   
                for (i = iPatchNum; i--;    )		{	fTempVal1_F[i] = 0.0;        fTempVal2_F[i] = 0.0;                 } 
                /*-------------------------------------------------------------------*/
	            for (i = iOFFSET[iRANK]; i--;    )
                {   if (iTDl_StabType[i] == 1) /* if it is a stick-slip patch/cell... =>apply the creep slip to load the  locked/strong patches	*/
    	            {	fTDl_CurStrssRateStk[i]  = fTDl_RefStrssRateStk[i] +fTDg_StrssRteChgStk[i+iSTARTPOS[iRANK]];			
    		            fTDl_CurStrssRateDip[i]  = fTDl_RefStrssRateDip[i] +fTDg_StrssRteChgDip[i+iSTARTPOS[iRANK]];	      
	                }   
	                else
	                {   fTDl_CurStrssRateStk[i]  = 0.0;
	                    fTDl_CurStrssRateDip[i]  = 0.0;             
		    }   }   }
		 	/*--------------------------------------------------- */	    
		 	for (i = iOFFSET[iRANK]; i--;   )
            {   if (iTDl_StabType[i] == 1) /* in that case (don't use post seismic and patch is stick-slip i.e., seismogenic, then I'll just consider past EQ stress to be new stress on fault */
                {   fTDl_StrssB4_H[i] = fTDl_LoclStrssH[i];                 fTDl_StrssB4_V[i] = fTDl_LoclStrssV[i];             fTDl_StrssB4_N[i] = fTDl_RefNormStrss[i];
                }
                else
                {   
                    fTemp1 = sqrtf(fTDl_LoclStrssH[i]*fTDl_LoclStrssH[i] + fTDl_LoclStrssV[i]*fTDl_LoclStrssV[i]);  /* compute the applied shear stress after EQ on fault patch*/
					fTemp2 = fTDl_CurrFric[i]*fTDl_RefNormStrss[i];  /* compute the patch'es strength*/  
                    if (iUsePostSeismic == 1)
		            {   /* here i would have to put the excess stress (same as below) and put into postseismic relaxation stressing */
		                fTDl_StrssB4_H[i]     = fTDl_LoclStrssH[i];         fTDl_StrssB4_V[i]    = fTDl_LoclStrssV[i];            fTDl_StrssB4_N[i]     = fTDl_RefNormStrss[i];
		                
		                fTemp = fTDl_PostStrssat_t0[i]*expf(-(fTimeStpInSecs*fTDl_PostStrssTStep[i])/fTDl_CharRelaxTime[i]);/* is going to be the post-seis "strength" at current time -relative to PostStrssTime*/
		                
		                if ((fTemp1 - fTemp2) > fTemp)
		                {   fTDl_PostStrssat_t0[i] = fTemp1 - fTemp2;  
		                    fTDl_PostStrssTStep[i]  = 1.0;              
        	        }   }   
                    else
                    {   if (fTemp1 > fTemp2) /* if that strength is below the currently applied stress, then lower the stress to be equal to static strength -would be ok to lower by some additional random fraction??*/
					    {   fTDl_StrssB4_H[i] = fTemp2/fTemp1 *fTDl_LoclStrssH[i];    fTDl_StrssB4_V[i] = fTemp2/fTemp1 *fTDl_LoclStrssV[i];    fTDl_StrssB4_N[i] = fTDl_RefNormStrss[i];
                        }
                        else
                        {   fTDl_StrssB4_H[i] = fTDl_LoclStrssH[i];                     fTDl_StrssB4_V[i] = fTDl_LoclStrssV[i];                 fTDl_StrssB4_N[i] = fTDl_RefNormStrss[i];
            }   }   }   }
            /*------------------------------------------------------------------ */
   	 	}
   	 	MPI_Barrier( MPI_COMM_WORLD );
   	}
	MPI_File_open(MPI_COMM_WORLD, cFileN_Out1, MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
	if (iRANK == 0)			{		MPI_File_write_at(fp_MPIout1, 0, &iEQcounter, 1, MPI_UNSIGNED,  &status);			}
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
	return ans;
}
