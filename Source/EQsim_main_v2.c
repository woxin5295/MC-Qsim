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
extern void LoadInputParameter(char **argv, unsigned int *iMDg_ChgFricBtwEQs, float *fMDg_AddNrmStrss,unsigned int *iSDg_FricLawUSED, unsigned int *iTDg_V1temp, unsigned int *iTDg_V2temp, unsigned int *iTDg_V3temp, unsigned int *iTDg_StabTypetemp, unsigned int *iTDg_SegIDtemp, float *fTDg_StrsRteStktemp, float *fTDg_StrsRteDiptemp, float *fMDg_Vp, float *fMDg_Poisson, float *fMDg_Lambda, float *fMDg_ShearMod, float *fMDg_MedDense,  float *fSDg_RefStatFric, float *fSDg_RefStatFr_vari, float *fSDg_RefDynFric, float *fSDg_RefDynFr_vari, float *fSDg_CritSlipDist, float *fSDg_CritSlipD_vari, float *fSDg_CritSlipVelo, float *fTDg_CentEpos, float *fTDg_CentNpos, float *fTDg_CentZpos, float *fTDg_StatFrictemp, float *fTDg_DynFrictemp, float *fVDg_Epostemp, float *fVDg_Npostemp, float *fVDg_Zpostemp, unsigned int iSegmentNum, unsigned int iPatchNum, unsigned int iVertexNum, unsigned int iUsdGrid, unsigned int iSrcGrid);								
extern void    DefineMoreParas(unsigned int *iTravTimes, unsigned int *iTDl_ttmax,  float *fTDl_Area, float *fTDl_RefNormStrss, float *fTDl_Curr_DcVal, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fVDg_Epostemp, const float *fVDg_Npostemp, const float *fVDg_Zpostemp, const unsigned int *iTDg_V1temp, const unsigned int *iTDg_V2temp, const unsigned int *iTDg_V3temp, const unsigned int *iTDl_SegID, const float *fSDg_CritSlipDist, const float *fSDg_CritSlipD_vari, const int *iSTARTPOS, const int *iOFFSET, int iRANK, float fMDg_Vp, float fMDg_MedDense, float fMDg_AddNrmStrss, float fIntSeisTimeStp,  float fMDg_g, float fdeltTincr, unsigned int iPatchNum, float *fRandVector, unsigned int *iRandPos, unsigned int iRandNumber);
extern void     Build_K_Matrix(float *fK_SS, float *iK_SD, float *iK_SO, float *iK_DS, float *iK_DD, float *iK_DO, unsigned int *iTDg_SelfIndLoctemp, const float *fVDg_Epostemp,  const float *fVDg_Npostemp, const float *fVDg_Zpostemp,  const unsigned int *iTDg_V1temp, const unsigned int *iTDg_V2temp, const unsigned int *iTDg_V3temp, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, float fMDg_UnitSlip, float fMDg_ShearMod, float fMDg_Lambda, unsigned int *iTDl_StabType, const float *fTDl_Curr_DcVal, const float *fTDl_StatFric, const float *fTDl_DynFric, const float *fTDl_RefNormStrss, const int *iSTARTPOS, const int *iOFFSET, int iRANK, unsigned int iPatchNum, const float *fTDl_Area);
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int main(int argc, char **argv)
{   if (argc != 10)	{	fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -n 4 ./Run_EQsim InputName RealizationNumber GridNumber(Faults) GridNumber(Loading) StoreRuptureModel4LargeEQs LargeEQmagnitude WhichDeltT2Use UsePostSeis EQrecordLength\n");		exit(10);			}
	/*-------------------------------------------------------------------*/
    unsigned int	i,				j,              k;
    int	    		iRANK,         	iSIZE;
    /*-------------------------------------------------------------------*/
    MPI_Init(&argc, &argv);               /* start MPI processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &iRANK);/* current proc. id    */
    MPI_Comm_size(MPI_COMM_WORLD, &iSIZE);/* # of processes      */ 

    MPI_Status status;
    MPI_Offset offset1;
    MPI_Offset offset2;
    /*-------------------------------------------------------------------*/
    int 	iOFFSET[iSIZE];		for (i = 0; i < iSIZE; i++)             		{       iOFFSET[i]      = 0;                              }
    int 	iSTARTPOS[iSIZE];   for (i = 0; i < iSIZE; i++)         			{       iSTARTPOS[i]    = 0;                              }
    /*-------------------------------------------------------------------*/
  	unsigned int 	iUsdGrid,			iSrcGrid,			iWrteKinMod,			iDeltT2Use,             iUsePostSeismic;
    float        	fKindModMinMag,     fRecLgth;
    char    *retch;
    /*----------------------------*/
  	sscanf(argv[3],"%d",  &iUsdGrid);    						sscanf(argv[4],"%d",  &iSrcGrid);    
    sscanf(argv[5],"%d",  &iWrteKinMod);						sscanf(argv[6],"%f",  &fKindModMinMag);  
    sscanf(argv[7],"%d",  &iDeltT2Use);	                        sscanf(argv[8],"%u",  &iUsePostSeismic);						
    sscanf(argv[9],"%f",  &fRecLgth); 		/*this is reclength in years*/	
	/*
      argv[1] is the file name -but only the "core" of the filename => like "AqabaFaults" for all models like "AqabaFaults_1_Roughn.dat" etc...
      argv[2] is the realization number (if I created 100 random models... -> which one to use..)
      argv[3] grid resolution to use (if more than one was used and stored in the binary files => which one to use
      argv[4] grid resolution of source (loading) grid => not used at the moment but needed for output
      argv[5] switch to decide whether to store large-earthquake kinematic models
      argv[6] this defines the magnitude for those "large-earthquakes" => EQs larger than this will be stored...
      argv[7] whether to use deltT based on v_p (== 1) or to use inf deltT (instantaneous signal distribution,  == 0)
      argv[8] this is length of record (in years)
	*/  
    /*-------------------------------------------------------------------*/
    unsigned int 	iSegmentNum,		iGridNum,			iPatchNum,			        iVertexNum; 
    float  			fMeanLegLgth,		fDummy;
  	char			ctempVals[512],		cFileN_In[512],		cFileN_Out1[512],			cFileN_Out2[512];   /* out1 is for the overall catalog; out2 is the stf for each patch that had slip => both are binaries!!!  	*/
	FILE    *fpIn;  
	/*----------------------------*/
    strcpy(cFileN_In,  argv[1]);              			   		strcat(cFileN_In,"_Summary.flt");
    strcpy(cFileN_Out1,argv[1]);              					strcat(cFileN_Out1,"_");				strcat(cFileN_Out1,argv[2]);		strcat(cFileN_Out1,"_Catalog.dat");
    strcpy(cFileN_Out2,argv[1]);              					strcat(cFileN_Out2,"_");				strcat(cFileN_Out2,argv[2]);		strcat(cFileN_Out2,"_STFs.dat");

    if ((fpIn = fopen(cFileN_In,"r")) == NULL)       		{   printf("Error -cant open *.flt file. in Main function -trying to open _Summary.flt here... file...\n");      exit(10);     }
	retch =fgets(ctempVals, 512, fpIn);							sscanf(ctempVals,"%*s %d",  &iSegmentNum); 
    retch =fgets(ctempVals, 512, fpIn);          				sscanf(ctempVals,"%*s %d",  &iGridNum);    
    /*----------------------------*/
    iUsdGrid = (iUsdGrid  < iGridNum) ? iUsdGrid : iGridNum;
    /*----------------------------*/
    for (i = 1; i < iUsdGrid; i++)			   				{   retch =fgets(ctempVals, 512, fpIn);   retch =fgets(ctempVals, 512, fpIn);}
    retch =fgets(ctempVals, 512, fpIn);							sscanf(ctempVals,"%*s %d  %d", &iPatchNum, &iVertexNum); 
    retch =fgets(ctempVals, 512, fpIn);          				sscanf(ctempVals,"%*s %e  %e", &fMeanLegLgth, &fDummy);
    fclose(fpIn);
    fMeanLegLgth = fMeanLegLgth*1000; /* now it is in meters */
    /*-------------------------------------------------------------------*/  
    unsigned int    iBASEelem;			iBASEelem   = (unsigned int)(iPatchNum/iSIZE); 
    unsigned int    iADDelem;			iADDelem    = (unsigned int)(iPatchNum%iSIZE);
    /*----------------------------*/
    for (i = 0; i < iSIZE;    i++)          				{ 	iOFFSET[i]     = iBASEelem;                       }
    for (i = 0; i < iADDelem; i++)       					{   iOFFSET[i]    += 1;                               }   
    for (i = 1; i < iSIZE;    i++)          				{   iSTARTPOS[i]   = iSTARTPOS[i-1] + iOFFSET[i-1];   }
    fprintf(stdout,"My Rank: %d  TriNum %d    BASEelem  %d   ADDelem  %d  STARTPOS  %d OFFSET   %d; start: %d     end: %d\n",iRANK, iPatchNum, iBASEelem, iADDelem, iSTARTPOS[iRANK], iOFFSET[iRANK], iSTARTPOS[iRANK],  (iSTARTPOS[iRANK]+iOFFSET[iRANK])  );
    /*-------------------------------------------------------------------*/ 
	unsigned int	*iMDg_ChgFricBtwEQs;        iMDg_ChgFricBtwEQs = (unsigned int *) calloc(1, sizeof(unsigned int));  
	float           *fMDg_AddNrmStrss;          fMDg_AddNrmStrss   = (float *)        calloc(1, sizeof(float));
    float           *fMDg_Vp;                   fMDg_Vp            = (float *)        calloc(1, sizeof(float));	                                
    float	        *fMDg_Poisson;              fMDg_Poisson       = (float *)        calloc(1, sizeof(float));						
    float	        *fMDg_Lambda;               fMDg_Lambda        = (float *)        calloc(1, sizeof(float));						            
    float	        *fMDg_ShearMod;             fMDg_ShearMod      = (float *)        calloc(1, sizeof(float));			
    float		    *fMDg_MedDense;             fMDg_MedDense      = (float *)        calloc(1, sizeof(float));		
	
    unsigned int 	*iSDg_FricLawUSED;          iSDg_FricLawUSED   = (unsigned int *) calloc(iSegmentNum, sizeof(unsigned int));	
    float           *fSDg_RefStatFric;          fSDg_RefStatFric   = (float *)        calloc(iSegmentNum, sizeof(float));                       
    float           *fSDg_RefStatFr_vari;       fSDg_RefStatFr_vari= (float *)        calloc(iSegmentNum, sizeof(float));		
    float 	        *fSDg_RefDynFric;           fSDg_RefDynFric    = (float *)        calloc(iSegmentNum, sizeof(float));                       
    float 	        *fSDg_RefDynFr_vari;        fSDg_RefDynFr_vari = (float *)        calloc(iSegmentNum, sizeof(float));		
	float 	        *fSDg_CritSlipDist;         fSDg_CritSlipDist  = (float *)        calloc(iSegmentNum, sizeof(float));                       
	float 	        *fSDg_CritSlipD_vari;       fSDg_CritSlipD_vari= (float *)        calloc(iSegmentNum, sizeof(float));		    				
	float 	        *fSDg_CritSlipVelo;         fSDg_CritSlipVelo  = (float *)        calloc(iSegmentNum, sizeof(float));
	      
    unsigned int	*iTDg_V1temp;               iTDg_V1temp	       = (unsigned int *) calloc(iPatchNum,   sizeof(unsigned int));         
    unsigned int	*iTDg_V2temp;               iTDg_V2temp        = (unsigned int *) calloc(iPatchNum,   sizeof(unsigned int));    	   
    unsigned int    *iTDg_V3temp;               iTDg_V3temp	 	   = (unsigned int *) calloc(iPatchNum,   sizeof(unsigned int));         
    unsigned int	*iTDg_StabTypetemp;         iTDg_StabTypetemp    = (unsigned int *) calloc(iPatchNum,   sizeof(unsigned int));    	   
	unsigned int	*iTDg_SegIDtemp;            iTDg_SegIDtemp     = (unsigned int *) calloc(iPatchNum,   sizeof(unsigned int));
	float           *fTDg_StrsRteStktemp;       fTDg_StrsRteStktemp= (float *)        calloc(iPatchNum,   sizeof(float));                  
	float           *fTDg_StrsRteDiptemp;       fTDg_StrsRteDiptemp= (float *)        calloc(iPatchNum,   sizeof(float));				
	float           *fTDg_CentEpos;             fTDg_CentEpos      = (float *)        calloc(iPatchNum,   sizeof(float));	                    
	float 	        *fTDg_CentNpos;             fTDg_CentNpos      = (float *)        calloc(iPatchNum,   sizeof(float));					
	float 	        *fTDg_CentZpos;             fTDg_CentZpos      = (float *)        calloc(iPatchNum,   sizeof(float));
	float 		    *fTDg_StatFrictemp;         fTDg_StatFrictemp  = (float *)        calloc(iPatchNum,   sizeof(float));                       
	float 	        *fTDg_DynFrictemp;          fTDg_DynFrictemp   = (float *)        calloc(iPatchNum,   sizeof(float));			
				
	float           *fVDg_Epostemp;             fVDg_Epostemp      = (float *)        calloc(iVertexNum,  sizeof(float));	                    
	float           *fVDg_Npostemp;             fVDg_Npostemp      = (float *)        calloc(iVertexNum,  sizeof(float));
	float 	        *fVDg_Zpostemp;             fVDg_Zpostemp      = (float *)        calloc(iVertexNum,  sizeof(float));
    /*-----------------------------------*/
	if (iRANK == 0) 
    {	
    
        LoadInputParameter(argv, iMDg_ChgFricBtwEQs, fMDg_AddNrmStrss, iSDg_FricLawUSED, iTDg_V1temp, iTDg_V2temp, iTDg_V3temp, iTDg_StabTypetemp, iTDg_SegIDtemp, fTDg_StrsRteStktemp, fTDg_StrsRteDiptemp, fMDg_Vp, fMDg_Poisson, fMDg_Lambda, fMDg_ShearMod, fMDg_MedDense,  fSDg_RefStatFric, fSDg_RefStatFr_vari, fSDg_RefDynFric, fSDg_RefDynFr_vari, fSDg_CritSlipDist, fSDg_CritSlipD_vari, fSDg_CritSlipVelo, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fTDg_StatFrictemp, fTDg_DynFrictemp, fVDg_Epostemp, fVDg_Npostemp, fVDg_Zpostemp, iSegmentNum, iPatchNum, iVertexNum, iUsdGrid, iSrcGrid);						
	
	}
    /*-----------------------------------*/
    float           fIntSeisTimeStp = 2.0/365.25; /* == 1 is one-day steps; is fraction of a year*/
    float           fTimeInSecs     = fIntSeisTimeStp*31536000.0; /* now the time is in seconds */
    unsigned int    *iTDl_StabType;             iTDl_StabType     = (unsigned int *) calloc(iOFFSET[iRANK],   sizeof(unsigned int));
    unsigned int    *iTDl_StabTypeChgd;         iTDl_StabTypeChgd = (unsigned int *) calloc(iOFFSET[iRANK],   sizeof(unsigned int)); 
    unsigned int    *iTDl_SegID;                iTDl_SegID        = (unsigned int *) calloc(iOFFSET[iRANK],   sizeof(unsigned int));
    float           *fTDl_RefStrssRateStk;      fTDl_RefStrssRateStk = (float *)     calloc(iOFFSET[iRANK],   sizeof(float));
    float           *fTDl_RefStrssRateDip;      fTDl_RefStrssRateDip = (float *)     calloc(iOFFSET[iRANK],   sizeof(float));
    float           *fTDl_CurStrssRateStk;      fTDl_CurStrssRateStk = (float *)     calloc(iOFFSET[iRANK],   sizeof(float));
    float           *fTDl_CurStrssRateDip;      fTDl_CurStrssRateDip = (float *)     calloc(iOFFSET[iRANK],   sizeof(float));  
    float           *fTDl_PostStrssat_t0;       fTDl_PostStrssat_t0  = (float *)     calloc(iOFFSET[iRANK],   sizeof(float));    
    float           *fTDl_PostStrssTime;        fTDl_PostStrssTime= (float *)        calloc(iOFFSET[iRANK],   sizeof(float));        
    float           *fTDl_StatFric;             fTDl_StatFric     = (float *)        calloc(iOFFSET[iRANK],   sizeof(float));
    float           *fTDl_CurrFric;             fTDl_CurrFric     = (float *)        calloc(iOFFSET[iRANK],   sizeof(float));
    float           *fTDl_DynFric;              fTDl_DynFric      = (float *)        calloc(iOFFSET[iRANK],   sizeof(float));
    /*----------------------------*/
    MPI_Bcast(iMDg_ChgFricBtwEQs,           1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(fMDg_AddNrmStrss,             1, MPI_FLOAT,    0, MPI_COMM_WORLD);	
	MPI_Bcast(iSDg_FricLawUSED,   iSegmentNum, MPI_UNSIGNED, 0, MPI_COMM_WORLD);			
	MPI_Bcast(iTDg_V1temp,          iPatchNum, MPI_UNSIGNED, 0, MPI_COMM_WORLD);			
	MPI_Bcast(iTDg_V2temp,          iPatchNum, MPI_UNSIGNED, 0, MPI_COMM_WORLD);			
    MPI_Bcast(iTDg_V3temp,          iPatchNum, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(fMDg_Vp,                      1, MPI_FLOAT, 0, MPI_COMM_WORLD);		            
	MPI_Bcast(fMDg_Poisson,                 1, MPI_FLOAT, 0, MPI_COMM_WORLD);	
    MPI_Bcast(fMDg_Lambda,                  1, MPI_FLOAT, 0, MPI_COMM_WORLD);                   
    MPI_Bcast(fMDg_ShearMod,                1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fMDg_MedDense,                1, MPI_FLOAT, 0, MPI_COMM_WORLD);                    		
	
    MPI_Bcast(fSDg_RefStatFric,   iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	                
    MPI_Bcast(fSDg_RefStatFr_vari,iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
    MPI_Bcast(fSDg_RefDynFric,    iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);                   
    MPI_Bcast(fSDg_RefDynFr_vari, iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fSDg_CritSlipDist,  iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);			        
    MPI_Bcast(fSDg_CritSlipD_vari,iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
    MPI_Bcast(fSDg_CritSlipVelo,  iSegmentNum, MPI_FLOAT, 0, MPI_COMM_WORLD);		
	   
	MPI_Bcast(fTDg_CentEpos,        iPatchNum, MPI_FLOAT, 0, MPI_COMM_WORLD);			        
	MPI_Bcast(fTDg_CentNpos,        iPatchNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(fTDg_CentZpos,        iPatchNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	
	
    MPI_Bcast(fVDg_Epostemp,       iVertexNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	                
	MPI_Bcast(fVDg_Npostemp,       iVertexNum, MPI_FLOAT, 0, MPI_COMM_WORLD);	        
	MPI_Bcast(fVDg_Zpostemp,       iVertexNum, MPI_FLOAT, 0, MPI_COMM_WORLD);			
	
    MPI_Scatterv(iTDg_SegIDtemp,      iOFFSET, iSTARTPOS, MPI_UNSIGNED, iTDl_SegID,        iOFFSET[iRANK], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Scatterv(iTDg_StabTypetemp,   iOFFSET, iSTARTPOS, MPI_UNSIGNED, iTDl_StabType,     iOFFSET[iRANK], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Scatterv(fTDg_StrsRteStktemp, iOFFSET, iSTARTPOS, MPI_FLOAT, fTDl_RefStrssRateStk, iOFFSET[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	MPI_Scatterv(fTDg_StrsRteDiptemp, iOFFSET, iSTARTPOS, MPI_FLOAT, fTDl_RefStrssRateDip, iOFFSET[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);

	MPI_Scatterv(fTDg_StatFrictemp,   iOFFSET, iSTARTPOS, MPI_FLOAT,    fTDl_StatFric,     iOFFSET[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD); 
	MPI_Scatterv(fTDg_StatFrictemp,   iOFFSET, iSTARTPOS, MPI_FLOAT,    fTDl_CurrFric,     iOFFSET[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD); 
	MPI_Scatterv(fTDg_DynFrictemp,    iOFFSET, iSTARTPOS, MPI_FLOAT,    fTDl_DynFric,      iOFFSET[iRANK], MPI_FLOAT,    0, MPI_COMM_WORLD);
	/*-------------------------------------------------------------------*/    
	float fdeltTincr; /* this one is for time steps during rupture => depends on distance of patches... but this needs calibration I believe*/
	  
	if (iDeltT2Use == 1)    {	fdeltTincr = 2.0*fMeanLegLgth/fMDg_Vp[0];           fdeltTincr = floor(fdeltTincr*10.0)/10.0;           }
	else                    {	fdeltTincr = 1.0E+09;                                                                                   }
	
	if (iRANK == 0)		{		fprintf(stdout,"delta T (coseismic) %f   fMeanLegLgth %f  RecLength %f\n",fdeltTincr,fMeanLegLgth,fRecLgth);			}
    /*-------------------------------------------------------------------*/
	unsigned int    iRandNumber = 1000000; 
	float           fMDg_g      = 9.81;
	float           fViscosity  = 1.0E+13; /* in Pa seconds; value doesn't seem right, but need something that corresponds to patch stiffness => to give proper char. time scale for post-seismic decay*/
	long            *iSeedVal;          iSeedVal          = (long *)         calloc(       1   ,              sizeof(long)); 
	unsigned int    *iTravTimes;        iTravTimes        = (unsigned int *) calloc(iOFFSET[iRANK]*iPatchNum, sizeof(unsigned int));
	unsigned int    *iTDl_ttmax;        iTDl_ttmax        = (unsigned int *) calloc(iOFFSET[iRANK],           sizeof(unsigned int));
	float           *fTDl_RefNormStrss; fTDl_RefNormStrss = (float *)        calloc(iOFFSET[iRANK],           sizeof(float));
    float           *fTDl_Area;         fTDl_Area         = (float *)        calloc(iOFFSET[iRANK],           sizeof(float));
    float           *fTDl_Curr_DcVal;   fTDl_Curr_DcVal   = (float *)        calloc(iOFFSET[iRANK],           sizeof(float));
	
	float           *fRandVector;       fRandVector       = (float *)        calloc(iRandNumber,              sizeof(float));
	unsigned int    *iRandPos;          iRandPos          = (unsigned int *) calloc(       1   ,              sizeof(unsigned int)); 
	
	time_t  t;
    srand((unsigned) time(&t)+iRANK*10);  /* initialize the random number generator; will generate 1E+6 random numbers => should be enough*/
	iSeedVal[0] = rand();
	for (i = 0; i < iRandNumber; i++ )       {           fRandVector[i] = ran0_inmain(iSeedVal);       } /* using a generator as described in "Numerical Recipes, page 279*/
	/*-------------------------------------------------------------------*/
	
		
    DefineMoreParas(iTravTimes, iTDl_ttmax, fTDl_Area, fTDl_RefNormStrss, fTDl_Curr_DcVal, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fVDg_Epostemp, fVDg_Npostemp, fVDg_Zpostemp, iTDg_V1temp, iTDg_V2temp, iTDg_V3temp, iTDl_SegID, fSDg_CritSlipDist, fSDg_CritSlipD_vari, iSTARTPOS, iOFFSET, iRANK, fMDg_Vp[0], fMDg_MedDense[0], fMDg_AddNrmStrss[0], fIntSeisTimeStp, fMDg_g, fdeltTincr, iPatchNum, fRandVector, iRandPos, iRandNumber);
	
	
	/*-------------------------------------------------------------------*/ 
	/* once I do include postseismic deformation or after slip, I'll have to add the additional patches here => iPatchNum +iBoundPatchNum */
	int          iVectPos;
	float        fMDg_UnitSlip = 1.0E-4*fMeanLegLgth; /* slip in meter => if meanleglength is 1000m, then a slip of .1m is used */
	float        *fK_SS;                fK_SS           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));	        
	float        *fK_SD;                fK_SD           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));					
    float        *fK_SO;                fK_SO           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));	        
    float        *fK_DD;                fK_DD           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));						
    float        *fK_DS;                fK_DS           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));	        
    float        *fK_DO;                fK_DO           = (float *)          calloc(iOFFSET[iRANK]*iPatchNum,sizeof(float));	
    unsigned int *iTempValues1;         iTempValues1    = (unsigned int *) calloc(iPatchNum,      sizeof(unsigned int)); 
    unsigned int *iTempValues2;         iTempValues2    = (unsigned int *) calloc(iPatchNum,      sizeof(unsigned int)); 
    unsigned int *iTDl_LocSelfInd;      iTDl_LocSelfInd = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int)); 		
    
    if (iRANK == 0)			{ 		fprintf(stdout,"Building K-matrix....\n");			}
	/*-------------------------------------------------------------------*/
	
	
    Build_K_Matrix(fK_SS, fK_SD, fK_SO, fK_DS, fK_DD, fK_DO, iTempValues1, fVDg_Epostemp, fVDg_Npostemp, fVDg_Zpostemp,  iTDg_V1temp, iTDg_V2temp, iTDg_V3temp, fTDg_CentEpos, fTDg_CentNpos, fTDg_CentZpos, fMDg_UnitSlip, fMDg_ShearMod[0], fMDg_Lambda[0], iTDl_StabType, fTDl_Curr_DcVal, fTDl_StatFric, fTDl_DynFric, fTDl_RefNormStrss, iSTARTPOS, iOFFSET, iRANK, iPatchNum,fTDl_Area);

	
	/*-------------------------------------------------------------------*/
	MPI_Reduce(iTempValues1, iTempValues2, iPatchNum, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Scatterv(iTempValues2, iOFFSET, iSTARTPOS, MPI_UNSIGNED, iTDl_LocSelfInd, iOFFSET[iRANK], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (iRANK == 0)			{ 		fprintf(stdout,"K-matrix done \n");			}
			
	/*-------------------------------------------------------------------*/  
	float *fTDg_StrssRteChgStk;     fTDg_StrssRteChgStk   = (float *) calloc(iPatchNum,      sizeof(float)); 
	float *fTDg_StrssRteChgDip;     fTDg_StrssRteChgDip   = (float *) calloc(iPatchNum,      sizeof(float)); 
	float *fTempValues1;            fTempValues1          = (float *) calloc(iPatchNum,      sizeof(float)); 
	float *fTempValues2;            fTempValues2          = (float *) calloc(iPatchNum,      sizeof(float)); 
	float *fTDl_CharRelaxTime;      fTDl_CharRelaxTime    = (float *) calloc(iOFFSET[iRANK], sizeof(float )); 		
	float fTempSlipStk,             fTempSlipDip;
	/*-------------------------------------------------------------------*/   
   	for (i = iOFFSET[iRANK]; i--;   )
    {		
    	fTDl_CharRelaxTime[i] = fabs(fViscosity/(fK_SS[iTDl_LocSelfInd[i]] > fK_DD[iTDl_LocSelfInd[i]] ? fK_SS[iTDl_LocSelfInd[i]] : fK_DD[iTDl_LocSelfInd[i]]));//(0.5*(fK_SS[iTDl_LocSelfInd[i]]+fK_DD[iTDl_LocSelfInd[i]]))); /* doesn't belong here but needs to be done at about that time and this loop is already here*/

    	if (iTDl_StabType[i] != 1) /* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
    	{	
    	    fTempSlipStk =  fTDl_RefStrssRateStk[i]/fK_SS[iTDl_LocSelfInd[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
    		fTempSlipDip =  fTDl_RefStrssRateDip[i]/fK_DD[iTDl_LocSelfInd[i]];
    
    		for (j = iPatchNum; j--;    )
    		{	if (iTDl_StabType[j] == 1) /* apply the creep slip to load the  unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */	
    			{   iVectPos         = i*iPatchNum + j;  				
    				fTempValues1[j] += fTempSlipStk*fK_SS[iVectPos] + fTempSlipDip*fK_DS[iVectPos];
    				fTempValues2[j] += fTempSlipStk*fK_SD[iVectPos] + fTempSlipDip*fK_DD[iVectPos];	
	}   }   }   }
    /*-------------------------------------------------------------------*/
    MPI_Allreduce(fTempValues1, fTDg_StrssRteChgStk, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
    MPI_Allreduce(fTempValues2, fTDg_StrssRteChgDip, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);    
    /*-------------------------------------------------------------------*/
	for (i = iOFFSET[iRANK]; i--;    )
    {	
        if (iTDl_StabType[i] == 1) /* if it is an unstable atch/cell... =>apply the creep slip to load */
    	{	fTDl_CurStrssRateStk[i]  = fTDl_RefStrssRateStk[i] +fTDg_StrssRteChgStk[i+iSTARTPOS[iRANK]];			
    		fTDl_CurStrssRateDip[i]  = fTDl_RefStrssRateDip[i] +fTDg_StrssRteChgDip[i+iSTARTPOS[iRANK]];	      
	    }   
	    else
	    {   fTDl_CurStrssRateStk[i]  = 0.0;
	        fTDl_CurStrssRateDip[i]  = 0.0;
	}   }
    /*-------------------------------------------------------------------*/
	/* now the main preparation work is done, define remaining variables before going into the EQ iteration loops*/    
	int             iTemp,          iTemp1; 
	float           fTimeYears;
	unsigned int    iEQcounter      = 0;
	unsigned int    iMaxMRFlength   = 10000;
    unsigned int    iMaxSTFlength   = 2000; 
    float           fMDg_CutStrss   = 2000.0; /* stress in Pa ==> earth tides change stress by about 5000Pa (= 0.005MPa); atmospheric pressure is ~100kPa => 2000Pa == 2% of atmospheric pressure change */
    float           fTemp,              fTemp1,                 fTemp2;
    float           fPostSeisIncrStk,   fPostSeisIncrDip;               
    /*---------------------------------------*/
	unsigned int *iMDl_EQongoing;     iMDl_EQongoing    = (unsigned int *) calloc(1 ,sizeof(unsigned int));	    				
	unsigned int *iMDg_EQongoing;     iMDg_EQongoing    = (unsigned int *) calloc(1 ,sizeof(unsigned int));		
	unsigned int *iMDl_WillContin;    iMDl_WillContin   = (unsigned int *) calloc(1 ,sizeof(unsigned int));	
	unsigned int *iMDg_WillContin;    iMDg_WillContin   = (unsigned int *) calloc(1 ,sizeof(unsigned int));	
	unsigned int *iMDl_ExitEQ;        iMDl_ExitEQ       = (unsigned int *) calloc(1, sizeof(unsigned int));						
    unsigned int *iMDg_ExitEQ;        iMDg_ExitEQ       = (unsigned int *) calloc(1, sizeof(unsigned int));
    unsigned int *iEQl_ActPtchNum;	  iEQl_ActPtchNum   = (unsigned int *) calloc(1, sizeof(unsigned int));						
	unsigned int *iEQg_ActPtchNum;	  iEQg_ActPtchNum   = (unsigned int *) calloc(1, sizeof(unsigned int));	
    unsigned int *iEQl_MRFlength;     iEQl_MRFlength    = (unsigned int *) calloc(1, sizeof(unsigned int));
	unsigned int *iEQg_MRFlength;     iEQg_MRFlength    = (unsigned int *) calloc(1, sizeof(unsigned int));  
    unsigned int *iMDl_TotalRuptT;    iMDl_TotalRuptT   = (unsigned int *) calloc(1, sizeof(unsigned int));	
    unsigned int *iMDl_ChangedStab;   iMDl_ChangedStab  = (unsigned int *) calloc(1, sizeof(unsigned int));
    unsigned int *iMDg_ChangedStab;   iMDg_ChangedStab  = (unsigned int *) calloc(1, sizeof(unsigned int));
    
    float        *fEQl_SeisPot;       fEQl_SeisPot      = (float *)        calloc(1, sizeof(float));	
	float        *fEQg_SeisPot;       fEQg_SeisPot      = (float *)        calloc(1, sizeof(float));
    
    unsigned int *iEQg_WrtStartPos;   iEQg_WrtStartPos  = (unsigned int *) calloc(iSIZE, sizeof(unsigned int));
     
    unsigned int *iTDl_WasActivated;  iTDl_WasActivated = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int)); 
	unsigned int *iTDl_t0;            iTDl_t0           = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int)); 	
	unsigned int *iTDl_nVal;          iTDl_nVal         = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int));		  
	unsigned int *iTDl_nPos;          iTDl_nPos         = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int)); 			
    unsigned int *iTDl_PosOfLastSlp;  iTDl_PosOfLastSlp = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int));  
	unsigned int *iEQl_ActPtchID;    iEQl_ActPtchID     = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int));					
	unsigned int *iEQl_t0ofPtch;     iEQl_t0ofPtch	    = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int));	
	unsigned int *iEQl_StabType;     iEQl_StabType      = (unsigned int *) calloc(iOFFSET[iRANK], sizeof(unsigned int));
	
	 
    float *fTDl_LoclStrssH;          fTDl_LoclStrssH   = (float *)         calloc(iOFFSET[iRANK], sizeof(float)); 				
	float *fTDl_LoclStrssV;          fTDl_LoclStrssV   = (float *)         calloc(iOFFSET[iRANK], sizeof(float));     
	float *fTDl_LoclStrssN;          fTDl_LoclStrssN   = (float *)         calloc(iOFFSET[iRANK], sizeof(float)); 
    float *fTDl_StrssB4_H;           fTDl_StrssB4_H    = (float *)         calloc(iOFFSET[iRANK], sizeof(float)); 				
	float *fTDl_StrssB4_V;           fTDl_StrssB4_V    = (float *)         calloc(iOFFSET[iRANK], sizeof(float));     
	float *fTDl_StrssB4_N;           fTDl_StrssB4_N    = (float *)         calloc(iOFFSET[iRANK], sizeof(float)); 	
	float *fTDl_EventSlipH;          fTDl_EventSlipH   = (float *)         calloc(iOFFSET[iRANK], sizeof(float));     
	float *fTDl_EventSlipV;          fTDl_EventSlipV   = (float *)         calloc(iOFFSET[iRANK], sizeof(float)); 
	float *fEQl_SlipHofPtch;         fEQl_SlipHofPtch  = (float *)         calloc(iOFFSET[iRANK], sizeof(float));    		
	float *fEQl_SlipVofPtch;         fEQl_SlipVofPtch  = (float *)         calloc(iOFFSET[iRANK], sizeof(float));
	float *fEQl_DeltTofPtch;         fEQl_DeltTofPtch  = (float *)         calloc(iOFFSET[iRANK], sizeof(float));
	
	float *fTDl_STF_H;               fTDl_STF_H        = (float *)         calloc(iOFFSET[iRANK]*iMaxSTFlength, sizeof(float)); 
	float *fTDl_STF_V;               fTDl_STF_V        = (float *)         calloc(iOFFSET[iRANK]*iMaxSTFlength, sizeof(float)); 
		
    float *fEQl_MRFvals;             fEQl_MRFvals      = (float *)         calloc(iMaxMRFlength, sizeof(float));
    float *fEQg_MRFvals;             fEQg_MRFvals      = (float *)         calloc(iMaxMRFlength, sizeof(float));			
		
	float *fTDg_EQStressH;           fTDg_EQStressH    = (float *)         calloc(iPatchNum, sizeof(float)); 				
	float *fTDg_EQStressV;           fTDg_EQStressV    = (float *)         calloc(iPatchNum, sizeof(float));     
	float *fTDg_EQStressN;           fTDg_EQStressN    = (float *)         calloc(iPatchNum, sizeof(float)); 
    float *fTDg_EQStressH2;          fTDg_EQStressH2   = (float *)         calloc(iPatchNum, sizeof(float)); 				
	float *fTDg_EQStressV2;          fTDg_EQStressV2   = (float *)         calloc(iPatchNum, sizeof(float));     
	float *fTDg_EQStressN2;          fTDg_EQStressN2   = (float *)         calloc(iPatchNum, sizeof(float)); 
    
    float *fTDg_PostEQStressH;       fTDg_PostEQStressH   = (float *)      calloc(iPatchNum, sizeof(float)); 				
	float *fTDg_PostEQStressV;       fTDg_PostEQStressV   = (float *)      calloc(iPatchNum, sizeof(float));     
	float *fTDg_PostEQStressN;       fTDg_PostEQStressN   = (float *)      calloc(iPatchNum, sizeof(float)); 
    float *fTDg_PostEQStressH2;      fTDg_PostEQStressH2  = (float *)      calloc(iPatchNum, sizeof(float)); 				
	float *fTDg_PostEQStressV2;      fTDg_PostEQStressV2  = (float *)      calloc(iPatchNum, sizeof(float));     
	float *fTDg_PostEQStressN2;      fTDg_PostEQStressN2  = (float *)      calloc(iPatchNum, sizeof(float)); 
    /*-------------------------------------------------------------------*/  
    free(fMDg_AddNrmStrss);     free(iTDg_V1temp);          free(iTDg_V2temp);          free(iTDg_V3temp); 
    free(iTDg_StabTypetemp);    free(iTDg_SegIDtemp);       free(fTDg_StrsRteStktemp);  free(fTDg_StrsRteDiptemp);
    free(fMDg_Vp);              free(fMDg_Poisson);         free(fMDg_MedDense);        free(fMDg_Lambda);
    free(fVDg_Epostemp);        free(fVDg_Npostemp);        free(fTDg_StatFrictemp);    free(fVDg_Zpostemp); 
    free(fTDg_DynFrictemp);     free(iTempValues1);         free(iTempValues2);         
    /*-------------------------------------------------------------------*/
    MPI_File fp_MPIout1;
	/*-------------------------------------------------------------------*/
    MPI_File_delete(cFileN_Out1, MPI_INFO_NULL);
	MPI_File_open(MPI_COMM_WORLD, cFileN_Out1, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
	if (iRANK == 0)			
	{	MPI_File_write(fp_MPIout1, &iEQcounter,   1, MPI_UNSIGNED, &status);
		MPI_File_write(fp_MPIout1, &iUsdGrid,     1, MPI_UNSIGNED, &status);	
		MPI_File_write(fp_MPIout1, &iSrcGrid,     1, MPI_UNSIGNED, &status);		
		MPI_File_write(fp_MPIout1, fMDg_ShearMod, 1, MPI_FLOAT,    &status);	
		MPI_File_write(fp_MPIout1, &fdeltTincr,   1, MPI_FLOAT,    &status);	
	}
	MPI_File_close(&fp_MPIout1);	
	/*-------------------------------------------------------------------*/
	for (j = iPatchNum; j--;    )
    {   fTempValues1[j]   = 0.0;                fTempValues2[j]   = 0.0;
	}	
	for (i = iOFFSET[iRANK]; i--;     ) 
    {	fTDl_StrssB4_H[i] = 0.0;			    fTDl_StrssB4_V[i] = 0.0;					fTDl_StrssB4_N[i] = fTDl_RefNormStrss[i]; /* STF start position (nVal) is always the same => define here...; then nPos is only the relative location to that point */ 
        iTDl_nPos[i]      = 0;			        iTDl_nVal[i]      = i*iMaxSTFlength;    /* is the same as corresponding function below //start position in matrix/vector  */    
    }    
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
    fTimeYears = 0.0;
    while (fTimeYears < fRecLgth)
    {	
		fTimeYears       += fIntSeisTimeStp;
		iMDl_EQongoing[0] = 0;
		/*-------------------------------------------------------------------*/
        if (iUsePostSeismic == 1) /*  https://en.wikipedia.org/wiki/Stress_relaxation    http://web.mit.edu/course/3/3.11/www/modules/visco.pdf      */
		{   for (j = iPatchNum; j--;    )
            {   fTempValues1[j]   = 0.0;                                                fTempValues2[j]   = 0.0;
	        } 
		    for (i = iOFFSET[iRANK]; i--;    ) 
            {	if ((iTDl_StabType[i] != 1) && (fTDl_PostStrssat_t0[i] > FLT_EPSILON)) /* don't use 0.0 but epsilon of float!*/
		        {   
		            fTemp1                 = sqrtf (fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);
    			    fTemp                  = fTemp1 - fTDl_CurrFric[i]*fTDl_StrssB4_N[i];
		            fTemp2                 = fTDl_PostStrssat_t0[i]*expf(-(fTimeInSecs*fTDl_PostStrssTime[i])/fTDl_CharRelaxTime[i]);
		            fTDl_PostStrssTime[i]  = fTDl_PostStrssTime[i] +1.0;
		              
		            if (fTemp > fTemp2)
		            {   fPostSeisIncrStk   = (fTemp-fTemp2)/fTemp *fTDl_StrssB4_H[i];
		                fPostSeisIncrDip   = (fTemp-fTemp2)/fTemp *fTDl_StrssB4_V[i];
		                fTemp1             = -1.0*fPostSeisIncrStk/fK_SS[iTDl_LocSelfInd[i]]; /* slip in horizontal */
		                fTemp2             = -1.0*fPostSeisIncrDip/fK_DD[iTDl_LocSelfInd[i]]; /* slip in vertical  */
                        
		                for (j = iPatchNum; j--;   ) /* here is then the relaxation i.e., stress redistribution according to the amount of loading that will be used this time step */
				        {   iVectPos      = i*iPatchNum + j;		
        	                fTempValues1[j] += fTemp1*fK_SS[iVectPos] +fTemp2*fK_DS[iVectPos];
		                    fTempValues2[j] += fTemp1*fK_SD[iVectPos] +fTemp2*fK_DD[iVectPos];
            }   }   }   } 
            MPI_Allreduce(fTempValues1, fTDg_EQStressH, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    MPI_Allreduce(fTempValues2, fTDg_EQStressV, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    	    
    	    for (i = iOFFSET[iRANK]; i--;    )
            {	fTDl_StrssB4_H[i] += fTDg_EQStressH[i+iSTARTPOS[iRANK]];                fTDl_StrssB4_V[i] += fTDg_EQStressV[i+iSTARTPOS[iRANK]];	
            }      
        }   
		/*-------------------------------------------------------------------*/
		for (i = iOFFSET[iRANK]; i--;    ) 
        {	if (iTDl_StabType[i] == 1 )
        	{   fTDl_StrssB4_H[i] += fTDl_CurStrssRateStk[i]*fIntSeisTimeStp;
    			fTDl_StrssB4_V[i] += fTDl_CurStrssRateDip[i]*fIntSeisTimeStp;
    			fTemp1             = sqrtf (fTDl_StrssB4_H[i]*fTDl_StrssB4_H[i] +fTDl_StrssB4_V[i]*fTDl_StrssB4_V[i]);
    			fTemp              = fTemp1 - fTDl_CurrFric[i]*fTDl_StrssB4_N[i]; 
    			
       			if (fTemp >= fMDg_CutStrss)
    			{	iMDl_EQongoing[0]     = 1;		iTDl_WasActivated[i]  = 1;          iTDl_t0[i]            = 0;	  
        }   }   }
        MPI_Allreduce(iMDl_EQongoing,  iMDg_EQongoing, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    	/*-------------------------------------------------------------------*/
    	if (iMDg_EQongoing[0] > 0)	/* EARTHQUAKE STARTS */	
    	{		
    		iMDg_ExitEQ[0]     = 0;             iMDl_TotalRuptT[0] = -1;                iEQcounter++;	
    	    iMDl_WillContin[0] = 1;
    		for (i = iOFFSET[iRANK]; i--;    ) /* each core goes through all patches and determines minimum time to failure */
            {   fTDl_LoclStrssH[i] = fTDl_StrssB4_H[i];         fTDl_LoclStrssV[i] = fTDl_StrssB4_V[i];             fTDl_LoclStrssN[i] = fTDl_StrssB4_N[i];	
    		}
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
    		while ((iMDg_EQongoing[0] == 1) && (iMDg_ExitEQ[0] == 0)) /* EARTHQUAKE ITERATION LOOP STARTS */
    		{   	
				iMDl_TotalRuptT[0]++;               iMDl_EQongoing[0]      = 0;         iMDl_ExitEQ[0]         = 0;          
    			
    			for (i = iPatchNum;  i--;  )									
    			{	fTDg_EQStressH[i]      = 0.0;	fTDg_EQStressV[i]      = 0.0;	    fTDg_EQStressN[i]      = 0.0;	
    			 //   fTDg_EQStressH2[i]     = 0.0;	fTDg_EQStressV2[i]     = 0.0;	    fTDg_EQStressN2[i]     = 0.0;						
    			}
    			for (i = iOFFSET[iRANK]; i--;    )
                {	/* determine new friction coefficient */
                    if (iTDl_WasActivated[i] == 1)
    				{	
    					if      (iSDg_FricLawUSED[ iTDl_SegID[ i ] ] == 1) /* classic static/dynamic friction */
    					{	
    					    fTDl_CurrFric[i] = fTDl_DynFric[i]; 
    					}
    					else if (iSDg_FricLawUSED[ iTDl_SegID[ i ] ] == 2) /* slip weakening -using accumulated slip of current event (on that patch) */
    					{   
    					    fTemp            = sqrtf(fTDl_EventSlipH[i]*fTDl_EventSlipH[i] +fTDl_EventSlipV[i]*fTDl_EventSlipV[i]);
    						fTDl_CurrFric[i] = fTDl_DynFric[i] + (fTDl_StatFric[i]- fTDl_DynFric[i])* (1.0 - fTemp/fTDl_Curr_DcVal[i]); /* the first part is the current friction then the second part defines delta_mue => change in friction coefficient; and the second one is the  */
    						if      ((iTDl_StabType[i] != 3) && (fTDl_CurrFric[i] < fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
    						else if ((iTDl_StabType[i] == 3) && (fTDl_CurrFric[i] > fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
    					}
    					else if (iSDg_FricLawUSED[ iTDl_SegID[ i ] ] == 3) /* velocity weakening -using slip increment from next iterative slip (with deltT that equates to velocity) */
    					{   
    					    fTemp            = sqrtf(fTDl_STF_H[iTDl_nVal[i]]*fTDl_STF_H[iTDl_nVal[i]] + fTDl_STF_V[iTDl_nVal[i]]*fTDl_STF_V[iTDl_nVal[i]]);
							fTDl_CurrFric[i] = fTDl_DynFric[i] + (fTDl_StatFric[i]- fTDl_DynFric[i])* (1.0 - fTemp/fTDl_Curr_DcVal[i]); 
    						if      ((iTDl_StabType[i] != 3) && (fTDl_CurrFric[i] < fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }
    						else if ((iTDl_StabType[i] == 3) && (fTDl_CurrFric[i] > fTDl_DynFric[i]))     {      fTDl_CurrFric[i] =  fTDl_DynFric[i];     }	
    				}   }
    				/*-------------------------------------------------------*/
    				/* determine amount of excess stress (if any available) */

					fTemp  = sqrtf(fTDl_LoclStrssH[i]*fTDl_LoclStrssH[i] + fTDl_LoclStrssV[i]*fTDl_LoclStrssV[i]);  /* this is the applied shear stress*/
					fTemp1 = fTemp - fTDl_CurrFric[i]*fTDl_LoclStrssN[i]; /* this is the amount of stress above current friction level in combination with normal stress at the location*/					
					fTemp2 = (fTDl_DynFric[i] *fTDl_LoclStrssN[i] - fTemp)/(fK_SS[iTDl_LocSelfInd[i]] > fK_DD[iTDl_LocSelfInd[i]] ? fK_SS[iTDl_LocSelfInd[i]] : fK_DD[iTDl_LocSelfInd[i]]);//(0.5*(fK_SS[iTDl_LocSelfInd[i]] + fK_DD[iTDl_LocSelfInd[i]])); /* this is releasable stress divided by self-stiffness => gives slip amount that would happen if patch fails "alone"....			take average of DD and SS stiffness */					

					if (((iTDl_WasActivated[i] == 1) && (fTemp1 >= fMDg_CutStrss)) || ((iTDl_WasActivated[i] == 0)  && (fTemp1 >= fMDg_CutStrss) && (fTemp2 >= fTDl_Curr_DcVal[i])))
					{	iMDl_WillContin[0]= 0;
					    iMDl_EQongoing[0] = 1;
						if (iTDl_WasActivated[i] == 0)		    {	iTDl_WasActivated[i] = 1;	        iTDl_t0[i] = iMDl_TotalRuptT[0];		        }						
						if (iTDl_nPos[i] < iMaxMRFlength)
						{	iTemp                = iTDl_nVal[i]+iTDl_nPos[i];
						    iTDl_PosOfLastSlp[i] = iTDl_nPos[i];
							fTDl_STF_H[iTemp]    = -1.0*(fTemp1/fTemp *fTDl_LoclStrssH[i]) /fK_SS[iTDl_LocSelfInd[i]]; /* horizontal slip at patch in current iteration (turn "n") */
							fTDl_STF_V[iTemp]    = -1.0*(fTemp1/fTemp *fTDl_LoclStrssV[i]) /fK_DD[iTDl_LocSelfInd[i]]; /* the "-1" is here because the stiffness matrix stuff gives the amount of slip nessessary to MAKE the observed stress; but I want to RELEASE it => opposite direction */
							fTDl_EventSlipH[i]  += fTDl_STF_H[iTDl_nVal[i]+iTDl_nPos[i]];				
							fTDl_EventSlipV[i]  += fTDl_STF_V[iTDl_nVal[i]+iTDl_nPos[i]];	
						}
						else
						{   iMDl_ExitEQ[0] = 1; /* define another exit variable => if I get here, "break" i.e., finish earthquake */
					}   }
					/*---------------------------------------------------------*/
					/* determine modified stress due to slip on fault patches */
					if (iTDl_WasActivated[i] == 1) 
    				{   
    				    for (j = iPatchNum; j--;   )
						{   iVectPos = i*iPatchNum + j;							
							iTemp    = iTDl_nPos[i] - iTravTimes[iVectPos];
							if (iTemp >= 0)
							{   iTemp  = (iTemp < iTDl_PosOfLastSlp[i]) ? iTemp : iTDl_PosOfLastSlp[i]; /* this is important! the STF has limited length => only use the part that contains actual information of current earthquake		*/						
								iTemp  = (iTemp < iMaxSTFlength)        ? iTemp : iMaxSTFlength;
								fTemp1 = 0.0;                             fTemp2 = 0.0; 
								
								for (k = (iTemp+1); k--;    )
                                {   fTemp1 += fTDl_STF_H[iTDl_nVal[i]+k];               fTemp2 += fTDl_STF_V[iTDl_nVal[i]+k]; 
								}
								fTDg_EQStressH[j] += fTemp1*fK_SS[iVectPos] +fTemp2*fK_DS[iVectPos];
								fTDg_EQStressV[j] += fTemp1*fK_SD[iVectPos] +fTemp2*fK_DD[iVectPos];
								fTDg_EQStressN[j] += fTemp1*fK_SO[iVectPos] +fTemp2*fK_DO[iVectPos];
						}	}   	
						iTDl_nPos[i]++;
    			}	}
    			/*---------------------------------------------------------------*/
    			/* last was done "locally" now transfer information to all patches*/
    			MPI_Allreduce(fTDg_EQStressH, fTDg_EQStressH2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(fTDg_EQStressV, fTDg_EQStressV2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(fTDg_EQStressN, fTDg_EQStressN2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			MPI_Allreduce(iMDl_EQongoing, iMDg_EQongoing,      1     , MPI_UNSIGNED,MPI_MAX, MPI_COMM_WORLD);
    			MPI_Allreduce(iMDl_WillContin,iMDg_WillContin,     1     , MPI_UNSIGNED,MPI_MAX, MPI_COMM_WORLD);
    			MPI_Allreduce(iMDl_ExitEQ,    iMDg_ExitEQ,         1     , MPI_UNSIGNED,MPI_MAX, MPI_COMM_WORLD);
    			/* combine "stress before EQ" and "EQ generated" changes to get "Currently applied stress" */	
    			for (i = iOFFSET[iRANK]; i--;    )
                {	
					fTDl_LoclStrssH[i] = fTDl_StrssB4_H[i] + fTDg_EQStressH2[i+iSTARTPOS[iRANK]];
        			fTDl_LoclStrssV[i] = fTDl_StrssB4_V[i] + fTDg_EQStressV2[i+iSTARTPOS[iRANK]];
        			fTDl_LoclStrssN[i] = fTDl_StrssB4_N[i] + fTDg_EQStressN2[i+iSTARTPOS[iRANK]];
                	fTDl_LoclStrssN[i] = (fTDl_LoclStrssN[i] > fTDl_RefNormStrss[i]*0.5) ? fTDl_LoclStrssN[i] : fTDl_RefNormStrss[i]*0.5;
                	// fTDl_LoclStrssN[i] = (fTDl_LoclStrssN[i] >             0.0)          ? fTDl_LoclStrssN[i] : 0.0;
                	/* the later is/are just a "safety measure" to ensure that normal stress is not getting too small... -but maybe it is not needed... => important is that stress is not turning negative*/
        		}
        		/*---------------------------------------------------------------*/
        		if ((iMDg_EQongoing[0] == 0) && (iMDg_WillContin[0] == 0) &&(fdeltTincr < 1.0E+9)) /* make sure that all the signal is out of the system => even if no slip occurred in last step on any patch; there may still be stress in the system that has not reached a receiver => wait until they all got their share */
        		{	float fTempH,       fTempV,         fTempN;
        		    
        		    for (i = iPatchNum;  i--;  )									
    		    	{	fTDg_PostEQStressH[i]  = 0.0;	fTDg_PostEQStressV[i]  = 0.0;		fTDg_PostEQStressN[i]  = 0.0;								
    			    }
        		    /* this here is done (with event slip etc) to figure out if rupture is actually able to continue i.e., if I have to keep working on the loop(s) above or can abort now... */        		    
        		    for (i = iOFFSET[iRANK]; i--;    )
                    {   if (iTDl_WasActivated[i] == 1) 
    				    {   for (j = iPatchNum; j--;   )
						    {   iVectPos               = i*iPatchNum + j;													    
								fTDg_PostEQStressH[j] += fTDl_EventSlipH[i]*fK_SS[iVectPos] +fTDl_EventSlipV[i]*fK_DS[iVectPos];
								fTDg_PostEQStressV[j] += fTDl_EventSlipH[i]*fK_SD[iVectPos] +fTDl_EventSlipV[i]*fK_DD[iVectPos];
								fTDg_PostEQStressN[j] += fTDl_EventSlipH[i]*fK_SO[iVectPos] +fTDl_EventSlipV[i]*fK_DO[iVectPos];
					}	}	}   	
					MPI_Allreduce(fTDg_PostEQStressH, fTDg_PostEQStressH2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			    MPI_Allreduce(fTDg_PostEQStressV, fTDg_PostEQStressV2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    			    MPI_Allreduce(fTDg_PostEQStressN, fTDg_PostEQStressN2, iPatchNum , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
					
					for (i = iOFFSET[iRANK]; i--;    )
                    {   fTempH  = fTDl_StrssB4_H[i] +fTDg_PostEQStressH2[i];
                        fTempV  = fTDl_StrssB4_V[i] +fTDg_PostEQStressV2[i];
                        fTempN  = fTDl_StrssB4_N[i] +fTDg_PostEQStressN2[i];
                        
                        fTemp   = sqrtf(fTempH*fTempH + fTempV*fTempV);
					    fTemp1  = fTemp - fTDl_CurrFric[i]*fTempN;
					    fTemp2  =(fTDl_DynFric[i]*fTempN - fTemp)/(fK_SS[iTDl_LocSelfInd[i]] > fK_DD[iTDl_LocSelfInd[i]] ? fK_SS[iTDl_LocSelfInd[i]] : fK_DD[iTDl_LocSelfInd[i]]);//(0.5*(fK_SS[iTDl_LocSelfInd[i]] + fK_DD[iTDl_LocSelfInd[i]])); /* this is releasable stress divided by self-stiffness => gives slip amount that would happen if patch fails "alone" .... */					
	 
					    if (((iTDl_WasActivated[i] == 1) && (fTemp1 >= fMDg_CutStrss)) || ((iTDl_WasActivated[i] == 0) && (fTemp2 >= fTDl_Curr_DcVal[i]))) /* if this if-statement is satisfied (there is sufficient stress still available) then I will continue to check */
	                    {	iMDl_EQongoing[0]  = 1;     
	                        iMDl_WillContin[0] = 1; 
					}   }	
					MPI_Allreduce(iMDl_EQongoing,  iMDg_EQongoing,      1  ,    MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);					
				    MPI_Allreduce(iMDl_WillContin, iMDg_WillContin,     1  ,    MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);					    
				}
    		} /* EARTHQUAKE ITERATION LOOP EXITED */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */
/*------------------------------------------------------------------ */   		
    		iEQl_ActPtchNum[0] = 0;
    		for (i = iOFFSET[iRANK]; i--;   )   /* now I define all the metrics of the event (its size, its MRF, etc...) */
            {   if (iTDl_WasActivated[i] == 1) 	
    			{	
    			    for (k = iTDl_PosOfLastSlp[i]; k--;    )       
    			    {	fEQl_MRFvals[iTDl_t0[i]+k] += sqrtf( fTDl_STF_H[iTDl_nVal[i]+k]*fTDl_STF_H[iTDl_nVal[i]+k] + fTDl_STF_V[iTDl_nVal[i]+k]*fTDl_STF_V[iTDl_nVal[i]+k]) *fTDl_Area[i];			}
    				
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
    		MPI_File_open(MPI_COMM_WORLD, cFileN_Out1,MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_MPIout1);
    		MPI_File_get_size(fp_MPIout1, &offset1);
    		MPI_Barrier( MPI_COMM_WORLD );
    			
			if (iRANK == 0)
			{	
			    fTemp2 = (log10f(fEQg_SeisPot[0]*fMDg_ShearMod[0])-9.1)/1.5; 
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
            iMDl_ChangedStab[0] = 0; /* this will record if stability type was changed*/
            for (i = iOFFSET[iRANK]; i--;   ) /* RESET ALL VALUES FOR PATCHES THAT MOVED */
            {
                /*------------------------------------------------------------- */
                if (iTDl_WasActivated[i] == 1) 	
    			{   if (iMDg_ChgFricBtwEQs[0] == 0)
		 		    {	fTDl_CurrFric[i]     = fTDl_StatFric[i];
		 		    }
		 		    else   
		 		    {	
		 		        fTDl_StatFric[i]   = fSDg_RefStatFric[iTDl_SegID[i]]  *(1.0 +fRandVector[iRandPos[0]] *fSDg_RefStatFr_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }       
                        fTDl_CurrFric[i]   = fTDl_StatFric[i];
                        fTDl_DynFric[i]    = fSDg_RefDynFric[iTDl_SegID[i]]   *(1.0 +fRandVector[iRandPos[0]] *fSDg_RefDynFr_vari[ iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }  
		 			    fTDl_Curr_DcVal[i] = fSDg_CritSlipDist[iTDl_SegID[i]] *(1.0 +fRandVector[iRandPos[0]] *fSDg_CritSlipD_vari[iTDl_SegID[i]]/100.0);       iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }      
		 		        /* test now again for stability type, while doing that check if any type was changed, if not then I don't have to recompute, if change is from cond stable to stable then I only change label, only if it goes from unstable to another state or vice versa, is it that I need to update loading*/

		 		        fTemp1 = (fTDl_DynFric[i] - fTDl_StatFric[i]) *fTDl_RefNormStrss[i]; /* this is the releasable stress amount*/
                        fTemp2 = fTemp1/(fK_SS[iTDl_LocSelfInd[i]] > fK_DD[iTDl_LocSelfInd[i]] ? fK_SS[iTDl_LocSelfInd[i]] : fK_DD[iTDl_LocSelfInd[i]]);//(0.5*(fK_SS[iTDl_LocSelfInd[i]] + fK_DD[iTDl_LocSelfInd[i]])); /* this is corresponding slip, sign is correct -> assuming that K-matrix at self-induced location is ALWAYS negative*/
        
                        if (fTemp2 > fTDl_Curr_DcVal[i]) /* have stress drop when going from stat to dyn friction and that drop i.e., corresponding slip is larger than Dc*/
                        {	if (iTDl_StabType[i] != 1)              {   iMDl_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
                            iTDl_StabType[i] = 1; /* patch is unstable*/
                        }
                        else /*if (fTemp2 <= fTDl_Curr_DcVal[i]) -if stress drop i.e., corresponding slip is smaller than Dc - */
                        {	if (fTemp1 < 0.0)  /*  is still weakening but just with slip that is lower than Dc => cond. stable*/
        	                {	if (iTDl_StabType[i] == 1)          {   iMDl_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
        	                    iTDl_StabType[i] = 2; /* patch is cond. stable*/
        	                }
        	                else /* no weakening but strengthening  => stable */
        	                {	if (iTDl_StabType[i] == 1)          {   iMDl_ChangedStab[0] = 1;        iTDl_StabTypeChgd[i]  = 1;      }
        	                    iTDl_StabType[i] = 3; /* patch is stable*/
                    }   }   }   
                }
                /*------------------------------------------------------------- */
		 		for (k = iMaxSTFlength; k--;   )     {	fTDl_STF_H[iTDl_nVal[i]+k] = 0.0;	    fTDl_STF_V[iTDl_nVal[i]+k] = 0.0;		} 	
		 		
		 		iEQl_ActPtchID[i]    = 0;               iEQl_t0ofPtch[i]     = 0;		        iTDl_WasActivated[i] = 0;				
		 		iTDl_t0[i]           = 0;               iTDl_nPos[i]         = 0;               iTDl_PosOfLastSlp[i] = 0;
		 		fEQl_DeltTofPtch[i]  = 0.0;             fEQl_SlipHofPtch[i]  = 0.0;			    fEQl_SlipVofPtch[i]  = 0.0;
    			fTDl_EventSlipH[i]   = 0.0;			    fTDl_EventSlipV[i]   = 0.0;
    		}
    		iEQl_ActPtchNum[0] 	= 0;		            iEQg_ActPtchNum[0] 	 = 0;               iEQl_MRFlength[0]    = 0;						
    		iEQg_MRFlength[0]   = 0;                    iMDl_ExitEQ[0]       = 0;		        iMDg_ExitEQ[0]       = 0;    
    		iMDl_EQongoing[0]   = 0;                    iMDg_EQongoing[0]    = 0;               iMDl_WillContin[0]   = 0;        
    		iMDg_WillContin[0]  = 0;                    fEQl_SeisPot[0]      = 0.0;		        fEQg_SeisPot[0]      = 0.0;
    		
    		for (i = iSIZE; i--;    )				{	iEQg_WrtStartPos[i] = 0;											    }
    		for (i = iMaxMRFlength; i--;    )	    {	fEQl_MRFvals[i]     = 0.0; 			    fEQg_MRFvals[i]      = 0.0;		}
            /*------------------------------------------------------------------------- */
            MPI_Allreduce(iMDl_ChangedStab, iMDg_ChangedStab,       1          ,    MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    		      	
		 	if (iMDg_ChangedStab[0] == 1) /* next, update the current stressing rate*/
		 	{   /*-------------------------------------------------------------------*/
	            for (j = iPatchNum; j--;    )
                {   fTempValues1[j]   = 0.0;            fTempValues2[j]   = 0.0;
	            }
	            for (i = iOFFSET[iRANK]; i--;   )
                {		
    	            if ((iTDl_StabType[i] != 1) && ( iTDl_StabTypeChgd[i] == 1))/* if the patch is NOT unstable but instead stable or conditionally stable ==> then it will creep*/
    	            {	fTempSlipStk =  fTDl_RefStrssRateStk[i]/fK_SS[iTDl_LocSelfInd[i]]; /* this is proposed creep (in X/yr) on the not-stick-slipping patches */
    		            fTempSlipDip =  fTDl_RefStrssRateDip[i]/fK_DD[iTDl_LocSelfInd[i]];
    
    		            for (j = iPatchNum; j--;    )
    		            {	if (iTDl_StabType[i] == 1) /* apply the creep slip to load the  unstable patches/cells; if static friction is larger than dyn -> is a stick-slip patch */	
    			            {   iVectPos                = i*iPatchNum + j;  				
    				            fTempValues1[j] += fTempSlipStk*fK_SS[iVectPos] + fTempSlipDip*fK_DS[iVectPos];
    				            fTempValues2[j] += fTempSlipStk*fK_SD[iVectPos] + fTempSlipDip*fK_DD[iVectPos];	    
	            }   }   }   }
                /*-------------------------------------------------------------------*/
                MPI_Allreduce(fTempValues1, fTDg_StrssRteChgStk, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); /* here, tempvalues contains the changes in stress on locked patches due to slip/creep on not-locked ones */
                MPI_Allreduce(fTempValues2, fTDg_StrssRteChgDip, iPatchNum, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);   
                for (i = iPatchNum; i--;    )		{	fTempValues1[i] = 0.0;        fTempValues2[i] = 0.0;                 } 
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
		                
		                fTemp = fTDl_PostStrssat_t0[i]*expf(-(fTimeInSecs*fTDl_PostStrssTime[i])/fTDl_CharRelaxTime[i]);/* is going to be the post-seis "strength" at current time -relative to PostStrssTime*/
		                
		                if ((fTemp1 - fTemp2) > fTemp)
		                {   fTDl_PostStrssat_t0[i] = fTemp1 - fTemp2;  
		                    fTDl_PostStrssTime[i]  = 1.0;              
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