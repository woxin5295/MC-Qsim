# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void LoadInputParameter(char **argv, unsigned int *iMDg_ChgFricBtwEQs, float *fMDg_AddNrmStrss, unsigned int *iSDg_FricLawUSED, unsigned int *iTDg_V1temp, unsigned int *iTDg_V2temp, unsigned int *iTDg_V3temp, unsigned int *iTDg_StabTypetemp, unsigned int *iTDg_SegIDtemp, float *fTDg_StrsRteStktemp, float *fTDg_StrsRteDiptemp, float *fMDg_Vp, float *fMDg_Poisson, float *fMDg_Lambda, float *fMDg_ShearMod, float *fMDg_MedDense,  float *fSDg_RefStatFric, float *fSDg_RefStatFr_vari, float *fSDg_RefDynFric, float *fSDg_RefDynFr_vari, float *fSDg_CritSlipDist, float *fSDg_CritSlipD_vari, float *fSDg_CritSlipVelo, float *fTDg_CentEpos, float *fTDg_CentNpos, float *fTDg_CentZpos, float *fTDg_StatFrictemp, float *fTDg_DynFrictemp, float *fVDg_Epostemp, float *fVDg_Npostemp, float *fVDg_Zpostemp, unsigned int iSegmentNum, unsigned int iPatchNum, unsigned int iVertexNum, unsigned int iUsdGrid, unsigned int iSrcGrid)							
{
	/*-------------------------------------------*/
	char 	cFileName1[512];		strcpy(cFileName1,argv[1]);       	strcat(cFileName1,"_Summary_dFricRoughn.nfo");
    char 	cFileName2[512];		strcpy(cFileName2,argv[1]);         strcat(cFileName2,"_");								strcat(cFileName2,argv[2]);			strcat(cFileName2,"_Roughn.dat");
    char 	cFileName3[512];		strcpy(cFileName3,argv[1]);         strcat(cFileName3,"_");								strcat(cFileName3,argv[2]);			strcat(cFileName3,"_Strgth.dat");
    char 	cFileName4[512];		strcpy(cFileName4,argv[1]);       	strcat(cFileName4,"_StressRate.dat");
	char    *retch;
	/*-------------------------------------------*/
	int		iTemp1,			iTemp2,			i,			j;
	size_t  ret;
	
	float   fTemp1;
	char    ctempVals[512];
	FILE    *fp1,			*fp2; 
	/*---------------------------------------------------------------------------------*/
	if ((fp1 = fopen(cFileName1,"r")) == NULL)       {   printf("Error -cant open *.flt file. LoadInputParameter function...\n");      exit(10);     }
    
    retch =fgets(ctempVals, 512, fp1);							sscanf(ctempVals,"%*s %e",&fMDg_MedDense[0]);	
    retch =fgets(ctempVals, 512, fp1);							sscanf(ctempVals,"%*s %e",&fTemp1);									fMDg_AddNrmStrss[0] = fTemp1*1.0E+6; /* convert from MPa to Pa*/
    retch =fgets(ctempVals, 512, fp1);							sscanf(ctempVals,"%*s %d",&iMDg_ChgFricBtwEQs[0]);		
	  
	for (i = 0; i < iSegmentNum; i++)
    {	for (j = 0; j < 19;	j++)					{	retch =fgets(ctempVals, 512, fp1);							}
    	retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %d", &iTemp1);								if (iTemp1 == 1)		{	iSDg_FricLawUSED[i] = 1;		}
    	retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_RefStatFric[i]);
    	retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_RefStatFr_vari[i]);
    	retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_RefDynFric[i]);
    	retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_RefDynFr_vari[i]);
        retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %d", &iTemp1);								if (iTemp1 == 1)		{	iSDg_FricLawUSED[i] = 2;		}
  		retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_CritSlipDist[i]);	
  		retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %e", &fSDg_CritSlipD_vari[i]);	
  		retch =fgets(ctempVals, 512, fp1);						sscanf(ctempVals,"%*s %d", &iTemp1);								if (iTemp1 == 1)		{	iSDg_FricLawUSED[i] = 3;		}
  		for (j = 0; j < 12;	j++)					{	retch =fgets(ctempVals, 512, fp1);							}
    }
	for (i = 0; i < iSegmentNum; i++)
    {	fSDg_CritSlipVelo[i] = 1.0; /* that is in m/s -don't know right now where this comes in; guess it is used when i want to do velocity weakening */
    }
	fclose(fp1);
	/*---------------------------------------------------------------------------------*/
	if ((fp1 = fopen(cFileName2,"rb")) == NULL)       {   printf("Error -cant open *_Roughn.dat file. LoadInputParameter function...\n");      exit(10);     }
	if ((fp2 = fopen(cFileName3,"rb")) == NULL)       {   printf("Error -cant open *_Strgth.dat file. LoadInputParameter function...\n");      exit(10);     }

	for (i = 1; i < iUsdGrid; i++)              /* all grids before the used one are tossed => read into "dummmy" */
	{	ret =fread(&iTemp1,   sizeof(  int),1,fp1);  /* this is currently read PatchNum  */
		ret =fread(&iTemp2,   sizeof(  int),1,fp1);  /* this is currently read VertexNum */
		
		fseek(fp1, (5L*sizeof(int)*(long)iTemp1 +5L*sizeof(float)*(long)iTemp2 +3L*sizeof(float)*(long)iTemp1), SEEK_CUR); /* skip over the values that I don't need to read */
		fseek(fp2, (1L*sizeof(int)*(long)1      +2L*sizeof(float)*(long)iTemp1 +1L*sizeof(int)*(long)iTemp1),   SEEK_CUR); /* skip over the values that I don't need to read */
	}
	ret =fread(&iTemp1, sizeof( int),1,fp1); 	/* this is currently read PatchNum  */
	ret =fread(&iTemp2, sizeof( int),1,fp1);		/* this is currently read VertexNum */
	if ((iPatchNum != iTemp1) || (iVertexNum != iTemp2))			{	fprintf(stdout,"Patch or vertex number not matching: patchum1 %d  patchum2 %d       vertexnum1 %d    vertexnum2 %d\n", iPatchNum, iTemp1, iVertexNum, iTemp2);		exit(10);		}
	/*-------------------------------------------*/
	ret =fread(iTDg_V1temp, sizeof(int),iPatchNum,fp1);
	ret =fread(iTDg_V2temp, sizeof(int),iPatchNum,fp1);
	ret =fread(iTDg_V3temp, sizeof(int),iPatchNum,fp1);
	/*-------------------------------------------*/
	fseek(fp1, 1L*sizeof(int)*(long)iPatchNum ,     SEEK_CUR); /* this is fault ID => but don't need it => skip over it */
	/*-------------------------------------------*/
	ret =fread(iTDg_SegIDtemp, sizeof(int),iPatchNum,fp1);	
	/*-------------------------------------------*/
	fseek(fp1, 2L*sizeof(float)*(long)iVertexNum,	SEEK_CUR); /* this would be local coordintates (from gridding) of the vertices..; not needed */
	/*-------------------------------------------*/
	ret =fread(fVDg_Epostemp, sizeof(float),iVertexNum,fp1);
	ret =fread(fVDg_Npostemp, sizeof(float),iVertexNum,fp1);
	ret =fread(fVDg_Zpostemp, sizeof(float),iVertexNum,fp1);
	/*-------------------------------------------*/
	ret =fread(fTDg_CentEpos, sizeof(float),iPatchNum,fp1);
	ret =fread(fTDg_CentNpos, sizeof(float),iPatchNum,fp1);	
	ret =fread(fTDg_CentZpos, sizeof(float),iPatchNum,fp1);
	
	/*-------------------------------------------*/
	fseek(fp2, 1L*sizeof(int),     SEEK_CUR); /* contains the patch number again -> skipped */
	/*-------------------------------------------*/
	ret =fread(fTDg_StatFrictemp,sizeof(float), iPatchNum,fp2); /* first realization of STATIC friction coefficient */
	ret =fread(fTDg_DynFrictemp, sizeof(float), iPatchNum,fp2); /* first realization of DYNAMIC friction coefficient */
	ret =fread(iTDg_StabTypetemp,  sizeof(int), iPatchNum,fp2); /* if patch is seismogenic or not (INCLUDES depth-dependence due to (a-b) and also fractal strength heterogeneity and unlocked patches !! */
	fclose(fp1);
	fclose(fp2);
	/*-------------------------------------------*/
	for (i = 0; i < iVertexNum; i++)			
	{	fVDg_Epostemp[i] *= 1000.0;				fVDg_Npostemp[i] *= 1000.0;				fVDg_Zpostemp[i] *= 1000.0;				}
	/*-------------------------------------------*/
	for (i = 0; i < iPatchNum; i++)			
	{	iTDg_V1temp[i]    -= 1;				iTDg_V2temp[i]   -= 1;			iTDg_V3temp[i]   -= 1;			
	    iTDg_SegIDtemp[i] -= 1; 		
		fTDg_CentEpos[i]  *= 1000.0;		fTDg_CentNpos[i] *= 1000.0;		fTDg_CentZpos[i] *= 1000.0;	
	}/* the minus 1 is to bring the index to be conform with C standard, starting at "0" */
	/*---------------------------------------------------------------------------------*/
	unsigned int *iUsedGridVect;
	/*---------------------------------------------------------------------------------*/
	if ((fp1 = fopen(cFileName4,"rb")) == NULL)       {   printf("Error -cant open *_StressRate.dat file. LoadInputParameter function...\n");      exit(10);     }
	ret =fread(fMDg_ShearMod, sizeof(float),1,fp1); 		/* this is parameter for elastic half-space "mu" == ShearModulus */
	ret =fread(fMDg_Lambda,   sizeof(float),1,fp1); 		/* this is parameter for elastic half-space "lambda" */
	
	fMDg_Poisson[0] = fMDg_Lambda[0]/(2.0*(fMDg_Lambda[0]+fMDg_ShearMod[0]));
	fTemp1          = (fMDg_MedDense[0] > 0.0) ? fMDg_MedDense[0] :  2700.0; /* if zero density is used for depth gradient then I still!! need to define a density for the wave propagation speed.... => is done here... */
	fMDg_Vp[0]      = sqrtf((fMDg_Lambda[0] +2.0*fMDg_ShearMod[0])/fTemp1);
	      
	ret =fread(&iTemp1, sizeof(int),1,fp1); /* number of source/loading grids */
	ret =fread(&iTemp2, sizeof(int),1,fp1); /* number of receiver/interacting faults/grids */
	if (iTemp1 < iSrcGrid) {   	printf("Error -source grid entered is larger than that in file. LoadInputParameter function...\n");        exit(10);     }
	if (iTemp2 < iUsdGrid) {	printf("Error -receiver grid entered is larger than that in file. LoadInputParameter function...\n");      exit(10);     }				

	iUsedGridVect  = (unsigned int *) calloc(iTemp2, sizeof(unsigned int));
	
	fseek(fp1, 1L*sizeof(int)*(long)iTemp1, SEEK_CUR); /* skip over the number of patches that are in the sourcegrid => not needed */
	ret =fread(iUsedGridVect, sizeof(int),iTemp2,fp1);
	
	for (i = 0; i < iTemp2; i++) 	 /* going through receivers */
	{	for (j = 0; j < iTemp1; j++) /* going through sources */
		{	
			if ( ((i+1) == iUsdGrid) && ((j+1) == iSrcGrid))
			{	ret =fread(fTDg_StrsRteStktemp, sizeof(float), iUsedGridVect[i],fp1);  /* this value iUsedGridVect[i] is the same as iPatchNum -i.e., has to be */
				ret =fread(fTDg_StrsRteDiptemp, sizeof(float), iUsedGridVect[i],fp1); 	
				break;
			}
			else
			{	fseek(fp1, 2L*sizeof(float)*(long)iUsedGridVect[i],     SEEK_CUR);
			}
	}	}

	fclose(fp1);
	/*---------------------------------------------------------------------------------*/
	free(iUsedGridVect);
	/*---------------------------------------------------------------------------------*/
    return;
}
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
void  DefineMoreParas(unsigned int *iTravTimes, unsigned int *iTDl_ttmax,  float *fTDl_Area, float *fTDl_RefNormStrss, float *fTDl_Curr_DcVal, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fVDg_Epostemp, const float *fVDg_Npostemp, const float *fVDg_Zpostemp, const unsigned int *iTDg_V1temp, const unsigned int *iTDg_V2temp, const unsigned int *iTDg_V3temp, const unsigned int *iTDl_SegID, const float *fSDg_CritSlipDist, const float *fSDg_CritSlipD_vari, const int *iSTARTPOS, const int *iOFFSET, int iRANK, float fMDg_Vp, float fMDg_MedDense, float fMDg_AddNrmStrss, float fIntSeisTimeStp,  float fMDg_g, float fdeltTincr, unsigned int iPatchNum, float *fRandVector, unsigned int *iRandPos, unsigned int iRandNumber)
{
	unsigned int   	i,				j,				iVectPos,		globi;
	float 			fTemp;
	float  			fP1[3],         fP2[3],    		fP3[3];
	float 			fP1P2[3],		fP1P3[3],		fP1P2crossP1P3[3];

	for (i = iOFFSET[iRANK]; i--;   )
	{	globi = i+iSTARTPOS[iRANK];
		/*---------------------------------------------------------------------------------*/
		fP1[0] = fVDg_Epostemp[iTDg_V1temp[globi]];        fP1[1] = fVDg_Npostemp[iTDg_V1temp[globi]];        fP1[2] = fVDg_Zpostemp[iTDg_V1temp[globi]];
        fP2[0] = fVDg_Epostemp[iTDg_V2temp[globi]];        fP2[1] = fVDg_Npostemp[iTDg_V2temp[globi]];        fP2[2] = fVDg_Zpostemp[iTDg_V2temp[globi]];
        fP3[0] = fVDg_Epostemp[iTDg_V3temp[globi]];        fP3[1] = fVDg_Npostemp[iTDg_V3temp[globi]];        fP3[2] = fVDg_Zpostemp[iTDg_V3temp[globi]];
       	
		fP1P2[0] = fP2[0] - fP1[0];			fP1P2[1] = fP2[1] - fP1[1];			fP1P2[2] = fP2[2] - fP1[2];
		fP1P3[0] = fP3[0] - fP1[0];			fP1P3[1] = fP3[1] - fP1[1];			fP1P3[2] = fP3[2] - fP1[2];
	
		fP1P2crossP1P3[0] = fP1P2[1]*fP1P3[2] - fP1P2[2]*fP1P3[1];
		fP1P2crossP1P3[1] = fP1P2[2]*fP1P3[0] - fP1P2[0]*fP1P3[2];
		fP1P2crossP1P3[2] = fP1P2[0]*fP1P3[1] - fP1P2[1]*fP1P3[0];

		fTDl_Area[i]      = 0.5*sqrtf( fP1P2crossP1P3[0]*fP1P2crossP1P3[0] +fP1P2crossP1P3[1]*fP1P2crossP1P3[1] +fP1P2crossP1P3[2]*fP1P2crossP1P3[2]);
		/*---------------------------------------------------------------------------------*/
		fTDl_RefNormStrss[i] = (int) (fMDg_MedDense *fMDg_g *fabs(fTDg_CentZpos[globi])) + fMDg_AddNrmStrss; /* all in Pa */
		/*---------------------------------------------------------------------------------*/
		fTDl_Curr_DcVal[i]   = fSDg_CritSlipDist[iTDl_SegID[i]] *(1.0 + fRandVector[iRandPos[0]]*fSDg_CritSlipD_vari[iTDl_SegID[i]]/100.0  );	
		
	
		iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % (iRandNumber-1);     }  
		    
		/*---------------------------------------------------------------------------------*/
		for (j = iPatchNum; j--;   )
		{	iVectPos = i*iPatchNum + j;
    		fTemp                = sqrtf( (fTDg_CentEpos[globi]-fTDg_CentEpos[j])*(fTDg_CentEpos[globi]-fTDg_CentEpos[j]) + (fTDg_CentNpos[globi]-fTDg_CentNpos[j])*(fTDg_CentNpos[globi]-fTDg_CentNpos[j]) + (fTDg_CentZpos[globi]-fTDg_CentZpos[j])*(fTDg_CentZpos[globi]-fTDg_CentZpos[j]) )/fMDg_Vp; /* the distance between both */
        	iTravTimes[iVectPos] = (unsigned int)(fTemp/fdeltTincr);
        	iTDl_ttmax[i]        = (iTDl_ttmax[i] > iTravTimes[iVectPos]) ? iTDl_ttmax[i] : iTravTimes[iVectPos];
        	
    	}
	}
	return;
}
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/