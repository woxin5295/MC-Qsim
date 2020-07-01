# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <float.h>
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void GetLocKOS_inLoadInput(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3]);

void LoadInputParameter(char **argv, const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int iRealizNum,const int iUsdGrid, const int iFltSegmNum, const int iFltPtchNum, const int iBndPtchNum, const int iFltVertNum, const int iBndVertNum, const float *fSD_StrssRate, const float *fSD_SlipRate, const float *fSD_SlipRake,  int *iMD_ChgFricBtwEQs, float *fMD_AddNrmStrss, float *fMD_Vp, float *fMD_Vs, float *fMD_Poisson, float *fMD_Lambda, float *fMD_ShearMod, float *fMD_MedDense,  int *iSD_FricLawUSED,  float *fSD_RefStatFric, float *fSD_RefStatFr_vari, float *fSD_RefDynFric, float *fSD_RefDynFr_vari, float *fSD_CritSlipDist, float *fSD_CritSlipD_vari, int *iTDg_V1,  int *iTDg_V2,  int *iTDg_V3, int *iTDl_SegID, float *fVDg_Epos, float *fVDg_Npos, float *fVDg_Zpos,float *fTDg_CentEpos, float *fTDg_CentNpos, float *fTDg_CentZpos, float *fTDl_RefStatFric, float *fTDl_RefDynFric, float *fTDl_StatFric, float *fTDl_DynFric,int   *iTDl_StabType, float *fTDl_RefStrssRateStk, float *fTDl_RefStrssRateDip, float *fTDl_SlipRate, float *fTDl_SlipRake, float *fTDl_CurrFric)
{
    /*-------------------------------------------*/
    char     cFileName1[512];       strcpy(cFileName1,argv[1]);         strcat(cFileName1,"_Summary_RoughStrength.txt");
    char     cFileName4[512];       strcpy(cFileName4,argv[1]);         strcat(cFileName4,"_BNDtrig.dat");
    char     cFileName3[512];       strcpy(cFileName3,argv[1]);         strcat(cFileName3,"_");                    strcat(cFileName3,argv[2]);            strcat(cFileName3,"_Strgth.dat");  
    char     cFileName2[512];                 
    if      (iRealizNum == 0)   {   strcpy(cFileName2,argv[1]);         strcat(cFileName2,"_FLTtrig.dat");         }
    else if (iRealizNum > 0)    {   strcpy(cFileName2,argv[1]);         strcat(cFileName2,"_");                    strcat(cFileName2,argv[2]);            strcat(cFileName2,"_Roughn.dat");     }

    char    *retch;
    /*-------------------------------------------*/
    int     iTemp1,             iTemp2,         i,            j;
    size_t  ret;    
    float   fTemp1;
    char    ctempVals[512];
    FILE    *fp1,               *fp2,           *fp3; 
    /*---------------------------------------------------------------------------------*/
    if ((fp1 = fopen(cFileName1,"r")) == NULL)         {          printf("Error -cant open *.flt file. LoadInputParameter function...\n");      exit(10);     }
    
    retch =fgets(ctempVals, 512, fp1);                            sscanf(ctempVals,"%*s %e",&fMD_MedDense[0]);    
    retch =fgets(ctempVals, 512, fp1);                            sscanf(ctempVals,"%*s %e",&fMD_AddNrmStrss[0]);                 fMD_AddNrmStrss[0] *= -1.0E+6; /* convert from MPa to Pa; use negative to convert from compressive=positive to compressive=negative*/
    retch =fgets(ctempVals, 512, fp1);                            sscanf(ctempVals,"%*s %e",&fMD_ShearMod[0]);                    fMD_ShearMod[0]    *=  1.0E+9; /* convert from GPa to Pa*/
    retch =fgets(ctempVals, 512, fp1);                            sscanf(ctempVals,"%*s %e",&fMD_Poisson[0]);
    retch =fgets(ctempVals, 512, fp1);                            sscanf(ctempVals,"%*s %d",&iMD_ChgFricBtwEQs[0]);        
     
    fMD_Lambda[0] = (2.0*fMD_ShearMod[0]*fMD_Poisson[0])/(1.0-2.0*fMD_Poisson[0]);
    fTemp1        = (fMD_MedDense[0] > 0.0) ? fMD_MedDense[0] :  2700.0; /* if zero density is used for depth gradient then I still!! need to define a density for the wave propagation speed.... => is done here... */
    fMD_Vp[0]     = sqrtf((fMD_Lambda[0] +2.0*fMD_ShearMod[0])/fTemp1); /* in m/s */
    fMD_Vs[0]     = sqrtf(fMD_ShearMod[0]/fTemp1); /* in m/s */
    /*-------------------------------------------*/    
    for (i = 0; i < iFltSegmNum; i++)
    {    for (j = 0; j < 19;    j++)                    {         retch =fgets(ctempVals, 512, fp1);                            }
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %d", &iSD_FricLawUSED[i]);                                
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_RefStatFric[i]);
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_RefStatFr_vari[i]);
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_RefDynFric[i]);
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_RefDynFr_vari[i]);    
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_CritSlipDist[i]);    
        retch =fgets(ctempVals, 512, fp1);                        sscanf(ctempVals,"%*s %e", &fSD_CritSlipD_vari[i]);            
        for (j = 0; j < 12;    j++)                     {         retch =fgets(ctempVals, 512, fp1);                            }
    }
    fclose(fp1);
    /*---------------------------------------------------------------------------------*/
    if ((fp1 = fopen(cFileName2,"rb"))     == NULL)   {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName2);      exit(10);     }
    if ((fp2 = fopen(cFileName3,"rb"))     == NULL)   {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName3);      exit(10);     }    
    if (iBndPtchNum > 0)
    {   if ((fp3 = fopen(cFileName4,"rb")) == NULL)   {   printf("Error -cant open  %s LoadInputParameter function...\n",cFileName4);     exit(10);     }
    }
    /*-------------------------------------------*/
    for (i = 1; i < iUsdGrid; i++)                   /* all grids before the used one are tossed => read into "dummmy" */
    {   ret =fread(&iTemp1,   sizeof(  int),1,fp1);  /* this is currently read PatchNum  */
        ret =fread(&iTemp2,   sizeof(  int),1,fp1);  /* this is currently read VertexNum */
        
        fseek(fp1, (5L*sizeof(int)*(long)iTemp1 +5L*sizeof(float)*(long)iTemp2 +3L*sizeof(float)*(long)iTemp1), SEEK_CUR); /* skip over the values that I don't need to read */
        fseek(fp2, (1L*sizeof(int)*(long)1      +2L*sizeof(float)*(long)iTemp1 +1L*sizeof(int)  *(long)iTemp1), SEEK_CUR); /* skip over the values that I don't need to read */
    }
    ret =fread(&iTemp1, sizeof( int),1,fp1);         /* this is currently read PatchNum  */
    ret =fread(&iTemp2, sizeof( int),1,fp1);        /* this is currently read VertexNum */
    if ((iFltPtchNum != iTemp1) || (iFltVertNum != iTemp2))            {    fprintf(stdout,"Patch or vertex number not matching: patchum1 %d  patchum2 %d       vertexnum1 %d    vertexnum2 %d\n", iFltPtchNum, iTemp1, iFltVertNum, iTemp2);        exit(10);        }
    /*-------------------------------------------*/
    ret =fread(iTDg_V1, sizeof(int),iFltPtchNum,fp1); /*might be ok.. but need to add an offset when loading the boundary box...*/
    ret =fread(iTDg_V2, sizeof(int),iFltPtchNum,fp1);
    ret =fread(iTDg_V3, sizeof(int),iFltPtchNum,fp1);
    /*-------------------------------------------*/
    fseek(fp1, 1L*sizeof(int)*(long)iFltPtchNum ,     SEEK_CUR); /* this is fault ID => but don't need it => skip over it */
    /*-------------------------------------------*/
    
    fseek(fp1, (1L*sizeof(int)*(long)(iSTARTPOS_F[iRANK])), SEEK_CUR); /*skip the part that is in front of segId for specific RANK*/
    ret =fread(iTDl_SegID, sizeof(int),iOFFSET_F[iRANK],fp1);   /*read the SegID for specific RANK*/
    fseek(fp1, (1L*sizeof(int)*(long)(iFltPtchNum -iSTARTPOS_F[iRANK]-iOFFSET_F[iRANK])), SEEK_CUR); /*skip over to end of SegIDs => in total i will have moved by iFltPtchNum*/
    
    /*-------------------------------------------*/
    fseek(fp1, 2L*sizeof(float)*(long)iFltVertNum,    SEEK_CUR); /* this would be local coordintates (from gridding) of the vertices..; not needed */
    /*-------------------------------------------*/
    ret =fread(fVDg_Epos, sizeof(float),iFltVertNum,fp1);
    ret =fread(fVDg_Npos, sizeof(float),iFltVertNum,fp1);
    ret =fread(fVDg_Zpos, sizeof(float),iFltVertNum,fp1);
    /*-------------------------------------------*/
    ret =fread(fTDg_CentEpos, sizeof(float),iFltPtchNum,fp1);
    ret =fread(fTDg_CentNpos, sizeof(float),iFltPtchNum,fp1);    
    ret =fread(fTDg_CentZpos, sizeof(float),iFltPtchNum,fp1);
    /*-------------------------------------------*/
    fseek(fp2, 1L*sizeof(int),     SEEK_CUR); /* contains the patch number again -> skipped */
    /*-------------------------------------------*/

    fseek(fp2, (1L*sizeof(float)*(long)(iSTARTPOS_F[iRANK])), SEEK_CUR); /*skip the part that is in front of segId for specific RANK*/
    ret =fread(fTDl_StatFric, sizeof(float),iOFFSET_F[iRANK],fp2);   /*read the SegID for specific RANK; first realization of STATIC friction coefficient */
    fseek(fp2, (1L*sizeof(float)*(long)(iFltPtchNum -iSTARTPOS_F[iRANK]-iOFFSET_F[iRANK])), SEEK_CUR); /*skip over to end of SegIDs => in total i will have moved by iFltPtchNum*/

    fseek(fp2, (1L*sizeof(float)*(long)(iSTARTPOS_F[iRANK])), SEEK_CUR); /*skip the part that is in front of segId for specific RANK*/
    ret =fread(fTDl_DynFric, sizeof(float),iOFFSET_F[iRANK],fp2);   /*read the SegID for specific RANK; first realization of DYNAMIC friction coefficient */
    fseek(fp2, (1L*sizeof(float)*(long)(iFltPtchNum -iSTARTPOS_F[iRANK]-iOFFSET_F[iRANK])), SEEK_CUR); /*skip over to end of SegIDs => in total i will have moved by iFltPtchNum*/

    fseek(fp2, (1L*sizeof(int)*(long)(iSTARTPOS_F[iRANK])), SEEK_CUR); /*skip the part that is in front of segId for specific RANK*/
    ret =fread(iTDl_StabType, sizeof(int),iOFFSET_F[iRANK],fp2);   /*read the SegID for specific RANK; if patch is seismogenic or not => is actually not used i.e., rewritten in K-matrix (INCLUDES depth-dependence due to (a-b) and also fractal strength heterogeneity and unlocked patches !! */
    fseek(fp2, (1L*sizeof(int)*(long)(iFltPtchNum -iSTARTPOS_F[iRANK]-iOFFSET_F[iRANK])), SEEK_CUR); /*skip over to end of SegIDs => in total i will have moved by iFltPtchNum*/
    
    fclose(fp1);
    fclose(fp2);
    /*-------------------------------------------*/
    if (iBndPtchNum > 0)
    {   ret =fread(&iTemp1, sizeof( int),1,fp3);         /* this is currently read PatchNum for boundary */
        ret =fread(&iTemp2, sizeof( int),1,fp3);        /* this is currently read VertexNum for boundary*/
        if ((iBndPtchNum != iTemp1) || (iBndVertNum != iTemp2))            {    fprintf(stdout,"Patch or vertex number not matching: patchum1 %d  patchum2 %d       vertexnum1 %d    vertexnum2 %d\n", iBndPtchNum, iTemp1, iBndVertNum, iTemp2);        exit(10);        }
        /*-------------------------------------------*/
        ret =fread(&iTDg_V1[iFltPtchNum], sizeof(int),iBndPtchNum,fp3); /*not sure if I need the "&" in front of the buffer,.. .the iFltPtchNum is the offset,...*/
        ret =fread(&iTDg_V2[iFltPtchNum], sizeof(int),iBndPtchNum,fp3); /*just tested it, worked on laptop (the offset stuff) -guess it will also work on HPC*/
        ret =fread(&iTDg_V3[iFltPtchNum], sizeof(int),iBndPtchNum,fp3); /*the "&" is required...*/
        /*-------------------------------------------*/
        fseek(fp3, 1L*sizeof(int)  *(long)iBndPtchNum ,   SEEK_CUR); /* this is fault section ID => but don't need it => skip over it */
        fseek(fp3, 2L*sizeof(float)*(long)iBndVertNum,    SEEK_CUR); /* this would be local coordintates (from gridding) of the vertices..; not needed */   
        /*-------------------------------------------*/
        ret =fread(&fVDg_Epos[iFltVertNum], sizeof(float),iBndVertNum,fp3);
        ret =fread(&fVDg_Npos[iFltVertNum], sizeof(float),iBndVertNum,fp3);
        ret =fread(&fVDg_Zpos[iFltVertNum], sizeof(float),iBndVertNum,fp3);
        /*-------------------------------------------*/
        ret =fread(&fTDg_CentEpos[iFltPtchNum], sizeof(float),iBndPtchNum,fp3);
        ret =fread(&fTDg_CentNpos[iFltPtchNum], sizeof(float),iBndPtchNum,fp3);    
        ret =fread(&fTDg_CentZpos[iFltPtchNum], sizeof(float),iBndPtchNum,fp3);
        /*-------------------------------------------*/    
        fclose(fp3);
        
        for (i = 0; i < iBndPtchNum; i++)
        {   iTDg_V1[iFltPtchNum+i]      += iFltVertNum; /*so that the vertices of the boundary-patches don't overlap with the fault! */
            iTDg_V2[iFltPtchNum+i]      += iFltVertNum;
            iTDg_V3[iFltPtchNum+i]      += iFltVertNum;       
    }   }
    /*---------------------------------------------------------------------------------*/
    for (i = 0; i < (iFltVertNum + iBndVertNum); i++)            
    {   fVDg_Epos[i]     *= 1000.0;             fVDg_Npos[i]     *= 1000.0;         fVDg_Zpos[i]     *= 1000.0; 
    }
    for (i = 0; i < (iFltPtchNum + iBndPtchNum); i++)            
    {          //   fprintf(stdout,"Z = %f \n",fVDg_Zpos[i]);     }
        iTDg_V1[i]        -= 1;                 iTDg_V2[i]       -= 1;              iTDg_V3[i]       -= 1;                  
        fTDg_CentEpos[i]  *= 1000.0;            fTDg_CentNpos[i] *= 1000.0;         fTDg_CentZpos[i] *= 1000.0;    
    }/* the minus 1 is to bring the index to be conform with C standard, starting at "0" */
    /*---------------------------------------------------------------------------------*/
    for (i = 0; i < iOFFSET_F[iRANK]; i++)
    {   iTDl_SegID[i]     -= 1; 
        fTDl_RefStrssRateStk[i] = cosf(fSD_SlipRake[iTDl_SegID[i]]) *fSD_StrssRate[iTDl_SegID[i]];
        fTDl_RefStrssRateDip[i] = sinf(fSD_SlipRake[iTDl_SegID[i]]) *fSD_StrssRate[iTDl_SegID[i]];
        fTDl_SlipRate[i]        = fSD_SlipRate[iTDl_SegID[i]];
        fTDl_SlipRake[i]        = fSD_SlipRake[iTDl_SegID[i]];
        fTDl_CurrFric[i]        = fTDl_StatFric[i];
        
        fTDl_RefStatFric[i]     = fTDl_StatFric[i];
        fTDl_RefDynFric[i]      = fTDl_DynFric[i];
    }
    /*---------------------------------------------------------------------------------*/
    return;
}
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
void  DefineMoreParas(const int iRANK, const int *iSTARTPOS_F, const int *iOFFSET_F, const int *iSTARTPOS_B, const int *iOFFSET_B, const  int iFltPtchNum, const  int iBndPtchNum, const  int *iTDg_V1, const  int *iTDg_V2, const  int *iTDg_V3, const  int *iTDl_SegID, const float *fTDg_CentEpos, const float *fTDg_CentNpos, const float *fTDg_CentZpos, const float *fVDg_Epos, const float *fVDg_Npos, const float *fVDg_Zpos, const float *fSD_CritSlipDist, const float *fSD_CritSlipD_vari, const float *fMD_Vp, const float *fMD_VpVsRatio, const float *fMD_MedDense, const float *fMD_AddNrmStrss, const float fMD_g, const float fdeltTincr, const float *fRandVector, const  int iRandNumber,  int *iRandPos, int *iTDlg_TravTimesP, int *iTDlg_TravTimesS, int *iMD_GlobTTmax, float *fTDl_Area, float *fTDl_RefNormStrss, float *fTDl_RefDcVal, float *fTDl_Curr_DcVal, float *fTDlg_LocSrcRcv_H, float *fTDlg_LocSrcRcv_V, float *fTDlg_LocSrcRcv_N)
{
    int     i,              j,              iVectPos,        globi;
    float   fTempP,         fTempS,         fTemp;
    float   fP1[3],         fP2[3],         fP3[3];
    float   fP1P2[3],       fP1P3[3],       fP1P2crossP1P3[3];
    float   fSrcRcvVect[3];         
    float   fvNrm[3],       fvStk[3],       fvDip[3];
        
    /*---------------------------------------------------------------------------------*/
    for (i = 0; i < iOFFSET_F[iRANK]; i++)
    {    globi  = i+iSTARTPOS_F[iRANK];
        /*---------------------------------------------------------------------------------*/
        fP1[0]   = fVDg_Epos[iTDg_V1[globi]];      fP1[1]   = fVDg_Npos[iTDg_V1[globi]];      fP1[2]   = fVDg_Zpos[iTDg_V1[globi]];
        fP2[0]   = fVDg_Epos[iTDg_V2[globi]];      fP2[1]   = fVDg_Npos[iTDg_V2[globi]];      fP2[2]   = fVDg_Zpos[iTDg_V2[globi]];
        fP3[0]   = fVDg_Epos[iTDg_V3[globi]];      fP3[1]   = fVDg_Npos[iTDg_V3[globi]];      fP3[2]   = fVDg_Zpos[iTDg_V3[globi]];
        
        fP1P2[0] = fP2[0] - fP1[0];                fP1P2[1] = fP2[1] - fP1[1];                fP1P2[2] = fP2[2] - fP1[2];
        fP1P3[0] = fP3[0] - fP1[0];                fP1P3[1] = fP3[1] - fP1[1];                fP1P3[2] = fP3[2] - fP1[2];
    
        fP1P2crossP1P3[0] = fP1P2[1]*fP1P3[2] - fP1P2[2]*fP1P3[1];
        fP1P2crossP1P3[1] = fP1P2[2]*fP1P3[0] - fP1P2[0]*fP1P3[2];
        fP1P2crossP1P3[2] = fP1P2[0]*fP1P3[1] - fP1P2[1]*fP1P3[0];
        /*---------------------------------------------------------------------------------*/
        GetLocKOS_inLoadInput(fvNrm, fvStk, fvDip, fP1, fP2, fP3);
        /*---------------------------------------------------------------------------------*/
        fTDl_RefNormStrss[i] = -1.0*(fMD_MedDense[0] *fMD_g *fabs(fTDg_CentZpos[globi])) + fMD_AddNrmStrss[0]; /* all in Pa; again -1 to convert to compressive=negative */
        /*---------------------------------------------------------------------------------*/
        fTDl_Area[i]         = 0.5*sqrtf( fP1P2crossP1P3[0]*fP1P2crossP1P3[0] +fP1P2crossP1P3[1]*fP1P2crossP1P3[1] +fP1P2crossP1P3[2]*fP1P2crossP1P3[2]);
        /*---------------------------------------------------------------------------------*/     
        fTDl_Curr_DcVal[i]   = fSD_CritSlipDist[iTDl_SegID[i]] *(1.0 + (fRandVector[iRandPos[0]]*2.0 -1.0)*fSD_CritSlipD_vari[iTDl_SegID[i]]/100.0  );  
        fTDl_RefDcVal[i]     = fTDl_Curr_DcVal[i];
        
        iRandPos[0]++;      if (iRandPos[0] >= iRandNumber) {   iRandPos[0] = rand() % iRandNumber;     } //this respawns new position of randPos somewhere in range of zero to randnumber-1     
                /*---------------------------------------------------------------------------------*/
        for (j = 0; j < iFltPtchNum;  j++)
        {   iVectPos                 = i*iFltPtchNum + j;
            fTempP                   = sqrtf( (fTDg_CentEpos[globi]-fTDg_CentEpos[j])*(fTDg_CentEpos[globi]-fTDg_CentEpos[j]) + (fTDg_CentNpos[globi]-fTDg_CentNpos[j])*(fTDg_CentNpos[globi]-fTDg_CentNpos[j]) + (fTDg_CentZpos[globi]-fTDg_CentZpos[j])*(fTDg_CentZpos[globi]-fTDg_CentZpos[j]) )/fMD_Vp[0]; /* the distance between both */
            fTempS                   = fTempP*fMD_VpVsRatio[0];
            /*---------------------------------------------------------------------------------*/
            iTDlg_TravTimesP[iVectPos] = (int)(fTempP/fdeltTincr);
            iTDlg_TravTimesS[iVectPos] = (int)(fTempS/fdeltTincr);
            /*---------------------------------------------------------------------------------*/
            iMD_GlobTTmax[0]         = (iMD_GlobTTmax[0] > iTDlg_TravTimesS[iVectPos]) ? iMD_GlobTTmax[0] : iTDlg_TravTimesS[iVectPos];     
            /*---------------------------------------------------------------------------------*/
            fSrcRcvVect[0]           = fTDg_CentEpos[j] - fTDg_CentEpos[globi]; /*this is vector from current patch (source) to receiver*/
            fSrcRcvVect[1]           = fTDg_CentNpos[j] - fTDg_CentNpos[globi]; 
            fSrcRcvVect[2]           = fTDg_CentZpos[j] - fTDg_CentZpos[globi];
            
            if (j == globi)
            {   fTDlg_LocSrcRcv_N[iVectPos]    = 0.0;       fTDlg_LocSrcRcv_H[iVectPos]    = 0.0;        fTDlg_LocSrcRcv_V[iVectPos]    = 0.0;
            }
            else
            {   fTemp                    = sqrtf(fSrcRcvVect[0]*fSrcRcvVect[0] +fSrcRcvVect[1]*fSrcRcvVect[1] +fSrcRcvVect[2]*fSrcRcvVect[2]);
                fSrcRcvVect[0]          /= fTemp;                fSrcRcvVect[1] /= fTemp;                fSrcRcvVect[2] /= fTemp; 
                /*---------------------------------------------------------------------------------*/                    
                fTDlg_LocSrcRcv_N[iVectPos]    = fvNrm[0]*fSrcRcvVect[0] + fvNrm[1]*fSrcRcvVect[1] + fvNrm[2]*fSrcRcvVect[2]; /*this rotates the vector from current source to receiver*/
                fTDlg_LocSrcRcv_H[iVectPos]    = fvStk[0]*fSrcRcvVect[0] + fvStk[1]*fSrcRcvVect[1] + fvStk[2]*fSrcRcvVect[2]; /*into local coordinate system (that is the idea...)*/
                fTDlg_LocSrcRcv_V[iVectPos]    = fvDip[0]*fSrcRcvVect[0] + fvDip[1]*fSrcRcvVect[1] + fvDip[2]*fSrcRcvVect[2];                
    }   }   }
    /*---------------------------------------------------------------------------------*/
    return;
}
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*/
void GetLocKOS_inLoadInput(float fvNrm[3], float fvStk[3], float fvDip[3], const float fP1[3], const float fP2[3], const float fP3[3])
{    /* this comes basically from Medhii's code... */
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