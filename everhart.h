#include <stdio.h>
//#include <conio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#define RADANVMAX 21

void FORCEE(double *Xc, double *Vc, double Tc, double *Fc);

int rada27e(double *X, double *V, double TI, double TF, double XL,
        int LL, int NV, int NI, int *NFout, int *NSout, int NCLASSin, int NOR,
        int FOLLOW)
{
/*ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ
³  E.EVERHART, PHISICS AND ASTRONOMY DEPARTMENT, UNIVERSITY OF DENVER 74
³  DENVER,COLORAD 80210. PHONE (303)-753-2238 OR 753-2362
³  PROGRAM ADAPTED FOR BESM-6 BY S.TARASEVICH, DEPT. OF ALGORITHMIZATION
³  INSTITUTE OF THEORETICAL ASTRONOMY, LENINGRAD 192187 USSR.  1975
³  INTEGRATOR FOR ORDERS 7,11,15,19,23,27.  NOR IS THE ORDER.
³  DOUBLE PRECISION VERSION FOR BESM-6  OF INTEGRATOR RADAU
³  NV IS THE NUMBER OF DEPENDENT VARIABLES
³  NCLASS IS 1 FOR 1ST-ORDER DIFF EQ,AND 2 FOR 2ND ORDER DIFF EQ.
³  IF FIRST DERIVATIVES ARE NOT PRESENT (CLASS IIS),THEN USE NCLASS=-2.
³  X(NV) IS THE INITIAL POSITION VECTOR.  IT RETURNS AS THE FINAL VALUE
³  V(NV) IS THE INITIAL VELOCITY VECTOP. IT RETURNS AS THE FINAL VALUE
³  IN THE CASE NCLASS IS UNITY,THEN V IS SIMPLY ZEROED
³  TI IS  INITIAL TIME, TF IS FINAL TIME, NF IS NUMBER OF FUNCTION EVALU
³  NS IS THE NUMBER OF SEQUENCES.
³  PROGRAM SET UP FOR A MAXIMUM OF 18 SIMULTANEOUS EQUATIONS.
³  LL CONTROLS SEQLENCE SIZE. THUS SS=10.**(-LL) IS DESIRED SIZE OF A TE
³  AS IN CONTROL SYSTEM I.
³  IF LL.LT.0, THEN XL IS THE SPECIFIED CONSTANT SEQUENCE SIZE
³  WILL INTEGRATE IN A NEGATIVE DIRECTION IF TF.LT.TI
*ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ*/
 int NF,NS,NCLASS;
 int KD,KD2,KE,KF,LA,N,K,L,J,M,JD,LC,LB,LD,LE,JDM,NL,NCOUNT,JR;

 static double 
        C[79],D[79],R[79],XI[79],
        H[15],W[14],U[14],
        ZERO=0.0,
        ONE=1.0,
        TM,T,S,Q,TEMP,RES,T2,W1,TP,SR;
 double VAL,TVAL,HSUM,SM,BDUBL,PW,TDIF,DIR,SS;
 static double 
               F1[RADANVMAX+1],FJ[RADANVMAX+1],Y[RADANVMAX+1],Z[RADANVMAX+1],
               BE[14][RADANVMAX+1],
               BT[14][RADANVMAX+1],B[14][RADANVMAX+1];
 int J2,NPQ,NSF,NPER,NCL,NES;
 int  NW[]={0,0,0,1,3,6,10,15,21,28,36,45,55,66,78};
 int  MC[]={0,1,13,24,34,43,51,58,64,69,73,76,78};
 int  NXI[]={0,2,3,4,5,6,7,8,9,10,11,12,13,3,6,10,15,21,28,36,45,55,66,
      78,4,10,20,35,56,84,120,165,220,286,5,15,35,70,126,210,330,495,
      715,6,21,56,126,252,462,792,1287,7,28,84,210,462,924,1716,8,36,
      120,330,792,1716,9,45,165,495,1287,10,55,220,715,11,66,286,12,
      78,13};
 double HH[]={0.0,
     .212340538239152 , .590533135559265, .911412040487296,
     .0985350857988264,.3045357266463639,.5620251897526139,
     .8019865821263918,.9601901429485313,.0562625605369221,
     .1802406917368924,.3526247171131696,.5471536263305554,
     .7342101772154105,.8853209468390958,.9775206135612875,
     .0362578128832095,.1180789787899987,.2371769848149604,
     .3818827653047060,.5380295989189891,.6903324200723622,
     .8238833438370047,.9256126102908040,.9855875903511235,
     .0252736203975203,.0830416134474051,.1691751003771814,
     .2777967151090321,.4015027202328608,.5318623869104160,
     .6599918420853348,.7771593929561621,.8753807748555569,
     .9479645488728194,.9899817195383196,.0186103650109879,
     .0614755408992690,.1263051786933106,.2098429717265625,
     .3078989982803983,.4155560359786595,.5274156139958823,
     .6378686027177612,.7413764592942375,.8327489886084423,
     .9074047753009974,.9616018612603216,.9926353489739107};

  NCLASS=NCLASSin;
      if(FOLLOW)
        {
          NPER=0;
          goto L2222;
        }
L2221:KD=(NOR-3)/2;
      KD2=KD/2;
      KE=KD+1;
      KF=KD+2;
      PW=ONE/(KD+3);
      NPER=0;
      NSF=0;
      NCL=NCLASS==1;
      NPQ=NCLASS< 2;
      SR=1.5;
      if(NV==1) SR=1.2;
      NES=LL<0;
      TDIF=TF-TI;
      DIR=TDIF/fabs(TDIF);
      if(NES) XL=fabs(XL)*DIR;
      NCLASS=abs(NCLASS);
      LA=KD2*KD2-1;

      for(N=2;N<=KF;N++)
      {
       LA=LA+1;
       H[N]=HH[LA];
       W[N-1]=ONE/(N+N*N*(NCLASS-1));
       U[N-1]=N+1;
      }
      for(K=1;K<=NV;K++)
      {
       if(NCL) V[K]=ZERO;
       for(L=1;L<=KE;L++)
       {
        BT[L][K]=ZERO;
        B[L][K]=ZERO;
       }
      }
      W1=ONE/NCLASS;
      for(J=1;J<=KD;J++)
      {
       M=MC[J];
       JD=J+1;
       for(L=JD;L<=KE;L++)
       {
        XI[M]=NXI[M]*W[J]/W[L];
        M=M+1;
       }
      }
      C[1]=-H[2]*W[1];
      D[1]=H[2]/W[2];
      R[1]=ONE/(H[3]-H[2]);
      LA=1;
      LC=1;
      for(K=3;K<=KE;K++)
      {
       LB=LA;
       LA=LC+1;
       LC=NW[K+1];
       JD=LC-LA;
       C[LA]=-H[K]*C[LB];
       C[LC]=(C[LA-1]/W[JD]-H[K])*W[JD+1];
       D[LA]=H[2]*D[LB]*W[K-1]/W[K];
       D[LC]=(D[LA-1]*W[K-1]+H[K])/W[K];
       R[LA]=ONE/(H[K+1]-H[2]);
       R[LC]=ONE/(H[K+1]-H[K]);
       if(K==3) goto L73;
       for(L=4;L<=K;L++)
       {
        LD=LA+L-3;
        LE=LB+L-4;
        JDM=LD-LA;
        C[LD]=W[JDM+1]*C[LE]/W[JDM]-H[K]*C[LE+1];
        D[LD]=(D[LE]+H[L-1]*D[LE+1])*W[K-1]/W[K];
        R[LD]=ONE/(H[K+1]-H[L-1]);
       }
       L73: continue;
      }
      SS=pow(10.0,(double)(-LL));
      NL=NI+30;
/*  SET IN A RESONABLE ESTIMATE TO T BASED ON EXPERIENCE.  (SAME SIGN AS */
      TP=( (NOR/11.0) * pow(0.5,0.4*LL) )*DIR*0.5;
      if(NES) TP=XL;
      if(TP/TDIF >  0.5) TP=0.5*TDIF;
      NF=0;
      NCOUNT=0;
L4000:NS=0;
      TM=TI;
      SM=10000.0;

      FORCEE(X,V,TM,F1);

      NF=NF+1;

/*  LOOP 58 FINDS THE BETA-VALUES FROM THE CORRECTED B-VALUES, USING D-C */

L722: for(K=1;K<=NV;K++)
      {
       BE[KE][K]=B[KE][K]/W[KE];
       for(J=1;J<=KD;J++)
       {
        JD=J+1;
        BE[J][K]=B[J][K]/W[J];
        for(L=JD;L<=KE;L++)
        {
          N=NW[L]+J;
          BE[J][K]=BE[J][K]+D[N]*B[L][K];
        }
       }
      }
      T=TP;
      TVAL=fabs(T);
      T2=pow(T,(double)(NCLASS));

     for(M=1;M<=NL;M++)
     {
      J2=1;
      for(J=2;J<=KF;J++)
      {
      JD=J-1;
      LA=NW[JD];
      JDM=J-2;
      S=H[J];
      Q=pow(S,(double)(NCLASS-1));
      if(NPQ) goto L5100;
      for(K=1;K<=NV;K++)
      {
       RES=B[KE][K];
       TEMP=RES*U[KE];
       for(L=1;L<=KD;L++)
       {
        JR=KE-L;
        RES=B[JR][K]+S*RES;
        TEMP=B[JR][K]*U[JR]+S*TEMP;
       }
       Y[K]=X[K]+Q*(T*V[K]+T2*S*(F1[K]*W1+S*RES));
       Z[K]=V[K]+S*T*(F1[K]+S*TEMP);
      }
      goto L5200;
L5100:for(K=1;K<=NV;K++)
      {
       RES=B[KE][K];
       for(L=1;L<=KD;L++)
       {
        JR=KE-L;
        RES=B[JR][K]+S*RES;
       }
       Y[K]=X[K]+Q*(T*V[K]+T2*S*(F1[K]*W1+S*RES));
      }
L5200:
      FORCEE(Y,Z,TM+S*T,FJ);
      NF=NF+1;
      if(J2) goto L702;
      for(K=1;K<=NV;K++)
      {
      TEMP=BE[JD][K];
      RES=(FJ[K]-F1[K])/S;
      N=LA;
      for(L=1;L<=JDM;L++)
      {
       N=N+1;
       RES=(RES-BE[L][K])*R[N];
      }
      BE[JD][K]=RES;
      TEMP=RES-TEMP;
      B[JD][K]=B[JD][K]+TEMP*W[JD];
      N=LA;
      for(L=1;L<=JDM;L++)
      {
       N=N+1;
       B[L][K]=B[L][K]+C[N]*TEMP;
      }
      }
      goto L174;
L702: J2=0;

      for(K=1;K<=NV;K++)
      {
       TEMP=BE[1][K];
       RES=(FJ[K]-F1[K])/S;
       BE[1][K]=RES;
       B[1][K]=B[1][K]+(RES-TEMP)*W[1];
      }
L174: continue;
      }
      if(M < NI) goto L175;
      HSUM=ZERO;
      VAL=pow(TVAL,(double)(-KE));
      for(K=1;K<=NV;K++)
      {
       BDUBL=B[KE][K];
       HSUM=HSUM+BDUBL*BDUBL;
      }
      HSUM=VAL*sqrt(HSUM);
      if(NSF) goto L175;
      if(fabs((HSUM-SM)/HSUM) < 0.01) break;   /*  GOTO L176; */
      SM=HSUM;
L175: continue;
     }

/* C  THIS NEXT PART THE PROPER STARTING VALUE FOR T */
L176: if(NSF) goto L180;
      TP=pow(SS/HSUM,(double)(PW))*DIR;
      if(NES) TP=XL;
      if(NES) goto L170;
      if(TP/T > ONE) goto L170;
      TP=0.8*TP;
      NCOUNT=NCOUNT+1;
      if(NCOUNT > 10)
        {
          *NSout=NS;
          *NFout=NF;
                      goto Lreturn;
        }
      goto L4000;
L170: NSF=1;
/*  FIND POSITION (AND VELOCITY FOR CLASS II AND IIS) AT THE END OF THE S */

L180: for(K=1;K<=NV;K++)
      {
      RES=B[KE][K];
      for(L=1;L<=KD;L++) RES=RES+B[L][K];
      X[K]=X[K]+V[K]*T+T2*(F1[K]*W1+RES);
      if(NCL) goto L35;
      RES=B[KE][K]*U[KE];
      for(L=1;L<=KD;L++) RES=RES+B[L][K]*U[L];
      V[K]=V[K]+T*(F1[K]+RES);
 L35: continue;
      }
      TM=TM+T;
      NS=NS+1;
      if(!NPER) goto L2222;
      *NSout=NS;
      *NFout=NF;
                      goto Lreturn;
                      
/*  CONTROL ON SIZE OF NEXT SEQUENCE AND RETURN WHEN TF IS REACHED
 */
L2222:FORCEE(X,V,TM,F1);
      NF=NF+1;
      if(NES) goto L341;
      TP=(pow(SS/HSUM,(double)(PW)))*DIR;
      if(TP/T > SR) TP=SR*T;
L341: if(NES) TP=XL;
      if(DIR*(TF-TM) > DIR*TP*1.001) goto L77;
      TP=TF-TM;
      NPER=1;
/*  PREDICT B-VALUES FOR NEXT SEQUENCE, */
 L77: Q=TP/T;
      for(K=1;K<=NV;K++)
      {
       RES=ONE;
      for(J=1;J<=KE;J++)
       {
      if(NS > 1) BT[J][K]=B[J][K]-BT[J][K];
      if(J == KE) goto L740;
      M=MC[J];
      JD=J+1;
      for(L=JD;L<=KE;L++)
      {
       B[J][K]=B[J][K]+XI[M]*B[L][K];
       M=M+1;
      }
L740: RES=RES*Q;
      TEMP=RES*B[J][K];
      B[J][K]=TEMP+BT[J][K];
      BT[J][K]=TEMP;
       }
      }
      NL=NI;
      goto L722;

Lreturn:
              return 0;
}