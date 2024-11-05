#include<stdio.h>
#include<stdlib.h>
#include<string.h>


// /0:C0/1:C1/2:Not/3:Copy/4:And2/5:Or2/6:Xor2/7:Ada2/8:DUMMY/
char fo1[] = 
  { 2,-1, -1,-1,  2,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1,  2,-1,  2,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
    1,-1,  0,-1,  2,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
    0, 2, -1,-1,  0, 2, -1,-1,  1, 2,  1, 2,  0, 2,  1, 2,  2, 2,
    0, 2,  1, 0,  0, 2, -1,-1,  1, 1,  1, 1,  0, 2,  1, 2,  2, 2,
    0, 2, -1,-1,  0, 0,  0, 1,  1, 2,  1, 2,  0, 2,  1, 2,  2, 2,
    0, 0,  1, 1,  0, 0,  0, 1,  1, 0,  0, 1,  0, 2,  1, 2,  2, 2,
    0, 0,  0, 1,  0, 2,  1, 1, -1,-1,  1, 1,  0, 0,  0, 1,  2, 2,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1};

char fo2[] =
  {-1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1,  1, 0, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1,  0, 1, -1,-1, -1,-1, -1,-1,
   -1,-1, -1,-1,  1, 1, -1,-1, -1,-1,  1, 0, -1,-1, -1,-1, -1,-1,
   -1,-1,  1, 0,  1, 0, -1,-1, -1,-1, -1,-1,  1, 1,  1, 0, -1,-1,
   -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1, -1,-1};

#define	MX_FB0N	(700)

#define	MX_N	(600)
#define	MX_NB	(3*MX_N+2)
#define	MX_NG	(MX_N*MX_N*8+MX_NB)

char *s0,*sbs,*sn0,*sbs1;
int *gts,*b1s,*b2s;
int *pt0,*pt0i,*pte,*ptei;

char *fns,*pts,*ptb;
int *ptf1,*ptf1i,*ptf2,*ptf2i;
int *ptb1,*ptb1i,*ptb2,*ptb2i;

//------------------------
int mallocs(){
  if((gts =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:gts\n");exit(0);}
  if((b1s =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:b1s\n");exit(0);}
  if((b2s =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:b2s\n");exit(0);}
  if((sbs =(char*)malloc((long long)sizeof(char)*MX_NB*MX_NG))==NULL){
	printf("E:sbs\n");exit(0);}
  if((fns =(char*)malloc((long)sizeof(char)*MX_NG*18*2))==NULL){printf("E:fns\n");exit(0);}
  if((pts =(char*)malloc((long)sizeof(char)*MX_NG*18*2))==NULL){printf("E:pts\n");exit(0);}
  if((ptf1=(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptf1i\n");exit(0);}
  if((ptf1i=(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptf1i\n");exit(0);}
  if((ptf2 =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptf2\n");exit(0);}
  if((ptf2i=(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptf2i\n");exit(0);}
  if((ptb1 =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptb1\n");exit(0);}
  if((ptb1i=(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptb1i\n");exit(0);}
  if((ptb2 =(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptb2\n");exit(0);}
  if((ptb2i=(int*)malloc((long)sizeof(int)*MX_NG))==NULL){printf("E:ptb2i\n");exit(0);}

  if((s0=(char*)malloc((long)sizeof(char)*MX_NB))==NULL){printf("E:s0\n");exit(0);}
  if((pt0=(int*)malloc((long)sizeof(int)*MX_NB))==NULL){printf("E:pt0\n");exit(0);}
  if((pt0i=(int*)malloc((long)sizeof(int)*MX_NB))==NULL){printf("E:pt0i\n");exit(0);}
  if((pte=(int*)malloc((long)sizeof(int)*MX_NB))==NULL){printf("E:pte\n");exit(0);}
  if((ptei=(int*)malloc((long)sizeof(int)*MX_NB))==NULL){printf("E:ptei\n");exit(0);}
  if((sn0=(char*)malloc((long)sizeof(char)*MX_NB))==NULL){printf("E:sn0\n");exit(0);}

  if((sbs1 =(char*)malloc((long long)sizeof(char)*MX_NB*MX_NG))==NULL){
	printf("E:sbs1\n");exit(0);}
  if((ptb =(char*)malloc((long)sizeof(char)*MX_NG*18*2))==NULL){printf("E:ptb\n");exit(0);}
}


// nxn2
// /0:C0/1:C1/2:Not/3:Copy/4:And2/5:Or2/6:Xor2/7:Ada2/
// ex. C0: gts=gts+[0]; b1s=b1s+[i1]; b2s=b2s+[-1]
// ex. And2: gts=gts+[4]; b1s=b1s+[i1]; b2s=b2s+[i2]
// (i1,i2,i3,i4,i5) -> (c,s,i3,i4),  c,s=i1+i2+(i3&i4), i5:Buffer
//
int appendAda42(int*gts,int*b1s,int*b2s,
   int p0,int p1,int p2,int p3,int p4){
  gts[0] = 3; b1s[0] = p4; b2s[0] = p2;
  gts[1] = 5; b1s[1] = p2; b2s[1] = p0;
  gts[2] = 7; b1s[2] = p1; b2s[2] = p2;
  gts[3] = 7; b1s[3] = p1; b2s[3] = p0;
  gts[4] = 4; b1s[4] = p3; b2s[4] = p2;
  gts[5] = 3; b1s[5] = p2; b2s[5] = p4;
}

// s0 = [0 for ii in range(nb1+2)] + i2v(nv,nb1+nb2)
// a1 = v2i(vs[0],0,nb1)
// a2 = v2i(vs[0],nb1*2+2,nb1*2+2+nb2)
//
// nb = 2*nb1+nb2+2
// ng = nb1+(6*nb1+2)*nb2
//
int nxn2(int*gts,int*b1s,int*b2s,int nb1,int nb2){
  int i1,i2,jj;

  jj=0;
  for(i1=nb2-1;i1>-1;i1--){
    gts[jj]=3; b1s[jj]=nb1+1; b2s[jj]=nb1*2+2+i1; jj++;
    for(i2=nb1-1;i2>-1;i2--){
      appendAda42(gts+jj,b1s+jj,b2s+jj,
        nb1+1, nb1+2+i1+i2, i2, nb1*2+2+i1, nb1);
      jj+=6;
    }
    gts[jj]=0; b1s[jj]=nb1+1; b2s[jj]=-1; jj++;
  }
  for(i1=0;i1<nb1;i1++){
    gts[jj]=0; b1s[jj]=nb1+2+i1; b2s[jj]=-1; jj++;
  }
  return(0);
}

//////

int trcSn0(int nb,int ng,char*fns,char*sn,int*pte,int*ptei){
  int ii,jj; long ngn;
  char *f;
  if(sn[0] < 0) return(-1);
  ngn = ng*18;
  for(ii=0;ii<nb;ii++){
    if(sn[ii] < 3){
      f = fns+pte[ii]*18;
      for(jj=0;jj<18;jj+=2)
        if(f[jj+ptei[ii]] != sn[ii]){ f[jj]=-1; f[jj+1]=-1; }
      f = f+ngn;
      for(jj=0;jj<18;jj+=2)
        if(f[jj+ptei[ii]] != sn[ii]){ f[jj]=-1; f[jj+1]=-1; }
  }  }
  return(0);
}

int trcFB0(int nb,int ng,char*fns,char*s0,char*sn,int*b1s,int*b2s,
  int*ptb1,int*ptb1i,int*ptb2,int*ptb2i,int n0,int ni){

  int ii,c[] = {1,2,3},trcF0(),trcB0(),countatbl();

  for(ii=0;ii<ni;ii++){
    if(trcF0(nb,ng,fns,s0,b1s,b2s,ptb1,ptb1i,ptb2,ptb2i,n0)==0) return(0);
    c[2]=c[1];c[1]=c[0];c[0]=countatbl(ng,fns,NULL);
    if(trcB0(nb,ng,fns,b1s,b2s,ptb1,ptb1i,ptb2,ptb2i,n0)==0) return(0);
    c[2]=c[1];c[1]=c[0];c[0]=countatbl(ng,fns,NULL);
    if(c[0]==c[1] && c[1]==c[2]) return(ii);
    /*if not fnschk0(fns): 
      if dsp: print("ZERO fns");
      return 0 */
  }
//  printf("OVER\n");
  return(-1);
}

int trcF0(int nb,int ng,char*fns,char*s0,int*b1s,int*b2s,
  int*ptb1,int*ptb1i,int*ptb2,int*ptb2i,int n0){ // n0=0

  int n,ngn,ii;
  char *f1,*f2,*fb1,*fb2;
  int ppi;
  
  ngn = ng*18;
  for(n=n0;n<ng;n++){
    f1 = fns+n*18; f2 = f1+ngn;
    if(ptb1[n] < n0){
      if(b1s[n] > -1){  //s0
        if(s0[b1s[n]]==0){
          f1[2]=-1;f1[8] =-1;f1[14]=-1; f2[2]=-1;f2[8] =-1;f2[14]=-1;
          f1[3]=-1;f1[9] =-1;f1[15]=-1; f2[3]=-1;f2[9] =-1;f2[15]=-1;
          f1[4]=-1;f1[10]=-1;f1[16]=-1; f2[4]=-1;f2[10]=-1;f2[16]=-1;
          f1[5]=-1;f1[11]=-1;f1[17]=-1; f2[5]=-1;f2[11]=-1;f2[17]=-1;
        }else if(s0[b1s[n]]==1){
          f1[0]=-1;f1[6] =-1;f1[12]=-1; f2[0]=-1;f2[6] =-1;f2[12]=-1;
          f1[1]=-1;f1[7] =-1;f1[13]=-1; f2[1]=-1;f2[7] =-1;f2[13]=-1;
          f1[4]=-1;f1[10]=-1;f1[16]=-1; f2[4]=-1;f2[10]=-1;f2[16]=-1;
          f1[5]=-1;f1[11]=-1;f1[17]=-1; f2[5]=-1;f2[11]=-1;f2[17]=-1;
        }else if(s0[b1s[n]]==2){
          f1[2]=-1;f1[8] =-1;f1[14]=-1; f2[2]=-1;f2[8] =-1;f2[14]=-1;
          f1[3]=-1;f1[9] =-1;f1[15]=-1; f2[3]=-1;f2[9] =-1;f2[15]=-1;
          f1[0]=-1;f1[6] =-1;f1[12]=-1; f2[0]=-1;f2[6] =-1;f2[12]=-1;
          f1[1]=-1;f1[7] =-1;f1[13]=-1; f2[1]=-1;f2[7] =-1;f2[13]=-1;
      } }
    }else{
      fb1 = fns+ptb1[n]*18; fb2 = fb1+ngn; ppi = ptb1i[n];
      for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==0 || fb2[ii+ppi]==0) break;
      if(ii==18){
        f1[0]=-1;f1[6] =-1;f1[12]=-1; f2[0]=-1;f2[6] =-1;f2[12]=-1;
        f1[1]=-1;f1[7] =-1;f1[13]=-1; f2[1]=-1;f2[7] =-1;f2[13]=-1;
      }
      for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==1 || fb2[ii+ppi]==1) break;
      if(ii==18){
        f1[2]=-1;f1[8] =-1;f1[14]=-1; f2[2]=-1;f2[8] =-1;f2[14]=-1;
        f1[3]=-1;f1[9] =-1;f1[15]=-1; f2[3]=-1;f2[9] =-1;f2[15]=-1;
      }
      for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==2 || fb2[ii+ppi]==2) break;
      if(ii==18){
        f1[4]=-1;f1[10]=-1;f1[16]=-1; f2[4]=-1;f2[10]=-1;f2[16]=-1;
        f1[5]=-1;f1[11]=-1;f1[17]=-1; f2[5]=-1;f2[11]=-1;f2[17]=-1;
    } }

    if(ptb2[n] < n0){
      if(b2s[n] > -1){  //s0
        if(s0[b2s[n]]==0){
          f1[6] =-1;f1[8] =-1;f1[10]=-1; f2[6] =-1;f2[8] =-1;f2[10]=-1;
          f1[7] =-1;f1[9] =-1;f1[11]=-1; f2[7] =-1;f2[9] =-1;f2[11]=-1;
          f1[12]=-1;f1[14]=-1;f1[16]=-1; f2[12]=-1;f2[14]=-1;f2[16]=-1;
          f1[13]=-1;f1[15]=-1;f1[17]=-1; f2[13]=-1;f2[15]=-1;f2[17]=-1;
        }else if(s0[b2s[n]]==1){
          f1[0] =-1;f1[2] =-1;f1[4] =-1; f2[0] =-1;f2[2] =-1;f2[4] =-1;
          f1[1] =-1;f1[3] =-1;f1[5] =-1; f2[1] =-1;f2[3] =-1;f2[5] =-1;
          f1[12]=-1;f1[14]=-1;f1[16]=-1; f2[12]=-1;f2[14]=-1;f2[16]=-1;
          f1[13]=-1;f1[15]=-1;f1[17]=-1; f2[13]=-1;f2[15]=-1;f2[17]=-1;
        }else if(s0[b2s[n]]==2){
          f1[6] =-1;f1[8] =-1;f1[10]=-1; f2[6] =-1;f2[8] =-1;f2[10]=-1;
          f1[7] =-1;f1[9] =-1;f1[11]=-1; f2[7] =-1;f2[9] =-1;f2[11]=-1;
          f1[0] =-1;f1[2] =-1;f1[4] =-1; f2[0] =-1;f2[2] =-1;f2[4] =-1;
          f1[1] =-1;f1[3] =-1;f1[5] =-1; f2[1] =-1;f2[3] =-1;f2[5] =-1;
      } }
    }else{
      if(b2s[n] > -1){ //単素子でない
        fb1 = fns+ptb2[n]*18; fb2 = fb1+ngn; ppi = ptb2i[n];
        for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==0 || fb2[ii+ppi]==0) break;
        if(ii==18){
          f1[0] =-1;f1[2] =-1;f1[4] =-1; f2[0] =-1;f2[2] =-1;f2[4] =-1;
          f1[1] =-1;f1[3] =-1;f1[5] =-1; f2[1] =-1;f2[3] =-1;f2[5] =-1;
        }
        for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==1 || fb2[ii+ppi]==1) break;
        if(ii==18){
          f1[6] =-1;f1[8] =-1;f1[10]=-1; f2[6] =-1;f2[8] =-1;f2[10]=-1;
          f1[7] =-1;f1[9] =-1;f1[11]=-1; f2[7] =-1;f2[9] =-1;f2[11]=-1;
        }
        for(ii=0;ii<18;ii+=2) if(fb1[ii+ppi]==2 || fb2[ii+ppi]==2) break;
        if(ii==18){
          f1[12]=-1;f1[14]=-1;f1[16]=-1; f2[12]=-1;f2[14]=-1;f2[16]=-1;
          f1[13]=-1;f1[15]=-1;f1[17]=-1; f2[13]=-1;f2[15]=-1;f2[17]=-1;
      } }
    }
    for(ii=0;ii<18;ii+=2) if(f1[ii]>-1 || f2[ii]>-1) break;  
    if(ii==18) return(0);
  }

  return(1);
}


int trcB0(int nb,int ng,char*fns,int*b1s,int*b2s,
  int*ptb1,int*ptb1i,int*ptb2,int*ptb2i,int n0){ //n0=-1

  int ii,jj,ngn;
  int pt1,pt2,pti1,pti2;
  char *f0,*f1,*f0b,*f1b;

  ngn = ng*18;
  for(ii=ng-1;ii>n0;ii--){
    pt1=ptb1[ii]; pti1=ptb1i[ii]; pt2=ptb2[ii]; pti2=ptb2i[ii];
    f0 = fns+ii*18; f1 = f0+ngn;
    if(pt1 > n0){
      f0b = fns+pt1*18; f1b = f0b+ngn;

      if(f0[0]<0 && f0[6]<0 && f0[12]<0 && f1[0]<0 && f1[6]<0 && f1[12]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti1]==0){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti1]==0){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      if(f0[2]<0 && f0[8]<0 && f0[14]<0 && f1[2]<0 && f1[8]<0 && f1[14]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti1]==1){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti1]==1){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      if(f0[4]<0 && f0[10]<0 && f0[16]<0 && f1[4]<0 && f1[10]<0 && f1[16]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti1]==2){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti1]==2){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      for(jj=0;jj<18;jj+=2) if(f0b[jj]>-1 || f1b[jj]>-1) break;
      if(jj==18) return(0);
    }

    if(pt2 > n0){
      f0b = fns+pt2*18; f1b = f0b+ngn;

      if(f0[0]<0 && f0[2]<0 && f0[4]<0 && f1[0]<0 && f1[2]<0 && f1[4]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti2]==0){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti2]==0){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      if(f0[6]<0 && f0[8]<0 && f0[10]<0 && f1[6]<0 && f1[8]<0 && f1[10]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti2]==1){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti2]==1){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      if(f0[12]<0 && f0[14]<0 && f0[16]<0 && f1[12]<0 && f1[14]<0 && f1[16]<0)
        for(jj=0;jj<18;jj+=2){
          if(f0b[jj+pti2]==2){ f0b[jj]=-1; f0b[jj+1]=-1;}
          if(f1b[jj+pti2]==2){ f1b[jj]=-1; f1b[jj+1]=-1;}
        }
      for(jj=0;jj<18;jj+=2) if(f0b[jj]>-1 || f1b[jj]>-1) break;
      if(jj==18) return(0);
    }
  }
  return(1);
}


int nextpts(int nb,int ng,int*b1s,int*b2s, int*pt0,int*pt0i,
  int*ptf1,int*ptf1i,int*ptf2,int*ptf2i,
  int*ptb1,int*ptb1i,int*ptb2,int*ptb2i,int*pte,int*ptei){

  int ii,jj,p0,p0i, f1,f1i,f2,f2i, b1,b1i,b2,b2i;
 
  // pt0
  for(ii=0;ii<nb;ii++){
    p0=-1; p0i=-1;
    for(jj=0;jj<ng;jj++){
      if(b1s[jj]==ii){ p0=jj; p0i=0; break;}
      if(b2s[jj]==ii){ p0=jj; p0i=1; break;}
    }
    pt0[ii]=p0; pt0i[ii]=p0i;
  }

  // ptf
  for(ii=0;ii<ng-1;ii++){
    f1=-1; f1i=-1;
    for(jj=ii+1;jj<ng;jj++){
      if(b1s[jj]==b1s[ii]){ f1=jj; f1i=0; break;}
      if(b2s[jj]==b1s[ii]){ f1=jj; f1i=1; break;}
    }
    ptf1[ii]=f1; ptf1i[ii]=f1i;
    f2=-1; f2i=-1;
    for(jj=ii+1;jj<ng;jj++){
      if(b1s[jj]==b2s[ii]){ f2=jj; f2i=0; break;}
      if(b2s[jj]==b2s[ii]){ f2=jj; f2i=1; break;}
    }
    ptf2[ii]=f2; ptf2i[ii]=f2i;
  }
  ptf1[ng-1]=-1; ptf1i[ng-1]=2; ptf2[ng-1]=-1; ptf2i[ng-1]=2;

  // ptb
  ptb1[0]=-1; ptb1i[0]=2; ptb2[0]=-1; ptb2i[0]=2;
  for(ii=1;ii<ng;ii++){
    b1=-1; b1i=-1;
    for(jj=ii-1;jj>-1;jj--){
      if(b1s[jj]==b1s[ii]){ b1=jj; b1i=0; break;}
      if(b2s[jj]==b1s[ii]){ b1=jj; b1i=1; break;}
    }
    ptb1[ii]=b1; ptb1i[ii]=b1i;
    b2=-1; b2i=-1;
    for(jj=ii-1;jj>-1;jj--){
      if(b1s[jj]==b2s[ii]){ b2=jj; b2i=0; break;}
      if(b2s[jj]==b2s[ii]){ b2=jj; b2i=1; break;}
    }
    ptb2[ii]=b2; ptb2i[ii]=b2i;
  }
  // pte
  for(ii=0;ii<nb;ii++){
    p0=-1; p0i=-1;
    for(jj=ng-1;jj>=0;jj--){
      if(b1s[jj]==ii){ p0=jj; p0i=0; break;}
      if(b2s[jj]==ii){ p0=jj; p0i=1; break;}
    }
    pte[ii]=p0; ptei[ii]=p0i;
  }

}

int gts2fns(char*fns,int* gts,int ng){
  int ii,ngi;
  ngi = ng*18;
  for(ii=0;ii<ng;ii++){
    memcpy(fns+18*ii,fo1+18*gts[ii],sizeof(char)*18);
    memcpy(fns+ngi+18*ii,fo2+18*gts[ii],sizeof(char)*18);
  }
}

int vs2vi(int* vo,char* vs){
  int ii;
  for(ii=0;ii<strlen(vs);ii++) vo[ii] = vs[ii]-48;
  vo[ii] = -1;
  return(strlen(vs));
}
char* vi2vs(char* vs,char* vi,int n){
  int ii;
  for(ii=0;ii<n;ii++) vs[ii] = vi[ii]+48;
  vs[ii] = '\0';
  return(vs);
}

int finput(char*fname,int*nb,int*ng,int*gts,int*b1s,int*b2s,char*s0,int*nb1,int*nb2,char*sn){
  FILE *fp;
  char buf[256]; int ii,nbb;

  if((fp=fopen(fname,"rt"))==NULL){
    printf("File error %s",fname); exit(0); }

  if(fscanf(fp,"%s %d",buf,&nbb)==EOF){printf("Error nb\n");exit(0);}
  if(strcmp(buf,"NB")==0){
    *nb = nbb; *nb1 = -1;
    if(fscanf(fp,"%s %d",buf,ng)==EOF){printf("Error ng\n");exit(0);}
    if(fscanf(fp,"%s",buf)==EOF){printf("Error GTS\n");exit(0);}
    for(ii=0;ii<*ng;ii++)
      if(fscanf(fp,"%d",&gts[ii])==EOF){printf("Error GTS %d\n",ii);exit(0);}
    if(fscanf(fp,"%s",buf)==EOF){printf("Error B1S\n");exit(0);}
    for(ii=0;ii<*ng;ii++)
      if(fscanf(fp,"%d",&b1s[ii])==EOF){printf("Error B1S %d\n",ii);exit(0);}
    if(fscanf(fp,"%s",buf)==EOF){printf("Error B2S\n");exit(0);}
    for(ii=0;ii<*ng;ii++)
      if(fscanf(fp,"%d",&b2s[ii])==EOF){printf("Error B2S %d\n",ii);exit(0);}
  }else{
    if(fscanf(fp,"%s %d",buf,nb2)==EOF){printf("Error nb2\n");exit(0);}
    *nb1 = nbb; *nb = 2*(*nb1)+*nb2+2;
//    ng = nb1+(6*nb1+2)*nb2
  }

  if(fscanf(fp,"%s",buf)==EOF){printf("Error S0\n");exit(0);}
  for(ii=0;ii<*nb;ii++)
    if(fscanf(fp,"%d",&s0[ii])==EOF){printf("Error S0 %d\n",ii);exit(0);}

  if(fscanf(fp,"%s",buf)==EOF){printf("Error Sn\n");exit(0);}
  sn[0] = -1;
  if(fscanf(fp,"%d",&sn[0])==EOF){printf("Error Sn0 %d\n",ii);exit(0);}
  if(sn[0] > -1)
    for(ii=1;ii<*nb;ii++)
      if(fscanf(fp,"%d",&sn[ii])==EOF){printf("Error Sn %d\n",ii);exit(0);}

  fclose(fp);
}

int foutput(char*fname,int nb,int ng,int*gts,int*b1s,int*b2s,char*s0){
  FILE *fp;
  int ii;

  if((fp=fopen(fname,"wt"))==NULL){
    printf("File error %s",fname); exit(0); }

  if(fprintf(fp,"NB %d\n",nb)==EOF){printf("Error nb\n");exit(0);}
  if(fprintf(fp,"NG %d\n",ng)==EOF){printf("Error ng\n");exit(0);}
  if(fprintf(fp,"GTS\n")==EOF){printf("Error GTS\n");exit(0);}
  for(ii=0;ii<ng;ii++)
    if(fprintf(fp,"%d ",gts[ii])==EOF){printf("Error GTS %d\n",ii);exit(0);}
  if(fprintf(fp,"\nB1S\n")==EOF){printf("Error B1S\n");exit(0);}
  for(ii=0;ii<ng;ii++)
    if(fprintf(fp,"%d ",b1s[ii])==EOF){printf("Error B1S %d\n",ii);exit(0);}
  if(fprintf(fp,"\nB2S\n")==EOF){printf("Error B2S\n");exit(0);}
  for(ii=0;ii<ng;ii++)
    if(fprintf(fp,"%d ",b2s[ii])==EOF){printf("Error B2S %d\n",ii);exit(0);}

  if(fprintf(fp,"\nS0\n")==EOF){printf("Error S0\n");exit(0);}
  for(ii=0;ii<nb;ii++)
    if(fprintf(fp,"%d",s0[ii])==EOF){printf("Error S0 %d\n",ii);exit(0);}
  if(fprintf(fp,"\n")==EOF){printf("Error S0\n");exit(0);}

  fclose(fp);
}

int fnsinput(char*fname,int*ng,char*fns){
  FILE *fp;
  char buf[256]; int ii,jj;
  char *fns1;

  if((fp=fopen(fname,"rt"))==NULL){
    printf("File error %s",fname); exit(0); }
  if(fscanf(fp,"%s %d",buf,ng)==EOF){printf("Error ng\n");exit(0);}
//  printf("Input fns: %d\n",*ng);
  for(ii=0;ii<*ng;ii++)
    for(jj=0;jj<18;jj++)
      if(fscanf(fp,"%d",fns+18*ii+jj)==EOF){printf("Error 0,%d,%d\n",ii);exit(0);}
  fns1 = fns+(*ng)*18;
  for(ii=0;ii<*ng;ii++)
    for(jj=0;jj<18;jj++)
      if(fscanf(fp,"%d",fns1+18*ii+jj)==EOF){printf("Error 1,%d,%d\n",ii);exit(0);}
  fclose(fp);
}

int fnsoutput(char*fname,int ng,char*fns){
  FILE *fp;
  int ii,jj;
  char *fns1;

  if((fp=fopen(fname,"wt"))==NULL){
    printf("File error %s",fname); exit(0); }
  printf("Output fns: %d\n",ng);
  if(fprintf(fp,"NG %d\n",ng)==EOF){printf("Error ng\n");exit(0);}
  for(ii=0;ii<ng;ii++){
    for(jj=0;jj<18;jj++)
      if(fprintf(fp,"%2d ",fns[18*ii+jj])==EOF){printf("Error 0,%d,%d\n",ii);exit(0);}
    fprintf(fp,"\n");
  }
  fns1 = fns+ng*18;
  for(ii=0;ii<ng;ii++){
    for(jj=0;jj<18;jj++)
      if(fprintf(fp,"%2d ",fns1[18*ii+jj])==EOF){printf("Error 1,%d,%d\n",ii);exit(0);}
    fprintf(fp,"\n");
  }
  fclose(fp);
}

int countatbl(int ng,char*fns,int*cnt2){
  int cnt=0,cnt1=0,ii,jj;
  for(ii=0;ii<ng;ii++)
    for(jj=0;jj<18;jj+=2){
      if(fns[ii*18+jj]>-1) cnt++;
      if(fns[(ii+ng)*18+jj]>-1) cnt++;
      if(fns[ii*18+jj]>-1 && fns[(ii+ng)*18+jj]>-1) cnt1++;
    }
  if(cnt2 != NULL) *cnt2=cnt1;
  return(cnt);
}


//char sbs[MX_NB*MX_NG];

char* trc0(int nb,char*s0,int ng,char*fns,int*b1s,int*b2s){
//  int bn;
  long long bn,ii;
  char *f1,*f2,*p1,*p2;
  char *s1;
  char buf[512];

  memcpy(sbs+0,s0,sizeof(char)*nb);
  memcpy(pts,fns,sizeof(char)*ng*18*2);

  ii = 0;
  while(ii < ng){
    bn = sbs[ii*nb+b1s[ii]];
    if(b2s[ii] > -1) bn += sbs[ii*nb+b2s[ii]]*3;
    p1 = pts+ii*18+bn*2; p2 = p1+ng*18;
    f1 = fns+ii*18+bn*2; f2 = f1+ng*18;
    s1 = sbs+(ii+1)*nb;
    if(p1[0] > -1){
      memcpy(s1,sbs+ii*nb,sizeof(char)*nb);
      s1[b1s[ii]] = f1[0];
      if(p1[1] > -1) s1[b2s[ii]] = f1[1];
      p1[0] = -1;
      ii += 1;
    }else if(p2[0] > -1){
      memcpy(s1,sbs+ii*nb,sizeof(char)*nb);
      s1[b1s[ii]] = f2[0];
      if(p2[1] > -1) s1[b2s[ii]] = f2[1];
      p2[0] = -1;
      ii += 1;
    }else{
      p1[0] = f1[0]; p1[1] = f1[1];
      p2[0] = f2[0]; p2[1] = f2[1];
      ii -= 1;
    }
//    printf("%d %d %d %d %s\n",ii,bn,b1s[ii],b2s[ii],vi2vs(buf,sbs+ii*nb,nb)); exit(0);
    if(ii < 0) return(NULL);
  }
  return(sbs+(long long)ng*nb);
}

//
// sbs[MX_NB*MX_NG],sbs1[MX_NB*MX_NG],pts[MX_NG*18*2]
// memcpy(sbs+0,s0,sizeof(char)*nb); ns0 = 1;
// default: n0=0; n1=ng
//
int trca0(int nb,int ng,char*fns,int*b1s,int*b2s,int ns0,int n0,int n1){
  int ns,nsb;
  long long bn,ii,jj;
  char *f1,*f2;
  char *s1;
  char buf[512];

  ns = ns0;
  for(ii=n0;ii<n1;ii++){
    memcpy(sbs1,sbs,sizeof(char)*nb*ns);
    nsb = 0;
    for(jj=0;jj<ns;jj++){
      bn = sbs1[jj*nb+b1s[ii]];
      if(b2s[ii] > -1) bn += sbs1[jj*nb+b2s[ii]]*3;
      f1 = fns+ii*18+bn*2; f2 = f1+ng*18;

      if(f1[0] > -1){
        s1 = sbs+nsb*nb;
        memcpy(s1,sbs1+jj*nb,sizeof(char)*nb);
        s1[b1s[ii]] = f1[0];
        if(f1[1] > -1) s1[b2s[ii]] = f1[1];
        nsb += 1;
      }
      if(f2[0] > -1){
        s1 = sbs+nsb*nb;
        memcpy(s1,sbs1+jj*nb,sizeof(char)*nb);
        s1[b1s[ii]] = f2[0];
        if(f2[1] > -1) s1[b2s[ii]] = f2[1];
        nsb += 1;
      }
    }
    ns = nsb;

    if(ii%100==1){
      printf("%d:%d\n",ii,ns);
//      printf("%d:%d  %s\n",ii,ns,vi2vs(buf,sbs+0*nb,nb));
    }
  }
//  printf("%d:%d  %s\n",ii,ns,vi2vs(buf,sbs+0*nb,nb));
  return(ns);
}


int addfns(int ng,char*pts,char*ptb){
  int ii,ngn;
  ngn = ng*18*2;
  for(ii=0;ii<ngn;ii+=2)
    if(ptb[ii] > -1){ pts[ii] = ptb[ii]; pts[ii+1] = ptb[ii+1];}
}

//
// sbs[MX_NB*MX_NG],sbs1[MX_NB*MX_NG],pts[MX_NG*18*2]
// memcpy(sbs+0,s0,sizeof(char)*nb); ns0 = 1;
// default: n0=0; n1=ng
//
// trcFB0(nb,ng,fns,s0,sn,b1s,b2s,ptb1,ptb1i,ptb2,ptb2i,n0,ni){
//
int trca10(int nb,int ng,char*fns,int*b1s,int*b2s,int*ptb1,int*ptb1i,int*ptb2,int*ptb2i){
  int ii,ns,nsb,jj;
  char *s0;
  int mx;
  mx = MX_FB0N;
  mx = 10;
  ns = 1;
  for(ii=0;ii<ng;ii++){
    ns = trca0(nb,ng,fns,b1s,b2s,ns,ii,ii+1);
    if(ii%100==50){
      memcpy(sbs1,sbs,sizeof(char)*nb*ns);
      memset(pts,0xff,sizeof(char)*ng*18*2);
      nsb=0;
      for(jj=0;jj<ns;jj++){
        memcpy(ptb,fns,sizeof(char)*ng*18*2);
        if(trcFB0(nb,ng,ptb,sbs1+jj*nb,NULL,b1s,b2s,ptb1,ptb1i,ptb2,ptb2i,ii,mx)!=0)
          memcpy(sbs+nsb*nb,sbs1+jj*nb,sizeof(char)*nb); nsb++;
        addfns(ng,pts,ptb);

      }
      ns = nsb;
      memcpy(fns,pts,sizeof(char)*ng*18*2);
    }
  }
  return(ns);
}

//open("a.txt","w").write(" ".join(map(str,[1,2,3])))
////
int main(int argc,char **argv){
  FILE *fp;
  char buf[65535],cmd[256];
  int nb,ng,ii,nb1,nb2,cnt,cnt2;
  char *sn;
  char ffnsout[] = "fnsout.txt";
  char ffnsin[] = "fnsin.txt";
  char fout[] = "out.txt";

  mallocs();

  if(argc < 2){printf("usage: til data [fns]\n"); exit(0); }

  finput(argv[1],&nb,&ng,gts,b1s,b2s,s0,&nb1,&nb2,sn0);
  if(nb1 > 0){
    ng = nb1+(6*nb1+2)*nb2;
    printf("Start.");
    nxn2(gts,b1s,b2s,nb1,nb2);
    printf("nxn2\n");
  }
  printf("nb:%ld ng:%ld nbgL,LL:%ld,%lld\n",nb,ng,(long)nb*ng,(long long)nb*ng); //exit(0);

  nextpts(nb,ng,b1s,b2s,pt0,pt0i,
    ptf1,ptf1i,ptf2,ptf2i,ptb1,ptb1i,ptb2,ptb2i,pte,ptei);

  printf("Sn0: %d\n",sn0[0]);

  if(argc >= 3) fnsinput(argv[2],&ng,fns);
  else gts2fns(fns,gts,ng);


  cnt=countatbl(ng,fns,&cnt2);printf("Link: %d, %d\n",cnt,cnt2);
  trcSn0(nb,ng,fns,sn0,pte,ptei);
  printf("  FB0: %d\n",trcFB0(nb,ng,fns,s0,NULL,b1s,b2s,ptb1,ptb1i,ptb2,ptb2i,0,MX_FB0N));
  cnt=countatbl(ng,fns,&cnt2);printf("Link: %d, %d\n",cnt,cnt2);

  sn = trc0(nb,s0,ng,fns,b1s,b2s);

  if(sn != NULL){
    if((fp=fopen(fout,"wt"))==NULL){ printf("ERROR fout\n"); exit(0);}
    fprintf(fp,"%s\n",vi2vs(buf,sn,nb));
    fclose(fp);
  }
  else printf("None\n");

  return(0);
}
