#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"

using std::cout;
using std::endl;
using std::cin;

///=============================================================================
const int    Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light
const double M_sun=1.98892*pow(10.,30); //in [kg].
const double vro_sun=226.0;
const double AU=1.495978707*pow(10.0,11.0);
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]
const double Mjupiter=1.898*pow(10,27); 
const double ddeg=0.5;
const double binary_fraction=double(2.0/3.0);
const int    GG=200;///number of bins of tE distribution
const double tE_min=0.0;
const double tE_max=100.0;  
const double stept=double((tE_max-tE_min)/GG); 
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double Avks=double(8.20922);

///=============================================================================

const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const int    N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo

///=============================================================================

///  SECTOR 12
const double Tobs=13.7;//days
const double cadence=double(30.0/60.0/24.0);//in days 
const int    M=2;///No. of filter G, T
const double detect[M]={21.0,18.0};
const double satu[M]=  {5.9, 1.0};
const double FWHM[M]=  {0.059*5.0, 21.0*3.0};//20 Sep 
const double AlAv[M]=  {1.13  , 0.8558}; ///https://iopscience.iop.org/article/10.3847/1538-3881/ab487f/pdf, 
const double sigma[M]= {0.017 , 0.02};//// G, T
const int ER=int(96); 
const double lcent=double(342.59208588 -360.0);//Sector 12
const double bcent=double(4.60638041); 
const int thre=5;  

///=============================================================================

struct source{
    int    nums,struc,cl, numt;
    double Ds,TET,FI,ratios;
    double nsbl[M], blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M];
    double lon, lat, Astar;
    double type, Tstar, logl, Rstar, mass;
    double od_disk, od_ThD, od_bulge, od_halo,opd;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart,nstarti;
    double nsdis[Num],nddis[Num];
    double ros,deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
};
struct lens{
    double Ml,Dl,ratiol,vl,vs,Vt,xls;
    double rhomaxl,tE,RE,t0, u0, tetE, mul;
    int    numl,struc;
};
struct detection{
   int    det,N_sim,N_det;
   double dchi,t,tmin,tmax;
   double ave_re,ave_vl,ave_vs,weight,deltam;
   double ave_opt,ave_opt2,ave_te,ave_u0,ave_Ds,ave_Dl;
   double ave_vt,ave_ml,ave_fred,ave_npsf,ave_tei; 
   double ave_aps,ave_ex,ave_dmag,ave_apbd,ave_bl,ave_apb;
   double ratio,nsdet,nstar,nevent,nfin; 
   double te[GG+1],effs[GG+1],effd[GG+1]; 
   double effs_com[GG+1], effd_com[GG+1];
   double tmag[ER], dmag[ER];  
   double timet;
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[2][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double Teff_b[N2],logl_b[N2],Mab_b[2][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[2][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logl_h[N4],Mab_h[2][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
};

///=============================================================================

int  TEdet(detection & d, double tE); 
int  Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s,CMD & cm, extinc & ex);
void func_lens(lens & l, source & s);
void vrel(source & s,lens & l);
void Disk_model(source & s, int , int , int);
void optical_depth(source & s);
void lensing(source & s, lens & l,detection & d,int, int, int);
double ErrorTESS(detection & d, double maga); 
double Interpol(double ds, extinc & ex);
double RandN(double sigma, double);
double RandR(double down, double up);

time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;
///===========================================================================//
///                                                                           //
///                  Main program                                             //
///                                                                           //
///===========================================================================//	
int main(){


   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   time(&_timeNow);
   printf("START time:   %s",ctime(&_timeNow));


    VBBinaryLensing vbb;
    vbb.Tol=1.e-4;
    //vbb.a1=0.0;  
    vbb.LoadESPLTable("./files/ESPL.tbl");
  
  
  
     source s;
     lens l;
     detection d;
     CMD cm;
     extinc ex;
     read_cmd(cm);
     
     
       
///=============================================================================     
     double PPM;
     FILE *err;
     err=fopen("./files/errorTESS2.txt", "r");
     if(!err){cout<<"cannot read ErrorTESS2.txt "<<endl; exit(0); }
     for(int i=0; i<ER; ++i){
     fscanf(err,"%lf     %lf   %lf\n", &d.tmag[i],&PPM, &d.dmag[i]);}
     fclose(err);
     
     
     
     
     char filename1[40];
     FILE* teefi;  FILE* fil1;   FILE* TEE;  
     FILE* result1;  
     FILE* result2; 
     result1=fopen("./files/TESS1.txt", "w");      
     result2=fopen("./files/TESS2.txt", "w");      
     fclose(result1); 
     fclose(result2);
     

///=============================================================================
     int nl=0, nb=0;
     int detectf,vv;
     double maga,effe; 
     double Nmicro=0.0,lonn,countt,countd, wei;   
     for(int i=0; i<(GG+1); ++i){d.effs_com[i]=0.0;   d.effd_com[i]=0.0;}
     double Dl[thre],Ds[thre],tE[thre],vt[thre],ml[thre], bl[thre],mb[thre],u0[thre],we[thre]; 
     double SD_Dl, SD_Ds, SD_tE, SD_vt, SD_ml, SD_bl, SD_mb, SD_u0;  
     
  
   //float(-12.0+bcent)
    for(s.lat= -0.89;    s.lat<=float(12.0+bcent);  s.lat+=ddeg){
    nb+=1;   nl=0;  
    for(lonn= float(-12.0+lcent) ;  lonn<=float(12.0+lcent);  lonn +=ddeg){      
    nl+=1;
    if(lonn<=0.0)   s.lon=360.0+lonn;
    else            s.lon=lonn;
    s.TET=double(360.0-s.lon)/RA;///radian
    s.FI= double(s.lat/RA);///radian
    Disk_model(s, 1, nb, nl);
    if(Extinction(ex,s)==1){
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>> NEW STEP <<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
    cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
    cout<<"nb:  "<<nb<<"\t nl:   "<<nl<<endl;
    
   
   
    d.ave_fred=d.ave_apb=d.ave_apbd=d.ave_dmag=d.ave_npsf=0.0; 
    d.ave_Ds=d.ave_Dl=d.ave_vt=d.ave_ml=d.ave_u0=d.ave_tei=0.0;
    d.ave_re=d.ave_vl=d.ave_vs=d.ave_te=d.nevent=0.0;
    d.ave_opt=d.ave_opt2=s.nstart=d.ratio=d.nfin=0.0;
    d.ave_bl=d.ave_aps=d.ave_ex=0.0; 
    
    for(int i=0;  i<Num;  ++i){s.nsdis[i]=0.0000145776534762; s.nddis[i]=0.0;}
    for(int i=0;  i<thre; ++i){Dl[i]=Ds[i]=tE[i]=vt[i]=ml[i]=bl[i]=mb[i]=u0[i]=we[i]=0.0;}
    for(int i=0; i<(GG+1);++i){d.effd[i]=0.0;  d.effs[i]=0.0000000754765;  d.te[i]=double(i*stept);}
    d.nstar=d.nsdet=d.nevent=0.00000645974;
    d.N_sim=d.N_det=0; 
    

    sprintf(filename1,"./files/distribution/%c%d%c%d.dat",'D',nb,'_',nl);
    fil1=fopen(filename1,"w");
    teefi=fopen("./files/tefzz.txt","w"); 
    fclose(teefi);
    
   


    do{
    do{
    func_source(s,cm,ex);
    func_lens(l,s);
    }while(l.tE>tE_max); 
    d.nstar+=1.0;
    optical_depth(s);
    s.nsdis[s.nums]+=1.0;
    maga=s.magb[1]-2.5*log10(s.blend[1] * s.Astar + 1.0-s.blend[1]);
    
    
    if(s.magb[1]>=satu[1] and maga<=detect[1]){
    d.N_sim+=1;
    d.weight=double(1.0/s.nsbl[1]); 
    d.nsdet   +=1.0*d.weight; 
    d.ave_aps += s.Map[1]*d.weight;
    d.ave_apb +=s.magb[1]*d.weight;
    d.ave_ex  +=s.Ai[1]*d.weight;
    d.ave_npsf+=s.nsbl[1]*d.weight;
    d.ave_opt +=s.opd*1.0e6*d.weight; 
    s.nddis[s.nums]+=1.0*d.weight; 
    detectf=0;  
    lensing(s,l,d, 0 ,  nb, nl);
    
    
    if(d.det==1 and d.dchi>=500.0 and d.timet<2.0 and d.timet>0.0){
    Dl[d.N_det]=l.Dl;        Ds[d.N_det]=s.Ds;    tE[d.N_det]=l.tE; 
    vt[d.N_det]=l.Vt;        ml[d.N_det]=l.Ml;    bl[d.N_det]=s.blend[1]; 
    mb[d.N_det]=s.magb[1];   u0[d.N_det]=l.u0;    we[d.N_det]= d.weight;  
    detectf=1; 
    d.N_det +=1; 
    d.nevent  +=1.0*d.weight; 
    d.ave_bl  +=s.blend[1]*d.weight; 
    d.ave_dmag+=d.deltam*d.weight;
    d.ave_apbd+=s.magb[1]*d.weight;
    d.ave_opt2+=s.opd*1.0e6*d.weight; 
    d.ave_te  +=l.tE*d.weight;
    d.ave_Ds  +=s.Ds*d.weight;
    d.ave_Dl  +=l.Dl*d.weight;
    d.ave_vt  +=l.Vt*d.weight;
    d.ave_vl  +=l.vl*d.weight;
    d.ave_vs  +=l.vs*d.weight;
    d.ave_ml  +=l.Ml*d.weight;
    d.ave_re  +=l.RE/AU*d.weight; 
    d.ave_u0  +=l.u0*d.weight;   
    if(s.cl<=4) d.ave_fred+=1.0*d.weight;}
   
   
   
    vv=int(TEdet(d,l.tE)); 
    d.effs[vv]    +=1.0;  
    d.effs_com[vv]+=1.0;
    if(detectf>0){d.effd[vv]+=1.0;  d.effd_com[vv]+=1.0;}
    teefi=fopen("./files/tefzz.txt","a+");
    fprintf(teefi,"%.7lf   %e   %d  %d\n",l.tE, d.weight, vv, detectf);
    fclose(teefi);


    if(detectf>0){
    if(d.N_det<=5 and d.N_det>0) lensing(s,l,d, int(d.N_det) , nb, nl);
    fprintf(fil1,
    "%d  %d  %d  %d  %.5lf  %.5lf  "///6
    "%d  %.5lf   %.5lf   %.5lf  "///10
    "%d   %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf   %.4lf   "   //23
    "%.5lf  %.5lf  %.8lf  %.8lf  %.2lf  %.2lf   %.4lf  %.4lf   "  //31
    "%.6lf  %.5lf  %.5lf  %.7lf  %.5lf  %.5lf  %.4lf  %.6lf  %.9lf  " ///40
    "%d  %.1lf  %d   %d  %.5lf\n", //44
    nb, nl, d.N_sim,d.N_det,s.lat, s.lon, //6
    l.struc, l.Ml, l.Dl, l.vl, //10
    s.struc, s.cl, s.mass, s.Ds, s.Tstar, s.Rstar, s.logl, s.type, l.vs, s.Mab[0], s.Mab[1], s.Map[0], s.Map[1],//23
    s.magb[0],s.magb[1], s.blend[0], s.blend[1], s.nsbl[0], s.nsbl[1], s.Ai[0], s.Ai[1], //31
    l.tE, l.RE/AU, l.t0, l.mul, l.Vt, l.u0, s.opd*1.0e6, s.ros, l.tetE,//40
    detectf,d.dchi,d.det,s.numt,d.timet);//44
    
    /*cout<<"**********************************************"<<endl;
    cout<<"N_sim: "<<d.N_sim<<"\t flag_detection: "<<detectf<<"\t d.N_det: "<<d.N_det<<endl;
    cout<<"det: "<<d.det<<"\t dci2: "<<d.dchi<<endl;
    cout<<"latitude: "<<s.lat<<"\t longtitude: "<<s.lon<<endl;
    cout<<"Ml (M_sun): "<<l.Ml<<"\t RE (AU): "<<l.RE/AU<<"\t Dl(Kpc): "<<l.Dl<<endl;
    cout<<"relative velosity (km/s): "<<l.Vt<<"\t Vs (km/s): "<<l.vs<<"\t Vl(km/s): "<<l.vl<<endl;
    cout<<"Ds: "<<s.Ds<<"\t tE (days): "<<l.tE<<"\t ros:  "<<s.ros<<endl;
    cout<<"app_mag[G-Gaia]: "<<s.Map[0]<<"\t app_mag[T-Tess]: "<<s.Map[1]<<endl;
    cout<<"app_mag_bg[G_Gaia]: "<<s.magb[0]<<"\t app_mag_bg[T_Tess]: "<<s.magb[1]<<endl;
    cout<<"belnding[G_Gaia]: "<<s.blend[0]<<"\t blending[T_Tess]: "<<s.blend[1]<<endl;
    cout<<"optical_depth(x10^6): "<<s.opd*1.0e6<<"\t u0: "<<l.u0<<endl;
    cout<<"N_blend[G_Gaia]: "<<s.nsbl[0]<<"\t N_blend[T_Tess]: "<<s.nsbl[0]<<endl; 
    cout<<"*********************************************************************"<<endl;*/
    }}
    
    }while(d.N_det< thre);
    fclose(fil1); 


///=============================================================================

     TEE=fopen("./files/tef_last.txt","w"); 
     for(int j=0;j<(GG+1); ++j){
     d.effd[j]=double(d.effd[j]/d.effs[j]);
     fprintf(TEE,"%d   %.3lf     %.2lf    %.2lf\n",j,d.te[j],d.effs_com[j],d.effd_com[j]);}     
     fclose(TEE);
     countd=0.0; d.ave_tei=0.0; 
     teefi=fopen("./files/tefzz.txt","r"); 
     if(!teefi){cout<<"cannot make file teefi:  "<<endl; exit(0);}
     while(!feof(teefi)){
     fscanf(teefi,"%lf   %lf  %d   %d\n",&l.tE,&wei,&vv,&detectf);
     if(detectf>0 and l.tE>0.0){
     countd += wei;  
     d.ave_tei+=fabs(d.effd[vv]/l.tE)*wei;}}
     d.ave_tei=double(d.ave_tei/countd); 
     fclose(teefi);



///=============================================================================
    d.ave_te= d.ave_te/d.nevent;
    d.ave_re= d.ave_re/d.nevent;
    d.ave_Ds= d.ave_Ds/d.nevent;
    d.ave_Dl= d.ave_Dl/d.nevent;
    d.ave_vt= d.ave_vt/d.nevent;
    d.ave_vl= d.ave_vl/d.nevent;
    d.ave_vs= d.ave_vs/d.nevent;
    d.ave_ml= d.ave_ml/d.nevent;
    d.ave_u0= d.ave_u0/d.nevent;
    d.ave_bl= d.ave_bl/d.nevent;
    d.ave_apbd=d.ave_apbd/d.nevent;
    d.ave_opt2=d.ave_opt2/d.nevent;
    d.ave_dmag=d.ave_dmag/d.nevent;
    d.ave_fred= d.ave_fred*100.0/d.nevent;
    d.ave_aps= d.ave_aps/d.nsdet;
    d.ave_apb= d.ave_apb/d.nsdet;
    d.ave_ex=   d.ave_ex/d.nsdet;
    d.ave_npsf=d.ave_npsf/d.nsdet;
    d.ave_opt= d.ave_opt/d.nsdet;
    
    s.nstart=0.0; 
    d.ratio= 2.0*d.ave_opt2*1.0e-6*d.ave_tei/M_PI;///this is changed 24 Sep 
    for(int i=1; i<Num; ++i){
    effe=double(s.nddis[i]/(s.nsdis[i]+0.0000064646745));
    s.nstart += s.Rostari[i]*(s.Nstart/s.Rostart)*effe;}
    d.nfin=d.ratio*Tobs*s.nstart;//  in degree-squared
    Nmicro+=d.nfin;
    result1=fopen("./files/TESS1.txt", "a+");
    fprintf(result1,
    "%d   %d   %.2lf  %.2lf  %.1lf  %d   %d   %.8lf  %.8lf  %.8lf  " //10
    "%.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  "//19
    "%.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  " //27
    "%e  %e  %e  %e  %e  %e  %e   %e   %e\n",  // 36
    nb, nl, s.lat,s.lon,d.nstar,d.N_sim,d.N_det,d.nevent,double(d.nevent*100.0/d.nsdet), double(d.nsdet*100.0/d.nstar),  //10
    d.ave_te,d.ave_re,d.ave_Ds,d.ave_Dl,d.ave_vt,d.ave_vl,d.ave_vs,d.ave_ml,d.ave_u0,//19
    d.ave_bl,d.ave_apbd, d.ave_dmag, d.ave_fred,d.ave_aps,d.ave_apb, d.ave_ex, d.ave_npsf,  //27
    s.Nstart,s.Rostart,d.ave_tei,d.ratio,d.ave_opt,d.nfin,Nmicro,s.nstart,d.ave_opt2);//36
    fclose(result1);
    
    
    /*
    SD_Dl=SD_Ds=SD_tE=SD_vt=SD_ml=SD_bl=SD_mb=SD_u0=0.0;  wei=0.0;
    for(int i=0; i<thre; ++i){ 
    wei+= we[i]; 
    SD_Dl+=pow(Dl[i]-d.ave_Dl,2.0)*we[i]; 
    SD_Ds+=pow(Ds[i]-d.ave_Ds,2.0)*we[i]; 
    SD_tE+=pow(tE[i]-d.ave_te,2.0)*we[i];    
    SD_vt+=pow(vt[i]-d.ave_vt,2.0)*we[i]; 
    SD_ml+=pow(ml[i]-d.ave_ml,2.0)*we[i];  
    SD_bl+=pow(bl[i]-d.ave_bl,2.0)*we[i];  
    SD_mb+=pow(mb[i]-d.ave_apbd,2.0)*we[i];  
    SD_u0+=pow(u0[i]-d.ave_u0,2.0)*we[i];}
    SD_Dl=sqrt(SD_Dl/(wei*(thre-1.0)/thre));
    SD_Ds=sqrt(SD_Ds/(wei*(thre-1.0)/thre));
    SD_tE=sqrt(SD_tE/(wei*(thre-1.0)/thre));
    SD_vt=sqrt(SD_vt/(wei*(thre-1.0)/thre));
    SD_ml=sqrt(SD_ml/(wei*(thre-1.0)/thre));
    SD_bl=sqrt(SD_bl/(wei*(thre-1.0)/thre));    
    SD_mb=sqrt(SD_mb/(wei*(thre-1.0)/thre));
    SD_u0=sqrt(SD_u0/(wei*(thre-1.0)/thre));  
    result2=fopen("./files/TESS2.txt", "a+");      
    fprintf(result2,"%d  %d  %.2lf %.2lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.5lf   %.5lf\n",
    nb, nl, s.lat,s.lon, SD_Dl, SD_Ds, SD_tE, SD_vt, SD_ml, SD_bl, SD_mb, SD_u0);       
    fclose(result2); 
    */
    
///=============================================================================

    cout<<"=============================================================="<<endl;
    cout<<" ******  latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
    cout<<"fraction_red_gients[%]: "<<d.ave_fred<<endl;
    cout<<"<Ds>: "<<d.ave_Ds<<"\t <Dl>: "<<d.ave_Dl<<endl;
    cout<<"<vl>: "<<d.ave_vl<<"\t <vs>: "<<d.ave_vs<<endl;
    cout<<"<Ml>: "<<d.ave_ml<<"\t <u0>: "<<d.ave_u0<<endl;
    cout<<"<Vt>: "<<d.ave_vt<<"\t <tE>: "<<d.ave_te<<endl;
    cout<<"<Map_sources>:"<<d.ave_aps<<"\t <EXtinction>: "<<d.ave_ex<<endl;
    cout<<"N_sim :  "<<d.N_sim<<"\t N_det:  "<<d.N_det<<"\t nstar:  "<<d.nstar<<endl;
    cout<<"<M_bg_detected>: "<<d.ave_apb<<"\t <M_bg_lensed>: "<<d.ave_apbd<<endl;
    cout<<"<effi_starDet[%]:  "<<double(d.nsdet*100.0/d.nstar)<<"\t effi_Microlensing[%]: "<<double(d.nevent*100.0/d.nsdet)<<endl;
    cout<<"No. lensing event[l,b]: "<<d.nfin<<"\t n_start[Omega_l] "<<s.nstart<<endl;
    cout<<"<opt>: "<<d.ave_opt<<"\t <opt2>: "<<d.ave_opt2<<"\t event_ra: "<<d.ratio<<endl;
    cout<<"=============================================================="<<endl;
     
     
     }//end of if extinction
     }}///end of subfield
     fclose(_randStream);
     return(0);
}
///==========================================================================//
///                                                                          //
///                   Lensing                                                //
///                                                                          //
///==========================================================================//
void lensing(source & s, lens & l,detection & d , int flagi, int nb, int nl)
{
    double maga,maga2,chi2,chi1,emt, u, As, sign, timee;
    VBBinaryLensing vbb;
    vbb.Tol=1.e-4;
    //vbb.a1=0.0;  
    vbb.LoadESPLTable("./files/ESPL.tbl");
  
    
    FILE *test;  
    FILE *test2; 
    char filenam1[40], filenam2[40];
    
    
    d.tmin=l.t0-1.5*l.tE;
    d.tmax=l.t0+1.5*l.tE;
    if(d.tmin> 0.0)    d.tmin=0.0; 
    if(d.tmax<Tobs)    d.tmax=Tobs;  
    double dt =double(l.tE/200.0);
    if(dt>cadence) dt=double(cadence/5.0); 
    int Nb=int(1.0+(d.tmax-d.tmin)/cadence)+2;
    int flag[Nb], save; 

    if(flagi>0){
    sprintf(filenam1,"./files/light/%c%d%c%d%c%d.dat",'L',nb,'_',nl,'_',flagi);
    test=fopen(filenam1,"w");
    sprintf(filenam2,"./files/light/%c%d%c%d%c%d.dat",'M',nb,'_',nl,'_',flagi);
    test2=fopen(filenam2,"w");}
    
    
    d.det=0; 
    s.numt=0;
    d.deltam=0.0; 
    chi1=chi2=d.dchi=0.0;
    for(int i=0; i<Nb; ++i) flag[i]=0;
    timee=0.0; 
    d.timet=-1.0; 
    double Asmax= vbb.ESPLMag2(l.u0, s.ros);
    double DEL= fabs( Asmax*s.blend[1]+1.0-s.blend[1] -1.0 )*0.1;
    
        
    for(d.t=d.tmin; d.t<=d.tmax; d.t=d.t+dt){
    timee+=dt; 
    u=sqrt(l.u0*l.u0+(d.t-l.t0)*(d.t-l.t0)/l.tE/l.tE);
    if(u==0.0) u=1.0e-50;
    if(u>float(s.ros*10.0)) As=(u*u+2.0)/(u*sqrt(u*u+4.0));
    else                    As=vbb.ESPLMag2(u, s.ros);   
    if(As<0.9){cout<<"ERROR AS: "<<As<<"\t u: "<<u<<"\t R:"<<s.ros<<endl; exit(0);}
    maga=s.magb[1]-2.5*log10(As*s.blend[1]+1.0-s.blend[1]);
    if(flagi>0)  fprintf(test2,"%.5lf  %.5lf  %.5lf\n",d.t, (d.t-l.t0)/l.tE, maga);
    
    
    if(fabs(As*s.blend[1]+1.0-s.blend[1]-1.0)>DEL and d.timet<0.0) d.timet=fabs(l.t0-d.t)/Tobs; 
    
    
    save=0; 
    if(timee>=cadence) {timee-=cadence; save=1;}
    if(maga>=satu[1] and maga<=detect[1] and d.t>0.0 and d.t<Tobs and save>0){
    emt= ErrorTESS(d, maga);  
    d.deltam+= emt; 
    sign=RandR(0.0,1.0);  
    if(sign<0.5)  sign=-1.0; 
    else          sign= 1.0;  
    maga2=maga + fabs(emt*RandN(1.0 , 2.5))*sign;
    chi1 += pow((maga2-s.magb[1])/emt,2.0);
    chi2 += pow((maga2-     maga)/emt,2.0);
    if(d.det==0){
    if(fabs(maga2-s.magb[1])<5.0*emt) flag[s.numt]=0;
    else{
    if(maga2<s.magb[1])  flag[s.numt]=+1; 
    else                 flag[s.numt]=-1;}
    if(s.numt>=3 and float(flag[s.numt]+flag[s.numt-1]+flag[s.numt-2]+flag[s.numt-3])>3.0) d.det=1;}
    s.numt+=1;
    if(flagi>0) fprintf(test,"%.5lf  %.5lf  %.5lf  %.5lf  %.1f\n",d.t,(d.t-l.t0)/l.tE,maga2, emt, sign);}
    }//time loop 
    
    d.dchi= fabs(chi2-chi1);
    if(flagi>0){fclose(test);    fclose(test2);  }
    d.deltam=double(d.deltam/s.numt/1.00);
    if(d.timet<0.0){
    cout<<"d.timet<0.0:    "<<d.timet<<endl;  
    cout<<"magb:  "<<s.magb[1]<<"\t t0:  "<<l.t0<<"\t tE:  "<<l.tE<<"\t ros: "<<s.ros<<endl;  
    cout<<"u0:  "<<l.u0<<"\t blend[1]:  "<<s.blend[1]<<"\t Asmax:  "<<Asmax<<"\t DEL:  "<<DEL<<endl; }
    //int uue; cin>>uue;  } 
}
///#############################################################################
double ErrorTESS(detection & d, double maga){  

    double emt=-1.0, m;      
    if(maga<d.tmag[0])   emt=d.dmag[0];  
    
    else if(maga>d.tmag[ER-1] or maga==d.tmag[ER-1]){
    m=  double(d.dmag[ER-1]-d.dmag[ER-2])/(d.tmag[ER-1]-d.tmag[ER-2]);  
    emt=double(d.dmag[ER-1]+m*(maga-d.tmag[ER-1])); }
    
    else{
    for(int i=0; i<int(ER-1); i++){
    if( double((maga-d.tmag[i])*(maga-d.tmag[i+1]))<0.0 or maga==d.tmag[i] ){
    m=  double(d.dmag[i+1]-d.dmag[i])/(d.tmag[i+1]-d.tmag[i]);  
    emt=double(d.dmag[i]+m*(maga-d.tmag[i]));
    break;}}}
    //if(emt<0.0 or emt<0.000001 or emt>0.5 or maga<0.0){
    //cout<<"Error emt:  "<<emt<<"\t maga:  "<<maga<<endl; }// exit(0);}
    return(emt); 
}
///#############################################################################
int TEdet(detection & d, double tE){
    int vf=-1; 
    for(int j=0; j<GG; ++j){
    if( (tE-d.te[j])*(tE-d.te[j+1])<0.0 or tE==d.te[j] ){vf=j; break; }}
    if(vf==-1 and (tE>tE_max or tE==tE_max)) vf=GG-1; 
    if(vf==-1 and (tE<tE_min or tE==tE_min)) vf=0;    
    if(vf<0 or tE<0.0 or vf>(GG-1)){cout<<"ERROR tE: "<<tE<<"\t vf: "<<vf<<endl;  exit(0);}
    return(vf); 
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] and ds<ex.dis[i+1]){
  F = ex.Extks[i]+(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0 or F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; exit(0); }
  return(F);
}
///==================================================================
double RandN(double sigma, double nn){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nn; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return(rr);
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye;
    double metal, age, gravity, MB, MV, MK, MI, G, T, GBP, GRP;
    char filename[40];
    FILE *fp2;
////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','i','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB,&MV, &MI, &MK, &cm.cl_d[j], &cm.type_d[j]);
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428- 0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                   +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430 + RandN(0.6,1.5);
    cm.Mab_d[0][j]= G;  
    cm.Mab_d[1][j]= T;
    if(cm.mass_d[j]<=0.0 or cm.Teff_d[j]<0.0 or metal>0.1 or age>10 or int(cm.cl_d[j])==6 or cm.type_d[j]>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<"\t metal:  "<<metal<<"\t age:  "<<age<<endl;
    cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;





////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','b','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,&MV,&MI,&MK,&cm.cl_b[j],&cm.type_b[j]);
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
    cm.Mab_b[0][j]= G;  
    cm.Mab_b[1][j]= T;
    if(cm.mass_b[j]<=0.0 or cm.Teff_b[j]<0.0 or age>10 or metal>0.9 or cm.cl_b[j]==6 or cm.type_b[j]>=8.0){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<"\t age:  "<<age<<"\t metal:  "<<metal<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<"\t Teff:  "<<cm.Teff_b[j]<<endl;  cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;






////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','k','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,&MV,&MI,&MK,&cm.cl_t[j],&cm.type_t[j]);
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
    cm.Mab_t[0][j]= G;  
    cm.Mab_t[1][j]= T;
    if(cm.mass_t[j]<=0.0 or cm.Teff_t[j]<0.0 or metal>0.2 or cm.cl_t[j]==6 or cm.type_t[j]>=8.0){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;




////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','h','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,&MV,&MI,&MK,&cm.cl_h[j],&cm.type_h[j]);
    
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
    cm.Mab_h[0][j]= G;  
    cm.Mab_h[1][j]= T;
    
    if(cm.mass_h[j]<=0.0 or age<0 or cm.cl_h[j]==6 or cm.Teff_h[j]<0.0 or metal>0.1 or cm.cl_h[j]>7 or cm.type_h[j]>9 or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"age:  "<<age<<"\t metal:  "<<metal<<"\t mass:  "<<cm.mass_h[j]<<endl; cin>>yye;}
    
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///==============================================================//
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD & cm, extinc & ex)
{
    int struc,nums, num;
    double rho,rf,Nblend[M], Mab[M], Map[M];
    double Ds,Ai[M],Av;
    double maxnb=0.0;  


    for(int i=0; i<M; ++i){
    Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    Nblend[i]=Nblend[i]+RandN(sqrt(Nblend[i]),1.5);
    if(Nblend[i]<=1.0)  Nblend[i]=1.0;
    if(Nblend[i]>maxnb) maxnb=Nblend[i];}
    
    for(int i=0; i<M; ++i){s.Fluxb[i]=0.0; s.magb[i]=0.0; s.Ai[i]=0.0; s.Map[i]=0.0; s.Mab[i]=0.0;}

   

    for(int k=1; k<=maxnb; ++k){
    
    do{
    nums=int(fabs((double)rand()/(double)(RAND_MAX+1.)*Num*1.0));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[nums] or nums<5);///distance larger than 50.0 
    Ds=double(nums*step);///in kpc
    //nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums; }
    //cout<<"Ds: "<<s.Ds<<"\t nums: "<<s.nums<<endl; 


    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
         if (rf<= s.rho_disk[nums]) struc=0;///thin disk
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
    else if (  rf<=s.Rostar0[nums]) struc=3;///halo
    if(k==1)    s.struc=struc;



    if(struc==0){///thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){
    s.Rstar=cm.Rs_d[num];
    s.type= cm.type_d[num];
    s.mass= cm.mass_d[num];
    s.Tstar=cm.Teff_d[num];
    s.logl= cm.logl_d[num];
    s.cl= cm.cl_d[num]; }
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])>5 or cm.Teff_d[num]<0.0 or cm.type_d[num]>8.0 or cm.mass_d[num]>1000.0){
    cout<<"Error(thin disk) temp: "<<cm.Teff_d[num]<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; exit(0); } }
 




    if(struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}//cout<<"Mab[i]: "<<s.Mab[i]<<endl; }
    if(k==1){
    s.Rstar=cm.Rs_b[num];
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Tstar=cm.Teff_b[num];
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])>5 or cm.Teff_b[num]<0.0 or cm.type_b[num]>8.0 or cm.mass_b[num]>10000.0){
    cout<<"Error(bulge) temp: "<<cm.Teff_b[num]<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   exit(0); }}




    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }//cout<<"Mab[i]: "<<s.Mab[i]<<endl;}
    if(k==1){
    s.Rstar= cm.Rs_t[num];
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Tstar=cm.Teff_t[num];
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])>5 or cm.Teff_t[num]<0.0 or cm.type_t[num]>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<cm.Teff_t[num]<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  exit(0);}}




    if(struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-1.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_h[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_h[num];
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Tstar=cm.Teff_h[num];
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])>5 or cm.Teff_h[num]<0.0 or cm.type_h[num]>8.0 or cm.mass_h[num]>10000.0){
    cout<<"Error(Galactic halo) temp: "<<cm.Teff_h[num]<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  exit(0);}}

   

    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av<0.0)    Av=0.0;
    for(int i=0; i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i] + 5.0*log10(Ds*100.0) + Ai[i];
    if(Nblend[i]>=k){s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));} }
    if(k==1){
    for(int i=0; i<M;  ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}}
    }///loop over the stars

   

    for(int i=0; i<M; ++i){
    s.magb[i]=-2.5*log10(fabs(s.Fluxb[i]));
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.nsbl[i]=double(Nblend[i]);
    if(int(s.nsbl[i])<1.0 or (s.nsbl[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or s.Fluxb[i]<0.0 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; exit(0); }}
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0 or s.Rstar<0.0 or s.mass>10000.0 or
    s.nums>Num or s.nums<=0 or Av<0.0 or s.cl<0 or s.cl>=6){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<"\t cl:  "<<s.cl<<endl; exit(0); }
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0)     Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lat>10.0) Lat=10.0; 
     if(Lon>100.0 and Lon<260.0) Lon=100.0; 
     if(Lon>360.000 or Lon<0.25 or fabs(Lat)>10.0 or (Lon>100 and Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");


     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }
     }}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{

    double f,test, t1, t2;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    for(int k=1; k<=s.nums; ++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}
    double mmin=double(13.0*Mjupiter/M_sun); 



    do{
    l.numl =int(RandR(1.0,s.nums-1.0));
    test =RandR(0.0,l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  exit(0); }
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);




   double  randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<=s.rho_disk[l.numl])    l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
  else    {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}




///#############################################################################

  if(l.struc==0){//thin disk
  t1=pow(2.0,-3.0);
  t2=pow(mmin,-0.7)*pow(0.08,0.7-1.6);
  do{
  l.Ml=RandR(mmin, 2.0);
  test=RandR(t1,t2);
  if(l.Ml<0.08)                     f=pow(l.Ml,-0.7)*pow(0.08,0.7-1.6);   
  else if(l.Ml<1.0 and l.Ml>=0.08)  f=pow(l.Ml,-1.6);
  else if(l.Ml>=1.0)                f=pow(l.Ml,-3.0);
  if(t1>t2 or f<t1 or f>t2){cout<<"Error_1 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}



  if(l.struc==1){///Galactic bulge
  t1=pow(1.4,-2.35)*pow(0.08,2.35-0.7); 
  t2=pow(mmin,-0.7);
  do{
  l.Ml=RandR(mmin, 1.4);
  test=RandR(t1,t2);
  if(l.Ml<0.08)         f=pow(l.Ml, -0.7); 
  else if(l.Ml>=0.08)   f=pow(l.Ml,-2.35)*pow(0.08,2.35-0.7);
  if(t1>t2 or f<t1 or f>t2 or l.Ml>1.4){cout<<"Error_2 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}



  if(l.struc==2  or l.struc==3){///thick disk  &&  Stellar Halo
  t1=pow(1.4,-0.5)*pow(0.08,0.5-0.7); 
  t2=pow(mmin,-0.7);
  do{
  l.Ml=RandR(mmin, 1.4); 
  test=RandR(t1,t2); 
  if(l.Ml<0.08)          f=pow(l.Ml,-0.7);   
  else if(l.Ml>=0.08)    f=pow(l.Ml,-0.5)*pow(0.08,0.5-0.7);
  if(t1>t2 or f<t1 or f>t2 or l.Ml>1.4){cout<<"Error_3 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}

///#############################################################################



  l.xls=l.Dl/s.Ds;
  l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
  l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
  vrel(s,l);
  l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
  s.ros=double(s.Rstar*Rsun*l.xls/l.RE);
  l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
  l.u0=fabs(RandR(0.0,1.0));
  if(l.u0==0.0) l.u0=1.0e-50; 
  l.tetE= double(l.RE/AU/l.Dl);///marcs
  l.t0=RandR(1.0 , Tobs-1.0);
  s.Astar= (l.u0*l.u0+2.0)/sqrt(l.u0*l.u0*(l.u0*l.u0+4.0)); 
    
  if(l.tE<0.0 or l.tE==0.0 or l.RE<0.0 or  l.Dl>s.Ds  or l.xls>=1.0 or s.ros<0.0 or l.t0>Tobs or l.Ml<0.0){ 
  cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<endl; 
  cout<<"xls:  "<<l.xls<<"\t ros:  "<<s.ros<<"\t t0:  "<<l.t0<<endl; exit(0);}    
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;


  double NN=2.5;
  double sigma_R_Disk, sigma_T_Disk, sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){ Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}



  for (int i=0;i<2; ++i){
  test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                           {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))                  {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))           {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))    {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4])){sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///======================================================================================
 if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}

///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= pi-fabs(tetd);
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
   
    if (l.Vt<0.0 or l.Vt>1.0e6 or l.Vt==0.0){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
//cout<<"Vt: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
}
///==================================================================
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int flag, int nb, int nl)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
   char filename[40];
   FILE *fill;

  
   if(flag>0){  
   sprintf(filename,"./files/density/%c%d%c%d.dat",'D',nb,'_', nl);
   fill=fopen(filename,"w");
   if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}}



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 
s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]
s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   if(flag>0)
   fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 
   if(flag>0)   fclose(fill);
 // fclose(fill);

 // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
 //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
