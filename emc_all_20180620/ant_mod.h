#ifndef _ANTMOD_H
#define _ANTMOD_H
#include <math.h>
#include "other_calc_mod.h"

#ifdef CPPDLL_EXPORTS
#define CPP_EXPORTS __declspec(dllexport)
#else
#define CPP_EXPORTS __declspec(dllimport)
#endif

struct para_1336sec
{
	double fai,fai3,thita,beta,ka,kp,kh,kv;
	int flag,flag_tilt;
};
struct CPP_EXPORTS para_imt5Gant
{//5G天线参数组
	//
	double phi, theta;//天线水平角，天线俯仰角（不需要用户输入）
	int rownum,colnum;//天线阵列行数、列数
	double G_Em, Am, SLA, phi_3db, theta_3db, theta_e,phi_e;
	//天线最大增益（dBi）,Am，SLA，天线半功率波束宽度（水平，度）
	//天线半功率波束宽度（水平，度）,天线半功率波束宽度（垂直，度）
	//水平电子倾角（度），俯仰电子倾角（度）（不需要用户输入）
	double freq;//频率（MHz）（不需要用户输入）

};

double Antepattern_1851_csc2(double theta_3,double theta_max,double floorlev,double theta)
{   
    double eps = 1e-20;
    double theta_1 = theta_3;
    double G = 0;
	double u;
	double G_theta1;
    if (theta >= -theta_3/0.88 && theta <= theta_3)
    {
		u = pi * 50.8 *sin(theta/180*pi)/theta_3;
        if (abs(u) < eps)
            G = 1;
        else
            G = sin(u)/u;
        G = 20*log10(G);
	}
    else if (theta > theta_3 && theta <= theta_max)
    {
		G_theta1 = sin( pi * 50.8 *sin(theta_1/180*pi)/theta_3 )/( pi * 50.8 *sin(theta_1/180*pi)/theta_3);
        G = G_theta1 * pow( csc(theta/180*pi)/csc(theta_1/180*pi) ,2.0);
        G = 20*log10(G);
	}
    else if (theta > theta_max && theta <= 90 )
        G = floorlev;
	return G;
}


void cal_gain_M1851(double theta,double theta_3dB,double FSSL,double Gout[3])
{    
    theta = mod(theta,360);
	
    if (theta > 180)
        theta = theta - 360;

    double theta_r = theta/180*pi;
    double eps = 1e-20;
	double u,G,G_peak,G_average;
	double Fu;
	if (FSSL >-13.2)
	{
		G = 0;
		G_peak = 0;
		G_average = 0;
	}
	else if (FSSL > -23)
	{
		u = pi *50.8*sin(theta_r)/theta_3dB;
		if (u == 0) 
			Fu = 1;
		else
			Fu = sin(u)/u;
		G = 20*log10(abs(Fu)); 
		if (G >= -5.75 && ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_peak = G;
		else 
		{
			G_peak = -8.584*log(2.876*abs(theta)/theta_3dB);
			if (G_peak < -30)
			  G_peak = -30;
		}
		  if (G >= -12.16 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_average = G;          
		  else 
		  {
			  G_average = -8.584*log(2.876*abs(theta)/theta_3dB) -3.72;
			  if (G_average < -30)
				  G_average = -30;
		  }
	}
	else if (FSSL > -32)
	{
		u = pi *68.8*sin(theta_r)/theta_3dB;
		if (abs(u-pi/2) < eps)
			Fu = 1/2; 
		else
			Fu = pi/2*cos(u)/( pi*pi/4 - u*u );

		G = 20*log10(abs(Fu)) - 20*log10(2/pi);
		if( G >= -14.4 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			G_peak = G;
		else 
		{
			G_peak = -17.51*log(2.33*abs(theta)/theta_3dB);
			if (G_peak < -50)
				G_peak = -50;
		}
		if (G >= -20.6 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			G_average = G;          
		else 
		{
			G_average = -17.51*log(2.33*abs(theta)/theta_3dB) -4.32;
			if (G_average < -50)
				G_average = -50;
		 }
	  }
	  else if (FSSL > -40)
	  {
		  u = pi *83.2*sin(theta_r)/theta_3dB;
		  if (u == 0)
			  Fu = 1/2;
		  else if (abs(u-pi)<eps)
			  Fu = pi/4;
		  else 
			  Fu = pi*pi/2/u*sin(u)/( pi*pi - u*u);
		  G = 20*log10(abs(Fu)) - 20*log10(1.0/2);
		  if (G >= -22.3 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_peak = G;
		  else 
		  {
			G_peak = -26.882*log(1.962*abs(theta)/theta_3dB);
			  if (G_peak < -60)
				  G_peak = -60;
		  }
		  if (G >= -29.0 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_average = G;          
		  else 
		  {
				G_average = -26.882*log(1.962*abs(theta)/theta_3dB) -4.6;
			  if( G_average < -60)
				  G_average = -60;
 		  }
	  }
	  else 
	  {
		  u = pi *95*sin(theta_r)/theta_3dB;
		  if( abs(u-pi/2)<eps)
			  Fu = 3/8;
		  else if (abs(u-pi*1.5)<eps)
			  Fu = pi/4;
		  else 
			  Fu = 3/8*pi*cos(u)*( 1/( pi*pi/4 - u*u ) - 1/( 9*pi*pi/4 - u*u) );

		  G = 20*log10(abs(Fu)) - 20*log10( 4/(3*pi ) );   
		   if( G >= -31.5 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_peak = G;
		  else 
			{
				G_peak = -35.84*log(1.756*abs(theta)/theta_3dB);
			  if (G_peak < -70)
				  G_peak = -70;
		   }
		  if (G >= -37.6 &&  ( mod(theta,360) < 90 || mod(theta,360) > 270 ))
			  G_average = G;          
		  else 
			{
				G_average = -35.84*log(1.756*abs(theta)/theta_3dB) -4.2;
			  if (G_average < -70)
				  G_average = -70;
  			 }
	  }
    Gout[0] = G;
	Gout[1] = G_average;
	Gout[2] = G_peak;
	return;
}
double A_H_Pattern(double theta)
{
	double A_H;
    if( abs(theta) < 0.5)
        A_H = 30;
    else if(abs(theta) < 2)
        A_H = 20;
    else if(abs(theta) < 5)
        A_H = 10;
    else if(abs(theta) < 40)
        A_H = 2.5;
    else 
        A_H = -5;
    return 30-A_H;
    
}

double A_V_Pattern(double phi)
{
	double A_V;
    if(phi >= 10)
        A_V = -5;
    else if(phi < 10 && phi >= 5)
        A_V = 10;
    else if(phi < 5 && phi >= -30)
        A_V = 30;
    else if(phi < -30 && phi >= -60)
        A_V = 2.5;
    else 
        A_V = -5;
	return 30-A_V;
}
double Antenna_1466_No3(double theta,double phi)
{    
    double A_H = A_H_Pattern(theta);
    double A_V = A_V_Pattern(phi);
    return 30 - ( A_H + A_V );
    
}
double G2D(double Freq_MHz,double G0,double Efficiency=0.6) 
{
	double C = 300000000,D;
	double lmd = C / Freq_MHz / 1000000;
	if(Efficiency == 0.6 )
		D = lmd * pow(10 , ((G0 - 7.7) / 20));
	else
		D = lmd * pow(10 , ((G0 - 20 * log10(sqrt(Efficiency) * pi)) / 20));
	return D;
}
double D2G(double Freq_MHz,double  Ant_Diameter ,double  Efficiency=0.6) 
{
	double C = 300000000;

	double lmd = C / Freq_MHz / 1000000;
	return 10 * log10(Efficiency * (pow(pi * Ant_Diameter / lmd,2) ));
}


double Gfai_F1245_D(double Freq_MHz,double Ant_Diameter,double AzimuthOffset,double Efficiency=0.6) 
{
	double C = 300000000;
	double d = Ant_Diameter;
	double lamda = C / Freq_MHz / 1000000;
	double fai = AzimuthOffset;
	double G;
	if(fai>180)
    fai = fai-360;
	if(fai<-180)
    fai = fai+360;
	double ita = Efficiency;
	double Gmax = D2G(Freq_MHz, d, ita);
	double G1 = 2 + 15 * log10(d / lamda);
	double faim = 20 * (lamda / d) * sqrt(Gmax - G1);
	double fair = 12.02 * (pow(d / lamda,-0.6));
	if(Freq_MHz < 1000 || Freq_MHz > 70000 )
		G=0;
    
	else if(d / lamda > 100 )
	{
		if(abs(fai) < faim )
			G = Gmax - 2.5 * 0.001 * pow(abs(fai) * d / lamda,2);
		else if(abs(fai) < (faim>fair?faim:fair) )
			G = G1;
		else if(abs(fai) < 48 )
			G = 29 - 25 * log10(abs(fai));
		else
			G = -13;
	}
	else
    {
		if(abs(fai) < faim )
			G = Gmax - 2.5 * 0.001 * pow(abs(fai) * d / lamda,2);
		else if(abs(fai) < 48 )
			G = 39 - 5 * log10(d / lamda) - 25 * log10(abs(fai));
		else
			G = -3 - 5 * log10(d / lamda);
	}
	return G;
}
double Gfai_F1245_G0(double Freq_MHz,double Gmainbeam ,double AzimuthOffset ,double Efficiency=0.6) 
{
	return Gfai_F1245_D(Freq_MHz, G2D(Freq_MHz, Gmainbeam,Efficiency), AzimuthOffset, Efficiency);
}
double Gfai_F699_D(double Freq_MHz,double Ant_Diameter,double AzimuthOffset,double Efficiency=0.6) 
{
	
	double C = 300000000;
	double d = Ant_Diameter;
	double lamda = C / Freq_MHz / 1000000;
	double fai = AzimuthOffset;
	double ita = Efficiency;
	double Gmax = D2G(Freq_MHz, d, ita);

	double G1 = 2 + 15 * log10(d / lamda);

	double faim = 20 * lamda / d * sqrt(Gmax - G1);

	double fair = 15.85 * (pow(d / lamda,-0.6));

	double fais = 144.5 * (pow(d / lamda,-0.2));
	double G;
	if (Freq_MHz < 100 || Freq_MHz > 70000 || d / lamda <= 0.63 )
		G = 0;
    
	else if(Freq_MHz < 1000 )
	{
		if( abs(fai) < faim )
			G = Gmax - 2.5 * 0.001 *pow (abs(fai) * d / lamda,2) ;
		else if( abs(fai) < (100 * lamda / d) )
			G = G1;
		else if( abs(fai) < fais )
			G = 52 - 10 * log10(d / lamda) - 25 * log10(abs(fai));
		else
			G = -2 - 5 * log10(d / lamda);
	}
	else if(d / lamda > 100 )
	{
		if (abs(fai) < faim )
			G = Gmax - 2.5 * 0.001 * pow(abs(fai) * d / lamda,2);
		else if (abs(fai) < fair )
			G = G1;
		else if (abs(fai) < 48 )
			G = 32 - 25 * log10(abs(fai));
		else
			G = -10;
	}
	else
	{    
		if( abs(fai) < faim )
			G = Gmax - 2.5 * 0.001 * pow(abs(fai) * d / lamda,2);
		else if( abs(fai) < (100 * lamda / d) )
			G = G1;
		else if( abs(fai) < 48 )
			G = 52 - 10 * log10(d / lamda) - 25 * log10(abs(fai));
		else
			G = 10 - 10 * log10(d / lamda);
	}
	return G;
}
double Gfai_F699_G0( double Freq_MHz ,double  Gmainbeam ,double  AzimuthOffset ,double Efficiency=0.6) 
{
	return Gfai_F699_D(Freq_MHz, G2D(Freq_MHz, Gmainbeam,Efficiency), AzimuthOffset, Efficiency);
}

double cal_ant_G0_m(double Freq_MHz,double Gmainbeam,double AzimuthOffset,double Efficiency=0.6,int flag=1)
{
	if(flag==1)
		return Gfai_F1245_G0(Freq_MHz,Gmainbeam,AzimuthOffset,Efficiency);
	else
		return Gfai_F699_G0(Freq_MHz,Gmainbeam,AzimuthOffset,Efficiency);
}

double singleBeam_S672(double  phi,double G_m,double phi_0,double L_s)
{

   double a = 2.58;
   double b = 6.32;
   double G;
    if( L_s <= -25 && L_s > -30)
	{	a = 2.88;b = 6.32;}
    else if(L_s <= -30)
    {   a = 3.16;b = 6.32;	}
    
    double phi_1 = phi_0*pow(10,(G_m + L_s + 20)/25);
    
    if(phi < a*phi_0)
        G = G_m - 3*pow(phi/phi_0,2) ;
    else if(phi < b*phi_0)
        G = G_m + L_s;
    else if(phi < phi_1)
        G = G_m + L_s + 20 - 25*log10(phi/phi_0);
    else 
        G = 0;
    return G;

}
double cal_thita_S672(double dG,double G_m,double phi_0,double L_s)
{
	double th_min = 180,th_max = 0,thita,th;
	double Gmin = singleBeam_S672( th_min,G_m,phi_0,L_s);
	double Gmax = singleBeam_S672( th_max,G_m,phi_0,L_s);
	double G = 10000, G0 = G_m - dG;
	if( G0<Gmin)
	{
		thita = th_min;
		//warning('待计算dG过小');
	
		return thita;
	}
	if (G0>Gmax)
	{	thita = th_max;
		//warning('待计算dG过大');
		return thita;
	}

	while(abs(G - G0)>0.01)
	{
		th = 0.5*(th_max+th_min);
		G = singleBeam_S672( th,G_m,phi_0,L_s);
		if( G<=G0 && G<Gmax)
		{
			th_min = th;
			Gmin = G;
		}
		else if( G>G0 && G> Gmin)
		{
			th_max = th;
			Gmax = G;
		}
	}
	thita = th;
	return thita;
}

double Ghr(double xh,double kh,double G180)
{
	double Ghr_out;
	double lamda_kh = 3*(1-pow(0.5,-1*kh));
	if( xh<=0.5)
		Ghr_out = -12*xh*xh;
	else
		Ghr_out = -12*pow(xh,2-kh) - lamda_kh;
	
	if( Ghr_out< G180)
		Ghr_out = G180;
	return Ghr_out;
}

double Gvr(double xv,double kv,double th3,double ka,double kp,int flag_bc)
{
	double Gvr_out;
	double C,xk,G180,lkv;
	if (flag_bc==1)
    {
		C = 10*log10(pow(180/th3,1.5)*(pow(4.0,-1.5)+kv)/(1+8*kp))/log10(22.5/th3);
        xk = sqrt(1 - 0.36*kv);
        G180 = -12+10*log10(1+8*kp)-15*log10(180/th3);
        lkv = 12 - C*log10(4.0) -10*log10(pow(4.0,-1.5)+kv);
        if( xv<xk)
            Gvr_out = -12*xv*xv;
        else if (xv<4)
            Gvr_out = -12+10*log10(pow(xv,-1.5)+kv);
        else if( xv<90/th3)
            Gvr_out = -1*lkv - C*log10(xv);
        else if (xv == 90/th3)
            Gvr_out = G180;
	}
    else if(flag_bc == 2)
    {
		C = 10*log10(pow(180/th3,1.5)*(pow(4.0,-1.5)+kv)/(1+8*ka))/log10(22.5/th3);
        xk = sqrt(1.33 - 0.33*kv);
        G180 = -15+10*log10(1+8*ka)-15*log10(180/th3);
        lkv = 12 - C*log10(4.0) -10*log10(pow(4.0,-1.5)+kv);
        if(xv<xk)
            Gvr_out = -12*xv*xv;
        else if (xv<4)
            Gvr_out = -15+10*log10(pow(xv,-1.5)+kv);
        else if( xv<90/th3)
            Gvr_out = -1*lkv -3 -C*log10(xv);
        else if (xv == 90/th3)
            Gvr_out = G180;
	}
	return Gvr_out;
}
double cal_G_F1336(double G0,para_1336sec parameter,double freq)
{
	double fai = parameter.fai;
	double fai3 = parameter.fai3;
	double thita = parameter.thita;//单位是度
	double beta = parameter.beta;
	double ka = parameter.ka;
	double kp = parameter.kp;
	double kh = parameter.kh;
	double kv = parameter.kv;
	int flag = parameter.flag;//1表示峰值增益，2平均增益；前者用于定点后者用于集总干扰兼容分析计算
	int flag_tilt = parameter.flag_tilt;//1表示天线选择机械下倾方式；0表示选择电子下倾方式
	double thitah,faih,xh,xv,G180,R,G,alfa,fai_th,fai_3m,psai;
	double x,psai_alfa;
	if(flag_tilt)
	{
		thitah = pi*thita/180;faih = pi*fai/180;
		beta = beta*pi/180;
		thita = asin(sin(thitah)*cos(beta)+cos(thitah)*cos(faih)*sin(beta));
		fai = acos((-sin(thitah)*sin(beta)+cos(thitah)*cos(faih)*cos(beta))/cos(thita));
		thita = thita*180/pi;beta = beta*180/pi;fai = fai*180/pi;
	}
	else
	{
		if(thita+beta>=0)
			thita = (thita+beta)*90/(90+beta);
		else
			thita = (thita+beta)*90/(90-beta);
	}
	double th3 = 31000*pow(10.0,-0.1*G0)/fai3;
	if(freq<=6e9)
	{
		if(flag==2)
		{
			xh = abs(fai)/fai3;
		   //     ffai = acos(cos(fai*pi/180).*cos(thita*pi/180))*180/pi;
			G180 = -15+10*log10(1+8*ka)-15*log10(180/th3);
			xv = abs(thita)/th3;
			R = (Ghr(xh,kh,G180)-Ghr(180/fai3,kh,G180))/(Ghr(0,kh,G180)-Ghr(180/fai3,kh,G180));
			G = G0 + Ghr(xh,kh,G180) + R*Gvr(xv,kv,th3,ka,kp,flag);
		}
		else if(flag == 1)
		{
			xh = abs(fai)/fai3;
			G180 = -12+10*log10(1+8*kp)-15*log10(180/th3);
			xv = abs(thita)/th3;
			R = (Ghr(xh,kh,G180)-Ghr(180/fai3,kh,G180))/(Ghr(0,kh,G180)-Ghr(180/fai3,kh,G180));
			G = G0+Ghr(xh,kh,G180)+R*Gvr(xv,kv,th3,ka,kp,flag);
		}
	}
	else
	{
		alfa = atan(tan(thita*pi/180)/sin(fai*pi/180));
		fai_th = fai3;
		if( abs(fai)<=fai_th)
			fai_3m = fai3;
		else if( abs(fai)<=180)
			fai_3m = 1/(sqrt(pow(cos(90*(abs(fai)-fai_th)/(180-fai_th))/fai3,2)
			+pow(sin(90*(abs(fai)-fai_th)/(180-fai_th))/th3,2)));
		psai = 180*acos(cos(fai*pi/180)*cos(thita*pi/180))/pi;
		if( psai<=90 && psai>=0)
			psai_alfa = 1/(sqrt(pow(cos(alfa)/fai3,2)+pow(sin(alfa)/th3,2)));
		else if(psai<=180)
			psai_alfa = 1/(sqrt(pow(cos(thita*pi/180)/fai_3m,2)+pow(sin(thita*pi/180)/th3,2)));

		x = psai/psai_alfa;
		if(flag==1)
		{
			if( x<1 && x>=0)
				G = G0-12*x*x;
			else if( x>=1)
				G = G0-12-15*log10(x);
		}
		else if( flag == 2)
		{
			if( x<1.152 && x>=0)
				G = G0-12*x*x;
			else if (x>=1.152)
				G = G0-15-15*log10(x);
		}
	}
	return G;
}
double cal_G_omni(double G0,double thita,double beta,int flag,double freq)
{
	double k;
	if(abs(thita)>90)
    return -1000;
	
	if( freq<3*1e9)
		k = 0.7;
	else
		k=0;

	if(thita+beta>=0)
		thita = 90*(thita+beta)/(90+beta);
	else
		thita = 90*(thita+beta)/(90-beta);

	double thita3 = 107.6*pow(10,-0.1*G0),thita4,G,thita5;
	if(flag == 1)
	{
		thita4 = thita3 * sqrt(1-(1/1.2)*log10(k+1));
		if(abs(thita)<thita4)
			G = G0-12*pow(thita/thita3,2);
		else if (abs(thita)<thita3)
			G = G0 - 12 + 10*log10(k+1);
		else
			G = G0 -12 +10*log10(pow(abs(thita)/thita3,-1.5)+k);
	}
	else
	{
		thita5 = thita3*sqrt(1.25-(1/1.2)*log10(k+1));
		if(abs(thita)<thita3)
			G = G0-12*pow(thita/thita3,2);
		else if (abs(thita)<thita5)
			G = G0 - 15 + 10*log10(k+1);
		else
			G = G0 -15 +10*log10(pow(abs(thita)/thita3,-1.5)+k);
    
	}
	return G;
}
double cal_ant_G0_mainbeamwidth(int flag_f6991245,double freq,double Gmainbeam,double d_G,double *Gave,double Efficiency=0.6)
{
	int k = 0;double thita_t=-180.0,thita;
	double Gfai,Gave_out=0,thita_out1,thita_out;
	if(*Gave == 1)
	{
		for (thita_t = -180;thita_t<=180;thita_t = thita_t+0.001)
		{
			if(flag_f6991245!=0)
				Gfai = cal_ant_G0_m(freq,Gmainbeam,thita_t,Efficiency);
			else
				Gfai = cal_ant_G0_m(freq,Gmainbeam,thita_t,Efficiency,2);
			if(Gfai>= Gmainbeam-d_G)
			{
				if(k==0)
					thita_out1 = thita_t;
				Gave_out = Gave_out+pow(10.0,Gfai/10);
				thita_out = thita_t;
				k = k+1;
			}
		}
		*Gave = 10*log10( Gave_out/k);
		thita = thita_out-thita_out1;
	}
	else
	{
		double th_min = 180, th_max = 0,th;
		double Gmin,Gmax;

		if(flag_f6991245 != 0)
		{	Gmin = cal_ant_G0_m(freq,Gmainbeam,th_min,Efficiency);Gmax = cal_ant_G0_m(freq,Gmainbeam,th_max,Efficiency);}
		else
		{	Gmin = cal_ant_G0_m(freq,Gmainbeam,th_min,Efficiency,2);Gmax = cal_ant_G0_m(freq,Gmainbeam,th_max,Efficiency,2);}

		double G = 10000;
		double G0 = Gmainbeam - d_G;
		if(G0<Gmin)
			return  2*th_min;

		if(G0>Gmax)
			return 2*th_max;
    
		while(abs(G - G0)>0.01)
		{
			th = 0.5*(th_max+th_min);
			if(flag_f6991245 != 0)
				G = cal_ant_G0_m(freq,Gmainbeam,th,Efficiency);
			else
				G = cal_ant_G0_m(freq,Gmainbeam,th,Efficiency,2);

			if (G<=G0 && G<Gmax)
			{	th_min = th;
				Gmin = G;
			}
			else if( G>G0 && G> Gmin)
			{
				th_max = th;
				Gmax = G;
			}
		}
		thita = 2*th;
	}
	return thita;
}
double single_pattern(double phi,double theta,double G_Em,double Am,double SLA,double phi_3db,double theta_3db);

// dhspace is d_H/wave_length
// dvspace is d_V/wave_length
// 其他含义与建议书一致
double imt_5G_antenna(para_imt5Gant para)
{
	double phi, theta, G_Em, Am, SLA, phi_3db, theta_3db, theta_e,phi_e,dhspace,dvspace;
	int rownum,colnum;
	phi = para.phi,theta = para.theta,G_Em = para.G_Em,Am = para.Am;
	SLA = para.SLA;
	phi_3db = para.phi_3db;theta_3db = para.theta_3db;
	theta_e = para.theta_e;	phi_e = para.phi_e;
	dhspace = 0.5*300/para.freq;dvspace =  0.5*300/para.freq;
	rownum = para.rownum;colnum = para.colnum;
	int m,n;
	double v_r,v_i,w_r,w_i,A_E,A_beam,temp;
	double sum_temp_r,sum_temp_i;
	sum_temp_r = 0;
	sum_temp_i = 0;
	for( m = 1;m < rownum; m ++ )
		for( n =1; n < colnum; n ++)
			{
			temp = 2*PI*( (n-1)*dvspace*cos(theta/180*PI) + (m-1)*dhspace*sin(theta/180*PI)*sin(phi/180*PI) );
			v_r = cos(temp);
			v_i = sin(temp);
			temp = 2*PI*( (n-1)*dvspace*sin(theta_e/180*PI) - (m-1)*dhspace*cos(theta_e/180*PI)*sin(phi_e/180*PI) ); 
			w_r = 1/sqrt(double(colnum*rownum))*cos(temp);
			w_i = 1/sqrt(double(colnum*rownum))*sin(temp);
			sum_temp_r = sum_temp_r + v_r*w_r - v_i*w_i;
			sum_temp_i = sum_temp_i + v_r*w_i + v_i*w_r;
            }
	A_E = single_pattern(phi,theta,G_Em,Am,SLA, phi_3db,theta_3db);
	A_beam = A_E + 10*log10(sum_temp_r*sum_temp_r + sum_temp_i*sum_temp_i);
	return A_beam;
}

double single_pattern(double phi,double theta,double G_Em,double Am,double SLA,double phi_3db,double theta_3db){

	double A_EH,A_EV,A_E,temp;

	temp = 12*pow((phi/phi_3db),2);
	if( temp <= Am )
		A_EH = -temp;
	else 
		A_EH = -Am;
	
	temp =  12*pow( ((theta-90)/theta_3db),2);
	if( temp <= SLA )
		A_EV = -temp;
	else 
		A_EV = -SLA;

	temp = -(A_EH + A_EV);
	if( temp <= Am )
		A_E = G_Em - temp;
	else
		A_E = G_Em - Am;
	
	return A_E;

}

#endif