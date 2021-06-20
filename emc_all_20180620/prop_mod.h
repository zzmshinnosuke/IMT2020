#ifndef _PROP_H
#define _PROP_H
#include <math.h>
#include "other_calc_mod.h"
#include "prop_mod_staticdata.h"
double cal_Lp_free(double freq,double distance)
{
	return 32.4+20*log10(freq*distance);
}
double F(double x)
{
	double out = 0.661*x+0.339*sqrt(pow(x,2)+5.51);
	return 1/out;
}
double cal_gammao_p676_11(double f,double p,double T,double rho)
{
	double thita = 300/T;
	double fi,a1,a2,a3,a4,a5,a6,e,df,deta,Si,Fi,d,NDf,Nf_oxy;
	double sum_NF=0;
	e = rho*T/216.7;
	for(int i=0;i<44;i++)
	{
		fi = para_o[i][0];
		a1 = para_o[i][1];
		a2 = para_o[i][2];
		a3 = para_o[i][3];
		a4 = para_o[i][4];
		a5 = para_o[i][5];
		a6 = para_o[i][6];
		df = a3*0.0001*(p*pow(thita,(0.8-a4))+1.1*e*thita);
		deta = (a5+a6*thita)*0.0001*(p+e)*pow(thita,0.8);
		Si = a1*pow(10.0,-7)*p*pow(thita,3)*exp(a2*(1-thita));
		Fi = f/fi*((df-deta*(fi-f))/(pow(fi-f,2)+df*df)+(df-deta*(fi+f))/(pow(fi+f,2)+df*df));
		sum_NF = sum_NF+Si*Fi;
	}
	d = 5.6*0.0001*(p+e)*pow(thita,0.8);
	NDf = f*p*thita*thita*((6.14*0.00001/d/(1+pow(f/d,2)))+(1.4*pow(10.0,-12)*p*pow(thita,1.5)/(1+1.9*0.00001*pow(f,1.5))));
	Nf_oxy = sum_NF+NDf;

	return 0.182*f*Nf_oxy;
}
double cal_gammaw_p676_11(double f,double p,double T,double rho)
{
	double thita = 300/T;
	double fi,b1,b2,b3,b4,b5,b6,e,df,Si,Fi;
	double sum_NF=0;
	for(int i=0;i<35;i++)
	{
		fi = para_w[i][0];
		b1 = para_w[i][1];
		b2 = para_w[i][2];
		b3 = para_w[i][3];
		b4 = para_w[i][4];
		b5 = para_w[i][5];
		b6 = para_w[i][6];
		e = rho*T/216.7;
		df = b3*0.0001*(p*pow(thita,b4)+b5*e*pow(thita,b6));


		Si = b1*0.1*e*pow(thita,3.5)*exp(b2*(1-thita));

		Fi = f/fi*(df/(pow(fi-f,2)+df*df)+df/(pow(fi+f,2)+df*df));

		sum_NF = sum_NF+Si*Fi;
	}

	return 0.182*f*sum_NF;
}
double cal_Ag_slant(double h1,double h2,double f,double dtot,double T = 288)
{
	double p = 1013.25;
	double h11 = h1<h2?h1:h2;
	double h20 = h2;
	h2 = h1>h2?h1:h2;
	h1 = h11;
	if (h2>20)
		h2 = 20;

	double rho = 7.5 ;
	double rho1 = rho*exp(-h2/2);
	double ee = rho1*T/216.7;
	double ptot = p*pow((1-0.02257*h2),5.256)+ee;

	double gammao = cal_gammao_p676_11(f,p,T,rho);
	double gammaw = cal_gammaw_p676_11(f,p,T,rho);

	double rp = ptot/1013.25;
	double t1 = 4.64*exp(-1*pow((f-59.7)/(2.87+12.4*exp(-7.9*rp)),2))/(1+0.066*pow(rp,-2.3));
	double t2 = 0.14*exp(2.12*rp)/(pow(f-118.75,2)+0.031*exp(2.2*rp));
	double t3 = 0.0114*f*(-0.0247+0.0001*f+1.61*pow(10.0,-6)*f*f)/(1-0.0169*f+
		4.1*0.00001*f*f+3.2*pow(10.0,-7)*pow(f,3))/(1+0.14*pow(rp,-2.6));

	double ho = 6.1*(1+t1+t2+t3)/(1+0.17*(pow(rp,-1.1)));
	ho = ho<10.7*pow(rp,0.3)?ho:10.7*pow(rp,0.3);
	double sigmaw = 1.013/(1+exp(-8.6*(rp-0.57)));
	double hw = 1.66*(1+1.39*sigmaw/(pow(f-22.235,2)+2.56*sigmaw)+3.37*sigmaw/(pow(f-183.31,2)+
		4.69*sigmaw)+1.58*sigmaw/(pow(f-325.1,2)+2.89*sigmaw));
	double Ao = gammao*ho;double Aw = gammaw*hw;
	double A = Ao+Aw;
	double thita1,thita2;
	double thita12[2]={0};
	cal_thita(dtot*1e3,h1*1e3,h20*1e3,thita12);
	thita1 = thita12[0];
	thita2 = thita12[1];

	thita1 = -thita1;
	double Re,x1,x11,x2,x22,tmp1,tmp2,tmp3,tmp4,Ag,hoo,hww;
	if(thita1>5)
	{
		hoo = ho*(exp(-1*h1/ho)-exp(-1*h2/ho));
		hww = hw*(exp(-1*h1/hw)-exp(-1*h2/hw));
		Ag = (gammao*hoo+gammaw*hww)/sin(thita1*pi/180);
	}
	else
	{
		if(thita1<0)
			thita1 = 0;

		Re = 8500;
		x1 = tan(thita1*pi/180)*sqrt((Re+h1)/ho);
		x11 = tan(thita1*pi/180)*sqrt((Re+h1)/hw);
		thita2 = acos((Re+h1)*cos(thita1*pi/180)/(Re+h2));
		x2 = tan(thita2*pi/180)*sqrt((Re+h2)/ho);
		x22 = tan(thita2*pi/180)*sqrt((Re+h2)/hw);
		tmp1 = sqrt(Re+h1)*F(x1)*exp(-1*h1/ho)/cos(thita1*pi/180);
		tmp2 = sqrt(Re+h2)*F(x2)*exp(-1*h2/ho)/cos(thita2*pi/180);
		tmp3 = sqrt(Re+h1)*F(x11)*exp(-1*h1/hw)/cos(thita1*pi/180);
		tmp4 = sqrt(Re+h2)*F(x22)*exp(-1*h2/hw)/cos(thita2*pi/180);
		Ag = gammao*sqrt(ho)*(tmp1-tmp2)+gammaw*sqrt(hw)*(tmp3-tmp4);
	}
	return Ag;
}
double cal_Lp_E2A_5G(double freq,double D,double p_l,double h1,double h2)
{
	double f = freq/1000;
	double lonlat0[2] = {0};
	double lonlat1[2] = {0};
	double p[2] = {D};
	double fai[2]={0};
	xy2lonlat(p,lonlat0,lonlat1);
	double D_r = D;
	double thita12[2] ={0};
double D_tr = cal_distance(lonlat0,lonlat1,fai,h1,h2);
cal_thita(D_r*1e3,h1,h2,thita12);
double thita1 = thita12[0],thita2 = thita12[1];
double Lp_free = cal_Lp_free(freq,D_tr);
double Lp_Ag = cal_Ag_slant(h1/1000,h2/1000,freq/1000,D,15+273);
double K1 = 93.0*pow(f,0.175),A1 = 0.05;
double thita = -thita1;
double Lces = pow(-K1*log(1-p_l/100)*1/(tan(A1*(1-thita/90)+pi*thita/180)),(0.5*(90-thita)/90))-1-0.6*Qi(p_l/100);
return Lp_free+Lp_Ag+Lces;

}
double cal_Lp_p528_rec(double Dis_km,double Freq_MHz,double H1,double H2,double Per)
{
	double t[5] = {1, 5, 10, 50, 95};
	double f[8] = {125, 300, 600, 1200, 2400, 5100, 9400, 15500};
	double D[1801] = {0};
	int t_len = 5,f_len = 8,D_len = 1801;
	int t_num,f_num,d_num;
	double f_inf,f_sup,t_inf,t_sup,d_inf,d_sup;
	int i= 1;
	for( i=0;i<=1800;i++)
		D[i] = i;
	int t_tmp,f_tmp,d_tmp;
	
	int t_inf_num,t_sup_num,f_inf_num,f_sup_num,d_inf_num,d_sup_num;
	//// step 1，寻找两个标称百分比
	if(Per<1 || Per>95)
		return -1;
	t_tmp = find_num(t,5,Per);
	if(t_tmp!=-1)
	{
		t_inf = Per;t_sup = Per;
		t_inf_num = t_tmp;t_sup_num = t_tmp;
	}
	else
	{
		t_num = find_num(t,5,Per,1);
		t_inf = t[t_num-1];t_sup = t[t_num];
		t_inf_num = t_num-1;t_sup_num = t_num;
	}
	//// step2，寻找两个标称频率
	
	if(Freq_MHz<125 || Freq_MHz>15500)
		return -1;
	f_tmp = find_num(f,8,Freq_MHz);
	if (f_tmp!=-1)
	{
		f_inf = Freq_MHz;f_sup = Freq_MHz;
		f_inf_num = f_tmp;f_sup_num = f_tmp;
	}
	else
	{
		f_num = find_num(f,8,Freq_MHz,1);
		f_inf = f[f_num-1];f_sup = f[f_num];
		f_inf_num = f_num-1;f_sup_num = f_num;
	}

	//// step3，寻找两个标称距离
	
	if(Dis_km<0 || Dis_km>1800)
		return -1;
	d_tmp = find_num(D,1801,Dis_km);
	if (d_tmp!=-1)
	{
		d_inf = Dis_km;d_sup = Dis_km;
		d_inf_num = d_tmp;d_sup_num = d_tmp;
	}
	else
	{
		d_num = find_num(D,1801,Dis_km,1);
		d_inf = D[d_num-1];d_sup = D[d_num];
		d_inf_num = d_num-1;d_sup_num = d_num;
	}

	//// step6，获得H1和H2对应的传输损耗
	//// step6.1获得H2的高低端标称值
	double h21_2[18] = {1000, 1000, 1000, 1000, 1000, 10000, 10000,
		10000, 10000, 10000, 10000, 20000, 20000, 20000, 20000, 20000, 20000, 20000};
	double h21_1[18] = {1.5, 15, 30, 60, 1000, 1.5, 15, 30, 60, 1000, 10000, 1.5, 15, 30, 60, 1000, 10000, 20000};
	double h2[3] = {1000, 10000, 20000};
	double h1[7] = {1.5, 15, 30, 60, 1000, 10000, 20000};
	double h1_inf,h1_sup,h2_inf,h2_sup;
	int h1_num,h2_num;
	if (H2<H1)
		return -1;
	if (H2<1000 || H2>20000)
		return -1;
	if (find_num(h2,3,H2)!=-1)
	{	h2_inf = H2;h2_sup = H2;}
	else
	{
		h2_num = find_num(h2,3,H2,1);
		h2_inf = h2[h2_num-1];h2_sup = h2[h2_num];
	}

	//// step6.2,获得H1的高低端标称值
	if (H1<1.5 || H1>20000)
		return -1;
	if (find_num(h1,7,H1)!=-1)
	{	h1_inf = H1;h1_sup = H1;}
	else
	{
		h1_num = find_num(h1,7,H1,1);
		h1_inf = h1[h1_num-1];h1_sup = h1[h1_num];
	}

	//// step6.6,获得dinf与dsup分别对应h1_inf,h2_inf的传播损耗
	int hnum;
	double min_value,max_value,x_min,x_max,x;
	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_inf);

	t_tmp = find_num(t,t_len,t_inf);f_tmp = find_num(f,f_len,f_inf);d_tmp = find_num(D,D_len,d_inf);
	double Lp_tinf_finf_dinf_h1inf_h2inf = data528[t_inf_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tinf_finf_dsup_h1inf_h2inf = data528[t_inf_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_finf_dinf_h1inf_h2inf;
	max_value = Lp_tinf_finf_dsup_h1inf_h2inf;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tinf_finf_d_h1inf_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_inf);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_inf));
	double Lp_tinf_finf_dinf_h1sup_h2inf = data528[t_inf_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tinf_finf_dsup_h1sup_h2inf = data528[t_inf_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_finf_dinf_h1sup_h2inf;
	max_value = Lp_tinf_finf_dsup_h1sup_h2inf;
	double Lp_tinf_finf_d_h1sup_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_finf_d_h1inf_h2inf;
	max_value = Lp_tinf_finf_d_h1sup_h2inf;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tinf_finf_d_h1_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_sup);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_sup));
	double Lp_tinf_finf_dinf_h1inf_h2sup = data528[t_inf_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tinf_finf_dsup_h1inf_h2sup = data528[t_inf_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_finf_dinf_h1inf_h2sup;
	max_value = Lp_tinf_finf_dsup_h1inf_h2sup;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tinf_finf_d_h1inf_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_sup);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_sup));
	double Lp_tinf_finf_dinf_h1sup_h2sup = data528[t_inf_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tinf_finf_dsup_h1sup_h2sup = data528[t_inf_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_finf_dinf_h1sup_h2sup;
	max_value = Lp_tinf_finf_dsup_h1sup_h2sup;
	double Lp_tinf_finf_d_h1sup_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_finf_d_h1inf_h2sup;
	max_value = Lp_tinf_finf_d_h1sup_h2sup;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tinf_finf_d_h1_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_finf_d_h1_h2inf;
	max_value = Lp_tinf_finf_d_h1_h2sup;
	x_min = log(h2_inf);x_max = log(h2_sup);x = log(H2);
	double Lp_tinf_finf_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_inf);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_inf));
	double Lp_tinf_fsup_dinf_h1inf_h2inf = data528[t_inf_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tinf_fsup_dsup_h1inf_h2inf = data528[t_inf_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_fsup_dinf_h1inf_h2inf;
	max_value = Lp_tinf_fsup_dsup_h1inf_h2inf;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tinf_fsup_d_h1inf_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_inf);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_inf));
	double Lp_tinf_fsup_dinf_h1sup_h2inf = data528[t_inf_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tinf_fsup_dsup_h1sup_h2inf = data528[t_inf_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_fsup_dinf_h1sup_h2inf;
	max_value = Lp_tinf_fsup_dsup_h1sup_h2inf;
	double Lp_tinf_fsup_d_h1sup_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_fsup_d_h1inf_h2inf;
	max_value = Lp_tinf_fsup_d_h1sup_h2inf;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tinf_fsup_d_h1_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_sup);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_sup));
	double Lp_tinf_fsup_dinf_h1inf_h2sup = data528[t_inf_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tinf_fsup_dsup_h1inf_h2sup = data528[t_inf_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_fsup_dinf_h1inf_h2sup;
	max_value = Lp_tinf_fsup_dsup_h1inf_h2sup;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tinf_fsup_d_h1inf_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_sup);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_sup));
	double Lp_tinf_fsup_dinf_h1sup_h2sup = data528[t_inf_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tinf_fsup_dsup_h1sup_h2sup = data528[t_inf_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tinf_fsup_dinf_h1sup_h2sup;
	max_value = Lp_tinf_fsup_dsup_h1sup_h2sup;
	double Lp_tinf_fsup_d_h1sup_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_fsup_d_h1inf_h2sup;
	max_value = Lp_tinf_fsup_d_h1sup_h2sup;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tinf_fsup_d_h1_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_fsup_d_h1_h2inf;
	max_value = Lp_tinf_fsup_d_h1_h2sup;
	x_min = log(h2_inf);x_max = log(h2_sup);x = log(H2);
	double Lp_tinf_fsup_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tinf_finf_d_h1_h2;
	max_value = Lp_tinf_fsup_d_h1_h2;
	x_min = log(f_inf);x_max = log(f_sup);x = log(Freq_MHz);
	double Lp_tinf_f_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_inf);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_inf));
	double Lp_tsup_finf_dinf_h1inf_h2inf = data528[t_sup_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tsup_finf_dsup_h1inf_h2inf = data528[t_sup_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_finf_dinf_h1inf_h2inf;
	max_value = Lp_tsup_finf_dsup_h1inf_h2inf;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tsup_finf_d_h1inf_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_inf);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_inf));
	double Lp_tsup_finf_dinf_h1sup_h2inf = data528[t_sup_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tsup_finf_dsup_h1sup_h2inf = data528[t_sup_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_finf_dinf_h1sup_h2inf;
	max_value = Lp_tsup_finf_dsup_h1sup_h2inf;
	double Lp_tsup_finf_d_h1sup_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_finf_d_h1inf_h2inf;
	max_value = Lp_tsup_finf_d_h1sup_h2inf;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tsup_finf_d_h1_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_sup);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_sup));
	double Lp_tsup_finf_dinf_h1inf_h2sup = data528[t_sup_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tsup_finf_dsup_h1inf_h2sup = data528[t_sup_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_finf_dinf_h1inf_h2sup;
	max_value = Lp_tsup_finf_dsup_h1inf_h2sup;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tsup_finf_d_h1inf_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_sup);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_sup));
	double Lp_tsup_finf_dinf_h1sup_h2sup = data528[t_sup_num][f_inf_num][d_inf_num][hnum+2];
	double Lp_tsup_finf_dsup_h1sup_h2sup = data528[t_sup_num][f_inf_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_finf_dinf_h1sup_h2sup;
	max_value = Lp_tsup_finf_dsup_h1sup_h2sup;
	double Lp_tsup_finf_d_h1sup_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_finf_d_h1inf_h2sup;
	max_value = Lp_tsup_finf_d_h1sup_h2sup;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tsup_finf_d_h1_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_finf_d_h1_h2inf;
	max_value = Lp_tsup_finf_d_h1_h2sup;
	x_min = log(h2_inf);x_max = log(h2_sup);x = log(H2);
	double Lp_tsup_finf_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_inf);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_inf));
	double Lp_tsup_fsup_dinf_h1inf_h2inf = data528[t_sup_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tsup_fsup_dsup_h1inf_h2inf = data528[t_sup_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_fsup_dinf_h1inf_h2inf;
	max_value = Lp_tsup_fsup_dsup_h1inf_h2inf;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tsup_fsup_d_h1inf_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_inf);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_inf));
	double Lp_tsup_fsup_dinf_h1sup_h2inf = data528[t_sup_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tsup_fsup_dsup_h1sup_h2inf = data528[t_sup_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_fsup_dinf_h1sup_h2inf;
	max_value = Lp_tsup_fsup_dsup_h1sup_h2inf;
	double Lp_tsup_fsup_d_h1sup_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_fsup_d_h1inf_h2inf;
	max_value = Lp_tsup_fsup_d_h1sup_h2inf;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tsup_fsup_d_h1_h2inf = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_inf,h2_sup);
	//hnum = find((h21(1,:)==h1_inf) & (h21(2,:)==h2_sup));
	double Lp_tsup_fsup_dinf_h1inf_h2sup = data528[t_sup_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tsup_fsup_dsup_h1inf_h2sup = data528[t_sup_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_fsup_dinf_h1inf_h2sup;
	max_value = Lp_tsup_fsup_dsup_h1inf_h2sup;
	x_min = log(d_inf);x_max = log(d_sup);x = log(Dis_km);
	double Lp_tsup_fsup_d_h1inf_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	hnum = find_num2(h21_1,h21_2,18,h1_sup,h2_sup);
	//hnum = find((h21(1,:)==h1_sup) & (h21(2,:)==h2_sup));
	double Lp_tsup_fsup_dinf_h1sup_h2sup = data528[t_sup_num][f_sup_num][d_inf_num][hnum+2];
	double Lp_tsup_fsup_dsup_h1sup_h2sup = data528[t_sup_num][f_sup_num][d_sup_num][hnum+2];
	min_value = Lp_tsup_fsup_dinf_h1sup_h2sup;
	max_value = Lp_tsup_fsup_dsup_h1sup_h2sup;
	double Lp_tsup_fsup_d_h1sup_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_fsup_d_h1inf_h2sup;
	max_value = Lp_tsup_fsup_d_h1sup_h2sup;
	x_min = log(h1_inf);x_max = log(h1_sup);x = log(H1);
	double Lp_tsup_fsup_d_h1_h2sup = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_fsup_d_h1_h2inf;
	max_value = Lp_tsup_fsup_d_h1_h2sup;
	x_min = log(h2_inf);x_max = log(h2_sup);x = log(H2);
	double Lp_tsup_fsup_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);

	min_value = Lp_tsup_finf_d_h1_h2;
	max_value = Lp_tsup_fsup_d_h1_h2;
	x_min = log(f_inf);x_max = log(f_sup);x = log(Freq_MHz);
	double Lp_tsup_f_d_h1_h2 = cal_medium(min_value,max_value,x_min,x_max,x);
	double Linf = Lp_tinf_f_d_h1_h2;
	double Lsup = Lp_tsup_f_d_h1_h2;
	double Qt = Qi(Per/100),Qinf = Qi(t_inf/100),Qsup = Qi(t_sup/100);
	double Lp_t_f_d_h1_h2;
	if(t_sup == t_inf)
		Lp_t_f_d_h1_h2 = Linf;
	else
		Lp_t_f_d_h1_h2 = Lsup*(Qinf - Qt)/(Qinf - Qsup) + Linf*(Qt - Qsup)/(Qinf - Qsup);
	
	return Lp_t_f_d_h1_h2;
}
double cal_dis_p528_rec(double Lp,double Freq,double H1,double H2,double Per)
{
	if (Freq>=15500 || Freq<=125)
		return -1;
	if (H1>H2)
		return -1;
	if (Per>95 || Per<1)
		return -1;
	int k = 1;
	double dis_min = 1,dis_max = 1000;
	double Lp_t= -100;
	double Lp_min = cal_Lp_p528_rec( dis_min,Freq, H1, H2, Per);
	double Lp_max = cal_Lp_p528_rec(dis_max,Freq,  H1, H2, Per);
	if(Lp<Lp_min)
		return dis_min;
	if(Lp>Lp_max)
		return dis_max;
	double dis_t;
	
	while (abs(Lp - Lp_t)>0.01)
	{
		dis_t = (dis_min+dis_max)/2;
		Lp_min = cal_Lp_p528_rec( dis_min,Freq, H1, H2, Per);
		Lp_max = cal_Lp_p528_rec(dis_max,Freq,  H1, H2, Per);
		Lp_t = cal_Lp_p528_rec(dis_t,Freq,  H1, H2, Per);
		if(Lp_t<Lp)
			dis_min = dis_t;

		if (Lp_t>Lp)
			dis_max = dis_t;
	}
	return dis_t;
}
double cal_dis_p525(double Lp,double freq)
{
	return pow( 10,((Lp-32.4)/20))/freq;
}

#endif