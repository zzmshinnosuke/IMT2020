#ifndef _OTHERCALC_H
#define _OTHERCALC_H

#include <iostream>
#include <math.h>
#define pi 3.141592653589793
#define PI 3.141592653589793
double min(double *xin,int len,int *minnum);
double max(double *xin,int len,int *maxnum);
double round(double x);
double sum(double *xin,int len);
int *sort(double *xin,int len,bool flag);
bool compare(double a,double b)
{return a>b;}
double round(double r)
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}        

double sum(double *xin,int len)
{
	double sum_out = 0;
	for(int i =0;i<len;i++)
		sum_out = sum_out+xin[i];
	return sum_out;
}
void swap(double *A,int i,int j)
{
	double temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}
void swapi(int *A,int i,int j)
{
	int temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}
int *sort(double *xin,int len,bool flag)
{
	int *sortnum = new int[len];
	for(int n = 0;n<len;n++)
		sortnum[n] = n;
	if(flag == true)
	{
		for(int j=0;j<len-1;j++)
		for(int i=0;i<len-1-j;i++)
			{
				if(xin[i]>xin[i+1])
				{
					swap(xin,i,i+1);
					swapi(sortnum,i,i+1);
				}
			}
	}
	else
	{
		for(int j=0;j<len-1;j++)
		for(int i=0;i<len-1-j;i++)
			{
				if(xin[i]<xin[i+1])
				{
					swap(xin,i,i+1);
					swapi(sortnum,i,i+1);
				}
			}
	}
	return sortnum;
}
void xy2lonlat(double pxy[2],double lonlat0[2],double lonlat[2])
{
	double x,y;

	x = pxy[0];
	y = pxy[1];
	int R = 6371;
	double lon0,lat0,lat,R_,lon;
	//double lonlat[2];
	lon0 = lonlat0[0]*pi/180;
	lat0 = lonlat0[1]*pi/180;

	lat = y/R + lat0;
	R_ = R*cos(lat);
	lon = x/R_ + lon0;
	lonlat[0] = lon*180/pi;
	lonlat[1] = lat*180/pi;
	//return lonlat;
}
void cal_thita(double d,double h1,double h2,double thita12[2])
{long R = 6371000;
//double thita12[2] = {0};
double tm = 0;
double tthita = 0;
        if(R+h1<=(R+h2)/cos(d/R))
		{
			tm = (R+h2)/(R+h1)-cos(d/R);
            tthita = sin(d/R)/tm;
            thita12[1] = atan(tthita);
            thita12[0] = pi - d/R - thita12[1];
		}
        else
        {
			tm = (R+h1)/(R+h2)-cos(d/R);
            tthita = sin(d/R)/tm;
            thita12[0] = atan(tthita);
            thita12[1] = pi - d/R - thita12[0];
		}
        if(d>=R*(acos(R/(R+h1))+acos(R/(R+h2))))
        {    thita12[0] = 180*acos(R/(R+h1))/pi;
		thita12[1] = 180*acos(R/(R+h2))/pi;}
        else
		{
            thita12[0] = 90-thita12[0]*180/pi;
            thita12[1] = 90-thita12[1]*180/pi;
		}
		//return thita12;
}

//  本函数用于计算地球上两点间的距离，输入为pA/pB（经纬坐标，单位度）,hA/B，单位米，输出dis，单位km
double cal_distance(double pA[2],double pB[2],double fai[2],double hA=0,double hB=0)
{
	double dis = 0;

	hA = hA/1000;hB = hB/1000;
	int R = 6371;
	double pA_long_lat[2],pB_long_lat[2];
	pA_long_lat[0] = pA[0]*pi/180;
	pA_long_lat[1] = pA[1]*pi/180;
	pB_long_lat[0] = pB[0]*pi/180;
	pB_long_lat[1] = pB[1]*pi/180;

	double radLat1 = pA_long_lat[1];
	double radLat2 = pB_long_lat[1];
	double a,b;
	b = pA_long_lat[0] - pB_long_lat[0];
	a = pA_long_lat[1] - pB_long_lat[1];
	double s = 2*asin(sqrt(sin(a/2)*sin(a/2)+cos(radLat1)*cos(radLat2)*sin(b/2)*sin(b/2)));
	dis = s*R;
	double d1,d2;
	double dis_t = R*(acos(R/(R+hA))+acos(R/(R+hB)));
	if (dis> dis_t) // 超出地理上的视距范围
		dis = sqrt(hA*hA+2*R*hA)+sqrt(hB*hB+2*R*hB)+dis-dis_t;
	else
	{
		d1 = R+hA;d2 = R+hB;
		dis = sqrt(d1*d1+d2*d2-2*d1*d2*cos(s));
	}

	// 计算两点间方位角
	double mLonA = pA_long_lat[0],mLonB = pB_long_lat[0];
	double mLatA = pA_long_lat[1],mLatB = pB_long_lat[1];
	double AngC = 0,Fai1 = 0,Fai2 = pi;
	double AngA = 0,AngB = 0,AngBeta = 0,TempX=0,DLon = 0,DLat = 0;
	if ((mLonA == mLonB) && (mLatA == mLatB))
	{
		AngC = 0;
		Fai1 = 0;
		Fai2 = pi;
	}
	else
	{
		AngA = pi * (90 - mLatB) / 180;
		AngB = pi * (90 - mLatA) / 180;
		AngBeta = pi * abs(mLonB - mLonA) / 180;
		TempX = cos(AngA) * cos(AngB) + sin(AngA) * sin(AngB) * cos(AngBeta);
		AngC = atan(sqrt(1 - TempX * TempX) / TempX);
		Fai1 = asin(sin(AngA) * sin(AngBeta) / sin(AngC));
		Fai2 = asin(sin(AngB) * sin(AngBeta) / sin(AngC));
    
		DLon = mLonB - mLonA;
		DLat = mLatB - mLatA;
		//根据象限调整角度值
		if (DLon >= 0)
			if (DLat >= 0)
			   //fai1不用调整            第一象限
				Fai2 = pi + Fai2;        //第三象限
			else
			{
				Fai1 = pi - Fai1;        //第二象限
				Fai2 = 2 * pi - Fai2 ;   //第四象限
			}
		else
			if (DLat >= 0)
			{
				Fai1 = 2 * pi - Fai1 ;   //第四象限
				Fai2 = pi - Fai2   ;     //第二象限
			}
			else
				Fai1 = pi + Fai1  ;      //第三象限
				//fai2不用调整            第一象限

	}
    

	if (Fai1>pi)
		Fai1 = Fai1-2*pi;

	if (Fai2>pi)
		Fai2 = Fai2-2*pi;
	fai[0] = Fai1*180/pi;
	fai[1] = Fai2*180/pi;
	return dis;
}
double cal_medium(double min_value,double max_value,double x_min,double x_max,double x)
{
	double mediumvalue=0;
	if(x_min == x_max && x_min == x)
		mediumvalue = 0.5*(min_value+max_value);
	else
		mediumvalue = min_value + (max_value-min_value)*(x-x_min)/(x_max-x_min);
	return mediumvalue;
}
double T(double y)
    {
		return sqrt(-2*log(y));     //(36c)
	}
double C(double z)
    {
		double C0 = 2.515516698;
        double C1 = 0.802853;
        double C2 = 0.010328;
        double D1 = 1.432788;
        double D2 = 0.189269;
        double D3 = 0.001308;
        return (((C2*T(z)+C1)*T(z))+C0)/(((D3*T(z)+D2)*T(z)+D1)*T(z)+1);//(36d)
	}
double Qi(double x)
{
	if( x<= 0.5)
		return T(x)-C(x);//(36a)
	else
		return -(T(1-x)-C(1-x)); //(36b)
}
//// 下述单位均为MHz
//// Pf表示干扰信号谱mask，格式为（20 -60;10 -40;……）
//// Rf表示受扰端选择性，格式为（20 -60;10 -40;……）
//// df表示两者频率间隔
//// inter_f表示拟合所需频率间隔，通常应为 各级带宽/2 的公约数
//// SP和SR分别表示干扰与受扰端的杂散发射电平(负数)，超出带宽标称值将以此值作为衰减值
double min(double *Pf,int len,int *min_i)
{
	int i = 0;
	double min_pf = Pf[0];
	for(i=0;i<len;i++) 	
		if( Pf[i]<=min_pf)
		{
			min_pf = Pf[i];*min_i = i;
		}
	return min_pf;
}
double min(double *Pf,int len)
{
	int i = 0;
	double min_pf = Pf[0];
	for(i=0;i<len;i++) 	
		if( Pf[i]<=min_pf)
		{
			min_pf = Pf[i];
		}
	return min_pf;
}
double max(double *Pf,int len)
{
	int i = 0;
	double max_pf = Pf[0];
	for(i=0;i<len;i++) 	
		if( Pf[i]>=max_pf)
		{ max_pf = Pf[i];}
	return max_pf;
}
double max(double *Pf,int len,int *max_i)
{
	int i = 0;
	double max_pf = Pf[0];
	for(i=0;i<len;i++) 	
		if( Pf[i]>=max_pf)
		{ max_pf = Pf[i];*max_i = i;}
	return max_pf;
}

///////////////////////////////////////////////////////////////////////
// find_num函数说明：x为输入数组，xlen为待寻找的长度，x_d为目标判断值
// flag为标志位，0表示x==x_d,1表示x>x_d,-1表示x<x_d,输出为从小到大第一个符合条件的位置。
int find_num(double x[],int xlen,double x_d,int flag = 0)
{
	int i=0;
	if (flag == 0)
	{
		for (;i<xlen;i++)
		{
			if(x[i] == x_d)
				return i;
		}
		if(i==xlen)
			return -1;
	}
	else if(flag >0)
	{
		for (;i<xlen;i++)
		{
			if(x[i] > x_d)
				return i;
		}
		if(i==xlen)
			return -1;
	}
	else
	{
		for (;i<xlen;i++)
		{
			if(x[i] < x_d)
				return i;
		}
		if(i==xlen)
			return -1;
	}


	return -1;
}
int find_num2(double x[],double y[],int len,double x_d,double y_d)
{
	int i=0;
	for(;i<len;i++)
	{
		if (x[i]==x_d && y[i]==y_d)
			return i;
	}
	return -1;
}
double cal_FDR_simple(double Bt,double Br)
{
	if (Bt<=Br)
		return 0;
	else
		return 10*log10(Bt/Br);
}		

double linear_point(double x1,double y1,double x2,double y2,double x)
{
	double y;
	if(x == x1) return y1;
	if(x == x2) return y2;
	y = y1+(x-x1)*(y2-y1)/(x2-x1);
	return y;
}
double mod(double x,double y)
{
	double n ;
	if(y==0)
		return x;
	else
	{
		n = floor(x/y);
		return x - n*y ;
	}
}
double csc(double x)
{
	return 1.0/sin(x);
}
double cal_FDR_inter(double *Pf,double *freq_P,double *Rf,double *freq_R,double df,double inter_f,double SP,double SR,int len_T,int len_R)
{
	double *Pf_t = new double[len_T+2];
	double *freq_Pt = new double[len_T+2];
	double *Rf_t = new double[len_T+2];
	double *freq_Rt = new double[len_T+2];

	for(int i=0;i<len_R;i++)
	{
		freq_Rt[i] = freq_R[i]+df;
		Rf_t[i] = Rf[i];
	}
	for(int i=0;i<len_T;i++)
	{
		freq_Pt[i] = freq_P[i];
		Pf_t[i] = Pf[i];
	}
	int *tempminmax = new int;
	double freq_minR = min(freq_R,len_R),freq_minP = min(freq_P,len_T);
	double freq_maxR = max(freq_R,len_R),freq_maxP = max(freq_P,len_T);
	double freq_min = freq_minR<freq_minP?freq_minR:freq_minP;
	double freq_max = freq_maxR<freq_maxP?freq_maxP:freq_maxR;
	if(freq_min<freq_minR)
	{
		freq_Rt[len_R] = freq_min;
		Rf_t[len_R] = SR;
		len_R++;
	}
	if( freq_min<freq_minP)
	{
		freq_Pt[len_T] = freq_min;
		Pf_t[len_T] = SP;
		//Pf = [Pf;freq_min SP];
		len_T++;
	}
	if (freq_max>freq_maxR)
	{
		freq_Rt[len_R] = freq_max;
		Rf_t[len_R] = SR;
		len_R++;
		//Rf = [Rf;freq_max SR];
	}
	if (freq_max>freq_maxP)
	{
		freq_Pt[len_T] = freq_max;
		Pf_t[len_T] = SP;
		len_T++;
		//Pf = [Pf;freq_max SP];
	}	

	//[len_Pf,tmp]=size(Pf);
	//[len_Rf,tmp] = size(Rf);
	int *tmp_T = new int[len_T];
	int *tmp_R = new int[len_R];
	double *Pf_d = new double[len_T];
	double *Rf_d = new double[len_R];
	
	tmp_T = sort(freq_Pt,len_T,true);
	for(int i =0 ;i<len_T;i++)
	{
		Pf_d[i] = Pf_t[tmp_T[i]];
	}

	tmp_R = sort(freq_Rt,len_R,true);
	for(int i =0 ;i<len_R;i++)
	{
		Rf_d[i] = Rf_t[tmp_R[i]];
	}

	// Pf_dr = 10.^(Pf_d/10);
	int len_freqmax = len_T+len_R+int(floor((freq_max-freq_min)/inter_f))+1;
	double *freq_squ = new double[len_freqmax];
	double freq_current = freq_min;
	int freqi = 0;//频率个数总数
	int Pi=0,Ri = 0;
	double *P_D = new double[len_freqmax];
	double *R_D = new double[len_freqmax];
	while(freq_current<freq_max)
	{
		if(freq_current>freq_Pt[Pi])
		{
			if(freq_current>freq_Rt[Ri])
			{
				if(freq_Pt[Pi]>freq_Rt[Ri])
				{
					freq_squ[freqi] = freq_Rt[Ri];
					Ri++;freqi++;
					freq_squ[freqi] = freq_Pt[Pi];
					Pi++;freqi++;
				}
				else if(freq_Pt[Pi]<freq_Rt[Ri])
				{
					freq_squ[freqi] = freq_Pt[Pi];
					Pi++;freqi++;
					freq_squ[freqi] = freq_Rt[Ri];
					Ri++;freqi++;
				}
				else
				{
					freq_squ[freqi]= freq_Pt[Pi];
					Pi++;Ri++;freqi++;
				}
			}
			else if (freq_current==freq_Rt[Ri])
			{
				freq_squ[freqi] = freq_Pt[Pi];
				Pi++;freqi++;
				freq_squ[freqi] = freq_current;
				Ri++;
			}
			else
			{
				freq_squ[freqi] = freq_Pt[Pi];
				Pi++;freqi++;
			}
		}
		else if(freq_current==freq_Pt[Pi])
		{
			Pi++;
			if(freq_current>freq_Rt[Ri])
			{
				freq_squ[freqi] = freq_Rt[Ri];
				Ri++;freqi++;
			}
			else if(freq_current==freq_Rt[Ri])
				Ri++;
		}
		else
		{
			if(freq_current>freq_Rt[Ri])
			{
				freq_squ[freqi] = freq_Rt[Ri];
				Ri++;freqi++;
			}
			else if(freq_current == freq_Rt[Ri])
			{
				freq_squ[freqi] = freq_Rt[Ri];
				Ri++;
			}
		}
		freq_squ[freqi] = freq_current;
		freq_current = freq_current+inter_f;
		freqi++;
	}
	if(freq_squ[freqi-1]<freq_max)
	{	freq_squ[freqi] = freq_max;freqi++;}
	Pi = 0;Ri = 0;
	int P1 = 0,P2 = 1;
	int R1 = 0,R2 = 1;
	for(int i = 0;i<freqi;i++)
	{
		if(freq_squ[i]>=freq_Pt[P2])
		{
			P1++;P2++;
			if(P2==len_T)
				P2 = len_T-1;
		}
		P_D[Pi] = linear_point(freq_Pt[P1],Pf_d[P1],freq_Pt[P2],Pf_d[P2],freq_squ[i]);
		Pi++;
		if(freq_squ[i]>=freq_Rt[R2])
		{
			R1++;R2++;
			if(R2==len_R)
				R2 = len_R-1;
		}
		R_D[Ri] = linear_point(freq_Rt[R1],Rf_d[R1],freq_Rt[R2],Rf_d[R2],freq_squ[i]);
		Ri++;
	}
	double A = 0;
	double pdi,pdi1;
	double rdi,rdi1;
	double Btmp=0,B=0;
	for(int i = 0;i<freqi-1;i++)
	{
		pdi = pow(10.0,P_D[i]/10.0);
		pdi1 = pow(10.0,P_D[i+1]/10.0);
		A = A+ 0.5*(pdi+pdi1)*(freq_squ[i+1] - freq_squ[i]);
		//Btmp = P_D[i]*R_D[i];
		rdi = pow(10.0,R_D[i]/10.0);
		rdi1 = pow(10.0,R_D[i+1]/10.0);
		B = B+0.5*(pdi*rdi+pdi1*rdi1)*(freq_squ[i+1] - freq_squ[i]);
	}
	double FDR = 10*log10(A/B);
	return FDR;
}
#endif