#ifndef _P452_H
#define _P452_H
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "other_calc_mod.h"
#include "prop_mod_staticdata.h"
using namespace std;
void HDR_REF(string Path,double AllAboutMap[33][12])
{
// global AllAboutMap;

	string HDR[33]=	{"W180N90.HDR",	"W140N90.HDR",	"W100N90.HDR",	"W060N90.HDR",
			"W020N90.HDR",	"E020N90.HDR",	"E060N90.HDR",	"E100N90.HDR",
			"E140N90.HDR",	"W180N40.HDR",	"W140N40.HDR",	"W100N40.HDR",
			"W060N40.HDR",	"W020N40.HDR",	"E020N40.HDR",	"E060N40.HDR",
			"E100N40.HDR",	"E140N40.HDR",	"W180S10.HDR",	"W140S10.HDR",
			"W100S10.HDR",	"W060S10.HDR",	"W020S10.HDR",	"E020S10.HDR",
			"E060S10.HDR",	"E100S10.HDR",	"E140S10.HDR",	"W180S60.HDR",
			"W120S60.HDR",	"W060S60.HDR",	"W000S60.HDR",	"E060S60.HDR",
			"E120S60.HDR"};

	// Path='..\Map\';
	string HDRfile[33];
	double R=6371;
	char temp[50]={0};
	FILE *f1 = NULL;
	const char *filename[33];
	// for part=1:1
	int part = 0;
	//float NROWS,NCOLS,NBANDS,NBITS,BANDROWBYTES,TOTALROWBYTES,BANDGAPBYTES,NODATA,ULXMAP,ULYMAP,XDIM,YDIM;
	for (part =0;part<33;part++)
	{
		HDRfile[part] = Path;
		HDRfile[part].append(HDR[part]);
		filename[part] = HDRfile[part].data(); 
		fopen_s(&f1,filename[part],"r");
		for (int i=0;i<41;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][0]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][1]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",& AllAboutMap[part][2]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][3]);
		for (int i=0;i<22;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",& AllAboutMap[part][4]);
		for (int i=0;i<22;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][5]);
		for (int i=0;i<22;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][6]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][7]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][8]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][9]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][10]);
		for (int i=0;i<15;i++)
			fscanf_s(f1,"%c",&temp[i]);
		fscanf_s(f1,"%lf",&AllAboutMap[part][11]);

		fclose(f1);

		//AllAboutMap[part][0]={NROWS,NCOLS,NBANDS,NBITS,BANDROWBYTES,TOTALROWBYTES,BANDGAPBYTES,NODATA,ULXMAP,ULYMAP,XDIM,YDIM};
	
	}
	return;
}
void pol2box(double p[2],double b[3]);
void box2pol(double bp[3],double pp[2]);
void DividePath(double p1[2],double p2[2],int n,double *pd1,double *pd2)
{
	double thita1=p1[0],	
	phai1=p1[1];

	double b1[3]={0},b2[3] = {0};
	pol2box(p1,b1);
	pol2box(p2,b2);

	double r1[3]={-1*cos(phai1+pi/2),-1*sin(phai1+pi/2),0};
	double r2[3];
	double tmp[2]={thita1+pi/2,phai1};
	pol2box(tmp,r2);
	double r3[3]={b1[0],b1[1],b1[2]};
	//r=[r1 r2 r3];

	double b2p[3]={0};
	double p2p[2] = {0};
	b2p[0] = r1[0]*b2[0]+r1[1]*b2[1]+r1[2]*b2[2];
	b2p[1] = r2[0]*b2[0]+r2[1]*b2[1]+r2[2]*b2[2];
	b2p[2] = r3[0]*b2[0]+r3[1]*b2[1]+r3[2]*b2[2];
	box2pol(b2p,p2p);
	double pdp[2]={0};
	double bdp[3]={0};
	double bd[3]={0};
	double pd[2]={0};
	for(int i=0;i<n;i++)
	{
		pdp[0] = (i+1.0)/(n+1)*p2p[0];
		pdp[1] = p2p[1];
		pol2box(pdp,bdp);
		bd[0] = r1[0]*bdp[0]+r2[0]*bdp[1]+r3[0]*bdp[2];
		bd[1] = r1[1]*bdp[0]+r2[1]*bdp[1]+r3[1]*bdp[2];
		bd[2] = r1[2]*bdp[0]+r2[2]*bdp[1]+r3[2]*bdp[2];
		box2pol(bd,pd);
		pd1[i] = pd[0];
		pd2[i] = pd[1];
	}	

	return ;
}
void pol2box(double p[2],double b[3])
{
	double thita=p[0],	phai=p[1];
	b[0]=cos(phai)*sin(thita);
	b[1]=sin(phai)*sin(thita);
	b[2]=cos(thita);
	return;
}
void box2pol(double bp[3],double pp[2])
{
	double x=bp[0],y=bp[1],z=bp[2];
	double thitap=acos(z);
	double phaip=acos(x/sin(thitap));
	phaip=phaip*(y>0)+(2*pi-phaip)*(y<=0);
	pp[0]=thitap;pp[1] = phaip;
	return ;
}
double sphere_trans(double longiA,double latiA,double longiB,double latiB)
{
	double R=6371;
	longiB=longiB*pi/180;
	longiA=longiA*pi/180;
	latiB=latiB*pi/180;
	latiA=latiA*pi/180;
	double a=latiB;
	double b=latiA;
	double belta=abs(longiA-longiB);
	double c=acos(sin(a)*sin(b)+cos(a)*cos(b)*cos(belta));
	double dis=abs(R*c);
	return dis;
}
void map(string Path,double longitudeA,double latitudeA,double longitudeB,double latitudeB,int n,double d[],double h[],int ch[])
{
	double AllAboutMap[33][12] = {0};
	HDR_REF(Path,AllAboutMap);
	double NROWS[33],NCOLS[33],NODATA[33],ULXMAP[33],ULYMAP[33],XDIM[33],YDIM[33];
    for(int i=0;i<33;i++)
	{
		NROWS[i]=AllAboutMap[i][1-1];
		NCOLS[i]=AllAboutMap[i][2-1];
		NODATA[i]=AllAboutMap[i][8-1];
		ULXMAP[i]=AllAboutMap[i][9-1];
		ULYMAP[i]=AllAboutMap[i][10-1];
		XDIM[i]=AllAboutMap[i][11-1];
		YDIM[i]=AllAboutMap[i][12-1];
	}
	
	string DEM[33]={   "W180N90.DEM ",	"W140N90.DEM ",	"W100N90.DEM ",	"W060N90.DEM ",
			"W020N90.DEM ",	"E020N90.DEM ",	"E060N90.DEM ",	"E100N90.DEM ",
			"E140N90.DEM ",	"W180N40.DEM ",	"W140N40.DEM ",	"W100N40.DEM ",
			"W060N40.DEM ",	"W020N40.DEM ",	"E020N40.DEM ",	"E060N40.DEM ",
			"E100N40.DEM ",	"E140N40.DEM ",	"W180S10.DEM ",	"W140S10.DEM ",
			"W100S10.DEM ",	"W060S10.DEM ",	"W020S10.DEM ",	"E020S10.DEM ",
			"E060S10.DEM ",	"E100S10.DEM ",	"E140S10.DEM ",	"W180S60.DEM ",
			"W120S60.DEM ",	"W060S60.DEM ",	"W000S60.DEM ",	"E060S60.DEM ",
			"E120S60.DEM "};
	//Path='..\Map\';

	double p1[2]={PI/2-latitudeA*PI/180,longitudeA*PI/180};
	double p2[2]={PI/2-latitudeB*PI/180,longitudeB*PI/180};
	double *pd1=new double[n];
	double *pd2 = new double[n];
	DividePath(p1,p2,n,pd1,pd2);
	for (int i=0;i<n+2;i++)
	{
		h[i]=0;
	}
	double R=6371;
	double longitude,latitude;
	double sea_symbol = -100000;
	int mapc = 0;
	string DEMfile[33];
	FILE *f2 = NULL;
	const char *filename[33];
	unsigned char a,b;
	double c;
	double thita,phai;
	double Col,Row;
	long mStartPos;
	int status;
	for (int ct=0;ct<=n+1;ct++)
	{
		c=0;
		if (ct==0)
		{    longitude=longitudeA;latitude=latitudeA;}
		else if (ct>0 && ct<=n)
		{
			thita=pd1[ct-1];
			phai=pd2[ct-1];
			longitude=(phai<=PI)*phai*180/PI+(phai>PI)*(-360+phai*180/PI);
			latitude=90-180*thita/PI;
		}
		else if (ct==n+1)
		{    longitude=longitudeB;latitude=latitudeB;}

		for (mapc=0;mapc<33;mapc++)
		{
			if ((longitude-ULXMAP[mapc])*(longitude-(NCOLS[mapc]*XDIM[mapc]+ULXMAP[mapc]))<=0 
					&& (latitude-ULYMAP[mapc])*(latitude-(-NROWS[mapc]*YDIM[mapc]+ULYMAP[mapc]))<=0)
			{
				DEMfile[mapc]=Path;
				DEMfile[mapc].append(DEM[mapc]);
				filename[mapc] = DEMfile[mapc].data();
				fopen_s(&f2,filename[mapc],"r");

				Col=abs(round((longitude-ULXMAP[mapc])/XDIM[mapc]+1));
				Row=abs(round((-latitude+ULYMAP[mapc])/YDIM[mapc]+1));
				if (Col>NCOLS[mapc]) 
                    Col = NCOLS[mapc];
                if(Row > NROWS[mapc])
                    Row = NROWS[mapc];
                
				mStartPos = long(2 * (Row - 1) * NCOLS[mapc] + 2 * Col - 2 );
				status=fseek(f2,mStartPos,0);
				fscanf_s(f2,"%c",&a);
				fscanf_s(f2,"%c",&b);
//				fread(a,8,1,f2);
//				fread(b,8,1,f2);
				c=a*256+b;
				if (c > (pow(2.0, 15) - 1))
					c= c-pow(2.0,16);
				if (c==NODATA[mapc])
	//                    c=0;
					sea_symbol = c;
				h[ct]=c;
				fclose(f2);
				break;
			}
		}
	}

	double dis=sphere_trans(longitudeA,latitudeA,longitudeB,latitudeB);
	d[0] = 0;
	for(int i=1;i<n+2;i++)
		d[i]=d[i-1]+dis/(n+1);
	double dB = 0;
	for (int i=1;i<=n+2;i++)        
	{
		if (h[i-1]!=0 && h[i-1]>-100)
			if (h[i-1]<100 && d[i-1] - dB <=50)
				ch[i-1]=1;
			else
				ch[i-1]=2;
          
		else
		{
			ch[i-1]=3;
			h[i-1] = 0;
			dB = d[i-1];
		}
	}
	dB = -1;
	for (int i=n+2;i>=1;i=i-1)
	{
		if (h[i-1]!=0 && h[i-1]>-100)
		{
			if (h[i-1]<100 && dB - d[i-1] <=50 && dB!=-1 && ch[i-1]!=1)
				ch[i-1]=1;
		}
		else
			dB = d[i-1];
	}
}
int find_intervals(int series[] ,int len_ser,int k1[], int k2[])
{
//find_intervals Find all intervals with consecutive 1's
//     [k1, k2] = find_intervals(series)
//     This function finds all 1's intervals, namely, the indices when the
//     intervals start and where they end
//
//     For example, for the input indices
//           0 0 1 1 1 1 0 0 0 1 1 0 0
//     this function will give back
//       k1 = 3, 10
//       k2 = 6, 11
//
//     Input arguments:
//     indices -   vector containing zeros and ones
//
//     Output arguments:
//     k1      -   vector of start-indices of the found intervals
//     k2      -   vector of end-indices of the found intervals
//
//     Example:
//     [k1, k2] = find_intervals(indices)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
	int count = 0;
	if (series[0] ==1 && series[1] == 1)
	{
		k1[count] = 0;
		count = count + 1;
	}

	for (int i = 1;i<len_ser-1;i++)
	{
		if ((series[i-1] == 0 ) && (series[i] == 1))
		{
			k1[count] = i;
			count = count + 1;
		}
	}
	count = 0;
	for (int i = 1;i<len_ser-1;i++)
	{
		if ((series[i] == 1) && (series[i+1] == 0 ))
		{
			k2[count] = i;
			count = count + 1;
		}
	}

	if (series[len_ser-1] ==1 && series[len_ser-2] == 1)
	{
		k2[count] = len_ser-1;
		count = count + 1;
	}
	return count;
}
double longest_cont_dist(double d[], int zone[], int zone_r,int len_d)
{
	//longest_cont_dist Longest continuous path belonging to the zone_r
//     dm = longest_cont_dist(d, zone, zone_r)
//     This function computes the longest continuous section of the
//     great-circle path (km) for a given zone_r
//
//     Input arguments:
//     d       -   vector of distances in the path profile
//     zone    -   vector of zones in the path profile
//     zone_r  -   reference zone for which the longest continuous section
//                 is computed
//
//     Output arguments:
//     dm      -   the longest continuous section of the great-circle path (km) for a given zone_r
//
//     Example:
//     dm = longest_cont_dist(d, zone, zone_r)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
//     v1    12FEB16     Ivica Stevanovic, OFCOM         included zone_r==12
	double dm = 0;
	int *start = new int[len_d];
	int *stop = new int[len_d];
	int *tmp = new int[len_d];
	int len_start = 0;
	double delta;
	if (zone_r  == 12)
	{
		for(int i=0;i<len_d;i++)
		{
			if(zone[i]==1 || zone[i]==2)
				tmp[i] = 1;
			else
				tmp[i]=0;
		}
		len_start = find_intervals(tmp,len_d,start,stop);
	}
	else
	{
		for(int i=0;i<len_d;i++)
		{
			if(zone[i]==zone_r)
				tmp[i] = 1;
			else
				tmp[i]=0;
		}
		len_start = find_intervals(tmp,len_d,start,stop);
	}

	int n = len_start;

	for (int i =0;i<n;i++)
	{
		delta = 0;
		if (d[stop[i]]<d[len_d-1])
			delta = delta + ( d[stop[i]+1]-d[stop[i]] )/2.0;

    
		if (d[start[i]]>0)
			delta = delta + ( d[stop[i]]-d[stop[i]-1] )/2.0;
		if(dm<d[stop[i]]-d[start[i]] + delta)
			dm = d[stop[i]]-d[start[i]] + delta;
	}
	return dm;
}
double beta0(double phi, double dtm, double dlm)
{
	////
//     This function computes the time percentage for which refractive index
//     lapse-rates exceeding 100 N-units/km can be expected in the first 100
//     m of the lower atmosphere
//     as defined in ITU-R P.452-16.
//
//     Input arguments:
//     phi     -   path centre latitude (deg)
//     dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
//     dlm     -   the longest continuous inland section of the great-circle path (km)
//
//     Output arguments:
//     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
//
//     Example:
//     b0 = beta0(phi, dtm, dlm)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab

	double b0,tau,mu1,mu4;
	tau = 1- exp(-(4.12*(1e-4)*pow(dlm,2.41)));       // (3a)

	mu1 = pow(pow(10.0,(-dtm/(16-6.6*tau))) + pow(10,(-5*(0.496 + 0.354*tau))),0.2); // (3)

	if (mu1 > 1)
		mu1 = 1;
   
	if (abs(phi) <= 70)
	{
		mu4 = pow(10.0, (-0.935 + 0.0176*abs(phi))*log10(mu1) );   // (4)
	    b0 = pow(10.0, -0.015*abs(phi) + 1.67 )*mu1*mu4;           // (2)  
	}
	else
	{
		mu4 = pow(10.0,0.3*log10(mu1));                            // (4)
	    b0 = 4.17*mu1*mu4;                                    // (2)
	}

	return b0;
}
double earth_rad_eff(double DN)
//earth_rad_eff Median value of the effective Earth radius
//     [ae, ab] = earth_rad_eff(DN)
//     This function computes the median value of the effective earth
//     radius, and the effective Earth radius exceeded for beta0// of time
//     as defined in ITU-R P.452-16.
//
//     Input arguments:
//     DN      -   the average radiorefractive index lapse-rate through the
//                 lowest 1 km of the atmosphere (N-units/km)
//
//     Output arguments:
//     ae      -   the median effective Earth radius (km)
//     ab      -   the effective Earth radius exceeded for beta0 // of time
//
//     Example:
//     [ae, ab] = earth_rad_eff(DN)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
{
	return  6371*157/(157-DN);     // (5)

}
double equiv_annual_percent(double pw, double theta, double omega)
{
// Function to compute equivalent annual time //p of interference exceedance (or // time that basic TX loss
// is NOT exceeded) for a given worst month time percentage (pw//)

// Reference: ITU-R P.452-15 Section 3.2.1

// Inputs:
// pw = the worst month time // (//)
// theta = the interference path center latitude (deg)
// omega = the fraction of the interference path over water (fraction) (See cell F13 on Step_2.)

// Outputs:
// p = the equivalent annual time percentage corresponding to the worst
// month time percentage (//)
	double deg = pi / 180;
	double p;
	double GL;
	// Equation (1a)
	if (abs(theta) <= 45) 
		GL = sqrt(1.1 + pow(abs(cos(2 * theta * deg)),0.7) );
	else
		GL = sqrt(1.1 - pow(abs(cos(2 * theta * deg)),0.7) );


	// Equation (1)
	p = pow(10.0,((log(pw) / log(10.0)) + (log(GL) / log(10.0)) - 0.186 * omega - 0.444) / (0.816 + 0.078 * omega));
	// See condition under Equation (1a)
	if (12 * p < pw) 
		p = pw / 12;
	return p;
}
double  path_fraction(double d[], int zone[], int zone_r,int len_d)
{
//path_fraction Path fraction belonging to a given zone_r
//     omega = path_fraction(d, zone, zone_r)
//     This function computes the path fraction belonging to a given zone_r
//     of the great-circle path (km) 
//
//     Input arguments:
//     d       -   vector of distances in the path profile
//     zone    -   vector of zones in the path profile
//     zone_r  -   reference zone for which the fraction is computed
//
//     Output arguments:
//     omega   -   path fraction belonging to the given zone_r
//
//     Example:
//     omega = path_fraction(d, zone, zone_r)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    02FEB16     Ivica Stevanovic, OFCOM         First implementation in matlab

	double dm = 0;
	int *start = new int[len_d];
	int *stop = new int[len_d];
	int *tmp = new int[len_d];
	for(int i=0;i<len_d;i++)
	{
		if(zone[i]==zone_r)
			tmp[i]=1;
		else
			tmp[i]=0;
	}
	double delta=0;
	int len_start = find_intervals(tmp,len_d,start,stop);
	for (int i = 0;i<len_start;i++)
	{	
		delta = 0;
		if (d[stop[i]]<d[len_d-1])
			delta = delta + ( d[stop[i]+1]-d[stop[i]] )/2.0;
		
    
		if (d[start[i]]>0)
			delta = delta + ( d[stop[i]]-d[stop[i]-1] )/2.0;
    
	   dm = dm + d[stop[i]]-d[start[i]] + delta;
   
	}

	double omega = dm/(d[len_d-1]-d[0]);

	return omega;
}
int closs_corr(double f, double d[], double h[], int zone[], double htg, double hrg, double ha_t, double ha_r, double dk_t, double dk_r,int len_d,
	double dc[], double hc[], int zonec[],double htrgc[2], double Ahtr[2])
{
//closs clutter loss correction according to P.452-16
//   function [dc,hc,zonec,htgc,htrc, Aht, Ahr] = closs_corr(d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
//
//   This function computes the height-gain correction as defined in ITU-R P.452-16 (Section 4.5.4)
//
//     Input parameters:
//     f       -   Frequency (GHz)
//     d       -   vector of distances di of the i-th profile point (km)
//     h       -   vector of heights hi of the i-th profile point (meters
//                 above mean sea level. Both vectors contain n+1 profile points
//     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
//     htg     -   Tx Antenna center heigth above ground level (m)
//     hrg     -   Rx Antenna center heigth above ground level (m)
//     ha_t    -   Nominal clutter height at the transmitting end (m, agl)
//     ha_r    -   Nominal clutter height at the receiving end (m, agl)
//     dk_t    -   distance from nominal clutter point to the Tx antenna (km)
//     dk_r    -   distance from nominal clutter point to the Rx antenna (km)
//
//     Output parameters:
//     dc      -   vector of distances in the height-gain model
//     hc      -   vector of heights in the height-gain model
//     zonec   -   Zone type: Coastal land (1), Inland (2) or Sea (3)in the height-gain model
//     htgc    -   Tx Antenna center heigth above ground level (m) in the height-gain model
//     hrgc    -   Rx Antenna center heigth above ground level (m) in the height-gain model
//     Aht     -   Additional losses to account for clutter shielding the
//     Ahr         transmitter and receiver. These should be set to zero if there is no
//                 such shielding
//
//     Example:
//     [dc,hc,zonec,htgc,hrgc, Aht, Ahr] = closs_corr(d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
//
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    09MAR16     Ivica Stevanovic, OFCOM         Initial version
//     v1    24NOV16     Ivica Stevanovic, OFCOM         introduced clutter nominal distance check
	int index1 = 0;
	int index2 = len_d-1;
	htrgc[0] = htg;
	htrgc[1] = hrg;
	Ahtr[0] = 0;
	Ahtr[1] = 0;

	double ha = ha_t;
	double dk = dk_t;
	double Ffc;
	int flagAht,flagAhr,kk=-1;
	if (ha > htg)
	{
		Ffc = 0.25+0.375*(1+tanh(7.5*(f-0.5)));  // (57a)
		Ahtr[0] = 10.25*Ffc*exp(-dk)*( 1- tanh(6*(htg/ha-0.625)) )-0.33; // (57)
		flagAht = 1;
		for(int i=0;i<len_d;i++)
		{
			if(d[i]>=dk)
			{
				kk = i;
				break;
			}
		}
		if (kk!=-1)
			index1 = kk;
		else
			index1 = len_d-1;
		htrgc[0] = ha_t;
	}
	ha = ha_r;
	dk = dk_r;
	kk = -1;
	if (ha > hrg)
	{
		Ffc = 0.25+0.375*(1+tanh(7.5*(f-0.5)));  // (57a)
		Ahtr[1] = 10.25*Ffc*exp(-dk)*( 1- tanh(6*(hrg/ha-0.625)) )-0.33;  // (57)
		flagAhr = 1;
		for(int i=0;i<len_d;i++)
		{
			if(d[i]<=d[len_d-1]-dk)
				kk = i;
		}

		//kk = find(d <= d(end)-dk);
		if(kk!=-1)
			index2 = kk;
		else
			index2 = 0;
    
		htrgc[1] = ha_r;
	}

	// Modify the path

	if (index2-index1 < 3) // at least two points between the clutter at Tx and Rx sides
		//error('tl_p452: closs_corr: the sum of clutter nominal distances is larger than the path length.');
		return -1;
    for(int i = 0;i<index2-index1+1;i++)
	{
		dc[i] = d[index1+i]-d[index1];
		hc[i] = h[index1+i];
		zonec[i] = zone[index1+i];
	}
	return index2-index1+1;
}
void smooth_earth_heights(double d[], double h[], double htg,double  hrg,double ae,double f,int len_d,//input
	double hstr[2], double hstrd[2], double htre[2], double hm[1], double dltr[2],double theta_tr[2], double theta_tot[1], int pathtype[1])//output
{
//smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.452-16
// [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f)
// This function derives smooth-Earth effective antenna heights according to
// Sections 4 and 5 of the Annex 2 of ITU-R P.452-16
//
// Input parameters:
// d         -   vector of terrain profile distances from Tx [0,dtot] (km)
// h         -   vector of terrain profile heigths amsl (m)
// htg, hrg  -   Tx and Rx antenna heights above ground level (m)
// ae        -   median effective Earth's radius (c.f. Eq (6a))
// f         -   frequency (GHz)
//
// Output parameters:
//
// hst, hsr     -   Tx and Rx antenna heigts of the smooth-Earth surface amsl (m)
// hstd, hsrd   -   Tx and Rx effective antenna heigts for the diffraction model (m)
// hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
// hm           -   The terrain roughness parameter (m)
// dlt          -   interfering antenna horizon distance (km)
// dlr          -   Interfered-with antenna horizon distance (km)
// theta_t      -   Interfering antenna horizon elevation angle (mrad)
// theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
// theta_tot    -   Angular distance (mrad)
// pathtype     -   1 = 'los', 2 = 'transhorizon'
//
// Rev   Date        Author                          Description
// -------------------------------------------------------------------------------
// v0    15JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
// v1    15JUN16     Ivica Stevanovic, OFCOM         Modifications related to LoS path

	int n = len_d;

	double dtot = d[n-1];

	//Tx and Rx antenna heights above mean sea level amsl (m)
	double hts = h[0] + htg;
	double hrs = h[n-1] + hrg;

	// Section 5.1.6.2

	double v1 = 0;
	for (int ii = 1;ii<n;ii++)
		v1 = v1 + (d[ii]-d[ii-1])*(h[ii]+h[ii-1]);  // Eq (161)
	
	double v2 = 0;
	for (int ii = 1;ii<n;ii++)
		v2 = v2 + (d[ii]-d[ii-1])*( h[ii]*( 2*d[ii] + d[ii-1] ) + h[ii-1] * ( d[ii] + 2*d[ii-1] ) );  // Eq (162)


	hstr[0] = (2*v1*dtot - v2)/dtot/dtot;       // Eq (163)
	hstr[1] = (v2- v1*dtot)/dtot/dtot;          // Eq (164)

	// Section 5.1.6.3

	double *HH = new double[len_d];
	for(int ii = 0;ii<len_d;ii++)
		HH[ii]= h[ii] - (hts*(dtot-d[ii]) + hrs*d[ii])/dtot;  // Eq (165d)

	double hobs = HH[1],alpha_obt=HH[1]/d[1],alpha_obr=HH[1]/(dtot-d[1]);
	for(int ii = 1;ii<n-1;ii++)
	{
		if(hobs<=HH[ii])
			hobs = HH[ii];
		if(alpha_obt<=HH[ii]/d[ii])
			alpha_obt = HH[ii]/d[ii];
		if(alpha_obr<=HH[ii]/(dtot-d[ii]))
			alpha_obr = HH[ii]/(dtot-d[ii]);
	}
	//hobs = max(HH(2:n-1));                 // Eq (165a)

	//alpha_obt = max( HH(2:n-1)./d(2:n-1) ); // Eq (165b)

	//alpha_obr = max( HH(2:n-1)./( dtot - d(2:n-1) ) ); // Eq (165c)

	// Calculate provisional values for the Tx and Rx smooth surface heights

	double gt = alpha_obt/(alpha_obt + alpha_obr);         // Eq (166e)
	double gr = alpha_obr/(alpha_obt + alpha_obr);         // Eq (166f)
	double hstp,hsrp;
	if (hobs <= 0)
	{
		hstp = hstr[0];                                 // Eq (166a)
		hsrp = hstr[1];                                 // Eq (166b)
	}
	else
	{
		hstp = hstr[0] - hobs*gt;                       // Eq (166c)
		hsrp = hstr[1] - hobs*gr;                       // Eq (166d)
	}

	// calculate the final values as required by the diffraction model

	if (hstp >= h[0])
		hstrd[0] = h[0];                                // Eq (167a)
	else
		hstrd[0] = hstp;                                // Eq (167b)

	if (hsrp > h[len_d-1])
		hstrd[1] = h[len_d-1];                              // Eq (167c)
	else
		hstrd[1] = hsrp;                                // Eq (167d)


	// Interfering antenna horizon elevation angle and distance
	double *theta = new double[n-2];
	for(int ii = 1;ii<n-1;ii++)
	{
		theta[ii-1] = 1000 * atan( (h[ii] - hts)/(1000 * d[ii] ) - d[ii]/(2*ae) );
	}

	//theta(theta < 0) = 0;  // condition below equation (152)

	double theta_t = max(theta,n-2);                           // Eq (154)
	
	double theta_td = 1000 * atan( (hrs - hts)/(1000 * dtot ) - dtot/(2*ae) );  // Eq (153)
	double theta_rd = 1000 * atan( (hts - hrs)/(1000 * dtot ) - dtot/(2*ae) );  // not defined in the recommendation

	if (theta_t > theta_td )  // Eq (150): test for the trans-horizon path
		pathtype[0] = 2; //transhorizon
	else
		pathtype[0] = 1; //los

	int kindex = -1;
	for (int ii=0;ii<n-2;ii++)
	{
		if(theta[ii]==theta_t)
			{kindex = ii;break;}
	}
	int lt = kindex+1;


	dltr[0] = d[lt];                             // Eq (155)

	// Interfered-with antenna horizon elevation angle and distance
	for(int ii = 1;ii<n-1;ii++)
	{
		theta[ii-1] = 1000 * atan( (h[ii] - hrs)/(1000 * (dtot - d[ii]) ) - (dtot - d[ii])/(2*ae) );
	}
	

	//theta(theta < 0) = 0;

	double theta_r = max(theta,n-2);

	for (int ii=0;ii<n-2;ii++)
	{
		if(theta[ii]==theta_r)
			{kindex = ii;}
	}
	int lr = kindex+1;

	dltr[1] = dtot - d[lr];                            // Eq (158)
	double lambda,Ce,numax;
	double *nu =new double[n-2];
	if (pathtype[0] == 1)
	{
		theta_t = theta_td;
		theta_r = theta_rd;
		lambda = 0.3/f;
		Ce = 1.0/ae;
		for(int ii=1;ii<n-1;ii++)
			nu[ii-1] = (h[ii] + 500*Ce*d[ii]*(dtot-d[ii])- (hts*(dtot- d[ii]) + hrs *d[ii])/dtot)* 
			 sqrt(0.002*dtot/(lambda*d[ii]*(dtot-d[ii])));
		numax = max(nu,n-2);
		for(int ii=0;ii<n-2;ii++)
			if(nu[ii]==numax)
				kindex = ii;
		//kindex = find(nu == numax);
		lt = kindex+1;  
		dltr[0] = d[lt];  
		dltr[1] = dtot - dltr[0];
		for(int ii=1;ii<n-1;ii++)
			if(dltr[1] <=dtot -d[ii])
				kindex = ii;
		lr = kindex+1;
	}

	// Angular distance

	theta_tot[0] = 1e3 * dtot/ae + theta_t + theta_r;         // Eq (159)

	// Section 5.1.6.4 Ducting/layer-reflection model

	// Calculate the smooth-Earth heights at transmitter and receiver as
	// required for the roughness factor
	if( hstr[0]>h[0])
		hstr[0] = h[0];                           // Eq (168a)
	if(hstr[1]>h[n-1])
		hstr[1] = h[n-1];                         // Eq (168b)

	// Slope of the smooth-Earth surface

	double m = (hstr[1] - hstr[0])/ dtot;                          // Eq (169)

	// The terminal effective heigts for the ducting/layer-reflection model

	htre[0] = htg + h[0] -   hstr[0];                       // Eq (170)
	htre[1] = hrg + h[n-1] - hstr[1];                       
	double *tmp_hm = new double[lr-lt+1];
	int kk = 0;
	for(int ii = lt-1;ii<lr;ii++)
	{
		tmp_hm[kk] = h[ii] - (hstr[0] + m*d[ii]);
		kk++;
	}
	hm[0] = max(tmp_hm,lr-lt+1);
	theta_tr[0] = theta_t;theta_tr[1] = theta_r;
	return;
}
double pl_los(double d, double f, double p, double b0, double w, double temp, double press, double dlt, double dlr,double Lb0pb[2])
{
//pl_los Line-of-sight transmission loss according to ITU-R P.452-16
//     This function computes line-of-sight transmission loss (including short-term effects)
//     as defined in ITU-R P.452-16.
//
//     Input parameters:
//     d       -   Great-circle path distance (km)
//     f       -   Frequency (GHz)
//     p       -   Required time percentage(s) for which the calculated basic
//                 transmission loss is not exceeded (//)
//     b0      -   Point incidence of anomalous propagation for the path
//                 central location (//)
//     w       -   Fraction of the total path over water (//)
//     temp    -   Temperature (degrees C)
//     press   -   Dry air pressure (hPa)
//     dlt     -   For a transhorizon path, distance from the transmit antenna to
//                 its horizon (km). For a LoS path, each is set to the distance
//                 from the terminal to the profile point identified as the Bullington
//                 point in the diffraction method for 50// time
//     dlr     -   For a transhorizon path, distance from the receive antenna to
//                 its horizon (km). The same note as for dlt applies here.
//
//     Output parameters:
//     Lbfsg   -   Basic transmission loss due to free-space propagation and
//                 attenuation by atmospheric gases
//     Lb0p    -   Basic transmission loss not exceeded for time percentage, p//, due to LoS propagation
//     Lb0b    -   Basic transmission loss not exceedd for time percentage, b0//, due to LoS propagation
//
//     Example:
//     [Lbfsg, Lb0p, Lb0b] = pl_los(d, f, p, b0, w, temp, press, dlt, dlr)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    04FEB14     Ivica Stevanovic, OFCOM         First implementation in matlab

	double T = temp + 273.15;

	// water vapor density
	double rho = 7.5 + 2.5 * w;  // (9a)

	// compute specific attenuation due to dry air and water vapor:

	double g_0 = cal_gammao_p676_11(f,press,T,rho);
	double g_w = cal_gammaw_p676_11(f,press,T,rho);
	double Ag = (g_0 + g_w) * d;  //(9)

	// Basic transmission loss due to free-space propagation and attenuation
	// by atmospheric gases
	double Lbfsg = 92.5 + 20.0*log10(f) + 20.0*log10(d) + Ag;  // (8)

	// Corrections for multipath and focusing effects at p and b0
	double Esp = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(p/50);   //(10a)
	double Esb = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(b0/50);  //(10b)

	// Basic transmission loss not exceeded for time percentage p// due to
	// LoS propagation
	double Lb0p = Lbfsg + Esp;    //(11)

	// Basic transmission loss not exceeded for time percentage b0// due to
	// LoS propagation
	double Lb0b = Lbfsg + Esb;    //(12)
	Lb0pb[0] = Lb0p;Lb0pb[1] = Lb0b;
	return Lbfsg;
}
double dl_bull(double d[], double h[], double hts,double  hrs,double  ap,double  f,int len_d)
{
//dl_bull Bullington part of the diffraction loss according to P.452-16
//   This function computes the Bullington part of the diffraction loss
//   as defined in ITU-R P.452-16 in 4.2.1
//
//     Input parameters:
//     d       -   vector of distances di of the i-th profile point (km)
//     h       -   vector of heights hi of the i-th profile point (meters
//     above mean sea level. Both vectors contain n+1 profile points
//     hts     -   transmitter antenna height in meters above sea level (i=0)
//     hrs     -   receiver antenna height in meters above sea level (i=n)
//     ap      -   the effective earth radius in kilometers
//     f       -   frequency expressed in GHz
//
//     Output parameters:
//     Lbull   -   Bullington diffraction loss for a given path
//
//     Example:
//     Lbull = dl_bull(d, h, hts, hrs, ap, f)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation in matlab


//// Body of function

// Effective Earth curvature Ce (km^-1)

	double Ce = 1/ap;

	// Wavelength in meters

	double lambda = 0.3/f;

	// Complete path length

	double dtot = d[len_d-1]-d[0];

	// Find the intermediate profile point with the highest slope of the line
	// from the transmitter to the point
	double *di =new double[len_d-2];
	double *hi = new double[len_d-2];
	double Stim = (h[1] + 500*Ce*d[1]*(dtot - d[1]) - hts)/d[1];
	for(int i=0;i<len_d-2;i++)
	{
		di[i] = d[i+1];
		hi[i] = h[i+1];
		if(Stim<(hi[i] + 500*Ce*di[i]*(dtot - di[i]) - hts)/di[i])
			Stim = (hi[i] + 500*Ce*di[i]*(dtot - di[i]) - hts)/di[i];
	}

	// Calculate the slope of the line from transmitter to receiver assuming a
	// LoS path

	double Str = (hrs - hts)/dtot;                                         // Eq (15)
	double numax,Luc,nub,Srim,dbp;
	numax = ( hi[0] + 500*Ce*di[0]*(dtot - di[0]) - ( hts*(dtot - di[0]) + hrs*di[0])/dtot ) *
					   sqrt(0.002*dtot/(lambda*di[0]*(dtot-di[0])));
	Srim = (hi[0] + 500*Ce*di[0]*(dtot-di[0])-hrs)/(dtot-di[0]);
	if (Stim < Str) // Case 1, Path is LoS
	{    
		// Find the intermediate profile point with the highest diffraction
		// parameter nu:
		for(int i = 0;i<len_d-2;i++)
			if(numax<( hi[i] + 500*Ce*di[i]*(dtot - di[i]) - ( hts*(dtot - di[i]) + hrs*di[i])/dtot ) *
					   sqrt(0.002*dtot/(lambda*di[i]*(dtot-di[i]))))
				numax =( hi[i] + 500*Ce*di[i]*(dtot - di[i]) - ( hts*(dtot - di[i]) + hrs*di[i])/dtot ) *
					 sqrt(0.002*dtot/(lambda*di[i]*(dtot-di[i])));
		Luc = 0;
		if (numax > -0.78)
			Luc = 6.9 + 20*log10(sqrt(pow(numax-0.1,2.0)+1) + numax - 0.1);   // Eq (13), (17)
	}
	else
	{    
		// Path is transhorizon
    
		// Find the intermediate profile point with the highest slope of the
		// line from the receiver to the point
		for(int i=0;i<len_d-2;i++)
			if(Srim<(hi[i] + 500*Ce*di[i]*(dtot-di[i])-hrs)/(dtot-di[i]))
				Srim = (hi[i] + 500*Ce*di[i]*(dtot-di[i])-hrs)/(dtot-di[i]);     // Eq (18)
    
		// Calculate the distance of the Bullington point from the transmitter:
    
		dbp = (hrs - hts + Srim*dtot)/(Stim + Srim);                // Eq (19)
    
		// Calculate the diffraction parameter, nub, for the Bullington point
    
		nub =  ( hts + Stim*dbp - ( hts*(dtot - dbp) + hrs*dbp)/dtot ) *
					   sqrt(0.002*dtot/(lambda*dbp*(dtot-dbp)));    // Eq (20)
    
		// The knife-edge loss for the Bullington point is given by
              
		Luc = 0;
		if (nub > -0.78)
			Luc = 6.9 + 20*log10(sqrt(pow(nub-0.1,2.0)+1) + nub - 0.1);   // Eq (13), (21)
	}

	// For Luc calculated using either (17) or (21), Bullington diffraction loss
	// for the path is given by

	double Lbull = Luc + (1 - exp(-Luc/6.0))*(10+0.02*dtot);         // Eq (22)
	return Lbull;
}
void dl_se_ft_inner(double epsr, double sigma,double d,double hte,double hre,double adft,double f,double Ldft[2])
{
////dl_se_ft_inner The inner routine of the first-term spherical diffraction loss
//   This function computes the first-term part of Spherical-Earth diffraction
//   loss exceeded for p// time for antenna heights
//   as defined in Sec. 4.2.2.1 of the ITU-R P.452-16, equations (30-37)
//
//     Input parameters:
//     epsr    -   Relative permittivity
//     sigma   -   Conductivity (S/m)
//     d       -   Great-circle path distance (km)
//     hte     -   Effective height of interfering antenna (m)
//     hre     -   Effective height of interfered-with antenna (m)
//     adft    -   effective Earth radius (km)
//     f       -   frequency (GHz)
//
//     Output parameters:
//     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p// time
//                implementing equations (30-37), Ldft(1) is for horizontal
//                and Ldft(2) for the vertical polarization
//
//     Example:
//     Ldft = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation in matlab


//// Body of the function

// Normalized factor for surface admittance for horizontal (1) and vertical
// (2) polarizations
	double K[2]={0};
	double beta_dft[2]={0};
	double X[2]={0},Yt[2]={0},Yr[2]={0},Fx[2]={0};
	
	K[0]= 0.036* pow(adft*f,(-1.0/3)) * pow( pow(epsr-1,2.0) + pow(18*sigma/f,2.0),(-1/4));   // Eq (30a)

	K[1] = K[0] * pow(pow(epsr,2.0) + pow(18*sigma/f,2.0),0.5);       // Eq (30b)

// Earth ground/polarization parameter
	for(int i=0;i<2;i++)
	{
		beta_dft[i] = (1 + 1.6*K[i]*K[i] + 0.67* pow(K[i],4.0))/( 1 + 4.5* K[i]*K[i] + 1.53* pow(K[i],4.0));  // Eq (31)
		X[i] = 21.88* beta_dft[i] * pow(f/ pow(adft,2.0),(1.0/3)) * d;          // Eq (32)
	// Normalized transmitter and receiver heights
		Yt[i] = 0.9575* beta_dft[i] * pow(pow(f,2.0) / adft,(1.0/3)) * hte;       // Eq (33a)
		Yr[i] = 0.9575* beta_dft[i] * pow(pow(f,2.0) / adft,(1.0/3)) * hre;       // Eq (33b)
	}
// Normalized distance
// Calculate the distance term given by:
	double Bt[2]={0},Br[2]={0},GYt[2],GYr[2];
for (int ii = 0;ii<2;ii++)
{
    if (X[ii] >= 1.6)
        Fx[ii] = 11 + 10*log10(X[ii]) - 17.6*X[ii];
    else
        Fx[ii] = -20*log10(X[ii]) - 5.6488* pow(X[ii],1.425);     // Eq (34)
	Bt[ii] = beta_dft[ii]*Yt[ii];
	Br[ii] = beta_dft[ii]*Yr[ii];
	if (Bt[ii]>2)
        GYt[ii] = 17.6*pow(Bt[ii] - 1.1,0.5) - 5*log10(Bt[ii] -1.1)-8;
    else
        GYt[ii] = 20*log10(Bt[ii] + 0.1* pow(Bt[ii],3.0));
    
    
    if (Br[ii]>2)
        GYr[ii] = 17.6*pow(Br[ii] - 1.1,0.5) - 5*log10(Br[ii] -1.1)-8;
    else
        GYr[ii] = 20*log10(Br[ii] + 0.1* pow(Br[ii],3.0));
    
    
    if (GYr[ii] < 2 + 20*log10(K[ii]))
        GYr[ii] = 2 + 20*log10(K[ii]);
    
    
    if (GYt[ii] < 2 + 20*log10(K[ii]))
        GYt[ii] = 2 + 20*log10(K[ii]);
	Ldft[ii] = -Fx[ii]-GYt[ii]-GYr[ii];
}


return ;
}
void dl_se_ft(double d,double hte,double hre,double adft,double f,double omega,double Ldft[2])
{
//dl_se_ft First-term part of spherical-Earth diffraction according to ITU-R P.452-16
//   This function computes the first-term part of Spherical-Earth diffraction
//   loss exceeded for p// time for antenna heights
//   as defined in Sec. 4.2.2.1 of the ITU-R P.452-16
//
//     Input parameters:
//     d       -   Great-circle path distance (km)
//     hte     -   Effective height of interfering antenna (m)
//     hre     -   Effective height of interfered-with antenna (m)
//     adft    -   effective Earth radius (km)
//     f       -   frequency (GHz)
//     omega   -   fraction of the path over sea
//
//     Output parameters:
//     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p// time
//                Ldft(1) is for the horizontal polarization
//                Ldft(2) is for the vertical polarization
//
//     Example:
//     Ldft = dl_se_ft(d, hte, hre, adft, f, omega)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation in matlab
//     v1    16DEC16     Ivica Stevanovic, OFCOM         corrected bug in dl_se_ft_inner function call where epsr and sigma order was interchanged


//// Body of function

// First-term part of the spherical-Earth diffraction loss over land

double epsr = 22;
double sigma = 0.003;

double Ldft_land[2]={0};
dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f,Ldft_land);

// First-term part of the spherical-Earth diffraction loss over sea

epsr = 80;
sigma = 5;
double Ldft_sea[2]={0};
dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f,Ldft_sea);


// First-term spherical diffraction loss 
for (int i=0;i<2;i++)
Ldft[i] = omega * Ldft_sea[i] + (1-omega)*Ldft_land[i];      // Eq (29)

return;
}

void dl_se(double d, double hte,double hre,double  ap,double  f,double  omega,double Ldsph[2])
{
//dl_se spherical-Earth diffraction loss exceeded for p// time according to ITU-R P.452-16
//   This function computes the Spherical-Earth diffraction loss exceeded
//   for p// time for antenna heights hte and hre (m)
//   as defined in Sec. 4.2.2 of the ITU-R P.452-16
//
//     Input parameters:
//     d       -   Great-circle path distance (km)
//     hte     -   Effective height of interfering antenna (m)
//     hre     -   Effective height of interfered-with antenna (m)
//     ap      -   the effective Earth radius in kilometers
//     f       -   frequency expressed in GHz
//     omega   -   the fraction of the path over sea
//
//     Output parameters:
//     Ldsph   -   The spherical-Earth diffraction loss not exceeded for p// time
//                 Ldsph(1) is for the horizontal polarization
//                 Ldsph(2) is for the vertical polarization
//
//     Example:
//     Ldsph = dl_se(d, hte, hre, ap, f, omega)
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    23DEC15     Ivica Stevanovic, OFCOM         Initial version
//     v1    01FEB16     Ivica Stevanovic, OFCOM         Introduced dl_se_ft



//// Body of function

// Wavelength in meters

	double lambda = 0.3/f;

	// Calculate the marginal LoS distance for a smooth path
	//double Ldsph[2]={0};
	double dlos = sqrt(2*ap) * (sqrt(0.001*hte) + sqrt(0.001*hre));    // Eq (23)
	double c,m,b,dse1,dse2,hse,hreq,aem;
	double Ldft[2]={0};
	if (d >= dlos)
	 {   // calculate diffraction loss Ldft using the method in Sec. 4.2.2.1 for 
		// adft = ap and set Ldsph to Ldft
    
		dl_se_ft(d, hte, hre, ap, f, omega,Ldsph);
		return ;
	}
	else
	{
		// calculate the smallest clearance between the curved-Earth path and
		// the ray between the antennas, hse
    
		c = (hte - hre)/(hte + hre);        // Eq (25d)
		m = 250*d*d/(ap*(hte +hre));        // eq (25e)
    
		b = 2*sqrt((m+1)/(3*m)) * cos( pi/3 + 1/3* acos( 3*c/2 * sqrt( 3*m/pow(m+1,3.0) ) ) );   // Eq (25c)
    
		dse1 = d/2*(1+b);           // Eq (25a)
		dse2 = d - dse1;            // Eq (25b) 
    
		hse = (hte - 500*dse1*dse1/ap)*dse2 + (hre - 500*dse2*dse2/ap)*dse1;
		hse = hse/d;                // Eq (24)
    
		// Calculate the required clearance for zero diffraction loss
    
		hreq = 17.456*sqrt(dse1 * dse2 * lambda/d);     // Eq (26)
    
		if (hse > hreq)
		{
			Ldsph[0] = 0;Ldsph[1] = 0;
			return;
		}
		else
		{
			// calculate the modified effective Earth radius aem, which gives
			// marginal LoS at distance d
        
			aem = 500*pow(d/( sqrt(hte) + sqrt(hre) ),2.0);     // Eq (27)
        
			// Use the method in Sec. 4.2.2.1 for adft ) aem to obtain Ldft
        
			 dl_se_ft(d, hte, hre, aem, f, omega,Ldft);
        
			if (Ldft < 0)
			{
				Ldsph[0] = 0;Ldsph[1] = 0;
				return;
			}
			else
			{
				Ldsph[0] = (1- hse/hreq)*Ldft[0];     // Eq (28)
				Ldsph[1] = (1- hse/hreq)*Ldft[1];     // Eq (28)
			}
			
		}
	}

	return;
}
void dl_delta_bull( double d[], double h[], double hts,double  hrs,double  hstd,double  hsrd,double  ap,double  f, double omega ,int len_d,double Ld[2])
{
//dl_delta_bull Complete 'delta-Bullington' diffraction loss model P.452-16
//   function Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega )
//
//   This function computes the complete 'delta-Bullington' diffraction loss
//   as defined in ITU-R P.452-16 (Section 4.2.3)
//
//     Input parameters:
//     d       -   vector of distances di of the i-th profile point (km)
//     h       -   vector of heights hi of the i-th profile point (meters
//     above mean sea level. Both vectors contain n+1 profile points
//     hts     -   transmitter antenna height in meters above sea level (i=0)
//     hrs     -   receiver antenna height in meters above sea level (i=n)
//     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
//     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
//     ap      -   the effective Earth radius in kilometers
//     f       -   frequency expressed in GHz
//     omega   -   the fraction of the path over sea
//
//     Output parameters:
//     Ld     -   diffraction loss for the general patha according to
//                Section 4.2.3 of ITU-R P.452-16. 
//                Ld(1) is for the horizontal polarization 
//                Ld(2) is for the vertical polarization
//
//     Example:
//     Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega )
//       
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


//// Body of function

// Use the method in 4.2.1 for the actual terrain profile and antenna
// heights. Set the resulting Bullington diffraction loss for the actual
// path to Lbulla

	double Lbulla = dl_bull(d, h, hts, hrs, ap, f,len_d);

	// Use the method in 4.2.1 for a second time, with all profile heights hi
	// set to zero and modified antenna heights given by

	double hts1 = hts - hstd;   // eq (38a)
	double hrs1 = hrs - hsrd;   // eq (38b)
	double *h1 = new double[len_d];
	for(int i=0;i<len_d;i++)
		h1[i]=0;

	// where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
	// resulting Bullington diffraction loss for this smooth path to Lbulls

	double Lbulls = dl_bull(d, h1, hts1, hrs1, ap, f,len_d);

	// Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
	// for the actual path length (dtot) with 

	double hte = hts1;             // eq (39a)
	double hre = hrs1;             // eq (39b)
	double dtot = d[len_d-1] - d[0];

	double Ldsph[2]={0};
	dl_se(dtot, hte, hre, ap, f, omega,Ldsph);

	// Diffraction loss for the general paht is now given by
	Ld[0]=Lbulla;Ld[1] = Lbulla;
	for(int i=0;i<2;i++)
		if(Ldsph[i] - Lbulls>0)
			Ld[i] = Ld[i]+Ldsph[i] - Lbulls;
	return;
}
double inv_cum_norm( double x )
{
	//inv_cum_norm approximation to the inverse cummulative normal distribution
//   I = inv_cum_norm( x )
//   This function implements an approximation to the inverse cummulative
//   normal distribution function for x <= 0.5 as defined in Attachment 3 to
//   Annex 1 of the ITU-R P.452-16
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version

	if (x < 0.000001)
		x = 0.000001;

	double tx = sqrt(-2*log(x));    // eq (172a)

	double C0 = 2.515516698;        // eq (172c)
	double C1 = 0.802853;           // eq (172d)
	double C2 = 0.010328;           // eq (172e)
	double D1 = 1.432788;           // eq (172f)
	double D2 = 0.189269;           // eq (172g)
	double D3 = 0.001308;           // eq (172h)

	double ksi = ( (C2*tx+C1)*tx + C0 )/ ( ((D3*tx + D2)*tx + D1)*tx + 1 );  // eq (172b)

	double I = ksi - tx;            // eq (172)

	return I;
}
void dl_p( double d[], double h[], double hts, double hrs,double hstd,double  hsrd,double f,double omega,double  p,
	double b0,double DN,int len_d ,double Ldp[2],double Ld50[2])
{
//dl_p Diffraction loss model not exceeded for p// of time according to P.452-16
//   function [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
//
//   This function computes the diffraction loss not exceeded for p// of time
//   as defined in ITU-R P.452-16 (Section 4.5.4)
//
//     Input parameters:
//     d       -   vector of distances di of the i-th profile point (km)
//     h       -   vector of heights hi of the i-th profile point (meters
//                 above mean sea level. Both vectors contain n+1 profile points
//     hts     -   transmitter antenna height in meters above sea level (i=0)
//     hrs     -   receiver antenna height in meters above sea level (i=n)
//     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
//     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
//     f       -   frequency expressed in GHz
//     omega   -   the fraction of the path over sea
//     p       -   percentage of time
//     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
//     DN      -   the average radio-refractive index lapse-rate through the
//                 lowest 1 km of the atmosphere. Note that DN is positive
//                 quantity in this procedure
//
//     Output parameters:
//     Ldp    -   diffraction loss for the general path not exceeded for p // of the time 
//                according to Section 4.2.4 of ITU-R P.452-16. 
//                Ldp(1) is for the horizontal polarization 
//                Ldp(2) is for the vertical polarization
//     Ld50   -   diffraction loss for p = 50//
//
//     Example:
//     [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
//       
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


//// Body of function

// Use the method in 4.2.3 to calculate diffractino loss Ld for effective 
// Earth radius ap = ae as given by equation (6a). Set median diffractino
// loss to Ldp50
	double ae = earth_rad_eff(DN);
	double ab = 6371*3;

	double ap = ae;

	dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega ,len_d,Ld50);

	if (p == 50)
	{
		Ldp[0] = Ld50[0];Ldp[1] = Ld50[1];
		return;
	}
	double Fi;
	double Ldb[2]={0};
	if (p < 50)
	{    
		// Use the method in 4.2.3 to calculate diffraction loss Ld for effective
		// Earth radius ap = abeta, as given in equation (6b). Set diffraction loss
		// not exceeded for beta0// time Ldb = Ld
    
		ap = ab;
    
		dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega ,len_d,Ldb);

		// Compute the interpolation factor Fi
    
		if (p > b0)
			Fi = inv_cum_norm(p/100) / inv_cum_norm(b0/100);   // eq (41a)
		else
			Fi = 1;
        
    
		// The diffraction loss Ldp not exceeded for p// of time is now given by
		for (int i=0;i<2;i++)
		Ldp[i] = Ld50[i] + Fi*(Ldb[i] - Ld50[i]);   // eq (42)
    
	}

	return;
}
double tl_anomalous(double dtot,double dlt, double dlr, double dct, double dcr, double dlm,double  hts, double hrs,double  hte, double hre, double hm,
	double theta_t,double  theta_r, double f, double p, double temp, double press, double omega, double ae, double b0)
{
//tl_anomalous Basic transmission loss due to anomalous propagation according to P.452-16
//   Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0)
//
//   This function computes the basic transmission loss occuring during
//   periods of anomalous propagation (ducting and layer reflection)
//   as defined in ITU-R P.452-16 (Section 4.4)
//
//     Input parameters:
//     dtot         -   Great-circle path distance (km)
//     dlt          -   interfering antenna horizon distance (km)
//     dlr          -   Interfered-with antenna horizon distance (km)
//     dct, dcr     -   Distance over land from the transmit and receive
//                      antennas tothe coast along the great-circle interference path (km).
//                      Set to zero for a terminal on a ship or sea platform
//     dlm          -   the longest continuous inland section of the great-circle path (km)
//     hts, hrs     -   Tx and Rx antenna heights aobe mean sea level amsl (m)
//     hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
//     hm           -   The terrain roughness parameter (m)
//     theta_t      -   Interfering antenna horizon elevation angle (mrad)
//     theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
//     f            -   frequency expressed in GHz
//     p            -   percentage of time
//     temp         -   Temperature (deg C)
//     press        -   Dry air pressure (hPa)
//     omega        -   fraction of the total path over water
//     ae           -   the median effective Earth radius (km)
//     b0           -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
//
//     Output parameters:
//     Lba    -   the basic transmission loss due to anomalous propagation
//               (ducting and layer reflection)
//
//     Example:
//     Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, b0)
//       
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version
//     v1    10NOV16     Ivica Stevanovic, OFCOM         added a line after eq. (54) to avoid division by zero
//     v2    13FEB17     Ivica Stevanovic, OFCOM         included lower limit for alpha and upper limit for mu2


//// Body of function

// empirical correction to account for the increasing attenuation with
// wavelength inducted propagation (47a)

	double Alf = 0;

	if (f < 0.5)
		Alf = 45.375 - 137.0*f + 92.5*f*f;

	// site-shielding diffraction losses for the interfering and interfered-with
	// stations (48)

	double theta_t1 = theta_t - 0.1*dlt;    // eq (48a)
	double theta_r1 = theta_r - 0.1*dlr;

	double Ast = 0;
	double Asr = 0;

	if (theta_t1 > 0)
		Ast = 20*log10(1 + 0.361*theta_t1*sqrt(f*dlt)) + 0.264*theta_t1*pow(f,(1.0/3));
	

	if (theta_r1 > 0)
		Asr = 20*log10(1 + 0.361*theta_r1*sqrt(f*dlr)) + 0.264*theta_r1*pow(f,1.0/3);
	

	// over-sea surface duct coupling correction for the interfering and
	// interfered-with stations (49) and (49a)

	double Act = 0;
	double Acr = 0;

	if (dct <= 5)
		if (dct <= dlt)
			if (omega >= 0.75)
				Act = -3*exp(-0.25*dct*dct)*(1+ tanh( 0.07*(50-hts) ));

	if (dcr <= 5)
		if (dcr <= dlr)
			if (omega >= 0.75)
				Acr = -3*exp(-0.25*dcr*dcr)*(1+ tanh( 0.07*(50-hrs) ));
			
	// specific attenuation (51)

	double gamma_d = 5e-5 * ae * pow(f,1.0/3);

	// angular distance (corrected where appropriate) (52-52a)

	theta_t1 = theta_t;
	theta_r1 = theta_r;

	if (theta_t> 0.1*dlt)
		theta_t1 = 0.1*dlt;
	if (theta_r > 0.1*dlr)
		theta_r1 = 0.1*dlr;

	double theta1 = 1e3*dtot/ae + theta_t1 + theta_r1;   
	double dI = (dtot - dlt - dlr)< 40?(dtot - dlt - dlr):40;   // eq (56a)
	double mu3 = 1;

	if (hm > 10)
		mu3 = exp( -4.6e-5 * (hm-10)*(43+6*dI) );  // eq (56)

	double tau = 1- exp(-(4.12e-4*pow(dlm,2.41)));       // eq (3a)
	double epsilon = 3.5;
	double alpha = -0.6 - epsilon*1e-9*pow(dtot,3.1)*tau;   // eq (55a)

	if (alpha < -3.4)
		alpha = -3.4;
	// correction for path geometry:

	double mu2 = pow( 500/ae * pow(dtot,2.0)/pow( sqrt(hte) + sqrt(hre),2.0 ),alpha);

	if (mu2 > 1)
		mu2 = 1;

	double beta = b0 * mu2 * mu3;      // eq (54)

	//beta = max(beta, eps);      // to avoid division by zero

	double Gamma = 1.076/pow(2.0058-log10(beta),1.012) * 
		exp( -( 9.51 - 4.8*log10(beta) + 0.198*pow(log10(beta),2.0))*1e-6*pow(dtot,1.13) );

	// time percentage variablity (cumulative distribution):

	double Ap = -12 + (1.2 + 3.7e-3*dtot)*log10(p/beta) + 12 * pow(p/beta,Gamma);  // eq (53)

	// time percentage and angular-distance dependent losses within the
	// anomalous propagation mechanism

	double Adp = gamma_d*theta1 + Ap;   // eq (50)

	// gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
	// whole path length

	// water vapor density

	double rho = 7.5 + 2.5*omega;

	double T = temp + 273.15;

	// compute specific attenuation due to dry air and water vapor:
	//[g_0, g_w] = p676d11_ga(f, press, rho, T);
	double g_0 = cal_gammao_p676_11(f,press,T,rho);
	double g_w = cal_gammaw_p676_11(f,press,T,rho);
	double	Ag = (g_0 + g_w) * dtot;  //(9)

	// total of fixed coupling losses (except for local clutter losses) between
	// the antennas and the anomalous propagation structure within the
	// atmosphere (47)

	double Af = 102.45 + 20*log10(f) + 20*log10(dlt + dlr) + Alf + Ast + Asr + Act + Acr;

	// total basic transmission loss occuring during periods of anomalaous
	// propagation 

	double Lba = Af + Adp + Ag;

	return Lba;
}
double tl_tropo(double dtot, double theta, double f, double p, double temp, double press, double N0, double Gt, double Gr )
//tl_tropo Basic transmission loss due to troposcatterer to P.452-16
//   Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
//
//   This function computes the basic transmission loss due to troposcatterer 
//   not exceeded for p// of time
//   as defined in ITU-R P.452-16 (Section 4.3)
//
//     Input parameters:
//     dtot    -   Great-circle path distance (km)
//     theta   -   Path angular distance (mrad)
//     f       -   frequency expressed in GHz
//     p       -   percentage of time
//     temp    -   Temperature (deg C)
//     press   -   Dry air pressure (hPa)
//     N0      -   path centre sea-level surface refractivity derived from Fig. 6
//     Gt,Gr   -   Antenna gain in the direction of the horizon along the
//                 great-circle interference path (dBi)
//
//     Output parameters:
//     Lbs    -   the basic transmission loss due to troposcatterer 
//                not exceeded for p// of time
//
//     Example:
//     Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
//       
//
//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


//// Body of function
{
	double T = temp + 273.15; 

	// Frequency dependent loss

	double Lf = 25*log10(f) - 2.5*pow(log10(f/2),2.0);    // eq (45a)

	// aperture to medium coupling loss (dB)

	double Lc = 0.051*exp(0.055*(Gt+Gr));             // eq (45b)

	// gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
	// whole path length

	double rho = 3.0;

	// compute specific attenuation due to dry air and water vapor:
	//[g_0, g_w] = p676d11_ga(f, press, rho, T);
		double g_0 = cal_gammao_p676_11(f,press,T,rho);
	double g_w = cal_gammaw_p676_11(f,press,T,rho);
	double Ag = (g_0 + g_w) * dtot;  //(9)

	// the basic transmission loss due to troposcatter not exceeded for any time
	// percentage p, below 50// is given

	double Lbs = 190 + Lf + 20*log10(dtot) + 0.573*theta - 0.15*N0 + Lc + Ag - 10.1*pow(-log10(p/50),0.7);

	return Lbs;
}
double tl_p452(double f, double p, double d[], double h[], int zone[],int worst_flag, double htg, double hrg,
	double phi_t, double phi_r,double Gt,double  Gr,int pol,double  dct,double dcr,double DN,double N0,
	double press,double temp,int len_d,double ha_t=0,double ha_r = 0,double dk_t = 0,double dk_r=0)
{
//tl_p452 basic transmission loss according to P.452-16
//   Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ha_t, ha_r, dk_t, dk_r )
//
//   This is the MAIN function that computes the basic transmission loss not exceeded for p// of time
//   as defined in ITU-R P.452-16 (Section 4.6). Other functions called from
//   this function are in ./src/ subfolder.
//
//     Input parameters:
//     f       -   Frequency (GHz)
//     p       -   Required time percentage for which the calculated basic
//                 transmission loss is not exceeded
//     d       -   vector of distances di of the i-th profile point (km)
//     h       -   vector of heights hi of the i-th profile point (meters
//                 above mean sea level. Both vectors contain n+1 profile points
//     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
//     htg     -   Tx Antenna center heigth above ground level (m)
//     hrg     -   Rx Antenna center heigth above ground level (m)
//     phi_t   -   Latitude of Tx station (degrees)
//     phi_r   -   Latitude of Rx station (degrees)
//     Gt, Gr  -   Antenna gain in the direction of the horizon along the
//                 great-circle interference path (dBi)
//     pol     -   polarization of the signal (1) horizontal, (2) vertical
//     dct     -   Distance over land from the transmit and receive
//     dcr         antennas to the coast along the great-circle interference path (km).
//                 Set to zero for a terminal on a ship or sea platform
//     DN      -   The average radio-refractive index lapse-rate through the
//                 lowest 1 km of the atmosphere (it is a positive quantity in this
//                 procedure) (N-units/km)
//     N0      -   The sea-level surface refractivity, is used only by the
//                 troposcatter model as a measure of location variability of the
//                 troposcatter mechanism. The correct values of DN and N0 are given by
//                 the path-centre values as derived from the appropriate
//                 maps (N-units)
//     press   -   Dry air pressure (hPa)
//     temp    -   Air temperature (degrees C)
//     ha_t    -   Clutter nominal height (m) at the Tx side
//     ha_r    -   Clutter nominal height (m) at the Rx side
//     dk_t    -   Clutter nominal distance (km) at the Tx side
//     dk_r    -   Clutter nominal distance (km) at the Rx side
//
//     Output parameters:
//     Lb     -   basic  transmission loss according to P.452-16
//
//     Example:
//     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp)
//     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ha_t, ha_r, dk_t, dk_r)

//     Rev   Date        Author                          Description
//     -------------------------------------------------------------------------------
//     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version
//     v1    10MAR16     Ivica Stevanovic, OFCOM         Introduced clutter losses
//     v2    15NOV16     Ivica Stevanovic, OFCOM         Corrected definition for Fj
//     v3    16NOV16     Ivica Stevanovic, OFCOM         Introduced the new version of ITU-R P.676-11 for gasseous attenuation calculation
//                                                       Corrected bug (epsr <-> sigma) in dl_se_ft
//                                                       Added a machine precision limit to avoid division by zero in tl_anomalous 
//                                                       Added a check for the number of points in the path profile 
//                                                       Added a check for the clutter loss nominal distances in cl_loss 
//     v4    10JAN17     Ivica Stevanovic, OFCOM         Corrected the Input parameters definition pol = 1(horizontal), pol = 2 (vertical) t
//     v5    13FEB17     Ivica Stevanovic, OFCOM         included lower limit for alpha and upper limit for mu2 in tl_anomalous

//   


// MATLAB Version 8.3.0.532 (R2014a) used in development of this code
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
//
// This function calls other functions that are placed in the ./src folder
// Test functions to verify/validate the current implementation are placed  in ./test folder

// s = pwd;
// if ~exist('p676d11_ga.m','file')
//     addpath([s '/src_452/'])
// end

// Read the input arguments 

// Compute the path profile parameters
// Path center latitude
	double phi_path = (phi_t + phi_r)/2;

	// Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
	int zone_r = 12;
	double dtm = longest_cont_dist(d,zone, zone_r,len_d);

	// Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
	zone_r = 2;
	double dlm = longest_cont_dist(d, zone, zone_r,len_d);

	// Compute b0
	double b0 = beta0(phi_path, dtm, dlm);

	double ae = earth_rad_eff(DN);
	double ab = 6371*3;
	// Compute the path fraction over see

	double omega = path_fraction(d, zone, 3,len_d);

	if (worst_flag)
		p = equiv_annual_percent(p, phi_path, omega);
	//     pp = pw2p(20,41,10,10);

	// Modify the path according to Section 4.5.4, Step 1 and compute clutter losses
	// only if not isempty ha_t and ha_r
	double *dc = new double[len_d];
	double *hc = new double[len_d];
	int *zonec = new int[len_d];
	double htrgc[2] = {0};
	double Ahtr[2]={0};
	int len_c = closs_corr(f, d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r,len_d,dc,hc,zonec,htrgc,Ahtr);

	//d = dc;
	//h = hc;
	//zone = zonec;
	//htg = htgc;
	//hrg = hrgc;
	double hstr[2],hstrd[2],htre[2],hm[1],dltr[2],theta_tr[2],theta[1];
	int pathtype[1];
	smooth_earth_heights(dc,hc,htrgc[0],htrgc[1],ae,f,len_c,hstr,hstrd,htre,hm,dltr,theta_tr,theta,pathtype);
	//[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

	double dtot = dc[len_c-1]-dc[0];

	//Tx and Rx antenna heights above mean sea level amsl (m)
	double hts = hc[0] + htrgc[0];
	double hrs = hc[len_c-1] + htrgc[1];

	// Effective Earth curvature Ce (km^-1)

	double Ce = 1/ae;

	// Wavelength in meters

	double lambda = 0.3/f;


	// Find the intermediate profile point with the highest slope of the line
	// from the transmitter to the point

	if (len_c<4)
		return -1234567;
	double *di = new double[len_c-2];
	double *hi=new double[len_c-2];
	double Stim = (hc[1] + 500*Ce*dc[1]*(dtot - dc[1]) - hts)/dc[1]; 
	for(int i = 0;i<len_c-2;i++)
	{
		di[i] = dc[i+1];
		hi[i] = hc[i+1];
		if(Stim<(hi[i] + 500*Ce*di[i]*(dtot - di[i]) - hts)/di[i])
			Stim = (hi[i] + 500*Ce*di[i]*(dtot - di[i]) - hts)/di[i];
	}

	// Calculate the slope of the line from transmitter to receiver assuming a
	// LoS path

	double Str = (hrs - hts)/dtot;                                         // Eq (15)

	// Calculate an interpolation factor Fj to take account of the path angular
	// distance (58)

	double THETA = 0.3;
	double KSI = 0.8;

	// changed the definition for Fj on 15DEC16.
	//Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (theta-THETA)/THETA) )
	double Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (Stim-Str)/THETA) );

	// Calculate an interpolation factor, Fk, to take account of the great
	// circle path distance:

	double dsw = 20;
	double kappa = 0.5;

	double Fk = 1.0 - 0.5*( 1.0 + tanh(3.0 * kappa * (dtot-dsw)/dsw) ); // eq (59)
	double Lb0pb[2];
	double Lbfsg = pl_los(dtot, f, p, b0, omega, temp, press, dltr[0], dltr[1],Lb0pb);
	//[Lbfsg, Lb0p, Lb0b] = pl_los(dtot, f, p, b0, omega, temp, press, dlt, dlr);

	double Ldp[2],Ld50[2];
	dl_p( dc, hc, hts, hrs, hstrd[0], hstrd[1], f, omega, p, b0, DN ,len_c,Ldp,Ld50);

	// The median basic transmission loss associated with diffraction Eq (43)
	double Lbd50[2];
	Lbd50[0] = Lbfsg + Ld50[0];Lbd50[1] = Lbfsg + Ld50[1];

	// The basic tranmission loss associated with diffraction not exceeded for
	// p// time Eq (44)

	double Lbd[2] = {Lb0pb[0] + Ldp[0],Lb0pb[0] + Ldp[1]};

	// A notional minimum basic transmission loss associated with LoS
	// propagation and over-sea sub-path diffraction

	double Lminb0p[2] = {Lb0pb[0] + (1-omega)*Ldp[0],Lb0pb[0] + (1-omega)*Ldp[1]};
	double Fi;
	if (p >= b0)
	{
	   Fi = inv_cum_norm(p/100)/inv_cum_norm(b0/100);   //eq (41a)
		for(int i=0;i<2;i++)
			Lminb0p[i] = Lbd50[i] + (Lb0pb[1] + (1-omega)*Ldp[i] - Lbd50[i])*Fi;   // eq (60)
   
	}

	// Calculate a notional minimum basic transmission loss associated with LoS
	// and transhorizon signal enhancements

	double eta = 2.5;

	double Lba = tl_anomalous(dtot, dltr[0], dltr[1], dct, dcr, dlm, hts, hrs, htre[0], htre[1], hm[0], theta_tr[0], theta_tr[1], f, p, temp, press, omega, ae, b0);

	double Lminbap = eta*log(exp(Lba/eta) + exp(Lb0pb[0]/eta));    // eq (61)

	// Calculate a notional basic transmission loss associated with diffraction
	// and LoS or ducting/layer reflection enhancements

	double Lbda[2] = {Lbd[0],Lbd[1]};
	double Lbam[2]={0};
	for (int i=0;i<2;i++)
	{
		if (Lminbap <= Lbd[i])
			Lbda[i] = Lminbap + (Lbd[i]-Lminbap)*Fk; 
		// Calculate a modified basic transmission loss, which takes diffraction and
		// LoS or ducting/layer-reflection enhancements into account

		Lbam[i] = Lbda[i] + (Lminb0p[i] - Lbda[i])*Fj;   // eq (63)
	}
	// Calculate the basic transmission loss due to troposcatter not exceeded
	// for any time percantage p 

	double Lbs = tl_tropo(dtot, theta[0], f, p, temp, press, N0, Gt, Gr );

	// Calculate the final transmission loss not exceeded for p// time
	double Lb_pol[2]={0};
	for (int i=0;i<2;i++)
	Lb_pol[i] = -5*log10(pow(10.0,-0.2*Lbs) + pow(10.0,-0.2*Lbam[i])) + Ahtr[0] + Ahtr[1];  // eq (64)

	double Lb = Lb_pol[pol-1];


	return Lb;
}
double cal_Lp_452(double lonlat_t[2],double lonlat_r[2],para_452 para)
{//// step1
	double Lp = 0;
	double freq = para.freq/1000; // GHz
	if (freq>50 || freq<0.1)
		return -1234567;
	double p_annual = para.p_annual;// p_annual,//

	double Lat_r = lonlat_r[1];// Rx纬度，单位度
	double Lon_r = lonlat_r[0];// Rx经度，单位度
	double htg = para.htg;//Tx天线地面上高度，单位m
	double hrg = para.hrg;//Rx地面上高度，单位m
	double Gt = para.Gt,Gr = para.Gr,press = para.press,temp = para.temp;
	int Pol = para.pol;
	if (press==-1)
		press = 1013.25;

	if (temp==-1)
		temp = 15;
	//// step2
	int worst_flag = para.worst_flag;

	//// clutter losses
	double ha_t = para.ha_t,dk_t = para.dk_t,ha_r = para.ha_r,dk_r = para.dk_r;

	//// step3
	double Delta_N = para.Delta_N,N_0 = para.N_0;

	//// step3bis
	double dct = para.dct,dcr = para.dcr;

	//// step4
	//MapPath = 'D:\Map\';
	// xmap = findstr(MapPath,'\');
	// MapPath = [MapPath(1:xmap(end-1)-1) '\Map\'];
	// MapPath='.\';
	string MapPath = "D:\\Map\\";
	int flag_cleansphere = para.flag_cleansphere,N_cal = para.N_cal;// 计算剖面时所需的点数
	int *ch = new int[N_cal+2];
	for(int i = 0;i<N_cal+2;i++)
		ch[i] = 2;

	double  Lat_t = lonlat_t[1];
	double  Lon_t = lonlat_t[0];
	double dis;
	double *D = new double[N_cal+2];
	double *h = new double[N_cal+2];
	D[0] = 0;h[0]=0;
	double fai[2]={0};
    if(flag_cleansphere)
    {
		dis=cal_distance(lonlat_t,lonlat_r,fai);
        for(int i = 1;i<N_cal+2;i++)
		{
			D[i] = D[i-1]+dis/(N_cal+1);
			h[i] = 0;
		}
	}
	else
       map(MapPath,Lon_t,Lat_t,Lon_r,Lat_r,N_cal,D,h,ch);

//     figure(1);
//     plot(D,h);pause;
//     kk = 1;
//     for ii = length(D):length(D)
//         [dlm,dtm] = dlmdtm(ii,D(1:ii),ch(1:ii,:));
//         PLoss(kk,pn) = basic_tx_loss3(prob,freq,Delta_N,htg,hrg,mean([Lat_t,Lat_r]),N_0,Gt(pn),Gr,omega,dtm,dlm,dct,dcr,0,D(1:ii),h,ha_t,dk_t,ha_r,dk_r,Pol,ch(1:ii,:),dtot);
// //         fprintf('在 //.2f km的位置上电波传播衰减为 //.2f dB\n',D(ii),PLoss(kk,pn));
//         kk = kk+1;
//     end
    Lp = tl_p452(freq, p_annual, D, h, ch,worst_flag, htg, hrg,Lat_t, Lat_r, Gt, Gr, Pol, dct, dcr, Delta_N, N_0, press,temp,N_cal+2, ha_t, ha_r, dk_t, dk_r);
	//Lp = 0;
	return Lp;
}

#endif