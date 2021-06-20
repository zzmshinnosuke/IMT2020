
#ifndef _IMT_ANTENNA_PATTERN
#define _IMT_ANTENNA_PATTERN

#include "math.h"
#include "stdio.h"

#define PI 3.14159265358979

double single_pattern(double phi,double theta,double G_Em,double Am,double SLA,double phi_3db,double theta_3db);

// dhspace is d_H/wave_length
// dvspace is d_V/wave_length
// 其他含义与建议书一致
double imt_antenna(double phi,double theta,int rownum,int colnum,double G_Em,double Am,double SLA,double phi_3db,double theta_3db,double theta_e,double phi_e,double dhspace=0.5,double dvspace=0.5)
{
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


#endif //_IMT_ANTENNA_PATTERN
