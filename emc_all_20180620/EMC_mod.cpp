#include "EMC_mod.h"

int deployment_imt(double pointx[],double pointy[],int stationnum,
	double Dp,double Dr,double radius=0,int flag=1,int flag_circle=1)
/*
  deployment_imt������������IMT��վ������pointx��pointyΪ��������ֱ꣨�ǣ���stationnumΪ���ɵĻ�վ������
  DpΪ�����뾶�������뾶�ڲ������ɻ�վ��DrΪ���ɰ뾶��radiusΪ���Ѱ뾶��flagΪ���������ʶλ��1�������2�ǹ���
  flag_circleΪ����Բ�������ʶλ��1��ʾԲ�ΰ뾶Dr��0��ʾ���α߳�Dr
*/
{
srand((unsigned int)time(NULL));
double x,y,R,Rinf,Rsup,rad,R2;
double randtmp = 0;
int num = 0;

if(flag==1)
{

	//pointxy = zeros(stationnum,2);

    while(num<stationnum)
	{
		randtmp = rand();
		x = Dr*(2*(randtmp/RAND_MAX)-1);
		randtmp = rand();
        y = Dr*(2*(randtmp/RAND_MAX)-1);
        if(x*x+y*y>= Dr*Dr || x*x+y*y<=Dp*Dp)
            continue;
        else
		{
            pointx[num] =x;pointy[num] =y;
			num = num+1;
		}
	}
	return num;
}
else if(flag==2)
{
	R = Dp+ Dr;
    Rinf = 0;Rsup= R;
    rad = radius;
	int k = 0,xi=1,yj=1;
    for(y = -Rsup,xi=1;y<=Rsup;y = y+0.75*sqrt(3.0)*rad,xi++)
	{
		for(x = Rsup,yj=1;x>=-Rsup;x = x-0.75*rad,yj++)
        {
			if((xi+yj)%2==0)
            {
				R2 = x*x+y*y;
				if(flag_circle)
				{
					if(R2<=Rsup*Rsup && R2>=Rinf*Rinf && R2>= Dp*Dp)
					{
						pointx[k]=x;
						pointy[k]=y;
						k = k+1;
					}
				}
				else
				{
					pointx[k]=x;
					pointy[k]=y;
					k = k+1;
				}
			}
		}
	}
	return k;

}
return 0;
}
void normrnd(double *out,int len,double mu=0,double sigma=1)
{
	srand((unsigned int)time(NULL));
	double x1,x2;
	for (int i = 0;i<len;i++)
	{
		x1 = rand()/RAND_MAX;
		x2 = rand()/RAND_MAX;
		out[i] = sqrt(-2*log(x1))*cos(x2*2*PI);
		out[i] = out[i]*sigma+mu;
	}
}
double normrnd(double mu=0,double sigma=1)
{
	srand((unsigned int)time(NULL));
	double x1,x2;
	x1 = rand()/RAND_MAX;
	x2 = rand()/RAND_MAX;
	double out = sqrt(-2*log(x1))*cos(x2*2*PI);
	out = out*sigma+mu;
	return out;
}

void raylrnd(double *out,double B,int len)
{
	if (B<0)
	{ out[0] = -1;return;}
	double x1,x2;
	for (int i=0;i<len;i++)
	{
		x1 = normrnd();x2 = normrnd();
		out[i] = sqrt(x1*x1+x2*x2)*B;
	}
}
void elec_angle_BS(int sta_num,double m_tilt,double *escan,double *etilt)
{
    double h_bs=6;    // antenna height of the base station
    double h_ue=1.5;  // antenna height of the user equipment
	double *dist_ue = new double[sta_num];
	//double *azim_ue = new double[sta_num];
    raylrnd(dist_ue,32,sta_num);     // the distance of the UE to base station is rayleigh distribution
    normrnd(escan,sta_num,0,30);   // the horizon angle is normal distribution
    //  the horizon angle is restricted  in [-60, 60] degree
    for (int i = 0;i<sta_num;i++)
	{
		if (escan[i]>60)
			escan[i] = 60;
		if (escan[i]<-60)
			escan[i] = -60;
		etilt[i] = atan( (h_bs-h_ue)/dist_ue[i] )*180/pi - m_tilt;
	}
}

double* IMT5G2NAV_32G_dyn(para_IMT5G para_IMT,paraNAV_32G para_NAV,para_EMCenv para_env,para_simu para_sim,MyFun myfun,double *I_IF_total_10)
{
//// ������ƶ�ͨ��ϵͳ���ź������ߵ絼��ϵͳ���ܷ������
// freq_IMT = 32000;

//// ****************************IMTϵͳ����**************************** ////
	double freq_IMT = para_IMT.freq;// IMTϵͳ����Ƶ�㣬MHz
	double Bt = para_IMT.B;// IMTϵͳ��200MHz����Ϊ�����з���
	double cell_dens[2];
	cell_dens[0] = para_IMT.cell_dens[0];//��վ�ܶȣ���λ��/km2
	cell_dens[1] = para_IMT.cell_dens[1];
	double h_IMT_BS = para_IMT.h;// IMT��վ���߸߶ȣ���
	double Ptmax_IMTBS = para_IMT.Pmax;// IMT��վ�����ֵ����Ϊ46dBm
	double Gtmax_IMTBS = para_IMT.Gmax;// IMT��վ���������������
	double Lossaver_BS = para_IMT.Lossaver_BS;// IMT�����ź�ƽ������˥�������縺����Ϊ50//����ƽ������ΪPtmax-Lossaver

	double TDDloss_BS = para_IMT.TDDloss;// TDD��Ծ�������
	double Lt_BS = para_IMT.L_Bs;// ���ŷ�������������ģ�dB
	//double fai_BS = para_IMT.phi_3db;// ���ŷ���վ����ˮƽ����빦�ʲ������
	double titl_BS = para_IMT.titl_BS;// ��������ǣ���λ�ȣ���Ӧ����������
	double Ra = para_IMT.Ra;// �ȵ�����ռ��������ٷֱ�

	//// *****************����ϵͳ����******************* ////
	int nav_type = para_NAV.ant_type;
	double Br = para_NAV.B;// �������ڽ��ջ�����MHz
	double Lr = para_NAV.Lfeeder;// ���ջ��������

	double Lcross = para_NAV.Lcross;// ����������IMT���߼�������
	double Grmax = para_NAV.Gmax;// ������������
	// Gr_sidelobe = -7;
	// Efficiency = 0.75;
	double NF_NAV = para_NAV.NF;// ����ϵ��
	double I_need = -120+NF_NAV;// ���ű�����׼��I0/N0 = -6
	// BIF_TT(:,1) = BIF_TT(:,1)/2;

	double h_NAV_sim = para_NAV.h;// ���������õķɻ����и߶�
	
	double lonlat_NAV[2]; lonlat_NAV[0] = para_NAV.lonlat[0];
	 lonlat_NAV[1] = para_NAV.lonlat[1];// �ɻ�λ�ó�ʼ��

	double speed_plane = para_NAV.speed; // �ɻ������ٶȣ��񺽷ɻ�ͨ��Ϊ900km/h
	double f_NAV = para_NAV.freq;// ����ϵͳ����Ƶ�㣬MHz
	double ant_scan_ele = para_NAV.ant_scan_ele;// ��������ָ�򣬵�λ��
	double ant_scan_azirpm = para_NAV.ant_scan_azirpm;// ���߷�λ��ת�٣���λ Ȧ/��

	//// ******************���ų�������******************* ////
	double lonlat_city[2];
	lonlat_city[0]= para_env.lonlat[0];//���������о�γ�������ʼ��
	lonlat_city[1]= para_env.lonlat[1];
	double Suburban_R = para_env.Suburban_R;// ��������Χ�뾶����λkm
	double Urban_R = para_env.Urban_R;
	double t_percent = para_env.t_percent;// ʱ��ٷֱ�
	double p_percent = para_env.p_percent;// ����ٷֱ�
	int stationnum_suburban = round(Ra*cell_dens[1]*(pi*Suburban_R*Suburban_R-pi*Urban_R*Urban_R));
	int stationnum_urban = round(Ra*cell_dens[0]*(pi*Urban_R*Urban_R));
	double *px_IMT = new double[stationnum_suburban+stationnum_urban];
	double *py_IMT = new double[stationnum_suburban+stationnum_urban];

	//int deployment_imt(double pointx[],double pointy[],int stationnum,
	//double Dp,double Dr,double radius=0,int flag=1,int flag_circle=1)
	int stationnum1 = deployment_imt(px_IMT,py_IMT,stationnum_suburban,Urban_R,Suburban_R);
	int stationnum2 = deployment_imt(px_IMT+stationnum_suburban,py_IMT+stationnum_suburban,
		stationnum_urban,Urban_R,Suburban_R);
	//[pxy_IMT_sub stationnum1]= deployment_imt(stationnum_suburban,Urban_R,Suburban_R);
	//[pxy_IMT_urb stationnum2]= deployment_imt(stationnum_urban,0,Urban_R);

	int stationnum = stationnum1+stationnum2;
	double *lon_IMTBS = new double[stationnum];
	double *lat_IMTBS = new double[stationnum];
	double lonlat_tmp[2] = {0};
	double pxy_tmp[2]={0};
	srand((unsigned int)time(NULL));
	double *azi_IMTBS = new double[stationnum];
	double rand_tmp = 0;
	para_imt5Gant para_5Gant = para_IMT.para_5Gant;
	para_5Gant.G_Em = Gtmax_IMTBS;
	para_5Gant.freq = freq_IMT;
			//para_IMTant = struct('phi',0.5,'theta',0.5,'phi_3dB',65,'theta_3dB',65,'A_m',30, 'SLA_v',30,'d_H',0.5*300/freq_IMT, ...
	//	'd_V',0.5*300/freq_IMT,'theta_tilt',10, 'phi_scan',0);
	for(int i=0;i<stationnum;i++)
	{
		pxy_tmp[0] = px_IMT[i];
		pxy_tmp[1] = py_IMT[i];
		xy2lonlat(pxy_tmp,lonlat_city,lonlat_tmp);
		lon_IMTBS[i] = lonlat_tmp[0];
		lat_IMTBS[i] = lonlat_tmp[1];
		rand_tmp = rand();
		azi_IMTBS[i] = rand_tmp/RAND_MAX*360-180;
	// Ϊÿ��IMT����Դ��ʼ����վ����ˮƽָ�򣬵�λ�ȣ�0�ȱ�ʾ������90�ȱ�ʾ������-90�ȱ�ʾ����
	}
	
	//para_IMTant = struct('phi',0.5,'theta',0.5,'phi_3dB',65,'theta_3dB',65,'A_m',30, 'SLA_v',30,'d_H',0.5*300/freq_IMT, ...
	//	'd_V',0.5*300/freq_IMT,'theta_tilt',10, 'phi_scan',0);
	//color_list = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];
	//load p676_paratable.mat;
	//// ********************** ��ż��ݷ���******************** ////
	double D_max = para_sim.Dmax;// �ɻ��Ӿ������D_max����ʼ������У�Զ�����D_max�����ټ��㣬��λkm
	double t_total = 2*D_max/speed_plane*3600;//����������ʱ��
	// Dp_total = [1:2:125];
	double t_per = para_sim.t_per;//����ʱ�䲽������λ��

	//double lonlat_NAV[2]={0};
	pxy_tmp[0] = D_max;pxy_tmp[1] = 0;
	xy2lonlat(pxy_tmp,lonlat_city,lonlat_NAV);// ����ɻ���ʼλ��
	double dis;
	double alfaIMT,aziNAV0;
	double alfatmp[2]={0};
	dis = cal_distance(lonlat_city,lonlat_NAV,alfatmp);// ����ɻ���ʼ�Ƕ�
	alfaIMT = alfatmp[0];aziNAV0 = alfatmp[1];
	double FDR_IMT2NAV_IF = 0;
	//  ant_scan_azi = ant_scan_azirpm;//ȡ����nav_type���״�ˮƽת������
	//	ant_ele = ant_scan_ele{nav_type};//ȡ����nav_type���״ﴹֱ��������
	//	h_NAV = h_NAV_sim{nav_type};
	//         for ele_no = 1:length(ant_ele)// ������nav_type���״ﴹֱ����
	//         for hnav = 1:length(h_NAV)
	//*******t_sim = [0:t_per:t_total];*****//����ʱ������
	double aziNAV,eleNAV,Dp;
	int t_no = 0;
	double *escan = new double[stationnum*para_sim.mont_num];
	double *etilt = new double[stationnum*para_sim.mont_num];
	double alfaNAV;
//	double alfatmp[2]={0};
	double dis_t;
	double lonlat_IMTBS[2]={0};
	double tmp_lonlat1[2] = {0};
	double tmp_lonlat2[2] = {0};
	double alfak,alfa2;
	double alfa_Gr;
	double Lp_10;
	double Gr;
	double phitmp;int mkl = 0;
	double Gt_IMT,I_IF0_10,I_IF_10;
	double mont_sum = 0,I_sumtmp = 0;
//	int t_no = 0;
	
	for (double t = 0;t<=t_total;t = t+t_per,t_no = t_no+1)// ��������ʱ��
	{
		aziNAV = aziNAV0+ t*ant_scan_azirpm*360/60;//���㵱ǰʱ���״﷽λ��ָ��
		eleNAV = ant_scan_ele;// ȡ����ǰʱ���״�����ָ��
		while (aziNAV>180 || aziNAV<-180) //��֤�״﷽λ�Ƿ�Χ�ڡ�-180,180��֮��
		{
			if (aziNAV>180)
				aziNAV = aziNAV-360;
			if (aziNAV<-180)
				aziNAV = aziNAV+360;
		}
		Dp = D_max-speed_plane/3600*t; //���㵱ǰʱ�̷ɻ������������λ�ã�����ʾ�Ѿ��ɹ��˳����Ͽ�
		//				D_sim(t_no) = Dp;
		pxy_tmp[0] = Dp;pxy_tmp[1] = 0;
		xy2lonlat(pxy_tmp,lonlat_city,lonlat_NAV);// ���㵼��ϵͳ�ľ�γ��
		elec_angle_BS(stationnum*para_sim.mont_num,titl_BS,escan,etilt);// ����IMT��վ���ߵ������
		//	escan = reshape(escan,stationnum,mont_num);
		//	etilt = reshape(etilt,stationnum,mont_num);
		mkl = 0;
		I_sumtmp = 0;
		for (int k =0;k<stationnum;k++) //����IMT���Ż�վ
		{
			lonlat_IMTBS[0] = lon_IMTBS[k];lonlat_IMTBS[1] = lat_IMTBS[k];
			dis = cal_distance(lonlat_IMTBS,lonlat_NAV,alfatmp);
			alfaIMT = alfatmp[0];alfaNAV = alfatmp[1];
			//����ÿ������վ������վ֮��Ĵ�Բ����
			
			cal_thita(dis*1000,h_IMT_BS,h_NAV_sim,alfatmp);
			alfak = alfatmp[0];alfa2 = alfatmp[1];
			tmp_lonlat1[0] = alfaNAV;tmp_lonlat1[1] =-1*alfa2;
			tmp_lonlat2[0] = aziNAV;tmp_lonlat2[1] = eleNAV;
            dis_t = cal_distance(tmp_lonlat1,tmp_lonlat2,alfatmp);
			alfa_Gr = 180*dis_t/6371/pi;
			if (nav_type==2)
				Gr = Gfai_F1245_G0(freq_IMT,Grmax,alfa_Gr,0.7);
			else if(nav_type == 3)
				Gr = Antenna_1466_No3(alfaNAV,-alfa2);
			else if(nav_type ==1)
				Gr = Gfai_F699_G0(freq_IMT,Grmax,alfa_Gr,0.7);
			
			phitmp = azi_IMTBS[k]+alfaIMT;
			while (phitmp>180 || phitmp<-180)
			{
				if (phitmp>180)
					phitmp = phitmp-360;
				
				if (phitmp<-180)
					phitmp = phitmp+360;
			}

			para_5Gant.phi = phitmp;
			para_5Gant.theta = 90+alfak-titl_BS;

			Lp_10 = cal_Lp_E2A_5G(freq_IMT,dis,p_percent,h_IMT_BS,h_NAV_sim);
			mont_sum =0;
			for (int mk=0;mk<para_sim.mont_num;mk++)
			{
				
				para_5Gant.theta_e =etilt[mkl];
				para_5Gant.phi_e = escan[mkl];
				mkl = mkl+1;
				Gt_IMT = imt_5G_antenna(para_5Gant);
                            
				I_IF0_10 = Ptmax_IMTBS-Lossaver_BS+Gt_IMT-Lt_BS-Lp_10-FDR_IMT2NAV_IF-Lcross+Gr-Lr;
				mont_sum = mont_sum+pow(10.0,I_IF0_10/10.0);
				//I_IF_10_mont(k,mk) = 10^(I_IF0_10/10);
			}
			I_IF_10 = mont_sum/para_sim.mont_num;
			I_sumtmp = I_sumtmp+I_IF_10;
		}
		I_IF_total_10[t_no] = 10*log10(I_sumtmp)-10*log10(Bt);
		myfun(t_total,t);
	}
	return I_IF_total_10;
}
