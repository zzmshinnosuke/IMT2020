#include "EMC_mod.h"

int deployment_imt(double pointx[],double pointy[],int stationnum,
	double Dp,double Dr,double radius=0,int flag=1,int flag_circle=1)
/*
  deployment_imt函数用于生成IMT基站，其中pointx、pointy为输出的坐标（直角）；stationnum为生成的基站数量；
  Dp为保护半径；保护半径内部不生成基站；Dr为生成半径；radius为蜂窝半径；flag为随机或规则标识位，1是随机、2是规则；
  flag_circle为生成圆形区域标识位，1表示圆形半径Dr，0表示方形边长Dr
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
//// 第五代移动通信系统干扰航空无线电导航系统集总仿真程序
// freq_IMT = 32000;

//// ****************************IMT系统参数**************************** ////
	double freq_IMT = para_IMT.freq;// IMT系统工作频点，MHz
	double Bt = para_IMT.B;// IMT系统以200MHz带宽为例进行分析
	double cell_dens[2];
	cell_dens[0] = para_IMT.cell_dens[0];//基站密度，单位个/km2
	cell_dens[1] = para_IMT.cell_dens[1];
	double h_IMT_BS = para_IMT.h;// IMT基站天线高度，米
	double Ptmax_IMTBS = para_IMT.Pmax;// IMT基站发射峰值功率为46dBm
	double Gtmax_IMTBS = para_IMT.Gmax;// IMT基站发射天线最大增益
	double Lossaver_BS = para_IMT.Lossaver_BS;// IMT发射信号平均功率衰减，网络负载率为50//，即平均功率为Ptmax-Lossaver

	double TDDloss_BS = para_IMT.TDDloss;// TDD活跃因子损耗
	double Lt_BS = para_IMT.L_Bs;// 干扰发射天线阵列损耗，dB
	//double fai_BS = para_IMT.phi_3db;// 干扰发射站天线水平方向半功率波束宽度
	double titl_BS = para_IMT.titl_BS;// 天线下倾角，单位度，对应郊区和市区
	double Ra = para_IMT.Ra;// 热点区域占城市面积百分比

	//// *****************导航系统参数******************* ////
	int nav_type = para_NAV.ant_type;
	double Br = para_NAV.B;// 单波束内接收机带宽，MHz
	double Lr = para_NAV.Lfeeder;// 接收机馈线损耗

	double Lcross = para_NAV.Lcross;// 接收天线与IMT天线极化隔离
	double Grmax = para_NAV.Gmax;// 天线主瓣增益
	// Gr_sidelobe = -7;
	// Efficiency = 0.75;
	double NF_NAV = para_NAV.NF;// 噪声系数
	double I_need = -120+NF_NAV;// 干扰保护标准，I0/N0 = -6
	// BIF_TT(:,1) = BIF_TT(:,1)/2;

	double h_NAV_sim = para_NAV.h;// 仿真所采用的飞机飞行高度
	
	double lonlat_NAV[2]; lonlat_NAV[0] = para_NAV.lonlat[0];
	 lonlat_NAV[1] = para_NAV.lonlat[1];// 飞机位置初始化

	double speed_plane = para_NAV.speed; // 飞机飞行速度，民航飞机通常为900km/h
	double f_NAV = para_NAV.freq;// 导航系统工作频点，MHz
	double ant_scan_ele = para_NAV.ant_scan_ele;// 天线仰角指向，单位度
	double ant_scan_azirpm = para_NAV.ant_scan_azirpm;// 天线方位角转速，单位 圈/分

	//// ******************干扰场景参数******************* ////
	double lonlat_city[2];
	lonlat_city[0]= para_env.lonlat[0];//待分析城市经纬度坐标初始化
	lonlat_city[1]= para_env.lonlat[1];
	double Suburban_R = para_env.Suburban_R;// 待分析范围半径，单位km
	double Urban_R = para_env.Urban_R;
	double t_percent = para_env.t_percent;// 时间百分比
	double p_percent = para_env.p_percent;// 地域百分比
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
	// 为每个IMT干扰源初始化基站天线水平指向，单位度，0度表示正北，90度表示正东，-90度表示正西
	}
	
	//para_IMTant = struct('phi',0.5,'theta',0.5,'phi_3dB',65,'theta_3dB',65,'A_m',30, 'SLA_v',30,'d_H',0.5*300/freq_IMT, ...
	//	'd_V',0.5*300/freq_IMT,'theta_tilt',10, 'phi_scan',0);
	//color_list = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];
	//load p676_paratable.mat;
	//// ********************** 电磁兼容分析******************** ////
	double D_max = para_sim.Dmax;// 飞机从距离城市D_max处开始飞向城市，远离城市D_max处不再计算，单位km
	double t_total = 2*D_max/speed_plane*3600;//仿真所需总时长
	// Dp_total = [1:2:125];
	double t_per = para_sim.t_per;//仿真时间步长，单位秒

	//double lonlat_NAV[2]={0};
	pxy_tmp[0] = D_max;pxy_tmp[1] = 0;
	xy2lonlat(pxy_tmp,lonlat_city,lonlat_NAV);// 计算飞机初始位置
	double dis;
	double alfaIMT,aziNAV0;
	double alfatmp[2]={0};
	dis = cal_distance(lonlat_city,lonlat_NAV,alfatmp);// 计算飞机初始角度
	alfaIMT = alfatmp[0];aziNAV0 = alfatmp[1];
	double FDR_IMT2NAV_IF = 0;
	//  ant_scan_azi = ant_scan_azirpm;//取出第nav_type类雷达水平转速数组
	//	ant_ele = ant_scan_ele{nav_type};//取出第nav_type类雷达垂直方向仰角
	//	h_NAV = h_NAV_sim{nav_type};
	//         for ele_no = 1:length(ant_ele)// 遍历第nav_type类雷达垂直仰角
	//         for hnav = 1:length(h_NAV)
	//*******t_sim = [0:t_per:t_total];*****//仿真时间数组
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
	
	for (double t = 0;t<=t_total;t = t+t_per,t_no = t_no+1)// 遍历仿真时间
	{
		aziNAV = aziNAV0+ t*ant_scan_azirpm*360/60;//计算当前时刻雷达方位角指向
		eleNAV = ant_scan_ele;// 取出当前时刻雷达仰角指向
		while (aziNAV>180 || aziNAV<-180) //保证雷达方位角范围在【-180,180】之间
		{
			if (aziNAV>180)
				aziNAV = aziNAV-360;
			if (aziNAV<-180)
				aziNAV = aziNAV+360;
		}
		Dp = D_max-speed_plane/3600*t; //计算当前时刻飞机距离城市中心位置，负表示已经飞过了城市上空
		//				D_sim(t_no) = Dp;
		pxy_tmp[0] = Dp;pxy_tmp[1] = 0;
		xy2lonlat(pxy_tmp,lonlat_city,lonlat_NAV);// 计算导航系统的经纬度
		elec_angle_BS(stationnum*para_sim.mont_num,titl_BS,escan,etilt);// 生成IMT基站天线电子倾角
		//	escan = reshape(escan,stationnum,mont_num);
		//	etilt = reshape(etilt,stationnum,mont_num);
		mkl = 0;
		I_sumtmp = 0;
		for (int k =0;k<stationnum;k++) //遍历IMT干扰基站
		{
			lonlat_IMTBS[0] = lon_IMTBS[k];lonlat_IMTBS[1] = lat_IMTBS[k];
			dis = cal_distance(lonlat_IMTBS,lonlat_NAV,alfatmp);
			alfaIMT = alfatmp[0];alfaNAV = alfatmp[1];
			//计算每个干扰站到受扰站之间的大圆距离
			
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
