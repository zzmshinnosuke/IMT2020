// EMCCLR.h

#pragma once
#include "../emc_all_20180620/EMC_mod.h"
#include "../emc_all_20180620/ant_mod.h"
using namespace System;

namespace EMCCLR {
	public ref class  para_imt5Gant1
	{//5G天线参数组
		//
	public: double phi, theta;//天线水平角，天线俯仰角（不需要用户输入）
		int rownum,colnum;//天线阵列行数、列数
		double G_Em, Am, SLA, phi_3db, theta_3db, theta_e,phi_e;
		//天线最大增益（dBi）,Am，SLA，天线半功率波束宽度（水平，度）
		//天线半功率波束宽度（水平，度）,天线半功率波束宽度（垂直，度）
		//水平电子倾角（度），俯仰电子倾角（度）（不需要用户输入）
		double freq;//频率（MHz）（不需要用户输入）
		void setss( int  a)
		{
			int rownum=a;
		}

	};
	public ref class  para_IMT5G1
	{	//5G系统参数
		//分别为：频率（MHz），带宽（MHz），高度（米），基站发射峰值功率（dBm），发射天线最大增益（dBi），TDD活跃因子损耗（dB）
		//天线阵列损耗（dB），天线下倾角（度），Ra（%），Rb（%），基站密度（市区、郊区，个/km2）
		// para_5Gant为5G天线参数组
	public:	double freq,B,h,Pmax,Gmax,Lossaver_BS,TDDloss,L_Bs,titl_BS,Ra,Rb;
			double cell_dens1,cell_dens2;
  

	};
	public ref struct paraNAV_32G1
	{	//5G机载雷达参数组
		// 分别为带宽（MHz），馈线损耗（dB），极化损耗（dB），天线主瓣增益（dBi），噪声系数（dB）
		// 接收机干噪比（dB），天线高度（米），经纬度（度），飞行速度（km/h），频率（MHz），天线仰角指向（度），天线方位角转速（圈/分）
	public:	double B,Lfeeder,Lcross,Gmax,NF,INR,h;
		double lonlat1,lonlat2;
		double speed,freq,ant_scan_ele,ant_scan_azirpm;
		Int32 ant_type;//表示天线类型标识，1对应F.699，2对应F.1245，3对应M.1466中第3类雷达
	};
	public ref struct para_EMCenv1
	{//电磁兼容分析环境参数组
		//分别为城市经纬度（度），郊区半径（km），市区半径（km），时间百分比（%），地域百分比（%）
	public:	double lonlat1,lonlat2;
		double Suburban_R,Urban_R,t_percent,p_percent;
	};
	public ref struct para_simu1
	{//仿真参数
		//仿真最大距离（km），时间采样间隔（秒）
	public:	double Dmax,t_per;
			int mont_num;
	};
	#pragma comment(lib,"emc_all_20180620.lib")
	#pragma managed
	public delegate int reprogess(int , int);
	
	public ref class Class1
	{
		// TODO: 在此处添加此类的方法。
	public:
		Class1();
		~Class1();
		void IMT5G2NAV_32G_dyn1(para_imt5Gant1^ %para_imt5gant,para_IMT5G1^ %para_IMT,paraNAV_32G1^ %para_NAV,para_EMCenv1^ %para_env,para_simu1^ %para_sim,array<double>^ I_IF_total_10,reprogess^ rp);
	};
}
