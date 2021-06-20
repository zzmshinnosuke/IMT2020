

#include <math.h>
#include <time.h>
#include "ant_mod.h"
#include "prop_mod.h"
#include "P1546.h"
#include "p452.h"
#include "prop_mod_staticdata.h"
#include "other_calc_mod.h"

#ifdef CPPDLL_EXPORTS
#define CPP_EXPORTS __declspec(dllexport)
#else
#define CPP_EXPORTS __declspec(dllimport)
#endif

struct CPP_EXPORTS para_IMT5G
{	//5G系统参数
	//分别为：频率（MHz），带宽（MHz），高度（米），基站发射峰值功率（dBm），发射天线最大增益（dBi），平均基站激活损耗（dB），TDD活跃因子损耗（dB）
	//天线阵列损耗（dB），天线下倾角（度），Ra（%），Rb（%），基站密度（市区、郊区，个/km2）
	// para_5Gant为5G天线参数组
	double freq,B,h,Pmax,Gmax,Lossaver_BS,TDDloss,L_Bs,titl_BS,Ra,Rb,cell_dens[2];
	para_imt5Gant para_5Gant;

};
struct CPP_EXPORTS paraNAV_32G
{	//5G机载雷达参数组
	// 分别为带宽（MHz），馈线损耗（dB），极化损耗（dB），天线主瓣增益（dBi），噪声系数（dB）
	// 接收机干噪比（dB），天线高度（米），经纬度（度），飞行速度（km/h），频率（MHz），天线仰角指向（度），天线方位角转速（圈/分）
	double B,Lfeeder,Lcross,Gmax,NF,INR,h,lonlat[2],speed,freq,ant_scan_ele,ant_scan_azirpm;
	int ant_type;//表示天线类型标识，1对应F.699，2对应F.1245，3对应M.1466中第3类雷达
};
struct CPP_EXPORTS para_EMCenv
{//电磁兼容分析环境参数组
	//分别为城市经纬度（度），郊区半径（km），市区半径（km），时间百分比（%），地域百分比（%）
	double lonlat[2],Suburban_R,Urban_R,t_percent,p_percent;
};
struct CPP_EXPORTS para_simu
{//仿真参数
	//仿真最大距离（km），时间采样间隔（秒）
	double Dmax,t_per;
	int mont_num;
};

extern "C" CPP_EXPORTS typedef int (__stdcall *MyFun)(int,int); //如果不加__stdcall，调用的时候一直报错，读取内存异常

extern "C" CPP_EXPORTS double* IMT5G2NAV_32G_dyn(para_IMT5G para_IMT,paraNAV_32G para_NAV,para_EMCenv para_env,para_simu para_sim,MyFun myfun,double *I_IF_total_10);
