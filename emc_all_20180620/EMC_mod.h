

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
{	//5Gϵͳ����
	//�ֱ�Ϊ��Ƶ�ʣ�MHz��������MHz�����߶ȣ��ף�����վ�����ֵ���ʣ�dBm������������������棨dBi����ƽ����վ������ģ�dB����TDD��Ծ������ģ�dB��
	//����������ģ�dB������������ǣ��ȣ���Ra��%����Rb��%������վ�ܶȣ���������������/km2��
	// para_5GantΪ5G���߲�����
	double freq,B,h,Pmax,Gmax,Lossaver_BS,TDDloss,L_Bs,titl_BS,Ra,Rb,cell_dens[2];
	para_imt5Gant para_5Gant;

};
struct CPP_EXPORTS paraNAV_32G
{	//5G�����״������
	// �ֱ�Ϊ����MHz����������ģ�dB����������ģ�dB���������������棨dBi��������ϵ����dB��
	// ���ջ�����ȣ�dB�������߸߶ȣ��ף�����γ�ȣ��ȣ��������ٶȣ�km/h����Ƶ�ʣ�MHz������������ָ�򣨶ȣ������߷�λ��ת�٣�Ȧ/�֣�
	double B,Lfeeder,Lcross,Gmax,NF,INR,h,lonlat[2],speed,freq,ant_scan_ele,ant_scan_azirpm;
	int ant_type;//��ʾ�������ͱ�ʶ��1��ӦF.699��2��ӦF.1245��3��ӦM.1466�е�3���״�
};
struct CPP_EXPORTS para_EMCenv
{//��ż��ݷ�������������
	//�ֱ�Ϊ���о�γ�ȣ��ȣ��������뾶��km���������뾶��km����ʱ��ٷֱȣ�%��������ٷֱȣ�%��
	double lonlat[2],Suburban_R,Urban_R,t_percent,p_percent;
};
struct CPP_EXPORTS para_simu
{//�������
	//���������루km����ʱ�����������룩
	double Dmax,t_per;
	int mont_num;
};

extern "C" CPP_EXPORTS typedef int (__stdcall *MyFun)(int,int); //�������__stdcall�����õ�ʱ��һֱ������ȡ�ڴ��쳣

extern "C" CPP_EXPORTS double* IMT5G2NAV_32G_dyn(para_IMT5G para_IMT,paraNAV_32G para_NAV,para_EMCenv para_env,para_simu para_sim,MyFun myfun,double *I_IF_total_10);
