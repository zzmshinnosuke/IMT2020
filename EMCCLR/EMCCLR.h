// EMCCLR.h

#pragma once
#include "../emc_all_20180620/EMC_mod.h"
#include "../emc_all_20180620/ant_mod.h"
using namespace System;

namespace EMCCLR {
	public ref class  para_imt5Gant1
	{//5G���߲�����
		//
	public: double phi, theta;//����ˮƽ�ǣ����߸����ǣ�����Ҫ�û����룩
		int rownum,colnum;//������������������
		double G_Em, Am, SLA, phi_3db, theta_3db, theta_e,phi_e;
		//����������棨dBi��,Am��SLA�����߰빦�ʲ�����ȣ�ˮƽ���ȣ�
		//���߰빦�ʲ�����ȣ�ˮƽ���ȣ�,���߰빦�ʲ�����ȣ���ֱ���ȣ�
		//ˮƽ������ǣ��ȣ�������������ǣ��ȣ�������Ҫ�û����룩
		double freq;//Ƶ�ʣ�MHz��������Ҫ�û����룩
		void setss( int  a)
		{
			int rownum=a;
		}

	};
	public ref class  para_IMT5G1
	{	//5Gϵͳ����
		//�ֱ�Ϊ��Ƶ�ʣ�MHz��������MHz�����߶ȣ��ף�����վ�����ֵ���ʣ�dBm������������������棨dBi����TDD��Ծ������ģ�dB��
		//����������ģ�dB������������ǣ��ȣ���Ra��%����Rb��%������վ�ܶȣ���������������/km2��
		// para_5GantΪ5G���߲�����
	public:	double freq,B,h,Pmax,Gmax,Lossaver_BS,TDDloss,L_Bs,titl_BS,Ra,Rb;
			double cell_dens1,cell_dens2;
  

	};
	public ref struct paraNAV_32G1
	{	//5G�����״������
		// �ֱ�Ϊ����MHz����������ģ�dB����������ģ�dB���������������棨dBi��������ϵ����dB��
		// ���ջ�����ȣ�dB�������߸߶ȣ��ף�����γ�ȣ��ȣ��������ٶȣ�km/h����Ƶ�ʣ�MHz������������ָ�򣨶ȣ������߷�λ��ת�٣�Ȧ/�֣�
	public:	double B,Lfeeder,Lcross,Gmax,NF,INR,h;
		double lonlat1,lonlat2;
		double speed,freq,ant_scan_ele,ant_scan_azirpm;
		Int32 ant_type;//��ʾ�������ͱ�ʶ��1��ӦF.699��2��ӦF.1245��3��ӦM.1466�е�3���״�
	};
	public ref struct para_EMCenv1
	{//��ż��ݷ�������������
		//�ֱ�Ϊ���о�γ�ȣ��ȣ��������뾶��km���������뾶��km����ʱ��ٷֱȣ�%��������ٷֱȣ�%��
	public:	double lonlat1,lonlat2;
		double Suburban_R,Urban_R,t_percent,p_percent;
	};
	public ref struct para_simu1
	{//�������
		//���������루km����ʱ�����������룩
	public:	double Dmax,t_per;
			int mont_num;
	};
	#pragma comment(lib,"emc_all_20180620.lib")
	#pragma managed
	public delegate int reprogess(int , int);
	
	public ref class Class1
	{
		// TODO: �ڴ˴���Ӵ���ķ�����
	public:
		Class1();
		~Class1();
		void IMT5G2NAV_32G_dyn1(para_imt5Gant1^ %para_imt5gant,para_IMT5G1^ %para_IMT,paraNAV_32G1^ %para_NAV,para_EMCenv1^ %para_env,para_simu1^ %para_sim,array<double>^ I_IF_total_10,reprogess^ rp);
	};
}
