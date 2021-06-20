// 这是主 DLL 文件。

#include "stdafx.h"
using namespace System;
using namespace System::Runtime::InteropServices;
#include "EMCCLR.h"

using namespace EMCCLR;

EMCCLR::Class1::Class1(){}

EMCCLR::Class1::~Class1(){}

void EMCCLR::Class1::IMT5G2NAV_32G_dyn1(para_imt5Gant1^ %para_imt5gant,para_IMT5G1^ %para_IMT,paraNAV_32G1^ %para_NAV,para_EMCenv1^ %para_env,para_simu1^ %para_sim,array<double>^ I_IF_total_10,reprogess^ rp)
{
	para_imt5Gant p0={para_imt5gant->phi,para_imt5gant->theta,para_imt5gant->rownum,para_imt5gant->colnum,para_imt5gant->G_Em,para_imt5gant->Am
		,para_imt5gant->SLA,para_imt5gant->phi_3db,para_imt5gant->theta_3db,para_imt5gant->theta_e,para_imt5gant->phi_e,para_imt5gant->freq};
	para_IMT5G p1={para_IMT->freq,para_IMT->B,para_IMT->h,para_IMT->Pmax,para_IMT->Gmax,para_IMT->Lossaver_BS,para_IMT->TDDloss,para_IMT->L_Bs,
		para_IMT->titl_BS,para_IMT->Ra,para_IMT->Rb,{para_IMT->cell_dens1,para_IMT->cell_dens2},p0};
	paraNAV_32G p2={para_NAV->B,para_NAV->Lfeeder,para_NAV->Lcross,para_NAV->Gmax,para_NAV->NF,para_NAV->INR,para_NAV->h,{para_NAV->lonlat1,
		para_NAV->lonlat2},para_NAV->speed,para_NAV->freq,para_NAV->ant_scan_ele,para_NAV->ant_scan_azirpm,para_NAV->ant_type
	};
	para_EMCenv p3={{para_env->lonlat1,para_env->lonlat2},para_env->Suburban_R,para_env->Urban_R,para_env->t_percent,para_env->p_percent};
	para_simu p4={para_sim->Dmax,para_sim->t_per,para_sim->mont_num};

	double* d=new double[I_IF_total_10->Length];
	IntPtr stub_rp=Marshal::GetFunctionPointerForDelegate(rp);
	double* result=IMT5G2NAV_32G_dyn( p1, p2, p3, p4,static_cast<MyFun>(stub_rp.ToPointer()), d);
	
	for(int i=0;i<I_IF_total_10->Length;i++)
	{
		I_IF_total_10[i]=result[i];
	}
}
