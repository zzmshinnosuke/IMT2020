using EMCCLR;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IMT2020
{
    public class Para
    {
        public Para()
        {
            ;
        }
        public Para(string name1,string time1,para_imt5Gant1 para_imt5gant1, para_IMT5G1 para_imt5g1, paraNAV_32G1 paraNav_32g1, para_EMCenv1 para_emcenv1, para_simu1 para_simu1)
        {
            this.name = name1;
            this.time = time1;
            this.para_imt5gant = para_imt5gant1;
            this.para_imt5g = para_imt5g1;
            this.paraNav_32g = paraNav_32g1;
            this.para_emcenv = para_emcenv1;
            this.para_simu = para_simu1;
        }
        public string name;
        public string time;
        public para_imt5Gant1 para_imt5gant;
        public para_IMT5G1 para_imt5g;
        public paraNAV_32G1 paraNav_32g;
        public para_EMCenv1 para_emcenv;
        public para_simu1 para_simu;
    }
}
