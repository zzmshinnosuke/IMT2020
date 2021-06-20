using EMCCLR;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Xml;

namespace IMT2020
{
    class LoadPara 
    {
        static String Path = Environment.CurrentDirectory;
        static String filePath = Path + "/Para.xml";

        public static bool savePara(Para para)
        {
            try
            {
                if (!File.Exists(filePath))  //如果不存在
                {
                    CreateXml();
                }
                //DeleteUser(nameStr);
                XmlDocument doc = new XmlDocument();
                doc.Load(filePath);
                XmlNode node = doc.SelectSingleNode("Para");
                XmlElement paranode = doc.CreateElement("P");
                paranode.SetAttribute("name", para.name);
                paranode.SetAttribute("time", getTimeDate2Long(DateTime.Now) + "");

                XmlElement imt5GantNode = doc.CreateElement("imt5Gant");
                imt5GantNode.SetAttribute("phi",para.para_imt5gant.phi+"");
                imt5GantNode.SetAttribute("theta", para.para_imt5gant.theta + "");
                imt5GantNode.SetAttribute("rownum", para.para_imt5gant.rownum + "");
                imt5GantNode.SetAttribute("colnum", para.para_imt5gant.colnum + "");
                imt5GantNode.SetAttribute("G_Em", para.para_imt5gant.G_Em + "");
                imt5GantNode.SetAttribute("Am", para.para_imt5gant.Am + "");
                imt5GantNode.SetAttribute("SLA", para.para_imt5gant.SLA + "");
                imt5GantNode.SetAttribute("phi_3db", para.para_imt5gant.phi_3db + "");
                imt5GantNode.SetAttribute("theta_3db", para.para_imt5gant.theta_3db + "");
                imt5GantNode.SetAttribute("theta_e", para.para_imt5gant.theta_e + "");
                imt5GantNode.SetAttribute("phi_e", para.para_imt5gant.phi_e + "");
                imt5GantNode.SetAttribute("freq", para.para_imt5gant.freq + "");

                XmlElement IMT5GNode = doc.CreateElement("IMT5G");
                IMT5GNode.SetAttribute("freq", para.para_imt5g.freq + "");
                IMT5GNode.SetAttribute("B", para.para_imt5g.B + "");
                IMT5GNode.SetAttribute("h", para.para_imt5g.h + "");
                IMT5GNode.SetAttribute("Pmax", para.para_imt5g.Pmax + "");
                IMT5GNode.SetAttribute("Gmax", para.para_imt5g.Gmax + "");
                IMT5GNode.SetAttribute("Lossaver_BS", para.para_imt5g.Lossaver_BS + "");
                IMT5GNode.SetAttribute("TDDloss", para.para_imt5g.TDDloss + "");
                IMT5GNode.SetAttribute("L_Bs", para.para_imt5g.L_Bs + "");
                IMT5GNode.SetAttribute("titl_BS", para.para_imt5g.titl_BS + "");
                IMT5GNode.SetAttribute("Ra", para.para_imt5g.Ra + "");
                IMT5GNode.SetAttribute("Rb", para.para_imt5g.Rb + "");
                IMT5GNode.SetAttribute("cell_dens1", para.para_imt5g.cell_dens1 + "");
                IMT5GNode.SetAttribute("cell_dens2", para.para_imt5g.cell_dens2 + "");

                XmlElement NAV_32GNode = doc.CreateElement("NAV_32G");
                NAV_32GNode.SetAttribute("B",para.paraNav_32g.B+"");
                NAV_32GNode.SetAttribute("Lfeeder",para.paraNav_32g.Lfeeder+"");
                NAV_32GNode.SetAttribute("Lcross",para.paraNav_32g.Lcross+"");
                NAV_32GNode.SetAttribute("Gmax",para.paraNav_32g.Gmax+"");
                NAV_32GNode.SetAttribute("NF",para.paraNav_32g.NF+"");
                NAV_32GNode.SetAttribute("INR",para.paraNav_32g.INR+"");
                NAV_32GNode.SetAttribute("h",para.paraNav_32g.h+"");
                NAV_32GNode.SetAttribute("lonlat1",para.paraNav_32g.lonlat1+"");
                NAV_32GNode.SetAttribute("lonlat2",para.paraNav_32g.lonlat2+"");
                NAV_32GNode.SetAttribute("speed",para.paraNav_32g.speed+"");
                NAV_32GNode.SetAttribute("freq",para.paraNav_32g.freq+"");
                NAV_32GNode.SetAttribute("ant_scan_ele",para.paraNav_32g.ant_scan_ele+"");
                NAV_32GNode.SetAttribute("ant_scan_azirpm",para.paraNav_32g.ant_scan_azirpm+"");
                NAV_32GNode.SetAttribute("ant_type",para.paraNav_32g.ant_type+"");

                XmlElement EMCenvNode = doc.CreateElement("EMCenv");
                EMCenvNode.SetAttribute("lonlat1",para.para_emcenv.lonlat1+"");
                EMCenvNode.SetAttribute("lonlat2",para.para_emcenv.lonlat2+"");
                EMCenvNode.SetAttribute("Suburban_R",para.para_emcenv.Suburban_R+"");
                EMCenvNode.SetAttribute("Urban_R",para.para_emcenv.Urban_R+"");
                EMCenvNode.SetAttribute("t_percent",para.para_emcenv.t_percent+"");
                EMCenvNode.SetAttribute("p_percent", para.para_emcenv.p_percent + "");

                XmlElement simuNode = doc.CreateElement("simu");
                simuNode.SetAttribute("Dmax", para.para_simu.Dmax + "");
                simuNode.SetAttribute("t_per", para.para_simu.t_per + "");
                simuNode.SetAttribute("mont_num", para.para_simu.mont_num + "");

                paranode.AppendChild(imt5GantNode);
                paranode.AppendChild(IMT5GNode);
                paranode.AppendChild(NAV_32GNode);
                paranode.AppendChild(EMCenvNode);
                paranode.AppendChild(simuNode);
                node.AppendChild(paranode);
                doc.Save(filePath);

            }
            catch (Exception e)
            {
                System.Console.Out.WriteLine(e.Message.ToString());
                return false;
            }
            return true;
        }

        public static List<Para> getAllPara()
        {
            List<Para> paraList = new List<Para>();
            try
            {
                if (!Directory.Exists(Path))
                {
                    new DirectoryInfo(Path).Create();

                }
                if (!File.Exists(filePath))  //如果不存在
                {
                    CreateXml();
                    return null;
                }
                else   //存在
                {
                    XmlDocument doc = new XmlDocument();
                    doc.Load(filePath);
                    XmlNodeList nodelist = doc.SelectSingleNode("Para").ChildNodes;
                    foreach (XmlNode xnode in nodelist)
                    {
                        if (xnode.Name == "P")
                        {
                                                                                                              
                            //天线参数，共十二个参数
                            XmlNode imt5Gantnode = xnode.SelectSingleNode("imt5Gant");
                            para_imt5Gant1 para_imt5gant = new para_imt5Gant1();
                            para_imt5gant.phi = double.Parse(imt5Gantnode.Attributes["phi"].Value.ToString());
                            para_imt5gant.theta = double.Parse(imt5Gantnode.Attributes["theta"].Value.ToString());
                            para_imt5gant.rownum = int.Parse(imt5Gantnode.Attributes["rownum"].Value.ToString());
                            para_imt5gant.colnum = int.Parse(imt5Gantnode.Attributes["colnum"].Value.ToString());
                            para_imt5gant.G_Em = double.Parse(imt5Gantnode.Attributes["G_Em"].Value.ToString());
                            para_imt5gant.Am = double.Parse(imt5Gantnode.Attributes["Am"].Value.ToString());
                            para_imt5gant.SLA = double.Parse(imt5Gantnode.Attributes["SLA"].Value.ToString());
                            para_imt5gant.phi_3db = double.Parse(imt5Gantnode.Attributes["phi_3db"].Value.ToString());
                            para_imt5gant.theta_3db = double.Parse(imt5Gantnode.Attributes["theta_3db"].Value.ToString());
                            para_imt5gant.theta_e = double.Parse(imt5Gantnode.Attributes["theta_e"].Value.ToString());
                            para_imt5gant.phi_e = double.Parse(imt5Gantnode.Attributes["phi_e"].Value.ToString());
                            para_imt5gant.freq = double.Parse(imt5Gantnode.Attributes["freq"].Value.ToString());

                            XmlNode IMT5Gnode = xnode.SelectSingleNode("IMT5G");
                            para_IMT5G1 para_imt5g = new para_IMT5G1();
                            //系统参数，共十三个参数
                            para_imt5g.freq = double.Parse(IMT5Gnode.Attributes["freq"].Value.ToString());
                            string temp = para_imt5g.freq.ToString();
                            System.Console.Out.WriteLine(temp);
                            para_imt5g.B = double.Parse(IMT5Gnode.Attributes["B"].Value.ToString());
                            para_imt5g.h = double.Parse(IMT5Gnode.Attributes["h"].Value.ToString());
                            para_imt5g.Pmax = double.Parse(IMT5Gnode.Attributes["Pmax"].Value.ToString());
                            para_imt5g.Gmax = double.Parse(IMT5Gnode.Attributes["Gmax"].Value.ToString());  //和天线参数中的天线最大增益是一个
                            para_imt5g.Lossaver_BS = double.Parse(IMT5Gnode.Attributes["Lossaver_BS"].Value.ToString());
                            para_imt5g.TDDloss = double.Parse(IMT5Gnode.Attributes["TDDloss"].Value.ToString());
                            para_imt5g.L_Bs = double.Parse(IMT5Gnode.Attributes["L_Bs"].Value.ToString());
                            para_imt5g.titl_BS = double.Parse(IMT5Gnode.Attributes["titl_BS"].Value.ToString());
                            para_imt5g.Ra = double.Parse(IMT5Gnode.Attributes["Ra"].Value.ToString());
                            para_imt5g.Rb = double.Parse(IMT5Gnode.Attributes["Rb"].Value.ToString());
                            para_imt5g.cell_dens1 = double.Parse(IMT5Gnode.Attributes["cell_dens1"].Value.ToString());
                            para_imt5g.cell_dens2 = double.Parse(IMT5Gnode.Attributes["cell_dens2"].Value.ToString());

                            XmlNode NAV_32Gnode = xnode.SelectSingleNode("NAV_32G");
                            paraNAV_32G1 paraNav_32g = new paraNAV_32G1();
                            //机载雷达参数，共十四个参数
                            paraNav_32g.B = double.Parse(NAV_32Gnode.Attributes["B"].Value.ToString());
                            paraNav_32g.Lfeeder = double.Parse(NAV_32Gnode.Attributes["Lfeeder"].Value.ToString());
                            paraNav_32g.Lcross = double.Parse(NAV_32Gnode.Attributes["Lcross"].Value.ToString());
                            paraNav_32g.Gmax = double.Parse(NAV_32Gnode.Attributes["Gmax"].Value.ToString());
                            paraNav_32g.NF = double.Parse(NAV_32Gnode.Attributes["NF"].Value.ToString());
                            paraNav_32g.INR = double.Parse(NAV_32Gnode.Attributes["INR"].Value.ToString());
                            paraNav_32g.h = double.Parse(NAV_32Gnode.Attributes["h"].Value.ToString());
                            paraNav_32g.lonlat1 = double.Parse(NAV_32Gnode.Attributes["lonlat1"].Value.ToString());
                            paraNav_32g.lonlat2 = double.Parse(NAV_32Gnode.Attributes["lonlat2"].Value.ToString());
                            paraNav_32g.speed = double.Parse(NAV_32Gnode.Attributes["speed"].Value.ToString());
                            paraNav_32g.freq = double.Parse(NAV_32Gnode.Attributes["freq"].Value.ToString());
                            paraNav_32g.ant_scan_ele = double.Parse(NAV_32Gnode.Attributes["ant_scan_ele"].Value.ToString());
                            paraNav_32g.ant_scan_azirpm = double.Parse(NAV_32Gnode.Attributes["ant_scan_azirpm"].Value.ToString());
                            paraNav_32g.ant_type = int.Parse(NAV_32Gnode.Attributes["ant_type"].Value.ToString());

                            XmlNode EMCenvnode = xnode.SelectSingleNode("EMCenv");
                            para_EMCenv1 para_emcenv = new para_EMCenv1();
                            //场景参数,共六个参数
                            para_emcenv.lonlat1 = double.Parse(EMCenvnode.Attributes["lonlat1"].Value.ToString());
                            para_emcenv.lonlat2 = double.Parse(EMCenvnode.Attributes["lonlat2"].Value.ToString());
                            para_emcenv.Suburban_R = double.Parse(EMCenvnode.Attributes["Suburban_R"].Value.ToString());
                            para_emcenv.Urban_R = double.Parse(EMCenvnode.Attributes["Urban_R"].Value.ToString());
                            para_emcenv.t_percent = double.Parse(EMCenvnode.Attributes["t_percent"].Value.ToString());
                            para_emcenv.p_percent = double.Parse(EMCenvnode.Attributes["p_percent"].Value.ToString());

                            XmlNode simunode = xnode.SelectSingleNode("simu");
                            para_simu1 para_simu = new para_simu1();
                            //仿真参数,共3个参数
                            para_simu.Dmax = double.Parse(simunode.Attributes["Dmax"].Value.ToString());
                            para_simu.t_per = double.Parse(simunode.Attributes["t_per"].Value.ToString());
                            para_simu.mont_num = int.Parse(simunode.Attributes["mont_num"].Value.ToString());

                            Para para = new Para(xnode.Attributes["name"].Value.ToString(), xnode.Attributes["time"].Value.ToString(), para_imt5gant, para_imt5g, paraNav_32g, para_emcenv, para_simu);
                            paraList.Add(para);
                        }
                    }
                }
            }
            catch (Exception e)
            {
                System.Windows.MessageBox.Show(e.Message.ToString());
                return null;
            }
            return paraList;
        }
        public static bool deletePara(Para para)
        {
            try
            {
                if (!File.Exists(filePath))  //如果不存在
                {
                    CreateXml();
                }

                XmlDocument doc = new XmlDocument();
                doc.Load(filePath);

                XmlNodeList nodelist = doc.SelectSingleNode("Para").ChildNodes;
                foreach (XmlNode xnode in nodelist)
                {
                    if (xnode.Name == "P" && xnode.Attributes["time"].Value.ToString() == para.time)
                    {
                        xnode.ParentNode.RemoveChild(xnode);
                    }
                }
                doc.Save(filePath);
            }
            catch (Exception)
            {
                return false;
            }
            return true;
        }
        


        /// <summary>
        /// 创建xml配置文件
        /// </summary>
        /// <returns></returns>
        private static bool CreateXml()
        {
            XmlDocument xmlDoc = new XmlDocument();
            XmlNode xmlnode = xmlDoc.CreateXmlDeclaration("1.0", "utf-8", "");
            xmlDoc.AppendChild(xmlnode);
            //创建根节点  
            XmlNode root = xmlDoc.CreateElement("Para");
            xmlDoc.AppendChild(root);
            try
            {
                xmlDoc.Save(filePath);
            }
            catch (Exception e)
            {
                //显示错误信息  
                Console.WriteLine(e.Message);
                return false;
            }
            return true;
        }

        /// <summary>
        /// 将日期转化为时间戳
        /// </summary>
        /// <param name="dt">日期</param>
        /// <returns>时间戳</returns>
        public static long getTimeDate2Long(DateTime dt)
        {
            DateTime startTime = TimeZone.CurrentTimeZone.ToLocalTime(new System.DateTime(1970, 1, 1));
            long timeStamp = (long)(dt - startTime).TotalSeconds;
            return timeStamp;
        }

    }
}
