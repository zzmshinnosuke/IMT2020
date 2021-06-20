using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Runtime.InteropServices;
using EMCCLR;
using System.Threading;

namespace IMT2020
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window
    {
        Thread objThread;  //定义一个线程，保证计算过程中，页面不会卡死
        para_imt5Gant1 para_imt5gant = new para_imt5Gant1();
        para_IMT5G1 para_imt5g = new para_IMT5G1();
        paraNAV_32G1 paraNav_32g = new paraNAV_32G1();
        para_EMCenv1 para_emcenv = new para_EMCenv1();
        para_simu1 para_simu = new para_simu1();
        public MainWindow()
        {
            InitializeComponent();
        }

        private void setPara()
        {
            this.Dispatcher.Invoke(new Action(delegate
           {
               //天线参数，共十二个参数
               para_imt5gant.phi = double.Parse(this.txt_phi.Text);
               para_imt5gant.theta = double.Parse(this.txt_theta.Text);
               para_imt5gant.rownum = int.Parse(this.txt_rownum.Text);
               para_imt5gant.colnum = int.Parse(this.txt_colnum.Text);
               para_imt5gant.G_Em = double.Parse(this.txt_G_Em.Text);
               para_imt5gant.Am = double.Parse(this.txt_Am.Text);
               para_imt5gant.SLA = double.Parse(this.txt_SLA.Text);
               para_imt5gant.phi_3db = double.Parse(this.txt_phi_3db.Text);
               para_imt5gant.theta_3db = double.Parse(this.txt_theta_3db.Text);
               para_imt5gant.theta_e = double.Parse(this.txt_theta_e.Text);
               para_imt5gant.phi_e = double.Parse(this.txt_phi_e.Text);
               para_imt5gant.freq = double.Parse(this.txt_freq_imt.Text);

               //系统参数，共十三个参数
               para_imt5g.freq = double.Parse(this.txt_freq.Text);
               para_imt5g.B = double.Parse(this.txt_B.Text);
               para_imt5g.h = double.Parse(this.txt_h.Text);
               para_imt5g.Pmax = double.Parse(this.txt_Pmax.Text);
               para_imt5g.Gmax = double.Parse(this.txt_G_Em.Text);  //和天线参数中的天线最大增益是一个
               para_imt5g.Lossaver_BS = double.Parse(this.txt_Lossaver_Bs.Text);
               para_imt5g.TDDloss = double.Parse(this.txt_TDDloss.Text);
               para_imt5g.L_Bs = double.Parse(this.txt_L_BS.Text);
               para_imt5g.titl_BS = double.Parse(this.txt_titl_BS.Text);
               para_imt5g.Ra = double.Parse(this.txt_Ra.Text) / 100;
               para_imt5g.Rb = double.Parse(this.txt_Rb.Text) / 100;
               para_imt5g.cell_dens1 = double.Parse(this.txt_cell_dens_1.Text);
               para_imt5g.cell_dens2 = double.Parse(this.txt_cell_dens_2.Text);

               //机载雷达参数，共十四个参数
               paraNav_32g.B = double.Parse(this.txt_B_NAN_32G.Text);
               paraNav_32g.Lfeeder = double.Parse(this.txt_Lfeeder.Text);
               paraNav_32g.Lcross = double.Parse(this.txt_Lcross.Text);
               paraNav_32g.Gmax = double.Parse(this.txt_Gmax_NAN_32G.Text);
               paraNav_32g.NF = double.Parse(this.txt_NF.Text);
               paraNav_32g.INR = double.Parse(this.txt_INR.Text);
               paraNav_32g.h = double.Parse(this.txt_h_NAN_32G.Text);
               paraNav_32g.lonlat1 = double.Parse(this.txt_lon_NAN_32G.Text);
               paraNav_32g.lonlat2 = double.Parse(this.txt_lat_NAN_32G.Text);
               paraNav_32g.speed = double.Parse(this.txt_speed.Text);
               paraNav_32g.freq = double.Parse(this.txt_freq_NAN_32G.Text);
               paraNav_32g.ant_scan_ele = double.Parse(this.txt_ant_scan_ele.Text);
               paraNav_32g.ant_scan_azirpm = double.Parse(this.txt_ant_scan_azirpm.Text);
               paraNav_32g.ant_type = this.cb_ant_type.SelectedIndex + 1;

               //场景参数,共六个参数
               para_emcenv.lonlat1 = double.Parse(this.txt_lon_EMCenv.Text);
               para_emcenv.lonlat2 = double.Parse(this.txt_lat_EMCenv.Text);
               para_emcenv.Suburban_R = double.Parse(this.txt_Suburban_R.Text);
               para_emcenv.Urban_R = double.Parse(this.txt_Urban_R.Text);
               para_emcenv.t_percent = double.Parse(this.txt_t_percent.Text);
               para_emcenv.p_percent = double.Parse(this.txt_p_percent.Text);

               //仿真参数,共3个参数
               para_simu.Dmax = double.Parse(this.txt_Dmax.Text);
               para_simu.t_per = double.Parse(this.txt_t_per.Text);
               para_simu.mont_num = int.Parse(this.txt_mont_num.Text);
           }));   
        }

        private void getPara(Para para_temp)
        {
            //天线参数，共十二个参数
            this.txt_phi.Text=para_temp.para_imt5gant.phi.ToString();
            this.txt_theta.Text=para_temp.para_imt5gant.theta.ToString();
            this.txt_rownum.Text=para_temp.para_imt5gant.rownum.ToString();
            this.txt_colnum.Text=para_temp.para_imt5gant.colnum.ToString();
            this.txt_G_Em.Text=para_temp.para_imt5gant.G_Em.ToString();
            this.txt_Am.Text=para_temp.para_imt5gant.Am.ToString();
            this.txt_SLA.Text=para_temp.para_imt5gant.SLA.ToString();
            this.txt_phi_3db.Text=para_temp.para_imt5gant.phi_3db.ToString();
            this.txt_theta_3db.Text=para_temp.para_imt5gant.theta_3db.ToString();
            this.txt_theta_e.Text=para_temp.para_imt5gant.theta_e.ToString();
            this.txt_phi_e.Text=para_temp.para_imt5gant.phi_e.ToString();
            this.txt_freq_imt.Text = para_temp.para_imt5gant.freq.ToString();

            //系统参数，共十三个参数
            this.txt_freq.Text = para_temp.para_imt5g.freq.ToString();
            this.txt_B.Text = para_temp.para_imt5g.B.ToString();
            this.txt_h.Text = para_temp.para_imt5g.h.ToString();
            this.txt_Pmax.Text = para_temp.para_imt5g.Pmax.ToString();
            this.txt_G_Em.Text = para_temp.para_imt5g.Gmax.ToString();
            this.txt_Lossaver_Bs.Text = para_temp.para_imt5g.Lossaver_BS.ToString();
            this.txt_TDDloss.Text = para_temp.para_imt5g.TDDloss.ToString();
            this.txt_L_BS.Text = para_temp.para_imt5g.L_Bs.ToString();
            this.txt_titl_BS.Text = para_temp.para_imt5g.titl_BS.ToString();
            this.txt_Ra.Text = (para_temp.para_imt5g.Ra*100).ToString();
            this.txt_Rb.Text = (para_temp.para_imt5g.Rb*100).ToString();
            this.txt_cell_dens_1.Text = para_temp.para_imt5g.cell_dens1.ToString();
            this.txt_cell_dens_2.Text = para_temp.para_imt5g.cell_dens2.ToString();

            //机载雷达参数，共十四个参数
            this.txt_B_NAN_32G.Text=para_temp.paraNav_32g.B.ToString();
            this.txt_Lfeeder.Text=para_temp.paraNav_32g.Lfeeder.ToString();
            this.txt_Lcross.Text = para_temp.paraNav_32g.Lcross.ToString();
            this.txt_Gmax_NAN_32G.Text = para_temp.paraNav_32g.Gmax.ToString();
            this.txt_NF.Text= para_temp.paraNav_32g.NF.ToString();
            this.txt_INR.Text= para_temp.paraNav_32g.INR.ToString();
            this.txt_h_NAN_32G.Text= para_temp.paraNav_32g.h.ToString();
            this.txt_lon_NAN_32G.Text = para_temp.paraNav_32g.lonlat1.ToString();
            this.txt_lat_NAN_32G.Text = para_temp.paraNav_32g.lonlat2.ToString();
            this.txt_speed.Text = para_temp.paraNav_32g.speed.ToString();
            this.txt_freq_NAN_32G.Text = para_temp.paraNav_32g.freq.ToString();
            this.txt_ant_scan_ele.Text = para_temp.paraNav_32g.ant_scan_ele.ToString();
            this.txt_ant_scan_azirpm.Text = para_temp.paraNav_32g.ant_scan_azirpm.ToString();
            this.cb_ant_type.SelectedIndex=para_temp.paraNav_32g.ant_type-1;

            //场景参数,共六个参数
            this.txt_lon_EMCenv.Text = para_temp.para_emcenv.lonlat1.ToString();
            this.txt_lat_EMCenv.Text = para_temp.para_emcenv.lonlat2.ToString();
            this.txt_Suburban_R.Text = para_temp.para_emcenv.Suburban_R.ToString();
            this.txt_Urban_R.Text = para_temp.para_emcenv.Urban_R.ToString();
            this.txt_t_percent.Text = para_temp.para_emcenv.t_percent.ToString();
            this.txt_p_percent.Text = para_temp.para_emcenv.p_percent.ToString();

            //仿真参数,共3个参数
            this.txt_Dmax.Text = para_temp.para_simu.Dmax.ToString();
            this.txt_t_per.Text=para_temp.para_simu.t_per.ToString();
            this.txt_mont_num.Text = para_temp.para_simu.mont_num.ToString();
        }

        public void compute()
        {
            double[] d = new double[(int)(2 * para_simu.Dmax / paraNav_32g.speed * 3600 / para_simu.t_per) + 1];
            Class1 cl = new Class1();
            cl.IMT5G2NAV_32G_dyn1(ref para_imt5gant, ref para_imt5g, ref paraNav_32g, ref para_emcenv, ref para_simu, d, reProgess);
            this.Dispatcher.Invoke(new Action(delegate
            {
                Graph graph = new Graph(d);
                graph.ShowDialog();
                this.btn_ok.IsEnabled = true;
            }));
        }

        public int reProgess(int all, int cur)
        {
            this.Dispatcher.Invoke(new Action(delegate
          {
              xbar.Maximum = all;
              xbar.Value = cur;
              this.lbl_percent.Content = (int)((double)cur / all * 100) + "%";
          }));
            return cur;
        }

        private void btn_ok_Click(object sender, RoutedEventArgs e)
        {
            setPara();
            objThread = new Thread(new ThreadStart(delegate
            {
                compute();
            }));
            objThread.Start();
            this.btn_ok.IsEnabled = false;
        }

        private void btn_cancle_Click(object sender, RoutedEventArgs e)
        {
            MessageBoxResult result = MessageBox.Show("是否取消仿真?", "", MessageBoxButton.YesNo);
            if (result == MessageBoxResult.Yes)
            {
                this.btn_ok.IsEnabled = true;
                xbar.Value = 0;
                this.lbl_percent.Content = "0%";
                objThread.Abort();
            }
        }

        private void btn_end_Click(object sender, RoutedEventArgs e)
        {
            MessageBoxResult result = MessageBox.Show("是否退出程序?", "", MessageBoxButton.YesNo);
            if (result == MessageBoxResult.Yes)
            {
                close();
            }
        }

        private void close()
        {
            System.Diagnostics.Process.GetCurrentProcess().Kill();
            this.Close();
        }

        private void Window_Closed_1(object sender, EventArgs e)
        {
            close();
        }

        private void btn_load_para_Click(object sender, RoutedEventArgs e)
        {
            LoadPara.getAllPara();
            ParaWindow paraWindow = new ParaWindow();
            paraWindow.RefreshLoadParaEvent += new RefreshLoadParaHandler(refreshMainPage);
            paraWindow.ShowDialog();
        }

        private void refreshMainPage(Para para)
        {
            getPara(para);
        }

        private void btn_save_para_Click(object sender, RoutedEventArgs e)
        {
            ParaNameWindow pnw = new ParaNameWindow();
            pnw.ShowDialog();
            setPara();
            if (pnw.getName() != "")
            {
                Para p = new Para(pnw.getName(), "", para_imt5gant, para_imt5g, paraNav_32g, para_emcenv, para_simu);
                LoadPara.savePara(p);
            }
        }
    }
}
