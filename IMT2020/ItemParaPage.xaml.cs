using System;
using System.Collections.Generic;
using System.Linq;
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

namespace IMT2020
{
    public delegate void RefreshParaHandler();
    public delegate void RefreshLoadParaHandler(Para para);
    /// <summary>
    /// ItemParaPage.xaml 的交互逻辑
    /// </summary>
    public partial class ItemParaPage : Page
    {
        Para para;
        public event RefreshParaHandler RefreshParaEvent;
        public event RefreshLoadParaHandler RefreshLoadParaEvent;
        public ItemParaPage(Para para)
        {
            InitializeComponent();
            this.lbl_title.Content = para.name;
            this.lbl_time.Content = getTimeLong2String(long.Parse(para.time));
            this.para = para;
        }

        public ItemParaPage()
        {
            InitializeComponent();
        }

        private void btn_load_Click(object sender, RoutedEventArgs e)
        {
            RefreshLoadParaEvent(para);
        }

        private void btn_delete_Click(object sender, RoutedEventArgs e)
        {
            LoadPara.deletePara(para);
            RefreshParaEvent();
        }

        /// <summary>
        /// 时间戳转为字符串
        /// </summary>
        /// <param name="t">时间戳</param>
        /// <returns>年月日时分秒字符串</returns>
        public static String getTimeLong2String(long t)
        {
            return new DateTime(1970, 01, 01, 08, 00, 00).AddSeconds(t).ToString("yyyy-MM-dd HH:mm:ss");
        }
    }
}
