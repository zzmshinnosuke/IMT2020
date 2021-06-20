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
using System.Windows.Shapes;

namespace IMT2020
{
    /// <summary>
    /// ParaNameWindow.xaml 的交互逻辑
    /// </summary>
    public partial class ParaNameWindow : Window
    {
        string name = "";
        public ParaNameWindow()
        {
            InitializeComponent();
            this.txt_name.Text = "para"+getTimeDate2Long(DateTime.Now) + "";
        }

        public string getName()
        {
            return this.name;
        }

        private void btn_ok_Click(object sender, RoutedEventArgs e)
        {
            name = this.txt_name.Text;
            if (name == "")
            {
                MessageBox.Show("名称不能为空");
                return;
            }
            this.Close();
        }

        /// <summary>
        /// 将日期转化为时间戳
        /// </summary>
        /// <param name="dt">日期</param>
        /// <returns>时间戳</returns>
        private static long getTimeDate2Long(DateTime dt)
        {
            DateTime startTime = TimeZone.CurrentTimeZone.ToLocalTime(new System.DateTime(1970, 1, 1));
            long timeStamp = (long)(dt - startTime).TotalSeconds;
            return timeStamp;
        }
    }
}
