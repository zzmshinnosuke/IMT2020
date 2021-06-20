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
    /// ParaWindow.xaml 的交互逻辑
    /// </summary>
    public partial class ParaWindow : Window
    {
        public event RefreshLoadParaHandler RefreshLoadParaEvent;
        public ParaWindow()
        {
            InitializeComponent();
            getAllPara();
        }
        private void getAllPara()
        {
            List<Para> paraList=LoadPara.getAllPara();
            for (int e = 0; e < paraList.Count; e++)
            {
                Para para = new Para();
                para = paraList[paraList.Count-e-1];
                Frame frame = new Frame();
                ItemParaPage ipp = new ItemParaPage(para);
                ipp.RefreshParaEvent += new RefreshParaHandler(refreshPage);
                ipp.RefreshLoadParaEvent += new RefreshLoadParaHandler(refreshMainPage);
                frame.Content = ipp;
                this.listBox.Items.Add(frame);
            }
        }

        private void refreshPage()
        {
            this.listBox.Items.Clear();
            getAllPara();
        }

        private void refreshMainPage(Para para)
        {
            RefreshLoadParaEvent(para);
            this.Close();
        }

    }
}
