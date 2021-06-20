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
using System.Windows.Threading;
using Visifire.Charts;

namespace IMT2020
{
    /// <summary>
    /// 绘制折线图的窗口
    /// </summary>
    public partial class Graph : Window
    {
        DispatcherTimer timer = new DispatcherTimer();
        public Graph(double[] y)
        {
            InitializeComponent();
            ct.Children.Clear();
            MinMax minmax = getMinMax(y);
            //创建一个图标
            Chart chart = new Chart();
            //设置图标的宽度和高度
            //chart.Width = 926;
            //chart.Height = 320;
            //是否启用打印和保存图片
            chart.ToolBarEnabled = true;
            //设置图标的属性
            chart.ScrollingEnabled = false;//是否启用或禁用滚动
            chart.View3D = false;//3D效果显示
            Axis yAxis = new Axis();
            //设置图标中Y轴的最小值永远为0           
           // yAxis.AxisMinimum = minmax.min;
            yAxis.AxisMaximum = (int)(minmax.max + minmax.scale * 0.1);           

            //设置图表中Y轴的后缀          
            chart.AxesY.Add(yAxis);
            // 创建一个新的数据线。               
            DataSeries dataSeries = new DataSeries();
            dataSeries.RenderAs = RenderAs.Line;//折线图

            // 设置数据点              
            DataPoint dataPoint;
            for (int i = 0; i < y.Length; i++)
            {
                // 创建一个数据点的实例。                   
                dataPoint = new DataPoint();
                // 设置X轴点                    
                dataPoint.AxisXLabel = i + "";
                //设置Y轴点                   
                dataPoint.YValue = y[i];
                dataPoint.MarkerSize = 0;
                //设置数据点颜色                  
                dataPoint.Color = new SolidColorBrush(Colors.SkyBlue);
                //添加数据点                   
                dataSeries.DataPoints.Add(dataPoint);
            }
            // 添加数据线到数据序列。                
            chart.Series.Add(dataSeries);
            chart.BorderThickness = new Thickness(0);
            ct.Children.Add(chart);
        }

        /// <summary>
        ///保存图表中的最大/最小值，和他们之间的距离 
        /// </summary>
        struct MinMax
        {
            public double min;
            public double max;
            public double scale;
        }
        private MinMax getMinMax(double[] y)
        {
            double min = y[0], max = y[0];
            MinMax minmax = new MinMax();
            foreach (double yy in y)
            {
                if (yy < min) min = yy;
                if (yy > max) max = yy;
            }
            minmax.min = min;
            minmax.max = max;
            minmax.scale = max - min;
            return minmax;
        }
    }
}
