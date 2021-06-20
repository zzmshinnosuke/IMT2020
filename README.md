# 5GSimulation
1.开发环境  win8.1（vmware虚拟机）  vs2012 .net4.0   
2.运行环境  安装.net4.0 ，运行直接把IMT2020下Release拷贝出来，将打包工具中的dll复制到运行目录下。（支持win7，win8，win10，winxp）
3.程序生成  首先生成 emc_all_20180620（右键生成），生产dll和lib在Release中，将lib复制到EMCCLR目录中，将lib和dll复制到IMT2020/Release目录下。然后EMCCLR（右键生成），生成的dll在Release中EMCCLR.dll。最后IMT2020，添加引用EMCCLR.dll和WPFVisifire.Charts.dll，执行Release，生产可执行文件。
4. emc_all_20180620 5G仿真计算模型
   EMCCLR 用于c#调用c++的公共运行库
   IMT2020 界面设计，以及仿真结果绘图
   打包需要 emc_all_20180620计算模型需要的v110运行库
   WPFVisifire.Charts.dll c#需要引入的绘折线图库

难点：c#调用c++中转接口的设计
