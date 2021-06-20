# 5GSimulation
5G仿真windows端应用程序
难点：c#调用c++中转接口的设计
wpf负责UI，c++负责仿真计算，中间件CLR负责两者之间的数据传递调用。

# Environments
win8.1（vmware虚拟机）  
vs2012   
.net4.0  
v110 运行库  

# Directory
emc_all_20180620 5G仿真计算模型  c++
EMCCLR 用于c#调用c++的公共运行库  c++
IMT2020 界面设计，以及仿真结果绘图 c#

# Requirments
WPFVisifire.Charts.dll  c#调用绘制折线图
mscvr110.dll和msvcp110.dll  在其他系统中运行可能需要，可以一起打包在一起，保证在多数的win系统中都可以运行。

# Run
（支持win7，win8，win10，winxp）


2.运行环境  安装.net4.0 ，运行直接把IMT2020下Release拷贝出来，将打包工具中的dll复制到运行目录下。
3.程序生成  首先生成 emc_all_20180620（右键生成），生产dll和lib在Release中，将lib复制到EMCCLR目录中，将lib和dll复制到IMT2020/Release目录下。然后EMCCLR（右键生成），生成的dll在Release中EMCCLR.dll。最后IMT2020，添加引用EMCCLR.dll和WPFVisifire.Charts.dll，执行Release，生产可执行文件。
4. 
   打包需要 emc_all_20180620计算模型需要的v110运行库
    c#需要引入的绘折线图库


