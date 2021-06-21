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
（支持win7，win8，win10，winxp） 需要安装.net4.0  
解决方案中存在3个项目，需要按对应顺序先生成动态链接库    
1. emc_all_20180620 右键“重新生成”，生成的链接库保存在./Release中  
2. EMCCLR 右键“重新生成”，生成的链接库保存在./Release中  
3. IMT2020 直接启动，生成的exe保存在了./Release中（IMT2020添加了引用EMCCLR.dll和WPFVisifire.Charts.dll）  
4. 最后生成的./Release 可以拷贝执行（或者打包），在win7或winxp上运行可以复制mscvr110.dll和msvcp110.dll两个链接库。



