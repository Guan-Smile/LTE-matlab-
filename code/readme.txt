

%%%SISO
使用方法：1.使用模式：Demo.m 运行即可。(可修改运行函数：[f_averge,BER]=mainFun_2(Snr,fd(i));%%%%%%%%%%%%%%%1:mainFun_2  2:mainFun_siso_noconvolution，分别对应有无卷积码的模式。)
               2.调试模式：main_ALL.m  运行即可。

实现了：时间同步、频率同步、信道估计（OFDM-15导频信息，SISO）、信道均衡（迫零均衡），未加卷积交织等，信道使用的是模拟EVA(Extended Vehicular A model,扩展车辆信道模型）的等效多径Rayleigh信道。

帧结构：（可调节，目前是：）

GUARD——PN——PN——PN——PN——OFDM——GUARD 

GUARD=20;
PN=255;
OFDM=1024+106(循环前缀)；
GUARD=20;
总帧长：datain_ALL ： 2190


%%MIMO
使用方法：1.使用模式：Demo_mimo.m 运行即可。
               2.调试模式：main_ALL_MIMO_test2.m  运行即可。




代码名及解释：

TimeLockFun时间同步函数
frequencLock频率同步函数
Canshu组帧函数
mse_G用于计算MSE
Demo系列：分为SISO,MIMO，用于最终结果分析，绘图

main_ALL：分为SISO,MIMO，全系统，用于调试测试。

mainFun：分为SISO,MIMO，main_ALL的函数版，用于demo中作函数调用。

