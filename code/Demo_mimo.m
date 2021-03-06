%% 数据分析用代码
close all
clear
clc

%% 测试结论：
%在仅有AWGN信道下，本系统可以达到无误码率。修改：mainFun_mimo中112--LTEChanOut = awgn(LTEChanOut,snr)
%修改为：LTEChanOut = awgn(channelInput,snr);%只通过awgn

%在LTE EVA标准信道模型（matlab自带）下（默认模式），即使使用最小的导频间隔只能使误码率在0.4左右。。。。

%% 参数设置
% 
% fd=[5 10 15];%导频间隔
fd=15;
F_averge_ALL=[];
BER_ALL=[];
MSE_ALL=[];
ALSUM=2;%%%%%%%%%%%数据循环次数1000
%% 频偏循环
for i=1:length(fd)

snr=-5:5:15;

logo=['导频间隔',num2str(fd(i)),'个子载波'];
% legg{i}=logo;
%% SNR循环

for Snr=snr
     f_averge_m=ones(1,ALSUM);
     BER_m=ones(2,ALSUM);
     Hmse_m=ones(4,ALSUM);
    for cont=1:ALSUM
   
    [Hmse,f_averge,BER]=mainFun_mimo_fini(Snr,fd(i),1);%%%%%%%%%%%%%%%1:RAY+AWGN  2:LTEMIMO
    f_averge_m(:,cont)=f_averge;
    BER_m(:,cont)=BER;
%     Hmse;
    Hmse_m(:,cont)=reshape(Hmse,[],1);
    end
    f_averge_m(find(f_averge_m==404))=[];
    M_FA=mean(f_averge_m);
    timeerror=find(BER_m==404);
    BER_m(timeerror)=[];
    BER_m=reshape(BER_m,2,[]);
    M_BER=mean(BER_m,2);
    F_averge_ALL=[F_averge_ALL,M_FA];%%%%%%%%%%%%%%%所有频偏HZ数据
    Hmse_m(find(Hmse_m==404))=[];
    Hmse_m=reshape(Hmse_m,4,[]);
    M_Hmse_m=mean(Hmse_m,2);
    BER_ALL=[BER_ALL,M_BER];%%%%%%%%%%%%%%%所有误码率数据
    MSE_ALL=[MSE_ALL,M_Hmse_m]
end


%% 绘制BER图像
figure()
% figure(8)
% hold on
semilogy(snr,BER_ALL(:,(i)*length(snr)-length(snr)+1:(i)*length(snr)),'-p')
xlabel('SNR')
ylabel('BER')
title('BER')

grid on
xlim([-5 15]);
% ylim([1e-4 0.5]);
legend('天线1','天线2');
% hold off
%% 绘制频率估计偏差图像
% figure(8)
% hold on
% plot(snr, F_averge_ALL((i)*length(snr)-length(snr)+1:(i)*length(snr)),'-p')
% xlabel('SNR')
% ylabel(' F_averge_ALL')
% title(' F_averge_ALL')
% % xlim([-5 15]);
% grid on
% % legend('导频间隔15');
% hold off
% legend(legg);
% 
% % ((i)*length(snr)-length(snr)+1:(i)*length(snr))
figure()
% figure(8)
% hold on
semilogy(snr,MSE_ALL(:,(i)*length(snr)-length(snr)+1:(i)*length(snr)),'-p')
xlabel('SNR')
ylabel('MSE')
title('MSE')

grid on
xlim([-5 15]);
% ylim([1e-4 0.5]);
legend('H11','H12','H21','H22');

end

