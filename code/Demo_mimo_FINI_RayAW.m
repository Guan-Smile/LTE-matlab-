%% ���ݷ����ô���
close all
clear
clc

%% ���Խ��ۣ�
%�ڽ���AWGN�ŵ��£���ϵͳ���Դﵽ�������ʡ��޸ģ�mainFun_mimo��112--LTEChanOut = awgn(LTEChanOut,snr)
%�޸�Ϊ��LTEChanOut = awgn(channelInput,snr);%ֻͨ��awgn

%��LTE EVA��׼�ŵ�ģ�ͣ�matlab�Դ����£�Ĭ��ģʽ������ʹʹ����С�ĵ�Ƶ���ֻ��ʹ��������0.4���ҡ�������

%% ��������
% 
% fd=[5 10 15];%��Ƶ���
fd=15;
F_averge_ALL=[];
BER_ALL=[];
ALSUM=1;%%%%%%%%%%%����ѭ������
%% Ƶƫѭ��
for i=1:length(fd)

snr=10:5:15;

logo=['��Ƶ���',num2str(fd(i)),'�����ز�'];
% legg{i}=logo;
%% SNRѭ��

for Snr=snr
     f_averge_m=ones(1,ALSUM);
     BER_m=ones(2,ALSUM);
    for cont=1:ALSUM
   
    [f_averge,BER]=mainFun_mimo_fini(Snr,fd(i),1);
    f_averge_m(:,cont)=f_averge;
    BER_m(:,cont)=BER;
    end
    M_FA=mean(f_averge_m);
    M_BER=mean(BER_m,2)
    F_averge_ALL=[F_averge_ALL,M_FA];%%%%%%%%%%%%%%%����ƵƫHZ����
    
    BER_ALL=[BER_ALL,M_BER];%%%%%%%%%%%%%%%��������������
    
end


%% ����BERͼ��
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
legend('����1','����2');
% hold off
%% ����Ƶ�ʹ���ƫ��ͼ��
% figure(8)
% hold on
% plot(snr, F_averge_ALL((i)*length(snr)-length(snr)+1:(i)*length(snr)),'-p')
% xlabel('SNR')
% ylabel(' F_averge_ALL')
% title(' F_averge_ALL')
% % xlim([-5 15]);
% grid on
% % legend('��Ƶ���15');
% hold off
% legend(legg);
% 
% % ((i)*length(snr)-length(snr)+1:(i)*length(snr))
end

