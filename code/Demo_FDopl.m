%% ���ݷ����ô���
close all
clear
clc
%% ��������

fd=[2,17,34,56,238];%H��������չ
% fd=100;
F_averge_ALL=[];
BER_ALL=[];
TIME_NUM=[];
ALSUM=1000;%%%%%%%%%%%%ѭ������
%% Ƶƫѭ��
for i=1:length(fd)

snr=-5:5:20;

% logo=['��������չ',num2str(fd(i)),'HZ'];
logo=['��������չ1HZ'];
legg{i}=logo;
%% SNRѭ��
for Snr=snr
     f_averge_m=ones(1,ALSUM);
     BER_m=ones(1,ALSUM);
    for cont=1:ALSUM
   
    [f_averge,BER]=mainFun_2(Snr,fd(i));%%%%%%%%%%%%%%%1:mainFun_2  2:mainFun_siso_noconvolution
    f_averge_m(:,cont)=f_averge;
    BER_m(:,cont)=BER;
    end
    f_averge_m(find(f_averge_m==404))=[];
    M_FA=mean(f_averge_m);
    timeerror=find(BER_m==404);
    BER_m(timeerror)=[];
    M_BER=mean(BER_m,2)
    F_averge_ALL=[F_averge_ALL,M_FA];%%%%%%%%%%%%%%%����ƵƫHZ����
    TIME_NUM=[TIME_NUM,length(timeerror)]
    BER_ALL=[BER_ALL,M_BER];%%%%%%%%%%%%%%%��������������
    
end

%% ����BERͼ��
figure()
% figure(8)
% hold on
semilogy(snr,BER_ALL((i)*length(snr)-length(snr)+1:(i)*length(snr)),'-p')
xlabel('SNR')
ylabel('BER')
title('BER')

grid on
xlim([-5 20]);
% ylim([1e-4 0.5]);
legend(legg(i));
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

% ((i)*length(snr)-length(snr)+1:(i)*length(snr))
end

