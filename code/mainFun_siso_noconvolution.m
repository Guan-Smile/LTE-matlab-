function [f_averge,BER]=mainFun(Snr,fD)%%%%
%% 1.����֡��

%��������
stp=14;%%%%��Ƶ�����15
N_subcarrier=1024;%���ز���128,1024
Npn=4;%%һ��N��pn����
global K

%% 1.1����ofdm����֡��
% ����ofdm���ƽ��ģ��
kk=1:stp+1:N_subcarrier;%(stp-1)/2:stp+1:1024;%��ǰ���15
mod = comm.OFDMModulator('NumGuardBandCarriers',[0;N_subcarrier-kk(end)],...
'PilotInputPort',true, ...%�Ƿ���뵼Ƶ
'FFTLength',N_subcarrier,...
'PilotCarrierIndices',[kk'], ...%���뵼Ƶ��λ������
'NumSymbols',1, ...
'CyclicPrefixLength',106,...%ѭ��ǰ׺�ĳ���
'InsertDCNull',false);  
modDim = info(mod); 
% showResourceMapping(mod)  %show
% ��������datain����Ƶ֡
% rng(22);
dain=randi([0 1],modDim.DataInputSize(1),1);% ��������֡
pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
dataIn = step(pskModulator,dain);% ����֡bpsk
% scatterplot(dataIn)%show
% Create a pilot signal that has the correct dimensions. 
pilotIn = complex(ones(modDim.PilotInputSize),ones(modDim.PilotInputSize)); % ���ɵ�Ƶ֡
% Apply OFDM modulation to the data and pilot signals. 
modData = step(mod,dataIn,pilotIn).*sqrt(N_subcarrier);%%%%%%%%%%ofdm������ɵ�����  ��128+106ѭ��ǰ׺��

%% 1.2�����PN���е�֡��
[PN,datain_ALL]=Canshu(Npn,modData);%%%%datain_ALL��ɵ�֡��

%% ����\���\��֯



%% ���ŵ�

% Rayleigh�ŵ�1
fs = 4e6;                                     % Hz
pathDelays = [0 3e-8 15e-8 31e-8 37e-8 71e-8 109e-8 173e-8 251e-8];    % sec
avgPathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];      % dB
% fD = 1;                                         % Hz
% % Rayleigh�ŵ�2
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-6];    % sec
% avgPathGains = [0 -10];      % dB
% fD = 1; 

% Create a Rayleigh channel System object
rchan1 = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',fD);%���ӻ�
% 'Visualization','Impulse and frequency responses')


after_Ray = rchan1(datain_ALL);
datain_ALL = awgn(after_Ray,Snr);%%%%%%%%aafterfm2:ͨ��Rayleigh�ŵ���aafterfm:��ͨ��Rayleigh�ŵ�

% datain_ALL = awgn(datain_ALL,Snr);%ֻͨ��awgn

%% ��Ƶƫ  fd
fd=100;%HZƵƫ
Rb=10e5;%%%%%%%%%��Դ��������
Ts=1./Rb;

cont=1:length(datain_ALL);%%%%%%%%%%ÿһ��Ķ���Ƶƫ
phase_pian = 2j*pi*fd.*Ts.*cont;%%��Ƶƫ���飺phase_pian = 2j*pi*fd.*Ts.*cont*0
datain_ALL=datain_ALL.*exp(phase_pian');

% scatterplot(datain_ALL);%show


%% ʱ��ͬ��
[Data_atertimelock,judg,Guard_atertimelock]=TimeLockFun(PN,datain_ALL,length(modData));
% if (judg==1)
%     !echo TIME LOCK success    
% else
%     !echo TIME LOCK Failed����
% end


%% Ƶƫ����
if ( Guard_atertimelock==404)
    BER=404
    f_averge=404;%%%%%%%%%%%��������
else
[f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

% [f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

%% ���

% ofdm���
demod = comm.OFDMDemodulator(mod);  
[dataOut, pilotOut] = step(demod,Data_atertFrelock);%%%%%%%%%%%ʹ��ʱ�滻ΪData_atertimelock!!!!!!Data_atertFrelock
dataOut = dataOut ./ sqrt(N_subcarrier);
pilotOut = pilotOut ./ sqrt(N_subcarrier);
%�Ƚ�
% isSame = (max(abs([dataIn(:) - dataOut(:); ...
%     pilotIn(:) - pilotOut(:)])) < 1e-10)

%% �ŵ�����

H_gu=(pilotOut)./(pilotIn); % lx: ���

 %����H����
 
 Hin=H_gu(:,1);
 for i=1:size(Hin)-1
    %for j=1:stp+1
    for k=1:stp
         ND(k+(i-1).*(stp+1))=Hin(i+1).*(k-1)./(stp+1)+Hin(i).*(stp+2-k)./(stp+1);
         co1 = (k)./(stp+1);
         co2 = (stp+1-k)./(stp+1);
         temp_lx(k+(i-1).*(stp))=Hin(i+1).*(k)./(stp+1)+Hin(i).*(stp+1-k)./(stp+1);

%           ND(j+(i-1)*(stp+1))=Hin(i);
    end    
 end
%  ND(stp+2+(size(Hin)-2).*(stp+1))=Hin(size(Hin));
temp_lx = temp_lx.';


%% �ŵ�����

ND_OUT= dataOut./ (temp_lx);%Ƶ��ԭ�źŹ���ֵ���Ƴ����ŵ�

%% ���ܷ���

%% ����������
pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);
data_FINI = step(pskDemodulator,ND_OUT);%%%%%%%%%%%���������о���BPSK�����dataOut(����ǰ)ND_OUT�����ƺ�

errorRate = comm.ErrorRate;
        errVec = errorRate(data_FINI,dain);   % data_FINI,dain;
%         ber=[ber;errVec(1)];
BER=errVec(1);
end