function [f_averge,BER]=mainFun_2(Snr,fD)%%%%

% clear
% clc
% close all
%% ��

%% 1.����֡��

%��������
stp=14;%%%%��Ƶ�����15
N_subcarrier=1024;%���ز���128,1024
Npn=4;%%һ��N��pn����
global K
snr=Snr;

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
dain=randi([0 1],modDim.DataInputSize(1)./2-6-8,1);% ��������֡  -8:crc^8   -4:crc^4
%CRC����
crcGen = comm.CRCGenerator([8 7 6 4 2 0]);
crcDet = comm.CRCDetector([8 7 6 4 2 0]);%%'z4+z3+z2+z+1'

encData = step(crcGen,dain);                % Append CRC bits

%��������
hConEnc = comm.ConvolutionalEncoder('TerminationMethod','Terminated');
hDec = comm.ViterbiDecoder('InputFormat','Hard','TerminationMethod','Terminated');

hError = comm.ErrorRate;
dain_JUAN = step(hConEnc, encData);%%%%%%%encData

% bu_zero=zeros(10-rem(length(dain_JUAN),10),1);
% dain_JUAN=[dain_JUAN;bu_zero];
%%��֯
JIA0_NUM=fD;%%%%%%%%%%%%%%%������

interleaver = comm.MatrixInterleaver('NumRows',JIA0_NUM,'NumColumns', length(dain_JUAN)/JIA0_NUM);
deinterleaver = comm.MatrixDeinterleaver('NumRows',JIA0_NUM,'NumColumns', length(dain_JUAN)/JIA0_NUM);

% dain_JIAO = dain_JUAN;
dain_JIAO = interleaver(dain_JUAN);

pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
dataIn = step(pskModulator,dain_JIAO);% ����֡bpsk
% scatterplot(dataIn)%show
% Create a pilot signal that has the correct dimensions. 
pilotIn = complex(ones(modDim.PilotInputSize),ones(modDim.PilotInputSize)); % ���ɵ�Ƶ֡
% Apply OFDM modulation to the data and pilot signals. 
modData = step(mod,dataIn,pilotIn).*sqrt(N_subcarrier);%%%%%%%%%%ofdm������ɵ�����  ��128+106ѭ��ǰ׺��

%% 1.2�����PN���е�֡��
[PN,datain_ALL]=Canshu(Npn,modData);%%%%datain_ALL��ɵ�֡��

%% ����



%% ���ŵ�

% Rayleigh�ŵ�1
fs = 4e6;                                     % Hz
pathDelays = [0 3e-8 15e-8 31e-8 37e-8 71e-8 109e-8 173e-8 251e-8];    % sec
avgPathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];      % dB
% fD = 1;                                         % Hz
% Rayleigh�ŵ�2
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-6];    % sec
% avgPathGains = [0 -10];      % dB
% fD = 1; 

% Create a Rayleigh channel System object
rchan1 = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',1);%���ӻ�
% 'Visualization','Impulse and frequency responses')


after_Ray = rchan1(datain_ALL);
datain_ALL = awgn(after_Ray,snr);%%%%%%%%aafterfm2:ͨ��Rayleigh�ŵ���aafterfm:��ͨ��Rayleigh�ŵ�

% datain_ALL = awgn(datain_ALL,snr);%ֻͨ��awgn

%% ��Ƶƫ
fd=0;%HZƵƫ
Rb=10e5;%%%%%%%%%��Դ��������
Ts=1./Rb;

cont=1:length(datain_ALL);%%%%%%%%%%ÿһ��Ķ���Ƶƫ
phase_pian = 2j*pi*fd.*Ts.*cont;%%��Ƶƫ���飺phase_pian = 2j*pi*fd.*Ts.*cont*0
datain_ALL=datain_ALL.*exp(phase_pian');

% scatterplot(datain_ALL)%show


%% ʱ��ͬ��
[Data_atertimelock,judg,Guard_atertimelock]=TimeLockFun(PN,datain_ALL,length(modData));
% if (judg==1)
%     !echo TIME LOCK success    
% else
%     !echo TIME LOCK Failed����
% end

%%%%%%%�����пɿ���TimeLockFun�еĻ�ͼ���֡�

%% Ƶƫ����
if ( Guard_atertimelock==404)
    BER=404;
    f_averge=404;%%%%%%%%%%%��������
else
[f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

%% ���

% ofdm���
demod = comm.OFDMDemodulator(mod);  
[dataOut, pilotOut] = step(demod,Data_atertimelock);%%%%%%%%%%%ʹ��ʱ�滻ΪData_atertimelock!!!!!!Data_atertFrelock
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


data_FINI = deinterleaver(data_FINI);%%%%%%%%%�⽻֯
data_FINI= step(hDec, data_FINI);%%%%%%%%ά�ر�����

af_crc=step(crcDet,data_FINI);
% errorRate = comm.ErrorRate;
% errors = step(hError, data, receivedBits);
        errVec = hError(dain,af_crc(1:length(dain)));   % data_FINI,dain;data_FINI(1:length(dain))
%         ber=[ber;errVec(1)];
BER=errVec(1);
end