function [f_averge,BER]=mainFun(Snr,fD)%%%%
%% 1.构造帧包

%参数设置
stp=14;%%%%导频间隔：15
N_subcarrier=1024;%子载波数128,1024
Npn=4;%%一共N段pn序列
global K

%% 1.1构造ofdm数据帧包
% 设置ofdm调制解调模块
kk=1:stp+1:N_subcarrier;%(stp-1)/2:stp+1:1024;%当前间隔15
mod = comm.OFDMModulator('NumGuardBandCarriers',[0;N_subcarrier-kk(end)],...
'PilotInputPort',true, ...%是否加入导频
'FFTLength',N_subcarrier,...
'PilotCarrierIndices',[kk'], ...%加入导频的位置序列
'NumSymbols',1, ...
'CyclicPrefixLength',106,...%循环前缀的长度
'InsertDCNull',false);  
modDim = info(mod); 
% showResourceMapping(mod)  %show
% 生成数据datain，导频帧
% rng(22);
dain=randi([0 1],modDim.DataInputSize(1),1);% 生成数据帧
pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
dataIn = step(pskModulator,dain);% 数据帧bpsk
% scatterplot(dataIn)%show
% Create a pilot signal that has the correct dimensions. 
pilotIn = complex(ones(modDim.PilotInputSize),ones(modDim.PilotInputSize)); % 生成导频帧
% Apply OFDM modulation to the data and pilot signals. 
modData = step(mod,dataIn,pilotIn).*sqrt(N_subcarrier);%%%%%%%%%%ofdm调制完成的数据  （128+106循环前缀）

%% 1.2构造加PN序列的帧包
[PN,datain_ALL]=Canshu(Npn,modData);%%%%datain_ALL完成的帧包

%% 调制\卷积\交织



%% 过信道

% Rayleigh信道1
fs = 4e6;                                     % Hz
pathDelays = [0 3e-8 15e-8 31e-8 37e-8 71e-8 109e-8 173e-8 251e-8];    % sec
avgPathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];      % dB
% fD = 1;                                         % Hz
% % Rayleigh信道2
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-6];    % sec
% avgPathGains = [0 -10];      % dB
% fD = 1; 

% Create a Rayleigh channel System object
rchan1 = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',fD);%可视化
% 'Visualization','Impulse and frequency responses')


after_Ray = rchan1(datain_ALL);
datain_ALL = awgn(after_Ray,Snr);%%%%%%%%aafterfm2:通过Rayleigh信道，aafterfm:不通过Rayleigh信道

% datain_ALL = awgn(datain_ALL,Snr);%只通过awgn

%% 加频偏  fd
fd=100;%HZ频偏
Rb=10e5;%%%%%%%%%信源比特速率
Ts=1./Rb;

cont=1:length(datain_ALL);%%%%%%%%%%每一项的都会频偏
phase_pian = 2j*pi*fd.*Ts.*cont;%%无频偏检验：phase_pian = 2j*pi*fd.*Ts.*cont*0
datain_ALL=datain_ALL.*exp(phase_pian');

% scatterplot(datain_ALL);%show


%% 时间同步
[Data_atertimelock,judg,Guard_atertimelock]=TimeLockFun(PN,datain_ALL,length(modData));
% if (judg==1)
%     !echo TIME LOCK success    
% else
%     !echo TIME LOCK Failed！！
% end


%% 频偏估计
if ( Guard_atertimelock==404)
    BER=404
    f_averge=404;%%%%%%%%%%%！！！！
else
[f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

% [f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

%% 解调

% ofdm解调
demod = comm.OFDMDemodulator(mod);  
[dataOut, pilotOut] = step(demod,Data_atertFrelock);%%%%%%%%%%%使用时替换为Data_atertimelock!!!!!!Data_atertFrelock
dataOut = dataOut ./ sqrt(N_subcarrier);
pilotOut = pilotOut ./ sqrt(N_subcarrier);
%比较
% isSame = (max(abs([dataIn(:) - dataOut(:); ...
%     pilotIn(:) - pilotOut(:)])) < 1e-10)

%% 信道估计

H_gu=(pilotOut)./(pilotIn); % lx: 点除

 %处理H矩阵
 
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


%% 信道均衡

ND_OUT= dataOut./ (temp_lx);%频域原信号过插值估计出的信道

%% 性能分析

%% 计算误码率
pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);
data_FINI = step(pskDemodulator,ND_OUT);%%%%%%%%%%%最终数据判决（BPSK解调）dataOut(估计前)ND_OUT（估计后）

errorRate = comm.ErrorRate;
        errVec = errorRate(data_FINI,dain);   % data_FINI,dain;
%         ber=[ber;errVec(1)];
BER=errVec(1);
end