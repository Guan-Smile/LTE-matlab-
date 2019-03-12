clear
clc
close all
%% 本代码独立于Demo.m，功能与mainFun.m完全等同，用于理解原理，调试使用。
%% 1.构造帧包

%参数设置
stp=14;%%%%导频间隔：14=>15
N_subcarrier=1024;%子载波数128,1024
Npn=4;%%一共N段pn序列
global K
snr=5;

%% 1.1构造非ofdm数据帧包
% 设置ofdm调制解调模块
kk=1:stp+1:N_subcarrier;%(stp-1)/2:stp+1:1024;%当前间隔15
kkk=[kk,kk+1];

mod = comm.OFDMModulator('NumGuardBandCarriers',[0;N_subcarrier-kk(end)-1],...
'PilotInputPort',true, ...%是否加入导频
'FFTLength',N_subcarrier,...
'PilotCarrierIndices',[sort(kkk)'], ...%加入导频的位置序列
'NumSymbols',1, ...
'CyclicPrefixLength',106,...%循环前缀的长度
'InsertDCNull',false);  
modDim = info(mod); 
% % showResourceMapping(mod)  %show
% 生成数据datain，导频帧
% 
% mod2 = comm.OFDMModulator('NumGuardBandCarriers',[1;N_subcarrier-kk(end)-1],...
% 'PilotInputPort',true, ...%是否加入导频
% 'FFTLength',N_subcarrier,...
% 'PilotCarrierIndices',[kk'+1], ...%加入导频的位置序列
% 'NumSymbols',1, ...
% 'CyclicPrefixLength',106,...%循环前缀的长度
% 'InsertDCNull',false);  
% modDim2 = info(mod2); 
% % %  showResourceMapping(mod2)  %show
% % 生成数据datain，导频帧

Fo_pilot=1:stp:modDim.DataInputSize(1)+length(kkk)/2-1;%%导频预留空位置序列
Fo_pilot2=stp:stp:modDim.DataInputSize(1)+length(kkk)/2;

%%% ofdm子载波正交，两天线发出的子载波正交(数据留出导频放置的位置)
% rng(22);
dain=randi([0 1],modDim.DataInputSize(1),1);

% OF_Data1=reshape(ANT1',1,[]);%%%%%%%%%%%% 生成数据帧

% rng(29);
dain2=randi([0 1],modDim.DataInputSize(1),1);
% ANT2=[nodain,dain2];
% OF_Data2=reshape(ANT2',1,[]);%%%%%%%%%%% 生成数据帧

pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
dataIn = step(pskModulator,dain);% 数据帧bpsk
dataIn2 = step(pskModulator,dain2);% 数据帧bpsk

% dataIn2=zeros(size(dataIn));

% dataIn(Fo_pilot)=0;%%%%%%%%%%%% 生成数据帧
% dataIn2(Fo_pilot2)=0;%%%%%%%%%%%% 生成数据帧

DATAIN=[dain,dain2];
DATAIN_AFBPSK=[dataIn,dataIn2];
% % OF_Data1=reshape(ANT1',1,[]);%%%%舍
% % OF_Data2=reshape(ANT2',1,[]);%%%舍
% scatterplot(dataIn)%show
% Create a pilot signal that has the correct dimensions. 
pilotIn_part = complex(ones(modDim.PilotInputSize(1)/2,1),ones(modDim.PilotInputSize(1)/2,1)); % 生成导频帧
nopilot=zeros(modDim.PilotInputSize(1)/2,1);
pilotIn=[pilotIn_part,nopilot];
pilotIn=reshape(pilotIn',1,[])';
pilotIn2=[nopilot,pilotIn_part];
pilotIn2=reshape(pilotIn2',1,[])';
% Apply OFDM modulation to the data and pilot signals. 
modData = step(mod,dataIn,pilotIn).*sqrt(N_subcarrier);%%%%%%%%%%ofdm调制完成的数据  （128+106循环前缀）
modData2 = step(mod,dataIn2,pilotIn2).*sqrt(N_subcarrier);%%%%%%%%%%ofdm调制完成的数据  （128+106循环前缀）
%% 1.2构造加PN序列的帧包
[PN,datain_ALL]=Canshu(Npn,modData);%%%%datain_ALL完成的帧包
[PN,datain_ALL2]=Canshu(Npn,modData2);%%%%datain_ALL完成的帧包
%  datain_ALL=zeros(size(datain_ALL2));
%  dataIn=zeros(size(dataIn2));
%  dain=zeros(size(dain));
channelInput=[datain_ALL,datain_ALL2];
% channelInput=[datain_ALL,zeros(size(datain_ALL2))];%%%等效单天线
DATAIN=[dain,dain2];
DATAIN_AFBPSK=[dataIn,dataIn2];
%% 调制


%%  单Rayleigh
fs = 4e6;                                     % Hz
pathDelays = [0 3e-6];    % sec
avgPathGains = [0 -10];      % dB
fD = 100;                                         % Hz
%  Rayleigh channel System object with the previously defined
% parameters and set the |Visualization| property to |Impulse and frequency
% responses| using name-value pairs.
% rchan1 = comm.RayleighChannel('SampleRate',fs, ...
%     'PathDelays',pathDelays, ...
%     'AveragePathGains',avgPathGains, ...
%     'MaximumDopplerShift',fD);%可视化
% rchan_out_11= step(rchan1,datain_ALL);
% rchan_out_21= step(rchan1,datain_ALL2);
% rchan_out_12= step(rchan1,datain_ALL);
% rchan_out_22= step(rchan1,datain_ALL2);
% 
% LTEChanOut=[rchan_out_21,0.1*rchan_out_12+rchan_out_22];
%% 过信道

% % Rayleigh信道eva
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-8 15e-8 31e-8 37e-8 71e-8 109e-8 173e-8 251e-8];    % sec
% avgPathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];      % dB
% fD = 1;                                         % Hz
% % Rayleigh信道2
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-6];    % sec
% avgPathGains = [0 -10];      % dB
% fD = 1; 

% Create an |LTEMIMOChannel| System object with a 2-by-2 antenna
% configuration and a medium correlation level.
lteChan = comm.LTEMIMOChannel(...
    'Profile',              'EVA 5Hz',...
    'AntennaConfiguration', '2x2',...
    'CorrelationLevel',     'Medium',...
    'AntennaSelection',     'Off',...
    'RandomStream',         'mt19937ar with seed',...
    'Seed',                 33,...
    'PathGainsOutputPort',  true);

mimoChan = comm.MIMOChannel('SampleRate',4e6, 'PathDelays',[0 1e-6], ...
    'AveragePathGains',[0 -4], 'NormalizePathGains',false, 'MaximumDopplerShift',1, ...
    'TransmitCorrelationMatrix',cat(3,eye(2),[1 0.1;0.1 1]), ...
    'ReceiveCorrelationMatrix',cat(3,[1 0.2;0.2 1],eye(2)), ...
    'RandomStream','mt19937ar with seed', 'Seed',33, 'PathGainsOutputPort',true);

% Filter the modulated data using the equivalent |mimoChannel| object.
[LTEChanOut,LTEPathGains] = lteChan(channelInput);%%%%%%%%%%%lteChan
% [ChanOut,pathGains] = mimoChan(channelInput);%%%%%%%%%%normal mimoChan
% LTEChanOut = ChanOut;
% LTEChanOut = awgn(channelInput,snr);%只通过awgn
AW_OUT = awgn(channelInput,snr);%只通过awgn
AW_OUT2 = awgn(channelInput,snr);%只通过awgn
AW_OUT3 = awgn(channelInput,snr);%只通过awgn
AW_OUT4 = awgn(channelInput,snr);%只通过awgn
% LTEChanOut = [0.2*AW_OUT(:,1)+ 0.4*AW_OUT(:,2),1*AW_OUT(:,1)+0.7*AW_OUT(:,2)];
% LTEChanOut = [0.1*AW_OUT(:,1)+1*AW_OUT2(:,2),1*AW_OUT3(:,1)+ 0.1*AW_OUT4(:,2)];

%% 加频偏
fd=0;%HZ频偏
Rb=10e5;%%%%%%%%%信源比特速率
Ts=1./Rb;

for lh=1:size(LTEChanOut,2)
cont=1:length(LTEChanOut(:,lh));%%%%%%%%%%每一项的都会频偏
phase_pian = 2j*pi*fd.*Ts.*cont;%%无频偏检验：phase_pian = 2j*pi*fd.*Ts.*cont*0
LTEChanOut(:,lh)=LTEChanOut(:,lh).*exp(phase_pian');

% scatterplot(LTEChanOut(:,lh))%show
end

%% 时间同步
Y=[];
H=[];

for lh=1:size(LTEChanOut,2)%%%%%%%%%%%大循环1.接受天线1；  2.接收天线2；
    
[Data_atertimelock,judg,Guard_atertimelock]=TimeLockFun(PN,LTEChanOut(:,lh),length(modData));
 % 频偏估计
[f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);
%数据重新存储


Data_AT_ALL(:,lh)=Data_atertimelock;%%%%%%%%%%%%%%%Data_atertimelock///Data_atertFrelock
Guard_AT_ALL(:,lh)=Guard_atertimelock;
% if (judg==1)
%     !echo TIME LOCK success    
% else
%     !echo TIME LOCK Failed！！
% end

%%%%%%%调试中可开启TimeLockFun中的绘图部分。

%% 频偏估计
% MIMO提前了：
% [f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

%% 解调

% ofdm解调h11,h12
demod = comm.OFDMDemodulator(mod);  
[dataOut_h11, pilotOut_h11] = step(demod,Data_atertFrelock);%%%%%%%%%%%使用时替换为Data_atertimelock!!!!!!Data_atertFrelock

dataOut_h11 = dataOut_h11 ./ sqrt(N_subcarrier);
H_real2(:,1)=dataOut_h11./DATAIN_AFBPSK(:,1);%%%%%%%%%h11
H_real2(:,2)=dataOut_h11./DATAIN_AFBPSK(:,2);%%%%%%%%%h21
pilotOut_h11 = pilotOut_h11 ./ sqrt(N_subcarrier);

%% 信道估计
H_gu1= pilotOut_h11./(pilotIn);%%%%%%h12
H_gu2= pilotOut_h11./(pilotIn2);%%%%%%h21

pilotOut_T1=reshape(H_gu1,2,[]);
pilotOut_T11=pilotOut_T1(1,:).';%%%%%%%%%%H11
% pilotOut_T21=pilotOut_T1(2,:)';%%%%%%%%%%inf

pilotOut_T1_2=reshape(H_gu2,2,[]);
pilotOut_T21=pilotOut_T1_2(2,:).';%%%%%%%%%%H21
% pilotOut_T21_2=pilotOut_T1_2(2,:)';%%%%%%%%%inf

%比较
% isSame = (max(abs([dataIn(:) - dataOut(:); ...
%     pilotIn(:) - pilotOut(:)])) < 1e-10)

% 
% H_gu_h11=(pilotOut_T11)./(pilotIn_part); % lx: 点除
% H_gu_h21=(pilotOut_T21)./(pilotIn_part); % lx: 点除

 %处理H矩阵
 
 Hin=pilotOut_T11(:,1);%%%%%%%线性内插求H11////%%%%%%%线性内插求H12
 for i=1:size(Hin)-1
    %for j=1:stp+1
    for k=1:stp
         ND(k+(i-1).*(stp+1))=Hin(i+1).*(k-1)./(stp+1)+Hin(i).*(stp+2-k)./(stp+1);
         co1 = (k)./(stp+1);
         co2 = (stp+1-k)./(stp+1);
         temp_x11(k+(i-1).*(stp))=Hin(i+1).*(k)./(stp+1)+Hin(i).*(stp+1-k)./(stp+1);

%           ND(j+(i-1)*(stp+1))=Hin(i);
    end    
 end
%  ND(stp+2+(size(Hin)-2).*(stp+1))=Hin(size(Hin));
temp_x11 = temp_x11.';
temp_x11(Fo_pilot)=[];

H=[H,temp_x11];
temp_x11=[];

 Hin=pilotOut_T21(:,1);%%%%%%%线性内插求H21////%%%%%%%线性内插求H22
 
 for i=1:size(Hin)-1
    %for j=1:stp+1
    for k=1:stp
         ND(k+(i-1).*(stp+1))=Hin(i+1).*(k-1)./(stp+1)+Hin(i).*(stp+2-k)./(stp+1);
         co1 = (k)./(stp+1);
         co2 = (stp+1-k)./(stp+1);
         temp_x21(k+(i-1).*(stp))=Hin(i+1).*(k)./(stp+1)+Hin(i).*(stp+1-k)./(stp+1);

%           ND(j+(i-1)*(stp+1))=Hin(i);
    end    
 end
%  ND(stp+2+(size(Hin)-2).*(stp+1))=Hin(size(Hin));
temp_x21 = temp_x21.';
temp_x21(Fo_pilot2)=[];
H=[H,temp_x21];%%%%%%%%%%h矩阵顺序：h11,h21,h12,h22
Y=[Y,dataOut_h11];%%%%%%%%矩阵顺序：h11,h21,h12,h22

temp_x21=[];

%% ceshi
figure()
 plot(real(H_real2));
 hold on;
 plot(real(H(:,2*(lh-1)+1:2*(lh))),'-p');
  title('对比真实信道和内插出来的信道')
  legend('真实信道1','真实信道2','内插后的信道1','内插后的信道2')
MSE(lh,:)=mse_G(H(:,2*(lh-1)+1:2*(lh)),H_real2)

end
% figure()
% plot(1:length(Y),Y,'-P');
%% 最大似然法
Get_DATA=[];
for j=1:length(H)
% Hg=reshape(H(j,:)',2,2);
Hg=[H(j,1),H(j,2);H(j,3),H(j,4)];
Y_USE=[Y(j,1);Y(j,2)];
HD1=Hg*[1;1];
HD2=Hg*[1;-1];
HD3=Hg*[-1;1];
HD4=Hg*[-1;-1];

SS=[[1,1];[1,-1];[-1,1];[-1,-1]];

D=[sum(abs(Y_USE-HD1).^2,1),sum(abs(Y_USE-HD2).^2,1),sum(abs(Y_USE-HD3).^2,1),sum(abs(Y_USE-HD4).^2,1)];



[VAL,PLACE]=min(D');

Get_DATA=[Get_DATA;SS(PLACE,:)];

% Get_DATA=reshape(Get_DATA,2,[]);
end
%% 性能分析

%% 计算误码率
pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);

for iwi=1:size(Get_DATA,2)
data_FINI = step(pskDemodulator,Get_DATA(:,iwi));%%%%%%%%%%%最终数据判决（BPSK解调）dataOut(估计前)ND_OUT（估计后）

errorRate = comm.ErrorRate;
        errVec = errorRate(data_FINI,DATAIN(:,iwi));   % data_FINI,dain;
%         ber=[ber;errVec(1)];
BER=errVec


figure()
stem(data_FINI/10,'P')
hold on
stem(DATAIN(:,iwi)/10,'o')
end

ND_OUT=  Y(:,1)./ H(:,1);%频域原信号过插值估计出的信道
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,1));   % data_FINI,dain;
BER_zf11=errVec(1)
ND_OUT=  Y(:,1)./ H(:,2);%频域原信号过插值估计出的信道
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,2));   % data_FINI,dain;
BER_zf21=errVec(1)

ND_OUT=  Y(:,2)./ H(:,3);%频域原信号过插值估计出的信道
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,1));   % data_FINI,dain;
BER_zf12=errVec(1)

ND_OUT=  Y(:,2)./ H(:,4);%频域原信号过插值估计出的信道
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,2));   % data_FINI,dain;
BER_zf22=errVec
%% zf BER
% for i=1:2
%     for j=2%1:4
%         ND_OUT=  Y(:,i)./ H(:,j);
%         data_FINI = step(pskDemodulator,ND_OUT);
%         reset(errorRate);
%         errVec = errorRate(data_FINI,DATAIN(:,1)); 
%         BER_zfd1=errVec(1)
%         reset(errorRate);
%         errVec = errorRate(data_FINI,DATAIN(:,2)); 
%         BER_zfd2=errVec(1)
%     end
% end
% % 
% figure()
% stem(data_FINI/10,'P')
% hold on
% stem(DATAIN(:,2)/10,'o')
oo=0;

% %% ceshi
% figure
%  plot(real(H_real2));
%  hold on;
%  plot(real(H(:,1)).*sqrt(N_subcarrier),'r-p');
%   title('对比真实信道和内插出来的信道')
%   legend('真实信道','内插后的信道')