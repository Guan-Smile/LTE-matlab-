clear
clc
close all
%% �����������Demo.m��������mainFun.m��ȫ��ͬ���������ԭ������ʹ�á�
%% 1.����֡��

%��������
stp=14;%%%%��Ƶ�����14=>15
N_subcarrier=1024;%���ز���128,1024
Npn=4;%%һ��N��pn����
global K
snr=5;

%% 1.1�����ofdm����֡��
% ����ofdm���ƽ��ģ��
kk=1:stp+1:N_subcarrier;%(stp-1)/2:stp+1:1024;%��ǰ���15
kkk=[kk,kk+1];

mod = comm.OFDMModulator('NumGuardBandCarriers',[0;N_subcarrier-kk(end)-1],...
'PilotInputPort',true, ...%�Ƿ���뵼Ƶ
'FFTLength',N_subcarrier,...
'PilotCarrierIndices',[sort(kkk)'], ...%���뵼Ƶ��λ������
'NumSymbols',1, ...
'CyclicPrefixLength',106,...%ѭ��ǰ׺�ĳ���
'InsertDCNull',false);  
modDim = info(mod); 
% % showResourceMapping(mod)  %show
% ��������datain����Ƶ֡
% 
% mod2 = comm.OFDMModulator('NumGuardBandCarriers',[1;N_subcarrier-kk(end)-1],...
% 'PilotInputPort',true, ...%�Ƿ���뵼Ƶ
% 'FFTLength',N_subcarrier,...
% 'PilotCarrierIndices',[kk'+1], ...%���뵼Ƶ��λ������
% 'NumSymbols',1, ...
% 'CyclicPrefixLength',106,...%ѭ��ǰ׺�ĳ���
% 'InsertDCNull',false);  
% modDim2 = info(mod2); 
% % %  showResourceMapping(mod2)  %show
% % ��������datain����Ƶ֡

Fo_pilot=1:stp:modDim.DataInputSize(1)+length(kkk)/2-1;%%��ƵԤ����λ������
Fo_pilot2=stp:stp:modDim.DataInputSize(1)+length(kkk)/2;

%%% ofdm���ز������������߷��������ز�����(����������Ƶ���õ�λ��)
% rng(22);
dain=randi([0 1],modDim.DataInputSize(1),1);

% OF_Data1=reshape(ANT1',1,[]);%%%%%%%%%%%% ��������֡

% rng(29);
dain2=randi([0 1],modDim.DataInputSize(1),1);
% ANT2=[nodain,dain2];
% OF_Data2=reshape(ANT2',1,[]);%%%%%%%%%%% ��������֡

pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
dataIn = step(pskModulator,dain);% ����֡bpsk
dataIn2 = step(pskModulator,dain2);% ����֡bpsk

% dataIn2=zeros(size(dataIn));

% dataIn(Fo_pilot)=0;%%%%%%%%%%%% ��������֡
% dataIn2(Fo_pilot2)=0;%%%%%%%%%%%% ��������֡

DATAIN=[dain,dain2];
DATAIN_AFBPSK=[dataIn,dataIn2];
% % OF_Data1=reshape(ANT1',1,[]);%%%%��
% % OF_Data2=reshape(ANT2',1,[]);%%%��
% scatterplot(dataIn)%show
% Create a pilot signal that has the correct dimensions. 
pilotIn_part = complex(ones(modDim.PilotInputSize(1)/2,1),ones(modDim.PilotInputSize(1)/2,1)); % ���ɵ�Ƶ֡
nopilot=zeros(modDim.PilotInputSize(1)/2,1);
pilotIn=[pilotIn_part,nopilot];
pilotIn=reshape(pilotIn',1,[])';
pilotIn2=[nopilot,pilotIn_part];
pilotIn2=reshape(pilotIn2',1,[])';
% Apply OFDM modulation to the data and pilot signals. 
modData = step(mod,dataIn,pilotIn).*sqrt(N_subcarrier);%%%%%%%%%%ofdm������ɵ�����  ��128+106ѭ��ǰ׺��
modData2 = step(mod,dataIn2,pilotIn2).*sqrt(N_subcarrier);%%%%%%%%%%ofdm������ɵ�����  ��128+106ѭ��ǰ׺��
%% 1.2�����PN���е�֡��
[PN,datain_ALL]=Canshu(Npn,modData);%%%%datain_ALL��ɵ�֡��
[PN,datain_ALL2]=Canshu(Npn,modData2);%%%%datain_ALL��ɵ�֡��
%  datain_ALL=zeros(size(datain_ALL2));
%  dataIn=zeros(size(dataIn2));
%  dain=zeros(size(dain));
channelInput=[datain_ALL,datain_ALL2];
% channelInput=[datain_ALL,zeros(size(datain_ALL2))];%%%��Ч������
DATAIN=[dain,dain2];
DATAIN_AFBPSK=[dataIn,dataIn2];
%% ����


%%  ��Rayleigh
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
%     'MaximumDopplerShift',fD);%���ӻ�
% rchan_out_11= step(rchan1,datain_ALL);
% rchan_out_21= step(rchan1,datain_ALL2);
% rchan_out_12= step(rchan1,datain_ALL);
% rchan_out_22= step(rchan1,datain_ALL2);
% 
% LTEChanOut=[rchan_out_21,0.1*rchan_out_12+rchan_out_22];
%% ���ŵ�

% % Rayleigh�ŵ�eva
% fs = 4e6;                                     % Hz
% pathDelays = [0 3e-8 15e-8 31e-8 37e-8 71e-8 109e-8 173e-8 251e-8];    % sec
% avgPathGains = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];      % dB
% fD = 1;                                         % Hz
% % Rayleigh�ŵ�2
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
% LTEChanOut = awgn(channelInput,snr);%ֻͨ��awgn
AW_OUT = awgn(channelInput,snr);%ֻͨ��awgn
AW_OUT2 = awgn(channelInput,snr);%ֻͨ��awgn
AW_OUT3 = awgn(channelInput,snr);%ֻͨ��awgn
AW_OUT4 = awgn(channelInput,snr);%ֻͨ��awgn
% LTEChanOut = [0.2*AW_OUT(:,1)+ 0.4*AW_OUT(:,2),1*AW_OUT(:,1)+0.7*AW_OUT(:,2)];
% LTEChanOut = [0.1*AW_OUT(:,1)+1*AW_OUT2(:,2),1*AW_OUT3(:,1)+ 0.1*AW_OUT4(:,2)];

%% ��Ƶƫ
fd=0;%HZƵƫ
Rb=10e5;%%%%%%%%%��Դ��������
Ts=1./Rb;

for lh=1:size(LTEChanOut,2)
cont=1:length(LTEChanOut(:,lh));%%%%%%%%%%ÿһ��Ķ���Ƶƫ
phase_pian = 2j*pi*fd.*Ts.*cont;%%��Ƶƫ���飺phase_pian = 2j*pi*fd.*Ts.*cont*0
LTEChanOut(:,lh)=LTEChanOut(:,lh).*exp(phase_pian');

% scatterplot(LTEChanOut(:,lh))%show
end

%% ʱ��ͬ��
Y=[];
H=[];

for lh=1:size(LTEChanOut,2)%%%%%%%%%%%��ѭ��1.��������1��  2.��������2��
    
[Data_atertimelock,judg,Guard_atertimelock]=TimeLockFun(PN,LTEChanOut(:,lh),length(modData));
 % Ƶƫ����
[f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);
%�������´洢


Data_AT_ALL(:,lh)=Data_atertimelock;%%%%%%%%%%%%%%%Data_atertimelock///Data_atertFrelock
Guard_AT_ALL(:,lh)=Guard_atertimelock;
% if (judg==1)
%     !echo TIME LOCK success    
% else
%     !echo TIME LOCK Failed����
% end

%%%%%%%�����пɿ���TimeLockFun�еĻ�ͼ���֡�

%% Ƶƫ����
% MIMO��ǰ�ˣ�
% [f_averge,Data_atertFrelock]=frequencLock(Data_atertimelock,Guard_atertimelock,Npn,K,Ts);

%% ���

% ofdm���h11,h12
demod = comm.OFDMDemodulator(mod);  
[dataOut_h11, pilotOut_h11] = step(demod,Data_atertFrelock);%%%%%%%%%%%ʹ��ʱ�滻ΪData_atertimelock!!!!!!Data_atertFrelock

dataOut_h11 = dataOut_h11 ./ sqrt(N_subcarrier);
H_real2(:,1)=dataOut_h11./DATAIN_AFBPSK(:,1);%%%%%%%%%h11
H_real2(:,2)=dataOut_h11./DATAIN_AFBPSK(:,2);%%%%%%%%%h21
pilotOut_h11 = pilotOut_h11 ./ sqrt(N_subcarrier);

%% �ŵ�����
H_gu1= pilotOut_h11./(pilotIn);%%%%%%h12
H_gu2= pilotOut_h11./(pilotIn2);%%%%%%h21

pilotOut_T1=reshape(H_gu1,2,[]);
pilotOut_T11=pilotOut_T1(1,:).';%%%%%%%%%%H11
% pilotOut_T21=pilotOut_T1(2,:)';%%%%%%%%%%inf

pilotOut_T1_2=reshape(H_gu2,2,[]);
pilotOut_T21=pilotOut_T1_2(2,:).';%%%%%%%%%%H21
% pilotOut_T21_2=pilotOut_T1_2(2,:)';%%%%%%%%%inf

%�Ƚ�
% isSame = (max(abs([dataIn(:) - dataOut(:); ...
%     pilotIn(:) - pilotOut(:)])) < 1e-10)

% 
% H_gu_h11=(pilotOut_T11)./(pilotIn_part); % lx: ���
% H_gu_h21=(pilotOut_T21)./(pilotIn_part); % lx: ���

 %����H����
 
 Hin=pilotOut_T11(:,1);%%%%%%%�����ڲ���H11////%%%%%%%�����ڲ���H12
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

 Hin=pilotOut_T21(:,1);%%%%%%%�����ڲ���H21////%%%%%%%�����ڲ���H22
 
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
H=[H,temp_x21];%%%%%%%%%%h����˳��h11,h21,h12,h22
Y=[Y,dataOut_h11];%%%%%%%%����˳��h11,h21,h12,h22

temp_x21=[];

%% ceshi
figure()
 plot(real(H_real2));
 hold on;
 plot(real(H(:,2*(lh-1)+1:2*(lh))),'-p');
  title('�Ա���ʵ�ŵ����ڲ�������ŵ�')
  legend('��ʵ�ŵ�1','��ʵ�ŵ�2','�ڲ����ŵ�1','�ڲ����ŵ�2')
MSE(lh,:)=mse_G(H(:,2*(lh-1)+1:2*(lh)),H_real2)

end
% figure()
% plot(1:length(Y),Y,'-P');
%% �����Ȼ��
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
%% ���ܷ���

%% ����������
pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);

for iwi=1:size(Get_DATA,2)
data_FINI = step(pskDemodulator,Get_DATA(:,iwi));%%%%%%%%%%%���������о���BPSK�����dataOut(����ǰ)ND_OUT�����ƺ�

errorRate = comm.ErrorRate;
        errVec = errorRate(data_FINI,DATAIN(:,iwi));   % data_FINI,dain;
%         ber=[ber;errVec(1)];
BER=errVec


figure()
stem(data_FINI/10,'P')
hold on
stem(DATAIN(:,iwi)/10,'o')
end

ND_OUT=  Y(:,1)./ H(:,1);%Ƶ��ԭ�źŹ���ֵ���Ƴ����ŵ�
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,1));   % data_FINI,dain;
BER_zf11=errVec(1)
ND_OUT=  Y(:,1)./ H(:,2);%Ƶ��ԭ�źŹ���ֵ���Ƴ����ŵ�
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,2));   % data_FINI,dain;
BER_zf21=errVec(1)

ND_OUT=  Y(:,2)./ H(:,3);%Ƶ��ԭ�źŹ���ֵ���Ƴ����ŵ�
data_FINI = step(pskDemodulator,ND_OUT);
  reset(errorRate)
 errVec = errorRate(data_FINI,DATAIN(:,1));   % data_FINI,dain;
BER_zf12=errVec(1)

ND_OUT=  Y(:,2)./ H(:,4);%Ƶ��ԭ�źŹ���ֵ���Ƴ����ŵ�
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
%   title('�Ա���ʵ�ŵ����ڲ�������ŵ�')
%   legend('��ʵ�ŵ�','�ڲ����ŵ�')