function [f_averge,dataout]=frequencLock(Data_atertimelock,Datain_G,Npn,K,Ts)%%%%传入数据：datain,频率同步使用PN序列个数
cont=1:length(Data_atertimelock);

% %% 参数设置
% fd=100;%HZ频偏
% K=255;%=K
% Rb=10e5
% Ts=1./Rb;
% Npn=20;%%多段频率同步pn序列
% guard_length=20;%保护序列长度
% 
% cont=1:K.*Npn+guard_length*2;
% 
% %% Guard保护序列
% guard=zeros(1,guard_length);
% 
% 
% %% PN序列生成
% h = commsrc.pn('GenPoly',[8 6 5 4 0],'NumBitsOut',255);%长度：255
% Hpn=generate(h);
% 
% datain=ones(Npn,1)*Hpn';
% 
% Datain=reshape(datain',[],1);
% 
% Datain_G= [guard';Datain;guard'];%%%加入保护序列
% 
% pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
% channelInput = step(pskModulator,Datain_G);
% % scatterplot(channelInput)%show
% 
% %% 加入频偏
% phase_pian = 2j*pi*fd.*Ts.*cont;%%无频偏检验：phase_pian = 2j*pi*fd.*Ts.*cont*0
% aafterfm=channelInput.*exp(phase_pian');
% 
% scatterplot(aafterfm)%show
% 
% %% 过信道
% 
% awchan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',-5,'RandomStream','mt19937ar with seed', ...
%     'Seed',98);
% aw_out=awchan(aafterfm);
% scatterplot(aw_out)


%% 相邻PN序列进行互相关

%  pure_pn=aw_out(guard_length+1:(end-guard_length));%除去保护序列提取pn序列
pure_PN=reshape(Datain_G,[],Npn);%重组为Npn个pn序列

pure_PN_c=conj(pure_PN);%共轭

for i=1:Npn-1
    P(:,i)=pure_PN_c(:,i).*pure_PN(:,i+1);
end
PP=sum(P);
T=K.*Ts;
fmax=1./(2*T);
f=(2*pi.*T).^(-1).*(atan(imag(PP)./real(PP)));

f_averge=mean(f);%%%估计的频偏

%% 消除频偏影响

phase_pian_xiao = 2j*pi*f_averge.*Ts.*cont;%%无频偏检验：phase_pian = 2j*pi*fd.*Ts.*cont*0
dataout=Data_atertimelock.*exp(phase_pian_xiao');
% scatterplot(dataout);
% 
% 
% %% 验证：BPSK解调，误码率
% pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);
% data_FINI = step(pskDemodulator,dataout);
% 
% errorRate = comm.ErrorRate;
% 
% %   rxData = pskDemodulator(rxSig);
%         % Collect the error statistics
%         errVec = errorRate(data_FINI,Datain_G);
% 


