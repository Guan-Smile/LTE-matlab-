function [f_averge,dataout]=frequencLock(Data_atertimelock,Datain_G,Npn,K,Ts)%%%%�������ݣ�datain,Ƶ��ͬ��ʹ��PN���и���
cont=1:length(Data_atertimelock);

% %% ��������
% fd=100;%HZƵƫ
% K=255;%=K
% Rb=10e5
% Ts=1./Rb;
% Npn=20;%%���Ƶ��ͬ��pn����
% guard_length=20;%�������г���
% 
% cont=1:K.*Npn+guard_length*2;
% 
% %% Guard��������
% guard=zeros(1,guard_length);
% 
% 
% %% PN��������
% h = commsrc.pn('GenPoly',[8 6 5 4 0],'NumBitsOut',255);%���ȣ�255
% Hpn=generate(h);
% 
% datain=ones(Npn,1)*Hpn';
% 
% Datain=reshape(datain',[],1);
% 
% Datain_G= [guard';Datain;guard'];%%%���뱣������
% 
% pskModulator = comm.PSKModulator('ModulationOrder',2,'PhaseOffset',0);
% channelInput = step(pskModulator,Datain_G);
% % scatterplot(channelInput)%show
% 
% %% ����Ƶƫ
% phase_pian = 2j*pi*fd.*Ts.*cont;%%��Ƶƫ���飺phase_pian = 2j*pi*fd.*Ts.*cont*0
% aafterfm=channelInput.*exp(phase_pian');
% 
% scatterplot(aafterfm)%show
% 
% %% ���ŵ�
% 
% awchan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',-5,'RandomStream','mt19937ar with seed', ...
%     'Seed',98);
% aw_out=awchan(aafterfm);
% scatterplot(aw_out)


%% ����PN���н��л����

%  pure_pn=aw_out(guard_length+1:(end-guard_length));%��ȥ����������ȡpn����
pure_PN=reshape(Datain_G,[],Npn);%����ΪNpn��pn����

pure_PN_c=conj(pure_PN);%����

for i=1:Npn-1
    P(:,i)=pure_PN_c(:,i).*pure_PN(:,i+1);
end
PP=sum(P);
T=K.*Ts;
fmax=1./(2*T);
f=(2*pi.*T).^(-1).*(atan(imag(PP)./real(PP)));

f_averge=mean(f);%%%���Ƶ�Ƶƫ

%% ����ƵƫӰ��

phase_pian_xiao = 2j*pi*f_averge.*Ts.*cont;%%��Ƶƫ���飺phase_pian = 2j*pi*fd.*Ts.*cont*0
dataout=Data_atertimelock.*exp(phase_pian_xiao');
% scatterplot(dataout);
% 
% 
% %% ��֤��BPSK�����������
% pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);
% data_FINI = step(pskDemodulator,dataout);
% 
% errorRate = comm.ErrorRate;
% 
% %   rxData = pskDemodulator(rxSig);
%         % Collect the error statistics
%         errVec = errorRate(data_FINI,Datain_G);
% 


