function [Data_use,judg,Guard_atertimelock]=TimeLockFun(PN,datain,lengh_data)
%%%%[数据帧包，时间同步是否成功(1:成功，0失败)]=TimeLockFun（PN序列，输入数据，数据帧长度）
%% 时间同步
%% 参数设置
K=510;%每个pn序列的长度
Npn=4;%%一共N段pn序列
% guard_length=20;%保护序列长度
% snr=-15
% global lengh_data%%%%%%%数据帧长度
% lengh_data=230
y=[];%门限
Tc=2;
men=[];

Tongbu_Num=Npn;%%同步个数

% %% 序列生成
%  [PN,datain]=Canshu(Npn);
% %  scatterplot(datain)%n段的pn序列+保护间隔
%  
%  
%  %% 过信道
%  awchan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',snr);
% aw_out=awchan(datain);



%% 时间同步：
tim=1:size(datain)-K+1;%%%
for kk=tim;
 use_relat=datain(kk:K-1+kk);
 dedsss_D=conj(use_relat).*PN;%乘pn序列后。双极性
 jude_sum=sum(dedsss_D);
 cd=sum(conj(use_relat).*use_relat);
 men=[men,Tc.*cd]; %门限
%  if (jude_sum>=b(kk))
%      record(kk)=jude_sum 
%  end
 y=[y,jude_sum];
end


dely_t=tim;
MK=y.*conj(y);
% figure(3)
% plot(dely_t,MK,'b-*')
% hold on 
% plot(dely_t,b,'k-*')
% grid on
% xlabel('采样时间k'); ylabel('M[k]'); 
% % legend('Rayleigh','AWGN+手动多径')%

%% 得出时间同步的值

for j=1:size(men,2)
%     MK(j);
%     men(j);
   if (MK(j)<men(j))
       MK(j)=0;    
   end
end
[dat,loc]=sort(MK);%%%最大值与所在位置(位置非正序)
%% 精同步：
true_loc=[];
ori=[];
for i=1:50

maxLoc1=sort(loc(end-Tongbu_Num:end));
% maxLoc2=[0,sort(loc(end-Tongbu_Num:end))]:

% dmaxLoc=maxLoc1-maxLoc2(1:size(maxLoc1,2));

b=bsxfun(@minus,maxLoc1,maxLoc1.');%求两两个之间的差值的结果
[x,y]=find(b==K);

true_loc_cont=unique([x,y]);%%%%正确的位置的序号
true_loc=maxLoc1(true_loc_cont);%%正确的位置横坐标
if ~isequal(ori,[])
true_loc=[ori,true_loc(Tongbu_Num:end)];
end

if (length(true_loc)>=4)
%     true_loc%%%%%%%%%%%%%%%%%show
%     !echo TIME LOCK test?
    JUDG2=diff(true_loc);
    [val,weizhi]=find(diff(true_loc)~=K);
    if isequal(weizhi,[1,2])
        true_loc(2)=[];                
    elseif isequal(weizhi,[1])
         true_loc(1)=[];        
    elseif isequal(weizhi,[2,3])
         true_loc(3)=[];       
    elseif isequal(weizhi,[3])
         true_loc(4)=[];        
    elseif isequal(weizhi,[1,2,3])
        true_loc(2)=[]; 
        true_loc(3)=[];
    end    
    if (length(true_loc)>=4)
%         true_loc;%%%%%%%%%%%%show
       break 
    end
end

Tongbu_Num=Tongbu_Num+1;
ori=true_loc;

end
    

%% 绘图：
% % figure()
% % plot(dely_t,MK,'b-*')
% % hold on 
% % plot(dely_t,men,'k-*')
% % grid on
% % xlabel('采样时间k'); ylabel('M[k]'); 
% % % legend('Rayleigh','AWGN+手动多径')%
% % plot(loc(end-Tongbu_Num:end),dat(end-Tongbu_Num:end),'rP')
if (length(true_loc)<3)
    Data_use=404;
else
    Data_use=datain(true_loc(2)+(Npn-fix(true_loc(2)/K))*K:true_loc(2)+(Npn-fix(true_loc(2)/K))*K+lengh_data-1);%%%%%%%%%%%%同步剪除后得到的数据帧部分
end
    % length(Data_use)show
if (length(true_loc)<3)
    Data_use=404;
else    
Guard_atertimelock=datain(true_loc(2)-fix(true_loc(2)/K)*K:true_loc(2)+(Npn-fix(true_loc(2)/K))*K-1);
end
% Guard_atertimelo=length(Guard_atertimelock)% show
if (length(Data_use)==lengh_data)

   pskDemodulator = comm.PSKDemodulator('ModulationOrder',2,'PhaseOffset',0);
   data_FINI = pskDemodulator(Data_use);
%     !echo TIME LOCK success
    judg=1;
    
else
    !echo TIME LOCK Failed！！
    judg=0;
    Data_use=[Data_use;zeros(lengh_data-length(Data_use),1)];
    Guard_atertimelock=404;
end

%% 绘图：
% figure(3)
% plot(dely_t,MK,'b-*')
% hold on 
% plot(dely_t,men,'k-*')
% grid on
% xlabel('采样时间k'); ylabel('M[k]'); 
% % legend('Rayleigh','AWGN+手动多径')%
% plot(loc(end-Tongbu_Num:end),dat(end-Tongbu_Num:end),'rP')