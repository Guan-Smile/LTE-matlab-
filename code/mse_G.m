function [MSE]=mse_G(f_gu,f_real)



MSE1=f_gu./f_real-1;
MSE2=abs(MSE1);
MSE=(mean(MSE2));
% ���ŵ��������ͼ��
% figure()
%  plot(real(f_real));
%  hold on;
%  plot(real(f_gu),'-p');
%   title('�Ա���ʵ�ŵ����ڲ�������ŵ�')
%   legend('��ʵ�ŵ�1','��ʵ�ŵ�2','�ڲ����ŵ�1','�ڲ����ŵ�2')
