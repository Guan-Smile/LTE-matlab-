function [MSE]=mse_G(f_gu,f_real)



MSE1=f_gu./f_real-1;
MSE2=abs(MSE1);
MSE=(mean(MSE2));
% 画信道估计情况图用
% figure()
%  plot(real(f_real));
%  hold on;
%  plot(real(f_gu),'-p');
%   title('对比真实信道和内插出来的信道')
%   legend('真实信道1','真实信道2','内插后的信道1','内插后的信道2')
