function [x,P]=ekf(func,x,P,func2,z,Q,R)
% Inputs:   f: function handle for f(x)
%           x: previous state estimate
%           P: previous estimated state covariance
%           h: function handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance

% Output:   x:  state estimate
%           P:  state covariance
[x1,A]=jacob(func,x);       %linearization of first function
P=A*P*A'+Q;                 
[z1,H]=jacob(func2,x1);     %linearization of second function
P12=P*H';                   %cross covariance
K=P12/(H*P12+R);            %Kalman filter gain
x=x1+K*(z-z1);              %state estimate
P=P-K*P12';                 %state covariance matrix


function [z,A]=jacob(fun,x)
% Jacobian of function fun
z=fun(x);
n=numel(x);
m=numel(z);
A=zeros(m,n);
h=n*eps;
for k=1:n
    x1=x;
    x1(k)=x1(k)+h*1i;
    A(:,k)=imag(fun(x1))/h;
end
