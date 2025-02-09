function J=Levy_Integral_full(x,V,ymin,ymax,k,S0,Kdisc)
% integral_ymin^ymax ...
%   (V(x+y)-V(x)-(exp(y)-1)V'(x))k(y)dy
N=length(x);
Nq=2*N; y=linspace(ymin,ymax,Nq);
dy=y(2)-y(1); w=ones(size(y))*dy; 
w(1)=w(1)/2; w(end)=w(end)/2;
w=w.*k(y);
dx=x(2)-x(1);
dV=(V(3:end)-V(1:end-2))./(2*dx);
J=zeros(N,1);
for i=2:N-1
    J(i)=sum( w.*...
        ( funV(x(i)+y,x,V,S0,Kdisc)...
          -V(i)-dV(i-1)*(exp(y)-1))  );
end

function Vf=funV(z,x,V,S0,Kdisc)
Vf=zeros(size(z));
% inside the grid
index=find( (z>x(1)).*(z<x(end)) );
Vf(index)=interp1(x,V,z(index));
% outside the grid --> BC
index=find( (z>=x(end)) );
% V= S0*exp(z)-K*exp(-r(T-t))
Vf(index)=S0*exp(z(index))-Kdisc;