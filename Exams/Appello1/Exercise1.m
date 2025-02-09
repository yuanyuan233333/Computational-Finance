clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT DELETE THE FOLLOWING STRING
% afmslafsdsdskdjdqwke3222312w3122
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input
S0=1; K=1; r=0.001; T=1;
% Levy triplet parameters
sigma=0.6;
lambda=3; muJ=-0.02; deltaJ=0.4;
k=@(y)lambda*exp(-(y-muJ).^2/(2*deltaJ^2))...
   /sqrt( 2*pi*deltaJ^2 ); 
% Discretization parameters
M=50; N=1000; % Unconditionally stable 
%% grids
dt=T/M; 
xmin=log( 0.1);%Smin=0.1S0=S0exp(xmin)
xmax=log( 3);%Smax=3S0=S0exp(xmax)
dx=(xmax-xmin)/N; 
x=linspace(xmin,xmax,N+1)';
[alpha,lambda_num,ymin,ymax]=...
    Levy_trunc(k,N)
%% Matrix
matA=sparse(N+1,N+1);
A=-(r-sigma^2/2-alpha)/(2*dx)+...
    sigma^2/(2*dx^2);
B=-1/dt-sigma^2/(dx^2)-(r+lambda);
C=(r-sigma^2/2-alpha)/(2*dx)+...
    sigma^2/(2*dx^2);
matA(1,1)=1; matA(end,end)=1;
for i=2:N
     matA(i, [i-1 i i+1])=[A B C];
end
%matA(2:end-1,2:end-1)=...
% spdiags([A*ones(N-1,1) B*ones(N-1,1), C*ones(N-1,1)],...
%            [-1 0 1],N-1,N-1);
%matA(2,1)=A; matA(end-1,end)=C;        
%% starting value --> payoff
V=max( S0*exp(x)-K, 0);
%% backward-in-time procedure
rhs=zeros(N+1,1); J=rhs;
for j=M-1:-1:0
    if lambda>0
        J=Levy_Integral(x,V,ymin,ymax,k,...
            S0,K*exp(-r*(T-(j+1)*dt)));
    end
    rhs(2:end-1)=-V(2:end-1)/dt;
    rhs(end)=S0*exp(xmax)-K*exp(-r*(T-j*dt));
    V=matA\(rhs-J);
end
S=S0*exp(x);
figure
plot(S,V); title('Call Price');
Price_FD=interp1( S,V, S0,'spline')        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT DELETE THE FOLLOWING STRING
% askdj2312312
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,lambda_num,ymin,ymax]=...
    Levy_trunc(k,N)
% Truncation of the integral domain
tol=1e-10;
ymin=-0.5;
while k(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.5;
while k(ymax)>tol
    ymax=ymax+0.5;
end
% Computation of alpha and lambda
Nq=2*N; y=linspace(ymin,ymax,Nq);
figure
plot(y,k(y))
lambda_num=trapz(y,k(y));
% dy=y(2)-y(1); w=ones(size(y))*dy; 
% w(1)=w(1)/2; w(end)=w(end)/2;
% lambda=sum(w.*k(y));
alpha=trapz(y,(exp(y)-1).*k(y));
end

function J=Levy_Integral(x,V,ymin,ymax,k,S0,Kdisc)
% integral_ymin^ymax V(x+y)k(y)dy
N=length(x);
Nq=2*N; y=linspace(ymin,ymax,Nq);
dy=y(2)-y(1); w=ones(size(y))*dy; 
w(1)=w(1)/2; w(end)=w(end)/2;
w=w.*k(y);
J=zeros(N,1);
for i=2:N-1
    J(i)=sum( w.*funV(x(i)+y,x,V,S0,Kdisc) );
end
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
end
