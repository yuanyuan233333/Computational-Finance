clear; close all;
% Price a European Call Option
% Merton process -> logprice PDE
% Finite Difference - Implicit Euler (Operator Splitting)
%--> We do not split the integral!!!

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
%% Matrix
matA=sparse(N+1,N+1);
A=-(r-sigma^2/2)/(2*dx)+...
    sigma^2/(2*dx^2);
B=-1/dt-sigma^2/(dx^2)-(r);
C=(r-sigma^2/2)/(2*dx)+...
    sigma^2/(2*dx^2);
matA(1,1)=1; matA(end,end)=1;
for i=2:N
     matA(i, [i-1 i i+1])=[A B C];
end
%% starting value --> payoff
V=max( S0*exp(x)-K, 0);
%% backward-in-time procedure
rhs=zeros(N+1,1); J=rhs;
for j=M-1:-1:0
    if lambda>0
        J=Levy_Integral_full(x,V,ymin,ymax,k,...
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



