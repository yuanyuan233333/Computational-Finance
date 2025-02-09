clear; close all;
% Price a Up and Out Put Option
% Merton process -> logprice PDE
% Finite Difference - Implicit Euler (Operator Splitting)
%--> We do not split the integral!!!

%% input
S0=95; K=95; r=0.0367; T=1; U = 125;
% Levy triplet parameters
sigma=0.126349; mu = -0.390078;
lambda=0.174814; muJ=-0.390078; deltaJ=0.338796;
k=@(y)lambda*exp(-(y-muJ).^2/(2*deltaJ^2))...
   /sqrt( 2*pi*deltaJ^2 ); 
% Discretization parameters
M=50; N=1000; % Unconditionally stable 
%% grids
dt=T/M; 
xmin=log( 0.1);
xmax=log( U/S0);%U in the domain
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
V=max( K-S0*exp(x), 0) ;
V(end) = 0; %>> U in the BC
%% backward-in-time procedure
rhs=zeros(N+1,1); J=rhs;
for j=M-1:-1:0
    if lambda>0
        J=Levy_Integral_full_P(x,V,ymin,ymax,k,...
            S0,K*exp(-r*(T-(j+1)*dt)));
    end
    rhs(2:end-1)=-V(2:end-1)/dt;
    rhs(1)=K*exp(-r*(T-j*dt)) - S0*exp(xmin);
    V=matA\(rhs-J);
end
S=S0*exp(x);
figure
plot(S,V); title('Put Price');
Price_FD=interp1( S,V, S0,'spline')        


