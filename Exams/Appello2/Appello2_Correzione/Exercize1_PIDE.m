clear; close all;


%% parameters
T = 2; K = 100; L = 80; U = 110;
S0 = 95; sigma = 0.4; r = 0.0367;
% ext VG parameters
sigmaVG = 0.12; thetaVG = 0.03; kVG = 0.20;
A = thetaVG/sigmaVG^2;
B = sqrt(thetaVG^2 + 2*sigmaVG^2/kVG)/sigmaVG^2;
k=@(y) 1/kVG./abs(y) .* exp(A.*y - B.*abs(y));

% Discretization parameters
M=50; N=2000; % Unconditionally stable 
%% grids
dt=T/M; 
xmin=log( L/S0);%%>>> L and U in the boundaries
xmax=log( U/S0);
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
V(1) = 0; V(end) = 0;%%>>> L and U in the payoff


%% backward-in-time procedure
rhs=zeros(N+1,1); J=rhs;
for j=M-1:-1:0
    J=Levy_Integral_full_KO(x,V,ymin,ymax,k,S0,K*exp(-r*(T-(j+1)*dt)));
    rhs(2:end-1)=-V(2:end-1)/dt;
    V=matA\(rhs-J);
    %%>>> zero BC
end
S=S0*exp(x);
figure
plot(S,V); title('Call Price');
Price_FD=interp1( S,V, S0,'spline')        


