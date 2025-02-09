clear
mu=0.03; sigma=0.5; S0=1;
T=10;
disp('Compute E[S(T)]');
Exact=S0*exp(mu*T)
disp('Exact Solution Simulation')
Nsim=1e6;
[Value,~,CI]=normfit( ...
    S0*exp((mu-sigma^2/2)*T+sigma*sqrt(T)*randn(Nsim,1)))
disp('Euler scheme')
M=1 
dt=T/M;
S=S0;
for i=1:M
    S=S.*(1+mu*dt+sigma*sqrt(dt)*randn(Nsim,1));
end
[Value,~,CI]=normfit( S )
M=10
dt=T/M;
S=S0;
for i=1:M
    S=S.*(1+mu*dt+sigma*sqrt(dt)*randn(Nsim,1));
end
[Value,~,CI]=normfit( S )
M=100
dt=T/M;
S=S0;
for i=1:M
    S=S.*(1+mu*dt+sigma*sqrt(dt)*randn(Nsim,1));
end
[Value,~,CI]=normfit( S )