clear; close all;
%% Simulation of the Merton model
mu=0.05; %drift
sigma=0.4; %volatility
lambda=2; %Poisson intensity
muJ=0.01; deltaJ=0.2; %Jumpsize parameters

T=2 % Maturity
M=100; % Number of steps in time
dt=T/M;
X=zeros(M+1,1);
Z=randn(M,1);
NT=poissrnd(lambda*T);
jumpT=sort(rand(1,NT)*T);
jumpSize=muJ+deltaJ*randn(NT,1);
for i=1:M
    X(i+1)=X(i)+mu*dt+sigma*sqrt(dt)*Z(i);
    % check if there are jumps in ((i-1)dt,idt]
    for j=1:NT
        if jumpT(j)>(i-1)*dt && jumpT(j)<=i*dt
            X(i+1)=X(i+1)+jumpSize(j);
        end
    end
end
plot(exp(X))

    
    
    
    