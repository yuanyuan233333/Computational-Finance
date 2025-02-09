clear; close all;
%% Simulation of the Kou model
mu=0.05; %drift
sigma=0.4; %volatility
lambda=2; %Poisson intensity
p=0.6; lambda_plus=10; lambda_minus=3; %Jumpsize parameters

T=2 % Maturity
M=100; % Number of steps in time
dt=T/M;
X=zeros(M+1,1);
Z=randn(M,1);
NT=poissrnd(lambda*T);
jumpT=sort(rand(1,NT)*T);
for i=1:M
    X(i+1)=X(i)+mu*dt+sigma*sqrt(dt)*Z(i);
    % check if there are jumps in ((i-1)dt,idt]
    for j=1:NT
        if jumpT(j)>(i-1)*dt && jumpT(j)<=i*dt
            u=rand;
            if p<u %positive jump
               jumpSize=exprnd(1/lambda_plus);
            else %negative jump
               jumpSize=-exprnd(1/lambda_minus);
            end               
            X(i+1)=X(i+1)+jumpSize;
        end
    end
end
plot(exp(X))

    
    
    
    