clear
lambda=2;
t=1;
Nsim=1e6;
%% Countdown Simulation
vN=zeros(Nsim,1);
for i=1:Nsim
    T=[]; % jump times
    T(1)=0;
    while T(end)<t
        % sample interarrival time
        tau=exprnd(1/lambda);
        T=[T, T(end)+tau];
    end
    T=T(2:end-1);
    Nt=length(T);
    vN(i)=Nt;
end
mean(vN)
var(vN)
%% Conditional Simulation
Nt=poissrnd(lambda*t,Nsim,1);
mean(Nt)
var(Nt)


