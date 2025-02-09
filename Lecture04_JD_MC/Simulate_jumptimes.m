clear
lambda=2;
t=1;
%% Countdown Simulation
T=[]; % jump times
T(1)=0;
while T(end)<t
    % sample interarrival time
    tau=exprnd(1/lambda);
    T=[T, T(end)+tau];
end
T=T(2:end-1)
Nt=length(T)
%% Conditional Simulation
Nt=poissrnd(lambda*t)
T=sort(rand(1,Nt)*t)


