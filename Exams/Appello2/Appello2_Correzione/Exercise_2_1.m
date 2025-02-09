clear all
clc
format long
T=1;
S0=95;
r=0.0367;
sigmaGBM=0.126349;
muJ=-0.390078;
deltaJ=0.338796;
lambda=0.174814;
Nsim=1e6;
Ndates=52;
disp('MC')
tic
[Path,~]=simulation_Merton_av(S0,T,Ndates,Nsim,sigmaGBM,muJ,deltaJ,lambda,r);

%% dmcheck
disp('risk-neutral check must be close to 0')
[check,~,CI]=normfit(Path(:,end)-S0*exp(r*T))
%%


disc_payoff=exp(-r*T)*max(Path(:,end)-mean(Path(:,[2:end]),2),0);
toc
[Price_MC,~,IC_MC]=normfit(disc_payoff)

disp('AV')
tic
[Path,Path_av]=simulation_Merton_av(S0,T,Ndates,Nsim,sigmaGBM,muJ,deltaJ,lambda,r);

disc_payoff=exp(-r*T)*max(Path(:,end)-mean(Path(:,[2:end]),2),0);
disc_payoff_av=exp(-r*T)*max(Path_av(:,end)-mean(Path_av(:,[2:end]),2),0);
toc
[Price_av,~,IC_av]=normfit((disc_payoff + disc_payoff_av )/2)
