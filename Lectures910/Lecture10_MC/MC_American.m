clear all
% Longstaff & Schwartz algorithm --------
%=== L=3 ================================
% Quadratic Polinomial

%==== INPUT =============================
S0=1; T=3; r=0.06; K=1.1;
M=3; dt=T/M;
N=8; % number of simulation
%--- From L&S article
S=[1.09 1.08 1.34;
   1.16 1.26 1.54;
   1.22 1.07 1.03;
   0.93 0.97 0.92;
   1.11 1.56 1.52;
   0.76 0.77 0.90;
   0.92 0.84 1.01;
   0.88 1.22 1.34];
%--- MC simulation 
% r=0.04; T=1; M=40; N=1e6; sigma=0.4; S0=1; K=1;
% S=AssetBS(r,sigma,S0,T,M,N);
% S=S(:,2:end); % we do not consider t=0 for early exercise
% dt=T/M;
%==== INITIALIZE =========================
Exercise_Time=M*ones(N,1);
Put=max(0,K-S(:,end)); %payoff

%==== BACKWARD-IN-TIME ===================
for j=M-1:-1:1
    Inmoney=find(S(:,j)<K);
    S_I=S(Inmoney,j);
    %-- Intrinsic Value
    IV=K-S_I;
    %-- Continuation Value 
    %- Regression
    A=[ones(length(S_I),1), S_I, S_I.^2];
    b=Put(Inmoney).*...
exp(-r*dt*(Exercise_Time(Inmoney)-j));
    alpha=A\b;
    %- Continuation Value 
    CV=A*alpha;
    %----------
    %== j is an exercise instant?
    Index=find(IV>CV);
    Early_Exercise=Inmoney(Index);
    % Update
    Put(Early_Exercise)=IV(Index);
    Exercise_Time(Early_Exercise)=j;
end
[price,~,CI]=normfit(Put.*exp(-r*dt*Exercise_Time))
%[~,EuPrice]=blsprice(S0,K,r,T,sigma)
