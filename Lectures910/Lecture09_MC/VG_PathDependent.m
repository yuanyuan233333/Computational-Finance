clear
%% Pricing Floating Stike Asian Call under VG
% Market/Contract Input
S0=100;
r=2/100;
K=105; T=0.5;
% monitoring -> daily
M=round(252*T);
% Model (VG) Input
theta=0.2; sigma=0.6; k=0.5;
% MC Input
Nsim=1e5;

%% Setup
dt=T/M;
char_exp=@(u) -log(1+u.^2*sigma^2*k/2-1i*theta*k*u)/k;
drift=r-char_exp(-1i);

%% 1. Standard MC
X=zeros(Nsim,M+1);
dS=k*icdf('gamma',rand(Nsim,M),dt/k,1);
Z=randn(Nsim,M);
for j=1:M
    X(:,j+1)=X(:,j)+drift*dt+...
        theta*dS(:,j)+sigma*sqrt(dS(:,j)).*Z(:,j);
end
S=S0*exp(X); clear X dS Z;
%disp('r-n check: must be zero')
%[check,~,CI_check]=normfit(S(:,end)-S0*exp(r*T))
discpayoff=exp(-r*T)*max( S(:,end)-mean(S(:,2:end),2),0);
[price,~,CI]=normfit(discpayoff)
%% AV MC
X=zeros(Nsim,M+1); XAV=zeros(Nsim,M+1);
dS=k*icdf('gamma',rand(Nsim,M),dt/k,1);
Z=randn(Nsim,M);
for j=1:M
    X(:,j+1)=X(:,j)+drift*dt+...
        theta*dS(:,j)+sigma*sqrt(dS(:,j)).*Z(:,j);
    XAV(:,j+1)=XAV(:,j)+drift*dt+...
        theta*dS(:,j)-sigma*sqrt(dS(:,j)).*Z(:,j);
end
S=S0*exp(X); SAV=S0*exp(XAV); clear X dS Z XAV;
%disp('r-n check: must be zero')
%[check,~,CI_check]=normfit(SAV(:,end)-S0*exp(r*T))
[price,~,CI]=normfit( exp(-r*T)*(...
    max( S(:,end)-mean(S(:,2:end),2),0) +...
    max( SAV(:,end)-mean(SAV(:,2:end),2),0) )/2 )
%% 3. Control Variable MC
Ef=S0*exp(r*T)-mean( S0*exp(r*dt*(1:M)) )
% Estimate Alpha
Nsim2=1000;
X=zeros(Nsim2,M+1);
dS=k*icdf('gamma',rand(Nsim2,M),dt/k,1);
Z=randn(Nsim2,M);
for j=1:M
    X(:,j+1)=X(:,j)+drift*dt+...
        theta*dS(:,j)+sigma*sqrt(dS(:,j)).*Z(:,j);
end
S=S0*exp(X); 
f=S(:,end)-mean( S(:,2:end),2 );
g=max(f,0)*exp(-r*T);
VC=cov(f,g);
alpha=-VC(1,2)/VC(1,1);
% Compute the price
X=zeros(Nsim,M+1);
dS=k*icdf('gamma',rand(Nsim,M),dt/k,1);
Z=randn(Nsim,M);
for j=1:M
    X(:,j+1)=X(:,j)+drift*dt+...
        theta*dS(:,j)+sigma*sqrt(dS(:,j)).*Z(:,j);
end
S=S0*exp(X); clear X dS Z;
f=S(:,end)-mean( S(:,2:end),2 ); clear S;
g=max(f,0)*exp(-r*T);
[price,~,CI]=normfit( g+alpha*(f-Ef) )


