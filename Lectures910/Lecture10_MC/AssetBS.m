function S=AssetBS(r, sigma, S0, T, M, Nsim)

X=zeros(Nsim,M+1);
dt=T/M;
for i=1:M
    X(:,i+1)=X(:,i)+(r-sigma^2/2)*dt+...
        sigma*sqrt(dt)*randn(Nsim,1);
end
S=S0*exp(X);