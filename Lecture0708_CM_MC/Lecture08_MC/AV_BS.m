clear
Nsim=1e6;
r=0.01; sigma=0.6; K=100; S0=96; T=1;
x1=randn(Nsim,1);
x2=-x1;
ST1=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*x1 );
ST2=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*x2 );
discpayoff1=exp(-r*T)*max( ST1-K,0);
discpayoff2=exp(-r*T)*max( ST2-K,0);
corr(x1,x2)
corr(ST1,ST2)
corr(discpayoff1,discpayoff2)