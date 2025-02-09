function [Path,Path_av]=simulation_Merton_av(S0,T,Ndates,Nsim,sigmaGBM,muJ,deltaJ,lambda,r)
dt=T/Ndates;
X=zeros(Nsim,Ndates+1);
X_av=zeros(Nsim,Ndates+1);
t=linspace(0,T,Ndates+1);
psi=@(u) -sigmaGBM^2/2*u.^2 ...
   +lambda*(exp(-deltaJ^2*u^2/2+1i*muJ*u)-1);
mu=r-psi(-1i);
Nt=icdf('Poisson',rand(Nsim,1),lambda*T);
for i=1:Nsim
    Tj=sort(T*rand(Nt(i),1));

    for j=1:Ndates
        Z_1=randn;
        X(i,j+1)=X(i,j) + mu*dt + sigmaGBM*sqrt(dt)*Z_1;
        X_av(i,j+1)=X(i,j) + mu*dt - sigmaGBM*sqrt(dt)*Z_1;
        for k=1:Nt(i)
            if (t(j)<Tj(k) && Tj(k)< t(j+1) )
                Z=randn;
                Y=muJ + deltaJ*Z;
                X(i,j+1)=X(i,j+1)+Y;
                Y_av=muJ - deltaJ*Z;
                X_av(i,j+1)=X_av(i,j+1)+Y_av;
            end
        end
    end

end



Path=S0*exp(X);
Path_av=S0*exp(X_av);









