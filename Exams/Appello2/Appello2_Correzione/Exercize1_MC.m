%%  MC ex1

%% parameters
T = 2; K = 100; L = 80; U = 110;
S0 = 95; sigmaBM = 0.4; r = 0.0367;
% ext VG parameters
sigmaVG = 0.12; thetaVG = 0.03; kVG = 0.20;
psi = @(u) -1/kVG * log(1 + u.^2*sigmaVG^2*kVG/2 - 1i*thetaVG*kVG.*u);
drift = r - sigmaBM^2/2 - psi(-1i);

%
Nsim = 1e6;
Ndates = round(12*T); %monthly monitoring
dt = T/Ndates;


%% STANDARD MC
%
% Asset simulation
% Initializations
X = zeros(Nsim,Ndates+1);

tic
for i = 1:Ndates
    %sim
    u = rand(Nsim,1);
    Z = randn(Nsim,1);
    W = randn(Nsim,1);
    %
    deltaS = kVG * icdf('Gamma',u,dt/kVG,1);
    %
    X(:,i+1) = X(:,i) + drift*dt + sigmaBM*sqrt(dt)*W + thetaVG*deltaS + sigmaVG*Z.*sqrt(deltaS);
end

S = S0*exp(X);
%% dmcheck
disp('risk-neutral check must be close to 0')
[check,~,CI]=normfit(S(:,end)-S0*exp(r*T))
%%
toc
% disc payoff
disc_payoff = exp(-r*T) .* max(S(:,end)-K,0) .* (min(S,[],2)>L) .* (max(S,[],2)<U);

%
[price,dummy,IC] = normfit(disc_payoff)




%% Antitetich variable


% Initializations
X = zeros(Nsim,Ndates+1);
X_AV = zeros(Nsim,Ndates+1);
tic
for i = 1:Ndates
    %sim
    u = rand(Nsim,1);
    Z = randn(Nsim,1);
    W = randn(Nsim,1);
    %
    deltaS = kVG * icdf('Gamma',u,dt/kVG,1);
    deltaS_AV = kVG * icdf('Gamma',1-u,dt/kVG,1);
    %
    X(:,i+1) = X(:,i) + drift*dt + sigmaBM*sqrt(dt)*W + thetaVG*deltaS + sigmaVG*Z.*sqrt(deltaS);
    X_AV(:,i+1) = X_AV(:,i) + drift*dt - sigmaBM*sqrt(dt)*W + thetaVG*deltaS_AV + sigmaVG*Z.*sqrt(deltaS_AV);
end

S = S0*exp(X);
S_AV = S0*exp(X);
toc

% disc payoff
disc_payoff = exp(-r*T) .* max(S(:,end)-K,0) .* (min(S,[],2)>L) .* (max(S,[],2)<U);
disc_payoff_AV = exp(-r*T) .* max(S_AV(:,end)-K,0) .* (min(S_AV,[],2)>L) .* (max(S_AV,[],2)<U);

% price
[price_AV,dummy,IC_AV] = normfit((disc_payoff+disc_payoff_AV)/2)






%% Control Variable

% Asset simulation with Nsim/100
% Initializations
X = zeros(Nsim/100,Ndates+1);
tic
for i = 1:Ndates
    %sim
    u = rand(Nsim/100,1);
    Z = randn(Nsim/100,1);
    W = randn(Nsim/100,1);
    %
    deltaS = kVG * icdf('Gamma',u,dt/kVG,1);
    %
    X(:,i+1) = X(:,i) + drift*dt + sigmaBM*sqrt(dt)*W + thetaVG*deltaS + sigmaVG*Z.*sqrt(deltaS);
end

S = S0*exp(X);

%
f = S(:,end);
g = exp(-r*T) .* max(f-K,0) .* (min(S,[],2)>L) .* (max(S,[],2)<U);
Ef = S0*exp(r*T);

VC = cov(f,g);
alpha = -VC(1,2)/VC(1,1);

% Asset simulation with Nsim
% Initializations
X = zeros(Nsim,Ndates+1);
tic
for i = 1:Ndates
    %sim
    u = rand(Nsim,1);
    Z = randn(Nsim,1);
    W = randn(Nsim,1);
    %
    deltaS = kVG * icdf('Gamma',u,dt/kVG,1);
    %
    X(:,i+1) = X(:,i) + drift*dt + sigmaBM*sqrt(dt)*W + thetaVG*deltaS + sigmaVG*Z.*sqrt(deltaS);
end

S = S0*exp(X);

%
f = S(:,end);
g = exp(-r*T) .* max(f-K,0) .* (min(S,[],2)>L) .* (max(S,[],2)<U);


[price_CV,dummy,IC_CV] = normfit(g+alpha*(f-Ef))
toc
