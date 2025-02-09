clear; close all;
Stock %Upload Data
Apple_MarketData=flipud(Apple_MarketData); %first element is the older
logreturn=log(Apple_MarketData(2:end)./...
    Apple_MarketData(1:end-1));
%logreturn=logreturn(end/2:end)
figure; hold on;
plot(logreturn);
title('Logreturn');
figure
hist(logreturn,40)
title('Logreturn histogram');

% Estimate B&S parameters from a historical time series
mean_value=mean(logreturn);
variance=var(logreturn);
dt=1/252;
sigma=sqrt( variance/dt )
mu=mean_value/dt+sigma^2/2
% Check Normality assumption
figure
qqplot(logreturn)
sk=skewness(logreturn)
ku=kurtosis(logreturn)





