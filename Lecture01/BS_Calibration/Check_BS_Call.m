clear; close all;
%>European Call Options
Maturity=4/12;
Data=[25.30	200 %Prices, Strikes
10.00	225
18.20	210
41.85	180
15.20	215
7.85	230
12.36	220
46.05	175
70.70	150
50.85	170
60.40	160
29.40	195];
%Data=[Strike,Time to Maturity,Market Price]
Data=[Data(:,2), Maturity*ones(size(Data(:,2))), Data(:,1)]
r=0.02; spot=218.75;
sigma=0.2;
options=optimoptions('fsolve','FunctionTolerance',1e-10,'StepTolerance',1e-8,'OptimalityTolerance',1e-8,'MaxIterations',800,'MaxFunctionEvaluations',500)
for i=1:size(Data,1)
    sigma=fsolve(@(vol) ...
        fun(vol,spot,Data(i,1),r,Data(i,2),Data(i,3)),...
        sigma,options);
    Vol(i)=sigma;
end
figure
plot(Data(:,1),Vol,'d'); 
xlabel('Strike'); ylabel('Impliedvol');










