%% Algorithm of Dassios-Zhao

rng(1929) %11
lambda=1; alpha=0.4; beta=0.6; %Setting of parameters % beta > alpha
lambda_0=3;  % lambda_0 > lambda
T=5;  % Time horizon for simulation 

Tj=0; %Vector of time jumps
lambda_tk=lambda_0;
Jump_time=0;

while (Jump_time<T)
    U_1=rand(1);
    U_2=rand(1);
    D_k=1+(beta*log(U_1))/(lambda_tk(end)-lambda);

    %Simulation of interarrival time
    if (D_k>0)
    S_1=(-1/beta)*log(D_k);
    S_2=(-1/lambda)*log(U_2);

    S_k=min(S_1,S_2);

    else 
     S_k=(-1/lambda)*log(U_2);   
    end

Jump_time=Tj(end)+S_k;
if(Jump_time<=T)
Tj=[Tj ; Jump_time];  %append new jump_time
intensity_value=lambda_0+sum(alpha*exp(-beta*Tj));

lambda_tk=[lambda_tk ; intensity_value]; %Update intensity
end
end


%creation of a vector which contains values for intensity in order to plot
%decaying
x=linspace(0,T,10000) ;
intensity_function=linspace(0,T,10000);
intensity_function(1)=lambda_0;
Tj=[Tj; T]; % useful for handling cycles
for i=2:length(x)
for j=2:length(Tj)
    if (Tj(j-1)<x(i) && Tj(j)>=x(i))
  intensity_function(i)= lambda_0+alpha*sum(exp(beta*Tj(2:j-1)))*exp(-beta*x(i));
    end
end
end
values=@(y) sum(heaviside(y-Tj));
%figure(1)
subplot(2,1,1)
plot(x,values(x),'-k','LineWidth',2);
grid on;
title('Hawkes process');
xlabel('t');
ylabel('N_t');

% figure(2)
subplot(2,1,2)
plot(x,intensity_function,'-k','LineWidth',2);
grid on;
title('Hawkes intensity');
xlabel('t');
ylabel('\lambda_t');