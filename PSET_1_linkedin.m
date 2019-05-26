clear;
clc;
load dataset.mat;
%% 1
figure
plot(dates,[OIL,SP500./10]); 
legend('OIL' , 'SP500');
% the sp500 is divided by ten to adjust the scale with the price of oil

%% 2
% Compute daily log changes in oil prices and in the stock market index (henceforth "returns")
%and report mean, standard deviation, skewness, kurtosis and 1rst order auto-correlation.

oil=price2ret(OIL);
sp=price2ret(SP500);

oil_mean=mean(oil);
oil_stdev=sqrt(var(oil));
oil_skewness=skewness(oil);
oil_kur=kurtosis(oil);

sp_mean=mean(sp);
sp_stdev=sqrt(var(sp));
sp_skewness=skewness(sp);
sp_kur=kurtosis(sp);

%% 3
%Compute the VaR and the ES with 1% and 5% confidence level for oil and stock market
%returns. Comment on the methodology you have used to compute the value-at-risk and the
%expected shortfall. 
n=size(oil,1);
n1=28;
n5=140;

oil_var1=max(oil(1:n1,1));
oil_var5=max(oil(1:n5,1));
oil_es1=mean(oil(1:n1,1));
oil_es5=mean(oil(1:n5,1));

sp_var1=max(sp(1:n1,1));
sp_var5=max(sp(1:n5,1));
sp_es1=mean(sp(1:n1,1));
sp_es5=mean(sp(1:n1,1));

%% 4
%Compute the unconditional correlation between the two return series.

correlation=corr(sp,oil);

%% 5
%Compute time-varying correlations between the two returns series using a rolling window of 20
%days. Plot your results and comment your findings. What is the average of the time-varying
%correlations?
corr20=zeros(n-20,1);
for i=21:n;
    corr20(i)=corr(sp(1:i),(oil(1:i)));
end
mean20=mean(corr20)*ones(n,1);
figure
plot(dates(2:end,1),[corr20,mean20]);
legend('Correlation on 20 days rolling window','Mean of correlation');

%% 6
er=price2ret(ER);
copper=price2ret(COPPER);
ust10y=price2ret(UST10Y);
er_11_14=er(1176:1958);
copper_11_14=copper(1176:1958);
ust10y_11_14=ust10y(1176:1958);
oil_11_14=oil(1176:1958);
n6=size(oil_11_14,1);
Y=oil_11_14; X=[ones(n6,1),copper_11_14,er_11_14,ust10y_11_14];
forecast_fit=ols(Y,X);  % i downloaded ols function for my econometric course
beta_hat1=forecast_fit.beta;
y_hat=forecast_fit.yhat;
eps=forecast_fit.resid;

%% 7
er_14_17=er(n6+1:end);
copper_14_17=copper(n6+1:end);
ust10y_14_17=ust10y(n6+1:end);
oil_14_17=oil(n6+1:end);
n7=size(er_14_17,1);
for_oil_14_17=beta_hat1(1)*ones(n7,1)+beta_hat1(2)*er_14_17+beta_hat1(3)*copper_14_17+beta_hat1(4)*ust10y_14_17;

figure
Y_14_17f=ret2price(for_oil_14_17);
Y_14_17r=ret2price(oil_14_17);
plot(dates(n6+1:end,1),[Y_14_17r,Y_14_17f]);
legend('Real Return on oil from 2014 to 2017','Forecasted Return on oil from 2014 to 2017');

%% 8 
%Compute the time-varying correlations (using the same 20-day rolling window) between stock
%market returns and demand-related changes in oil prices (i.e., the fitted values from (1)) for
%the sample 6/30/2011 to latest available observation. Plot your results
sp_14_17=sp(n6+1:end);
corr20_8=zeros(n7-20,1);
for i=21:n7;
    corr20_8(i)=corr(sp_14_17(1:i),(for_oil_14_17(1:i)));
end
mean20_8=mean(corr20_8)*ones(n7,1);
figure
plot(dates(end+1-n7:end,1),[corr20_8,mean20_8]);
legend('Correlation on 20 days rolling window','Mean of correlation');

%% 9
%Compute the time-varying correlations (using the same 20-day rolling window) between stock
%market returns and supply-related changes in oil prices (i.e., the residuals values from (1))
%for the sample 6/30/2011 to the latest available observation. Plot your results results and
%comment briefly.

corr20_9=zeros(n6-20,1);
for i=21:n6;
    corr20_9(i)=corr(sp_14_17(1:i),(eps(1:i)));
end
mean20_9=mean(corr20_9)*ones(n6,1);
figure
plot(dates(end-n6+1:end,1),[corr20_9,mean20_9]);
legend('Correlation on 20 days rolling window','Mean of correlation');

%% 10
n10=size(dates(1177:end),1);
vix_11_17=price2ret(VIX(1176:end))+1;
sp_11_17=sp(1176:end);
er_11_17=er(1176:end);
copper_11_17=copper(1176:end);
ust10y_11_17=ust10y(1176:end);
oil_11_17=oil(1176:end);

Y10=oil_11_17; X10=[ones(n10,1),er_11_17,copper_11_17,ust10y_11_17,vix_11_17];

ols_10=ols(Y10,X10);

fit_10=ols_10.yhat;
eps_10=ols_10.resid;


corr20_10_1=zeros(n10-20,1);
for i=21:n10;
    corr20_10_1(i)=corr(sp_11_17(1:i),(fit_10(1:i)));
end
mean20_10_1=mean(corr20_10_1)*ones(n10,1);
figure
plot(dates(end+1-n10:end,1),[corr20_10_1,mean20_10_1]);
legend('Correlation on 20 days rolling window','Mean of correlation');

corr20_10_2=zeros(n10-20,1);
for i=21:n10;
    corr20_10_2(i)=corr(sp_11_17(1:i),(eps_10(1:i)));
end
mean20_10_2=mean(corr20_10_2)*ones(n10,1);
figure
plot(dates(end+1-n10:end,1),[corr20_10_2,mean20_10_2]);
legend('Correlation on 20 days rolling window','Mean of correlation');
