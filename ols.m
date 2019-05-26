function results=ols(y,x,p)
% PURPOSE: least-squares regression 
%---------------------------------------------------
% USAGE: results = ols(y,x)
% where: y = dependent variable vector    (nobs x 1)
%        x = independent variables matrix (nobs x nvar)
%---------------------------------------------------
% RETURNS: a structure
%        results.meth  = 'ols'
%        results.beta  = bhat     (nvar x 1)
%        results.tstat = t-stats  (nvar x 1)
%        results.bstd  = std deviations for bhat (nvar x 1)
%        results.yhat  = yhat     (nobs x 1)
%        results.resid = residuals (nobs x 1)
%        results.sige  = e'*e/(n-k)   scalar
%        results.rsqr  = rsquared     scalar
%        results.rbar  = rbar-squared scalar
%        results.dw    = Durbin-Watson Statistic
%        results.nobs  = nobs
%        results.nvar  = nvars
%        results.y     = y data vector (nobs x 1)
%        results.bint  = (nvar x2 ) vector with 95% confidence intervals on beta
%---------------------------------------------------
% SEE ALSO: prt(results), plt(results)
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Barry Dillon (CICG Equity)
% added the 95% confidence intervals on bhat

if (nargin > 3); error('Wrong # of arguments to ols'); 
else
 [nobs nvar] = size(x); [nobs2 junk] = size(y);
 if (nobs ~= nobs2); error('x and y must have same # obs in ols'); 
 end;
end;

results.meth = 'ols';
results.y = y;
results.nobs = nobs;
results.nvar = nvar;


xpxi = (x'*x)\eye(nvar);

T=size(x,1);
results.beta = xpxi*(x'*y);
results.yhat = x*results.beta;
results.resid = y - results.yhat;
sigu = results.resid'*results.resid;
results.sige = sigu/(nobs-nvar);
tmp = (results.sige)*(diag(xpxi));
results.cov=(results.sige)*(xpxi);
sigb=sqrt(tmp);
results.bstd = sigb;
results.tstat = results.beta./(sqrt(tmp));
results.pvalue  = 2*(1-tcdf( abs(results.tstat), T-size(results.beta,1) ));
ym = y - mean(y);
rsqr1 = sigu;
rsqr2 = ym'*ym;
results.rsqr = 1.0 - rsqr1/rsqr2; % r-squared
rsqr1 = rsqr1/(nobs-nvar);
rsqr2 = rsqr2/(nobs-1.0);
if rsqr2 ~= 0
    results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared
else
    results.rbar = results.rsqr;
end;
ediff = results.resid(2:nobs) - results.resid(1:nobs-1);
results.dw = (ediff'*ediff)/sigu; % durbin-watson


%% RESET TEST
% Default, powers 2 and 3.
if nargin==2;
    p=3;
end;
pw=[2:p];

Z=results.yhat.^pw;


W=[x,Z]; % extended set of regressors

K=size(x,2); % number of original regressors

L=size(W,2)-K; % number of powers of y_hat

R=[zeros(L,K), eye(L)]; % generate the selection matrix of restrictions

c=zeros(L,1);  % vector of constants 

gamma_hat=inv(W'*W)*W'*y; % estimates on the auxiliary extended regression

resid_aux=y-W*gamma_hat;  % residuals auxiliary component

s2_aux=resid_aux'*resid_aux/(T-K-L);

results.RESET_test=(R*gamma_hat-c)'*inv(R*inv(W'*W)*R')*(R*gamma_hat-c)/(L*s2_aux); % RESET test-statistic

results.pvalue_RESET_test=1-fcdf(results.RESET_test,L,T-K-L);  % p-value (all probability of the F on the right tail) - pvalue>0.05 : the linear specification is ok
%% Print on screen

parno   = (1:size(results.beta,1))';

Res     = [ parno results.beta results.bstd results.tstat results.pvalue];

 disp('OLS estimation');
 
fprintf('\n\n\n **********************************************************************\n');
if T-size(results.beta,1)<=0;
    fprintf('\nWarning\n')
    fprintf('Model contains more parameters than observations \n')
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
end
fprintf('Number of observations: %12.4f\n',T);
fprintf('Number of parameters    %12.4f\n',size(results.beta,1));
fprintf(' **********************************************************************\n');
fprintf('       parameter       beta        stderr    t-student      p-value\n');
fprintf('  %12.0f %12.4f  %12.4f %12.4f %12.4f\n', Res' );