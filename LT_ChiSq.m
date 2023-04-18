function LT = LT_ChiSq(df,ncp,coef,type)

% Laplace transform of linear combination of independent RVs from Chi-Square distribution
% 
% BASIC THEORY:
% if X_i ~ N(mu_i,1) then Y_i=(X_i)^2 ~ ChiSquare(df_i,ncp_i)
% ncp_i=(mu_i)^2 and ncp = sum_{i=1}^N ncp_i
% hence if X_i ~ N(0,1) then Y_i=(X_i)^2 ~ ChiSquare(df_i=1,ncp_i=0)
%
% LT_ChiSq evaluates the Laplace transform of  of  Z =  sum_{i=1}^N Y_i
% The Laplace transform of a sum of independent random variables is equal to the
% product of the individual Laplace transforms
%
% Laplace transform of Z is defined by
% L(s) = exp((-ncp*s)/(1+2*sigma^2*s))/((1+2*sigma^2*s)^(df/2))
% under the condition Re(s)>-1/(2*sigma^2)
%
% SYNTAX:
% LT = LT_ChiSq(t,df,ncp,type)
%
% INPUTS:
% df    - vector of the degrees of freedom of the the chi-squared random
%         variables.  If df is scalar, it is assumed that all degrees of
%         freedom are equal. If empty, default value is df = 1.
% ncp   - vector of the non-centrality parameters of the the chi-squared
%         random variables. If ncp is scalar, it is assumed that all
%         non-centrality parameters are equal. If empty, default value is
%         ncp = 0.
% coef  - vector of the coefficients of the linear combination of the
%          chi-squared random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
% type  - type of the otput function (default value is PDF, however the CDF
%         is also possible)

%% check the input parameters
narginchk(1,4);

if nargin < 4, type  = "pdf"; end
if nargin < 3, coef  = 1; end
if nargin < 2, ncp  = 0; end

%% Laplace transform

if type == "pdf"
    LT = @(t) (2*t+1).^(-df/2);
elseif type == "cdf"
    LT = @(t) (1./t).*(2*t + 1).^(-df/2);
end

end