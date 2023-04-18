function LT = LT_exp(lambda,type)
% Laplace transform of PDF or CDF of Exponential distribution with
% parameter lambda.
% 
% SYNTAX:
% LT = LT_exp(lambda,type)
%
% INPUTS:
% lambda - Parameter of exponential distribution
% type   - type of the otput function (default value is PDF, however the CDF
%          is also possible)


%% check the input parameters
narginchk(1,2);

if nargin < 2, type  = "pdf"; end

%% Laplace transform

if type == "pdf"
    LT = @(t) lambda./(lambda + t);
elseif type == "cdf"
    LT = @(t) (1./t).*(lambda./(lambda + t));
end

end