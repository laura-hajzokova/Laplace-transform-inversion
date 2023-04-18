function LT = LT_FisherSnedecor(d1,d2,type)
% Laplace transform of PDF or CDF of F-distribution with d1 and d2 degrees
% of freedom, i.e. F(d1,d2).
% 
% SYNTAX:
% LT = LT_FisherSnedecor(d1,d2,type)
%
% INPUTS:
% d1, d2 - degrees of freedom of F-distribution F(d1,d2)
% type   - type of the otput function (default value is PDF, however the CDF
%          is also possible)

%% check the input parameters
narginchk(2,3);

if nargin < 3, type  = "pdf"; end

%% Laplace transform

if type == "pdf"
    LT = @(t) gamma((d1+d2)/2)/gamma(d2).*HypergeomU(d1/2,1-d2/2,d2/d1*t);
    %LT = @(t) gamma(d1/2)/beta(d1/2,d2/2).*HypergeomU(d1/2,1-d2/2,d2/d1*t);
elseif type == "cdf"
    LT = @(t) (1./t) .* gamma((d1+d2)/2)/gamma(d2) .*HypergeomU(d1/2,1-d2/2,d2/d1 *t);
end

end