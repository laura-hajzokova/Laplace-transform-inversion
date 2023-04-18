function LT = LT_gamma(a,b,type)
% Laplace transform of PDF or CDF of Gamma distribution with parameters
% a,b, whereas PDF of the distribution Gamma(a,b) is defined as 
% f(t) = b^a*t.^(a-1).*exp(-b*t)/gamma(a)
% 
% SYNTAX:
% LT = LT_gamma(a,b,type)
%
% INPUTS:
% a      - parameter of Gamma distribution
% b      - parameter of Gamma distribution
% type   - type of the otput function (default value is PDF, however the CDF
%          is also possible)

%% check the input parameters
narginchk(2,3)

if nargin<3, type="pdf"; end

%% Laplace transform

if type == "pdf"
    LT = @(t) b^a./(b+t).^a;
elseif type == "cdf"
    LT = @(t) (1./t).*(b^a./(b+t).^a);
end

end