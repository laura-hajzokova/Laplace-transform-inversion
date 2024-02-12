function [x,fun]=InvLT(LTfun,xmin,xmax,n,a,ns,nd)
% INVLAP – Numerical Inversion of the Laplace Transform LTfun(s). 
%  InvLT is a general-purpose algorithm that can be used to invert any
%  arbitrary Laplace transform (LT), our specific objective in this context
%  is to apply it to the LT of probability distributions that can be
%  expressed in terms of their LT.
% 
% SYNTAX:
%    [x,fun]=InvLT(LTfun,xmin,xmax,n,a,ns,nd)
%
% INPUTS:
%  LTfun - function handle which evaluates the Laplace transform for any
%          array of complex values as a string 
%  xmin  - lower limit of x-values where the inverted function is evaluated
%  xmax  - upper limit of x-values where the inverted function is evaluated
%  n     - total number of values where the inverted function is evaluated
%  a     - parameter of the method (default value is a = 6)
%  ns    - parameter of the method (default value is ns = 20)
%  nd    - parameter of the method (default value is ns = 19)
% 
% EXAMPLE1
%  [x,fun]=InvLT(@(s)s./(s.^2+4*pi^2),0,10,201);
%  plot(x,fun,'.-')
% 
% EXAMPLE2
%  [x,fun]=InvLT(@(s)tanh(s)./s,0,10,201);
%  plot(x,fun,'.-')
% 
% EXAMPLE3
%  [x,fun]=InvLT(@(s)1./(sqrt(s).*s),0,100,201);
%  plot(x,fun,'.-')
%
% EXAMPLE4
%  [x,fun]=InvLT(@(s) 1./(1+2*s).^(5/2),0,20,201);
%  plot(x,fun,'.-')
%
% EXAMPLE5
%  df = 5;
%  [x,fun]=InvLT(@(s) 1./(1+2*s).^(5/2),0,20,201);
%  plot(x,fun,'.-')
%  grid on, hold on
%  plot(x,chi2pdf(x,df),'--')
%
% EXAMPLE6
%  df = 5;
%  [x,fun]=InvLT(@(s) 1./(s.*(1+2*s).^(df/2)),0,20,201);
%  plot(x,fun,'.-')
%  grid on, hold on
%  plot(x,chi2cdf(x,df),'--')
%
% EXAMPLE7 (PDF of a linear combination of independent RVs)
%  LTgen = @(s,df) 1./(1+2*s).^(df/2);
%  df = [4 5 7];
%  LTfun = @(s) LTgen(s/3,df(1)) .* LTgen(s/3,df(2)) .* LTgen(s/3,df(3));
%  [x,fun]=InvLT(LTfun,0,15,201);
%  plot(x,fun,'.-'), grid on
%  xlabel('x')
%  ylabel('pdf')
%  title('PDF of a linear combination of independent RVs')
%
% EXAMPLE8 (CDF of a linear combination of independent RVs)
%  LTgen = @(s,df) 1./(1+2*s).^(df/2);
%  df = [4 5 7];
%  LTfun = @(s) (1./s) .* LTgen(s/3,df(1)) .* LTgen(s/3,df(2)) .* LTgen(s/3,df(3));
%  [x,fun]=InvLT(LTfun,0,15,201);
%  plot(x,fun,'.-'), grid on
%  xlabel('x')
%  ylabel('cdf')
%  title('CDF of a linear combination of independent RVs')
%
% InvLT is a version of the algorithm INVLAP proposed by 
%  Juraj Valsa (2023). Numerical Inversion of Laplace Transforms in Matlab
%  (https://www.mathworks.com/matlabcentral/fileexchange/ ...
%  32824-numerical-inversion-of-laplace-transforms-in-matlab),
%  MATLAB Central File Exchange. Retrieved March 3, 2023.  
%
% REMARK
%  Inversion of Laplace transforms is a very important procedure used in
%  solution of complex linear systems. The function f(t)=INVLAP(F(s))
%  offers a simple, effective and reasonably accurate way to achieve the
%  result. It is based on the paper:
%
%   J. Valsa and L. Brancik: Approximate Formulae for Numerical Inversion
%   of Laplace Transforms, Int. Journal of Numerical Modelling: Electronic
%   Networks, Devices and Fields, Vol. 11, (1998), pp. 153-166 
%  
%  The transform LTfun(s) may be any reasonable function of complex
%  variable s^α, where α is an integer or non-integer real exponent. Thus,
%  the function INVLAP can solve even fractional problems and invert
%  functions LTfun(s) containing rational, irrational or transcendental
%  expressions. The function does not require to compute poles nor zeroes
%  of LTfun(s). It is based on values of LTfun(s) for selected complex
%  values of the independent variable s. The resultant computational error
%  can be held arbitrarily low at the cost of CPU time. With the today’s
%  computers and their speed this does not present any serious limitation
%  (see Examples).
% 
% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '03-Mar-2023 17:15:47'

%% CHECK THE INPUTS AND OUTPUTS
narginchk(1, 7);
if nargin < 7, nd = []; end
if nargin < 6, ns = []; end
if nargin < 5, a = []; end
if nargin < 4, n = []; end
if nargin < 3, xmax = []; end
if nargin < 2, xmin = []; end

if isempty(xmin), xmin = 0; end
if isempty(xmax), xmax = 0; end
if isempty(n), n = 100; end
if isempty(a), a = 6; end
if isempty(ns), ns = 20; end
if isempty(nd), nd = 19; end

% t=0 is not allowed
if xmin == 0  
    xmin = eps^(2/3);  
end

% time vector
x = linspace(xmin,xmax,n)'; 

% prepare necessary coefficients
ni  = 1:ns+1+nd;
aa  = a + (ni-1) .* pi .* 1i;
bb  = -exp(a) .* (-1).^ni;
n   = 1:nd;

bdif = fliplr(cumsum(gamma(nd+1)./gamma(nd+2-n)./gamma(n)))./2^nd;
bb(ns+2:ns+1+nd) = bb(ns+2:ns+1+nd).*bdif;
bb(1) = bb(1)/2;

% complex frequency s
s   = aa ./ x;                 
bt  = bb ./ x;
szs = size(s);

% functional value F(s)
btF = reshape(bt(:) .* LTfun(s(:)),szs);          

% original f(tt)
fun  = sum(real(btF),2);  

