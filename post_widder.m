function [fun,time] = post_widder(f,t)
%  Numerical Inversion of the Laplace Transform f(s). 
%  post_widder is an algorithm based on Post-Widder inversion formula
%  that can be used to invert any arbitrary Laplace transform (LT).
%  Our specific objective in this context is to apply it to the LT 
%  of probability distributions that can be expressed in terms of their LT.
% 
% SYNTAX:
%    [fun,time] = post_widder(f,t)
%
% INPUTS:
%  f   - function handle which evaluates the Laplace transform for any
%          array of complex values as a string 
%  t   - vector of values where the inverted function is evaluated
%
% InvLT is a version of the algorithm INVLAP proposed by 
%  Juraj Valsa (2023). Numerical Inversion of Laplace Transforms in Matlab
%  (https://www.mathworks.com/matlabcentral/fileexchange/ ...
%  32824-numerical-inversion-of-laplace-transforms-in-matlab),
%  MATLAB Central File Exchange. Retrieved March 3, 2023.  
%
% REMARK
%  Inversion of Laplace transforms is a very important procedure used in
%  solution of complex linear systems. The function offers a way to achieve
%  the result. It is based on the paper:
%
%   Abate J., Whitt, W. (1995). Numerical Inversion of Laplace Transforms 
%   of Probability Distributions, ORSA Journal on Computing, 7 (1), 36-43
%  
% Laura Hajzokova (laura.hajzokova@savba.sk)

tic

m = 6;
g = [];
t = t';

for i=1:m
    n = 10*i;
    e = 8;
    r = 1/(10^(e/(2*n)));
    u = (n+1)./(t*2*n*r^n);
    h = pi/n;

    sum1 = 0;

    for j=1:n-1
        s = (n+1)*(1 - r*exp(1i*h*j))./t;
        sum1 = sum1 + (-1)^j * real(f(s));
        
    end
    sum1 = real(f((n+1)*(1-r)./t) + 2*sum1 + (-1)^n*f((n+1)*(1+r)./t));
    g = [g;u.*sum1];
end

w = [];
for k=1:m
    w = [w;(-1)^(m-k) * k^m /(factorial(k)*factorial(m-k))];
end

fun = w'*g;
fun = fun';

time = toc;

end