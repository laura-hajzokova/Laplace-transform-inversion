function [ilt,time] = improved_talbot(F,t,N,contour)
%  Numerical Inversion of the Laplace Transform F(s). 
%  improved_talbot is an algorithm based on evaluation of Bromwich integral
%  that can be used to invert any arbitrary Laplace transform (LT).
%  Our specific objective in this context is to apply it to the LT 
%  of probability distributions that can be expressed in terms of their LT.
% 
% SYNTAX:
%    [ilt,time] = improved_talbot(F,t,N,contour)
%
% INPUTS:
%  F        - function handle which evaluates the Laplace transform for any
%             array of complex values as a string 
%  t        - vector of values where the inverted function is evaluated
%  N        - number of nodes used in quadrature (the default value is N=32)
%  contour  - the contour used in contour transformation in Bromwich
%             integral (default value is "talbot", other possible contours are
%             "parabola" and "hyperbola")
%
% REMARK
%  Inversion of Laplace transforms is a very important procedure used in
%  solution of complex linear systems. The function offers a way to achieve
%  the result. It is based on the papers:
%
% Dingfelder, B., Weideman, J.A.C. (2015). An improved Talbot method for 
% numerical Laplace transform inversion, Numerical Algorithms, 68, 167-183
% Trefethen, L. N., Weideman, J. A. C., Schmelzer, T. (2006). Talbot 
% quadratures and rational approximations, BIT Numerical Mathematics, 46, 653–670
%  
% Laura Hajzokova (laura.hajzokova@savba.sk)

%% PARAMETERS

narginchk(2,4);

if nargin<4, contour="talbot"; end
if nargin<3, N=32; end

a = 0.6407;
m = 0.5017;
s = 0.6122;
n = 0.2645;

tic
t = t';
k = 1:(N-1); % k=1:N
theta = -pi + 2*pi/N*(k-1/2);

% meshgrid
[T,THETA] = meshgrid(t,theta);

if contour=="talbot"
    z = (-s + m.*THETA.*cot(a.*THETA) + n*1i.*THETA);
    dz = (m.*cot(a.*THETA) - a*m.*THETA.*csc(a*THETA).^2 + 1i*n);
elseif contour=="hyperbola"
    z = 2.246*(1-sin(1.1721-0.3443*1i.*THETA));
    dz = 2.246*cos(1.1721 - 0.3443*1i.*THETA)*0.3443*1i;
elseif contour=="parabola"
    z = (0.1309 - 0.1194.*THETA.^2 + 0.25*1i*THETA);
    dz = (- 2*0.1194.*THETA + 0.25*1i);
end

% inverse transform
ilt = real((1./(t*1i)).*sum(exp(N*z).*F(N*z./T).*dz));
ilt = ilt';

time = toc;
end