clc;
clear all;
close all;

% INPUT
x = linspace(0,5,100)';
d1 = 4;
d2 = 2;
c = (d2/d1)^(d2/2)/beta(d1/2,d2/2);
a = d2/d1;
m = d1/2 - 1;
n = -(d1+d2)/2;
N = 32;
M = 10;

time_total = zeros(5,1);

%% PDF
for i=1:M
    
pdf_F = @(t) c.*t.^m.*(t+a).^n;
pdf_lt_F = LT_FisherSnedecor(d1,d2);

[pdf_ilt_P,pdf_t_P] = improved_talbot(pdf_lt_F,x,N,'parabola');
pdf_ilt_P(1) = 0;
[pdf_ilt_H,pdf_t_H] = improved_talbot(pdf_lt_F,x,N,'hyperbola');
pdf_ilt_H(1) = 0;
[pdf_ilt_T,pdf_t_T] = improved_talbot(pdf_lt_F,x,N,'talbot');
pdf_ilt_T(1) = 0;
[pdf_ilt_PV,pdf_t_PV] = post_widder(pdf_lt_F,x);
pdf_ilt_PV(1) = 0;
tic
[xx1,pdf_ilt_V] = InvLT(pdf_lt_F,x(1),x(end),length(x));
pdf_t_V = toc;
pdf_ilt_V(1) = 0;

time_total = time_total + [pdf_t_P;pdf_t_H;pdf_t_T;pdf_t_PV;pdf_t_V];
end

%% RESULTS - PDF

contour = ["Parabola        ";"Hyperbola       ";"Talbot cotangent"];
time = time_total/M;
error = [sum(abs(pdf_F(x)-pdf_ilt_P));sum(abs(pdf_F(x)-pdf_ilt_H));sum(abs(pdf_F(x)-pdf_ilt_T));sum(abs(pdf_F(x)-pdf_ilt_PV));sum(abs(pdf_F(xx1)-pdf_ilt_V))];

fprintf('============================================================================== \n');
fprintf('                                 RESULTS (PDF) \n');
fprintf('------------------------------------------------------------------------------ \n');
fprintf('     method          |    contour            |    cpu time   |    error    \n');
fprintf('------------------------------------------------------------------------------ \n');
fprintf('     talbot          |    %s   |    %4.6f   |    %4.4e   \n' ,[contour';time(1:3)';error(1:3)']);
fprintf('------------------------------------------------------------------------------ \n');
fprintf('     post widder     |    %s   |    %4.6f   |    %4.4e   \n' ,["                ";time(4)';error(4)']);
fprintf('------------------------------------------------------------------------------ \n');
fprintf('     hosono          |    %s   |    %4.6f   |    %4.4e   \n' ,["                ";time(5)';error(5)']);
fprintf('------------------------------------------------------------------------------ \n');


%% CDF
% cdf_F = @(t) double(betainc(d1*t./(d1*t + d2),d1/2,d2/2));
% cdf_lt_F = LT_FisherSnedecor(d1,d2,"cdf");
% 
% [cdf_ilt_P,cdf_t_P] = improved_talbot(pdf_lt_F,x,N,'parabola');
% cdf_ilt_P(1) = 0;
% [cdf_ilt_H,cdf_t_H] = improved_talbot(pdf_lt_F,x,N,'hyperbola');
% cdf_ilt_H(1) = 0;
% [cdf_ilt_T,cdf_t_T] = improved_talbot(pdf_lt_F,x,N,'talbot');
% cdf_ilt_T(1) = 0;
% [cdf_ilt_PV,cdf_t_PV] = post_widder(pdf_lt_F,x);
% cdf_ilt_PV(1) = 0;
% tic
% [xx12,cdf_vilt_F] = InvLT(cdf_lt_F,x(1),x(end),length(x));
% cdf_vt_F = toc;
% cdf_vilt_F(1) = 0;

%% RESULTS - CDFs

% time = [cdf_t_chisq;cdf_t_exp;cdf_t_gamma;cdf_t_F];
% error = [sum(abs(cdf_chisq(x)-cdf_ilt_chisq));sum(abs(cdf_exp(x)-cdf_ilt_exp));sum(abs(cdf_gamma(x)-cdf_ilt_gamma));sum(abs(cdf_F(x)-cdf_ilt_F))];
% ptime = [cdf_pt_chisq;cdf_pt_exp;cdf_pt_gamma;cdf_pt_F];
% perror = [sum(abs(cdf_chisq(x)-cdf_pilt_chisq));sum(abs(cdf_exp(x)-cdf_pilt_exp));sum(abs(cdf_gamma(x)-cdf_pilt_gamma));sum(abs(cdf_F(x)-cdf_pilt_F))];
% vtime = [cdf_vt_chisq;cdf_vt_exp;cdf_vt_gamma;cdf_vt_F];
% verror = [sum(abs(cdf_chisq(xx02)-cdf_vilt_chisq));sum(abs(cdf_exp(xx12)-cdf_vilt_exp));sum(abs(cdf_gamma(xx22)-cdf_vilt_gamma));sum(abs(cdf_F(xx32)-cdf_vilt_F))];
% 
% fprintf('============================================================================== \n');
% fprintf('                                 RESULTS (CDFs) \n');
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('        method       |    distribution     |    time       |    error    \n');
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('        talbot       |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';time';error']);
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('     post widder     |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';ptime';perror']);
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('        hosono       |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';vtime';verror']);
% fprintf('------------------------------------------------------------------------------ \n');
