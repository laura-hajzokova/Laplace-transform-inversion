clc;
clear all;
close all;

% INPUT
x = linspace(0,5,100)';
timeChi = zeros(3,1);
timeExp = zeros(3,1);
timeGamma = zeros(3,1);
timeF = zeros(3,1);
M = 10;

for i=1:M
    
%% CHI-SQUARE DISTRIBUTION
df = 5;

pdf_chisq = @(t) 1/(2^(df/2)*gamma(df/2)) * t.^(df/2 - 1).*exp(-t./2);
pdf_lt_chisq = LT_ChiSq(df);
[pdf_ilt_chisq,pdf_t_chisq] = improved_talbot(pdf_lt_chisq,x);
pdf_ilt_chisq(1) = 0;
[pdf_pilt_chisq,pdf_pt_chisq] = post_widder(pdf_lt_chisq,x);
pdf_pilt_chisq(1) = 0;
tic; 
[xx01,pdf_vilt_chisq] = InvLT(pdf_lt_chisq,x(1),x(end),length(x)); 
pdf_vt_chisq = toc;
pdf_vilt_chisq(1) = 0;

timeChi = timeChi + [pdf_t_chisq;pdf_pt_chisq;pdf_vt_chisq];

% cdf_chisq = @(x) chi2cdf(x,df);
% cdf_lt_chisq = LT_ChiSq(df,0,1,"cdf");
% [cdf_ilt_chisq,cdf_t_chisq] = improved_talbot(cdf_lt_chisq,x);
% cdf_ilt_chisq(1) = 0;
% [cdf_pilt_chisq,cdf_pt_chisq] = post_widder(cdf_lt_chisq,x);
% cdf_pilt_chisq(1) = 0;
% tic
% [xx02,cdf_vilt_chisq] = InvLT(cdf_lt_chisq,x(1),x(end),length(x));
% cdf_vt_chisq = toc;
% cdf_vilt_chisq(1) = 0;

%% F DISTRIBUTION
d1 = 4;
d2 = 2;
c = (d2/d1)^(d2/2)/beta(d1/2,d2/2);
a = d2/d1;
m = d1/2 - 1;
n = -(d1+d2)/2;
 
pdf_F = @(t) c.*t.^m.*(t+a).^n;
pdf_lt_F = LT_FisherSnedecor(d1,d2);
[pdf_ilt_F,pdf_t_F] = improved_talbot(pdf_lt_F,x);
pdf_ilt_F(1) = 0;
[pdf_pilt_F,pdf_pt_F] = post_widder(pdf_lt_F,x);
pdf_pilt_F(1) = 0;
tic
[xx11,pdf_vilt_F] = InvLT(pdf_lt_F,x(1),x(end),length(x));
pdf_vt_F = toc;
pdf_vilt_F(1) = 0;

timeF = timeF + [pdf_t_F;pdf_pt_F;pdf_vt_F];

% cdf_F = @(t) double(betainc(d1*t./(d1*t + d2),d1/2,d2/2));
% cdf_lt_F = LT_FisherSnedecor(d1,d2,"cdf");
% [cdf_ilt_F,cdf_t_F] = improved_talbot(cdf_lt_F,x);
% cdf_ilt_F(1) = 0;
% [cdf_pilt_F,cdf_pt_F] = post_widder(cdf_lt_F,x);
% cdf_pilt_F(1) = 0;
% tic
% [xx12,cdf_vilt_F] = InvLT(cdf_lt_F,x(1),x(end),length(x));
% cdf_vt_F = toc;
% cdf_vilt_F(1) = 0;
% 

%% EXPONENTIAL DISTRIBUTION
lambda = 1;

pdf_exp = @(t) lambda*exp(-lambda*t);
pdf_lt_exp = LT_exp(lambda);
[pdf_ilt_exp,pdf_t_exp] = improved_talbot(pdf_lt_exp,x);
pdf_ilt_exp(1) = 1;
[pdf_pilt_exp,pdf_pt_exp] = post_widder(pdf_lt_exp,x);
pdf_pilt_exp(1) = 1;
tic
[xx21,pdf_vilt_exp] = InvLT(pdf_lt_exp,x(1),x(end),length(x));
pdf_vt_exp = toc;
pdf_vilt_exp(1) = 1;

timeExp = timeExp + [pdf_t_exp;pdf_pt_exp;pdf_vt_exp];

% cdf_exp = @(t) 1-exp(-lambda*t);
% cdf_lt_exp = LT_exp(lambda,"cdf");
% [cdf_ilt_exp,cdf_t_exp] = improved_talbot(cdf_lt_exp,x);
% cdf_ilt_exp(1) = 0;
% [cdf_pilt_exp,cdf_pt_exp] = post_widder(cdf_lt_exp,x);
% cdf_pilt_exp(1) = 0;
% tic
% [xx22,cdf_vilt_exp] = InvLT(cdf_lt_exp,x(1),x(end),length(x));
% cdf_vt_exp = toc;
% cdf_vilt_exp(1) = 0;

%% GAMMA DISTRIBUTION
a = 2;
b = 1/2;

pdf_gamma = @(t) b^a*t.^(a-1).*exp(-b*t)/gamma(a);
pdf_lt_gamma = LT_gamma(a,b);
[pdf_ilt_gamma,pdf_t_gamma] = improved_talbot(pdf_lt_gamma,x);
pdf_ilt_gamma(1) = 0;
[pdf_pilt_gamma,pdf_pt_gamma] = post_widder(pdf_lt_gamma,x);
pdf_pilt_gamma(1) = 0;
tic
[xx31,pdf_vilt_gamma] = InvLT(pdf_lt_gamma,x(1),x(end),length(x));
pdf_vt_gamma = toc;
pdf_vilt_gamma(1) = 0;

timeGamma = timeGamma + [pdf_t_gamma;pdf_pt_gamma;pdf_vt_gamma];

% cdf_gamma = @(x) gamcdf(x,a,1/b);
% cdf_lt_gamma = LT_gamma(a,b,"cdf");
% [cdf_ilt_gamma,cdf_t_gamma] = improved_talbot(cdf_lt_gamma,x);
% cdf_ilt_gamma(1) = 0;
% [cdf_pilt_gamma,cdf_pt_gamma] = post_widder(cdf_lt_gamma,x);
% cdf_pilt_gamma(1) = 0;
% tic
% [xx32,cdf_vilt_gamma] = InvLT(cdf_lt_gamma,x(1),x(end),length(x));
% cdf_vt_gamma = toc;
% cdf_vilt_gamma(1) = 0;
end

timeChi = timeChi/10;
timeExp = timeExp/10;
timeGamma = timeGamma/10;
timeF = timeF/10;

%% RESULTS - PDFs

distribution = ["Chi-Square    ";"Exponential   ";"Gamma   	    ";"FisherSnedecor"];

error = [sum(abs(pdf_chisq(x)-pdf_ilt_chisq));sum(abs(pdf_exp(x)-pdf_ilt_exp));sum(abs(pdf_gamma(x)-pdf_ilt_gamma));sum(abs(pdf_F(x)-pdf_ilt_F))];
perror = [sum(abs(pdf_chisq(x)-pdf_pilt_chisq));sum(abs(pdf_exp(x)-pdf_pilt_exp));sum(abs(pdf_gamma(x)-pdf_pilt_gamma));sum(abs(pdf_F(x)-pdf_pilt_F))];
verror = [sum(abs(pdf_chisq(xx01)-pdf_vilt_chisq));sum(abs(pdf_exp(xx11)-pdf_vilt_exp));sum(abs(pdf_gamma(xx21)-pdf_vilt_gamma));sum(abs(pdf_F(xx31)-pdf_vilt_F))];

time = [timeChi(1);timeExp(1);timeGamma(1);timeF(1)];
ptime = [timeChi(2);timeExp(2);timeGamma(2);timeF(2)];
vtime = [timeChi(3);timeExp(3);timeGamma(3);timeF(3)];

fprintf('============================================================================== \n');
fprintf('                                 RESULTS (PDFs) \n');
fprintf('------------------------------------------------------------------------------ \n');
fprintf('        method       |    distribution     |    cpu time   |    error    \n');
fprintf('------------------------------------------------------------------------------ \n');
fprintf('        talbot       |    %s   |    %4.6f   |    %4.4e   \n' ,[distribution';time';error']);
fprintf('------------------------------------------------------------------------------ \n');
fprintf('     post widder     |    %s   |    %4.6f   |    %4.4e   \n' ,[distribution';ptime';perror']);
fprintf('------------------------------------------------------------------------------ \n');
fprintf('        hosono       |    %s   |    %4.6f   |    %4.4e   \n' ,[distribution';vtime';verror']);
%fprintf('------------------------------------------------------------------------------ \n');

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
% fprintf('        method       |    distribution     |    cpu time   |    error    \n');
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('        talbot       |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';time';error']);
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('     post widder     |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';ptime';perror']);
% fprintf('------------------------------------------------------------------------------ \n');
% fprintf('        hosono       |    %s   |    %4.6f   |    %4.8e   \n' ,[distribution';vtime';verror']);
% fprintf('------------------------------------------------------------------------------ \n');

%% Plot

% f = figure;
% f.Position = [350 50 850 700];
% 
% figure(f)
% subplot(3,2,1)
% plot(x,pdf_lt_chisq(x),'b-')
% hold on
% plot(x,pdf_ilt_chisq,'r-')
% title('Laplace transform and inverse - PDF of Chi-Square')
% legend('LT','iLT')
% 
% subplot(3,2,2)
% plot(x,cdf_lt_chisq(x),'b-')
% hold on
% plot(x,cdf_ilt_chisq,'r-')
% title('Laplace transform and inverse - CDF of Chi-Square')
% legend('LT','iLT')
% ylim([0 2])
% 
% subplot(3,2,3)
% plot(x,pdf_lt_exp(x),'b-')
% hold on
% plot(x,pdf_ilt_exp,'r-')
% title('Laplace transform and inverse - PDF of Exponential')
% legend('LT','iLT')
% 
% subplot(3,2,4)
% plot(x,cdf_lt_exp(x),'b-')
% hold on
% plot(x,cdf_ilt_exp,'r-')
% title('Laplace transform and inverse - CDF of Exponential')
% legend('LT','iLT')
% ylim([0 2])
% 
% subplot(3,2,5)
% plot(x,pdf_lt_gamma(x),'b-')
% hold on
% plot(x,pdf_ilt_gamma,'r-')
% title('Laplace transform and inverse - PDF of Gamma')
% legend('LT','iLT')
% 
% subplot(3,2,6)
% plot(x,cdf_lt_gamma(x),'b-')
% hold on
% plot(x,cdf_ilt_gamma,'r-')
% title('Laplace transform and inverse - CDF of Gamma')
% legend('LT','iLT')
% ylim([0 2])
