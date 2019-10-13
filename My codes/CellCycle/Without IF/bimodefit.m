%% Fitting bimodality to the DNA Content.

% This does not check if the distribution is indeed bimdodal - there are
% ways to check bimodality of a distribution. You may want to ememd the
% code to your taste.

% clear all;
% close all;
% clc;

function [G1_mean, var1, G2_mean, var2] = bimodefit(filename)
A = load(filename);
x = A(:,1);
bimode = @(x,p,mu1,mu2,sigma1,sigma2) p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);

pStart = .7;
% muStart = [40000000, 80000000];
muStart = quantile(x,[.25 .75]);
sigmaStart = sqrt(var(x) - .25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];


lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];
paramEsts = mle(x, 'pdf',bimode, 'start',start, 'lower',lb, 'upper',ub)
                      
statset('mlecustom')

options = statset('MaxIter',10000, 'MaxFunEvals',6000);
paramEsts = mle(x, 'pdf',bimode, 'start',start, 'lower',lb, 'upper',ub, 'options',options)
                      
% bins = 0:.5:7.5;
% h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
% h = histogram(x,80);
% 
% h.FaceColor = [.9 .9 .9];
% xgrid = linspace(1.1*min(x),1.1*max(x),200);
% pdfgrid = bimode(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
% hold on
% plot(xgrid,pdfgrid,'-'); 
% hold off
% xlabel('x')
% ylabel('Probability Density')

G1_mean = paramEsts(1,2);
var1 = paramEsts(1,3);
G2_mean = paramEsts(1,4);
var2 = paramEsts(1,5);
f = fopen('G1_peak.dat','w');  
    fprintf(f,'%f\t%f\t%f\t%f', G1_mean, var1, G2_mean, var2);

end
