%% Fitting Gaussian to the S phase cells.

% This does not check if the distribution is indeed bimdodal - there are
% ways to check bimodality of a distribution. You may want to ememd the
% code to your taste.

% clear all;
% close all;
% clc;

function [G2_mean, std2] = SCellFit(filename)

if exist('G2_M.dat', 'file')

A = load('G2_M.dat');
x = A(:,1);

else
    CellCycleStage(filename);
end
    
unimode = @(x,mu1,sigma1) normpdf(x,mu1,sigma1);

muStart = mean(x);
sigmaStart = std(x);
start = [muStart sigmaStart];


lb = [0 0];
ub = [Inf Inf];
                      
statset('mlecustom')

options = statset('MaxIter',20000, 'MaxFunEvals',40000);
paramEsts = mle(x, 'pdf',unimode, 'start',start, 'lower',lb, 'upper',ub, 'options',options);
                      

G2_mean = paramEsts(1,1);
std2 = paramEsts(1,2);

    

    f = fopen('G2_peak.dat','w');  
    fprintf(f,'%f\t%f\n', G2_mean, std2);

end