function [mean,var] = GausPDFProd(meanData,varData,meanPrior,varPrior)
% Calculates product of two Gaussian PDFs, namely 'data' and 'prior'
%
% INPUTS
% meanData, varData     - Mean and variance on data
% meanPrior, varPrior   - Mean and variance based on prior
%
% OUTPUTS
% mean                  - Mean of posterior
% var                   - Variance of posterior


ssvar = varData.^2 + varPrior.^2;
var = varData^2.*varPrior^2./ssvar;

A = varData.^(-2).*meanData + varPrior.^(-2).*meanPrior;
B = varData.^(-2)+varPrior.^(-2);
mean = A./B;