% Lithium fluoride temperature dependent conductivity.
% In this MATLAB script we plotted and fit data for cooling region, 
% fit to find activation energy (Resistance with temperature). 
% 10 October 2021
% Written by Najwa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
x=xlsread('datasheet1.xlsx');
onT=x(:,1);              %1/Temperature values
logonv=x(:,2);           %log 1/v(voltage) values
scatter(onT,logonv,'.')  %plot log 1/v vs1/T
xlabel('1/T' ,'FontSize',15);
ylabel('Log (R)', 'FontSize', 15);
% Constants
B = 8.617333262145e-5;    % Boltzmann constant in eV/K
% Fits polynomial of order 1 (i.e. y=mx+c)
hold on
P = polyfit(onT,logonv,1)
yfit = P(1)*onT+P(2);
% Plot the fit on top of the data
plot(onT,yfit,'red');
text(2e-3,-13,['logonv=' num2str(P(1)) 'onT' num2str(P(2))],'color','black', 'FontSize', 15);
hold off
% The best fit line that can be obtained (R^2=0.9963).