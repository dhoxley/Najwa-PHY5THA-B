% This MATLAB script for Data set (3) analysis.
% This file contains tow columns , one for time and the other 
% temperature.
% 22 October 2021
% Written by Najwa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Percent (Importdata) implies that the data is read from a file and separated into numbers and text.

clear 
clc
% There are 3 cool matrix containing (resistance , time and temperature).
% The numbers are in an array B.data.
B=importdata('datasheet3_20oct.xlsx');
cool=B.data;
temp=cool(:,2);
time=cool(:,1);
figure('Name','Raw Data');
plot(time,temp,'.')

% convertthe temp to Kelvin
ktemp=temp+273
figure('Name','cooling curve in Kelvin')
plot(time,ktemp,'.')
xlabel('Time in Seconds');
ylabel('Temperature in Kelvin')
title('Plot for the temperature vs time (cooling curve)');

% Environment temperature in kelvin.
rel_ktemp=ktemp-300
figure('Name','T-Tenv(kelvin)vstime (sec)');
plot(time,rel_ktemp,'.')
log_rel_ktemp=log(rel_ktemp)
figure('Name','log(T-Tenv)vs time');
plot(time,log_rel_ktemp,'.')
P=polyfit(time,log_rel_ktemp,1)
disp('gradient is')
gradient=P(1)
disp('time constant is')
tau=-1/P(1)
yfit=P(1)*time+P(2)
hold on 
plot(time,yfit,'r-.')
R=corrcoef(time,log_rel_ktemp)
disp('PearsonRsquared value is')
rsq=R(1,2)*R(1,2)
str1=sprintf('%.5f',rsq);
disp(rsq)
txt1={'R Squared=',str1};
text(1.04e4,5,txt1)
tr2=sprintf('%.5f',tau);
% Newton's law of cooling is 
% Log(T-Tenv)=log(T0-Tenv)-t/tau
% Tau=1/gradient of graph
% Calculate and display time constant
disp(tau)
txt2={'Time constant(s):',tr2};
text(1.01e4,3.5,txt2)
xlabel('Time in Seconds');
ylabel('Log(temperature)')
title('Plot for log (temperature) vs time (cooling curve) ');