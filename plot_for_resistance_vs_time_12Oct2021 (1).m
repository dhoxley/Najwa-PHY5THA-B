% Lithium fluoride temperature dependent conductivity 
% In this MATLAB script we plotted (Resistance vs time )from data set 1.
% 12 October 2021
% Written by Najwa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
% This file contains Time in (YYYY-MM-DD HH-MM-SS)format in first column.
% Convert time to time_elapsed. 
% Also we have resistance in column 3.
B=importdata('dataset1Oct.txt')
cool=B.data;
time=cool(:,1);
resistance=cool(:,3);
figure('Name','Raw Data');
plot(time,resistance,'.')
% time in seconds since first measurement.
time_elapsed=(time-time(1));
figure('Name','Elapsed time')
plot(time_elapsed,resistance,'.')
xlabel('Time in second');
ylabel('Resistance in Ohms');
title( 'Plot for resistance vs time');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%