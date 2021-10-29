% This MATLAB script to analyse data set 5 .
% First step to plot of resistance of heated LiF vs time . 
% Next step to plot Log (R) vs 1/T plot and fitting for cooling region.
% Then For cooling region we try to finding the gradient value. 
% and activation energies for all the cooling regions.
% 29 October 2021
% Written by Najwa
% Modified by David Hoxley 21 October 2021 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all   % clear varilables in workspace in case
close all   % close all figures
% time constant in seconds obtained from Data Set 
tau=186;
temp_env=300; % environment temp in kelvin
temp_init=450+273; % initial temp of cooling curve in Kelvin- guess from max temp of data set #1
Boltz=8.617e-5  % Boltzmann const in eV/Kelvin 
 % importdata() reads from the file and separates numbers from text
% The numbers are in an array B.data
B= importdata('dataset5rawdata.txt')
% separate out data from header text 
% cool is a 3 col matrix containing the timestamp, elapsedtime, resistance
cool=B.data;       
conductance=cool(:,3); % Separate out into conductance vs time
time=0.1*cool(:,2); % correction factor to align elapsedtime with timestamp 

% plot raw conductance data
figure('Name','Conductance vs time'); 
plot(time,conductance,'.')
xlabel('time (s)')
ylabel('conductance (Ohm-1.m-1)')

% plot conductance vs array position to make breaking up easier
figure('Name','Conductance vs position'); 
plot(conductance, '.')
 
% separate out cooling curve #1
conductance1=[conductance(1:2e4)];
time1=[time(1:2e4)];

% plot conductance vs time
figure('Name','Conductance1 vs time1'); 
plot(time1,conductance1,'.')
xlabel('time (s)')
ylabel('conductance (Ohm-1.m-1)')

% convert time into temperature 
% use Newton's law of cooling
temp1=(temp_init-temp_env)*exp(-time1/tau) + temp_env ;

% % plot temp vs time
% figure('Name','temp1 vs time1'); 
% plot(time1,temp1,'.')
% xlabel('time (s)')  
% ylabel('Temperature (K)')

% plot conductance vs temp for region 1
figure('Name','conductance1 vs temp1'); 
plot(temp1,conductance1,'.')
xlabel('Temperature (K)')
ylabel('conductance (arb units)')
% Separate into intrinsic and extrinsic segments
% intrinsic is hot where conduction dominated by high density of defects with high activation energy
% extrinsic is colder where conduction dominated by fewer impurity defects
% of lower activation energy
intrinsic1=[conductance1(1:100)];
temp1_intrinsic=[temp1(1:100)];
extrinsic1=[conductance1(100:500)];
temp1_extrinsic=[temp1(100:500)];

% plot intrinsic region to check
figure('Name','intrinsic1 vs temp'); 
plot(temp1_intrinsic,intrinsic1,'.')
 
% plot extrinsic region to check
figure('Name','extrinsic1 vs temp'); 
plot(temp1_extrinsic,extrinsic1,'.')

%  prepare for Arrhenius plot intrinsic+extrinsic regions
logG1=log(conductance1(1:500));
invT1=1./temp1(1:500);

% % % make Arrhenius plot whole range
% figure('Name','log G1 vs 1/temp1'); 
% plot(invT1,logG1,'.')

%  prepare for Arrhenius plot intrinsic
logG1_intrinsic=log(intrinsic1);
invT1_intrinsic=1./temp1_intrinsic;

% fit intrinsic Arrhenius
p=polyfit(invT1_intrinsic,logG1_intrinsic,1);
yfit_in=p(1)*invT1_intrinsic+ p(2);
slope_in=p(1);
intercept_in=p(2);

% % make Arrhenius plot intrinsic with fit
% figure('Name','log G1_intrinsic vs 1/temp1'); 
% plot(invT1_intrinsic,logG1_intrinsic,'.')
% hold on
% plot(invT1_intrinsic,yfit_in,'red')
% text(1.38e-3,0.33,['logG_IN=' num2str(p(1)) '* invT + ' num2str(p(2))],'color','red');
% hold off

% prepare for Arrhenius plot extrinsic
logG1_extrinsic=log(extrinsic1);
invT1_extrinsic=1./temp1_extrinsic;

% fit extrinsic Arrhenius
p=polyfit(invT1_extrinsic,logG1_extrinsic,1);
yfit_ex=p(1)*invT1_extrinsic+ p(2);
slope_ex=p(1);
intercept_ex=p(2);

% % make Arrhenius plot extrinsic
% figure('Name','log G1_extrinsic vs 1/temp1'); 
% plot(invT1_extrinsic,logG1_extrinsic,'.')
% hold on
% plot(invT1_extrinsic,yfit_ex,'red')
% text(1.8e-3,0.34,['logG_EX=' num2str(p(1)) '* invT + ' num2str(p(2))],'color','red');
% hold off

% make Arrhenius plot whole range again
% this time with fit overplot
figure('Name','log G1 vs 1/temp1'); 
hold on
plot(invT1,logG1,'.')
plot(invT1_extrinsic,yfit_ex,'red')
% text(1.8e-3,0.33,['logG_EX=' num2str(slope_ex) '* invT + ' num2str(intercept_ex)],'color','red');
plot(invT1_intrinsic,yfit_in,'blue')
% text(1.8e-3,0.34,['logG_IN=' num2str(slope_in) '* invT + ' num2str(intercept_in)],'color','blue');
xlabel( '1/temperature (K^{-1})'  )
ylabel('log(conductance) ')
legend('LiF cooling data', 'fit to extrinsic region', 'fit to intrinsic region')
hold off

disp('Intrincic slope is')
disp(slope_in)
disp('Extrinsic slope is')
disp(slope_ex)

ActivationEnergy_in=-slope_in*Boltz
ActivationEnergy_ex=-slope_ex*Boltz


% check fits and find uncertainties
mdl_in=fitlm(invT1_intrinsic,logG1_intrinsic)
disp('fitlm slope is')
check_slope_in=mdl_in.Coefficients.Estimate(2)
disp('fitlm Standard error is')
SE_in=mdl_in.Coefficients.SE(2)
mdl_ex=fitlm(invT1_extrinsic,logG1_extrinsic)
disp('fitlm slope is')
check_slope_ex=mdl_ex.Coefficients.Estimate(2)
disp('fitlm Standard error is')
SE_ex=mdl_ex.Coefficients.SE(2)


% Display final activation energies with uncertainties
disp('INtrincic Activation in eV is')
disp(ActivationEnergy_in)
disp('absolute uncertainty in INtrinsic activation energy')
uncert_in=ActivationEnergy_in*SE_in/check_slope_in
%2sf_uncert_in=round(uncert_in,2,'significant')
disp('EXtrinsic Activation in eV is')
disp(ActivationEnergy_ex)
disp('absolute uncertainty in EXtrinsic activation energy')
uncert_ex=ActivationEnergy_ex*SE_ex/check_slope_ex
%2sf_uncert_ex=round(uncert_ex,2,'significant')

% y = round(x,2,'significant')