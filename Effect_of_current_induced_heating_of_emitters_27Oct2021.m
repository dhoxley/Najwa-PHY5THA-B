% This MATLAB script builds a simple model of current-induced heating. 
% a: Starting with derive a formula for the resistance of a cylinder of radius( R) 
% and length( Z )made of material of conductivity sigma.
% b: A current( I )passing through this cylinder will dissipated power in Joules/second 
% (P = I^2R ). 
% c: combine this with the formula from part a above to derive a formula expressing P 
% in terms of I, R, Z and sigma the specific heat capacity C of a material is defined as the amount of energy required to raise the temperature of 1 kg of the material by one degree. 
% Thus the temperature rise of an object deltaT is related to the amount of heat Q pumped into it by Q=mass*C*deltaT. use this and your answers to a and b above to derive a formula for the rise in temperature as a function of time
% for a cylinder having a current pass through it.
% Then I plotted  this as a graph, using values for a 1mm diameter and 10 mm long cylinder of LiF. Take the specific heat capacity of LiF as 1562 J Kg-1
% Rise in temperature as a function of time
% 27 October 2021
% Written by Najwa


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimensions of the cylinder.
% r= radius and z= length of the cylinder in m.
% These dimensions are chosen because they are the same as the dimensions of the pillars.
r=5e-4; z=0.01;          
% Cross-sectional area of cylinder in order to find resistance that depends
% on resistivity, length and area.
A=pi*r^2;

I=0.1;                      % Current in Amps 
t=0:1000                    % Time in Seconds
% Resistivity in ohm.m for lithium fluoride.
roh_liF=10               
% Resistance in Ohm. for LiF.
R_liF=roh_liF* z/A
% Power dissipated in Joules/second. 
P_LiF=I*I*R_liF;
% Mass of LiF in Kg.
m_liF=2.06e-5    
% Specific heat capacity of LiF.
c_liF=1562 % J/(kg K)
% Change in temperature depends on the heat supply by the following formula.  
% Q=m*c*deltaT; 
% And power dissipated is defined as heat supply per unit time; 
% i.e P=Q/t.
% Now rise in temperature as a function of time for a cylinder.
deltaT_LiF=P_LiF*t/(m_liF*c_liF) 
plot(t,deltaT_LiF,'red',"LineWidth",2)
hold on
x1=0; x2=1000;y1=3.95e7; y2=3.95e7;
plot([x1,x2,],[y1,y2],'--', 'color', 'red');
xlabel('Timerange in second ');
ylabel('Rise in temperature in Kelvin')
legend('rise in temperature','melting point');
title('Plot for rise in temperature as a function of time for lithium fluoride');
hold off
% Melting point of LiF =276.9 K % melting point means evaluate how long the
% cylinder takes to melt.
% Now repeating the above calculation for Si
% Resistivity in ohm.m for silicon.
roh_Si=6.40e2 
% Mass on silicon in Kg.
m_Si=1.83e-5
% Specific heat capacity of silicon in J/(kg K)
c_Si=710 ;
% Resistance of silicon in ohm.
R_Si=roh_Si* z/A
% Power dissipated in joule/second for silicon
P_Si=I*I*R_Si;
% Change in temperature depends on the heat supply by the following formula.  
% Q=m*c*deltaT; 
% And power dissipated is defined as heat supply per unit time; 
% i.e P=Q/t.
% Now rise in temperature.
deltaT_Si=P_Si*t/(m_Si*c_Si) 
plot(t,deltaT_Si,'blue',"LineWidth",2)
hold on
x1=0; x2=1000;y1=6.3e9; y2=6.3e9;
plot([x1,x2,],[y1,y2],'--', 'color', 'blue');
xlabel('Timerange');
ylabel('Rise in temperature')
legend('rise in temperature','melting point');
hold off
xlabel('Timerange in second ');
ylabel('Rise in temperature in Kelvin')
title('Plot for rise in temperature as a function of time for silicon');
% Melting point of silicon=279.3 K % melting point means evaluate how long the
% cylinder takes to melt.
% Now repeating the above calculation for copper.
roh_Cu=1.68e-8      % Resistivity in ohm.m for copper.
m_Cu=  7.03e-5      % kg 
c_Cu=385            % J/(kg K)
R_Cu=roh_Cu* z/A
P_Cu=I*I*R_Cu;
% Q=m*c*deltaT; % m= mass and c=specific heat capacity.
% Now rise in temperature.
deltaT_Cu=P_Cu*t/(m_Cu*c_Cu) 
plot(t,deltaT_Cu,'black',"LineWidth",2)
hold on
x1=0; x2=1000;y1=0.079; y2=0.079;
plot([x1,x2,],[y1,y2],'--', 'color', 'black');
xlabel('Timerange in second');
ylabel('Rise in temperature in Kelvin')
legend('rise in temperature','melting point');
hold off
xlabel('Timerange in second ');
ylabel('Rise in temperature in Kelvin')
%title('Plot for rise in temperature as a function of time for copper');

% Melting point of copper=273.08 K % melting point means evaluate how long the
% cylinder takes to melt.
% Now repeating the above calculation for tungsten.
roh_w=5.6e-8        % Resistivity in ohm.m for tungsten.
m_w= 1.54e-4        % kg 
c_w=133.976         % J/(kg K)
R_w=roh_w* z/A
P_w=I*I*R_w;
% Q=m*c*deltaT; % m= mass and c=specific heat capacity.
% Now raise in temperature.
deltaT_w=P_w*t/(m_w*c_w) 
plot(t,deltaT_w,'magenta',"LineWidth",2)
hold on
x1=0; x2=1000;y1=0.345; y2=0.345;
plot([x1,x2,],[y1,y2],'--', 'color', 'magenta');
xlabel('Timerange');
ylabel('Rise in temperature')
legend('rise in temperature','melting point');
title('Plot for rise in temperature as a function of time for lithium fluoride');
hold off
xlabel('Timerange');
ylabel('Rise in temperature')
title('Plot for rise in temperature as a function of time for tungsten');
% Melting point of tungsten =353 K % melting point means evaluate how long the
% cylinder takes to melt.
% Plot all four on same axes for comparison
figure()
plot(t,deltaT_LiF,'red',"LineWidth",2)
hold on
plot(t,deltaT_Si,'blue',"LineWidth",2)
hold on
plot(t,deltaT_Cu,'black',"LineWidth",2)
hold on
plot(t,deltaT_w,'magenta',"LineWidth",2)
hold off
xlabel('Timerange in Seconds');
ylabel('Rise in temperature in K')
title('Plot for rise in temperature (time) for LiF, silicon, copper and tungsten');
legend('red for LiF','blue for Si','black for Cu','green for w')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this part of the MATLAB script we are looking for the
% effect of Current-induced heating of emitters.

% i. Temperature as function of time T(t) for cylindrical pillar of same geometry as micropillars, assuming conductivity constant with temperature
% ii. T(t) for cylindrical pillar assuming conductivity decreases linearly with temperature

% iii. T(t) for cylinder where conductivity increases linearly with time

% 1. Using Conductivity Values for intrinsic, doped Si

% 2. Using Conductivity Values of LiF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% Dimensions of the cylinder.
% r= radius and z= length of the cylinder in m.
% These dimensions are chosen because they are the same as the dimensions of the pillars.
r=5e-4; z=0.01;          
% Cross-sectional area of cylinder in order to find resistance that depends
% on resistivity,length and area.
A=pi*r^2;

I=0.1;                      % Current in Amps 
T0=0                       %initial temperature in K
t=0:10                    % Time in Seconds
% Resistivity in ohm.m for lithium fluoride.
roh_liF=10               
% Resistance in Ohm. for LiF.
R0_liF=roh_liF* z/A
R1_liF=R0_liF*15e-3; %Random value when i have tried to plot this i found that
% with increasing value  of R1 the temperature rise curve fall down.
% Mass of LiF in Kg.
m_liF=2.06e-5;   
% Specific heat capacity of LiF.
c_liF=1562 % J/(kg K)
% Change in temperature depends on the heat supply by the following formula.  
% Q=m*c*deltaT; 
% And power dissipated is defined as heat supply per unit time; 
% i.e P=Q/t.
% Now rise in temperature as a function of time for a cylinder.
T_LiF=T0*exp(I^2*R1_liF*t)  +  I^2*R0_liF*t
scatter(t,T_LiF,'.red',"LineWidth",2)
xlabel('Timerange in Seconds');
ylabel('Rise in temperature in Kelvin')
% title('Plot for rise in temperature as a function of time for lithium fluoride');
% In case of LiF, rise in temperature increases directly with decreasing
% conductivity.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now repeating the above calculation for Si
% Resistivity in ohm.m for silicon.
roh_Si=6.40e2 
% Mass on silicon in Kg.
m_Si=1.83e-5
% Specific heat capacity of silicon in J/(kg K)
c_Si=710 ;
% Resistance of silicon in ohm.
R0_Si=roh_Si* z/A;
R1_Si=2;
% Change in temperature depends on the heat supply by the following formula.  
% Q=m*c*deltaT; 
% And power dissipated is defined as heat supply per unit time; 
% i.e P=Q/t.
% Now rise in temperature.
T_Si=T0*exp(I^2*R1_Si*t)  +  I^2*R0_Si*t
scatter(t,T_Si,'.blue',"LineWidth",2)
xlabel('Timerange in Seconds');
ylabel('Rise in temperature in Kelvin')
%title('Plot for rise in temperature as a function of time for silicon');
% In case of pure semiconductor like silicon, rise in temperature decrease exponentially with
% decreasing  conductivity. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now repeating the above calculation for copper.
roh_Cu=1.68e-8      % Resistivity in ohm.m for copper.
m_Cu=  7.03e-5      % kg 
c_Cu=385            % J/(kg K)
R0_Cu=roh_Cu* z/A;
R1_Cu=R0_Cu*15e-3;

% Q=m*c*deltaT; % m= mass and c=specific heat capacity.
% Now rise in temperature.
T_Cu=T0*exp(I^2*R1_Cu*t)  +  I^2*R0_Cu*t
scatter(t,T_Cu,'.black',"LineWidth",2)
xlabel('Timerange in Seconds');
ylabel('Rise in temperature in Kelvin')
%title('Plot for rise in temperature as a function of time for copper');
% In case of metal rise in temperature increases linearly with increasing
% conductivity. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%