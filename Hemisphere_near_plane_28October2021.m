
% In this MATLAB script, I plot J(E) and FN curves for 
% work functions of 5, 2, 1, and 0.1 eV for a 50 micron radius hemisphere and a 100 um gap,
% and then adjust the voltage from 1 V to 10 kV,
% and then plot J(E) and FN curves for a 50 micron radius hemisphere and a 100 um gap. 
% Plot all four curves on the same axis, and then fit a straight line to the FN plots to demonstrate that I have recovered the work function in the manner anticipated.
% Then, in the following phase, I added the noise to the data; in practise, I have a noise floor of 10 nA.
% Then, in the following step, I converted the J(E) plots into I(V) plots using the geometry, for example, for (a hemisphere on a plane), and then I added noise of 10 nA to the plots.
% Then define a threshold voltage as the voltage at which the current increases above the level of the background noise.
% Fowler-Nordheim- type- calculations-equation
% Hemisphere near plane
% 28 October 2021
% Witten by Najwa 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/Incorporate realistic values for field enhancement.
clc
clear

a= 1.541434e-6;         % correction factor in A ev V^-2
b=1 ;                   % b is a constant arising from the derivation of the Fowler-Nordheim equation.
phi= 1;                 % work Function,in eV.
vf=1;                   % constant
V=0.005:0.005:0.5;
% a's are the separation between the top of the sphere.
% r's are the radius of the sphere.
% Gamma is field enhancement factor.
% F is the electric field.

a1=1e-3; a2=1e-3;a3=1e-3; r1=1e-3;r2=1e-4;r3=1e-5;
gamma1=0.9*((r1+a1)/r1)
gamma2=0.9*((r2+a2)/r2)
gamma3=0.9*((r3+a3)/r3)
F1=V/a1; F2=V/a2; F3=V/a3;
% Use Fowler-Nordheim equation to calculate current density in Amp/m^2
% F_max is the maximum electric field after applying field enhancement factor.
F1_max=  gamma1*F1; F2_max=gamma2*F2; F3_max=gamma3*F3;
X1=exp(-vf*b*phi.^(3/2)./F1_max);
X2=exp(-vf*b*phi.^(3/2)./F2_max);
X3=exp(-vf*b*phi.^(3/2)./F3_max);
% Calculate the field emitted current density J1,J2,J3 
% in Amps per square metre.
J1=a*phi^(-1)*F1_max.^(2).*X1 ;% current density in Amp/m^2
J2=a*phi^(-1)*F2_max.^(2).*X2;
J3=a*phi^(-1)*F3_max.^(2).*X3;
% Current density (J) against electric field (F).
figure()
plot(F1,J1,'red',"LineWidth",2);
hold on
plot(F2,J2,'blue',"LineWidth",2);
plot(F3,J3,'black',"LineWidth",2);
hold off
xlabel('F(V/m)');
ylabel('J(amp/m^2)');
legend('F1,J1 with radius of 1x10^-3 m','F2,J2 with radius of 1x10^-4 m','F3,J3 with radius of 1x10^-5 m');
title('Current density vs electric field');

% For each of current J vs field F , I plotted ln(J/F^2) vs 1/F,and I got a
% straight line.
% Each plot of J(F) is accompanied by a plot of plot of ln(J/F^2) versus 1/F,
% fitting 
f1= 1./F1;
f2=1./F2;
f3=1./F3;                % m/V
j1=log(J1./F1_max.^(2)); % Amp.m/V
j2=log(J2./F2_max.^(2));
j3=log(J3./F3_max.^(2));
plot(f1,j1','O');
P=polyfit(f1,j1,1);
yfit=P(1)*f1+P(2);
hold on
plot(f1,yfit,'red','linewidth',2);
text(0.02,-13.48,['j1=' num2str(P(1)) 'f1' num2str(P(2))],'color','red');
plot(f2,j2,'o');
P=polyfit(f2,j2,1);
yfit=P(1)*f2+P(2);
plot(f2,yfit,'blue','linewidth',2);
text(0.02,-13.47,['j2=' num2str(P(1)) 'f2' num2str(P(2))],'color','blue');
plot(f3,j3,'o');
P=polyfit(f3,j3,1);
yfit=P(1)*f3+P(2);
plot(f3,yfit,'black','linewidth',2);
text(0.02,-13.46,['j3=' num2str(P(1)) 'f3' num2str(P(2))],'color','black');
xlabel('1/F');
ylabel('Log J/F^2 ');
title('Plot of ln(J/F^2 )versus 1/F');
legend('f1,j1','gradient 1','f2,j2','gradient 2','f3,j3','gradient 3')
hold off
% Gradient(slope) = vf*b* (workfunction)^1.5/gamma(enhancement factor)
% verified for each curves. e.g, 0.56=1*1*(1)^1.5/1.8
% i. e. 0.56=0.56 ###
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2/noise calculation
% Y Noisy= YOriginal+ amplitude*rand(1,length(YOriginal));
% plot(t,[x,y])
plot(F1,J1,'.','color', 'red');
hold on
noise_floor= 1e-8; % in A/m2
noisy_J1=J1+noise_floor*(rand(1, length(J1)));
plot(F1,noisy_J1,'-',"Color",'red');
xlabel('F (V/m)');
ylabel('J(amps/m^2)');
%title('Curent density vs electric field with noise');
hold off

% 3/Convert J(E) plot to I(V),then add noise of 10 nA,
% then define a threshold voltage as the point 
% at which the current rises above the noise.
 
I1=3*pi*r1.^(2)*J1;
I2=3*pi*r2.^(2)*J2;
I3=3*pi*r3.^(2)*J3;
plot(V,I1, 'red','linewidth',2);
hold on
plot(V,I2,'blue','linewidth',2);
hold on
plot(V,I3,'black','linewidth',2);
xlim([0,0.1]);
ylim([0, 1.5e-7]);
xlabel('voltage V(V)');
ylabel('current I(A)');
legend('V,I1','V,I2','V,I3')
legend('V,I1 with radius of 1x10^-3 m','V,I2 with radius of 2x10^-3m','V,I3 with radius of 3x10^-3m');
%title('Current Vs Voltage graph');
hold off

% Adding noise at noise floor 10 nA.

noise_floor= 1e-8 ;
noisy_I1=I1+noise_floor*(rand(1, length(I1)));
plot(V,I1,'.','Color','red');
hold on
plot(V,noisy_I1,'-','color','red');
noisy_I2=I2+noise_floor*(rand(1, length(I2)));
plot(V,I2,'.','color','blue');
plot(V,noisy_I2,'-','color','blue');
noisy_I3=I3+noise_floor*(rand(1, length(I3)));
plot(V,I3,'.','color','black');
plot(V,noisy_I3,'-','color','black');
x1=0; x2=0.5; y1=1e-8;y2=1e-8;
plot([x1,x2],[y1,y2]);
xlim([0,0.1]);
ylim([0, 1.5e-7]);
idxmin = find(I1== 1e-8);
idxmax1=find(V==0.015);
idxmax2=find(V==0.025);
idxmax3=find(V==0.03);

plot(V,I1,'+r','LineWidth',4,'MarkerIndices',[idxmax1, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I2,'+b','LineWidth',4,'MarkerIndices',[idxmax2, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I3,'+k','LineWidth',4,'MarkerIndices',[idxmax3, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
ylabel('Current I(A)');
xlabel('Voltage V(V)');
legend('VI1,with radius 1x10^-3m','VI1 with noise','VI2 with radius 2x10^-3m','V,I2 with noise','VI3 with radius 3x10^-3 m','V,I3 with noise');
%title('Current versus voltage graph with noise');
grid on
hold off


% Threshold voltage In V

 Vth1=0.015 ; Vth2=0.025; Vth3=0.03; 
Vth=[Vth1,Vth2,Vth3];
r=[r1,r2,r3];
scatter(r,Vth,'o')
hold on
plot(r,Vth,'red');
xlabel('Radius r(m)');
ylabel('Threshold (V)');
%title('Thresold voltage vs radius of hemisphere at constant distance of 100micron');
%hold off

clc
clear
% Fowler-Nordheim- type- calculations-equation
a= 1.541434e-6;        % correction factor in A ev V^-2
b=1 ;                  % b is a constant arising from the derivation of the Fowler-Nordheim equation.
phi= 1;                % work Function,in eV.
vf=1;                  % constant
V=0.005:0.005:0.5;

a1=1e-3; a2=1e-4;a3=1e-5; r1=1e-5;r2=1e-5;r3=1e-5;
gamma1=0.9*((r1+a1)/r1)
gamma2=0.9*((r2+a2)/r2)
gamma3=0.9*((r3+a3)/r3)
F1=V/a1; F2=V/a2; F3=V/a3;

F1_max=  gamma1*F1; F2_max=gamma2*F2; F3_max=gamma3*F3;
X1=exp(-vf*b*phi.^(3/2)./F1_max);
X2=exp(-vf*b*phi.^(3/2)./F2_max);
X3=exp(-vf*b*phi.^(3/2)./F3_max);

% Calculate the field emitted current density J1,J2,J3 
% in Amps per square metre.
J1=a*phi^(-1)*F1_max.^(2).*X1 ;
J2=a*phi^(-1)*F2_max.^(2).*X2;
J3=a*phi^(-1)*F3_max.^(2).*X3;
% Current density (J) against electric field (F).
figure()
plot(F1,J1,'red','linewidth',2);
hold on
plot(F2,J2,'blue','linewidth',2);
plot(F3,J3,'black','linewidth',2);
hold off
xlabel('F(V/m)');
ylabel('J(amp/m^2)');
xlim([0,1e4]);
ylim([0,3000]);
legend('F1,J1','F2,J2','F3,J3')
legend('F1,J1 with separation of 1x10^-3 m','F2,J2 with separation of 2x10^-4 m','F3,J3 with separation of 3x10^-5 m');
%title('Current density vs electric field');

% For each of current J vs field F , I plotted ln(J/F^2) vs 1/F,and I got a
% straight line.
% Each plot of J(F) is accompanied by a plot of plot of ln(J/F^2) versus 1/F,
% fitting 
f1= 1./F1;
f2=1./F2;
f3=1./F3;                % m/V
j1=log(J1./F1_max.^(2)); % Amp.m/V
j2=log(J2./F2_max.^(2));
j3=log(J3./F3_max.^(2));
plot(f1,j1','O');
hold on
P=polyfit(f1,j1,1);
yfit=P(1)*f1+P(2);
plot(f1,yfit,'red','linewidth',2);
text(0.025,-13.3842,['j1=' num2str(P(1)) 'f1' num2str(P(2))],'color','red');
plot(f2,j2,'o');
P=polyfit(f2,j2,1);
yfit=P(1)*f2+P(2);
plot(f2,yfit,'blue','linewidth',2);
text(0.025,-13.3844,['j2=' num2str(P(1)) 'f2' num2str(P(2))],'color','blue');
plot(f3,j3,'o');
P=polyfit(f3,j3,1);
yfit=P(1)*f3+P(2);
plot(f3,yfit,'black','linewidth',2);
text(0.025,-13.3846,['j3=' num2str(P(1)) 'f3' num2str(P(2))],'color','black');

xlabel('1/F');
ylabel('Log J/F^2 ');
legend('f1,j1','gradient 1','f2,j2','gradient 2','f3,j3','gradient 3')
%title('Plot of ln(J/F^2) versus 1/F');
hold off
% Gradient(slope) = vf*b* (workfunction)^1.5/gamma(enhancement factor)
% verified for each curves. e.g, 0.56=1*1*(1)^1.5/1.8
% i. e. 0.56=0.56 ###

% 2/noise calculation
% Y Noisy= YOriginal+ amplitude*rand(1,length(YOriginal));
% plot(t,[x,y])
plot(F1,J1,'red');
hold on
noise_floor= 1e-8; % in A/m^2
noisy_J1=J1+noise_floor*(rand(1, length(J1)));
plot(F1,noisy_J1,'-',"Color",'red');
xlabel('F (V/m)');
ylabel('J(amps/m^2)');
title('Curent density vs electric field with noise');
hold off

% 3/Convert J(E) plot to I(V),then add noise of 10 nA,
% then define a threshold voltage as the point 
% at which the current rises above the noise.
 
I1=3*pi*r1.^(2)*J1;
I2=3*pi*r2.^(2)*J2;
I3=3*pi*r3.^(2)*J3;
plot(V,I1, 'red','linewidth',2);
hold on
plot(V,I2,'blue','linewidth',2);
plot(V,I3,'black','linewidth',2);
xlim([0,0.1]);
ylim([0,1.5e-7]);
xlabel('Voltage V(V)');
ylabel('Current I(A)');
legend('V,I1','V,I2','V,I3')
title('Current vs voltage graph');
hold off

% Adding noise at noise floor 10 nA.

noise_floor= 1e-8 ;
noisy_I1=I1+noise_floor*(rand(1, length(I1)));
plot(V,I1,'.','Color','red');
hold on
plot(V,noisy_I1,'-','color','red');
noisy_I2=I2+noise_floor*(rand(1, length(I2)));
plot(V,I2,'.','color','blue');
plot(V,noisy_I2,'-','color','blue');
noisy_I3=I3+noise_floor*(rand(1, length(I3)));
plot(V,I3,'.','color','black');
plot(V,noisy_I3,'-','color','black');
 x1=0; x2=0.5; y1=1e-8;y2=1e-8;
 plot([x1,x2],[y1,y2]);
 xlim([0,0.1]);
ylim([0, 1.5e-7]);
idxmin = find(I1== 1e-8);
idxmax1=find(V==0.02);
idxmax2=find(V==0.025);
idxmax3=find(V==0.01);

plot(V,I1,'+r','LineWidth',4,'MarkerIndices',[idxmax1, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I2,'+b','LineWidth',4,'MarkerIndices',[idxmax2, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I3,'+k','LineWidth',4,'MarkerIndices',[idxmax3, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
ylabel('Current I(A)');
xlabel('Voltage V(V)');
legend('V,I1','noise1','V,I2','noise2','V,I3','noise3')
title('Current versus voltage graph with noise');
grid on
hold off


% Threshold voltage In V

Vth1=0.02 ; Vth2=0.025; Vth3=0.1; 
Vth=[Vth1,Vth2,Vth3];
a=[a1,a2,a3];
scatter(a,Vth,'o');
hold on
plot(a,Vth,'red');
xlabel('Distance a(m)');
ylabel('Threshold voltage(V)');
title('Threshold voltage vs anode distance at constant hemisphere radius of 10 micron');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
% Fowler-Nordheim- type- calculations-equation
a= 1.541434e-6;             % correction factor in A ev V^-2
b=1 ;                       % b is a constant arising from the derivation of the Fowler-Nordheim equation.
%  Phi is work Function,in eV.
vf=1;                       % constant
V=0.005:0.005:0.5;

a1=1e-3; a2=1e-3;a3=1e-3; r1=1e-5;r2=1e-5;r3=1e-5;
gamma1=0.9*((r1+a1)/r1)
gamma2=0.9*((r2+a2)/r2);gamma3=0.9*((r3+a3)/r3);
F1=V/a1; F2=V/a2; F3=V/a3;

% Phi is a work function in V/m % random value
phi1=5; phi2=2; phi3=1;
F1_max=  gamma1*F1; F2_max=gamma2*F2; F3_max=gamma3*F3;
X1=exp(-vf*b*phi1.^(3/2)./F1_max);
X2=exp(-vf*b*phi2.^(3/2)./F2_max);
X3=exp(-vf*b*phi3.^(3/2)./F3_max);

% Calculate the field emitted current density J1,J2,J3 
% in Amps per square metre.

J1=a*phi1^(-1)*F1_max.^(2).*X1; 
J2=a*phi2^(-1)*F2_max.^(2).*X2;
J3=a*phi3^(-1)*F3_max.^(2).*X3;
% Current density (J) against electric field (F).
figure()
plot(F1,J1,'red','linewidth',2);
hold on
plot(F2,J2,'blue','linewidth',2);
plot(F3,J3,'black','linewidth',2);
hold off
xlabel('F(V/m)');
ylabel('J(amp/m^2');
legend('F1,J1','F2,J2','F3,J3')
title('Current density vs electric field');

% For each of current J vs field F , I plotted ln(J/F^2) vs 1/F,and I got a
% straight line.
% Each plot of J(F) is accompanied by a plot of plot of ln(J/F^2) versus 1/F,
% fitting 
f1= 1./F1;
f2=1./F2;
f3=1./F3;                % m/V
j1=log(J1./F1_max.^(2)); % Amp.m/V
j2=log(J2./F2_max.^(2));
j3=log(J3./F3_max.^(2));
plot(f1,j1','O');
hold on
P=polyfit(f1,j1,1);
yfit=P(1)*f1+P(2);
plot(f1,yfit,'red','linewidth',2);
text(0.05,-5.8,['j1=' num2str(P(1)) 'f1' num2str(P(2))],'color','red');
plot(f2,j2,'o');
P=polyfit(f2,j2,1);
yfit=P(1)*f2+P(2);
plot(f2,yfit,'blue','linewidth',2);
text(0.05,-4.4,['j2=' num2str(P(1)) 'f2' num2str(P(2))],'color','blue');
plot(f3,j3,'o');
P=polyfit(f3,j3,1);
yfit=P(1)*f3+P(2);
plot(f3,yfit,'black','linewidth',2);
text(0.05,-5,['j3=' num2str(P(1)) 'f3' num2str(P(2))],'color','black');
xlabel('1/F');
ylabel('LogJ/F^2 ');
title('Plot of ln(J/F^2) versus 1/F');
hold off
% Gradient(slope) = vf*b* (workfunction)^1.5/gamma(enhancement factor)
% verified for each curves. e.g, 0.12=1*1*(1)^1.5/90.90
% i. e. 0.12=0.11 ###

% 2/noise calculation
% Y Noisy= YOriginal+ amplitude*rand(1,length(YOriginal));
% plot(t,[x,y])
plot(F1,J1,'red');
hold on
noise_floor= 1e-8;       % in A/m^2
noisy_J1=J1+noise_floor*(rand(1, length(J1)));
plot(F1,noisy_J1,'-',"Color",'red');
xlabel('F (V/m)');
ylabel('J(amps/m^2)');
title('Curent density vs electric field with noise');
hold off

% 3/Convert J(E) plot to I(V),then add noise of 10 nA,
% then define a threshold voltage as the point 
% at which the current rises above the noise.
 
I1=3*pi*r1.^(2)*J1;
I2=3*pi*r2.^(2)*J2;
I3=3*pi*r3.^(2)*J3;
plot(V,I1, 'red','linewidth',2);
hold on
plot(V,I2,'blue','linewidth',2);
plot(V,I3,'black','linewidth',2);
xlim([0,0.1]);
ylim([0, 1.5e-7]);
xlabel('Voltage V(V)');
ylabel('Current I(A)');
legend('V,I1 with phi1=5 m','V,I2 with phi2=2m','V,I3 with phi3=1m');
title('Current vs voltage graph')
hold off

% Adding noise at noise floor 10 nA. 

noise_floor= 1e-8; 
noisy_I1=I1+noise_floor*(rand(1, length(I1)));
plot(V,I1,'.','Color','red');
hold on
plot(V,noisy_I1,'-','color','red');
noisy_I2=I2+noise_floor*(rand(1, length(I2)));
plot(V,I2,'.','color','blue');
plot(V,noisy_I2,'-','color','blue');
noisy_I3=I3+noise_floor*(rand(1, length(I3)));
plot(V,I3,'.','color','black');
plot(V,noisy_I3,'-','color','black');
x1=0; x2=0.5; y1=1e-8;y2=1e-8;
plot([x1,x2],[y1,y2]);
xlim([0,0.1]);
ylim([0, 1.5e-7]);
idxmin = find(I1== 1e-8);
idxmax1=find(V==0.065);
idxmax2=find(V==0.04);
idxmax3=find(V==0.03);

plot(V,I1,'+r','LineWidth',4,'MarkerIndices',[idxmax1, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I2,'+b','LineWidth',4,'MarkerIndices',[idxmax2, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I3,'+k','LineWidth',4,'MarkerIndices',[idxmax3, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
ylabel('Current I(A)');
xlabel('Voltage V(V)');
legend('VI1,with phi1=5m','VI1 with noise','VI2 with phi2=2 m','V,I2 with noise','VI3 with phi3=1 m','V,I3 with noise');
title('Current versus voltage graph with noise')
grid on
hold off

% Threshold voltage In V

Vth1=0.065 ; Vth2=0.04; Vth3=0.03; 
Vth=[Vth1,Vth2,Vth3];
phi=[phi1,phi2,phi3];
scatter(phi,Vth,'o')
hold on
plot(phi,Vth,'red');
xlabel('Workfunction phi(eV)');
ylabel('Vth(V)');
title('Threshold voltage vs work function at constant a=100 micron and r=10 micron');