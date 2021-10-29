% In this MATLAB script to calculate the field emission current density J as function of applied electric field F for metal cathodes of differing geometry.
% i. Parallel plate
% 1. J(F) and F-N plots for parallel plate capacitor of appropriate dimensions
% 2. I(V)  plots of parallel plate capacitor with noise and threshold voltage Vth
% 3. Vth as function of surface barrier
% in this part I do not used field enhancement factor (gamma)!
% Fowler-Nordheim- type- calculations-equation
% 27 October 2021
% Witten by Najwa    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. J(F) and F-N plots for parallel plate capacitor of appropriate dimensions

clc
clear
a= 1.541434e-6;    % correction factor in A ev V^-2
b=1;               % b is a constant arising from the derivation of the Fowler-Nordheim equation
phi= 1;            % work Function,in eV.
vf=1;              % constant
V=0.005:0.005:0.5;        % V is the range of voltages applied to the sample.
% a's are the separation between the top of the sphere.
% r's are the radius of the sphere.
% F is the electric field.
% x is a separation between two parallel plate in m
x1=1e-3;x2=2e-3;x3=3e-3; 
% L is a length and W is a width of the plate in m
L=1e-3; W=1e-3;
% F is a  electric field strength in V/m
F1=V/x1; F2=V/x2; F3=V/x3;
% Use Fowler-Nordheim equation to calculate current density in Amp/m^2.


X1=exp(-vf*b*phi.^(3/2)./F1);
X2=exp(-vf*b*phi.^(3/2)./F2);
X3=exp(-vf*b*phi.^(3/2)./F3);


% Calculate the field emitted current density J1,J2,J3 
% in Amps per square metre.

J1=a*phi^(-1)*F1.^(2).*X1 ;
J2=a*phi^(-1)*F2.^(2).*X2;
J3=a*phi^(-1)*F3.^(2).*X3;
% Plot J(F) for parallel plate.
figure()
plot(F1,J1,'.r'); 
hold on
plot(F2,J2,'.k')
plot(F3,J3,'.b');
xlabel('F(V/m)');
ylabel('J(amp/m^2)');
legend('F1,J1 with separation of 1x10^-3 m','F2,J2 with separation of 2x10^-3 m','F3,J3 with separation of 3x10^-3 m');
title('Current density vs electric field strength',"FontSize",15);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Insert code to  plot the F-N plot for the three emitters.
% Make fowler-Nordheim plot i.e. 'plot of ln(J/F^2) versus 1/F'
% Fit F-N curves to extract gradient and check that it is what we
% expect.
f1= 1./F1;
f2=1./F2;
f3=1./F3;            % m/V
j1=log(J1./F1.^(2)); % Amp.m/V
j2=log(J2./F2.^(2));
j3=log(J3./F3.^(2));
plot(f1,j1','O');
P=polyfit(f1,j1,1);
yfit=P(1)*f1+P(2);
hold on
plot(f1,yfit,'red','linewidth',2);
text(0.25,-13.35,['j1=' num2str(P(1)) 'f1' num2str(P(2))],'color','red');
plot(f2,j2,'o');
P=polyfit(f2,j2,1);
yfit=P(1)*f2+P(2);
plot(f2,yfit,'black','linewidth',2);
text(0.25,-13.50,['j2=' num2str(P(1)) 'f2' num2str(P(2))],'color','black');
plot(f3,j3,'o');
P=polyfit(f3,j3,1);
yfit=P(1)*f3+P(2);
plot(f3,yfit,'blue','linewidth',2);
text(0.25,-13.60,['j3=' num2str(P(1)) 'f3' num2str(P(2))],'color','blue');
xlabel('1/F');
ylabel('LogJ/F2 ');
legend('f1,j1','gradient 1','f2,j2','gradient 2','f3,j3','gradient 3')
title('Plot of ln(J/F^2) versus 1/F',"FontSize",15);
hold off
%% because the value of phi for each curve is same we got same gradient for each curve.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. I(V)  plots of parallel plate capacitor with noise and threshold voltage Vth
% convert current density to current using dimensions of plates
I1=L*W*J1;
I2=L*W*J2;
I3=L*W*J3;
% plotI(V) curves

plot(V,I1, 'red',"LineWidth",2 );
hold on
plot(V,I2,'black','LineWidth',2);
hold on
plot(V,I3,'blue',"LineWidth",2);
xlim([0,0.5]);
ylim([0, 1.5e-7]);
xlabel('Voltage V(V)');
ylabel('Current I(A)')
legend('V,I1 with separation of 1x10^-3 m','V,I2 with separation of 2x1-^-3 m','V,I3 with separation of 3x10^-3 m')
title('Current vs voltage graph',"FontSize",15);
hold off
% Plot log(1/V^2) vs 1/V.
v=1./V;
logon1=log(I1./V.^2);
logon2=log(I2./V.^2);
logon3=log(I3./V.^2);
plot(v,logon1,'.r')
plot(v,logon1','O');
P=polyfit(v,logon1,1);
yfit=P(1)*v+P(2);
hold on
plot(v,yfit,'red','linewidth',2);
text(50,-13.7,['logon1=' num2str(P(1)) 'v' num2str(P(2))],'color','red');
plot(v,logon2,'o');
P=polyfit(v,logon2,1);
yfit=P(1)*v+P(2);
plot(v,yfit,'black','linewidth',2);
text(50,-15.3,['logon2=' num2str(P(1)) 'v' num2str(P(2))],'color','black');
plot(v,logon3,'o');
P=polyfit(v,logon3,1);
yfit=P(1)*v+P(2);
plot(v,yfit,'blue','linewidth',2);
text(50,-16,['logon3=' num2str(P(1)) 'v' num2str(P(2))],'color','blue');
xlabel('1/v');
ylabel('Log (I/V^2) ');
legend('v,logon1','gradient 1','v,logon2','gradient 2','v,lohon3','gradient 3')
title('Plot of ln(I/V^2) versus 1/V',"FontSize",15);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adding noise to I(V) curves.
% Set noise_floor at 10 nA i.e. 1e-8 A

noise_floor= 1e-8;
noisy_I1=I1+noise_floor*(rand(1, length(I1)));
plot(V,I1,'.','Color','red');
hold on
plot(V,noisy_I1,'-','color','red');
noisy_I2=I2+noise_floor*(rand(1, length(I2)));
plot(V,I2,'.','color','black');
plot(V,noisy_I2,'-','color','black');
noisy_I3=I3+noise_floor*(rand(1, length(I3)));
plot(V,I3,'.','color','blue');
plot(V,noisy_I3,'-','color','blue');
x1=0; x2=0.5; y1=1e-8; y2=1e-8;
plot([x1,x2],[y1,y2]);
xlim([0,0.5]);
ylim([0, 1.5e-7]);
idxmin = find(I1== 1e-8);
idxmax1=find(V==0.08);
idxmax2=find(V==0.16);
idxmax3=find(V==0.25);

plot(V,I1,'+r','LineWidth',4,'MarkerIndices',[idxmax1, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I2,'+k','LineWidth',4,'MarkerIndices',[idxmax2, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I3,'+b','LineWidth',4,'MarkerIndices',[idxmax3, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)


ylabel('Current I(A)');
xlabel('Voltage V(V)');
legend('V,I1 with separation of 1x10^-3 m','V,I1 with noise','V,I2 with separation of 2x10^-3 m','V,I2 with noise','V,I3 with separation of 3x10^-3 m','V,I3 with noise');
title('Current versus voltage graph with noise',"FontSize",10);
grid on
hold off
% Threshold voltage in(V)
 
Vth1=0.08 ; Vth2=0.16; Vth3=0.25; 
%=====================================
clc
clear
% 3. threshold against surface barrier
a= 1.541434e-6;      % correction factor in A ev V^-2
b=1 ;                % eV ^(-3/2)V nm^-1 b it is related to the model used in deriving the expressions.
% phi is work Function,in eV.
vf=1;                % constant
V=0.005:0.005:0.5;
% a1,a2 and a3 are radius of hemisphere.
% length and width of the plate in m
x1=1e-3;x2=1e-3;x3=1e-3; 
L=1e-3;W=1e-3;
phi1=5; phi2=2; phi3=1;

F1=V/x1; F2=V/x2; F3=V/x3;
X1=exp(-vf*b*phi1.^(3/2)./F1);
X2=exp(-vf*b*phi2.^(3/2)./F2);
X3=exp(-vf*b*phi3.^(3/2)./F3);
% current density in Amp/m squared.
J1=a*phi1^(-1)*F1.^(2).*X1 ;
J2=a*phi2^(-1)*F2.^(2).*X2;
J3=a*phi3^(-1)*F3.^(2).*X3;
I1=L*W*J1;
I2=L*W*J2;
I3=L*W*J3
plot(V,I1, 'red','linewidth',2);
hold on
plot(V,I2,'black','linewidth',2);
hold on
plot(V,I3,'blue','linewidth',2);
xlim([0,0.5]);
ylim([0, 1.5e-7]);
ylabel('current I(A)');
xlabel('voltage V(V)');
legend('V,I1','V,I2','V,I3');
title('Current versus voltage graph',"FontSize",15);

hold off

% Adding noise at noise floor of 10 nA .

noise_floor= 1e-8;
noisy_I1=I1+noise_floor*(rand(1, length(I1)))
plot(V,I1,'.','Color','red')
hold on
plot(V,noisy_I1,'-','color','red')
hold off
hold on
noisy_I2=I2+noise_floor*(rand(1, length(I2)))
plot(V,I2,'.','color','black')
hold on
plot(V,noisy_I2,'-','color','black')
hold off
hold on
noisy_I3=I3+noise_floor*(rand(1, length(I3)))
plot(V,I3,'.','color','blue')
hold on
plot(V,noisy_I3,'-','color','blue')
x1=0; x2=0.5;y1=1e-8, y2=1e-8;
plot([x1,x2,],[y1,y2])
grid on
xlim([0,0.5]);
ylim([0, 1.5e-7]);
idxmin = find(I1== 1e-8);
idxmax1=find(V==0.19);
idxmax2=find(V==0.13);
idxmax3=find(V==0.08);

plot(V,I1,'+r','LineWidth',4,'MarkerIndices',[idxmax1, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I2,'+k','LineWidth',4,'MarkerIndices',[idxmax2, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
plot(V,I3,'+b','LineWidth',4,'MarkerIndices',[idxmax3, idxmin],...
    'MarkerFaceColor','red',...
    'MarkerSize',15)
xlabel('Volage in V')
ylabel('Current in nAmps')
legend('V,I1 for phi=5 ev','noiseI1','V,I2 for phi=2 eV','noiseI2','V,I3 for phi=1 eV','noiseI3')
hold off
title('Adding noise to the current vs voltage graph',"FontSize",15)
hold off
hold off

% Threshold voltage in(V)
% Vth against surface barrier 
Vth1=0.19 ; Vth2=0.13; Vth3=0.08; 
phi=[phi1,phi2,phi3]
Vth=[Vth1,Vth2,Vth3] 
scatter(phi,Vth,'O')
hold on
plot(phi,Vth,'red');
xlabel('Surface barrier(eV)');
ylabel('Threshold (V)');
title('Threshold voltage against surface barrier',"FontSize",15)
hold off
% vth is proportional to work function .