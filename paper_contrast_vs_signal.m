r1 = 4.5; % (1/s)(1/mM)

T10_blood = 1.4; % s
T10_kidney = 1.2; % s

T10 = T10_blood;

R10 = 1/T10;

Cgd = 0:0.01:1; % mM

R1_t = R10 + r1*Cgd;

% figure
% plot(Cgd,1./R1_t)
% xlabel('Concentration [mM]')
% ylabel('T1 [s]')


FA = 12; % degrees
TR = 4e-3; % s
theta = FA/180*pi;

baseVal = 5959;
S0 = baseVal/(1-exp(-R10*TR))*sin(theta)./(1-exp(-R10*TR)*cos(theta));
% S0 = 10000;
S_t = S0*(1-exp(-R1_t*TR))*sin(theta)./(1-exp(-R1_t*TR)*cos(theta));


figure
plot(Cgd,S_t)
xlabel('Concentration [mM]')
ylabel('Signal')


dS_t = (S_t - S_t(1))/S_t(1);
figure
plot(dS_t,Cgd)
ylabel('Concentration [mM]')
xlabel('Signal Change')

stepCgd = (Cgd(2)-Cgd(1))/(dS_t(2)-dS_t(1));
linCgd = dS_t*stepCgd;
figure
plot(dS_t,Cgd,'b-')
hold on
plot(dS_t,linCgd,'r-')
ylabel('Concentration [mM]')
xlabel('Signal Change')

figure
plot(dS_t,(linCgd-Cgd)./Cgd*100,'r-')
ylabel('Error (%)')
xlabel('Signal Change')

% figure
% plot(Cgd,dS_t)
% xlabel('Concentration [mM]')
% ylabel('Signal Change')