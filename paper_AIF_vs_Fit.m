addpath(genpath('./Source'));

% This script generates the concentration curves for the renal cortex and
% medulla using the 2-compartment model and "trueAIF".

%% Load Global Parameters
paper_loadParameters;
% baseFolderName
% tmax
% hct
% GFR

%% Load Base ROI
load([baseFolderName '/' 'baseROI.mat'])
Vvox = baseROI.Vvox;
ABase = baseROI.ABase;
RCBase = baseROI.RCBase;
RMBase = baseROI.RMBase;
RCSBase = baseROI.RCSBase;
LCBase = baseROI.LCBase;
LMBase = baseROI.LMBase;
LCSBase = baseROI.LCSBase;

%% Model Parameters 
RCVoxCnt = size(RCBase,1);
LCVoxCnt = size(LCBase,1);
k21Vec_R = GFR/(RCVoxCnt*Vvox); % Ktrans
% k21Vec_L = GFR/(size(LCBase,1)*Vvox); % Ktrans
k21Vec_L = k21Vec_R; % Left is half the size of right (50/50 looks bad)
k12 = 1.50; 
vb = 0.35;
vb_med = 0.20;
d = 1.8;
tau = 0;


%% Generate AIF
AIFScalingFactor = 1;
x = [1.459;2.500;0.30046;0.695;0.1200;0.212;4.650;0.0745;13.078;0.953;20/60];
% x = [1.959;2.500;0.30046;0.695;0.1000;0.212;4.650;0.0745;13.078;0.953;20/60]; % Sharper Peak
t_AIF = (0:1:tmax)';
trueAIF = AIFScalingFactor*ParametricAIF(x,t_AIF').';



% load('htrAIF.mat');
% t_AIF = (0:1:tmax)';
% trueAIF = interp1(htrAIF.t,htrAIF.aorta,t_AIF+80,'linear','extrap')/(1-hct);



figure
plot(trueAIF)

xdata = [trueAIF,t_AIF];

% fprintf('Saving true AIF ...\n');
% save('trueAIF.mat', 'trueAIF','-v7.3');

%% Generate Renal Curves
T = length(trueAIF);
artComp = zeros(T,length(GFR));
tubComp = zeros(T,length(GFR));
for kk=1:length(GFR)
    k21 = k21Vec_R(kk);
    xFit = [k21 k12 vb d tau];
    [Ckidney,Cart,Ctub] = ThreeCompartment(xFit,xdata);
    artComp(:,kk) = Cart(:);
    tubComp(:,kk) = Ctub(:);
end

Ckidney = artComp + tubComp;

figure
plot(Ckidney)

figure
plot(tubComp)


%% Save targetPSC
truePSC.Ckidney = Ckidney;
truePSC.trueAIF = trueAIF;
    
% fprintf('Saving true percent signal change curves ...\n');
% save('truePSC.mat', 'truePSC','-v7.3');

%% Analyze with 3 compartment model

AIFScalingFactor = 1;
x = [1.459;2.500;0.30046;0.695;0.1200;0.212;4.650;0.0745;13.078;0.953;20/60];
% x = [1.959;2.500;0.30046;0.695;0.1000;0.212;4.650;0.0745;13.078;0.953;20/60]; % Sharper Peak
x = [1.459;2.500;0.30046;0.695;0.1800;0.212;4.650;0.0745;13.078;0.953;20/60];
t_AIF = (0:1:tmax)';
newAIF = AIFScalingFactor*ParametricAIF(x,t_AIF').';
t_newAIF = 0:T-1;

% newAIF = fitInfo.Cp_aortaHTR;
% t_newAIF = fitInfo.t_aortaHTR;

% newAIF = fitInfo.Cp_aortaVS;
% t_newAIF = fitInfo.t_aortaVS;


G = length(GFR);
kTransVec = zeros(G,1);
for gg = 1:G
    [xFit,resnorm] = FitThreeCompartmentAsymmetric(Ckidney(:,gg),0:T-1,newAIF,t_newAIF);
    kTransVec(gg) = xFit(1);
end

GFR_Est = kTransVec.*(RCVoxCnt*Vvox);

figure
subplot(1,2,1)
plot(t_newAIF,newAIF,'b-')
hold on
plot(t_AIF,trueAIF,'r--')
hold off
subplot(1,2,2)
plot(GFR,GFR_Est,'b-')
hold on
plot(GFR,GFR,'r--')
hold off