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
AIFScalingFactor = 0.5; %0.5
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

fprintf('Saving true AIF ...\n');
save('trueAIF.mat', 'trueAIF','-v7.3');

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
    
fprintf('Saving true percent signal change curves ...\n');
save('truePSC.mat', 'truePSC','-v7.3');


%%
% Set ktrans to GFR/Vcx. Normal GFR is 125/2 for a single kidney.
% Sweep GFR from 10 to 80 [ml/min]

%% Generate enhanced ROIs
T = length(trueAIF);
ABaseline = mean(ABase(:,1));
ATarget = ((1-hct)*trueAIF*ABaseline + ABaseline);


RCBaseline = mean(RCBase(:,1));
RMBaseline = mean(RMBase(:,1));
LCBaseline = mean(LCBase(:,1));
LMBaseline = mean(LMBase(:,1));
RCTarget = zeros(T,length(GFR),class(ATarget));
RMTarget = zeros(T,length(GFR),class(ATarget));
LCTarget = zeros(T,length(GFR),class(ATarget));
LMTarget = zeros(T,length(GFR),class(ATarget));

CMIntersect = 80; % seconds
for kk=1:length(GFR)
    k21_R = k21Vec_R(kk);
    k21_L = k21Vec_L(kk);
    xFit_R = [k21_R k12 vb d tau];
    xFit_L = [k21_L k12 vb d tau];
    [Ckidney_R,Cart_R,Ctub_R] = ThreeCompartment(xFit_R,xdata);
    [Ckidney_L,Cart_L,Ctub_L] = ThreeCompartment(xFit_L,xdata);
    Cmed_R = vb_med*Cart_R + (1-vb_med)*Ctub_R;
    Cmed_L = vb_med*Cart_L + (1-vb_med)*Ctub_L;
    RCTarget(:,kk) = Ckidney_R*RCBaseline + RCBaseline;
    RMTarget(:,kk) = (RCTarget(:,kk)-RMBaseline).*Cmed_R/Cmed_R(CMIntersect)+RMBaseline;
    LCTarget(:,kk) = Ckidney_L*LCBaseline + LCBaseline;
    LMTarget(:,kk) = (LCTarget(:,kk)-LMBaseline).*Cmed_L/Cmed_L(CMIntersect)+LMBaseline;
    
%     plot(LCTarget(:,kk),'k-')
%     hold on
%     plot(LMTarget(:,kk),'b-')
%     plot(ATarget,'r-')
%     hold off
%     pause
end


%% Save targetROI

targetROI.ATarget = ATarget;
targetROI.RCTarget = RCTarget;
targetROI.RMTarget = RMTarget;
targetROI.LCTarget = LCTarget;
targetROI.LMTarget = LMTarget;
targetROI.GFR = GFR;
targetROI.Vvox = Vvox;
targetROI.RCVoxCnt = RCVoxCnt;
targetROI.LCVoxCnt = LCVoxCnt;

fprintf('Saving target ROIs ...\n');
save('targetROI.mat', 'targetROI','-v7.3');

