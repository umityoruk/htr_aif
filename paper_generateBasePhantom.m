addpath(genpath('./Source'));

%% Load Global Parameters
paper_loadParameters;

% tmax

%% Load mru dataset
[mruDataset, phantomMask] = prepDataset4;

fprintf('Saving phantom masks ...\n');
save('phantomMask.mat', 'phantomMask','-v7.3');


mruIm = single(mruDataset.mruIm);

baselineTime = 6;

% Fix the initial time points
mruIm(:,:,:,1:baselineTime) = repmat(mruIm(:,:,:,baselineTime),[1 1 1 baselineTime]);


%% Interpolate mru dataset to create the base phantom
tmin = 30; % Baseline was too long (~50s)
t_mruIm = (0:size(mruIm,4)-1)*mruDataset.TRes;
t_new = tmin:tmin+tmax;
% t_new = 0:t_mruIm(end);


createBasePhantom(mruIm,t_mruIm,t_new) % Saves as "basePhantom.mat"


%% Generate base roi vectors
load('basePhantom.mat')
T = size(basePhantom,4);
baseROI.ABase = zeros(sum(sum(sum(phantomMask.AMask))),T,class(basePhantom));
baseROI.RCBase = zeros(sum(sum(sum(phantomMask.RCMask))),T,class(basePhantom));
baseROI.RMBase = zeros(sum(sum(sum(phantomMask.RMMask))),T,class(basePhantom));
baseROI.RCSBase = zeros(sum(sum(sum(phantomMask.RCSMask))),T,class(basePhantom));
baseROI.LCBase = zeros(sum(sum(sum(phantomMask.LCMask))),T,class(basePhantom));
baseROI.LMBase = zeros(sum(sum(sum(phantomMask.LMMask))),T,class(basePhantom));
baseROI.LCSBase = zeros(sum(sum(sum(phantomMask.LCSMask))),T,class(basePhantom));

for tt=1:T
    currentIm = basePhantom(:,:,:,tt);
    baseROI.ABase(:,tt) = currentIm(phantomMask.AMask);
    baseROI.RCBase(:,tt) = currentIm(phantomMask.RCMask);
    baseROI.RMBase(:,tt) = currentIm(phantomMask.RMMask);
    baseROI.RCSBase(:,tt) = currentIm(phantomMask.RCSMask);
    baseROI.LCBase(:,tt) = currentIm(phantomMask.LCMask);
    baseROI.LMBase(:,tt) = currentIm(phantomMask.LMMask);
    baseROI.LCSBase(:,tt) = currentIm(phantomMask.LCSMask);
end
baseROI.Vvox = mruDataset.Vvox;

fprintf('Saving base ROIs ...\n');
save('baseROI.mat', 'baseROI','-v7.3');