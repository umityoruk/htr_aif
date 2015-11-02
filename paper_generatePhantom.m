load('trueAIF.mat');
Ktrans = 0.28;
createPhantomWithKtrans(trueAIF, t_AIF, Ktrans);


%% Load offline recon (dicoms)

% Link or copy the outphase folder to OfflineRecon in the filesystem
mruImFilename = ['./OfflineRecon/'];

loadData=true;
if (loadData)
    [mruIm,TRes,TR,TE,FA,Vvox,info]=loadDicomFolder(fullfile(pwd, mruImFilename));
    mruIm=double(abs(mruIm));
end

simple4DViewer(mruIm);

t_mruIm = (0:size(mruIm,4)-1)*TRes*1e-3;



% Load sampling info (copy it from one of the mru cases)
load('samplingInfo.mat');
t_phantom = samplingInfo.t;


%% Create base phantom
HTR_createBasePhantom(mruIm,t_mruIm,t_phantom);



%% Now export aortaMask
aortaMaskFilename = ['./AortaDicom/'];

loadData=true;
if (loadData)
    [aortaMaskRaw]=loadDicomFolder(fullfile(pwd, aortaMaskFilename));
    aortaMaskRaw=double(abs(aortaMaskRaw));
end

threshold = 1024
aortaMask = double(aortaMaskRaw > threshold);

% Need to adjust the number of slices because stupid Osirix changes the
% number of slices during reslicing!!
aortaMask = aortaMask(:,:,1:3:end);
aortaMask = aortaMask(:,:,1:end-1);

% Check if you cut the right slices out.
simple4DViewer(mruIm,aortaMask);

save('aortaMask.mat','aortaMask');




