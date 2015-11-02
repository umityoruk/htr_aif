% Load the dataset for analysis
[mruDataset,mruMaskset] = prepDataset_mrucase2;

mruIm = mruDataset.mruIm;
TRes = mruDataset.TRes;
TR = mruDataset.TR;
FA = mruDataset.FA;
AMask = mruMaskset.AMask;
RCMask = mruMaskset.RCMask;
LCMask = mruMaskset.LCMask;
RMMask = mruMaskset.RMMask;
LMMask = mruMaskset.LMMask;
RCSMask = mruMaskset.RCSMask;
LCSMask = mruMaskset.LCSMask;
RPMask = mruMaskset.RPMask;
LPMask = mruMaskset.LPMask;
Vvox = mruDataset.Vvox;

sliceSelect = 13;

%%
RCMask = squeeze(RCMask(:,:,sliceSelect));
LCMask = squeeze(LCMask(:,:,sliceSelect));
RMMask = squeeze(RMMask(:,:,sliceSelect));
LMMask = squeeze(LMMask(:,:,sliceSelect));
RCSMask = squeeze(RCSMask(:,:,sliceSelect));
LCSMask = squeeze(LCSMask(:,:,sliceSelect));
%%
simple4DViewer(mruIm)

%%
mip_ss = 13:15;

mipIm = max(mruIm(:,:,mip_ss,:),[],3);

% simple4DViewer(mipIm)

figure
imagesc(squeeze(mipIm(:,:,:,11)))
colormap('gray')


%%

bg = squeeze(mruIm(:,:,sliceSelect,end));

figure
imagesc(bg)
colormap('gray')
hold on
emptyL = zeros(size(bg));

colorCx = [0.0 0.0 1.0];
colorMd = [0.0 1.0 0.0];
colorCs = [0.7 0.0 0.0];

RC = cat(3,colorCx(1)*RCMask,colorCx(2)*RCMask,colorCx(3)*RCMask);
h=imagesc(RC);
set(h,'AlphaData',RCMask);

RM = cat(3,colorMd(1)*RMMask,colorMd(2)*RMMask,colorMd(3)*RMMask);
h=imagesc(RM);
set(h,'AlphaData',RMMask);

RCS = cat(3,colorCs(1)*RCSMask,colorCs(2)*RCSMask,colorCs(3)*RCSMask);
h=imagesc(RCS);
set(h,'AlphaData',RCSMask);


LC = cat(3,colorCx(1)*LCMask,colorCx(2)*LCMask,colorCx(3)*LCMask);
h=imagesc(LC);
set(h,'AlphaData',LCMask);

LM = cat(3,colorMd(1)*LMMask,colorMd(2)*LMMask,colorMd(3)*LMMask);
h=imagesc(LM);
set(h,'AlphaData',LMMask);

LCS = cat(3,colorCs(1)*LCSMask,colorCs(2)*LCSMask,colorCs(3)*LCSMask);
h=imagesc(LCS);
set(h,'AlphaData',LCSMask);

axis('image')