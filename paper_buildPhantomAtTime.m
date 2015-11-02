function phantom3d = paper_buildPhantomAtTime(basePhantom,phantomMask,targetROI,gfrInd,timeInd)

phantom3d = basePhantom(:,:,:,timeInd);

AMask = phantomMask.AMask;
RCMask = phantomMask.RCMask;
RMMask = phantomMask.RMMask;
LCMask = phantomMask.LCMask;
LMMask = phantomMask.LMMask;

phantom3d(AMask) = phantom3d(AMask)*targetROI.ATarget(timeInd)/mean(phantom3d(AMask));
phantom3d(RCMask) = phantom3d(RCMask)*targetROI.RCTarget(timeInd,gfrInd)/mean(phantom3d(RCMask));
phantom3d(RMMask) = phantom3d(RMMask)*targetROI.RMTarget(timeInd,gfrInd)/mean(phantom3d(RMMask));
phantom3d(LCMask) = phantom3d(LCMask)*targetROI.LCTarget(timeInd,gfrInd)/mean(phantom3d(LCMask));
phantom3d(LMMask) = phantom3d(LMMask)*targetROI.LMTarget(timeInd,gfrInd)/mean(phantom3d(LMMask));

end