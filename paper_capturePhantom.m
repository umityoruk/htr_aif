
sliceSelect = 17;
im0 = newPhantom(:,:,sliceSelect,20);
im1 = newPhantom(:,:,sliceSelect,40);
im2 = newPhantom(:,:,sliceSelect,120);


% Normalize the images
maxVal = max([max(max(im0)); max(max(im1)); max(max(im2))]);

brightness = 1.0;
im0 = im0/maxVal*brightness;
im1 = im1/maxVal*brightness;
im2 = im2/maxVal*brightness;

pullDown = 0.00;
im0 = im0 - pullDown; im0(im0<0) = 0;
im1 = im1 - pullDown; im1(im1<0) = 0;
im2 = im2 - pullDown; im2(im2<0) = 0;


cutImages = true;
if (cutImages)
    leftCut = 30;
    rightCut = 25;
    topCut = 0;
    bottomCut = 100;
    im0 = im0(topCut+1:end-bottomCut,leftCut+1:end-rightCut);
    im1 = im1(topCut+1:end-bottomCut,leftCut+1:end-rightCut);
    im2 = im2(topCut+1:end-bottomCut,leftCut+1:end-rightCut);
end

imwrite(im0,'phantomIm0.png')
imwrite(im1,'phantomIm1.png')
imwrite(im2,'phantomIm2.png')