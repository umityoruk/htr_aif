function [mruImNew] = paper_removeFlicker(mruIm)

Nx = size(mruIm,1);
Ny = size(mruIm,2);
Nz = size(mruIm,3);
T = size(mruIm,4);

mruImAve = squeeze(mean(mruIm,4));

mruImNew = zeros(size(mruIm),class(mruIm));
figure
for zz = 1:Nz
    for yy = 1:Ny
        for xx = 1:Nx
%             sig = squeeze(mruIm(xx,yy,zz,:));
%             smoothSig = smooth(sig,0.1,'lowess');
            if (xx >= 40 && yy >= 40 && zz >= 17)
                sig = squeeze(mruIm(xx,yy,zz,:));
                smoothSig = smooth(sig,0.1,'lowess');
                subplot(1,2,1)
                imagesc(mruImAve(:,:,zz))
                colormap('gray')
                hold on
                plot(xx,yy,'r+', 'MarkerSize', 20,'LineWidth',2);
                hold off
                subplot(1,2,2)
                plot(sig,'k-')
                hold on
                plot(smoothSig,'r-','LineWidth',2);
                hold off
                title(['(' num2str(xx) ' ' ...
                           num2str(yy) ' ' ...
                           num2str(zz) ')']);
                pause
            end
%             mruImNew(xx,yy,zz,:) = smoothSig;
        end
    end
end

end