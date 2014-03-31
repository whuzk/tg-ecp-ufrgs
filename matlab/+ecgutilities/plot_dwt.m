function h = plot_dwt(c,l,lev,len,h)

nbcol = 64;
if ~iscell(c)
    c = detcoef(c,l,1:lev);
end
cfd = dwt_matrix(c,lev,len);
cfd = wcodemat(cfd,nbcol,'row');

if isempty(h)
    colormap(pink(nbcol));
    h = image(cfd);
    tics = 1:lev;
    labs = int2str(tics');
    setup_plot(tics,labs);
else
    set(h, 'cdata', cfd);
end


function cfd = dwt_matrix(d,lev,len)
cfd = zeros(lev,len);
for k = 1:lev
    dk = d{k}(:)';
    dk = dk(ones(1,2^k),:);
    cfd(k,:) = wkeep1(dk(:)',len);
end
cfd =  cfd(:);
I = find(abs(cfd)<sqrt(eps));
cfd(I) = zeros(size(I));
cfd = reshape(cfd,lev,len);

function setup_plot(tics,labs)
set(gca,...
    'YTickLabelMode','manual','YDir','normal', ...
    'Box','On','YTick',tics,'YTickLabel',labs);
title('Discrete Wavelet Transform, Absolute Coefficients.');
xlabel('Time (or Space)')
ylabel('Level');