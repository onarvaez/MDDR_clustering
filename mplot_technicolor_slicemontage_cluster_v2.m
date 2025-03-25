function mplot_technicolor_slicemontage_cluster(method, dps, slices_path, clim, opt,slmoNrow,slmol_v,slmob_v)
% function mplot_technicolor_slicemontage(method, dps, slices_path, clim, opt)

if nargin < 5
    opt = mdm_opt();
end
opt = mplot_opt(opt);

%Plot parameter maps

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;

if isfield(opt,'k_range')
    k_range = opt.k_range;
else
    k_range = 1:sz(3); 
end


% if isfield(opt,'k_range')
%     k_range = opt.k_range;
% else
%     k_range = 1:sz(3); 
% end
% if isfield(opt,'j_range')
%     j_range = opt.j_range;
% else
%     j_range = 1:sz(2); 
% end
% if isfield(opt,'i_range')
%     i_range = opt.i_range;
% else
%     i_range = 1:sz(1); 
% end


s0_nk = dps.s0(:,:,k_range);
s_b2000_nk = dps.s_b2000(:,:,k_range);
smax = quantile(s0_nk(s0_nk>0),.999,'all');
clim.s2000 = clim.s0;
clim.s0 = smax*clim.s0;
clim.s2000 = quantile(s_b2000_nk(s0_nk>0),.999,'all')*clim.s2000;

mask = dps.s0 > clim.mask_threshold*smax;

pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

Nbins = numel(dps.bin);

if strcmp(method,'dtod')
    plotfields.gray = {'s0';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'; 'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta';'dmdisodnu';'dmsddeltadnu'};
    
    Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
    Ncolumns = 7;
    Nrows = 8;
    papersize = [Ncolumns/opt.mplot.zoomx Nrows*pixaspect*imaspect/opt.mplot.zoomy];
    papersize = papersize/papersize(1)*17.56;

    dbottom = 1/Nrows;
    dleft = 1/Ncolumns;

    position.height = 1.01*dbottom; 
    position.width = 1.01*dleft;
    
%     position.left_v =   [0    1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 1.5 2.5 3.5 4.5 6   7   6   7   6   1.5 1.5 1.5 2.5 2.5 2.5 3.5 3.5 3.5 4.5 4.5 4.5 5.5 5.5 5.5 6.5 6.5 6.5 0]/Ncolumns;
%     position.bottom_v = [3.75 5.5 5.5 5.5 5.5 4.5 4.5 4.5 4.5 3.5 3.5 3.5 3.5 5.5 5.5 4.5 4.5 3.5 2   1   0   2   1   0   2   1   0   2   1   0   2   1   0   2   1   0   1.75]/Nrows;
    position.left_v =   [0    1.5 2.5 1.5 2.5 2   5   6   5   6   5.5   ...
        1.5 1.5 1.5 2.5 2.5 2.5 5   5   5   6   6   6   3.5 3.5 3.5 0]/Ncolumns;
    
    position.bottom_v = [3.75 5.5 5.5 4.5 4.5 3.5 5.5 5.5 4.5 4.5 3.5 ...
        2   1   0   2   1   0   2   1   0   2   1   0   2   1   0   1.75]/Nrows;

elseif strcmp(method,'dtor1r2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvr1r2'; 'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
%   plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2';'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta';'mr1';'mr2';'dmdisodnu';'dmsddeltadnu'};


%6clust
    Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
    Ncolumns = 9;
    Nrows = slmoNrow;
    papersize = [Ncolumns/opt.mplot.zoomx Nrows*pixaspect*imaspect/opt.mplot.zoomy];
    papersize = papersize/papersize(1)*17.56;

    dbottom = 1/Nrows;
    dleft = 1/Ncolumns;

    position.height = 1.01*dbottom; 
    position.width = 1.01*dleft;
    
    position.left_v =   slmol_v/Ncolumns;
    
    position.bottom_v = slmob_v/Nrows;
      
end


msf_mkdir(fullfile(slices_path));

for nk = k_range
% for nk = 2
    figure(2), clf
    map_count = 1;
    axh_v = [];
    for c = 1:numel(plotfields.gray)
        im3d = double(dps.(plotfields.gray{c}));
        im2d = mask(:,:,nk).*im3d(:,:,nk);
        im2d(im2d>.99*max(clim.(plotfields.gray{c}))) = .99*max(clim.(plotfields.gray{c})); %Factor .99 to avoid black speckles in pdf
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);
        imagesc(im2d')
        colormap(axh,gray(256))
        set(axh,'CLim',clim.(plotfields.gray{c}))
        axh_v = [axh_v; axh];
        map_count = map_count+1;
    end

    for c = 1:numel(plotfields.hotcold)
        im3d = double(dps.(plotfields.hotcold{c}));
        im2d = mask(:,:,nk).*im3d(:,:,nk);
        im2d(im2d>.99*max(clim.(plotfields.hotcold{c}))) = .99*max(clim.(plotfields.hotcold{c})); %Factor .99 to avoid black speckles in pdf
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);

        imagesc(im2d')
        colormap(axh,mplot_cmaphotcold(128))
        set(axh,'CLim',clim.(plotfields.hotcold{c}))   
        axh_v = [axh_v; axh];
        map_count = map_count+1;

    end

    
    for c = 1:numel(plotfields.bin)
        for nbin = 1:Nbins
            clear im3d
            cind = (dps.bin{nbin}.(plotfields.bin{c})-min(clim.(plotfields.bin{c})))...
                /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
            im3d = dist_cind2rgb_jet(cind);
            im3d.bright = mask.*(dps.bin{nbin}.f);
            im2d = zeros(sz(2),sz(1),3);
            im2d(:,:,1) = im3d.r(:,:,nk)';
            im2d(:,:,2) = im3d.g(:,:,nk)';
            im2d(:,:,3) = im3d.b(:,:,nk)';
            im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]);
            position.left = position.left_v(map_count);
            position.bottom = position.bottom_v(map_count);
            axh = axes('position',[position.left position.bottom position.width position.height]);
            
            imagesc(im2d)
            set(axh,'CLim',[0 1])   
            axh_v = [axh_v; axh];
            map_count = map_count+1;
        end
    end

    for nbin = 1:Nbins
        clear im3d
        im3d.r = dps.bin{nbin}.mdxx;
        im3d.g = dps.bin{nbin}.mdyy;
        im3d.b = dps.bin{nbin}.mdzz;
        im3d.bright = mask.*(dps.bin{nbin}.f);
        im2d = zeros(sz(2),sz(1),3);
        im2d(:,:,1) = im3d.r(:,:,nk)';
        im2d(:,:,2) = im3d.g(:,:,nk)';
        im2d(:,:,3) = im3d.b(:,:,nk)';
        im2d = .99*im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]); %factor .99 added to avoid black voxels in pdf
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);

        imagesc(im2d)
        set(axh,'CLim',[0 1])   
        axh_v = [axh_v; axh];
        map_count = map_count+1;
    end

    clear im3d
    im3d.r = (dps.bin{5}.f);
    im3d.g = (dps.bin{1}.f);
    im3d.b = (dps.bin{4}.f);
    %im3d.bright = mask.*((dps.bin{1}.f) + (dps.bin{2}.f) + (dps.bin{3}.f));
    %im3d.bright = mask;
%     im3d.bright
    im2d = zeros(sz(2),sz(1),3);
    im2d(:,:,1) = im3d.r(:,:,nk)';
    im2d(:,:,2) = im3d.g(:,:,nk)';
    im2d(:,:,3) = im3d.b(:,:,nk)';
    %im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]);
    position.left = position.left_v(map_count);
    position.bottom = position.bottom_v(map_count);
    axh = axes('position',[position.left position.bottom position.width position.height]);

    imagesc(im2d)
    set(axh,'CLim',[0 1])   
    axh_v = [axh_v; axh];
    map_count = map_count+1;

    if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
        set(axh_v,'YDir','reverse')
    elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
        set(axh_v,'XDir','reverse','YDir','normal')
    else
        set(axh_v,'YDir','normal')
    end
    axis(axh_v,'tight','off')
    set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])

    % Temporary grid lines for guiding placement of labels
    if (0)
        axh_grid = axes('position',[0 0 1 1]);
        xgrid_x = repmat([0; Ncolumns], [1 4*Nrows+1]);
        xgrid_y = repmat((0:.25:Nrows), [2 1]);
        ygrid_x = repmat((0:.25:Ncolumns), [2 1]);
        ygrid_y = repmat([0; Nrows], [1 4*Ncolumns+1])
        plot(axh_grid,xgrid_x, xgrid_y, 'g-') 
        hold(axh_grid,'on')
        plot(axh_grid,ygrid_x, ygrid_y, 'g-') 
        axis(axh_grid,'tight','off')
    end

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
    %fig_name = ['slice' num2str(slice_select)];
    fig_name = ['slice' num2str(nk)];
    print(fullfile(slices_path,fig_name),'-loose','-dsvg','-r1200')
    print(fullfile(slices_path,fig_name),'-loose','-dpdf')
    save(fullfile(slices_path,'clim.mat'),'clim')
    %print(fullfile(slices_path,fig_name),'-loose','-dpng')
    %savefig(gcf,fullfile(slices_path,fig_name),'compact')
end