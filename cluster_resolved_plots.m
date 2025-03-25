%% Prepare paths
clearvars; close all; clc;


%bootstraps path
bs_path = '.../dtor1r2d/bootstraps';

%select frequency range
method = 'dtor1r2d'; omega_v = [1 100]*2*pi;
clust_num = [6];

%for slice_select = 1:numel(slices)
%------------
tic;
pdata_path = fileparts(fileparts(bs_path));
roi_select = load(fullfile(bs_path,'1/mfs.mat'));

roi = roi_select.mfs.mask; 

clust_dir = fullfile(pdata_path,'clusters/');


load(fullfile(clust_dir, 'my_clusters_v3.mat'));


cl_mask = clust_stru{1,clust_num-1}.clmask;
num_clust = clust_stru{1,clust_num-1}.Nclusters-1;
slmoNrow = clust_stru{1,clust_num-1}.slmo_grid.Nrow;
slmol_v = clust_stru{1,clust_num-1}.slmo_grid.left_v;
slmob_v = clust_stru{1,clust_num-1}.slmo_grid.bottom_v;
clust_order = clust_stru{1,clust_num-1}.clorder;
opt = mdm_opt();
opt.dtod.maps_omega = omega_v;
opt.(method).maps_omega = omega_v;

cls_path = cell(0,0);  
for cldata = 1:numel(clust_order)
    
       
    cls_path{1+numel(cls_path)} = clust_stru{1,clust_num-1}.clustFracN(:,:,:,clust_order(cldata)); 

end  


pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(pdata_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(pdata_path,method,'bootstraps');        
end


%% Make figures

%Set all the color limits 

disp('Setting color limits . . . .')

%rat invivo 
clim.s0 = 1*[0 1]; % Multiplied with s0 below
clim.mdiso = 4e-9*[0 1]; %1e-9*[0 1]
clim.msddelta = 1*[0 1];
clim.mr1 = 0.8*[0 1]; %1 %0.6 in first figure
clim.mr2 = 35*[0 1];
clim.vdiso = 1e-18*[0 1];
clim.vsddelta = .10*[0 1];
clim.vr1 = 0.2*[0 1];
clim.vr2 = 200*[0 1];


clim.cvdisosddelta = sqrt(max(clim.vdiso)*max(clim.vsddelta))*[-1 1]; %clim.cvdisosddelta = .2e-9*[-1 1];
clim.cvdisor1 = sqrt(max(clim.vdiso)*max(clim.vr1))*[-1 1];
clim.cvdisor2 = sqrt(max(clim.vdiso)*max(clim.vr2))*[-1 1];
clim.cvsddeltar1 = sqrt(max(clim.vsddelta)*max(clim.vr1))*[-1 1];
clim.cvsddeltar2 = sqrt(max(clim.vsddelta)*max(clim.vr2))*[-1 1];
clim.cvr1r2 = sqrt(max(clim.vr1)*max(clim.vr2))*[-1 1];
clim.mask_threshold = eps; %clim.mask_threshold = 0.01;


%rat
clim.dmdisodnu = 1.5e-3*max(clim.mdiso)*[-1 1];
clim.dmsddeltadnu = 1.5e-3*max(clim.msddelta)*[-1 1];
clim.dvdisodnu = 1e-3*max(clim.vdiso)*[-1 1];
clim.dvsddeltadnu = 1e-3*max(clim.vsddelta)*[-1 1];
clim.dcvdisosddeltadnu = 1e-2*max(clim.cvdisosddelta)*[-1 1];

%%
disp('Loading bootstraps to get the cluster-resolved means . . . .')

% Loop over datasets
for Nclust = 1:clust_num    
    
    %cl_or = clust_order(Nclust);
    bs_dps = mdm_dps_collectbs_cluster(method, bs_paths{1}, opt, cl_mask,Nclust);
    %bs_dps = mdm_dps_collectbs_cluster(method, bs_paths{1}, opt, cl_mask,slice_select,Nclust);

        
    if ~all(cellfun('isempty',bs_dps))
       % median_dps = mdm_dps_median_cluster(bs_dps,1);
        median_dps = mdm_dps_median(bs_dps);
        clear bs_dps
 
    end
    
    all_dps.bin{Nclust} = median_dps;  
end

for Nclust = 1:numel(clust_order)
select = Nclust;

clust_select = clust_order(select);
cluster_fraction = clust_stru{1,clust_num-1}.clustFracN(:,:,:,clust_select);
%cluster_fraction = clust_stru{1,clust_num-1}.clustFrac(:,:,:,clust_select);
dps_cf.f{Nclust} = cluster_fraction;
end


disp('Loading bootstraps one more time (sorry memory) to get the per-voxel means . . . .')

bs_dps = mdm_dps_collectbs_global_v2(method, bs_paths{1}, opt, cls_path);
        
 if ~all(cellfun('isempty',bs_dps))
        %median_dps = mdm_dps_median_cluster(bs_dps,slice_select);
        median_dps = mdm_dps_median(bs_dps);

        clear bs_dps
 

end



median_dps.bin = all_dps.bin;

for cf = 1:clust_num
median_dps.bin{cf}.f = dps_cf.f{1,cf};

end
%%
disp('Making figures . . . .')

%figures and maps
%mplot_technicolor_nii(method, median_dps, pmaps_paths{ndata}, clim, opt)
slice_select = 1;
mplot_technicolor_nii_cluster_v2(method, median_dps, pmaps_paths{ndata}, clim, opt,slice_select)
mplot_technicolor_slicemontage_cluster_v2(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),['slicemontage_cluster_' int2str(clust_num)]), clim, opt,slmoNrow,slmol_v,slmob_v)
mplot_globalstats(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),[method '_globalstats']), clim)
%     
t1 = toc;
disp(['Done. ' 'Total time: ' int2str(t1) ' seconds'])

  
