%% Set the bootstrap path and options
clearvars; close all; clc;
%path to bootstraps
param.bs_path = '..../dtor1r2d/bootstraps';

%%Parameters
param.slices = 1; % Select the slice of interest (select 1 if it is single slice)
param.num_clust = [2 3 4 5 6 7 8 9 10]; % Select number of clusters
param.clust_replicates = 3; % Selecet number of replicates
% change omega values depending of protocol
param.method = 'dtor1r2d'; param.omega_10_90 = [1 100]*2*pi; param.omega_50 = [50]*2*pi; 
% save Path
param.clust_sv_path = [fileparts(fileparts(param.bs_path)) filesep 'clusters' filesep];

%% Load Boostrap
%[boot_stru,~] = my_load_bootstrap_function(bs_path,method);
[boot_stru,~] = my_load_bootstrap_function(param.bs_path,param.method);

param.Size_img = [size(boot_stru.w_4d,1) size(boot_stru.w_4d,2) size(boot_stru.w_4d,3)];
% Get Mask
param.Mask = logical(sum(boot_stru.w_4d,4)~=0);
param.Mask = repmat(param.Mask,1,1,1,size(boot_stru.w_4d,4));

%% Compute omega dependant parameters
[dparo_4d_low,dperpo_4d_low] = my_dtor1r2d_dist2parso(boot_stru.dpar_4d,boot_stru.dperp_4d,boot_stru.rpar_4d,boot_stru.rperp_4d,boot_stru.d0_4d,param.omega_10_90(1,1));
[disoo_4d_low,~,~,ddeltao_4d_low,~,~] = my_dtd_pars2dpars(dparo_4d_low,dperpo_4d_low);
[dparo_4d_high,dperpo_4d_high] = my_dtor1r2d_dist2parso(boot_stru.dpar_4d,boot_stru.dperp_4d,boot_stru.rpar_4d,boot_stru.rperp_4d,boot_stru.d0_4d,param.omega_10_90(1,2));
[disoo_4d_high,~,~,ddeltao_4d_high,~,~] = my_dtd_pars2dpars(dparo_4d_high,dperpo_4d_high);

boot_stru.Diff_disoo_4d = disoo_4d_high-disoo_4d_low;
boot_stru.Diff_ddeltao_4d = ddeltao_4d_high-ddeltao_4d_low;
clearvars dparo_4d_low dperpo_4d_low disoo_4d_low sddeltao_4d_low
clearvars dparo_4d_high dperpo_4d_high disoo_4d_high sddeltao_4d_high

[boot_stru.dparo_4d,boot_stru.dperpo_4d] = my_dtor1r2d_dist2parso(boot_stru.dpar_4d,boot_stru.dperp_4d,boot_stru.rpar_4d,boot_stru.rperp_4d,boot_stru.d0_4d,param.omega_50);
[boot_stru.disoo_4d,~,~,boot_stru.ddeltao_4d,~,boot_stru.sddeltao_4d] = my_dtd_pars2dpars(boot_stru.dparo_4d,boot_stru.dperpo_4d);

%% Prepare clust_stru : get mask bootstrap only, mean and std normalisation
%param.clust_variables = {'Diff_disoo'}; % 'w' ,'Diff_ddeltao'
param.clust_variables = {'disoo','ddeltao','r1','r2','Diff_disoo'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'disoo','ddeltao','r1','r2','Diff_disoo','Diff_ddeltao'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'disoo','sddeltao','r1','r2','Diff_disoo'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'disoo','sddeltao','r1','r2','Diff_disoo','Diff_ddeltao'}; % 'w' 
%param.clust_variables = {'disoo','sddeltao','Diff_disoo'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'disoo','ddeltao'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'r1','r2'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'disoo','ddeltao','r1','r2'}; % 'w' ,'Diff_ddeltao'
%param.clust_variables = {'dparo_4d_low','dparo_4d_high','dperpo_4d_low','dperpo_4d_high'}; % 'w' ,'Diff_ddeltao'

param.clust_Nvariables = size(param.clust_variables,2);
clustData = single(zeros(sum(param.Mask(:)),param.clust_Nvariables));
for ind_var = 1:param.clust_Nvariables
    % Select mask bootstrap only
    clustData(:,ind_var) = single(boot_stru.([param.clust_variables{1,ind_var} '_4d'])(param.Mask));    %clustData(:,ind_var) = single(boot_stru.([param.clust_variables{1,ind_var}])(param.Mask));

    % Remove mean value
    %clustData(:,ind_var) = clustData(:,ind_var)- mean(clustData(:,ind_var));
    % Divide by std value
    clustData(:,ind_var) = clustData(:,ind_var)/std(clustData(:,ind_var));
end

%% weight the cluster data
weights = single(boot_stru.w_4d(param.Mask));
max_weights = quantile(weights,0.999).*1;
weights = weights./max_weights;
weights(weights>1) = 0;
weights = round(weights.*100);
% hist(weights,1000);

clustDataW = single(zeros([100 size(clustData)]));
log_mask = false(size(clustDataW));
for ind_param = 1:size(clustData,2)
for ind = 1:size(weights,1)
    if weights(ind,1) ~=0
    clustDataW(1:weights(ind,1),ind,ind_param) = repmat(clustData(ind,ind_param),weights(ind,1),1);
    log_mask(1:weights(ind,1),ind,ind_param) = true(weights(ind,1),1);
    end
end
end
clustDataW = clustDataW(log_mask);
clustDataW = reshape(clustDataW,size(clustDataW,1)./param.clust_Nvariables,param.clust_Nvariables);
clearvars log_mask

%% clustering
if ~exist([param.clust_sv_path 'my_clusters_v3.mat'],'file')
    for ind_clu = 1:size(param.num_clust,2)
    %for ind_clu = 5
         disp('Clustering . . . .')
         disp(ind_clu + 1)
        tic
        param.opts = statset('Display','final','MaxIter',500);
        clust_stru{ind_clu}.GMModel = fitgmdist(clustDataW,param.num_clust(1,ind_clu),'regularization',0.01,'Replicates',param.clust_replicates,'Options',param.opts);
        [idx,nlog,p,logpdf] = cluster(clust_stru{ind_clu}.GMModel,clustData);
        toc
 
        %%save probabilities
        %prob = sum(exp(logpdf),2);
        %% reshape cluster indices to 4D matrix
        indices_4d = find(param.Mask);
        [dim1,dim2,dim3,dim4] = ind2sub(size(param.Mask),indices_4d);
        clust_stru{ind_clu}.idx_4d = zeros(size(param.Mask));
        for ind_clust = 1:numel(idx)
            clust_stru{ind_clu}.idx_4d(dim1(ind_clust),dim2(ind_clust),dim3(ind_clust),dim4(ind_clust)) = idx(ind_clust);
            %clust_stru{ind_clu}.p_4d(dim1(ind_clust),dim2(ind_clust),dim3(ind_clust),dim4(ind_clust)) = prob(ind_clust);
        end
        clust_stru{ind_clu}.Nclusters = param.num_clust(1,ind_clu);
        clearvars idx dim1 dim2 dim3 dim4 indices_4d ind_clust;
    end
    save([param.clust_sv_path 'my_clusters_v3.mat'],'clust_stru','-v7.3')
else
    load([param.clust_sv_path 'my_clusters_v3.mat']);
end

%% Compute cluster fractions

for ind_Nclu = 1:size(param.num_clust,2)
    for ind_clu = 1:clust_stru{ind_Nclu}.Nclusters
        idx_4d_tmp = zeros(size(clust_stru{ind_Nclu}.idx_4d));
        idx_4d_tmp(clust_stru{ind_Nclu}.idx_4d==ind_clu) =1;
        clust_stru{ind_Nclu}.clustFrac(:,:,:,ind_clu) = sum(idx_4d_tmp,4)./size(clust_stru{ind_Nclu}.idx_4d,4);
        k_mclu = clust_stru{ind_Nclu}.clustFrac(:,:,:,ind_clu);
        clust_stru{ind_Nclu}.clustFracN(:,:,:,ind_clu) = (k_mclu-min(k_mclu(:))) ./ (max(k_mclu(:)-min(k_mclu(:))));
    end
end


for ind_Nclu = 1:size(param.num_clust,2)
    for ind_clu = 1:clust_stru{ind_Nclu}.Nclusters-1
        clust_stru{ind_Nclu}.slmo_grid = slicemontage_grid(clust_stru{ind_clu}.Nclusters);
    end
end

%% Make and save bootstraps mask 
for ind_Nclu = 1:size(param.num_clust,2)

% figure(1)
% for ind_Ncluster = 1:clust_stru{1,ind_Nclu}.Nclusters
%     subplot(1,clust_stru{1,ind_Nclu}.Nclusters,ind_Ncluster)
%     imagesc(clust_stru{1,ind_Nclu}.clustFrac(:,:,param.slices,ind_Ncluster));
%     axis image
%     xticks([]); yticks([]);
% end
% prompt = {'Enter the order of the clusters (only numbers, comma separated) acording to the displayed figure (first three: wm, gm, csf/pbs):'};
% dlgtitle = 'Clusters';
% wdims = [1 70];
% answer = inputdlg(prompt,dlgtitle,wdims);
[prop,clust_order] = sort(clust_stru{ind_Nclu}.GMModel.ComponentProportion,'descend');

%clust_order = [str2num(answer{1,1})];

for group = 1:numel(1:clust_stru{1,ind_Nclu}.Nclusters)

    gru = clust_order(group);
    %mapa = reshape(clust_stru{1,1}.idx_4d,size(disoo,1),[]);
    mapa = clust_stru{1,ind_Nclu}.idx_4d;
    %separate each cluster and make it binary
    mapa(mapa~=gru)= 0;
    mapa(mapa==gru)=1;
    %multiply each cluster mask to each parameter_bootstrap solution
    
    boot_mask = mapa;
    mask_bps.clust{gru} = boot_mask;
   
end


% Reorder the whole bootstrap solutions into separated 10 solutions that 
%correspond to the 10 solutions per bootstrap folder 

clust_lo = mask_bps.clust;
b_no= size(mapa,4)/10;
rep_boot = 10;




sep_vec = 1:10:b_no*rep_boot;

for Nclust = 1:clust_stru{1,ind_Nclu}.Nclusters
    
    clust = clust_lo{1,clust_order(Nclust)};
    
    clust_folder = zeros((size(clust,1)),(size(clust,2)),1,b_no);
    
    for Nfolder = 1:b_no

        clust_folder_temp = clust(:,:,:,sep_vec(Nfolder):(sep_vec(Nfolder)+9));
        
        cluster_fl.folder{Nfolder} = clust_folder_temp;
    end

cl_ten.cf{Nclust} = cluster_fl.folder;


end

% save the mask and other options for the next step
clust_stru{ind_Nclu}.clmask = cl_ten.cf;
clust_stru{ind_Nclu}.clorder = clust_order;

end

save([param.clust_sv_path 'my_clusters_v3.mat'],'clust_stru','-v7.3')
%% display STD decay curves
for ind_Nclu = 1:size(param.num_clust,2)
    idx_lin = clust_stru{ind_Nclu}.idx_4d(param.Mask);
    [gc,~] = groupcounts(idx_lin);
    for ind_clu = 1:clust_stru{ind_Nclu}.Nclusters
        std_clu_tmp(ind_clu,:) = std(clustData(idx_lin==ind_clu,:),1);
    end
    std_clu(ind_Nclu,:) = (sum(std_clu_tmp.*gc,1))./sum(gc);
end
std_clu = cat(1,ones(1,size(std_clu,2)),std_clu);

figure(1)
set(gcf,'color','w')
title('Standard deviation curves')
plot([1 param.num_clust], std_clu)
ylabel('Standard deviation')
xlabel('N cluster')
legend(param.clust_variables)

%% Display cluster fractions images
%param.slices = 2;
figure(2)
set(gcf,'color','w')
set(gcf,'Position',[744,49.8,937.8,1020.8])
MaxClusters = 0;
for ind_Nclu = 1:size(param.num_clust,2)
    MaxClusters = max([MaxClusters clust_stru{ind_Nclu}.Nclusters]);
end
for ind_Nclu = 1:size(param.num_clust,2)
    [~,sort_order_clu] = sort(clust_stru{ind_Nclu}.GMModel.ComponentProportion,'descend');
    sort_prop = sort(clust_stru{ind_Nclu}.GMModel.ComponentProportion,'descend');
    for ind_clu = 1:clust_stru{ind_Nclu}.Nclusters
        subplot(MaxClusters,size(param.num_clust,2),(ind_clu-1)*(size(param.num_clust,2))+(ind_Nclu))
        imagesc(flip(flip(clust_stru{ind_Nclu}.clustFracN(:,:,param.slices,sort_order_clu(ind_clu))',1),2));
        xL=xlim;
        yL=ylim;
        text(0.99*xL(2),min(0.99*yL(1)),num2str(round(sort_prop(ind_clu).*100)),'HorizontalAlignment','right','VerticalAlignment','top','Color','white','FontSize',7)
        colormap jet
        axis image
        if ind_clu ==1
            title(num2str(clust_stru{ind_Nclu}.Nclusters));
        end
        xticks([]); yticks([]);
    end
end

%% Display maps and distributions
% parameters
param.Nbins = 128;
param.edgesDiso = 0:5/(param.Nbins-1):5;
%param.edgesDiso = 0:2/(param.Nbins-1):2;
param.edgesDdelta = -0.5:1.5/(param.Nbins-1):1;
%param.edgesDdelta = 0:1/(param.Nbins-1):1;
%cx
%param.edgesR1 = 0:.6/(param.Nbins-1):.6;
%wm
param.edgesR1 = 0:1.5/(param.Nbins-1):1.5;
%wm
param.edgesR2 = 0:80/(param.Nbins-1):80;
%cx
%param.edgesR2 = 0:60/(param.Nbins-1):60;
param.edgesDiffDiso = -0.1:0.3/(param.Nbins-1):0.3;
%param.edgesDiffDiso = -0.1:0.1/(param.Nbins-1):0.1;
param.edgesDiffDdelta = -0.5:1/(param.Nbins-1):0.5;
param.Variable_to_display = 6;

% Figures
for ind_Nclu = 1:size(param.num_clust,2)
    figure('Name',['Images and distribution for ' num2str(clust_stru{ind_Nclu}.Nclusters) ' clusters'])
    set(gcf,'color','w')
    [~,sort_order_clu] = sort(clust_stru{ind_Nclu}.GMModel.ComponentProportion,'descend');
    % display fractions images
    for ind_clu = 1:clust_stru{ind_Nclu}.Nclusters
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,1+(ind_clu-1)*(param.Variable_to_display+1))
        imagesc(flip(flip(clust_stru{ind_Nclu}.clustFracN(:,:,param.slices,sort_order_clu(ind_clu))',1),2));
        colormap jet
        axis image
        if ind_clu ==1
            title(num2str(clust_stru{ind_Nclu}.Nclusters));
        end
        xticks([]); yticks([]);
        %% Display Component distribution
        Mask_cmp = zeros(size(clust_stru{ind_Nclu}.idx_4d));
        Mask_cmp(clust_stru{ind_Nclu}.idx_4d==sort_order_clu(ind_clu)) =1;
        Mask_cmp = logical(Mask_cmp);
        
        % Diso
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,2+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesDiso,dist2bar(param.edgesDiso,boot_stru.disoo_4d./1e-9,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesDiso,dist2bar(param.edgesDiso,boot_stru.disoo_4d(Mask_cmp)./1e-9,boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1 + max(dist2bar(param.edgesDiso,boot_stru.disoo_4d(Mask_cmp)./1e-9,boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('D_{iso}/10^{-9}','FontSize', 8);
        end
        
        % Ddelta²
%         subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,3+(ind_clu-1)*(param.Variable_to_display+1))
%         hold on
%             bar(param.edgesDdelta,dist2bar(param.edgesDdelta,boot_stru.sddeltao_4d,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
%             bar(param.edgesDdelta,dist2bar(param.edgesDdelta,boot_stru.sddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
%         yticks([]);
%         ylim([0 1+max(dist2bar(param.edgesDdelta,boot_stru.sddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)))*1.2])
%         if ind_clu ==1
%             title('D_{\Delta}^{2}','FontSize', 8);
%         end
        
        % Ddelta
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,3+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesDdelta,dist2bar(param.edgesDdelta,boot_stru.ddeltao_4d,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesDdelta,dist2bar(param.edgesDdelta,boot_stru.ddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1+max(dist2bar(param.edgesDdelta,boot_stru.ddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('D_{\Delta}','FontSize', 8);
        end
        
        % R1
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,4+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesR1,dist2bar(param.edgesR1,boot_stru.r1_4d,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesR1,dist2bar(param.edgesR1,boot_stru.r1_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1+max(dist2bar(param.edgesR1,boot_stru.r1_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('R1 (s-1)','FontSize', 8);
        end
        
        % R2
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,5+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesR2,dist2bar(param.edgesR2,boot_stru.r2_4d,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesR2,dist2bar(param.edgesR2,boot_stru.r2_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1+max(dist2bar(param.edgesR2,boot_stru.r2_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('R2 (s-1)','FontSize', 8);
        end
        
        % Delta Diso
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,6+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesDiffDiso,dist2bar(param.edgesDiffDiso,boot_stru.Diff_disoo_4d./1e-9,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesDiffDiso,dist2bar(param.edgesDiffDiso,boot_stru.Diff_disoo_4d(Mask_cmp)./1e-9,boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1+max(dist2bar(param.edgesR2,boot_stru.Diff_disoo_4d(Mask_cmp)./1e-9,boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('{\Delta}D_{iso}','FontSize', 8);
        end
        
        % Delta Ddelta
        subplot(clust_stru{ind_Nclu}.Nclusters,param.Variable_to_display+1,7+(ind_clu-1)*(param.Variable_to_display+1))
        hold on
            bar(param.edgesDiffDdelta,dist2bar(param.edgesDiffDdelta,boot_stru.Diff_ddeltao_4d,boot_stru.w_4d),'FaceColor',[0 0 0],'EdgeColor','none','BarWidth', 1)
            bar(param.edgesDiffDdelta,dist2bar(param.edgesDiffDdelta,boot_stru.Diff_ddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)),'FaceColor',[1 0 0],'EdgeColor','none','BarWidth', 1)
        yticks([]);
        ylim([0 1+max(dist2bar(param.edgesR2,boot_stru.Diff_ddeltao_4d(Mask_cmp),boot_stru.w_4d(Mask_cmp)))*1.2])
        if ind_clu ==1
            title('{\Delta}D_{\Delta}','FontSize', 8);
        end
    end
end




%%%%%%%   DTDOR1R2 functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[boot_stru,m_bsmerge] = my_load_bootstrap_function(bs_path,method)
bsno = msf_getdirno(bs_path);
m_bsmerge = [];
for nbs = 1:numel(bsno)
    mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
    if exist(mfs_fn,'file')==2
        mfs = mdm_mfs_load(mfs_fn);
        m_bsmerge = cat(4,m_bsmerge,mfs.m);
        if nbs ==1
            size_bsn = size(mfs.m);
            if strcmp(method,'dtor1r2d')
                size_bsn(1,4) = (size_bsn(1,4)-1)/10;
                [boot_stru.dpar_4d,boot_stru.dperp_4d,boot_stru.theta_4d,boot_stru.phi_4d,boot_stru.d0_4d,boot_stru.rpar_4d,boot_stru.rperp_4d,boot_stru.r1_4d,boot_stru.r2_4d,boot_stru.w_4d] = my_dtor1r2d_bsmerge_initializepars(size_bsn,numel(bsno));
            end
        end
        if strcmp(method,'dtor1r2d')
            [dpar_temp,dperp_temp,theta_temp,phi_temp,d0_temp,rpar_temp,rperp_temp,r1_temp,r2_temp,w_temp] = dtor1r2d_4d_m2pars(mfs.m);
            boot_stru.dpar_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(dpar_temp);
            boot_stru.dperp_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(dperp_temp);
            boot_stru.theta_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(theta_temp);
            boot_stru.phi_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(phi_temp);
            boot_stru.d0_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(d0_temp);
            boot_stru.rpar_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(rpar_temp);
            boot_stru.rperp_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(rperp_temp);
            boot_stru.r1_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(r1_temp);
            boot_stru.r2_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(r2_temp);
            boot_stru.w_4d(:,:,:,1+(nbs-1)*size_bsn(1,4):nbs*size_bsn(1,4)) = single(w_temp);
        end
    end
end
end

function [dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2,w] = my_dtor1r2d_bsmerge_initializepars(size_bsn,Nbsn)
dpar = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
dperp = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
theta = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
phi = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
d0 = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
rpar = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
rperp = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
r1 = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
r2 = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
w = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
end

function [dparo,dperpo] = my_dtor1r2d_dist2parso(par,perp,rpar,rperp,d0,omega_v)
dparo = single(d0 - (d0 - par)./(1 + omega_v.^2./rpar.^2));
dperpo = single(d0 - (d0 - perp)./(1 + omega_v.^2./rperp.^2));
end

function [disoo,danisoo,dratioo,ddeltao,sdanisoo,sddeltao] = my_dtd_pars2dpars(dparo,dperpo)
disoo = single((dparo + 2*dperpo)/3);
danisoo = single((dparo - dperpo)/3);
dratioo = single(msf_notfinite2zero(dparo./dperpo));
ddeltao = single(msf_notfinite2zero(danisoo./disoo));
sdanisoo = single(danisoo.^2);
sddeltao = single(msf_notfinite2zero(sdanisoo./disoo.^2));
end

function [bar_hist] = dist2bar(edges_hist,parameter_data,weights_data)
bar_hist = zeros(size(edges_hist));
for ind_hist = 1:numel(bar_hist)-1
    bar_hist(1,ind_hist) = sum(weights_data(parameter_data >= edges_hist(1,ind_hist) & parameter_data < edges_hist(1,ind_hist+1)));
end
end
