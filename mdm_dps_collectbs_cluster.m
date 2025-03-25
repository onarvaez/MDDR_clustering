function bs_dps = mdm_dps_collectbs_cluster(method, bs_path, opt, norm_clust,Nclust)
%function bs_dps = mdm_dps_collectbs_cluster(method, bs_path, opt, norm_clust,slice_select,Nclust)



bsno = msf_getdirno(bs_path);        
bs_dps = cell(numel(bsno),1);
load_mfs_success = ones(numel(bsno),1);
%parfor nbs = 1:numel(bsno)
for nbs = 1:numel(bsno)
    cls_mask= norm_clust{1,Nclust}{nbs};
    mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
    if exist(mfs_fn,'file')==2

        try
            
            mfs = mdm_mfs_load(mfs_fn);
            m = double(mfs.m);

            %bs_dps{nbs} = feval([method '_4d_fit2param_cluster'], m, [], opt,cls_mask,slice_select);
            bs_dps{nbs} = feval([method '_4d_fit2param_cluster'], m, [], opt,cls_mask);


%             if strcmp(method,'dtr2d')
%                 bs_dps{nbs} = dtr2d_4d_fit2param(m, [], opt);
%             elseif strcmp(method,'dtr1d')
%                 bs_dps{nbs} = dtr1d_4d_fit2param(m, [], opt);
%             elseif strcmp(method,'dtr1r2d')
%                 bs_dps{nbs} = dtr1r2d_4d_fit2param(m, [], opt);
%             elseif strcmp(method,'dtd')
%                 bs_dps{nbs} = dtd_4d_fit2param(m, [], opt);
%             end

            bs_dps{nbs}.nii_h = mfs.nii_h;
        
        catch
            load_mfs_success(nbs) = 0;
            warning(['mdm_dps_collectbs failed loading ' mfs_fn])
        end
    end
end


bs_dps = bs_dps(logical(load_mfs_success));



