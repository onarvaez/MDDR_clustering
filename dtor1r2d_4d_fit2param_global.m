function dps = dtor1r2d_4d_fit2param_global(mfs_fn, dps_fn, opt,norm_clust)
% function dps = dtor1r2d_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtor1r2d_opt(opt);
%dps = mdm_mfs_load(mfs_fn);
%k_clust = (norm_clust);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end

% create parameter maps and save them
%mask = dps.m;
m = dps.m;
dps = rmfield(dps,'m');

omega = opt.dtor1r2d.maps_omega;
Nomega = numel(omega);
for nomega = 1:Nomega
    %[dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtor1r2d_4d_m2pars(m);
    [dpar,dperp,theta,phi,r1,r2,w] = dtor1r2d_4d_m2parso(m,omega(nomega));

    sz = size(m);
    nn = size(dpar,4);

    %Calculate derived parameters
    [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
    [diso,daniso,dratio,ddelta,sdaniso,sddelta] = dtd_pars2dpars(dpar,dperp);

    dtor1r2ds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
        'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r1',r1,'r2',r2);

    dps = dtor1r2d_dtor1r2ds2dps(dps, dtor1r2ds);

    % reshape help functions
    sz_reshape  = msf_size(m(:,:,:,1), 3);
    g_reshape = @(a,n) reshape(a, prod(sz_reshape(1:3)), n);
    f_reshape = @(a,n) reshape(a, sz_reshape(1), sz_reshape(2), sz_reshape(3), n);
    dt = cat(4,dps.mdxx,dps.mdyy,dps.mdzz,dps.mdxy,dps.mdxz,dps.mdyz);
    dps = tm_dt_to_dps(g_reshape(dt, 6)*1e9, dps, f_reshape, 0.0001);
%     [min(diso(:)) max(diso(:))]
%     [min(sddelta(:))  max(sddelta(:))]


    %Per-cluster statistical measures
    for nbin = 1:numel(norm_clust)   

        load(norm_clust{nbin});
          dps_bin.no = nbin;
          
          
%         diso_pre = diso .* norm_k_mclu;
%         dratio_pre = dratio .* norm_k_mclu;
%         sddelta_pre = sddelta .* norm_k_mclu;
%         r1_pre = r1 .* norm_k_mclu;
%         r2_pre = r2 .* norm_k_mclu;
%         
        dtor1r2ds_temp = dtor1r2ds;
        
%         dtor1r2ds_temp.diso = diso_pre;
%         dtor1r2ds_temp.dratio = dratio_pre;
%         dtor1r2ds_temp.sddelta = sddelta_pre;
%         dtor1r2ds_temp.r1 = r1_pre;
%         dtor1r2ds_temp.r2 = r2_pre;
%         
        dtor1r2ds_temp.diso = diso;
        dtor1r2ds_temp.dratio = dratio;
        dtor1r2ds_temp.sddelta = sddelta;
        dtor1r2ds_temp.r1 = r1;
        dtor1r2ds_temp.r2 = r2;
        
        
        dtor1r2ds_temp.w = dtor1r2ds.w;
        dps_bin = dtor1r2d_dtor1r2ds2dps(dps_bin, dtor1r2ds_temp);
        dps_bin.f = dps_bin.s0./dps.s0;
        dps.bin{nbin} = dps_bin;
    end
        dps.omega{nomega} = dps;

end


%Rate of change with nu = omega/(2*pi);
fields = {'mdiso'; 'msddelta'; 'vdiso'; 'vsddelta'; 'cvdisosddelta'};
for nfield = 1:numel(fields)
    field = fields{nfield};
    parnam = ['d' field 'dnu'];
    dps.(parnam) = (dps.omega{end}.(field) - dps.omega{1}.(field))./(omega(end)-omega(1))*2*pi;
    parnam2 = ['nd' field 'dnu'];
    dps.(parnam2) = (dps.omega{end}.(field) - dps.omega{1}.(field))./(dps.omega{end}.(field) + dps.omega{1}.(field))*2./(omega(end)-omega(1))*2*pi;
end

fields = {'mdiso'; 'msddelta'};
for nbin = 1:numel(norm_clust)   
    for nfield = 1:numel(fields)
        field = fields{nfield};
        parnam = ['d' field 'dnu'];
        dps.bin{nbin}.(parnam) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(omega(end)-omega(1))*2*pi;
        parnam2 = ['nd' field 'dnu'];
        dps.bin{nbin}.(parnam2) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(dps.omega{end}.bin{nbin}.(field) + dps.omega{1}.bin{nbin}.(field))*2./(omega(end)-omega(1))*2*pi;
    end
end

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end
