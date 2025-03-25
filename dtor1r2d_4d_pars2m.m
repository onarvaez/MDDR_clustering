function m1 = dtor1r2d_4d_pars2m(dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2,w)

% Get size of the input variables
sz = size(dpar); % Expected: [122, 128, 1, 1000]

% Initialize the full-sized output matrix
m1 = zeros(sz(1), sz(2), sz(3), 10000); % Size: [122, 128, 1, 101]

% Create logical indexing pattern
ind = false(10000,1);
ind(1:10:end) = 1; % Matches original function

% Assign each variable back to its corresponding place in 'm'
m1(:,:,:,circshift(ind,  0,1)) = dpar;
m1(:,:,:,circshift(ind,  1,1)) = dperp;
m1(:,:,:,circshift(ind,  2,1)) = theta;
m1(:,:,:,circshift(ind,  3,1)) = phi;
m1(:,:,:,circshift(ind,  4,1)) = d0;
m1(:,:,:,circshift(ind,  5,1)) = rpar;
m1(:,:,:,circshift(ind,  6,1)) = rperp;
m1(:,:,:,circshift(ind,  7,1)) = r1;
m1(:,:,:,circshift(ind,  8,1)) = r2;
m1(:,:,:,circshift(ind,  9,1)) = w;

end
