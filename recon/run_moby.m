clear all

addpath ../utils;
addpath ./funcs;

% recon parameters
n_iter = 10;
recon_window = [-1,15]; % match with the data file
corrt_window = [-20, -5];
imsize = [62 62 112];
vsize = [0.8008 0.8008 1.6021];
data_folder = 'D:\HPC_backup\MOBY_water';
mu = 0; % detector timing offset in ns

% triples to be reconstructed
lmdata = reshape(touch(fullfile(data_folder, 'final_simple_tlim-1_15/pr/-1.0_15.0.lm'), '*int16'), 5, []);
dt = touch(fullfile(data_folder, 'final_simple_tlim-1_15/pr/w1/-1.0_15.0.mul_fac'), 'float32');
          
% total doubles
lmdata2 = reshape(touch(fullfile(data_folder, 'final_simple_true_doubles/trues.lm'), '*int16'), 5, []);
num_doubles = size(lmdata2, 2);

% count number of randoms in the correction window
s = dir(fullfile(data_folder, 'final_simple_tlim-20_-5/pr/-20.0_-5.0.lm'));
num_randoms = s.bytes / 10;
gamma = num_randoms * diff(recon_window) / diff(corrt_window) / num_doubles;

% get psf kernel
kernel = get_psf_kernel([9 9 13], [0.8 0.8 1.6], vsize); 

% sensitivity
sen = touch('../scanner/uih3_sensimage_wAC_water_62s112-voxelsize_0.8008mm_1.6021-brd0_5-single.smap', 'float32');
sen = psf_blur(sen, kernel, imsize); % blur the sensitivity map

% scanner
tc = reshape(touch('../scanner/NX_nonDOI_xtal_pos_xy','single'), 2, []);
to = touch('../scanner/NX_nonDOI_xtal_pos_z', 'single'); 
tof_info = [250, 39];

%% SIMPLE reconstruction
disp('---Doing MLEM reconstruction...')
[~, at] = mlem(lmdata, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], [], sen); 

disp('---Doing weighted MLEM reconstruction...')
[~, w1] = mlem(lmdata, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], dt, sen);

disp('---Doing MLEM reconstruction...') % can be slow
[~, db] = mlem(lmdata2, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], [], sen);

% division
disp('---Doing lifetime calculation...')
m1 = zeros([imsize, n_iter]);
for ii = 1:n_iter
    m1 = (w1 - gamma * mean(recon_window) * db) ./ (at - gamma * db) - mu;
end

% visualize
m1_slice = reshape(m1(:, 5), [62 62 112]);
m1_slice = squeeze(m1_slice(22, :, :));

figure; imshow_zj([0,50], [0,180], rot90(m1_slice, 1), [0.7025,1.1525]); c2=colorbar;c2.Label.String = 'Lifetime (ns)'; c2.FontSize=18; axis off; 
 

function kernel = initialize_psf(psf_size, psf_fwhm, vsize)
%   psf_size: [x y z] in voxels
%   psf_fwhm: [x y z] in mm
    kernel_xy = get_xy_kernel(psf_size(1:2), psf_fwhm(1:2), vsize(1:2));
    kernel_z  = get_z_kernel(psf_size(3), psf_fwhm(3), vsize(3));
    kernel = kernel_xy .* reshape(kernel_z,1,1,[]);
end

