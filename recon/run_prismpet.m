clear all

addpath ../utils;
addpath ./funcs;

% recon parameters
n_iter = 10;
recon_window = [-1,15]; % match with the data file
corrt_window = [-15, -1];
imsize = [128 128 64];
vsize = [1 1 1];
data_folder = 'D:\Prism-PET\lm_for_SIMPLE';
mu = 0.239; % detector timing offset in ns

% triples to be reconstructed
lmdata = reshape(touch(fullfile(data_folder, 'meat_-1_15_en25000_35000_40000_DC_eID.lm'), '*int16'), 5, []);
dt = touch(fullfile(data_folder, 'meat_-1_15_en25000_35000_40000_DC_eID.mul_fac'), 'float32');
          
% total doubles
lmdata2 = reshape(touch(fullfile(data_folder, 'meat_doubles_eID.lm'), '*int16'), 5, []);
num_doubles = size(lmdata2, 2);

% count number of randoms in the correction window
s = dir(fullfile(data_folder, 'meat_-15_-1_en25000_35000_40000_DC_eID.lm'));
num_randoms = s.bytes / 10; % every event has 10 bytes
gamma = num_randoms * diff(recon_window) / diff(corrt_window) / num_doubles;

% get psf kernel
kernel = get_psf_kernel([11 11 13], [0.8 0.8 0.8], vsize); 

% sensitivity
sen = touch('../scanner/GeoSen_imsize128_64_vsize1_1.smap', 'float32');
sen = psf_blur(sen, kernel, imsize); % blur the sensitivity map

% scanner
tc = reshape(touch('../scanner/prismpet_det_pos_xy', 'single'), 2, []);
to = touch('../scanner/prismpet_det_pos_z', 'single'); 
tof_info = [270, 5];

%% SIMPLE reconstruction
disp('---Doing MLEM reconstruction...')
[~, at] = mlem(lmdata, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], [], sen); 

disp('---Doing weighted MLEM reconstruction...')
[~, w1] = mlem(lmdata, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], dt, sen);

disp('---Doing MLEM reconstruction...')
[~, db] = mlem(lmdata2, n_iter, tc, to, tof_info, imsize, vsize, kernel, [], [], [], sen);

% division
disp('---Doing lifetime calculation...')
m1 = zeros([imsize, n_iter]);
for ii = 1:n_iter
    m1 = (w1 - gamma * mean(recon_window) * db) ./ (at - gamma * db) - mu;
end

% visualize
it = 8; % iteration to show
m1_slice = reshape(m1(:, it), imsize);
m1_slice = squeeze(m1_slice(73, 41:88, 23:45));

figure, imagesc(m1_slice, [0.35,1.4]), c2=colorbar; c2.Label.String = 'Lifetime (ns)'; c2.FontSize=18; axis equal, axis off; 
 



