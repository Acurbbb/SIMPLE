function kernel = get_psf_kernel(psf_size, psf_fwhm, vsize)
%   psf_size: [x y z] in voxels
%   psf_fwhm: [x y z] in mm
%   vsize: [x, y, z] in mm
kernel_xy = get_xy_kernel(psf_size(1:2), psf_fwhm(1:2), vsize(1:2));
kernel_z  = get_z_kernel(psf_size(3), psf_fwhm(3), vsize(3));
kernel = kernel_xy .* reshape(kernel_z,1,1,[]);
end


function kernel_xy = get_xy_kernel(ksize, FWHM, voxsize)
%   ksize: [x, y], kernel size in voxel
%   FWHM: FWHM of kernel in mm
%   voxsize: voxel size in mm

sig_x = FWHM(1)/2/sqrt(2*log(2));
sig_y = FWHM(2)/2/sqrt(2*log(2));

x = (-(ksize(1)-1)/2:(ksize(1)-1)/2) * voxsize(1);
y = (-(ksize(2)-1)/2:(ksize(2)-1)/2) * voxsize(2);

[X, Y] = meshgrid(x, y);

kernel_xy = exp( -( X.^2 / sig_x^2 + Y.^2 / sig_y^2 )/2 );

kernel_xy = kernel_xy/sum(kernel_xy, 'all');
end


function kernel_z = get_z_kernel(ksize, FWHM, voxsize)
%   ksize: scalar, kernel size in voxel
%   FWHM: FWHM of kernel in mm
%   voxsize: voxel size in mm

sig = FWHM/2/sqrt(2*log(2));

z = (-(ksize-1)/2:(ksize-1)/2) * voxsize;

kernel_z = exp(-(z.^2)/2/sig^2);
kernel_z = kernel_z/sum(kernel_z, 'all');
end
