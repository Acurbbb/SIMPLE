function [x, xs] = mlem(lm, n_iter, tc, to, tof_info, imsize, vsize, kernel, x, af, mf, simg) 

    if nargin < 12 || isempty(simg)
        simg = ones(imsize);
    end

    if nargin < 11 || isempty(mf)
        mf = ones(size(lm,2), 1);
    end

    if nargin < 10 || isempty(af)
        af = zeros(size(lm,2), 1);
    end

    if nargin < 9 || isempty(x)
        x = 1e-5 * ones(prod(imsize), 1); % initialize with small values
    end
    
    if nargin < 8 || isempty(kernel)
        kernel = 1;
    end

    if nargin < 7 
        error("Not enough input!");
    end

    % list-mode EM recon          
    if nargout > 1
        xs = zeros([prod(imsize), n_iter]);
    end

    for it = 1:n_iter    
        xb = psf_blur(x, kernel, imsize);
        px = fproj_tof_mt(xb, imsize, vsize, tc, to, lm, tof_info) + af;

        xem = bproj_tof_mt(mf(:) ./ (px(:)+1e-10), imsize, vsize, tc, to, lm, tof_info);
        xem = psf_blur(xem, kernel, imsize);

        x(simg>0) = x(simg>0) .* xem(simg>0) ./ simg(simg>0);
        x(simg<=0) = 0;

        if nargout > 1
            xs(:,it) = x;
        end

        disp(['   OSEM iteration ', num2str(it),' finished...']);
    end
end

