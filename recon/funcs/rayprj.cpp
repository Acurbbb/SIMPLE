#include "mex.h"
#include "rtracer.h"
#include <sys/time.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs == 0)
        return;     

#if USE_FPROJ
#if USE_BPROJ
    mexErrMsgTxt("USE_FPROJ and USE_BPROJ can't be together!");
    return;
#endif
#endif 
    
#if USE_BPROJ
#if USE_FPROJ
    mexErrMsgTxt("USE_FPROJ and USE_BPROJ can't be together!");
    return;
#endif
#endif
    
#if USE_FPROJ

    // forward projector:
    //
    //  proj = rayfp(image, imgsize, voxsize, xtal_xy_position, ring_z_location, lmdata, tof_info)
    //
    // backprojector:
    //
    //  image = raybp(prjs, imgsize, voxsize, xtal_xy_position, ring_z_location, lmdata, tof_info)
    //
    //

    // read in image content
    IMAGE_DATA_TYPE* image = 0;
    if (nrhs < 1) {
        mexErrMsgTxt("need an image!");
        return;
    } else {
        image = (IMAGE_DATA_TYPE*)mxGetPr(prhs[0]);
        if (image == NULL) {
            mexErrMsgTxt("invalid image!");
            return;
        }        
    }
    
    // size based on image content
    INT n0 = mxGetN(prhs[0]);
    INT n1 = mxGetM(prhs[0]);
    INT nv = n0 * n1;

#endif
    
#if USE_BPROJ
    
    REAL* proj = 0;
    if (nrhs < 1)  {
        mexErrMsgTxt("need projections!");
        return;
    } else {
        proj = mxGetPr(prhs[0]);
        if (proj == NULL) {
            //mexWarnMsgTxt("no proj found!");
        }
    }

    INT num_of_projs = mxGetN(prhs[0]) * mxGetM(prhs[0]);
#endif
    
    // read in image dimension
    REAL* imgsize = 0; 
    if (nrhs < 2) {
        mexErrMsgTxt("not specify image dimension!");
        return;
    } else {
        imgsize = mxGetPr(prhs[1]);
        if ((imgsize == NULL) ||
            (imgsize[0] <= 0) || 
            (imgsize[1] <= 0) || 
            (imgsize[2] <= 0)) {
            mexErrMsgTxt("invalid image size!");
            return;
        }
    }     

    // convert to integers
    INT ni = INT(imgsize[0]);
    INT nj = INT(imgsize[1]);
    INT nk = INT(imgsize[2]);    
    INT nsize = ni * nj * nk;

#if USE_FPROJ    
    // compare to image vector size
    if (nsize != nv) {
        mexErrMsgTxt("invalid image vector or image size!");
        return;
    }
#endif

    // read in voxel size
    REAL* voxsize = 0;
    if (nrhs < 3) {
        mexErrMsgTxt("not specify voxel size!");    
        return;
    } else {
        voxsize = mxGetPr(prhs[2]);
        if ((voxsize == NULL) || 
            (voxsize[0] <= 0) || 
            (voxsize[1] <= 0) || 
            (voxsize[2] <= 0)) {
            mexErrMsgTxt("invalid voxel size!");
            return;
        }
    }
             
    // xtal transaxial position (2d) 
    XTAL_XY_POSITION* xtal_xy_positions = 0;
    if (nrhs < 4) {
        mexErrMsgTxt("not specify crystal transaxial coordinates!");
        return;
    } else {
        xtal_xy_positions = (XTAL_XY_POSITION*)mxGetPr(prhs[3]);
        if ((xtal_xy_positions == NULL) || 
            (mxGetM(prhs[3]) != 2)) {
            mexErrMsgTxt("invalid crystal positions!");
            return;
        }
    }
    INT num_of_xtals = mxGetN(prhs[3]);
        
    // ring z locations
    REAL* ring_z_locations = 0;
    if (nrhs < 5) {
        mexErrMsgTxt("not specify ring axial coordinates!");
        return;
    } else {
        ring_z_locations = mxGetPr(prhs[4]);
        if ((ring_z_locations == NULL)) {
            mexErrMsgTxt("invalid ring axial coordinates!");
            return;
        }
    }
    INT nring = mxGetM(prhs[4]) * mxGetN(prhs[4]);
    if (nring == 0) {
        mexErrMsgTxt("invalid ring axial coordinates!");
        return;
    }

    if (nrhs < 6) {
        mexErrMsgTxt("no listmode data found!");
        return;
    } else {
        
        if ((mxGetClassID(prhs[5]) != mxINT16_CLASS) ||
            (mxGetPr(prhs[5]) == NULL)) {
            mexErrMsgTxt("invalid listmode data, must be INT16!");
            return;
        }
        
        if (mxGetM(prhs[5]) != 5) {
            mexErrMsgTxt("invalid listmode format, must be 5 x num_of_prompts!");
            return;
        }
    }
    LMEVENT* lmdata = (LMEVENT*)mxGetPr(prhs[5]);
    INT num_of_prompts = mxGetN(prhs[5]);

#if USE_BPROJ
    
    if ((num_of_prompts != num_of_projs) && (proj != 0)) {
        mexErrMsgTxt("wrong number of projections or prompts!");
        return;
    }

#endif
    
    REAL tw_resolution = -1;
    REAL tw_spacing = -1;

#if USE_TOF

    if (nrhs < 7) {
        mexErrMsgTxt("not specify timing resoluton and window spacing!");
        return;
    } else {
        REAL* twinfo = mxGetPr(prhs[6]);
        if ((twinfo == 0) || (mxGetN(prhs[6]) * mxGetM(prhs[6]) < 2)) {
            mexErrMsgTxt("invalid timing info!");
            return;
        }
        
        tw_resolution = twinfo[0];
        tw_spacing = twinfo[1];
    }
#endif

    // show something
#if USE_FPROJ
#if DEBUG
    mexPrintf("forward projecting %d prompts ... ", num_of_prompts);
#endif    
#endif
#if USE_BPROJ
#if DEBUG
    mexPrintf("backprojecting %d prompts ... ", num_of_prompts);
#endif    
#endif
#if DEBUG
    mexPrintf("\tImage size: %d,%d,%d\n"
              "\tVoxel size: %g,%g,%g\n"
              "\tNumber of crystals on a ring: %d\n"
              "\tNumber of rings: %d\n"
              "\tNumber of prompts: %d\n"
              "\tTOF info: %g, %g\n",
              ni, nj, nk, voxsize[0], voxsize[1], voxsize[2],
              num_of_xtals, nring, num_of_prompts, tw_resolution, tw_spacing);
#endif

    // create raytracer
    ImageRayTracer::NUMBER_OF_SIGMA = 3;
    ImageRayTracer::GWSAMPLESIZE = 10240;
    ImageRayTracer raytracer(ni, nj, nk, 
                             voxsize[0], voxsize[1], voxsize[2],
                             tw_resolution * 0.15 / (2*sqrt(2*log(2))),
                             tw_spacing * 0.15);


	// create a timer
	timeval hres_t0;
	timeval hres_t1;
	gettimeofday(&hres_t0, NULL);

#if USE_FPROJ    
    // create projection buffer
    plhs[0] = mxCreateDoubleMatrix(num_of_prompts, 1, mxREAL);
    REAL* proj = mxGetPr(plhs[0]);
    
    INT n;
    LMEVENT e;
#if USE_OMP
    #pragma omp parallel for private(n, e)
#endif
    for (n = 0; n < num_of_prompts; n ++) {
        
        e = lmdata[n];
#if USE_BRESENHAM
        proj[n] = raytracer.fprojBresenham(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        image, e.tbin_id);
#else        

#if USE_LINEAR_INTERP
        proj[n] = raytracer.fprojLinterp(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        image, e.tbin_id);

#else
        proj[n] = raytracer.fproj(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        image, e.tbin_id);

#endif
#endif                        
    }
#endif
    
#if USE_BPROJ
    plhs[0] = mxCreateDoubleMatrix(nsize, 1, mxREAL);
    REAL* bpimg = mxGetPr(plhs[0]);
    
    INT n;
#if USE_OMP    
    #pragma omp parallel for private(n)
#endif    
    for (n = 0; n < num_of_prompts; n ++) {
        
        LMEVENT e = lmdata[n];
        REAL weight = (proj == 0) ? 1.0 : proj[n]; 
//        if (weight > 0) {
#if USE_SQ_LOR   
        raytracer.bproj_sq(xtal_xy_positions[e.xtal_id_1].x,
                            xtal_xy_positions[e.xtal_id_1].y,
                            ring_z_locations[e.ring_id_1],
                            xtal_xy_positions[e.xtal_id_2].x,
                            xtal_xy_positions[e.xtal_id_2].y,
                            ring_z_locations[e.ring_id_2], 
                            weight, bpimg, e.tbin_id);
#else     

#if USE_BRESENHAM
#if 0
        mexPrintf("#%d: [%f %f] [%f %f]\n", n + 1, 
                    xtal_xy_positions[e.xtal_id_1].x,
                    xtal_xy_positions[e.xtal_id_1].y,
                    xtal_xy_positions[e.xtal_id_2].x,
                    xtal_xy_positions[e.xtal_id_2].y);
#endif			
        raytracer.bprojBresenham(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        weight, bpimg, e.tbin_id);
#else

#if USE_LINEAR_INTERP // linear interp
#if 0
        mexPrintf("#%d: [%f %f] [%f %f]\n", n + 1, 
                    xtal_xy_positions[e.xtal_id_1].x,
                    xtal_xy_positions[e.xtal_id_1].y,
                    xtal_xy_positions[e.xtal_id_2].x,
                    xtal_xy_positions[e.xtal_id_2].y);
#endif			
        raytracer.bprojLinterp(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        weight, bpimg, e.tbin_id);
        
#else
        raytracer.bproj(xtal_xy_positions[e.xtal_id_1].x,
                        xtal_xy_positions[e.xtal_id_1].y,
                        ring_z_locations[e.ring_id_1],
                        xtal_xy_positions[e.xtal_id_2].x,
                        xtal_xy_positions[e.xtal_id_2].y,
                        ring_z_locations[e.ring_id_2], 
                        weight, bpimg, e.tbin_id);
#endif

#endif
                            
#endif
//        }
    }
    
#endif

    
    gettimeofday(&hres_t1, NULL);
    double et = (hres_t1.tv_sec - hres_t0.tv_sec) * 1000.0;
    et += (hres_t1.tv_usec - hres_t0.tv_usec) / 1000.0;
#if DEBUG    
    mexPrintf("OK, done! Elapsed time is %.2lf ms.\n", et);
#endif

    return;
    
}
