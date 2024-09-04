function [m_xtal_pos_xy, m_ring_z] = get_crystal_position()
    num_of_rings = 156; % m_sysparms.crystal_array_size.a * m_sysparms.number_of_detector_modules.a
    ring_pitch =  3.2043 + 0; % m_sysparms.crystal_size.axial + m_sysparms.crystal_gap_size.axial
    num_of_blocks_t = 20; % number_of_detector_modules
    num_of_xtals_d_zmx = 1; % crystal_array_size.d
    num_of_xtals_t = 48; % crystal_array_size.t
    ring_diameter = 498.6; % detector_ring_diameter
    xtal_pitch_t = 1.6017 + 0; % m_sysparms.crystal_size.transaxial_front + m_sysparms.crystal_gap_size.transaxial_front
    xtal_length = 18.1; % crystal_size.depth 
    img_size_ij = 62; % image_size.i 
    img_size_k = 112; % image_size.k
    vox_size_ij = 0.8008; % voxel_size_i
    vox_size_k = 1.6021; % voxel_size_k
    tw_resolution = 250; % tof_info.resolution 
    tw_spacing = 39;% tof_info.bin_size
    half_ang = false;


    xb = num_of_xtals_d_zmx*xtal_length/2;
    yb = num_of_xtals_t*xtal_pitch_t/2;
    zb = num_of_rings*ring_pitch/2;
    blockangle = 2*pi/num_of_blocks_t;    
    ringradius = ring_diameter/2;

    ang_step = 360.0 / num_of_blocks_t;
    if half_ang
        ang_init = 0.5 * ang_step;
    else
        ang_init = 0;
    end

    m_xtal_pos_xy = zeros(2, num_of_xtals_d_zmx * num_of_blocks_t * num_of_xtals_t); 

    for b = 0:(num_of_xtals_d_zmx-1)

        zmx_xl = b*xtal_length*2;

        for a = 0:(num_of_blocks_t-1)

            angle = ang_init + ang_step * a;
            bx = (ring_diameter + xtal_length + zmx_xl) * 0.5 * cos(angle * pi / 180.0);
            by = (ring_diameter + xtal_length + zmx_xl) * 0.5 * sin(angle * pi / 180.0);

            for i = 0:(num_of_xtals_t-1)

                xtalxy_temp = zeros(2,1);
                x = xtal_pitch_t * (-num_of_xtals_t * 0.5 + 0.5 + i);
                xtalxy_temp(1) = bx + x * cos((angle + 90.0) * pi / 180.0);
                xtalxy_temp(2) = by + x * sin((angle + 90.0) * pi / 180.0);
                m_xtal_pos_xy(:, b*num_of_blocks_t*num_of_xtals_t+a*num_of_xtals_t+i+1) = xtalxy_temp;

            end
        end
    end

    m_ring_z = zeros(1, num_of_rings);
    for i = 0 : (num_of_rings-1)
        m_ring_z(1, i+1) = (-num_of_rings * 0.5 + i + 0.5)*ring_pitch;
    end