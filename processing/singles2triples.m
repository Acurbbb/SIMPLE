% This script aims to convert listmode singles data to triples
% Major operations
% 1) do time shifting to simulate lifetime delay
% 2) add energy blur
% 3) add gap to the axial ID

clear all;

addpath ./mex;
tsr = 1; % time scaling, if higher activity is desired 
is_rej = true; % geometric rejection

% path of source file
src_folder = 'D:\HPC_backup\MOBY_water\MOBY_20kBq%cc_0000.0_1800.0_0_0_0\raw_singles'
src_files = dir(fullfile(src_folder, '*_singles.lm'));

% mask of lesion and geometry
lsmask = reshape(touch('mouse_lesion.img', 'int16'), 100, 100, 360); 
lssize = size(lsmask);
vsize = [0.5 0.5 0.5];
sPos = [0, 0, 0]; % shift in source, default is center
origin = -lssize.*vsize/2; % where the first voxel begins

fov = [-25,25;-25,25;-90,90]'; %[x;y;z] in mm

% get a look-up table of crystal position, for geometric rejection
[m_xtal_pos_xy, m_ring_z] = get_crystal_position();

% input data, indexes of the rows, refer to GATE for detailed definition
nrow = 9;      % # of entries of an event
id.evt = 1;    % event ID of singles
id.sposx = 2;  % source position x
id.sposy = 3;  % source position y
id.sposz = 4;  % source position z
id.time = 5;   % time (sec)
id.en = 6;     % energy (MeV)
id.tx = 7;     % transaxial ID
id.ax = 8;     % axial ID
id.cp = 9;     % compton scatter
nTx = 960; % number of transaxial xtals for the NX scanner
minTxDiff_p_a = 20; % triples with prompt gamma and 511 within that amount of transaxial distance in crystals will be rejected
maxTxDiff_p_a = nTx - minTxDiff_p_a; 
minAxDiff_p_a = 10; % triples with prompt gamma and 511 within that amount of axial distance in crystals will be rejected
minTxDiff_a_a = 0;  % triples with two 511s within that amount of transaxial distance in crystals will be rejected
maxTxDiff_a_a = nTx - minTxDiff_a_a;

% time shifting options
opt_time.oPs_bg = 2.5*10^-9; % lifetime of background
opt_time.oPs_le = 2.0*10^-9; % lifetime of lesion
opt_time.frac_Ps = 0.4;      % intensity of positronium formation [0, 1]
opt_time.pPs = 0.125*10^-9;  % 0.125 ns of p-Ps lifetime 
opt_time.direct = 0.4*10^-9; % 0.4 ns of direct annihilation
shift_en = 0.53;             % above which the prompt gamma times are shifted. Note this cut is applied before energy blur

% time windows
pGammaFuture = 20e-9; % prompt gamma time window in the future, corresponding to negative lifetime measurement
pGammaPast = 20e-9;   % prompt gamma time window in the past, corresponding to positive lifetime measurement
coincTime = 1e-9;     % coincidence time window, same as the standard PET 
delayTime = 1.5e-7;   % a delay time 
% delayTime_co_pgr = delayTime / 3; % not used

% energy windows for sorting (applied on blurred energy)
annihEn = [0.43, 0.63];
pGammaEn = [0.7, Inf];

% open output file
output_folder = 'D:\HPC_backup\MOBY_water'; % specify a folder for output
save_dir = fullfile(output_folder, ['tric_SIMPLE_rej_pa_no_aa_tsr', num2str(tsr),'x_en', num2str(annihEn(1)),'-',num2str(annihEn(2)),'_',num2str(pGammaEn(1)),...
                    '_Ah',num2str(coincTime*1e9),'_PG',num2str(-pGammaFuture*1e9),'_',num2str(pGammaPast*1e9),'_dl',num2str(delayTime*1e9)]);

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

mkdir(fullfile(save_dir, 'pr'));
% mkdir(fullfile(save_dir, 'dl_co'));
% mkdir(fullfile(save_dir, 'dl_pg'));
% mkdir(fullfile(save_dir, 'dl_co_pg'));
% mkdir(fullfile(save_dir, 'dl_co_pgr'));


%% triples grouping
tic

% Group the small files into batches. 
% Outer loop over batches
% Inner parallel loop over files in a batch 
max_job = 96; % num of files per batch
n_it = ceil(length(src_files) / max_job); % num of batches
elapsed_prev = 0;
for it = 1:n_it
    
    src_files_temp = src_files((max_job*(it-1)+1):min(max_job*it,length(src_files)));
    disp(['processing ', src_files_temp(1).name,' to ', src_files_temp(end).name,'...']);

    % parallel process a batch of files
    parfor iw = 1:length(src_files_temp)
        
        % read file. Make sure the time order is ascending
        sgl = read_sgl_struct(fullfile(src_files_temp(iw).folder, src_files_temp(iw).name));
        if sum(diff(sgl(1,:))<0) ~= 0
            error('EventID not ascending!');
        end

        % time scaling. Scale the time to simulate a higher activity level.
        if tsr ~= 1
            sgl = scale_time(sgl, tsr, id.evt, id.time);
        end

        % shifting
        sgl = add_gap_axID_NX(sgl, id.ax); % particular to this simulation
        sgl = add_lifetime(sgl, shift_en, id, lsmask, sPos, origin, vsize, opt_time);  
        sgl = add_energy_blur(sgl, 0.13, 0.511, id.en);

        % sort it to be chronological
        [~, sort_idx] = sort(sgl(id.time,:));
        sgl = sgl(:, sort_idx);
                       
        % triples pairing
        tric_pr = pair_triples_takeAllGood(sgl, coincTime, 0, [-pGammaPast,pGammaFuture], 0, annihEn, pGammaEn, minTxDiff_p_a, maxTxDiff_p_a, ...
            minAxDiff_p_a, minTxDiff_a_a, maxTxDiff_a_a, m_xtal_pos_xy, m_ring_z, fov);
        
        % tric_dl_co = pair_triples_takeAllGood(sgl, coincTime, -delayTime, [-pGammaPast,pGammaFuture], 0, annihEn, pGammaEn, minTxDiff_p_a, maxTxDiff_p_a, ...
        %     minAxDiff_p_a, minTxDiff_a_a, maxTxDiff_a_a, m_xtal_pos_xy, m_ring_z, fov);

        % tric_dl_pg = pair_triples_takeAllGood(sgl, coincTime, 0, [-pGammaPast,pGammaFuture], -delayTime, annihEn, pGammaEn, minTxDiff_p_a, maxTxDiff_p_a, ...
        %     minAxDiff_p_a, minTxDiff_a_a, maxTxDiff_a_a, m_xtal_pos_xy, m_ring_z, fov);

        % tric_dl_co_pg = pair_triples_takeAllGood(sgl, coincTime, -delayTime, [-pGammaPast,pGammaFuture], -delayTime, annihEn, pGammaEn, minTxDiff_p_a, maxTxDiff_p_a, ...
        %     minAxDiff_p_a, minTxDiff_a_a, maxTxDiff_a_a, m_xtal_pos_xy, m_ring_z, fov);

        % tric_dl_co_pgr = pair_triples_takeAllGood(sgl, coincTime, -delayTime_co_pgr, [-pGammaPast,pGammaFuture], -delayTime, annihEn, pGammaEn, minTxDiff_p_a, maxTxDiff_p_a, ...
        %     minAxDiff_p_a, minTxDiff_a_a, maxTxDiff_a_a, m_xtal_pos_xy, m_ring_z, fov);

        % travel difference correction
        tric_pr = correct_tof_lifetime(tric_pr, id, m_xtal_pos_xy, m_ring_z);
        % tric_dl_co = correct_tof_lifetime(tric_dl_co, id, m_xtal_pos_xy, m_ring_z);
        % tric_dl_pg = correct_tof_lifetime(tric_dl_pg, id, m_xtal_pos_xy, m_ring_z);
        % tric_dl_co_pg = correct_tof_lifetime(tric_dl_co_pg, id, m_xtal_pos_xy, m_ring_z);
        % tric_dl_co_pgr = correct_tof_lifetime(tric_dl_co_pgr, id, m_xtal_pos_xy, m_ring_z);

        % write
        save_to_lm5_t(tric_pr, id, 0, 0, fullfile(save_dir, 'pr'), src_files_temp(iw).name(1:13), 39.0625*10^-3);
        % save_to_lm5_t(tric_dl_co, id, delayTime, 0, fullfile(save_dir, 'dl_co'), src_files_temp(iw).name(1:13), 39.0625*10^-3);
        % save_to_lm5_t(tric_dl_pg, id, 0, delayTime, fullfile(save_dir, 'dl_pg'), src_files_temp(iw).name(1:13), 39.0625*10^-3);
        % save_to_lm5_t(tric_dl_co_pg, id, delayTime, delayTime, fullfile(save_dir, 'dl_co_pg'), src_files_temp(iw).name(1:13), 39.0625*10^-3);
        % save_to_lm5_t(tric_dl_co_pgr, id, delayTime_co_pgr, delayTime, fullfile(save_dir, 'dl_co_pgr'), src_files_temp(iw).name(1:13), 39.0625*10^-3);
    end

    elapsed = toc; 
    fprintf('Elapsed time of the last iteration is %.6f s\n', elapsed-elapsed_prev);
    elapsed_prev = elapsed;
    if mod(it,2) == 0
        fprintf('Estimated time remained: %.6f min\n', (n_it-it)*elapsed/it/60);
    end
end
fclose all;
toc   

function sgl = scale_time(sgl, tsr, idevt, idt)
% eventID should be ascending
    diff_eID = sgl(idevt,2:end) - sgl(idevt, 1:end-1); 
    id_start = [1, find(diff_eID > 0) + 1]; % find the idx of the first single in a decay
    id_start_with_end = [id_start, size(sgl,2)+1]; % add an end for looping
    
    for ii = 1 : (length(id_start_with_end)-1)      
        t_temp = sgl(idt, id_start_with_end(ii):(id_start_with_end(ii+1)-1));
        
        diff_to_first =  t_temp - t_temp(1);
        
        sgl(idt, id_start_with_end(ii)) = sgl(idt, id_start_with_end(ii)) / tsr; % scale the first single
        sgl(idt, id_start_with_end(ii):(id_start_with_end(ii+1)-1)) = sgl(idt, id_start_with_end(ii)) + diff_to_first; % keep the difference 
    end
end

function in = add_lifetime(in, shift_en, id, lemask, sPos, origin, vsize, opt_time)
    for ii = 1:size(in,2)
        if in(id.en,ii) > shift_en
            in(id.time,ii) = max(0, in(id.time,ii) - generate_time_with_mask(in(:,ii), id, lemask, sPos, origin, vsize, opt_time));
        end
    end
end

function lifetime = generate_time_with_mask(prompt_temp, id, lemask, sPos, origin, vsize, opt_t)
    sidx = ceil((prompt_temp(id.sposx:id.sposz,1)'-sPos-origin)./vsize);
    sidx = min(sidx, size(lemask));
    sidx = max(sidx, 1);
    if lemask(sidx(1),sidx(2),sidx(3)) ~= 0
        lifetime = get_lifetime(opt_t.oPs_le, opt_t.pPs, opt_t.direct, opt_t.frac_Ps);
    else
        lifetime = get_lifetime(opt_t.oPs_bg, opt_t.pPs, opt_t.direct, opt_t.frac_Ps);
    end

end

function time = get_lifetime(t_o, t_p, t_d, frac_ps)
    channel_decider = rand(1);
    if channel_decider <= 1-frac_ps
        time = min(1000*t_d, randexp(t_d));
    elseif channel_decider <= 1-1/4*frac_ps
        time = min(1000*t_o, randexp(t_o));
    else
        time = min(1000*t_p, randexp(t_p));
    end
end

function time = randexp(lt)
    x = rand(1);
    time = -lt*log(-x+1);
end

function singles = add_energy_blur(singles, percent, reference, id_en)
    FWHM_0 = reference*percent;
    FWHMs = FWHM_0 .* sqrt(reference./max(singles(id_en,:),0.2));
    sigmas = FWHMs/2/sqrt(2*log(2));
    % cut the blur range
    bound = 6*sigmas;
    en_blur = min(bound, randn(1,size(singles,2)).*sigmas);
    en_blur = max(-bound, en_blur);
    singles(id_en,:) = singles(id_en, :) + en_blur;
end

function trip = correct_tof_lifetime(trip, id, m_xtal_pos_xy, m_ring_z)
    % travel time correction
    c = 2.998*10^11; % mm/s
    for ii = 1:(size(trip,2)/3)
        trip_temp = trip(:, (3*(ii-1)+1):(3*(ii-1)+3));
        an0 = [m_xtal_pos_xy(:,trip_temp(id.tx, 1)+1); m_ring_z(trip_temp(id.ax, 1)+1)];
        an1 = [m_xtal_pos_xy(:,trip_temp(id.tx, 2)+1); m_ring_z(trip_temp(id.ax, 2)+1)];
        pt = [m_xtal_pos_xy(:,trip_temp(id.tx, 3)+1); m_ring_z(trip_temp(id.ax, 3)+1)];
        
        mid = (an0 + an1)/2;
        
        relative_coin_dist_diff = (-trip_temp(id.time, 1) + trip_temp(id.time, 2)) * c / sqrt(sum((an0-an1).^2));
        an_pos = mid + relative_coin_dist_diff / 2 * (an0 - an1);

        dist_pt = sqrt(sum((pt- an_pos).^2));
        dist_an0 = sqrt(sum((an0- an_pos).^2));
        dist_an1 = sqrt(sum((an1- an_pos).^2));
        
        trip(id.time, 3*(ii-1)+3) = trip_temp(id.time, 3) - (dist_pt/c-0.5*(dist_an0/c+dist_an1/c));
    end
end

function save_to_lm5_t(trip, id, delayCoinc, delayPG, folder, basename, tof_bin)
% lm5 + time measurement in two files
    data = zeros(5, 1/3*size(trip,2), 'int16');
    time = zeros(1, size(data,2), 'single');
    for ii = 1:(size(trip,2)/3)
        trip_temp = trip(:, (3*(ii-1)+1):(3*(ii-1)+3));
        lifetime = 1e9*(0.5*(trip_temp(id.time, 1)+trip_temp(id.time, 2)+delayCoinc) - (trip_temp(id.time, 3)+delayPG)) ;
        tof = 1e9*(trip_temp(id.time, 1) - trip_temp(id.time,2) - delayCoinc);
        data(:, ii) = [int16(trip_temp(id.tx,1)); int16(trip_temp(id.ax,1)); ...
                       int16(trip_temp(id.tx,2)); int16(trip_temp(id.ax,2)); ...
                       int16(tof / tof_bin)]; 
        time(ii) = lifetime;
    end
    dump(data, fullfile(folder, [basename,'.lm']), 'int16');
    dump(time, fullfile(folder, [basename,'.float']), 'float32');
end

function save_etype(trip, id, delayCoinc, delayPG, folder, basename, dtype)
% show the event type in prompt pairing 
% lm6: tx1, ax1, tx2, ax2, tof, lifetime, event type
% 0: true, 1: type I, 2: type II.a, 3: type II.b, 4: type III
    ide = id.evt;
    coin = int8(zeros(1, 1/3*size(trip,2)));
    for ii = 1:(size(trip,2)/3)   
        trip_tmp = trip(:, (3*(ii-1)+1):(3*(ii-1)+3));
        if trip_tmp(ide,1) == trip_tmp(ide,2) && trip_tmp(ide,2) == trip_tmp(ide,3)
            coin(ii) = 0;
        elseif trip_tmp(ide,1) == trip_tmp(ide,2)
            coin(ii) = 1;
        elseif trip_tmp(ide,1) == trip_tmp(ide,3)
            coin(ii) = 2;
        elseif trip_tmp(ide,2) == trip_tmp(ide,3)
            coin(ii) = 3;
        else
            coin(ii) = 4;
        end
    end
    dump(coin, fullfile(folder, [basename,'.etype_int8']), dtype);
end
      
function out = add_gap_axID_NX(in, id_ax)
    in(id_ax, :) = in(id_ax, :) + 1;
    in(id_ax, in(id_ax, :)>24) = in(id_ax, in(id_ax, :)>24) + 2;
    in(id_ax, in(id_ax, :)>50) = in(id_ax, in(id_ax, :)>50) + 2;
    in(id_ax, in(id_ax, :)>76) = in(id_ax, in(id_ax, :)>76) + 2;
    in(id_ax, in(id_ax, :)>102) = in(id_ax, in(id_ax, :)>102) + 2;
    in(id_ax, in(id_ax, :)>128) = in(id_ax, in(id_ax, :)>128) + 2;
    out = in;
end


