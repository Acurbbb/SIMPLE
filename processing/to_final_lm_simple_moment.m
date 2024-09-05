% INPUT: int16 list mode data + float32 time measurements
% OUTPUT: combined int16 list mode data and the correponding time
% measurements in two files. The input events are thresholded by the
% reconstruction time window and the correction time window
clear all; 
fclose all;

tof_bin = 39.0625*10^-3; 

mode = {'recon' ,'correct'};

parent_folder = 'D:\HPC_backup\MOBY_water';
src_folder = 'tric_SPLIT_rej_pa_no_aa_tsr1x_en0.43-0.63_0.7_Ah1_PG-25_85_dl150';

for ii = 1:length(mode)
    if strcmp(mode{ii}, 'recon')
        tmin = -1;
        tmax = 15;
    elseif strcmp(mode{ii}, 'correct')
        tmin = -20;
        tmax = -5;
    end


    save_folder = ['final_simple_tlim', num2str(tmin),'_',num2str(tmax)];

    src_dir_pr = fullfile(parent_folder, src_folder, 'pr');
    % src_dir_dl_pg = fullfile(parent_folder, src_folder, 'dl_pg');
    % src_dir_dl_co = fullfile(parent_folder, src_folder, 'dl_co');
    % src_dir_dl_co_pg = fullfile(parent_folder, src_folder, 'dl_co_pg');
    % src_dir_dl_co_pgr = fullfile(parent_folder, src_folder, 'dl_co_pgr');
    % src_dir_true = fullfile(parent_folder, src_folder, 'true');

    save_dir_pr = fullfile(parent_folder, save_folder, 'pr');  
    % save_dir_dl_pg = fullfile(parent_folder, save_folder, 'dl_pg');
    % save_dir_dl_co = fullfile(parent_folder, save_folder, 'dl_co');
    % save_dir_dl_co_pg = fullfile(parent_folder, save_folder, 'dl_co_pg');
    % save_dir_dl_co_pgr = fullfile(parent_folder, save_folder, 'dl_co_pgr');
    % save_dir_true = fullfile(parent_folder, save_folder, 'true'); 

    mkdir_ws(save_dir_pr);
    % mkdir_ws(save_dir_dl_pg);
    % mkdir_ws(save_dir_dl_co);
    % mkdir_ws(save_dir_dl_co_pg);
    % mkdir_ws(save_dir_dl_co_pgr);

    tic
    consolidate_and_save(src_dir_pr, save_dir_pr, tmin, tmax);
    toc
    % consolidate_and_save(src_dir_dl_pg, save_dir_dl_pg, 1, tmin, tmax);
    % toc
    % consolidate_and_save(src_dir_dl_co, save_dir_dl_co, 1, tmin, tmax);
    % toc
    % consolidate_and_save(src_dir_dl_co_pg, save_dir_dl_co_pg, 1, tmin, tmax);
    % toc
    % consolidate_and_save(src_dir_dl_co_pgr, save_dir_dl_co_pgr, 1, tmin, tmax);
    % toc
    % consolidate_and_save(src_dir_true, save_dir_true, fullfile(save_dir_true,'w1'), 1, tmin, tmax);
    % toc
end



function mkdir_ws(folder)
mkdir(fullfile(folder, 'w1'));
end

function consolidate_and_save(src_dir, save_dir_lm, tmin, tmax)
% input: lm5 and float lifetime, output: concolidated and thresholded data

outLMFile = fullfile(save_dir_lm, sprintf('%04.1f_%04.1f.lm',tmin,tmax));
outTimeFile_w1 = fullfile(save_dir_lm, 'w1', sprintf('%04.1f_%04.1f.mul_fac',tmin,tmax));
fid1 = fopen(outLMFile, 'wb');
fid2 = fopen(outTimeFile_w1, 'wb'); 

src_files = dir(fullfile(src_dir, '*lm'));
for fi = 1:length(src_files)
    % read data
    data = int16(reshape(touch(fullfile(src_files(fi).folder, src_files(fi).name), '*int16'), 5, []));
    time = single(reshape(touch(fullfile(src_files(fi).folder, strrep(src_files(fi).name, '.lm', '.float')), '*single'), 1, []));
    
    % do thresholding
    data = data(:, time > tmin & time < tmax);
    time = time(time > tmin & time < tmax);
    
    fwrite(fid1, data, 'int16');
    fwrite(fid2, time, 'single');
end

fclose(fid1); fclose(fid2); 
end


