% Function to process MEPs and rename trials based on VR txt files.
% It also extracts MEPs amplitudes and latencies for each block.
% Developed for Lieve Filbrich experiment : using virtual reality, TMS
% during sensitization using capsaicin.
%
% Folder and subject naming must be:
% i.e.: sub-00$_group-ctrl or sub-00$_group-caps
% 
% The main_folder must contain:
%   - raw_data folder (must contain all participants data folders),
%       - the subject folder must contain:
%               - 'MRI' folder
%               - 'Sessions' folder
%               - the files from unity (.txt, .meta, .csv)
%
%   - pre-processed data folder (created by the function, to save letswave data files)
%   - results folder  (created by the function, to save .mat MEP structure and .csv files)
%   - scripts_toolboxes folder (with function and toolboxes)
% 
% BLOCK 1 (BLK1) is the PRE capsaicin block
% BLOCKS 2 to 5 (BLK2 to BLK5) are the POST capsaicin blocks
% 
%
% Required toolboxes:
%   natsortfiles (2022, Stephen Cobeldick)
%   letswave6
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Cédric Lenoir, NeTMeD, IoNS, UCLouvain, February 2026

function mep_process

%% 1) LOAD DATA AND SET A FEW PARAMETERS

% set the paths and initialize letswave
% Select the general data folder which contains all subjects folders
main_folder = uigetdir('C:\','Select the main folder');
cd(main_folder)
raw_data_folder = fullfile(main_folder,'raw_data');
results_folder = fullfile(main_folder,'results');
pre_proc_folder = fullfile(main_folder,'pre-processed');

if ~isfolder(results_folder); mkdir(results_folder); end
if ~isfolder(pre_proc_folder); mkdir(pre_proc_folder); end

addpath(genpath(main_folder))
close all hidden

% initialize letswave 6
letswave()
clc

% dialog box
prompt = {'\fontsize{12} Subject ID (3 digits, i.e.,:001)? :','\fontsize{12} Group? 0->CONTROL / 1->CAPSAICIN :',...
    '\fontsize{12} Sensitized arm (L/R): ','\fontsize{12} Stimulated hemisphere (L/R):',...
    '\fontsize{12} PRE-CAPSAICIN BLOCK1 TMS start index? : ','\fontsize{12} PRE-CAPSAICIN BLOCK1 TMS stop index? : ',...
    '\fontsize{12} POST-CAPSAICIN BLOCK2 TMS start index? : ','\fontsize{12} POST-CAPSAICIN BLOCK2 TMS stop index? : ',...
    '\fontsize{12} POST-CAPSAICIN BLOCK3 TMS start index? : ','\fontsize{12} POST-CAPSAICIN BLOCK3 TMS stop index? : ',...
    '\fontsize{12} POST-CAPSAICIN BLOCK4 TMS start index? : ','\fontsize{12} POST-CAPSAICIN BLOCK4 TMS stop index? : ',...
    '\fontsize{12} POST-CAPSAICIN BLOCK5 TMS start index? : ','\fontsize{12} POST-CAPSAICIN BLOCK5 TMS stop index? : ',...
    '\fontsize{12} Baseline RMS threshold (µV): ','\fontsize{12} Index of TMS event to exclude (comma separated): ',};
dlgtitle = 'LOAD SUBJECT DATA';
opts.Interpreter = 'tex';
dims = repmat([1 80],16,1);
definput = {'001','1','R','L','114','137','138','164','165','192','193','220','221','246','15',''};
info = inputdlg(prompt,dlgtitle,dims,definput,opts);
% make sure that subject ID number is 3 digits
if size(char(info(1)),2) == 3
    subject_id = char(info(1));
else
    subject_id = sprintf('%.3d',str2double(info(1)));
end
group_cond = str2double(info(2));
if group_cond == 0
    group_str = 'ctrl';
elseif group_cond == 1
    group_str = 'caps';
else
    disp('Group is not correct! try again')
end
sensi_arm = char(info(3));
stim_hemi = char(info(4));
ixd_tms_pre_start = str2double(info(5));
ixd_tms_pre_stop = str2double(info(6));
ixd_tms_pst_start = str2double([info(7) info(9) info(11) info(13)]);
ixd_tms_pst_stop = str2double([info(8) info(10) info(12) info(14)]);
abs_thrshld = str2double(info(15));
tms_bad = char(info(16));
if ~isempty(tms_bad)
    tms_bad = regexp(tms_bad, '\d+', 'match');
    for ibad = 1:size(tms_bad,2)
        tms_idx_bad(ibad,1) = str2double(tms_bad{1,ibad});
    end
else
    tms_idx_bad = [];
end

block_num = 5;

% check if the number of blocks is 5 and that we have TMS indices for each block
if sum(isnan([ixd_tms_pre_start ixd_tms_pre_stop ixd_tms_pst_start ixd_tms_pst_stop])) ~= 0
    disp('At least 1 TMS index or BLOCK is missing ! Check.')
    close all hidden
    return
else
end

% store info in structure
sub_info = struct;
sub_info.sub_ID = subject_id;
sub_info.group = group_str;
sub_info.sensitized_arm = sensi_arm;
sub_info.stimulated_hemisph = stim_hemi;
sub_info.TMS_triggers.pre = [ixd_tms_pre_start;ixd_tms_pre_stop];
sub_info.TMS_triggers.pst = [ixd_tms_pst_start;ixd_tms_pst_stop];
sub_info.bad_TMS = tms_bad;

% folder name of Visor EMG data
gen_subject_folder = fullfile(raw_data_folder,strcat('sub-',subject_id,'_group-',group_str));
subject_folder = fullfile(raw_data_folder,strcat('sub-',subject_id,'_group-',group_str),'Sessions');
sub_pre_proc_folder = fullfile(pre_proc_folder,sprintf('sub-%s_group-%s',subject_id,group_str));
sub_results_folder = fullfile(results_folder,sprintf('sub-%s_group-%s',subject_id,group_str));
if ~isfolder(sub_pre_proc_folder); mkdir(sub_pre_proc_folder); end
if ~isfolder(sub_results_folder); mkdir(sub_results_folder); end

% list sessions
session_list = dir(subject_folder);
session_list = session_list(~ismember({session_list.name}, {'.', '..'}));

% list EMG CNT files
file_list = dir(fullfile(session_list(end).folder,session_list(end).name,'*emg.cnt*'));

% import CNT file
[out_data,~] = RLW_import_CNT(fullfile(session_list.folder,session_list.name,file_list.name));
filename = file_list.name(1:end-8);
savename = strcat('sub-',subject_id,'_group-',group_str,'_',filename);
out_data.header.name = savename;
CLW_save(sub_pre_proc_folder,out_data.header,out_data.data);
clear out_data
[header, data] = CLW_load(fullfile(sub_pre_proc_folder,savename));
sr = 1/header.xstep;


%% 2) Pre-processing steps in letswave

% DC removal
[dc_header, dc_data] = RLW_dc_removal(header,data,'linear_detrend',0);

% high pass filter Butterworth 4 Hz; order 4
low_cutoff = 4;
order = 2;
[filt_header, filt_data] = RLW_butterworth_filter(dc_header,dc_data,'filter_type','highpass','low_cutoff',low_cutoff,'filter_order',order);

% segmentation
xstart = -0.2;
xduration = 0.7;
[seg_header, seg_data] = RLW_segmentation(filt_header, filt_data, {'1'},'x_start',xstart,'x_duration',xduration);
seg_header.chanlocs.labels = 'EMG1';
seg_header.name = strcat(header.name,' DC HPfilt ep');

% save dataset
CLW_save(sub_pre_proc_folder,seg_header, seg_data);


%% 3) Remove irrelevant TMS events

% keep TMS events after rMT and remove bad TMS events if there are any
% arrange epochs of valid TMS triggers for each TMS block
valid_tms_idx = cell(block_num,1);
for iblock = 1:block_num
    if iblock == 1
        valid_tms_idx{iblock,1} = setdiff((ixd_tms_pre_start:ixd_tms_pre_stop),tms_idx_bad);
    else
        valid_tms_idx{iblock,1} = setdiff((ixd_tms_pst_start(iblock-1):ixd_tms_pst_stop(iblock-1)),tms_idx_bad);
    end
end

% save valid epochs in separate lw files for each block
for iblock = 1:block_num
    [block_header{iblock,1}, block_data{iblock,1}] = RLW_arrange_epochs(seg_header, seg_data,valid_tms_idx{iblock,1});
    block_header{iblock,1}.name = strcat(seg_header.name,' BLK ',num2str(iblock));
    CLW_save(sub_pre_proc_folder,block_header{iblock,1}, block_data{iblock,1});
end


%% 4) Sort MEPs according to trial condition (from VR files) per block and save lw file for each event

event_labels ={'base','sensi','ctrl'};
% load the text files from the VR software
% get the names of the files and get the date and time
txt_files = dir(fullfile(gen_subject_folder,'*TMSSensitized*.csv'));
for ifile = 1:size(txt_files,1)
    idx_date = strfind(txt_files(ifile).name,'TMSSensitized_');
    txt_date_files{ifile,1} = txt_files(ifile).name(idx_date:end);
end
% sort files by date and time
[~, idx, ~] = natsort(txt_date_files,'\d+');

% read the .csv file to get the events codes for BLK 1 to 5
for iblock = 1:block_num
    temp = readtable(fullfile(gen_subject_folder,txt_files(idx(iblock)).name));
    for itrial = 1:size(temp,1)
        idx_base(itrial,1) = strfind(temp.TYPE(itrial),'BASE');
        idx_sensi(itrial,1) = strfind(temp.TYPE(itrial),'WASPsensi');
        idx_ctrl(itrial,1) = strfind(temp.TYPE(itrial),'WASPcontrol');
    end

    % rename the events
    for itrial = 1:size(temp,1)
        if ~isempty(idx_base{itrial,1})
            block_header{iblock,1}.events(itrial).code = event_labels{1,1};
        else
        end
        if ~isempty(idx_sensi{itrial,1})
            block_header{iblock,1}.events(itrial).code = event_labels{1,2};
        else
        end
        if ~isempty(idx_ctrl{itrial,1})
            block_header{iblock,1}.events(itrial).code = event_labels{1,3};
        else
        end
    end
    % overwrite the datasets with renamed label events
    CLW_save(sub_pre_proc_folder,block_header{iblock,1}, block_data{iblock,1});
    clear temp idx_date idx_ctrl idx_sensi idx_base
end


%% 5) Discard MEPs if baseline activity

% define baseline window to check for baseline activity 200 ms before TMS pulse
bsln_time_window = [-0.2 0];
bsln_sample_window = [1 round((abs(bsln_time_window(1))-abs(bsln_time_window(2)))*sr+1)];

% First check: for baseline single trial EMG activity (RMS) per block
% threshold is mean(RMS) + 2.5*std(RMS)
% (as in Sulcova D. et al. bioRxiv 2022; Grandjean and Duque NIMG 2020)

% prepare x values for plotting
xval = -0.2:1/sr:0.5;

% loop until no outliers is identified + automatic plot of the discarded trial
disp(' ')
disp('Iterative trial RMS for baseline above mean RMS')
disp(' ')
for iblock = 1:block_num

    run_loop = 1;
    iloop = 1;
    temp_data{iblock,iloop} = squeeze(block_data{iblock,1});

    while run_loop

           % compute RMS for each trial
        for itrial = 1:size(temp_data{iblock,iloop},1)
            bsln_rms{iblock,1}(itrial,1) = rms(temp_data{iblock,iloop}(itrial,bsln_sample_window(1):bsln_sample_window(2)));
        end

        % threshold mean RMS +/- 2.5*SD (permissive)
        avg_bsln_rms(iblock,iloop) = mean(bsln_rms{iblock,1});
        sd_bsln_rms(iblock,iloop) = std(bsln_rms{iblock,1});
        thrshld(iblock,iloop) = avg_bsln_rms(iblock,iloop)+2.5*sd_bsln_rms(iblock,iloop);

        idx_val = 1;
        idx_exc = 1;
        idx_exclud{iblock,iloop} = [];
        for itrial = 1:size(temp_data{iblock,iloop},1)
            if bsln_rms{iblock,1}(itrial,1) < thrshld(iblock,iloop)
                valid_data{iblock,1}(idx_val,:) = temp_data{iblock,iloop}(itrial,:);
                idx_val = idx_val+1;
            else
                idx_exclud{iblock,iloop}(idx_exc,1) = itrial;
                idx_exc = idx_exc+1;
                % plot discarded trial
                f = figure('Color','w','Position',[0 0 1500 700]);
                figure(f);
                temp_plot_data{iblock,1} = squeeze(block_data{iblock,1});
                plot(xval,temp_plot_data{iblock,1}(idx_exclud{iblock,iloop}(end),:),'k','LineWidth',1)
                xline(0,'--r','TMS','LabelOrientation','horizontal')
                yline(0,'-k')
                yline(bsln_rms{iblock,1}(idx_exclud{iblock,iloop}(end),1),'--k','RMS','LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left')
                yline(avg_bsln_rms(iblock,iloop),'--m','Average RMS','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right')
                ax = gca;
                ax.TickDir = 'out';
                ax.XLabel.String = 'time (s)';
                ax.YLabel.String = 'amplitude (µV)';
                ax.YLim = [-200 200];
                ax.Box = 'off';
                title({'First check iterative RMS'},{strcat(['sub-',subject_id,' MEP#',num2str(idx_exclud{iblock,iloop}(end)),' in block ',num2str(iblock)])})
                % pause()
                close (f)
            end
        end
        disp(strcat(['Block-',num2str(iblock),': ',num2str(length(idx_exclud{iblock,iloop})),' MEP discarded after iteration-',num2str(iloop)]))
        if isempty(idx_exclud{iblock,iloop})
            run_loop = 0;
        else
            iloop = iloop+1;
            temp_data{iblock,iloop} = valid_data{iblock,1};
            valid_data{iblock,1} = [];
            bsln_rms{iblock,1} = [];
            idx_exclud{iblock,iloop} = [];
        end
    end
end

% store indices of discarded trials
all_exclud = cell(block_num,1);
for iblock = 1:block_num
    for icleaning = 1:iloop
        all_exclud{iblock,1} = [all_exclud{iblock,1}; idx_exclud{iblock,icleaning}];
    end
end

% Second check: threshold on baseline RMS exceeding +/-15 µV
% (as in Sulcova et al. BioRxiv 2022; Morozova et al. Sci Reports 2024)
idx_exclud{1,size(idx_exclud,2)+1} = [];

for iblock = 1:block_num

    idx_val = 1;
    idx_exc = 1;
    for itrial = 1:size(valid_data{iblock,1},1)
        if bsln_rms{iblock,1}(itrial,1) < abs_thrshld
            valid_data{iblock,2}(idx_val,:) = valid_data{iblock,1}(itrial,:);
            idx_val = idx_val+1;
        else
            idx_exclud{iblock,end}(idx_exc,1) = itrial;
            idx_exc = idx_exc+1;
        end
    end
end
% Warning: display the number of MEP excluded for each block
disp(' ')
disp(['Threshold for baseline at ',num2str(abs_thrshld),' µV'])
disp(' ')
for iblock = 1:size(idx_exclud,1)
    if isempty(idx_exclud{iblock,end})
        disp(strcat(['Block-',num2str(iblock),': 0 MEP discarded']))
    else
        disp(strcat(['Block-',num2str(iblock),': ',num2str(size(idx_exclud{iblock,end},1)),' MEP discarded']))
    end
end

% plot the discarded trials
for iblock = 1:block_num
    for itrial = 1:size(idx_exclud{iblock,end},1)
        f = figure('Color','w','Position',[0 0 1500 700]);
        figure(f);

        temp_plot_data{iblock,1} = squeeze(block_data{iblock,1});
        plot(xval,temp_plot_data{iblock,1}(idx_exclud{iblock,end}(itrial,1),:),'k','LineWidth',1)
        
        xline(0,'--r','TMS','LabelOrientation','horizontal')
        yline(0,'-k')
        yline(bsln_rms{iblock,1}(idx_exclud{iblock,end}(itrial,1),1),'--k','RMS','LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left')
        yline(abs_thrshld,'--m','absolute threshold','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right')
        ax = gca;
        ax.TickDir = 'out';
        ax.XLabel.String = 'time (s)';
        ax.YLabel.String = 'amplitude (µV)';
        ax.YLim = [-200 200];
        ax.Box = 'off';
        title({'Second check RMS threshold'},{strcat(['sub-',subject_id,' MEP#',num2str(idx_exclud{iblock,end}(itrial,1)),' in block ',num2str(iblock)])})
        % pause()
        close(f)
    end
end

% merge the indices of all excluded trials
for iblock = 1:block_num
    all_exclud{iblock,1} = [all_exclud{iblock,1}; idx_exclud{iblock,end}];
end

% save the valid MEPs for each block
for iblock = 1:block_num
    clean_header{iblock,1} = block_header{iblock,1};
    clean_header{iblock,1}.name = strcat([clean_header{iblock,1}.name,' clean']);
    % discard events or trials from all_exclud
    clean_header{iblock,1}.events(all_exclud{iblock,1}) = [];
    % fix epochs numbering
    for ievent = 1:size(clean_header{iblock,1}.events,2)
        clean_header{iblock,1}.events(ievent).epoch = ievent;
    end
    clean_data{iblock,1}(:,1,1,1,1,:) = valid_data{iblock,2};
    clean_header{iblock,1}.datasize = size(clean_data{iblock,1});
    CLW_save(sub_pre_proc_folder,clean_header{iblock,1},clean_data{iblock,1});
end


%% 6) Split the valid MEP trials in separate files per event types

for iblock = 1:block_num
    out_datasets = RLW_segmentation2(clean_header{iblock,1}, clean_data{iblock,1},event_labels,'x_start',xstart,'x_duration',xduration);
    for ilabel = 1:size(out_datasets,2)
        CLW_save(sub_pre_proc_folder,out_datasets(ilabel).header, out_datasets(ilabel).data);
    end
    clear out_datasets
end


%% 7) extract MEPs peak-to-peak amplitude in specific window between 10 ms and 50 ms after TMS

response_window = round([0.210*sr+1 0.250*sr+1]);
% do it for each trial type sorted file, first list them
list_base = dir(fullfile(sub_pre_proc_folder,'*base *.mat'));
list_sensi = dir(fullfile(sub_pre_proc_folder,'*sensi *.mat'));
list_ctrl = dir(fullfile(sub_pre_proc_folder,'*ctrl *.mat'));

% for base trials
base_idx_bad = 1;
base_bad_mep = cell(size(list_base,1),1);
base_answ = cell(size(list_base,1),1);
for ifile = 1:size(list_base,1)
    [base_header{ifile,1}, base_data{ifile,1}] = CLW_load(fullfile(sub_pre_proc_folder,list_base(ifile).name));
    base_data{ifile,1} = squeeze(base_data{ifile,1});
    for itrial = 1:size(base_data{ifile,1},1)
        [base_max_p{ifile,itrial}, base_max_lat{ifile,itrial}] = max(base_data{ifile,1}(itrial,response_window(1):response_window(2)));
        [base_min_p{ifile,itrial}, base_min_lat{ifile,itrial}] = min(base_data{ifile,1}(itrial,response_window(1):response_window(2)));

        % check if time inteval between positive and negative maxima is ok
        if abs(base_max_lat{ifile,itrial} - base_min_lat{ifile,itrial}) > round(0.015*sr,1)
            base_bad_mep{ifile,1}(base_idx_bad,1) = itrial;
            disp(' ')
            disp(strcat(['Time interval between peak is too long check trial#',num2str(itrial),' in: ',list_base(ifile).name]));
            disp(' ')
            figure, plot(base_data{ifile,1}(itrial,response_window(1):response_window(2))), hold on, plot(base_max_lat{ifile,itrial},base_max_p{ifile,itrial},'*r','MarkerSize',10)
            plot(base_min_lat{ifile,itrial},base_min_p{ifile,itrial},'*r','MarkerSize',10), title({strcat('trial#',num2str(itrial))},{list_base(ifile).name})
            base_answ{ifile,1}(base_idx_bad,1) = questdlg('Discard this MEP?','P-2-P?','yes','no','yes');
            base_idx_bad = base_idx_bad+1;
        else
        end
    end
end
% for ctrl trials
ctrl_idx_bad = 1;
ctrl_bad_mep = cell(size(list_ctrl,1),1);
ctrl_answ = cell(size(list_ctrl,1),1);
for ifile = 1:size(list_ctrl,1)
    [ctrl_header{ifile,1}, ctrl_data{ifile,1}] = CLW_load(fullfile(sub_pre_proc_folder,list_ctrl(ifile).name));
    ctrl_data{ifile,1} = squeeze(ctrl_data{ifile,1});
    for itrial = 1:size(ctrl_data{ifile,1},1)
        [ctrl_max_p{ifile,itrial}, ctrl_max_lat{ifile,itrial}] = max(ctrl_data{ifile,1}(itrial,response_window(1):response_window(2)));
        [ctrl_min_p{ifile,itrial}, ctrl_min_lat{ifile,itrial}] = min(ctrl_data{ifile,1}(itrial,response_window(1):response_window(2)));

        % check if time inteval between positive and negative maxima is ok
        if abs(ctrl_max_lat{ifile,itrial} - ctrl_min_lat{ifile,itrial}) > round(0.015*sr,1)
            ctrl_bad_mep{ifile,1}(ctrl_idx_bad,1) = itrial;
            disp(' ')
            disp(strcat(['Time interval between peak is too long check trial#',num2str(itrial),' in: ',list_ctrl(ifile).name]));
            disp(' ')
            figure, plot(ctrl_data{ifile,1}(itrial,response_window(1):response_window(2))), hold on, plot(ctrl_max_lat{ifile,itrial},ctrl_max_p{ifile,itrial},'*r','MarkerSize',10)
            plot(ctrl_min_lat{ifile,itrial},ctrl_min_p{ifile,itrial},'*r','MarkerSize',10), title({strcat('trial#',num2str(itrial))},{list_ctrl(ifile).name})
            ctrl_answ{ifile,1}(ctrl_idx_bad,1) = questdlg('Discard this MEP?','P-2-P?','yes','no','yes');
            ctrl_idx_bad = ctrl_idx_bad+1;
        else
        end
    end
end

% for sensi trials
sensi_idx_bad = 1;
sensi_bad_mep = cell(size(list_sensi,1),1);
sensi_answ = cell(size(list_sensi,1),1);
for ifile = 1:size(list_sensi,1)
    [sensi_header{ifile,1}, sensi_data{ifile,1}] = CLW_load(fullfile(sub_pre_proc_folder,list_sensi(ifile).name));
    sensi_data{ifile,1} = squeeze(sensi_data{ifile,1});
    for itrial = 1:size(sensi_data{ifile,1},1)
        [sensi_max_p{ifile,itrial}, sensi_max_lat{ifile,itrial}] = max(sensi_data{ifile,1}(itrial,response_window(1):response_window(2)));
        [sensi_min_p{ifile,itrial}, sensi_min_lat{ifile,itrial}] = min(sensi_data{ifile,1}(itrial,response_window(1):response_window(2)));

        % check if time inteval between positive and negative maxima is ok
        if abs(sensi_max_lat{ifile,itrial} - sensi_min_lat{ifile,itrial}) > round(0.015*sr,1)
            sensi_bad_mep{ifile,1}(sensi_idx_bad,1) = itrial;
            disp(' ')
            disp(strcat(['Time interval between peak is too long check trial#',num2str(itrial),' in: ',list_sensi(ifile).name]));
            disp(' ')
            figure, plot(sensi_data{ifile,1}(itrial,response_window(1):response_window(2))), hold on, plot(sensi_max_lat{ifile,itrial},sensi_max_p{ifile,itrial},'*r','MarkerSize',10)
            plot(sensi_min_lat{ifile,itrial},sensi_min_p{ifile,itrial},'*r','MarkerSize',10), title({strcat('trial#',num2str(itrial))},{list_sensi(ifile).name})
            sensi_answ{ifile,1}(sensi_idx_bad,1) = questdlg('Discard this MEP?','P-2-P?','y','n','y');
            sensi_idx_bad = sensi_idx_bad+1;
        else
        end
    end
end


%% 8) Store and save MEPs amplitudes in a single .mat file

MEP = struct;
MEP.sub_info = sub_info;
% for base trials
for ifile = 1:size(list_base,1)
    if contains(list_base(ifile).name,'BLK1')
        for itrial = 1:size(base_data{ifile,1},1)
            MEP.base.amp_PRE{ifile,1}(itrial,1) = (abs(base_max_p{ifile,itrial}) + abs(base_min_p{ifile,itrial}));
            if base_max_lat{ifile,itrial} < base_min_lat{ifile,itrial}
                MEP.base.lat_PRE{ifile,1}(itrial,1) = (base_max_lat{ifile,itrial}+10)/sr;
            else
                MEP.base.lat_PRE{ifile,1}(itrial,1) = (base_min_lat{ifile,itrial}+10)/sr;
            end
        end
    else
        for itrial = 1:size(base_data{ifile,1},1)
            MEP.base.amp_PST{ifile,1}(itrial,1) = (abs(base_max_p{ifile,itrial}) + abs(base_min_p{ifile,itrial}));
            if base_max_lat{ifile,itrial} < base_min_lat{ifile,itrial}
                MEP.base.lat_PST{ifile,1}(itrial,1) = (base_max_lat{ifile,itrial}+10)/sr;
            else
                MEP.base.lat_PST{ifile,1}(itrial,1) = (base_min_lat{ifile,itrial}+10)/sr;
            end
        end
    end
end

% for control trials
for ifile = 1:size(list_ctrl,1)
    for itrial = 1:size(ctrl_data{ifile,1},1)
        MEP.ctrl.amp_PST{ifile,1}(itrial,1) = (abs(ctrl_max_p{ifile,itrial}) + abs(ctrl_min_p{ifile,itrial}));
        if ctrl_max_lat{ifile,itrial} < ctrl_min_lat{ifile,itrial}
            MEP.ctrl.lat_PST{ifile,1}(itrial,1) = (ctrl_max_lat{ifile,itrial}+10)/sr;
        else
            MEP.ctrl.lat_PST{ifile,1}(itrial,1) = (ctrl_min_lat{ifile,itrial}+10)/sr;
        end
    end
end

% for sensisitized trials
for ifile = 1:size(list_sensi,1)
    for itrial = 1:size(sensi_data{ifile,1},1)
        MEP.sensi.amp_PST{ifile,1}(itrial,1) = (abs(sensi_max_p{ifile,itrial}) + abs(sensi_min_p{ifile,itrial}));
        if sensi_max_lat{ifile,itrial} < sensi_min_lat{ifile,itrial}
            MEP.sensi.lat_PST{ifile,1}(itrial,1) = (sensi_max_lat{ifile,itrial}+10)/sr;
        else
            MEP.sensi.lat_PST{ifile,1}(itrial,1) = (sensi_min_lat{ifile,itrial}+10)/sr;
        end
    end
end


% update and save number of discarded MEPs
for ifile = 1:size(base_answ,1)
    if ~isempty(base_answ{ifile,1})
        MEP.base.warning_exclud = base_answ{ifile,1};
    else
        MEP.base.warning_exclud = [];
    end
end
for ifile = 1:size(ctrl_answ,1)
    if ~isempty(ctrl_answ{ifile,1})
        MEP.ctrl.warning_exclud = ctrl_answ{ifile,1};
    else
        MEP.ctrl.warning_exclud = [];
    end
end
for ifile = 1:size(sensi_answ,1)
    if ~isempty(sensi_answ{ifile,1})
        MEP.sensi.warning_exclud = sensi_answ{ifile,1};
    else
        MEP.sensi.warning_exclud = [];
    end
end
MEP.excludedPerBlock = all_exclud;

% save indices of suspicious MEPs
MEP.base.warningMEP = base_bad_mep;
MEP.ctrl.warningMEP = ctrl_bad_mep;
MEP.sensi.warningMEP = sensi_bad_mep;

save(fullfile(sub_results_folder,strcat(savename,'_meps')),'MEP')


%% 9) Exclude suspicious MEPs and concatenate MEPs amplitudes across PRE and POST capsaicin blocks 

temp_base_PRE = [];
temp_base_PST = [];
for ifile = 1:size(list_base,1)
    if contains(list_base(ifile).name,'BLK1')
        if ~isempty(base_answ{ifile,1})
            MEP.base.amp_PRE{ifile,1}(base_bad_mep{ifile,1}) = [];
        else
        end
        temp_base_PRE = MEP.base.amp_PRE{ifile,1};
    else
        if ~isempty(base_answ{ifile,1})
            MEP.base.amp_PST{ifile,1}(base_bad_mep{ifile,1}) = [];
        else
        end
        temp_base_PST = [temp_base_PST;MEP.base.amp_PST{ifile,1}];
    end
end

temp_ctrl_PST = [];
for ifile = 1:size(list_ctrl,1)
    if ~isempty(ctrl_answ{ifile,1})
        MEP.ctrl.amp_PST{ifile,1}(ctrl_bad_mep{ifile,1}) = [];
    else
    end
    temp_ctrl_PST = [temp_ctrl_PST;MEP.ctrl.amp_PST{ifile,1}];
end

temp_sensi_PST = [];
for ifile = 1:size(list_sensi,1)
    if ~isempty(sensi_answ{ifile,1})
        MEP.sensi.amp_PST{ifile,1}(sensi_bad_mep{ifile,1}) = [];
    else
    end
    temp_sensi_PST = [temp_sensi_PST;MEP.sensi.amp_PST{ifile,1}];
end

% save MEP structure
save(fullfile(sub_results_folder,strcat(savename,'_meps')),'MEP')


%% 10) Export the MEPs amplitude data for the current subject per trial type in .csv files

% baseline
[max_size, idx_max] = max([size(temp_base_PRE,1) size(temp_base_PST,1)]);
if idx_max == 1
    temp_base_PST(end+1:max_size,1) = NaN;
elseif idx_max == 2
    temp_base_PRE(end+1:max_size,1) = NaN;
end
base_table = array2table([temp_base_PRE,temp_base_PST],'VariableNames',{'base_amp_PRE','base_amp_PST'});
base_csv_name = strcat(savename,'_meps_base.csv');
writetable(base_table,fullfile(sub_results_folder,base_csv_name),'Delimiter', ',')

% WASP control
ctrl_table = array2table(temp_ctrl_PST,'VariableNames',{'ctrl_amp_PST'});
ctrl_csv_name = strcat(savename,'_meps_ctrl.csv');
writetable(ctrl_table,fullfile(sub_results_folder,ctrl_csv_name),'Delimiter', ',')

% WASP sensitized
sensi_table = array2table(temp_sensi_PST,'VariableNames',{'sensi_amp_PST'});
sensi_csv_name = strcat(savename,'_meps_sensi.csv');
writetable(sensi_table,fullfile(sub_results_folder,sensi_csv_name),'Delimiter', ',')

end