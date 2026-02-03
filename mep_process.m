% function to process MEPs and rename trials based on VR txt files
% it also extract MEPs amplitudes and latencies for each block.
% 
% TO DO:
% decide on folder and subject level naming of the acquisitions
% sub-00$
% main_folder contains the raw_data folder, pre-processed data folder,
% results folder, (scripts folder?)
% 
% Block 1 (BLK1) is the block pre-HFS
% Block 2 to 5 (BLK2 to BLK5) are the blocks post-HFS
% 
% (sub-001 is test_giulia)
% (sub-002 is Laurie)
%
% the fucntion requires:
%   the package natsortfiles 
%   letswave6
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Cédric Lenoir, NeTMeD, IoNS, UCLouvain, January 2026

function mep_process


%% 1) LOAD DATA AND SET A FEW PARAMETERS

% set the paths and initialize letswave
% Select the general data folder which contains all subjects folders
main_folder = uigetdir('C:\','Select main folder');
cd(main_folder)
raw_data_folder = fullfile(main_folder,'raw_data');
results_folder = fullfile(main_folder,'results');
pre_proc_folder = fullfile(main_folder,'pre-processed');

if ~isfolder(raw_data_folder); mkdir(raw_data_folder); end
if ~isfolder(results_folder); mkdir(results_folder); end
if ~isfolder(pre_proc_folder); mkdir(pre_proc_folder); end

addpath(genpath(main_folder))

% initialize letswave 6
% check if lw is already on the path if not select the folder of the toolobx
if contains(path, 'C:\Users\cedlenoir\Documents\MATLAB\letswave6-master')
    letswave();
    clc
else
    lw_path = uigetdir('C:\Users\cedlenoir\Documents\MATLAB','Select Letswave 6 folder');
    addpath(genpath(lw_path))
    letswave();
    clc
end
% check if natsorfiles toolbox is on the path
if contains(path, 'C:\Users\cedlenoir\Documents\MATLAB\natsortfiles')
else
    natsortfiles_path = uigetdir('C:\Users\cedlenoir\Documents\MATLAB','Select natsortfiles folder');
    addpath(genpath(natsortfiles_path))
    clc
end
% dialog box
prompt = {'\fontsize{12} Subject ID? :','\fontsize{12} Sensitized arm (L/R): ','\fontsize{12} Stimulated hemisphere (L/R):',...
    '\fontsize{12} TMS PRE-CAPS BLK1 start index? : ','\fontsize{12} TMS PRE-CAPS BLK1 stop index? : ',...
    '\fontsize{12} TMS POST-CAPS BLK2 start index? : ','\fontsize{12} TMS POST-CAPS BLK2 stop index? : ',...
    '\fontsize{12} TMS POST-CAPS BLK3 start index? : ','\fontsize{12} TMS POST-CAPS BLK3 stop index? : ',...
    '\fontsize{12} TMS POST-CAPS BLK4 start index? : ','\fontsize{12} TMS POST-CAPS BLK4 stop index? : ',...
    '\fontsize{12} TMS POST-CAPS BLK5 start index? : ','\fontsize{12} TMS POST-CAPS BLK5 stop index? : ',...
    '\fontsize{12} EMG trial plots? : YES -> 1 OR NO -> 0','\fontsize{12} Comments: ',};
dlgtitle = 'LOAD SUBJECT DATA';
opts.Interpreter = 'tex';
dims = repmat([1 80],15,1);
definput = {'002','R','L','114','137','138','164','165','192','193','220','221','246','1',''};
info = inputdlg(prompt,dlgtitle,dims,definput,opts);
subject_id = char(info(1));
sensi_arm = char(info(2));
stim_hemi = char(info(3));
ixd_tms_pre_start = str2double(info(4));
ixd_tms_pre_stop = str2double(info(5));
ixd_tms_pst_start = str2double([info(6) info(8) info(10) info(12)]);
ixd_tms_pst_stop = str2double([info(7) info(9) info(11) info(13)]);
plot_opt_EMG = str2double(info(14));
notes = str2double(info(15));
block_num = 5;
% store info in structure
sub_info = struct;
sub_info.sub_ID = subject_id;
sub_info.sensitized_arm = sensi_arm;
sub_info.stimulated_hemisph = stim_hemi;
sub_info.TMS_triggers.pre = [ixd_tms_pre_start;ixd_tms_pre_stop];
sub_info.TMS_triggers.pst = [ixd_tms_pst_start;ixd_tms_pst_stop];
sub_info.comments = notes;

% folder name of Visor EMG data
gen_subject_folder = fullfile(raw_data_folder,strcat('sub-',subject_id));
subject_folder = fullfile(raw_data_folder,strcat('sub-',subject_id),'Sessions');
sub_pre_proc_folder = fullfile(pre_proc_folder,sprintf('sub-%s',subject_id));
if ~isfolder(sub_pre_proc_folder); mkdir(sub_pre_proc_folder); end

% list sessions
session_list = dir(subject_folder);
session_list = session_list(~ismember({session_list.name}, {'.', '..'}));

% list EMG CNT files
file_list = dir(fullfile(session_list.folder,session_list.name,'*emg.cnt*'));

% import CNT file
[out_data,~] = RLW_import_CNT(fullfile(session_list.folder,session_list.name,file_list.name));
filename = file_list.name(1:end-8);
savename = strcat('sub-',subject_id,'_',filename);
out_data.header.name = savename;
CLW_save(sub_pre_proc_folder,out_data.header,out_data.data);
clear out_data
[header, data] = CLW_load(fullfile(sub_pre_proc_folder,savename));
sr = 1/header.xstep;

% pre-processign steps in lw
% DC removal
[dc_header, dc_data] = RLW_dc_removal(header,data,'linear_detrend',0);

% high pass filter Butterworth 4 Hz; order 4
low_cutoff = 4;
order = 4;
[filt_header, filt_data] = RLW_butterworth_filter(dc_header,dc_data,'filter_type','highpass','low_cutoff',low_cutoff,'filter_order',order);

% segmentation
xstart = -0.2;
xduration = 0.7;
[seg_header, seg_data] = RLW_segmentation(filt_header, filt_data, {'1'},'x_start',xstart,'x_duration',xduration);
seg_header.chanlocs.labels = 'EMG1';
seg_header.name = strcat(header.name,' DC HPfilt ep');

% save dataset
CLW_save(sub_pre_proc_folder,seg_header, seg_data);

% keep TMS events after rMT
% arrange epochs of valid TMS triggers for each TMS block
valid_tms_idx = cell(block_num,1);
for iblock = 1:block_num
    if iblock == 1
        valid_tms_idx{iblock,1} = (ixd_tms_pre_start:ixd_tms_pre_stop);
    else
        valid_tms_idx{iblock,1} = (ixd_tms_pst_start(iblock-1):ixd_tms_pst_stop(iblock-1));
    end
end
for iblock = 1:block_num
    [block_header{iblock,1}, block_data{iblock,1}] = RLW_arrange_epochs(seg_header, seg_data,valid_tms_idx{iblock,1});
    block_header{iblock,1}.name = strcat(seg_header.name,' BLK ',num2str(iblock));
    CLW_save(sub_pre_proc_folder,block_header{iblock,1}, block_data{iblock,1});
end

% Sort MEPs according to trial condition per block ans save lw file for each events
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

% define baseline window to check for baseline activity
bsln_time_window = [-0.2 -0.01];
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
for iblock = 1:5

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
                title(strcat(['sub-',subject_id,' MEP#',num2str(idx_exclud{iblock,iloop}(end)),' in block ',num2str(iblock)]))
                pause()
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
all_exclud = cell(5,1);
for iblock = 1:block_num
    for icleaning = 1:iloop
        all_exclud{iblock,1} = [all_exclud{iblock,1}, idx_exclud{iblock,icleaning}];
    end
end

% Second check: threshold on baseline RMS exceeding +/-15 µV
% (as in Sulcova et al. BioRxiv 2022; Morozova et al. Sci Reports 2024)
abs_thrshld = 15;
idx_exclud{1,size(idx_exclud,2)+1} = [];

for iblock = 1:5

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
for iblock = 1:5
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
        title(strcat(['sub-',subject_id,' MEP#',num2str(idx_exclud{iblock,end}(itrial,1)),' in block ',num2str(iblock)]))
        pause()
        close(f)
    end
end

% merge the indices of all excluded trials
for iblock = 1:block_num
    all_exclud{iblock,1} = sort([all_exclud{iblock,1}; idx_exclud{iblock,end}]);
end

% save the valid MEPs for each block
for iblock = 1:size(valid_data,1)
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

% slipt the valid MEP trials in separate files per event types
for iblock = 1:block_num
    out_datasets = RLW_segmentation2(clean_header{iblock,1}, clean_data{iblock,1},event_labels,'x_start',xstart,'x_duration',xduration);
    for ilabel = 1:size(out_datasets,2)
        CLW_save(sub_pre_proc_folder,out_datasets(ilabel).header, out_datasets(ilabel).data);
    end
    clear out_datasets
end

% define MEP peak-to-peak amplitude criterion > 50 µV in specific window (s)
mep_threshold = 50;
response_window = round([0.205*sr+1 0.265*sr+1]);
for iblock = 1:size(valid_data,1)
    idx_bad = 1;
    bad_mep{iblock,idx_bad} = [];
    for itrial = 1:size(valid_data{iblock,2},1)
        [max_p{iblock,itrial}, max_lat{iblock,itrial}] = max(valid_data{iblock,2}(itrial,response_window(1):response_window(2)));
        [min_p{iblock,itrial}, min_lat{iblock,itrial}] = min(valid_data{iblock,2}(itrial,response_window(1):response_window(2)));

        % check if time between positive and negative maxima is ok OR if MEP amplitude >= 50 µV
        if abs(max_lat{iblock,itrial} - min_lat{iblock,itrial}) > round(0.015*sr,1) || (max_p{iblock,itrial} + abs(min_p{iblock,itrial})) < mep_threshold
            bad_mep{iblock,1}(idx_bad,1) = itrial;
            idx_bad = idx_bad+1;
        else
        end
    end
end

% store MEPs amplitude
for iblock = 1:size(valid_data,1)
    for itrial = 1:size(valid_data{iblock,2},1)
        if isempty(bad_mep{iblock,1})
            mep_amp{iblock,1}(itrial,1) = (max_p{iblock,itrial} + abs(min_p{iblock,itrial}));
            if max_lat{iblock,itrial} < min_lat{iblock,itrial}
                mep_lat{iblock,1}(itrial,1) = max_lat{iblock,itrial};
            else
                mep_lat{iblock,1}(itrial,1) = min_lat{iblock,itrial};
            end
        elseif ~isempty(bad_mep{iblock,1})
            

    end
end

% extract MEP latencies




end