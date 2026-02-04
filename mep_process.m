% function to process MEPs and rename trials based on VR txt files
% it also extracts MEPs amplitudes and latencies for each block.
% 
% Folder and subject names naming are build whatever names are encoded in
% the VR software
% sub-00$
% main_folder contains the raw_data folder, pre-processed data folder,
% results folder, scripts folder (with toolboxes)
% 
% Block 1 (BLK1) is the pre-HFS block
% Block 2 to 5 (BLK2 to BLK5) are the post-HFS blocks
% 
% 
% 
% 
%
%
% Required toolboxes:
%   natsortfiles (2022, Stephen Cobeldick)
%   letswave 6
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
close all hidden

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
    '\fontsize{12} Baseline RMS threshold (µV): ','\fontsize{12} Comments: ',};
dlgtitle = 'LOAD SUBJECT DATA';
opts.Interpreter = 'tex';
dims = repmat([1 80],15,1);
definput = {'002','R','L','114','137','138','164','165','192','193','220','221','246','15',''};
info = inputdlg(prompt,dlgtitle,dims,definput,opts);
subject_id = char(info(1));
sensi_arm = char(info(2));
stim_hemi = char(info(3));
ixd_tms_pre_start = str2double(info(4));
ixd_tms_pre_stop = str2double(info(5));
ixd_tms_pst_start = str2double([info(6) info(8) info(10) info(12)]);
ixd_tms_pst_stop = str2double([info(7) info(9) info(11) info(13)]);
abs_thrshld = str2double(info(14));
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
sub_results_folder = fullfile(results_folder,sprintf('sub-%s',subject_id));
if ~isfolder(sub_pre_proc_folder); mkdir(sub_pre_proc_folder); end
if ~isfolder(sub_results_folder); mkdir(sub_results_folder); end

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


%% 2) Pre-processing steps in letswave

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


%% 3) Remove irrelevant TMS events

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
all_exclud = cell(block_num,1);
for iblock = 1:block_num
    for icleaning = 1:iloop
        all_exclud{iblock,1} = [all_exclud{iblock,1}, idx_exclud{iblock,icleaning}];
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
        pause()
        close(f)
    end
end

% merge the indices of all excluded trials
for iblock = 1:block_num
    all_exclud{iblock,1} = sort([all_exclud{iblock,1}; idx_exclud{iblock,end}]);
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
% for base trials
for ifile = 1:size(list_base,1)
    for itrial = 1:size(base_data{ifile,1},1)
        MEP.base.amp{ifile,1}(itrial,1) = (abs(base_max_p{ifile,itrial}) + abs(base_min_p{ifile,itrial}));
        if base_max_lat{ifile,itrial} < base_min_lat{ifile,itrial}
            MEP.base.lat{ifile,1}(itrial,1) = (base_max_lat{ifile,itrial}+10)/sr;
        else
            MEP.base.lat{ifile,1}(itrial,1) = (base_min_lat{ifile,itrial}+10)/sr;
        end
    end
end

% for control trials
for ifile = 1:size(list_ctrl,1)
    for itrial = 1:size(ctrl_data{ifile,1},1)
        MEP.ctrl.amp{ifile,1}(itrial,1) = (abs(ctrl_max_p{ifile,itrial}) + abs(ctrl_min_p{ifile,itrial}));
        if ctrl_max_lat{ifile,itrial} < ctrl_min_lat{ifile,itrial}
            MEP.ctrl.lat{ifile,1}(itrial,1) = (ctrl_max_lat{ifile,itrial}+10)/sr;
        else
            MEP.ctrl.lat{ifile,1}(itrial,1) = (ctrl_min_lat{ifile,itrial}+10)/sr;
        end
    end
end

% for sensisitized trials
for ifile = 1:size(list_sensi,1)
    for itrial = 1:size(sensi_data{ifile,1},1)
        MEP.sensi.amp{ifile,1}(itrial,1) = (abs(sensi_max_p{ifile,itrial}) + abs(sensi_min_p{ifile,itrial}));
        if sensi_max_lat{ifile,itrial} < sensi_min_lat{ifile,itrial}
            MEP.sensi.lat{ifile,1}(itrial,1) = (sensi_max_lat{ifile,itrial}+10)/sr;
        else
            MEP.sensi.lat{ifile,1}(itrial,1) = (sensi_min_lat{ifile,itrial}+10)/sr;
        end
    end
end
% save number of discarded MEPs
MEP.excludedPerBlock = all_exclud;
% save indices of suspicious MEPs
MEP.base.warningMEP = base_bad_mep;
MEP.base.warningMEP = ctrl_bad_mep;
MEP.base.warningMEP = sensi_bad_mep;
save(fullfile(sub_results_folder,strcat(savename,'_meps')),'MEP')


%% 9) Concatenate MEPs amplitudes across blocks and exclude suspicious MEPs
variable_names = {'base', 'ctrl', 'sensi'};
temp_base = [];
for ifile = 1:size(MEP.base.amp,1)
    if ~isempty(base_answ{ifile,1})
        MEP.base.amp{ifile,1}(base_bad_mep{ifile,1}) = [];
        temp_base = [temp_base;MEP.base.amp{ifile,1}];
    else
        temp_base = [temp_base;MEP.base.amp{ifile,1}];
    end
end

temp_ctrl = [];
for ifile = 1:size(MEP.ctrl.amp,1)
    if ~isempty(ctrl_answ{ifile,1})
        MEP.ctrl.amp{ifile,1}(ctrl_bad_mep{ifile,1}) = [];
        temp_ctrl = [temp_ctrl;MEP.ctrl.amp{ifile,1}];
    else
        temp_ctrl = [temp_ctrl;MEP.ctrl.amp{ifile,1}];
    end
end

temp_sensi = [];
for ifile = 1:size(MEP.sensi.amp,1)
    if ~isempty(sensi_answ{ifile,1})
        MEP.sensi.amp{ifile,1}(sensi_bad_mep{ifile,1}) = [];
        temp_sensi = [temp_sensi;MEP.sensi.amp{ifile,1}];
    else
        temp_sensi = [temp_sensi;MEP.sensi.amp{ifile,1}];
    end
end


%% 10) Export the MEPs amplitude data for the current subject per trial type in .csv files
% baseline
base_table = array2table(temp_base,'VariableNames',{'base_amp'});
base_csv_name = strcat(savename,'_meps_base.csv');
writetable(base_table,fullfile(sub_results_folder,base_csv_name),'Delimiter', ',')
% WASP control
ctrl_table = array2table(temp_ctrl,'VariableNames',{'ctrl_amp'});
ctrl_csv_name = strcat(savename,'_meps_ctrl.csv');
writetable(ctrl_table,fullfile(sub_results_folder,ctrl_csv_name),'Delimiter', ',')
% WASP sensitized
sensi_table = array2table(temp_sensi,'VariableNames',{'sensi_amp'});
sensi_csv_name = strcat(savename,'_meps_sensi.csv');
writetable(sensi_table,fullfile(sub_results_folder,sensi_csv_name),'Delimiter', ',')

end