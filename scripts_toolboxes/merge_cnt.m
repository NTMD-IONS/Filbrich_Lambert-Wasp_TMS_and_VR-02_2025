% merge and concatenate 2 .cnt files from EMG to be used in
% mep_process.m script

% first import manually using letswave the 2 .cnt files from the different sessions in the
% first subfolder session that has been crerated

% choose the Sessions subfolder
session_folder = uigetdir('F:\','Select the main folder');
cd(session_folder)

% list the imported lw files
list_lw = dir(fullfile(session_folder,'*emg.mat'));

% load the lw files of the imported datasets
for ifile = 1:size(list_lw,1)
    [header{ifile,1}, data{ifile,1}] = CLW_load(list_lw(ifile).name);
    temp_data{ifile,1} = squeeze(data{ifile,1});
end

% concatenate the datasets
concat_data(1,1,1,1,1,:) = [temp_data{1,1}; temp_data{2,1}];

% concatenate headers
concat_header = header{1,1};
concat_header.name = strcat('concat_',concat_header.name);
concat_header.datasize(6) = size(concat_data,6);
icount = 1;
for ievent = 1:size(header{2,1}.events,2)
    if strcmp(header{2,1}.events(ievent).code,'1')
        concat_header.events(size(header{1,1}.events,2)+icount) = header{2,1}.events(ievent);
        concat_header.events(size(header{1,1}.events,2)+icount).latency = header{2,1}.events(ievent).latency + header{1,1}.datasize(6)*header{1,1}.xstep;
        icount = icount+1;
    else
    end
end

% save lw file
CLW_save(session_folder,concat_header,concat_data)