% Convert Biopac data file from .acq to .mat vector to be used with DRIFTER
% Revis data
%
% Hanna Halme 02/2016
% hanna.halme@hus.fi

clear all
close all

sessions = ['1.sessio/';'2.sessio/';'3.sessio/'];

for s=1:size(sessions,1)
    session = sessions(s,:);
    
    % Path to Biopac files
    path = ['/opt/MR_data/Silja_Raty/Revis_0007_rs/Biopac/',session];
    fname = dir(path);
    file = ([path, fname(3).name]);

    % load signal file (might take a while)
    if exist(file, 'file') == 2
        acq=load_acq(file);
    end

    % cut the whole data before first fMRI series
    mri_start = find(acq.data(:,6)==max(acq.data(:,6)),1);
    data = acq.data(mri_start:end,:);
    data = data(2200000:end,:);

    % extract pulse, respiration and time stamps
    timestamp = data(:,6);
    pulse = data(:,2);
    resp = data(:,1);

    fs = 1000; %sampling rate in Hz
    tr = 2; % TR in s
    vol = 184; %number of volumes in rs-fMRI

    % extract signals for fMRI1
    start1 = find(timestamp==max(timestamp),1);
    end1 = start1 + vol*tr*fs;
    pulse1 = pulse(start1:end1);
    resp1 = resp(start1:end1);
    pulse(1:end1) = [];
    resp(1:end1) = [];
    timestamp(1:end1) = [];

    % extract signals for fMRI2
    start2 = find(timestamp==max(timestamp),1);
    end2 = start2 + vol*tr*fs;
    pulse2 = pulse(start2:end2);
    resp2 = resp(start2:end2);
    
    

    % save signals as .mat 1xN vectors
    disp('Saving signal files...')
    save([path, 'pulse_fMRI1'],'pulse1');
    save([path, 'resp_fMRI1'],'resp1');
    save([path, 'pulse_fMRI2'],'pulse2');
    save([path, 'resp_fMRI2'],'resp2');

end

disp('ALL FILES READY!');