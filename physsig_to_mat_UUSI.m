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
    data = acq.data;

    % extract pulse, respiration and time stamps
    timestamp = data(:,6);
    pulse = data(:,2);
    resp = data(:,1);
       
    figure;plot(timestamp);title('original timestamp');

    % change the timestamp vector so that it is constant (=5) during EPI
    time_id = find(timestamp==max(timestamp));
    for i=1:length(time_id)-1
        if time_id(i+1)-time_id(i) <= 2000
            timestamp(time_id(i):time_id(i+1)) = max(timestamp);
        end
    end

    figure(2);plot(timestamp);title('constant timestamp');

    % find EPI series matching the length of rs-fMRI
    vec=[0; ones([360000,1])*max(timestamp)];
    rsepi_id=strfind(timestamp',vec'); %id's corresponding to starting point of rs-epi

    %%
    fs = 1000; %sampling rate in Hz
    tr = 2; % TR in s
    vol = 184; %number of volumes in rs-fMRI

    % extract signals for fMRI1
    start1 = rsepi_id(1);
    end1 = start1 + vol*tr*fs;
    pulse1 = pulse(start1:end1);
    resp1 = resp(start1:end1);

    % extract signals for fMRI2
    start2 = rsepi_id(2);
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