function Data = run_drifter_noSPM_kopio(fMRI,ref_data_puls,ref_data_resp)

    addpath('/opt/MATLAB/NIfTI_20140122/')
          
    fmri = load_untouch_nii(fMRI);
    obs_data = double(fmri.img);
  
    % Provide two structures: "data" and "refdata" such that
    data.data = obs_data;
    data.dt = fmri.hdr.dime.pixdim(5)*1000; %tr in milliseconds

    refdata{1}.dt = 1/50;
    refdata{1}.freqlist = 50:120; % Vector of possible frequencies in bpm
    refdata{1}.N = 1; % 
    refdata{1}.Nimm = 2;
    refdata{1}.downdt = 0.1;

    refdata{2}.dt = 1/50;
    refdata{2}.freqlist = 6:15; % Vector of possible frequencies in bpm 
    refdata{2}.N = 2;
    refdata{2}.Nimm = 2; 
    refdata{2}.downdt = 0.1;

    if exist(ref_data_puls, 'file') == 2
        puls = load(ref_data_puls);
        puls = puls.sig;
    	refdata{1}.data = puls;
    	%refdata{1}.poverall = 0.1;
    else
	disp('MISSING DATA!!!')
    end
        
    if exist(ref_data_resp, 'file') == 2
        resp = load(ref_data_resp);
        resp = resp.sig;
    	refdata{2}.data = resp; 
    	%refdata{2}.poverall = 0.1;
    else
	disp('MISSING DATA!!!')
    end

  
    % Run DRIFTER if at least one phsyiological signal exists
    if exist(ref_data_puls, 'file') == 2 | exist(ref_data_resp, 'file') == 2
    	[new_data,refdata] = drifter(data,refdata);
	save_nii(make_nii(new_data.estimate, [], [], 16),'drifter_corrected.nii.gz');
    	%save_nii(make_nii(new_data.noise, [], [], 16),'drifter_noise.nii.gz');
 	Data = new_data.estimate;
    else
	save_nii(make_nii(data.data, [], [], 16),'drifter_corrected.nii.gz');
	Data = data.data;
    end
   

end
