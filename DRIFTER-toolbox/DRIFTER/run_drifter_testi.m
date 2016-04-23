function output_list = run_drifter_testi(fMRI,ref_data_puls,ref_data_resp)
  addpath /usr/local/MATLAB/R2015a/spm8

  ref_data_puls = load(ref_data_puls);
  ref_data_puls = ref_data_puls.A;
  ref_data_resp = load(ref_data_resp);
  ref_data_resp = ref_data_resp.A;

  spm('defaults','fmri');
  spm_jobman('initcfg'); % useful in SPM8 only 

% METHOD SPECIFICATION

  % Choose in which mode to run the method (0 = DRIFTER with all noise 
  % removed, 1 = DRIFTER with only physiological noise removed, 
  % 2 = a RETROICOR implementation)
  matlabbatch{1}.spm.tools.drifter.mode = 1;
  
  % Choose with which prefix to output the data
  matlabbatch{1}.spm.tools.drifter.prefix = 'n';
      
  % Assign the EPI fMRI data files (as we have not done anything prior to 
  % this, we use the original data files loaded into f)
  fmri = load_untouch_nii(fMRI);
  matlabbatch{1}.spm.tools.drifter.epidata.files = fmri.img;

    % The  Interscan interval, TR, (specified in milliseconds)
  matlabbatch{1}.spm.tools.drifter.epidata.tr = fmri.hdr.dime.pixdim(5)*1000;
  
% REFERENCE 1: The cardiac signal
  matlabbatch{1}.spm.tools.drifter.refdata(1).name     = 'Cardiac Signal';
  matlabbatch{1}.spm.tools.drifter.refdata(1).data     = ref_data_puls;
  matlabbatch{1}.spm.tools.drifter.refdata(1).dt       = 1/50;
  matlabbatch{1}.spm.tools.drifter.refdata(1).downdt   = 1/2;
  matlabbatch{1}.spm.tools.drifter.refdata(1).freqlist = 60:120;
  matlabbatch{1}.spm.tools.drifter.refdata(1).N        = 3;
  matlabbatch{1}.spm.tools.drifter.refdata(1).Nimm     = 1;
  
% REFERENCE 2: The respiration signal (see above for details)
  matlabbatch{1}.spm.tools.drifter.refdata(2).name     = 'Respiratory Signal';
  matlabbatch{1}.spm.tools.drifter.refdata(2).data     = ref_data_resp;
  matlabbatch{1}.spm.tools.drifter.refdata(2).dt       = 1/50;
  matlabbatch{1}.spm.tools.drifter.refdata(2).downdt   = 1/2;
  matlabbatch{1}.spm.tools.drifter.refdata(2).freqlist = 10:20;
  matlabbatch{1}.spm.tools.drifter.refdata(2).N        = 4;
  matlabbatch{1}.spm.tools.drifter.refdata(2).Nimm     = 1;
     
% RUN THE JOB (either in interactive mode or as a batch job)
  %data = spm_jobman('interactive',matlabbatch);
  output_list = spm_jobman('run',matlabbatch);
  
end
