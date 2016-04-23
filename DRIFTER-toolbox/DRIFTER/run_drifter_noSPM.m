
  %fMRI = '/home/hannahalme/Desktop/MTBI/results/20141108_1/motion_corrected.nii.gz'
  %ref_data_puls =  '/home/hannahalme/Desktop/MTBI/data/image_data/meilahti/20141108_1/physsig/puls1.puls.mat'
  %ref_data_resp = '/home/hannahalme/Desktop/MTBI/data/image_data/meilahti/20141108_1/physsig/resp1.resp.mat'
  ref_data_puls = load(ref_data_puls);
  ref_data_puls = ref_data_puls.A;
  ref_data_resp = load(ref_data_resp);
  ref_data_resp = ref_data_resp.A;
          
  fmri = load_untouch_nii(fMRI);
  obs_data = double(fmri.img);
  
  % Provide two structures: "data" and "refdata" such that
  data.data = obs_data;
  data.dt = fmri.hdr.dime.pixdim(5)*1000;

  % And for each reference signal
  refdata{1}.dt = 1/50;
  refdata{1}.freqlist = 60:90; % Vector of possible frequencies in bpm
  refdata{1}.data = ref_data_puls;
  refdata{1}.N = 0;
  refdata{1}.Nimm = 0;
  refdata{1}.downdt = 1.0;

  refdata{2}.dt = 1/50;
  refdata{2}.freqlist = 10:20; % Vector of possible frequencies in bpm
  refdata{2}.data = ref_data_resp;  
  refdata{2}.N = 0;
  refdata{2}.Nimm = 0; 
  refdata{2}.downdt = 1.0;
  
  % Run DRIFTER
  [data,refdata] = drifter(data,refdata);
  save_nii(make_nii(data.data),'/home/hannahalme/Desktop/MTBI/results/20141108_1/drifter_corrected1.nii.gz')
  save_nii(make_nii(data.noise),'/home/hannahalme/Desktop/MTBI/results/20141108_1/drifter_noise1.nii.gz')
 % save_nii(make_nii(data.data),out_fmri)
 % save_nii(make_nii(data.noise),out_noise)

