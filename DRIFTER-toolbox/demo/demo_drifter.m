%% DEMO_drifter - Example usage of the DRIFTER toolbox
%
% Description:
%   This file demonstrates the use of the DRIFTER method in SPM. It 
%   uses data files that are available for download on the toolbox 
%   web page. Refer to the toolbox documentation for further details.
%
% See also:
%   http://becs.aalto.fi/en/research/bayes/drifter/
%
% Version:
%   updated 2012-04-17
%
% Copyright:
%   Arno Solin, 2011

%% Set up

  % Make sure SPM is in your path (Must be changed!!!)
  addpath /usr/local/MATLAB/R2015a/spm8
  
  % Initialize SPM defaults
  spm('Defaults','fMRI');
  spm_jobman('initcfg'); % useful in SPM8 only
  
  % Load filenames using SPM built-in functions
  f = spm_select('FPListRec','demodata/','^d.*\.img$');
  
  % Load reference signal data (respiration, cardiac), 
  % both sampled at dt=0.001 s
  load references.mat
  
  
%% Preliminary stuff

  % Here you can do some preliminary fixing of the data. Such as
  % - REALIGN
  % - COREGISTER
  % - SEGMENT
  % - ... (not necessarily prior to using DRIFTER)

  % Help on these can be found e.g. in the SPM documentation file.
  
  
%% Use the new toolbox

% METHOD SPECIFICATION

  % Choose in which mode to run the method (0 = DRIFTER with all noise 
  % removed, 1 = DRIFTER with only physiological noise removed, 
  % 2 = a RETROICOR implementation)
  jobs{1}.spm.tools.drifter.mode = 1;
  
  % Choose with which prefix to output the data
  jobs{1}.spm.tools.drifter.prefix = 'n';
      
  % Assign the EPI fMRI data files (as we have not done anything prior to 
  % this, we use the original data files loaded into f)
  jobs{1}.spm.tools.drifter.epidata.files = f;
  
  % The  Interscan interval, TR, (specified in milliseconds)
  jobs{1}.spm.tools.drifter.epidata.tr = 100;
  
% REFERENCE 1: The cardiac signal
    
  % Signal name (for debugging)
  jobs{1}.spm.tools.drifter.refdata(1).name     = 'Cardiac Signal';
    
  % Assign the reference data that was loaded earlier
  jobs{1}.spm.tools.drifter.refdata(1).data     = cardiac;
  
  % Set the sampling interval of the reference signal (dt=0.001)
  jobs{1}.spm.tools.drifter.refdata(1).dt       = 1/1000;

  % To speed up the estimation we downsample the signal
  jobs{1}.spm.tools.drifter.refdata(1).downdt   = 1/10;
  
  % List of possible frequencies for this phenomenon (in beats per min)
  jobs{1}.spm.tools.drifter.refdata(1).freqlist = 60:120;

  % Number of periodics to estimate (fundamental + number of harmonics)
  jobs{1}.spm.tools.drifter.refdata(1).N        = 3;

  % There is no need to estimate as many periodics while finding the
  % frequency. Therefore we use only the fundamental here.
  jobs{1}.spm.tools.drifter.refdata(1).Nimm     = 1;
  
% REFERENCE 2: The respiration signal (see above for details)
  jobs{1}.spm.tools.drifter.refdata(2).name     = 'Respiratory Signal';
  jobs{1}.spm.tools.drifter.refdata(2).data     = respiration;
  jobs{1}.spm.tools.drifter.refdata(2).dt       = 1/1000;
  jobs{1}.spm.tools.drifter.refdata(2).downdt   = 1/10;
  jobs{1}.spm.tools.drifter.refdata(2).freqlist = 10:70;
  jobs{1}.spm.tools.drifter.refdata(2).N        = 4;
  jobs{1}.spm.tools.drifter.refdata(2).Nimm     = 1;
     
% RUN THE JOB (either in interactive mode or as a batch job)
  %spm_jobman('interactive',jobs);
  spm_jobman('run',jobs);
  
  
%% Further analysis

  % Now you have the cleaned data stored on disk with the prefix specified
  % above. Hereafter the data can be analyzed in various ways with e.g. the
  % default functionalities provided by SPM.
  % See the SPM documentation for help.

  