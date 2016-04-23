fprintf(1,'Executing %s at %s:\n',mfilename,datestr(now));
ver,
try,run_drifter_noSPM_kopio('/home/hannahalme/Desktop/MTBI/results/20141108_1/motion_corrected.nii.gz','/home/hannahalme/Desktop/MTBI/data/image_data/meilahti/20141108_1/physsig/puls1.puls.mat', '/home/hannahalme/Desktop/MTBI/data/image_data/meilahti/20141108_1/physsig/resp1.resp.mat')
,catch ME,
fprintf(2,'MATLAB code threw an exception:\n');
fprintf(2,'%s\n',ME.message);
if length(ME.stack) ~= 0, fprintf(2,'File:%s\nName:%s\nLine:%d\n',ME.stack.file,ME.stack.name,ME.stack.line);, end;
end;