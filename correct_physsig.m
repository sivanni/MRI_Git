%poistetaan 5000-koodilla merkityt synkkapulssien ajoitukset fysiologisista
%signaaleista...

clear all
close all

data_dir = '/home/hannahalme/Desktop/MTBI/data/image_data/meilahti/';
a = dir(fullfile(data_dir));

for i=3:size(a,1)
    subj = a(i,1).name;
    subj_dir = [data_dir,subj, '/physsig/'];
    filenames = ['puls1.puls.mat'; 'puls2.puls.mat'; 'resp1.resp.mat'; 'resp2.resp.mat'];
    for j=1:size(filenames,1)
        fname = [subj_dir,filenames(j,:)];
        if exist(fname, 'file') == 2
            sig = load(fname);
            sig = sig.A;
            first_pulse = find(sig==5000,1);
            %remove data until first pulse
            sig(1:first_pulse) = [];
            %remove pulse timings
            sig(find(sig==5000)) = [];
            save([subj_dir, filenames(j,1:5),'.mat'], 'sig')
        end
    end
end
