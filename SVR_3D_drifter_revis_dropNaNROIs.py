""" 

Resting state fMRI connectivity analysis with Support Vector Regression (SVR)

This script implements the method described in [1].

Overview of the method:
1. Parcellate the brain into ROIs according to Harvard-Oxford atlas
2. On each iteration, remove one ROI from data 1 and generate a linear
model based on timeseries of the remaining voxels
3. Predict the timeseries of the removed ROI with the linear model
4. Calculate correlation between the predicted timeseries and true timeseries in data 2
5. Repeat 1-4 using data 2 for training and data 1 for testing
6. Calculate correlation between model 1 and model 2 estimated from independent datasets
7. Map correlation and reprodicibility coefficients to ROIs in standard brain anatomy


Hanna Halme // 2015
Questions and comments: hanna.halme@hus.fi

[1] Craddock, R.C., Milham, M.P., & LaConte, S.M. (2013). Predicting intrinsic
brain activity. Neuroimage, 82, 127-136.

"""
import nipype.interfaces.fsl as fsl	    
import os
import numpy as np
import sklearn
import nibabel
from sklearn.svm import SVR 
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# Results from preprocessing workflows
# data_path = '/opt/MR_data/Silja_Raty/Revis_0007_rs/Results_nodeoblique/Functional_nosearch/'
# results_path = '/opt/MR_data/Silja_Raty/Revis_0007_rs/Results_nodeoblique/Connectivity_analysis/'

subj_dirs = ['Revis_0022_rs']  # os.listdir(data_path)

# Process all patients in a directory
for k in range(0, len(subj_dirs)):
    for s in (['1', '2', '3']):
        subj = subj_dirs[k]
        data_path = '/opt/MR_data/Silja_Raty/'+subj+'/Results_nodeoblique/Functional_nosearch/data'+s+'/'
        results_path = '/opt/MR_data/Silja_Raty/'+subj+'/Results_nodeoblique/Connectivity_analysis_noNaN/data'+s+'/'
        print "\n\nCurrently processing subject: ", subj

        # Preprocessed fMRI data
        func1 = data_path + 'func1_2/warpall/drifter_corrected_calc_trim_intnorm_warp.nii.gz'
        func2 = data_path + 'func2_2/warpall/drifter_corrected_calc_trim_intnorm_warp.nii.gz'

        mask = nibabel.load('/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz').get_data()

        print "Collecting timeseries data..."
        func1_data = nibabel.load(func1).get_data()
        func2_data = nibabel.load(func2).get_data()
        volnum = np.shape(func1_data)[-1]

        print "Running support vector regression..."
        # Number of ROIs in Harvard-Oxford atlas
        roin = 96

        # Linear kernel & other parameters according to Craddock's paper
        svr = SVR(kernel="linear", C=100.0, epsilon=0.1, max_iter=30000)
        scaler = StandardScaler()
        pipe = Pipeline([('scaler', scaler), ('clf', svr)])

        # These were np.empty([roin]) -- potentially dangerous, sets existing memory value to the array
        corrcoef1 = np.zeros([roin])
        corrcoef2 = np.zeros([roin])

        reprod_coef = np.zeros([roin])

        # Create array of ROIs excluding NaN ROIs in any subject
        all_roi_indeces = range(0, roin)
        nanROIs = np.array([15, 16, 17, 21, 22, 27, 28, 29, 30, 49, 50, 53, 54, 65, 66, 67, 68, 73, 74, 76])
        roi_indeces_with_no_nans = np.delete(all_roi_indeces, nanROIs - 1)

        # for roi_index in range(0,roin):
        for roi_index in roi_indeces_with_no_nans:
            # Reload ROI data on every iteration
            ROI_data = nibabel.load('/usr/share/fsl/5.0/data/atlases/'
                                    'HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz').get_data()

            # Define voxels belonging to a current ROI and set to zero
            ids = np.where(ROI_data == roi_index+1)
            ROI_data[ids[0], ids[1], ids[2]] = 0

            # Find the voxels belonging to ROIs containing NaN values, and put them to zero
            for roi_to_remove in nanROIs:
                ROI_data[ROI_data == roi_to_remove] = 0

            # Find indices for voxels in ROIs to be actually used to predict the target ROI time sequence
            ids2 = np.where(ROI_data > 0)

            roi_timeseries1 = np.mean(func1_data[ids[0], ids[1], ids[2], :], axis=0)
            roi_timeseries2 = np.mean(func2_data[ids[0], ids[1], ids[2], :], axis=0)

            # Remove background and ROI voxels from training data
            brain_timeseries1 = func1_data[ids2[0], ids2[1], ids2[2], :]
            brain_timeseries2 = func2_data[ids2[0], ids2[1], ids2[2], :]

            # Remaining ROI timeseries in fMRI2 = true timeseries
            roi_timeseries1 = (roi_timeseries1-np.mean(roi_timeseries1)) / np.std(roi_timeseries1)
            roi_timeseries2 = (roi_timeseries2-np.mean(roi_timeseries2)) / np.std(roi_timeseries2)

            # Predict timeseries with SVR
            pred_labels1 = pipe.fit(brain_timeseries1.T, roi_timeseries1).predict(brain_timeseries2.T)
            model1 = svr.coef_
            pred_labels1 = (pred_labels1-np.mean(pred_labels1))/np.std(pred_labels1)

            # Same thing again, now use the other data for training...
            pred_labels2 = pipe.fit(brain_timeseries2.T, roi_timeseries2).predict(brain_timeseries1.T)
            model2 = svr.coef_
            pred_labels2 = (pred_labels2-np.mean(pred_labels2))/np.std(pred_labels2)

            # Add to predicted timeseries 1
            corrcoef1[roi_index] = np.corrcoef(pred_labels1, roi_timeseries2)[0, 1]
            print roi_index, '1: ', corrcoef1[roi_index]
            corrcoef2[roi_index] = np.corrcoef(pred_labels2, roi_timeseries1)[0, 1]
            print roi_index, '2: ', corrcoef2[roi_index]
            reprod_coef[roi_index] = np.corrcoef(model1, model2)[0, 1]
            print reprod_coef[roi_index]

            model1, model2, ROI_data = [], [], []

        # Save results

        # Map correlation coefficients to a brain atlas
        print "Done! Mapping correlation coefficient to a brain volume..."
        ROI_atlas = nibabel.load('/usr/share/fsl/5.0/data/atlases/HarvardOxford/'
                                 'HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz')
        ROI_data = ROI_atlas.get_data()
        affine = ROI_atlas.get_affine()
        corr_coef_map1 = np.zeros(np.shape(ROI_data))
        corr_coef_map2 = np.zeros(np.shape(ROI_data))
        reprod_coef_map = np.zeros(np.shape(ROI_data))

        for j in range(0, roin-1):
            ids = np.where(ROI_data == j + 1)
            corr_coef_map1[ids] = corrcoef1[j]
            corr_coef_map2[ids] = corrcoef2[j]
            reprod_coef_map[ids] = reprod_coef[j]

        # print "Saving images..."
        # corr_img1 = nibabel.Nifti1Image(corr_coef_map1, affine)
        # corr_img2 = nibabel.Nifti1Image(corr_coef_map2, affine)
        # reprod_img = nibabel.Nifti1Image(reprod_coef_map, affine)
        # nibabel.save(corr_img1,results_path + 'subj' + subj + '_corrcoef1.nii.gz')
        # nibabel.save(corr_img2,results_path + 'subj' + subj + '_corrcoef2.nii.gz')
        # nibabel.save(reprod_img, results_path + 'subj' + subj +'_reprod_coef.nii.gz')

        # Save correlations as mat files
        from scipy.io import savemat
        savemat(results_path + 'subj' + subj + '_corrcoef_fmri1', {'corrcoef_fmri1': corrcoef1})
        savemat(results_path + 'subj' + subj + '_corrcoef_fmri2', {'corrcoef_fmri2': corrcoef2})
        savemat(results_path + 'subj' + subj + '_reprod_coef', {'reprod_coef': reprod_coef})

        # NaN-alueet 15, 16, 17, 21, 22, 27, 28, 29, 30, 49, 50, 53, 54, 65, 66, 67, 68, 73, 74, 76