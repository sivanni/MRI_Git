""" 

Resting state fMRI connectivity analysis with Support Vector Regression (SVR)

This script implements the method described in [1].

Overview of the method:
1. Parcellate the brain into ROIs according to Harvard-Oxford atlas
2. On each iteration, remove one ROI from the data 1 and generate a linear
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
data_path = '/opt/Laskenta/Control_Room/Biomedicum/results/'

# change the subject directory!!
subj_dirs = ['20150407']

for k in range(0,len(subj_dirs)):
	subj = subj_dirs[k]

	print "\n\nCurrently processing subject: ", subj

	# Preprocessed fMRI data 
	func1 = data_path + subj + '/func1_2/warpall/drifter_corrected_calc_trim_intnorm_warp.nii.gz'
	func2 = data_path + subj + '/func2_2/warpall/drifter_corrected_calc_trim_intnorm_warp.nii.gz'


	mask = nibabel.load('/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz').get_data()

	print "Collecting timeseries data..."
	func1_data = nibabel.load(func1).get_data()
	func2_data = nibabel.load(func2).get_data()
	volnum = np.shape(func1_data)[-1]
	
	print "Running support vector regression..."
	# Number of ROIs in Harvard-Oxford atlas
	roin=96

	#linear kernel & other parameters according to Craddock's paper
	svr = SVR(kernel="linear", C=100.0, epsilon=0.1, max_iter=30000)
	scaler = StandardScaler()
	pipe = Pipeline([('scaler', scaler),('clf',svr)])

	corrcoef1 = np.empty([roin])
	corrcoef2 = np.empty([roin])

	reprod_coef = np.empty([roin])

	for i in range(0,roin):
		# Reload ROI data on every iteration
		ROI_data = nibabel.load('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz').get_data()

		# Define voxels belonging to a current ROI and set to zero
		ids = np.where(ROI_data==i+1)
		ROI_data[ids[0],ids[1], ids[2]] = 0
		ids2 = np.where(ROI_data>0)
	
		roi_timeseries1 = np.mean(func1_data[ids[0],ids[1], ids[2], :], axis=0)
		roi_timeseries2 = np.mean(func2_data[ids[0],ids[1], ids[2], :], axis=0)
	
		# Remove background and ROI voxels from training data
		brain_timeseries1 = func1_data[ids2[0], ids2[1], ids2[2], :]
		brain_timeseries2 = func2_data[ids2[0], ids2[1], ids2[2], :]

		# Remaining ROI timeseries in fMRI2 = true timeseries	
		roi_timeseries1 =  (roi_timeseries1-np.mean(roi_timeseries1))/np.std(roi_timeseries1)
		roi_timeseries2 =  (roi_timeseries2-np.mean(roi_timeseries2))/np.std(roi_timeseries2)

		# Predict timeseries with SVR	
		pred_labels1 = pipe.fit(brain_timeseries1.T, roi_timeseries1).predict(brain_timeseries2.T)
		model1 = svr.coef_
		pred_labels1 =  (pred_labels1-np.mean(pred_labels1))/np.std(pred_labels1)

		# Same thing again, now use the other data for training...
		pred_labels2 = pipe.fit(brain_timeseries2.T, roi_timeseries2).predict(brain_timeseries1.T)
		model2 = svr.coef_
		pred_labels2 =  (pred_labels2-np.mean(pred_labels2))/np.std(pred_labels2)

		# Add to predicted timeseries 1
		corrcoef1[i] = np.corrcoef(pred_labels1, roi_timeseries2)[0,1]
		print i, '1: ', corrcoef1[i]
		corrcoef2[i] = np.corrcoef(pred_labels2, roi_timeseries1)[0,1]
		print i, '2: ', corrcoef2[i]
		reprod_coef[i] = np.corrcoef(model1, model2)[0,1]	
		print reprod_coef[i]

		model1, model2, ROI_data = [], [], []

	# Save results

	# Map correlation coefficients to a brain atlas
	print "Done! Mapping correlation coefficient to a brain volume..."
	ROI_atlas = nibabel.load('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz')
	ROI_data = ROI_atlas.get_data()
	affine = ROI_atlas.get_affine()
	corr_coef_map1 = np.zeros(np.shape(ROI_data))
	corr_coef_map2 = np.zeros(np.shape(ROI_data))
	reprod_coef_map = np.zeros(np.shape(ROI_data))

	for j in range(0,roin-1):
		ids = np.where(ROI_data==j+1)
		corr_coef_map1[ids] = corrcoef1[j]
		corr_coef_map2[ids] = corrcoef2[j]
		reprod_coef_map[ids] = reprod_coef[j]

	print "Saving images..."
	corr_img1 = nibabel.Nifti1Image(corr_coef_map1, affine)
	corr_img2 = nibabel.Nifti1Image(corr_coef_map2, affine)
	reprod_img = nibabel.Nifti1Image(reprod_coef_map, affine)
	nibabel.save(corr_img1,'results/subj' + subj + '_corrcoef1.nii.gz')
	nibabel.save(corr_img2,'results/subj' + subj + '_corrcoef2.nii.gz')
	nibabel.save(reprod_img,'results/subj' + subj +'_reprod_coef.nii.gz')

	# Save correlations as mat files
	from scipy.io import savemat	
	savemat('results/subj' + subj + '_corrcoef_fmri1', {'corrcoef_fmri1':corrcoef1})
	savemat('results/subj' + subj + '_corrcoef_fmri2', {'corrcoef_fmri2':corrcoef2})
	savemat('results/subj' + subj + '_reprod_coef', {'reprod_coef':reprod_coef})
