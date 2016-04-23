""" 

Compare prediction accuracy and reproducibility of a single
patient to a distribution of healthy controls and detect
outlier ROIs

The outliers are presented in a box plot, in which ROIs
showing accuracy or reproducibility < 3 std below the normal
distribution are marked with a star and a red text.

In addition, the ROIs are mapped to the patient's individual
anatomy (.nii.gz file), which can be visualised with e.g. FSLView

Hanna Halme // 2015
Questions and comments: hanna.halme@hus.fi

"""
import nipype.interfaces.io as nio           
import os
import numpy as np
import nibabel
import nipype.interfaces.fsl as fsl	    
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt
import os

# TODO: muuta tahan polut joissa potilaiden ja verrokkien analyysien tulokset
# esikasitelty potilasdata
patient_data = '/media/HannaHalmeLevy/results/Patient_data/Preprocessed_DRIFTER/'
# potilaiden SVR-analyysin tulokset
patient_results = '/media/HannaHalmeLevy/results/Patient_data/DRIFTER_SVR/'
# esikasitelty verrokkidata
pilot_data =  '/media/HannaHalmeLevy/results/Pilot_data/'
# verrokkien SVR-analyysin tuloksey
pilot_results = '/media/HannaHalmeLevy/results/Pilot_data/Preprocessed_DRIFTER/'

pilot_dir = os.listdir(pilot_results)

# Load patient data
subjects = os.listdir(patient_data)

for i in range(0, len(subjects)):

	test_subj = subjects[i]

	# Harvard-Oxford lateralized cortical atlas
	cort_labels = np.loadtxt('HO_labels.txt', delimiter=',', usecols=(4,),dtype=str)

	print "Testing outliers for subject", test_subj
	ROI_corr1 = loadmat(patient_results + 'subj' + test_subj + '_corrcoef_fmri1.mat')['corrcoef_fmri1'][0]
	ROI_corr2 = loadmat(patient_results + 'subj' + test_subj + '_corrcoef_fmri2.mat')['corrcoef_fmri2'][0]
	reprod = loadmat(patient_results  + 'subj' + test_subj +'_reprod_coef.mat')['reprod_coef'][0]

	# Average prediction accuracies for fMRI1 and fMRI2
	test_accuracy = (ROI_corr1 + ROI_corr2) / 2.0
	test_reproducibility = reprod
	anat = patient_data + test_subj + '/anat/skullstrip/t1_mpr_sag_p2_iso_warp_reoriented_skullstrip.nii.gz'

	# Number of healthy control subjects
	n=len(pilot_dir)
	roin=len(cort_labels)
	cort_ROI_list1 = np.zeros([roin,n])
	cort_ROI_list2 = np.zeros([roin,n])

	# Distributions of prediction accuracy & reproducibility for healthy subjects
	for i in range(0,n):

		subj = pilot_dir[i]
		print "Currently processing pilot subject: ", subj

		ROI_corr1 = loadmat(pilot_data + 'DRIFTER_SVR/subj' + subj + '_corrcoef_fmri1.mat')['corrcoef_fmri1'][0]
		ROI_corr2 = loadmat(pilot_data + 'DRIFTER_SVR/subj' + subj + '_corrcoef_fmri2.mat')['corrcoef_fmri2'][0]
		reprod = loadmat(pilot_data + 'DRIFTER_SVR/subj' + subj +'_reprod_coef.mat')['reprod_coef'][0]

		# Accuracy & reproducibility
		cort_ROI_list1[:,i] = (ROI_corr1 + ROI_corr2) / 2.0
		cort_ROI_list2[:,i] = reprod

	# sort from largest to smallest
	accuracy = np.nan_to_num(cort_ROI_list1)
	ids1 = np.argsort(np.mean(accuracy,axis=1))
	accuracy, cort_labels = accuracy[ids1], cort_labels[ids1]
	reproducibility = cort_ROI_list2
	reproducibility = reproducibility[ids1]
	test_accuracy = test_accuracy[ids1]
	test_reproducibility = test_reproducibility[ids1]

	# Calculate & detect outliers
	diff_acc = (test_accuracy-np.mean(accuracy,axis=1))/np.std(accuracy,axis=1)
	diff_rep = (test_reproducibility-np.mean(reproducibility,axis=1))/np.std(reproducibility,axis=1)

	# Make box & whiskers plots for each Harvard-Oxford region
	plt.figure(figsize=(10,20))
	plt.boxplot(accuracy.T,vert=0,showmeans=True)
	ticks, labels = plt.yticks(np.arange(1,np.shape(cort_labels)[0]+1), cort_labels, fontsize=12)
	for i in range(0,roin-1):
	    if diff_acc[i] < -3.0:
	    	plt.scatter(test_accuracy[i], i, marker="*")
		labels[i-1].set_color('red') 
	plt.xlabel('Prediction Accuracy')
	plt.title('Accuracy for patient '+test_subj)
	plt.gcf().subplots_adjust(left=0.50)

	# Save figure as pdf
	plt.savefig(patient_results + test_subj + '_accuracy.pdf')

	# Plot outlier regions to a brain map in the patient's anatomical space
	print "Plotting accuracy outlier ROIs to a brain map..."
	ROI_atlas = nibabel.load('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz')
	orig_ROI_data = ROI_atlas.get_data()
	affine = ROI_atlas.get_affine()
	acc_roi_data = np.zeros(np.shape(orig_ROI_data))
	for i in range(1,roin):
		ids = np.where(orig_ROI_data==i)
		acc_roi_data[ids[0],ids[1], ids[2]] = diff_acc[i]
	acc_roi_data[np.where(acc_roi_data>0)] = 1 #remove positive values
	acc_img = nibabel.Nifti1Image(-acc_roi_data, affine)#change sign to positive
	nibabel.save(acc_img, patient_results + 'subj' + test_subj + '_diff_accuracy.nii.gz')

	# Transform to native anatomical space
	print "Transforming to native space..."
	flirt = fsl.FLIRT(bins=800, cost='mutualinfo', reference=anat)
	flirt.inputs.in_file = patient_results + 'subj' + test_subj + '_diff_accuracy.nii.gz'
	flirt.inputs.out_file = patient_results + 'subj' + test_subj + '_diff_accuracy_transformed.nii.gz'
	flirt.run()

	# Make box & whiskers plots for each Harvard-Oxford region
	plt.figure(figsize=(10,20))
	plt.boxplot(reproducibility.T,vert=0,showmeans=True)
	ticks, labels = plt.yticks(np.arange(1,np.shape(cort_labels)[0]+1), cort_labels, fontsize=12)
	for i in range(0,roin-1):
	    if diff_rep[i] < -3.0:
	    	plt.scatter(test_reproducibility[i], i, marker="*")
		labels[i-1].set_color('red')
	plt.xlabel('Reproducibility')
	plt.title('Reproducibility for patient '+test_subj)
	plt.gcf().subplots_adjust(left=0.50)

	# Save figure as pdf
	plt.savefig(patient_results +  test_subj + '_reproducibility.pdf')

	# Plot outlier regions to a brain map in the patient's anatomical space
	print "Plotting reproducibility outlier ROIs to a brain map..."
	rep_roi_data = np.zeros(np.shape(orig_ROI_data))
	for i in range(1,roin):
		ids = np.where(orig_ROI_data==i)
		rep_roi_data[ids[0],ids[1], ids[2]] = diff_rep[i]
	rep_roi_data[np.where(rep_roi_data>0)] = 1 #remove positive values
	rep_img = nibabel.Nifti1Image(-rep_roi_data, affine)#change sign to positive
	nibabel.save(rep_img, patient_results + 'subj' + test_subj + '_diff_reproducibility.nii.gz')

	# Transform to native anatomical space
	print "Transforming to native space..."
	flirt.inputs.in_file = patient_results + 'subj' + test_subj + '_diff_reproducibility.nii.gz'
	flirt.inputs.out_file = patient_results + 'subj' + test_subj + '_diff_reproducibility_transformed.nii.gz'
	flirt.run()


	print "DONE!"
