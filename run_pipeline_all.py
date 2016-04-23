import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # Utilities
import nipype.interfaces.fsl as fsl	     # fsl
import nipype.interfaces.afni as afni	     # afni
import nipype.pipeline.engine as pe          # pypeline engine
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.interfaces.matlab import MatlabInputSpec, MatlabCommand    # matlab
from nipype.interfaces.nipy.preprocess import Trim
import sklearn
import nibabel
from sklearn.svm import SVR 
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import os
import shutil
import copy
import numpy as np

#TODO: change path
data_path = '/media/Arkisto/hannahalme/data/image_data/potilaat/'
results_path = '/home/hannahalme/Desktop/MTBI/scripts/Biomedicum/results/'

subj_dirs = os.listdir(data_path)

#for k in range(0,len(subj_dirs)):
for k in range(2,3):
	subj = subj_dirs[k]
	try:
	    os.stat(results_path+subj)
	except:
	    os.mkdir(results_path+subj)

	os.chdir(results_path+subj)
	print "Currently processing subject: ", subj

	### ANATOMICAL PREPROCESSING
	"""Nipype-pipeline: anatomical preprocessing

	Anatomical data preprocessing steps:

	1. Transform slices from oblique to axial orientation and to FSL std orientation
	2. Skull strip
	(3. Tissue segmentation, comment the lines out if you wish to do this)
	4. Register to MNI152 standard template with linear and nonlinea registration
	"""

	anat = data_path + subj + '/anat/t1_mpr_sag_p2_iso.nii.gz'

	# Initialize workflow_anat
	workflow_anat = pe.Workflow(name='anat')
	workflow_anat.base_dir = '.'

	#TODO: tarkista polut
	ref_brain = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
	ref_mask = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
	reference_skull = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm.nii.gz'
	fnirt_config = '/usr/share/fsl/5.0/etc/flirtsch/T1_2_MNI152_2mm.cnf'

	# Reorient to FSL standard orientation
	deoblique = pe.Node(interface=afni.Warp(in_file=anat, deoblique=True, outputtype='NIFTI_GZ'), name='deoblique')
	reorient = pe.Node(interface=fsl.Reorient2Std(output_type='NIFTI_GZ'), name='reorient')
	workflow_anat.connect(deoblique, 'out_file', reorient, 'in_file')

	# AFNI skullstrip
	skullstrip = pe.Node(interface=afni.SkullStrip(args='-push_to_edge -ld 30', outputtype='NIFTI_GZ'), name='skullstrip')
	workflow_anat.connect(reorient, 'out_file', skullstrip, 'in_file')

	# Segment with FSL FAST
	#tissue priors
	#tissue_path = '/usr/share/fsl/5.0/data/standard/tissuepriors/2mm/'
	#csf_prior = tissue_path + 'avg152T1_csf_bin.nii.gz'
	#white_prior = tissue_path + 'avg152T1_white_bin.nii.gz'
	#gray_prior = tissue_path + 'avg152T1_gray_bin.nii.gz'

	#segmentation = pe.Node(interface=fsl.FAST(number_classes=3, use_priors=True, img_type=1), name='segmentation')
	#workflow_anat.connect(skullstrip, 'out_file', segmentation, 'in_files')

	# Register to standard MNI template
	#1. linear
	linear_reg = pe.Node(interface=fsl.FLIRT(cost='corratio', reference=ref_brain), name='linear_reg')

	#2.nonlinear
	nonlinear_reg = pe.Node(interface=fsl.FNIRT(fieldcoeff_file=True, jacobian_file=True, ref_file=ref_brain, refmask_file=ref_mask), name='nonlinear_reg')

	inv_flirt_xfm = pe.Node(interface=fsl.utils.ConvertXFM(invert_xfm=True), name='inv_linear_xfm')

	workflow_anat.connect(skullstrip, 'out_file', linear_reg, 'in_file')
	workflow_anat.connect(linear_reg, 'out_matrix_file', nonlinear_reg, 'affine_file')
	workflow_anat.connect(skullstrip, 'out_file', nonlinear_reg, 'in_file')
	workflow_anat.connect(linear_reg, 'out_matrix_file', inv_flirt_xfm, 'in_file')

	# Run workflow_anat
	workflow_anat.write_graph()
	workflow_anat.run()
	print "ANATOMICAL PREPROCESSING DONE! Results in ", results_path+subj


	### FUNCTIONAL PREPROCESSING
	"""
	Nipype-pipeline: fMRI preprocessing

	Functional data preprocessing steps:

	1. 4D realignment: motion and slice timing correction
	2. Reorientation to FSL standard orientation
	3. DRIFTER: physiological signal filtering
	4. Skull stripping
	5. Removing of 3 first volumes
	6. Intensity normalization
	7. Registration to MNI152 standard space

	Workflow1 (result folder func1_1 and func2_1) = steps 1-3
	Workflow2 (result folder func1_2 and func2_2) = steps 4-7
	"""

	for data in (['func1','func2']):
		subj = subj_dirs[k]
		try:
		    os.stat(results_path+subj)
		except:
		    os.mkdir(results_path+subj)

		os.chdir(results_path+subj)
		print "Currently processing subject: ", subj

		# Define input files: 2xfMRI + 1xMPRAGE
		func1 = data_path + subj + '/func/restingstatefMRI1.nii.gz'
		func2 = data_path + subj + '/func/restingstatefMRI2.nii.gz'
		anat = results_path + subj + '/anat/skullstrip/t1_mpr_sag_p2_iso_warp_reoriented_skullstrip.nii.gz'


		# Physiological signal files
		physsig1 = data_path + subj + '/physsig/pulse_fMRI1.mat'
		physsig2 = data_path + subj + '/physsig/resp_fMRI1.mat'
		physsig3 = data_path + subj + '/physsig/pulse_fMRI2.mat'
		physsig4 = data_path + subj + '/physsig/resp_fMRI2.mat'


		#Initialize workflows
		workflow1 = pe.Workflow(name = data + '_1')

		workflow1.base_dir = '.'
		inputnode = pe.Node(interface=util.IdentityInterface(fields=['source_file']), name='inputspec')
		outputnode = pe.Node(interface=util.IdentityInterface(fields=['mean_image', 'move_par']), name='outputspec')

		if data == 'func1':
			inputnode.inputs.source_file = func1
		else:
			inputnode.inputs.source_file = func2		

		# Motion correction + slice timing correction
		realign4d = pe.Node(interface=SpaceTimeRealigner(), name='realign4d')
		realign4d.inputs.ignore_exception=True
		realign4d.inputs.slice_times='asc_alt_siemens'
		realign4d.inputs.slice_info = 2
		realign4d.inputs.tr = 2.00

		workflow1.connect(inputnode, 'source_file', realign4d, 'in_file')
		workflow1.connect(realign4d, 'par_file', outputnode, 'move_par')

		# Reorient
		deoblique = pe.Node(interface=afni.Warp(deoblique=True, outputtype='NIFTI_GZ'), name='deoblique')
		workflow1.connect(realign4d, 'out_file', deoblique, 'in_file')
		reorient = pe.Node(interface=fsl.Reorient2Std(output_type='NIFTI_GZ'), name='reorient')
		workflow1.connect(deoblique, 'out_file', reorient, 'in_file')
		workflow1.connect(reorient, 'out_file', outputnode, 'out_file')
		
		# Run workflow1
		workflow1.write_graph()
		workflow1.run()
		
		# Initialize DRIFTER
		if data=='func1':
			#check is physiological signals exist, if not, then None
			infile = results_path+subj+'/' + data + '_1/reorient/corr_restingstatefMRI1_warp_reoriented.nii.gz'
			physsig_a = physsig1
			if not os.path.isfile(physsig1):
			    	print "MISSING PULS DATA!"
			physsig_b = physsig2
			if not os.path.isfile(physsig2):
				print "MISSING RESP DATA!"

		else:
			infile = results_path+subj+'/' + data + '_1/reorient/corr_restingstatefMRI2_warp_reoriented.nii.gz'
			physsig_a = physsig3
			if not os.path.isfile(physsig3):
			    	print "MISSING PULS DATA!"
			physsig_b = physsig4
			if not os.path.isfile(physsig4):
				print "MISSING RESP DATA!"

		print "Running DRIFTER..."
		script="run_drifter_noSPM_kopio('%s','%s','%s')" %  (infile, physsig_a, physsig_b)
		# TODO: muuta seuraavat polut oikeiksi!!
		# tahan Matlabin asennuspolku
		MatlabCommand.set_default_paths(['/usr/local/MATLAB/R2015a/spm8/toolbox/DRIFTER-toolbox/DRIFTER', '/usr/local/MATLAB/R2015a/spm8/','/usr/local/MATLAB/NIFTI20130306'])
		# TODO: tahan DRIFTER-toolboxin polku
		mlab = MatlabCommand(script=script, mfile=True, paths='/media/HannaHalmeLevy/scripts/DRIFTER-toolbox/DRIFTER/', terminal_output = "stream")

		try:
			os.stat(results_path+subj+'/' + data + '_1/drifter/drifter_corrected.nii.gz')
		except:
			drifter_result = mlab.run()
			os.mkdir(results_path+subj+'/' + data + '_1/drifter/')
			os.rename(results_path+subj+'/drifter_corrected.nii.gz', results_path+subj+'/' + data + '_1/drifter/drifter_corrected.nii.gz')
		
		# Initialize workflow2
		workflow2 = pe.Workflow(name=data + '_2')
		workflow2.base_dir = '.'

		inputnode2 = pe.Node(interface=util.IdentityInterface(fields=['drifter_result']), name='inputspec2')
		outputnode2 = pe.Node(interface=util.IdentityInterface(fields=['result_func']), name='outputspec2')

		inputnode2.inputs.drifter_result=results_path+subj+'/' + data + '_1/drifter/drifter_corrected.nii.gz'

		# AFNI skullstrip and mean image skullstrip
		tstat1 = pe.Node(interface=afni.TStat(args='-mean',outputtype="NIFTI_GZ"), name='tstat1')
		automask = pe.Node(interface=afni.Automask(dilate=1,outputtype="NIFTI_GZ"), name='automask')
		skullstrip = pe.Node(interface=afni.Calc(expr = 'a*b',outputtype="NIFTI_GZ"), name='skullstrip')
		tstat2 = pe.Node(interface=afni.TStat(args='-mean',outputtype="NIFTI_GZ"), name='tstat2')

		workflow2.connect(inputnode2, 'drifter_result', tstat1,'in_file')
		workflow2.connect(tstat1, 'out_file', automask, 'in_file')
		workflow2.connect(automask, 'out_file', skullstrip, 'in_file_b')
		workflow2.connect(inputnode2, 'drifter_result', skullstrip, 'in_file_a')
		workflow2.connect(skullstrip, 'out_file', tstat2, 'in_file')

		# Remove n (3) first volumes
		trim = pe.Node(interface=Trim(begin_index=3), name='trim')
		workflow2.connect(skullstrip, 'out_file', trim, 'in_file')

		# Spatial smoothing, kernel sigma 2.00 mm (5 mm is too much)
		smooth = pe.Node(interface=fsl.maths.SpatialFilter(operation='mean', 			terminal_output='stream', kernel_shape='gauss', kernel_size=1.5, 			nan2zeros=True), name='smooth')

		workflow2.connect(trim,'out_file', smooth, 'in_file')

		# Normalize the median value of each run to 10000
		intnorm = pe.Node(interface=fsl.ImageMaths(op_string='-ing 10000', suffix='_intnorm'), name='intnorm')
		workflow2.connect(smooth, 'out_file', intnorm, 'in_file')

		# Register to standard space
		mean2anat = pe.Node(fsl.FLIRT(bins=40, cost='normmi', dof=12, interp='nearestneighbour'), name='mean2anat')
		mean2anat.inputs.reference = anat
		workflow2.connect(tstat2, 'out_file', mean2anat, 'in_file')

		# Transform mean functional image
		warpmean = pe.Node(interface=fsl.ApplyWarp(), name='warpmean')
		warpmean.inputs.ref_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
		warpmean.inputs.field_file = results_path+subj+'/anat/nonlinear_reg/t1_mpr_sag_p2_iso_warp_reoriented_skullstrip_fieldwarp.nii.gz'
		workflow2.connect(tstat2, 'out_file', warpmean, 'in_file')
		workflow2.connect(mean2anat, 'out_matrix_file', warpmean, 'premat')

		# Transform all images
		warpall = pe.Node(interface=fsl.ApplyWarp(), name='warpall')
		warpall.inputs.ref_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
		warpall.inputs.field_file = results_path+subj+'/anat/nonlinear_reg/t1_mpr_sag_p2_iso_warp_reoriented_skullstrip_fieldwarp.nii.gz'
		workflow2.connect(intnorm, 'out_file', warpall, 'in_file')
		workflow2.connect(mean2anat, 'out_matrix_file', warpall, 'premat')
		workflow2.connect(warpall, 'out_file', outputnode2, 'result_func')

		# Run workflow2
		workflow2.write_graph()
		workflow2.run()

	print "FUNCTIONAL PREPROCESSING DONE! Results in ", results_path+subj


	### SVR 3D
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
	""" 

	print "\n\nSVR analysis, currently processing subject: ", subj

	# Preprocessed fMRI data 
	func1 = data_path + subj + '/func1_2/warpall/drifter_corrected_calc_trim_filt_intnorm_warp.nii.gz'
	func2 = data_path + subj + '/func2_2/warpall/drifter_corrected_calc_trim_filt_intnorm_warp.nii.gz'


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

	#print "Saving images..."
	#corr_img1 = nibabel.Nifti1Image(corr_coef_map1, affine)
	#corr_img2 = nibabel.Nifti1Image(corr_coef_map2, affine)
	#reprod_img = nibabel.Nifti1Image(reprod_coef_map, affine)
	#nibabel.save(corr_img1,results_path + 'subj' + subj + '_corrcoef1.nii.gz')
	#nibabel.save(corr_img2,results_path + 'subj' + subj + '_corrcoef2.nii.gz')
	#nibabel.save(reprod_img, results_path + 'subj' + subj +'_reprod_coef.nii.gz')

	# Save correlations as mat files
	from scipy.io import savemat	
	savemat(results_path + 'subj' + subj + '_corrcoef_fmri1', {'corrcoef_fmri1':corrcoef1})
	savemat(results_path + 'subj' + subj + '_corrcoef_fmri2', {'corrcoef_fmri2':corrcoef2})
	savemat(results_path + 'subj' + subj + '_reprod_coef', {'reprod_coef':reprod_coef})

	
	### CALCULATE STATISTICS
	""" 
	Compare prediction accuracy and reproducibility of a single
	patient to a distribution of healthy controls and detect
	outlier ROIs

	The outliers are presented in a box plot, in which ROIs
	showing accuracy or reproducibility < 3 std below the normal
	distribution are marked with a star and a red text.

	In addition, the ROIs are mapped to the patient's individual
	anatomy (.nii.gz file), which can be visualised with e.g. FSLView
	"""

	pilot_data = '/home/hannahalme/Desktop/MTBI/scripts/Biomedicum/Pilot_data/Preprocessed/'
	pilot_dir = os.listdir(pilot_data)

	# Harvard-Oxford lateralized cortical atlas
	cort_labels = np.loadtxt('HO_labels.txt', delimiter=',', usecols=(4,),dtype=str)

	print "Testing outliers for subject", subj
	ROI_corr1 = loadmat(results_path+'subj' + subj + '_corrcoef_fmri1.mat')['corrcoef_fmri1'][0]
	ROI_corr2 = loadmat(results_path+'subj' + subj + '_corrcoef_fmri2.mat')['corrcoef_fmri2'][0]
	reprod = loadmat(results_path+'subj' + subj +'_reprod_coef.mat')['reprod_coef'][0]

	# Average prediction accuracies for fMRI1 and fMRI2
	test_accuracy = (ROI_corr1 + ROI_corr2) / 2.0
	test_reproducibility = reprod
	anat = patient_results + subj + '/anat/skullstrip/t1_mpr_sag_p2_iso_warp_reoriented_skullstrip.nii.gz'

	# Number of healthy control subjects
	n=len(pilot_dir)
	roin=len(cort_labels)
	cort_ROI_list1 = np.zeros([roin,n])
	cort_ROI_list2 = np.zeros([roin,n])

	# Distributions of prediction accuracy & reproducibility for healthy subjects
	for i in range(0,n):
		subj = pilot_dir[i]
		print "Currently processing pilot subject: ", subj

		ROI_corr1 = loadmat(results_path+'Drifter_SVR/subj' + subj + '_corrcoef_fmri1.mat')['corrcoef_fmri1'][0]
		ROI_corr2 = loadmat(results_path+'Drifter_SVR/subj' + subj + '_corrcoef_fmri2.mat')['corrcoef_fmri2'][0]
		reprod = loadmat(results_path+'Drifter_SVR/subj' + subj +'_reprod_coef.mat')['reprod_coef'][0]

		# Accuracy & reproducibility
		cort_ROI_list1[:,i] = (ROI_corr1 + ROI_corr2) / 2.0
		cort_ROI_list2[:,i] = reprod

	# sort from largest to smallest
	accuracy = np.nan_to_num(cort_ROI_list1)
	ids1 = np.argsort(np.mean(accuracy,axis=1))
	accuracy, cort_labels = accuracy[ids1], cort_labels[ids1]
	reproducibility = cort_ROI_list2
	reproducibility = reproducibility[ids1]

	# Calculate & detect outliers
	diff_acc = (test_accuracy-np.mean(accuracy,axis=1))/np.std(accuracy,axis=1)
	diff_rep = (test_reproducibility-np.mean(reproducibility,axis=1))/np.std(reproducibility,axis=1)

	# Make box & whiskers plots for each Harvard-Oxford region
	plt.figure(figsize=(10,20))
	plt.boxplot(accuracy.T,vert=0,showmeans=True)
	ticks, labels = plt.yticks(np.arange(0,np.shape(cort_labels)[0]), cort_labels, fontsize=12)
	for i in range(0,roin-1):
	    if diff_acc[i] < -3.0:
	    	plt.scatter(test_accuracy[i], i, marker="*")
		labels[i].set_color('red') 
	plt.xlabel('Prediction Accuracy')
	plt.title('Accuracy for patient '+subj)
	plt.gcf().subplots_adjust(left=0.50)

	# Save figure as pdf
	plt.savefig(results_path + subj + '_accuracy.pdf')

	# Plot outlier regions to a brain map in the patient's anatomical space
	print "Plotting accuracy outlier ROIs to a brain map..."
	#TODO: tarkista polku
	ROI_atlas = nibabel.load('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr25-2mm.nii.gz')
	orig_ROI_data = ROI_atlas.get_data()
	affine = ROI_atlas.get_affine()
	acc_roi_data = np.zeros(np.shape(orig_ROI_data))
	for i in range(1,roin):
		ids = np.where(orig_ROI_data==i)
		acc_roi_data[ids[0],ids[1], ids[2]] = diff_acc[i]
	acc_roi_data[np.where(acc_roi_data>0)] = 1 #remove positive values
	acc_img = nibabel.Nifti1Image(-acc_roi_data, affine)#change sign to positive
	nibabel.save(acc_img,results_path+'subj' + subj + '_diff_accuracy.nii.gz')

	# Transform to native anatomical space
	print "Transforming to native space..."
	flirt = fsl.FLIRT(bins=800, cost='mutualinfo', reference=anat)
	flirt.inputs.in_file = results_path+'subj' + subj + '_diff_accuracy.nii.gz'
	flirt.inputs.out_file = results_path+'subj' + subj + '_diff_accuracy_transformed.nii.gz'
	flirt.run()

	# Make box & whiskers plots for each Harvard-Oxford region
	plt.figure(figsize=(10,20))
	plt.boxplot(reproducibility.T,vert=0,showmeans=True)
	ticks, labels = plt.yticks(np.arange(0,np.shape(cort_labels)[0]), cort_labels, fontsize=12)
	for i in range(0,roin-1):
	    if diff_rep[i] < -3.0:
	    	plt.scatter(test_reproducibility[i], i, marker="*")
		labels[i].set_color('red')
	plt.xlabel('Reproducibility')
	plt.title('Reproducibility for patient '+subj)
	plt.gcf().subplots_adjust(left=0.50)

	# Save figure as pdf
	plt.savefig('results/'+ subj + '_reproducibility.pdf')

	# Plot outlier regions to a brain map in the patient's anatomical space
	print "Plotting reproducibility outlier ROIs to a brain map..."
	rep_roi_data = np.zeros(np.shape(orig_ROI_data))
	for i in range(1,roin):
		ids = np.where(orig_ROI_data==i)
		rep_roi_data[ids[0],ids[1], ids[2]] = diff_rep[i]
	rep_roi_data[np.where(rep_roi_data>0)] = 1 #remove positive values
	rep_img = nibabel.Nifti1Image(-rep_roi_data, affine)#change sign to positive
	nibabel.save(rep_img,results_path+'subj' + subj + '_diff_reproducibility.nii.gz')

	# Transform to native anatomical space
	print "Transforming to native space..."
	flirt.inputs.in_file = results_path+'subj' + subj + '_diff_reproducibility.nii.gz'
	flirt.inputs.out_file = results_path+'subj' + subj + '_diff_reproducibility_transformed.nii.gz'
	flirt.run()

	# Remove FLIRT mat files
	os.remove('subj' + subj + '_diff_accuracy_flirt.mat')
	os.remove('subj' + subj + '_diff_reproducibility_flirt.mat')

print "ALL DONE!"
