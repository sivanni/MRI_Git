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

Running in terminal:
ipython
(in python terminal):
%run preprocess_functional_data

Hanna Halme // 2015
Questions and comments: hanna.halme@hus.fi

"""

import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # Utilities
import nipype.interfaces.fsl as fsl	     # fsl
import nipype.interfaces.afni as afni	     # afni
import nipype.pipeline.engine as pe          # pypeline engine
import nipy
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.interfaces.matlab import MatlabInputSpec, MatlabCommand    # matlab
from nipype.interfaces.nipy.preprocess import Trim
import os
import shutil
import copy
import numpy as np

########################################################################################
data_path = '/opt/dicom_images/'
results_path = '/opt/Laskenta/Control_Room/Biomedicum/results/'

#TODO: selvita fysiologisten signaalien muuntaminen!!!

subj_dirs = os.listdir(data_path)
for k in range(0,len(subj_dirs)):

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

		#topup files: b0 images with opposing polarities & encoding txt file
		#b0_LR = data_path + subj + '/func/b0_LR.nii.gz'
		#b0_RL = data_path + subj + '/func/b0_RL.nii.gz'
		#enc_file = '/home/hannahalme/Desktop/MTBI/data/acqparam.txt'
		#enc_file2 = '/home/hannahalme/Desktop/MTBI/data/acqparam2.txt'

		# Physiological signal files
		physsig1 = data_path + subj + '/physsig/puls1.mat'
		physsig2 = data_path + subj + '/physsig/resp1.mat'
		physsig3 = data_path + subj + '/physsig/puls2.mat'
		physsig4 = data_path + subj + '/physsig/resp2.mat'


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
		realign4d.inputs.tr = 2.000

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
			if os.path.isfile(physsig1):
			    	physsig_a = physsig1
			else:
				physsig_a = None
			if os.path.isfile(physsig2):
			    	physsig_b = physsig2
			else:
				physsig_b = None

		else:
			infile = results_path+subj+'/' + data + '_1/reorient/corr_restingstatefMRI2_warp_reoriented.nii.gz'
			if os.path.isfile(physsig3):
			    	physsig_a = physsig3
			else:
				physsig_a = None
			if os.path.isfile(physsig4):
			    	physsig_b = physsig4
			else:
				physsig_b = None

		print "Running DRIFTER..."
		script="run_drifter_noSPM_kopio('%s','%s','%s')" %  (infile, physsig_a, physsig_b)
		MatlabCommand.set_default_paths(['/usr/share/spm8/','/opt2/MATLAB/NIFTI20140122/'])
		mlab = MatlabCommand(script=script, mfile=True, paths='/opt/Laskenta/Control_Room/Biomedicum/DRIFTER-toolbox/DRIFTER/', terminal_output = "stream")

		try:
			os.stat(results_path+subj+'/' + data + '_1/drifter/drifter_corrected.nii.gz')
		except:
			drifter_result = mlab.run()
			os.mkdir(results_path+subj+'/' + data + '_1/drifter/')
			os.rename(results_path+subj+'/drifter_corrected.nii.gz', results_path+subj+'/' + data + '_1/drifter/drifter_corrected.nii.gz')
			os.rename(results_path+subj+'/drifter_noise.nii.gz', results_path+subj+'/' + data + '_1/drifter/drifter_noise.nii.gz')
		
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

		workflow2.connect(inputnode2, 'drifter_result', tstat1, 'in_file')
		workflow2.connect(tstat1, 'out_file', automask, 'in_file')
		workflow2.connect(automask, 'out_file', skullstrip, 'in_file_b')
		workflow2.connect(inputnode2, 'drifter_result', skullstrip, 'in_file_a')
		workflow2.connect(skullstrip, 'out_file', tstat2, 'in_file')

		# Remove n (3) first volumes
		trim = pe.Node(interface=Trim(begin_index=3), name='trim')
		workflow2.connect(skullstrip, 'out_file', trim, 'in_file')

		# Normalize the median value of each run to 10000
		intnorm = pe.Node(interface=fsl.ImageMaths(op_string='-ing 10000', suffix='_intnorm'), name='intnorm')
		workflow2.connect(trim, 'out_file', intnorm, 'in_file')

		# Register to standard space
		mean2anat = pe.Node(fsl.FLIRT(bins=40, interp='nearestneighbour', cost='normmi', dof=12), name='mean2anat')
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
