"""
Nipype-pipeline: fMRI preprocessing

Functional data preprocessing steps:

1. 4D realignment: motion and slice timing correction
2. Reorientation to FSL standard orientation (no deoblique)
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
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.interfaces.matlab import MatlabInputSpec, MatlabCommand    # matlab
from nipype.interfaces.nipy.preprocess import remove3vol
import os
import shutil
import copy
import numpy as np
import glob
from subprocess import call

#########################################################################################

data_path = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Functional_scans/'
results_path = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Results_nodeoblique/Functional_nosearch/data2/'

subj_dirs = os.listdir(data_path)

for data in (['func1','func2']):

	os.chdir(results_path)

	# Define input files: 2xfMRI + 1xMPRAGE
	func1 = data_path +  '/sessio2_a/epi.nii.gz'
	func2 = data_path + '/sessio2_b/epi.nii.gz'
	anat = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Results_nodeoblique/00002_session2/anat/reorient/anat_LR_brain_reoriented.nii.gz'
	field_file = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Results_nodeoblique/00002_session2/anat/nonlinear_reg/anat_LR_brain_reoriented_fieldwarp.nii.gz'

	# Physiological signal files
	physsig1 = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Biopac/2.sessio/pulse_fMRI1.mat'
	physsig2 = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Biopac/2.sessio/resp_fMRI1.mat'
	physsig3 = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Biopac/2.sessio/pulse_fMRI2.mat'
	physsig4 = '/opt/MR_data/Silja_Raty/Revis_0017_rs/Biopac/2.sessio/resp_fMRI2.mat'


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
	#deoblique = pe.Node(interface=afni.Warp(deoblique=True, outputtype='NIFTI_GZ'), name='deoblique')
	#workflow1.connect(realign4d, 'out_file', deoblique, 'in_file')
	reorient = pe.Node(interface=fsl.Reorient2Std(output_type='NIFTI_GZ'), name='reorient')
	workflow1.connect(realign4d, 'out_file', reorient, 'in_file')
	workflow1.connect(reorient, 'out_file', outputnode, 'out_file')
	
	# Run workflow1
	workflow1.write_graph()
	workflow1.run()
	
	# Initialize DRIFTER
	if data=='func1':
		#check is physiological signals exist, if not, then None
		infile = results_path+'/' + data + '_1/reorient/corr_epi_reoriented.nii.gz'
		physsig_a = physsig1
		if not os.path.isfile(physsig1):
		    	print "MISSING PULS DATA!"
		physsig_b = physsig2
		if not os.path.isfile(physsig2):
			print "MISSING RESP DATA!"

	else:
		infile = results_path+'/' + data + '_1/reorient/corr_epi_reoriented.nii.gz'
		physsig_a = physsig3
		if not os.path.isfile(physsig3):
		    	print "MISSING PULS DATA!"
		physsig_b = physsig4
		if not os.path.isfile(physsig4):
			print "MISSING RESP DATA!"

	print "Running DRIFTER..."
	script="run_drifter_noSPM_kopio('%s','%s','%s')" %  (infile, physsig_a, physsig_b)
	MatlabCommand.set_default_paths(['/opt/Laskenta/Control_Room/Biomedicum/DRIFTER-toolbox/DRIFTER/', '/opt/MATLAB/R2015a/spm8/','/opt/MATLAB/NIfTI_20140122/'])
	mlab = MatlabCommand(script=script, mfile=True, paths=['/opt/Laskenta/Control_Room/Biomedicum/DRIFTER-toolbox/DRIFTER/', '/opt/MATLAB/R2015a/spm8/','/opt/MATLAB/NIfTI_20140122/'], terminal_output = "stream")

	try:
		os.stat(results_path + data + '_1/drifter/drifter_corrected.nii.gz')
	except:
		drifter_result = mlab.run()
		os.mkdir(results_path + data + '_1/drifter/')
		os.rename(results_path+'/drifter_corrected.nii.gz', results_path+'/' + data + '_1/drifter/drifter_corrected.nii.gz')
	
	# Initialize workflow2
	workflow2 = pe.Workflow(name=data + '_2')
	workflow2.base_dir = '.'

	inputnode2 = pe.Node(interface=util.IdentityInterface(fields=['drifter_result']), name='inputspec2')
	outputnode2 = pe.Node(interface=util.IdentityInterface(fields=['result_func']), name='outputspec2')

	inputnode2.inputs.drifter_result=results_path+'/' + data + '_1/drifter/drifter_corrected.nii.gz'


	# Call fslcpgeom source dest, source is reorient output nii.gz file and dest is drifter folder nii.gz file
	reoriented_file=results_path+'/' + data + '_1/reorient/corr_epi_reoriented.nii.gz' 
	drifted_file=results_path+'/' + data + '_1/drifter/drifter_corrected.nii.gz'
	call(["fslcpgeom", reoriented_file, drifted_file])

	# AFNI skullstrip and mean image skullstrip
	mean_epi = pe.Node(interface=afni.TStat(args='-mean',outputtype="NIFTI_GZ"), name='mean_epi')
	skullstrip3D = pe.Node(interface=afni.Automask(dilate=1,outputtype="NIFTI_GZ"), name='skullstrip3D')
	skullstrip4D = pe.Node(interface=afni.Calc(expr = 'a*b',outputtype="NIFTI_GZ"), name='skullstrip4D’)
	mean_epi_brain = pe.Node(interface=afni.TStat(args='-mean',outputtype="NIFTI_GZ"), name='mean_epi_brain')

	workflow2.connect(inputnode2, 'drifter_result', mean_epi,'in_file')
	workflow2.connect(mean_epi, 'out_file', skullstrip3D, 'in_file')
	workflow2.connect(skullstrip3D, 'out_file', skullstrip, 'in_file_b')
	workflow2.connect(inputnode2, 'drifter_result', skullstrip, 'in_file_a')
	workflow2.connect(skullstrip, 'out_file', mean_epi_brain, 'in_file')

	# Remove n (3) first volumes
	remove3vol = pe.Node(interface=remove3vol(begin_index=3), name='remove3vol')
	workflow2.connect(skullstrip, 'out_file', remove3vol, 'in_file')

	# Spatial smoothing, kernel sigma 2.00 mm (5 mm is too much)
	#smooth = pe.Node(interface=fsl.maths.SpatialFilter(operation='mean', terminal_output='stream', kernel_shape='gauss', kernel_size=1.5, 			nan2zeros=True), name='smooth')

	#workflow2.connect(remove3vol,'out_file', smooth, 'in_file')

	# Normalize the median value of each run to 10000
	intnorm = pe.Node(interface=fsl.ImageMaths(op_string='-ing 10000', suffix='_intnorm'), name='intnorm')
	workflow2.connect(remove3vol, 'out_file', intnorm, 'in_file')

	# Register to standard space
	#mean2anat = pe.Node(fsl.FLIRT(bins=40, cost='normmi', dof=12, interp='nearestneighbour'), name='mean2anat')
	mean2anat = pe.Node(fsl.FLIRT(no_search=True), name='mean2anat')
	mean2anat.inputs.reference = anat
	workflow2.connect(mean_epi_brain, 'out_file', mean2anat, 'in_file')

	# Transform mean functional image
	warpmean = pe.Node(interface=fsl.ApplyWarp(), name='warpmean')
	warpmean.inputs.ref_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
	warpmean.inputs.field_file = field_file
	workflow2.connect(mean_epi_brain, 'out_file', warpmean, 'in_file')
	workflow2.connect(mean2anat, 'out_matrix_file', warpmean, 'premat')

	# Transform all images
	warpall = pe.Node(interface=fsl.ApplyWarp(), name='warpall')
	warpall.inputs.ref_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
	warpall.inputs.field_file = field_file
	workflow2.connect(intnorm, 'out_file', warpall, 'in_file')
	workflow2.connect(mean2anat, 'out_matrix_file', warpall, 'premat')
	workflow2.connect(warpall, 'out_file', outputnode2, 'result_func')

	# Run workflow2
	workflow2.write_graph()
	workflow2.run()

	print "FUNCTIONAL PREPROCESSING DONE! Results in ", results_path
