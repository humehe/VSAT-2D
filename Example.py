import sys, os
import numpy as np
from numpy import mean,median
from progressbar import *
from termcolor import colored

from Fnc_Stk_Dir_2D import *
from Fnc_Stk_Utl_2D import *
from Fnc_Stk_Mth_2D import *
from Fnc_Stk_Tbl_2D import *
from Fnc_Stk_Spc_2D import *
from Fnc_Stk_Stk_2D import *
from Fnc_Stk_Plt_2D import *



#####################################1################################################
###################################STACKS#############################################
immoments       = True
continuum       = True
stack_light     = True #SUM MED AVG

if line == '13CO':
	restframe_frequency      =   110.20137E9           
elif line == '12CO':
	restframe_frequency      =   115.271208E9
elif line == '18CO':
	restframe_frequency      =   109.78217340E9
for element in itlpd(channel_width,sbsmn,sbsms):

	#Results
	stk_hme_dir = home + 'Stack_Results-'+ line +'-2D/'
	img_dir_res = stk_hme_dir + 'IMAGES/'
	ifc_dir_img = img_dir_res + 'COLLAPSED/' 
	stp_dir_res = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_res = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stk_dir_res = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'

	DIR_CAT_IPT  = [cats_dir]
	DIR_SPC_IPT  = [img_dir]
	DIR_RES      = [stk_hme_dir,img_dir_res,ifc_dir_img,stp_dir_res,tbl_dir_res,plt_dir_res,stk_dir_res]

	cat_ipt_tbl  =   cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) 
	Check_directories(cat_ipt_tbl,cat_parent,DIR_RES=DIR_RES)
	###
	cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

	print
	print colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'green')
	print

	Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
	id_glx        = Cat_Ipt_Tbl[1]
	fits          = Cat_Ipt_Tbl[2]
	delta_nu      = Cat_Ipt_Tbl[4]
	z             = Cat_Ipt_Tbl[8]  #[5]
	Lfir          = Cat_Ipt_Tbl[11] #[6]
	nu            = Cat_Ipt_Tbl[13] #[14]
	vel           = Cat_Ipt_Tbl[14] #[15]
	num_obj       = len(Cat_Ipt_Tbl[0])

	############################################################
	#######Correcting by number of channels used in split#######
	############################################################
	Line_Widths_Msr = cat_dir + line + '-stacker-gauss' + tbl_ext_ipt
	print
	print colored('Correction of flux by number of channels used in task split/immoments!','yellow')
	print colored('Info Line Width from table: ' + str(Line_Widths_Msr) ,'green')
	print
	Cat_Ipt_Lnw_Wdth_Tbl   =  Table_Read_Width(Line_Widths_Msr,tbl_format_ipt)
	id_glx_ref   = Cat_Ipt_Lnw_Wdth_Tbl[1]
	line_width   = Cat_Ipt_Lnw_Wdth_Tbl[23]

	indx_arr       = []
	id_arr         = []
	line_width_arr = [] 

	[indx_arr.append(np.where(galaxy_identification==id_glx_ref)[0][0]) for galaxy_identification in id_glx]
	[id_arr.append((int(id_glx_ref[indexing]))) for indexing in indx_arr]

	[line_width_arr.append(line_width[indexing]) for indexing in indx_arr]

	print
	print len(id_glx)
	print len(indx_arr)
	print len(line_width_arr)
	print 
	print 'This sample : '
	print [j for j in id_glx]
	print 'Reference   : '
	print id_arr
	print

	############################################################
	#######Correcting by number of channels used in split#######
	############################################################

	nchan         = (2*subcube_width /element[0])+1
	slice_nmb     = int(np.ceil(((nchan-1)/2)))#-1 #19 

	if immoments == True and continuum == False:
		clp_sfx = '.IM0'
		stk_sfx = str(element[0]) + 'kms.IM0'
	###CORRECTION FOR CONT####
	elif immoments == True and continuum == True:
		clp_sfx = '.cont'
		stk_sfx = str(element[0]) + 'kms.IM0.CNT'
	else:
		pass
		clp_sfx = ''
		stk_sfx = str(element[0]) + 'kms'
	###CORRECTION FOR CONT####

	fits_new = []
	if continuum == True:
		[fits_new.append(image_fits + clp_sfx+'.fits') for image_fits in (fits)]
	else:
		[fits_new.append(image_fits + '.' + str(element[0]) + 'kms'+clp_sfx+'.fits') for image_fits in (fits)]

	print
	print colored('Reading files as : '+str(fits_new[0]),'cyan')
	print

	cubetoread     = [ifc_dir_img + str(file) for file in fits_new]
	print
	print colored('Reading files as : '+str(cubetoread[0]),'cyan')
	print

	img_stack =[]
	[img_stack.append(apfts.getdata(img,memmap=False)) for img in cubetoread]
	print
	print "\n".join([str(np.asarray(spec).shape) for spec in img_stack])
	print
	weights        = np.arange(0,len(fits),1)
	weights        = np.asarray(Lfir)

	stk_ofn_prfx = cat_parent + '-' + str(element[1])
	stk_ofn_prfx = cat_parent + '-' + str(element[2]) + '-' +str(element[1])

	Stack_Res     = Cube_Stack_2D(cubetoread,stk_ofn_prfx,weights,
					sig_clp    = False      ,sufix=stk_sfx,freq_obs_f=restframe_frequency,
					stack_lite = stack_light,
					cp_bs_hdrs  = True,
					stt_var     = True,
					spc_wdt_dir = element[0],
					stt_mst_tbl = Cat_Ipt_Tbl,stt_hdr=element[2],
					lnw_wdt     = line_width_arr)
	#############################Add Headers to Stack Results##############################
	name = cat_parent + '-' + str(element[1])
	name = cat_parent + '-' + str(element[2]) + '-' +str(element[1])
	bs_func = ''
	sufix = element[0]
	spc_dir_dst = stk_dir_res

	print
	print "\n".join([file for file in Stack_Res])
	print
	###############################Add Headers to Stack Results##############################

#########################################STACKS#############################################
###########################################1################################################
###########################################2################################################
###################################FIT###############################################

circular_gaus   = True
#################################EXTRACT-SQUARE-REGION###################################
for element in itlpd(channel_width,sbsmn,sbsms,func):
	#Results
	stk_hme_dir = home + 'Stack_Results-'+ line +'-2D/'
	img_dir_res = stk_hme_dir + 'IMAGES/'
	ifc_dir_img = img_dir_res + 'COLLAPSED/' 
	stp_dir_res = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_tbl = tbl_dir_res + 'PLOTS/' 
	plt_dir_res = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stm_dir_plt = plt_dir_res + 'STAMPS/'
	ana_dir_plt = plt_dir_res + 'ANALYSIS/'
	mcm_dir_plt = plt_dir_res + 'MCMC/'    
	res_dir_plt = plt_dir_res + 'RESULTS/' 
	stk_dir_res = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'
	DIR_CAT_IPT = [cats_dir]
	DIR_SPC_IPT = [img_dir]
	DIR_RES     = [
					stk_hme_dir,stp_dir_res,tbl_dir_res,
					plt_dir_res,
					stm_dir_plt,ana_dir_plt,mcm_dir_plt,res_dir_plt,
					stk_dir_res
				]

	cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

	print colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan')
	print

	Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
	id_glx        = Cat_Ipt_Tbl[1]
	fits          = Cat_Ipt_Tbl[2]
	delta_nu      = Cat_Ipt_Tbl[4]
	z             = Cat_Ipt_Tbl[8]  #[5]
	Lfir          = Cat_Ipt_Tbl[11] #[6]
	nu            = Cat_Ipt_Tbl[13] #[14]
	vel           = Cat_Ipt_Tbl[14] #[15]
	num_obj       = len(Cat_Ipt_Tbl[0])

	############################################################
	#######Correcting by number of channels used in split#######
	############################################################
	Line_Widths_Msr = cat_dir + line + '-stacker-gauss' + tbl_ext_ipt
	Cat_Ipt_Lnw_Wdth_Tbl   =  Table_Read_Width(Line_Widths_Msr,tbl_format_ipt)
	id_glx_ref   = Cat_Ipt_Lnw_Wdth_Tbl[1]
	line_width   = Cat_Ipt_Lnw_Wdth_Tbl[23]

	indx_arr       = []
	id_arr         = []
	line_width_arr = [] 

	[indx_arr.append(np.where(galaxy_identification==id_glx_ref)[0][0]) for galaxy_identification in id_glx]
	[id_arr.append((int(id_glx_ref[indexing]))) for indexing in indx_arr]
	[line_width_arr.append(line_width[indexing]) for indexing in indx_arr]

	print
	print len(id_glx)
	print len(indx_arr)
	print len(line_width_arr)
	print 
	print 'This sample : '
	print [j for j in id_glx]
	print 'Reference   : '
	print id_arr
	print
	print colored('Correction of flux by number of channels used in task split/immoments!','yellow')
	print colored('Info Line Width from table: ' + str(Line_Widths_Msr) ,'green')
	print
	############################################################
	#######Correcting by number of channels used in split#######
	############################################################

	if immoments == True and continuum == False:
		clp_sfx = '.IM0'
		stk_sfx = str(element[0]) + 'kms.IM0'
	###CORRECTION FOR CONT####
	elif immoments == True and continuum == True:
		clp_sfx = '.cont'
		stk_sfx = str(element[0]) + 'kms.IM0.CNT'
	else:
		pass
		clp_sfx = ''
		stk_sfx = str(element[0]) + 'kms'
	###CORRECTION FOR CONT####


	#cube2bExt         = stk_dir_res+ 'CII_ALPINE-'+ str(sbsms) + '-' +str(sbsmn) + '-stk-'+ func +'-'+str(channel_width)+'kms.fits'
	cube2bExt         = stk_dir_res+ 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-'+ str(element[3])+'-'+stk_sfx+'.fits'
	print cube2bExt
	Header_Get_Add(cube2bExt,'RESTFRQ',restframe_frequency)
	scale_deg         = Header_Get_Add(cube2bExt,'CDELT2',0.000138888888889)
	scale_arcsec      = scale_deg*3600#0.00027777778
	reg_wdth_x        = 25 #arcsec 8
	reg_wdth_x        = reg_wdth_x/scale_arcsec #pixels
	reg_wdth_rad_x    = int(np.ceil(reg_wdth_x/2)) #pixels

	reg_wdth_y        = 25 #arcsec 5
	reg_wdth_y        = reg_wdth_y/scale_arcsec #pixels
	reg_wdth_rad_y    = int(np.ceil(reg_wdth_y/2)) #pixels

	X_center,Y_Center = 128,128#128,133

	########################################MASK-CIRCULAR-REGION##################################
	for mask_radi_as_ms,mask_radi_as_in,mask_radi_as_ot in zip(apertures_measu,apertures_inner,apertures_outer):
		#cube2bplot       = stk_dir_res + 'CII_ALPINE-'+ str(sbsms) + '-' +str(sbsmn) +'-stk-'+ func +'-'+str(channel_width)+'kms.fits'
		cube2bplot       = stk_dir_res + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3]+'-'+stk_sfx+'.fits'
		scale_deg        = Header_Get(cube2bplot,'CDELT2')
		scale_arcsec     = scale_deg*3600#0.00027777778
		mask_radi_px_ms  = mask_radi_as_ms / scale_arcsec #pixels
		#mask_radi_as_in  = 5 #arcsec 8 5
		mask_radi_px_in  = mask_radi_as_in / scale_arcsec #pixels
		#mask_radi_as_ot  = mask_radi_as_in * 2
		mask_radi_px_ot  = mask_radi_as_ot / scale_arcsec

		X0_F,Y0_F        = 128,128
		X0_F_0,Y0_F_0    = 128,128
		print
		print cube2bplot
		print 'Masking circular aperture.'
		print 'Center           : ',X0_F_0,Y0_F_0
		print 'Radii int[arcsec]: ',mask_radi_as_in
		print 'Radii int[pixels]: ',mask_radi_px_in
		print 'Radii out[arcsec]: ',mask_radi_as_ot
		print 'Radii out[pixels]: ',mask_radi_px_ot
		print 'Radii msr[arcsec]: ',mask_radi_as_ms
		print 'Radii msr[pixels]: ',mask_radi_px_ms
		Slices_Files = CCube_Spatial_Extract_Circular(cube2bplot,
										X0_F,Y0_F,
										mask_radi_px_in,mask_radi_as_in,
										mask_radi_px_ot,mask_radi_as_ot,
										mask_radi_px_ms,mask_radi_as_ms,
										x_ref=X0_F_0,y_ref=Y0_F_0,verbose=True,plt_slcs=True,frq_r=restframe_frequency, prefix=prefix_line,
										Splt_Hdr_Cmt_cp=element[2],
										dest_dir_stp = stp_dir_res)
		Plot_CCube_Slices(Slices_Files[0],  #CSEC_ofn_c_in,
							Slices_Files[1],#CSEC_ofn_d_in,
							Slices_Files[2],#CSEC_ofn_m_in,
							Slices_Files[3],#CSEC_ofn_c_ot,
							Slices_Files[4],#CSEC_ofn_d_ot,
							Slices_Files[5],#CSEC_ofn_m_ot,
							frq_r=restframe_frequency,
							prefix=prefix_line)
		########################################MASK-CIRCULAR-REGION###################################
		cat_tbl       = cat_dir + 'CII_Sources_ALPINE-' + str(sbsmn) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_ALPINE-' + line + '-' + str(sbsmn) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_ALPINE-' + line + '-' + str(sbsms) + '-' +str(sbsmn) + tbl_ext_ipt

		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + str(element[1]) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[1]) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

		print colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan')
		print

		Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
		fits          = Cat_Ipt_Tbl[2]
		delta_nu      = Cat_Ipt_Tbl[4]
		z             = Cat_Ipt_Tbl[8]  #[5]
		Lfir          = Cat_Ipt_Tbl[11] #[6]
		nu            = Cat_Ipt_Tbl[13] #[14]
		vel           = Cat_Ipt_Tbl[14] #[15]
		num_obj       = len(Cat_Ipt_Tbl[0])

		z_sample_avg  = np.mean(z)
		z_sample_med  = np.median(z)
		print
		print 'Redshift (avg): ',z_sample_avg
		print 'Redshift (med): ',z_sample_med
		print 'subcube_width : ',subcube_width
		print
		#########################################PLOT-FREQ PROFILE ON REGION############################
		cube2bplot1  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_ot) + 'as_msk_ot.fits'
		cube2bplot2  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_in) + 'as_msk_in.fits'
		cube2bplot3  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_ms) + 'as_msk_ms.fits'

		cube2bplot4  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_ot) + 'as_dta_ot.fits'
		cube2bplot5  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_in) + 'as_dta_in.fits'
		cube2bplot6  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + stk_sfx +'-crc-' + str(mask_radi_as_ms) + 'as_dta_ms.fits'

		Header_Get_Add(cube2bplot3,'CDELT2',0.000138888888889)
		Header_Get_Add(cube2bplot6,'CDELT2',0.000138888888889)

		###############################FIT 2D GAUSSIAN ON COLLAPSED CUBE################################
		print
		print ana_dir_plt
		print
		for function in ['sum']:#,'med','avg']:
			CCube_fit_2D_Gaussian(cube2bplot3,slc_nmb=None,
						clp_fnc      = function ,
						x_ref        = X0_F_0,y_ref=Y0_F_0,
						circular     = circular_gaus,
						SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm),SIGMAY_f2DG=fwhm2sigma(9*tms_sgm),
						displ_s_f    = True,verbose=True,						
						rst_frq      = restframe_frequency,
						frq_r        = restframe_frequency,prefix=prefix_line,sgm_fnc=None,
						sgm_wgth_tms = 'slice_1fw',src_sze_fxd = fixed_size,
						dest_dir_plt = ana_dir_plt)#1sgm-2,3,4,5sgm,slice_1fw,slice_1fw
			CCube_fit_2D_Gaussian(cube2bplot6,slc_nmb=None,
						clp_fnc      = function ,
						x_ref        = X0_F_0,y_ref=Y0_F_0,
						circular     = circular_gaus,
						SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm),SIGMAY_f2DG=fwhm2sigma(9*tms_sgm),
						displ_s_f    = True,verbose=True,						
						rst_frq      = restframe_frequency,
						frq_r        = restframe_frequency,prefix=prefix_line,sgm_fnc=None,
						sgm_wgth_tms = 'slice_1fw',src_sze_fxd = fixed_size,
						dest_dir_plt = ana_dir_plt)#1sgm-2,3,4,5sgm,slice_1fw,slice_1fw
		###############################FIT 2D GAUSSIAN ON COLLAPSED CUBE################################
#quit()
############################################2#################################################
###########################################FIT################################################
#########################################Cube Stat############################################
############################################3#################################################
func            = ['avg','med']   #sum med avg
cw              = [20]
immoments       = True
continuum       = True
for element in itlpd(channel_width,sbsmn,sbsms,func):	

	cat_ipt_tbl    = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1])
	cat_tbl        = cat_ipt_tbl + tbl_ext_ipt#

	print
	print colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan')
	print

	Cat_Ipt_Tbl = Table_Read(cat_tbl,tbl_format_ipt)
	fits        = Cat_Ipt_Tbl[2]
	delta_nu    = Cat_Ipt_Tbl[4]
	z           = Cat_Ipt_Tbl[8]
	Lfir        = Cat_Ipt_Tbl[11]
	nu          = Cat_Ipt_Tbl[13]
	vel         = Cat_Ipt_Tbl[14]
	num_obj     = len(Cat_Ipt_Tbl[0])
	#
	z_sample_avg  = np.mean(z)
	z_sample_med  = np.median(z)
	print
	print 'Redshift (avg): ',z_sample_avg
	print 'Redshift (med): ',z_sample_med
	print 'subcube_width : ',subcube_width

	#Results
	stk_hme_dir = home + 'Stack_Results-'+ line +'-2D/'
	img_dir_res = stk_hme_dir + 'IMAGES/'
	ifc_dir_img = img_dir_res + 'COLLAPSED/'  
	stp_dir_res = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_tbl = tbl_dir_res + 'PLOTS/' 
	plt_dir_res = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stk_dir_res = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'
	DIR_CAT_IPT = [cats_dir]
	DIR_SPC_IPT = [img_dir]
	DIR_RES     = [stk_hme_dir,stp_dir_res,tbl_dir_res,plt_dir_res,stk_dir_res,plt_dir_tbl]

	Check_directories(cat_ipt_tbl,cat_parent)

	apertures_measu = [15]
	apertures_inner = [10]
	apertures_outer = [20]

	if immoments == True and continuum == False:
		clp_sfx = '.IM0'
		stk_sfx = str(element[0]) + 'kms.IM0'
	###CORRECTION FOR CONT####
	elif immoments == True and continuum == True:
		clp_sfx = '.cont'
		stk_sfx = str(element[0]) + 'kms.IM0.CNT'
	else:
		pass
		clp_sfx = ''
		stk_sfx = str(element[0]) + 'kms'
	###CORRECTION FOR CONT####

	for mask_radi_as_ms,mask_radi_as_in,mask_radi_as_ot in zip(apertures_measu,apertures_inner,apertures_outer):
		cube2bplot_in  = stp_dir_res + line + '-CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+stk_sfx+'-crc-'+str(mask_radi_as_in)+'as_msk_in.fits'
		cube2bplot_ot  = stp_dir_res + line + '-CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+stk_sfx+'-crc-'+str(mask_radi_as_ot)+'as_msk_ot.fits'
		cube2bplot_ms  = stp_dir_res + line + '-CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+stk_sfx+'-crc-'+str(mask_radi_as_ms)+'as_msk_ms.fits'
		slc_nmb1 = 0#int(Header_Get(cube2bplot_in,'MAX_SNA'))
		Cube_Stat(cube2bplot_in,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency)
		Cube_Stat(cube2bplot_ot,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency)
		Cube_Stat(cube2bplot_ms,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency)
#quit()
#########################################Cube Stat############################################
############################################3#################################################

############################################5#################################################
#######################################Luminosity Plot########################################
####################################LAST VERSION WORKING #####################################
method          = 3                      #1-2DGFXSGM-AVG 2-2DGFX2FWHM-AVG(SPL) 3-2DGF (IM0)
error           = 1                      #N sigma confidence intervals
plt_scl         = None                   #Z, Y, both None
log_lm          = 'error'                #both,value,error, None
func1           = 'avg'
func2           = 'med'
cw              = 20
sbsmn           = [0]                 #A, L, H
sbsms           = 'CNT_B'
mask_radi_as_ms = 15
continuum       = False

ERR_MC_ERR_CMP  = False
iterations_mc   = 10000
plot_dist_hist  = True
ERR_MC_ERR_PLT  = True
var_ctr_val_avg = True #Use average as central value for a given property, False uses median

lit_res         = True   #Plot literature Results
ind_res         = True   #Plot Individual Measurements
immoments       = True

MCMC_generator(iterations_mc,line1,line2,method,error,'CNT_B',sbsmn,spc_wdt_dir=cw,mask_radi_as_ms=mask_radi_as_ms)
quit()
