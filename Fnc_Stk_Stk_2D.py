import bottleneck as bn
from astropy import stats as apsts

import scipy.integrate as integrate

from Fnc_Stk_Dir_2D import *
from Fnc_Stk_Spc_2D import *
from Fnc_Stk_Tbl_2D import *

####Fnc_Stk_Stk_2D####
def Cube_Stack_2D(CCubes2bStacked,name,wght_img_2bstack,sig_clp,*args, **kwargs):
	wrt_fits         = kwargs.get('wrt_fits'       ,True)
	pst_msk          = kwargs.get('pst_msk'        ,False)
	pst_smt          = kwargs.get('pst_smt'        ,False)
	pst_cnt          = kwargs.get('pst_cnt'        ,False)
	stack_ext        = kwargs.get('stack_ext'      ,None)
	new_CRVAL1_head  = kwargs.get('new_CRVAL1_head',None)
	new_CDELT1_head  = kwargs.get('new_CDELT1_head',None)
	smt_spc_pst      = kwargs.get('smt_spc_pst'    ,False)
	smooth_shape     = kwargs.get('smooth_shape'   ,'gaussian')
	wght_type        = kwargs.get('wght_type'      ,None)
	wcs              = kwargs.get('wcs'            ,None)
	sufix            = kwargs.get('sufix'          ,'')
	freq_obs_f       = kwargs.get('freq_obs_f'     ,99999)

	stack_lite       = kwargs.get('stack_lite'     ,True)

	spc_wdt_dir      = kwargs.get('spc_wdt_dir'    ,500)

	cp_bs_hdrs       = kwargs.get('cp_bs_hdrs'     ,False)

	stt_var			 = kwargs.get('stt_var',False)
	stt_mst_tbl		 = kwargs.get('stt_mst_tbl',None)
	stt_hdr			 = kwargs.get('stt_hdr',None)

	lnw_wdt			 = kwargs.get('lnw_wdt',None)

	img_2bstack    = [apgtdt(img,memmap=False) for img in CCubes2bStacked]
	wcs            = kwargs.get('wcs'            ,apwcs(CCubes2bStacked[0]))
	try:
		wcs       = wcs.dropaxis(3) 
	except IndexError:
		pass

	print
	print 'Number of galaxies to be stacked (histogram): ',len(img_2bstack)

	if sig_clp == True:
		img_flt       = astropy.stats.sigma_clip(img_2bstack,sigma=sigma_cut,axis=0,iters=None,cenfunc=sigma_cen_fct, copy=True)

		print
		print colored('Sigma-clipping for stacking!','yellow')
		print colored('Sigma Cut                    : ' + str(sigma_cut),'yellow')
		print colored('Central function             : ' + str(sigma_cen_fct), 'yellow')
		print colored('Central Value for clipping   : ' + str(sigma_cen_fct),'yellow')

		img_flt.set_fill_value(sigma_msk_fill_val)
		img_flt_filled = img_flt.filled()
		img_stat       = img_flt_filled
	elif sig_clp == False:
		img_stat   = img_2bstack

	print 
	print np.asarray(img_stat).shape
	img_stat = np.squeeze(img_stat) 
	print
	print np.asarray(img_stat).shape

	wght_img_copy = wght_img_2bstack
	wght_img_stat = wght_img_2bstack 
	wght_img_stat = np.asarray(wght_img_stat)
	img_staw      = []

	img_stat_smw_f = []
	[img_staw.append(np.asarray(img_stat)[j]*np.asarray(wght_img_stat)[j]) for j in range(len(wght_img_stat))]
	img_staw      = np.asarray(img_staw)
	[img_stat_smw_f.append(np.divide(np.asarray(img_staw)[j],np.asarray(img_stat)[j])) for j in range(len(img_stat))]

	print
	print colored('Original shape                                               : '+str(np.asarray(img_stat).shape),'cyan')
	img_stat = np.squeeze(img_stat) 
	N,l,m  = np.asarray(img_stat).shape 
	print colored('Squeezed useless extra dimensions                            : '+str(np.asarray(img_stat).shape),'cyan')
	print colored('Dimension Numer of Cubes, X size, Y size : '+str(N)+', '+str(l)+', '+str(m),'cyan')
	print

	img_res_sum = bn.nansum(np.array(img_stat)             , axis=0)
	img_res_avg = bn.nanmean(np.array(img_stat)            , axis=0)
	img_res_med = bn.nanmedian(np.array(img_stat)          , axis=0)

	print
	print colored('Sum, Mean, Median : Stacked data cubes OK','yellow')
	print 

	if stack_lite == False:

		img_stat_hst_y = []
		img_stat_hsw_y = []

		pb = ProgressBar(l)
		for y_dim in range(l):
			pb.update()
			Y_ROW = np.asarray(img_stat)[:,y_dim,:]
			Transpose  = np.asarray(Y_ROW).T
			Transposw  = np.asarray(Y_ROW).T
			img_stat_hst_x = []
			img_stat_hsw_x = []

			for x_dim in range(len(Transpose)):
				if np.isnan(sigma_msk_fill_val) == True:
					non_msk_num = int(np.count_nonzero(~np.isnan(Transpose[x_dim])))
					msk_num     = int(np.count_nonzero(np.isnan(Transpose[x_dim])))
					img_stat_hst_x.append(float(non_msk_num))

					non_msk_num_wghts = int(np.count_nonzero(~np.isnan(Transposw[x_dim])))
					msk_num_wghts     = int(np.count_nonzero(np.isnan(Transposw[x_dim])))
					img_stat_hsw_x.append(float(non_msk_num_wghts))

				elif np.isnan(sigma_msk_fill_val) == False:
					pass
					non_msk_num = int(np.count_nonzero(Transpose[x_dim]!=sigma_msk_fill_val))
					img_stat_hst_x.append(float(non_msk_num))
					non_msk_num_wghts = int(np.count_nonzero(Transposw[x_dim]!=sigma_msk_fill_val))
					img_stat_hsw_x.append(float(non_msk_num_wghts))
				else:
					pass
			
			img_stat_hst_x = np.reshape(img_stat_hst_x,(m))
			img_stat_hsw_x = np.reshape(img_stat_hsw_x,(m))

			img_stat_hst_y.append(img_stat_hst_x)
			img_stat_hsw_y.append(img_stat_hsw_x)
			#ENDS HISTO

		img_sts_hst = np.asarray(img_stat_hst_y)
		img_res_std = bn.nanstd(np.array(img_stat), axis=0)

		print
		print colored('Histogram, Std: Stacked data cubes OK','yellow')
		print 

		img_res_suw_pre = np.asarray(bn.nansum(np.array(img_staw)                , axis=0))
		img_sts_wsu_pre = np.asarray(bn.nansum(np.array(img_stat_smw_f)          , axis=0))

		img_sts_wsu_pre = np.squeeze(img_sts_wsu_pre)
		img_res_suw_pre = np.squeeze(img_res_suw_pre)

		print
		print colored('Weights Sum Weighted Sum pre computations: OK','yellow')
		print
		print img_res_suw_pre.shape

		img_sts_hsw = data=np.asarray(img_stat_hsw_y)                           #histogram of weighted cubes
		img_sts_wsu = data=img_sts_wsu_pre                                      #sum of weights
		img_res_suw = data=img_res_suw_pre                                      #weighted sum
		img_res_avw = data=img_res_suw_pre.astype(float)/img_sts_wsu_pre.astype(float) 
		
		print
		print colored('SW Histogram, Sum of weights, Weighted Sum: Stacked data cubes OK','yellow')
		print

		img_res_1sl = np.nanpercentile(np.array(img_stat), 15.9, axis=0)
		img_res_1sh = np.nanpercentile(np.array(img_stat), 84.1, axis=0)
		img_res_2sl = np.nanpercentile(np.array(img_stat), 2.30, axis=0)
		img_res_2sh = np.nanpercentile(np.array(img_stat), 97.7, axis=0)
		img_res_3sl = np.nanpercentile(np.array(img_stat), 0.20, axis=0)
		img_res_3sh = np.nanpercentile(np.array(img_stat), 99.8, axis=0)
		img_res_p25 = np.nanpercentile(np.array(img_stat), 25.0, axis=0)
		img_res_p75 = np.nanpercentile(np.array(img_stat), 75.0, axis=0)


		print 'Stacked images through : sum, mean, median, and percentiles: '
		print '17., 83.0, (1 sigma)'
		print '2.5, 97.5, (2 sigma)'
		print '0.5, 99.5, (3 sigma)'
		print '25., 75.0, (interquantile)'
		print		
		print colored('Percentiles: Stacked data cubes OK','yellow')
		print 
	elif stack_lite == True:
		pass


	bs_func = kwargs.get('bs_func','')

	if wrt_fits==True:
		if  '-BS-' in name:
			print (colored(name,'yellow'))
			spc_dir_dst = str_bst_stk + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass
		elif  '-BS_MST' in name:
			print (colored(name,'yellow'))
			spc_dir_dst = stt_bst_stk + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass
		else:
			spc_dir_dst = stk_dir_res + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass

		spec_file_sum_ofn = spc_dir_dst + str(name) + bs_func + '-stk-sum-' + str(sufix) + '.fits'
		spec_file_avg_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avg-' + str(sufix) + '.fits'
		spec_file_med_ofn = spc_dir_dst + str(name) + bs_func + '-stk-med-' + str(sufix) + '.fits'
		spec_file_hst_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hst-' + str(sufix) + '.fits'
		spec_file_std_ofn = spc_dir_dst + str(name) + bs_func + '-stk-std-' + str(sufix) + '.fits'
		spec_file_p25_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p25-' + str(sufix) + '.fits'
		spec_file_p75_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p75-' + str(sufix) + '.fits'
		spec_file_1sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sl-' + str(sufix) + '.fits'
		spec_file_1sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sh-' + str(sufix) + '.fits'
		spec_file_2sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sl-' + str(sufix) + '.fits'
		spec_file_2sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sh-' + str(sufix) + '.fits'
		spec_file_3sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sl-' + str(sufix) + '.fits'
		spec_file_3sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sh-' + str(sufix) + '.fits'
	

		spec_file_hsw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hsw-' + str(sufix) + '.fits'
		spec_file_wsu_ofn = spc_dir_dst + str(name) + bs_func + '-stk-wsu-' + str(sufix) + '.fits'
		spec_file_suw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-suw-' + str(sufix) + '.fits'
		spec_file_avw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avw-' + str(sufix) + '.fits'

		spec_file_sum     = Wrt_FITS_File(img_res_sum,spec_file_sum_ofn)
		spec_file_avg     = Wrt_FITS_File(img_res_avg,spec_file_avg_ofn)
		spec_file_med     = Wrt_FITS_File(img_res_med,spec_file_med_ofn)

		if stack_lite == False:
			spec_file_hst     = Wrt_FITS_File(img_sts_hst,spec_file_hst_ofn)
			spec_file_std     = Wrt_FITS_File(img_res_std,spec_file_std_ofn)

			spec_file_p25     = Wrt_FITS_File(img_res_p25,spec_file_p25_ofn)
			spec_file_p75     = Wrt_FITS_File(img_res_p75,spec_file_p75_ofn)
			spec_file_1sl     = Wrt_FITS_File(img_res_1sl,spec_file_1sl_ofn)
			spec_file_1sh     = Wrt_FITS_File(img_res_1sh,spec_file_1sh_ofn)
			spec_file_2sl     = Wrt_FITS_File(img_res_2sl,spec_file_2sl_ofn)
			spec_file_2sh     = Wrt_FITS_File(img_res_2sh,spec_file_2sh_ofn)
			spec_file_3sl     = Wrt_FITS_File(img_res_3sl,spec_file_3sl_ofn)
			spec_file_3sh     = Wrt_FITS_File(img_res_3sh,spec_file_3sh_ofn)

			spec_file_hsw     = Wrt_FITS_File(img_sts_hsw,spec_file_hsw_ofn)
			spec_file_wsu     = Wrt_FITS_File(img_sts_wsu,spec_file_wsu_ofn)
			spec_file_suw     = Wrt_FITS_File(img_res_suw,spec_file_suw_ofn)
			spec_file_avw     = Wrt_FITS_File(img_res_avw,spec_file_avw_ofn)

			OPT_STCK_FLS = [spec_file_sum_ofn,spec_file_avg_ofn,spec_file_med_ofn,spec_file_hst_ofn,
			spec_file_std_ofn,
			spec_file_p25_ofn,spec_file_p75_ofn,
			spec_file_1sl_ofn,spec_file_1sh_ofn,
			spec_file_2sl_ofn,spec_file_2sh_ofn,
			spec_file_3sl_ofn,spec_file_3sh_ofn,
			spec_file_hsw_ofn,spec_file_wsu_ofn,spec_file_suw_ofn,spec_file_avw_ofn]
		elif stack_lite == True:
			OPT_STCK_FLS = [spec_file_sum_ofn,spec_file_avg_ofn,spec_file_med_ofn]
		[Header_Updt(spec_sts_res,'STK_NUM' ,len(img_2bstack), header_comment = 'Number of galaxies used for Stack') for spec_sts_res in OPT_STCK_FLS]
	else:
		pass


	print 'Images Stacked files names: '
	print
	print colored(spec_file_sum_ofn,'cyan')
	print colored(spec_file_avg_ofn,'cyan')
	print colored(spec_file_med_ofn,'cyan')
	if stack_lite == False:
		print colored(spec_file_hst_ofn,'cyan')
		print colored(spec_file_std_ofn,'cyan')
		print colored(spec_file_p25_ofn,'cyan')
		print colored(spec_file_p75_ofn,'cyan')
		print colored(spec_file_1sl_ofn,'cyan')
		print colored(spec_file_1sh_ofn,'cyan')
		print colored(spec_file_2sl_ofn,'cyan')
		print colored(spec_file_2sh_ofn,'cyan')
		print colored(spec_file_3sl_ofn,'cyan')
		print colored(spec_file_3sh_ofn,'cyan')

		print colored(spec_file_hsw_ofn,'yellow')
		print colored(spec_file_avw_ofn,'yellow')
		print colored(spec_file_suw_ofn,'yellow')

	elif stack_lite == True:
		pass

	if stack_lite == True:
		FNL_SPEC_RES = [spec_file_med,spec_file_avg,spec_file_sum]
	elif stack_lite == False:
		FNL_SPEC_RES = [
					spec_file_med,spec_file_avg,spec_file_sum,spec_file_std,
					spec_file_hst,
					spec_file_1sl,spec_file_1sh,
					spec_file_2sl,spec_file_2sh,
					spec_file_3sl,spec_file_3sh,
					spec_file_p25,spec_file_p75,
					spec_file_hsw,spec_file_wsu,spec_file_suw,spec_file_avw]

	lnewdth_sample_avg  = np.mean(lnw_wdt)
	lnewdth_sample_med  = np.median(lnw_wdt)
	lnewdth_sample_1sl  = np.nanpercentile(lnw_wdt, 15.9)
	lnewdth_sample_1sh  = np.nanpercentile(lnw_wdt, 84.1)
	lnewdth_sample_2sl  = np.nanpercentile(lnw_wdt, 2.30)
	lnewdth_sample_2sh  = np.nanpercentile(lnw_wdt, 97.7)
	lnewdth_sample_3sl  = np.nanpercentile(lnw_wdt, 0.20)
	lnewdth_sample_3sh  = np.nanpercentile(lnw_wdt, 99.8)
	lnewdth_sample_p25  = np.nanpercentile(lnw_wdt, 25.0)
	lnewdth_sample_p75  = np.nanpercentile(lnw_wdt, 75.0)

	[Header_Get_Add(Stacked_Cube,'STW_AVG',lnewdth_sample_avg,header_comment='Line Width [GHz] Average')               for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_MED',lnewdth_sample_med,header_comment='Line Width [GHz] Median')                for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_1SL',lnewdth_sample_1sl,header_comment='Line Width [GHz] 1 sgm lw lmt 15.9 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_1SH',lnewdth_sample_1sh,header_comment='Line Width [GHz] 1 sgm hg lmt 84.1 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_2SL',lnewdth_sample_2sl,header_comment='Line Width [GHz] 2 sgm lw lmt 2.30 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_2SH',lnewdth_sample_2sh,header_comment='Line Width [GHz] 2 sgm hg lmt 97.7 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_3SL',lnewdth_sample_3sl,header_comment='Line Width [GHz] 3 sgm lw lmt 0.20 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_3SH',lnewdth_sample_3sh,header_comment='Line Width [GHz] 3 sgm hg lmt 99.8 pct') for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_P25',lnewdth_sample_p25,header_comment='Line Width [GHz] 25 pct')                for Stacked_Cube in OPT_STCK_FLS]
	[Header_Get_Add(Stacked_Cube,'STW_P75',lnewdth_sample_p75,header_comment='Line Width [GHz] 75 pct')                for Stacked_Cube in OPT_STCK_FLS]

	if cp_bs_hdrs == True:
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BSCALE')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BZERO')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BMAJ')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BMIN')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BPA')      for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BTYPE')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'EQUINOX')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'RADESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'BUNIT')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'RADESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'LONPOLE')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'LATPOLE')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC1_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC2_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC3_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC1_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC2_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC3_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC1_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC2_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PC3_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CTYPE1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRVAL1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CDELT1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRPIX1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CUNIT1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CTYPE2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRVAL2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CDELT2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRPIX2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CUNIT2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CTYPE3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRVAL3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CDELT3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CRPIX3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'CUNIT3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PV2_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'PV2_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'RESTFRQ')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'SPECSYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'ALTRVAL')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'ALTRPIX')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'VELREF')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'TELESCOP') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSERVER') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'DATE-OBS') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'TIMESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSRA')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSDEC')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSGEO-X') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSGEO-Y') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'OBSGEO-Z') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'DATE')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,CCubes2bStacked[0],'ORIGIN')   for stk_res_flr in OPT_STCK_FLS]
	else:
		pass
	if stt_var == True:
		print
		print (colored('Adding stat to fits headers!','yellow'))
		print
		tbl_sts = Table_Ipt_Cat_Stats(stt_mst_tbl,stt_hdr)
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][0] ,tbl_sts[1][0] ,header_comment='Redshift Average')                         for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][1] ,tbl_sts[1][1] ,header_comment='Redshift Median')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][2] ,tbl_sts[1][2] ,header_comment='Redshift 1 sgm lw lmt 15.9 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][3] ,tbl_sts[1][3] ,header_comment='Redshift 1 sgm hg lmt 84.1 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][4] ,tbl_sts[1][4] ,header_comment='Redshift 2 sgm lw lmt 2.30 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][5] ,tbl_sts[1][5] ,header_comment='Redshift 2 sgm hg lmt 97.7 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][6] ,tbl_sts[1][6] ,header_comment='Redshift 3 sgm lw lmt 0.20 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][7] ,tbl_sts[1][7] ,header_comment='Redshift 3 sgm hg lmt 99.8 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][8] ,tbl_sts[1][8] ,header_comment='Redshift 25 pct')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][9] ,tbl_sts[1][9] ,header_comment='Redshift 75 pct')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][10],tbl_sts[1][10],header_comment=str(tbl_sts[2]) + ' Average')               for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][11],tbl_sts[1][11],header_comment=str(tbl_sts[2]) + ' Median')                for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][12],tbl_sts[1][12],header_comment=str(tbl_sts[2]) + ' 1 sgm lw lmt 15.9 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][13],tbl_sts[1][13],header_comment=str(tbl_sts[2]) + ' 1 sgm hg lmt 84.1 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][14],tbl_sts[1][14],header_comment=str(tbl_sts[2]) + ' 2 sgm lw lmt 2.30 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][15],tbl_sts[1][15],header_comment=str(tbl_sts[2]) + ' 2 sgm hg lmt 97.7 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][16],tbl_sts[1][16],header_comment=str(tbl_sts[2]) + ' 3 sgm lw lmt 0.20 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][17],tbl_sts[1][17],header_comment=str(tbl_sts[2]) + ' 3 sgm hg lmt 99.8 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][18],tbl_sts[1][18],header_comment=str(tbl_sts[2]) + ' 25 pct')                for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][19],tbl_sts[1][19],header_comment=str(tbl_sts[2]) + ' 75 pct')                for Stacked_Cube in OPT_STCK_FLS]		
	else:
		pass
	return OPT_STCK_FLS
####Fnc_Stk_Stk_2D####