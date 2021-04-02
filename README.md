# VSAT-2D
Valparaíso Stacking Analysis Tool 2D is part of the Valparaíso Stacking Analysis Tool (VSAT), and provide a series of tools for selecting, stacking, and analyzing _moment-0_ intensity maps from interferometric datasets . It is intended for stacking samples of datacubes belonging to large extragalactic catalogs by selecting subsamples of galaxies defined by their available properties (_e.g. redshift, stellar mass, star formation rate_) being possible to generate diverse (_e.g. median, average, weighted average, histogram_) composite spectra. However, it is possible to also use VSAT on smaller datasets containing any type of astronomical object and any type of 2D image.

![Alt text](./Figures-IM/Scheme2.jpg?raw=true "3D moment-0 map Scheme.")

## Content

1. Fnc_Stk_Dir.py:
   - Location of the input catalogue and spectral data. 
   - Parameters for selecting subsamples of galaxies according to their physical properties.
   - MCMC parameters.
   - Location of the resulting products of the stacking analyses _e.g. stamps, tables, plots,and stacked spectra_.

2. Fnc_Stk_Mth.py:
   - Math functions (e.g. cosmological constants, gaussian profiles for line fitting) needed throughout the stacking analysis.

3. Fnc_Stk_Spc.py 
   - Tools for modyfing datacubes including _e.g. masking, adding missing frequencies, eextract regions etc_

4. Fnc_Stk_Stt.py 
   - Statistical funtions for datacubes.

5. Fnc_Stk_Plt.py
   - Plot templates used throughout the stacking analysis. 

6. Fnc_Stk_Stk.py
   - Core of the 3D stacking tool.

7. Fnc_Stk_Fts.py
   - Funtions to access and modify (add, modify, delete) fits headers

8. Fnc_Stk_Tbl.py
   - Functions to read, write and modify different tables.
 
 9. Fnc_Stk_Utl.py
   - Auxiliary functions for the stacking analysis.

## Parameters
VSAT-2D generates composite 2D images coming from _moment-0_ intensity maps generated with CASA from ALMA interferometric observation. However it is possible to combine any other set of 2D images. After the composite datacubes are generated, it is possible to measure the flux of a source through a gaussian model. 

###### "Stacking"
There are two different options to use the stacking procedure a _lite_ version (```stack_lite=True```) which will generate _sum, median and average_ compositte datacubes and a _full_ version (```stack_lite=False```) which additionally will create _histograms, weighted average, percentiles_ composite datacubes. By default the lite version is defined. Through ```sigma_clipping=True```it is possible to exlude outliers that exceed n-times (```sigma_cut```) the mean/median ``` sigma_cen_fct ``` of the stacked pixels. 

###### "Stamps"
To measure the flux it is possible to create smaller datacubes ("_stamps_") around any partticular _ra, dec_ position withiin a circular region. ```apertures_measu``` defines the flux measurement regioin, while ```apertures_inner``` and ```apertures_outer```define an outter ring useful for noise estimations.  

###### "Fitting"
The flux estimation is computed analytically through a 3D-gauussian model. First the spectrral/velociity location of the maximum flux emission is determined through a 1D gaussian model, althoough it is possible to fix the channel at which the peak is located. Then the flux contained in a region previously defined by ```apertures_measu```is computed through a 2D gaussian profile to obtain the size ($\sigma_{x,y}$) and the amplitude (_A_).

###### "MCMC"
To compute the Confident Intervals (CIs) of the flux measurments it is possible to run Monte Carlo simulations defined by the flux measurements previously computed and by the statistical properties of the used sample. ```iterations_mc``` define the nuumer of repetititions, ```plot_dist_hist=True``` will create hiistograms of the simulations if the lines defined by ```line1```and ```line2```.

## Example

The Exaample.py script contains an example too stack a sample of 27 galaxies belonging to the Valpara\'iso ALMA/APEX Line Emission Survey(VALES). The sample of spectra can be downloaded from the [zenodo repository](). Then by simple running ```python Example.py``` will complete all the following steps below. The following  snippets are extracts contained in the Example.py file and will guide you through the file. 

###### "Stacking"
The following snippet will stack the galaxies.

```python
Stack_Res     = Cube_Stack_2D(cubetoread,stk_ofn_prfx,weights,
				sig_clp     = False      ,
				sufix       = stk_sfx.   ,
				freq_obs_f  = restframe_frequency,
				stack_lite  = stack_light,
				cp_bs_hdrs  = True,
				stt_var     = True,
				spc_wdt_dir = element[0],
				stt_mst_tbl = Cat_Ipt_Tbl,stt_hdr=element[2],
				lnw_wdt     = line_width_arr)
```
This will generate the following fits files in the results directory (```~/Example/Stack_Results-12CO-2D/STACKS/20/```):

```
 - CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-med-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-sum-20kms.IM0.CNT.fits
```

If ```stack_lite = False``` additional composite spectra will be gnerated:

```
 - CII_HATLAS-CNT_B-0-stk-hst-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-std-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-hsw-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-suw-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-wsu-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-avw-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-1sl-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-1sh-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-2sl-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-2sh-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-3sl-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-3sh-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-p25-20kms.IM0.CNT.fits
 - CII_HATLAS-CNT_B-0-stk-p75-20kms.IM0.CNT.fits
```


###### "Stamps"
To measure the source's flux, stamps can be created. The following snippet will create a 15'', 10'' and 20'' stamps considering ```CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT.fits``` as inputt image. In this example these regions use the image center (```X0_F, Y0_F```) as reference but this can be defined with the ```X_C,Y_C```parameters.

```
python
Slices_Files = Cube_Spatial_Extract_Circular(cube2bplot,
						X0_F,Y0_F,
						mask_radi_px_in,mask_radi_as_in,
						mask_radi_px_ot,mask_radi_as_ot,
						mask_radi_px_ms,mask_radi_as_ms,		
						x_ref=X0_F_0,y_ref=Y0_F_0,
						verbose=True,
						frq_r=restframe_frequency, prefix=prefix_line,
						Splt_Hdr_Cmt_cp=element[2],
						dest_dir_stp = stp_dir_res)
```
The stamps will be located in the ```~/Example/Stack_Results-12CO-2D/STAMPS/20/``` directory, and correspond to three different regions: measurement (ms), inner (in) and outter (ot) witth their corresponding fits files: circunscribed (crc), data (dta), and masked (msk) .

```
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-10as_crc_in.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-10as_dta_in.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-10as_msk_in.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-20as_crc_ot.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-20as_dta_ot.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-20as_msk_ot.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-15as_crc_ms.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-15as_dta_ms.fits
 - 12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-15as_msk_ms.fits
```

Which can be plotted to visualize the stamps.

```
python
Plot_CCube_Slices(Slices_Files[0],#CSEC_ofn_c_in,
		Slices_Files[1],
		Slices_Files[2],
		Slices_Files[3],
		Slices_Files[4],
		Slices_Files[5],
		frq_r=restframe_frequency,
		prefix=prefix_line)
```

![Alt text](./Figures-IM/12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-10as_crc_in-slices.jpg?raw=true "3D moment-0 map Scheme.")

###### "2D-Gaussian Fit"


Finally a 2D gaussian fit can be performed. 

```
python
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
```

This will generate a figure with three panels including the image, the moodel and the residuals. Additioinally model and resiidual fits files will be created in the ```~/Example/Stack_Results-12CO-2D/STAMPS/``` directory.

```
- 12CO-CII_HATLAS-CNT_B-0-stk-med-20kms.IM0.CNT-crc-15as_dta_ms-2DCGF-sum-MDL.fits
- 12CO-CII_HATLAS-CNT_B-0-stk-med-20kms.IM0.CNT-crc-15as_dta_ms-2DCGF-sum-RSD.fits
```

![Alt text](./Figures-IM/12CO-CII_HATLAS-CNT_B-0-stk-avg-20kms.IM0.CNT-crc-15as_dta_ms-2DCGF-sum-RSD.jpg?raw=true "Stacked spectra COSMOS field.")


###### "Stats"
Stats on the stamps fits files can be obtained through:

```python
Cube_Stat(cube2bplot_in,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency)
```

This will generate asciii and csv tables in the ```~/Example/Stack_Results-12CO-2D/TABLES/``` directory.

```
12CO-CII_HATLAS-RDS_B-0-stk-med-250kms-crc-15as_msk_ms-stt.dat
12CO-CII_HATLAS-RDS_B-0-stk-med-250kms-crc-15as_msk_ms-stt.csv
```

###### "MCMC"
It is possible to compute the Confidence Inteervals (CIs) of the flux measurements throough the 2D gauussian profiile on the composite stacked images through Monte Carlo process.

```
python
MCMC_generator(iterations_mc,line1,line2,method,error,'CNT_B',sbsmn,spc_wdt_dir=cw,mask_radi_as_ms=mask_radi_as_ms)
```

This will generate a series of tables in the ```~/Example/Stack_Results-12CO-2D/PLOTS/20/MCMC/``` directory containg the MCMC statistics:

```
- Sources-MC-10000-HATLAS-12CO-13CO-CNT_B-0-M3.dat
- CII_HATLAS-CNT_B-MS-3-Z-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-CNT_B-MS-3-FLX-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-CNT_B-MS-3-LUM-LOG-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-CNT_B-MS-3-CNT_B-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-FLX-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-LUM-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-LUM-LOG-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-CNT_B-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-CNT_B-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-Z-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-Z-0-0-stk-20kms-crc-15as_msk_ms-stt.dat
```

Tables foor plotting purposes in the ``` Stack_Results-12CO-2D/TABLES/PLOTS ``` directory:

```
- CII_HATLAS-12CO-13CO-CNT_B-MS-3-MC-10000-3-FLX-0-0-stk-20kms-crc-15as_msk_ms-stt-PLT.dat
```


And a plot containing the MCMC results.

![Alt text](./Figures-IM/Sources-MC-10000-HATLAS-12CO-13CO-CNT_B-0-M3.jpg?raw=true "Stacked spectra VALES field.")	
S
## Dependencies
Currently VSAT works only with astropy 2.0 as it relies on pyraf continuum task for continuum normalization. However a new version will be released dropping this dependency.
 - [astropy](https://www.astropy.org)
 - [bottleneck](https://pypi.org/project/Bottleneck/)
 - [pandas](https://pandas.pydata.org)
 - [scipy](https://www.scipy.org)
 - [numpy](https://numpy.org)
 - [lmfit](https://lmfit.github.io/lmfit-py/)
 - [matplotlib](https://matplotlib.org)
 - [termcolor](https://pypi.org/project/termcolor/)
 - [progressbar](https://pypi.org/project/progressbar2/)
## License

BSD 3-Clause License

Copyright (c) 2021, VSAT-1D developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

