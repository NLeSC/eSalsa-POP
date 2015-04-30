eSalsa POP GPU Version
======================

This is the user guide of GPU-POP (the GPU-enabled version of the Parallel Ocean Program).

Unfortunately, only the computation of vertical mixing coefficients happens on the GPU as of now. Hopefully, we 
will be able to port more of POP in the future.

How to get it?
--------------

You can get the GPU version using the following command:

git clone git@github.com:NLeSC/eSalsa-POP.git --branch gpu

This will create a directory called eSalsa-POP that will become the POPDIR for the GPU version. There are other ways to download the code, for example directly from 
<https://github.com/NLeSC/eSalsa-POP/tree/gpu> But using the git command above, it is easier to stay up to date with the latest version of the GPU version.


Required settings
-----------------

It's **extremely important** to realize that the GPU version will only produce correct results for 0.1 degree resolution runs with the following namelist settings:

&vmix_kpp_nml  
   Prandtl         = 10.0  
   rich_mix        = 50.0  
   lrich           = .true.  
   ldbl_diff	   = .true.  
   lshort_wave     = .true.  
   lcheckekmo	   = .false.  
   lhoriz_varying_bckgrnd = .false.  
   llangmuir              = .false.  
   linertial              = .false.  
   num_v_smooth_Ri = 1  
/  
(the false values in vmix_kpp_nml are actually default values so can be left out of the namelist)

&grid_nml  
  partial_bottom_cells = .true.  
/  

&tidal_nml  
  ltidal_mixing = .false.  
/  
(this is the default value, namelist can be left empty)

&sw_absorption_nml  
  sw_absorption_type   = 'jerlov'  
   jerlov_water_type    =    3  
/  
(these are the default values, namelist can be left empty)

&state_nml  
   state_choice = 'mwjf'  
   state_range_opt = 'enforce'  
/  

There is no fundamental reason why only this specific list of settings is supported. It was mainly done to save development time. If you really need an option different from what is listed here, 
let me know.


Configuring the GPU Version
---------------------------

There are a couple of things that are specific to the GPU version:

First, in pop_in you should have a namelist like this:

&gpu_mod_nml  
  use_gpu = .true.  
/  

And secondly, if you want to use a block size other than 60 by 60 you should also change the gpu_domain.h file in the directory 'source'. It's important that nx_block and ny_block are set to +4 of 
the values you use in POP_DomainSizeMod.F90. I know this is not very user friendly, but there is no simple way to prevent this without introducing yet another settings file.

Actually, I recommend to not change the block size to anything other than 60x60 for the 0.1 degree resolution. Larger blocks will allow for less land-only blocks to be 
removed, smaller blocks will increase the total amount of work because the increase in removed land-only blocks is outweighed by the increase in the number of halo cells around blocks. 60x60 is the 
optimal setting for the 0.1 degree resolution and present-day landmask.


Building GPU-POP
----------------

The build file used on Cartesius is called bull-gpu-intel.gnu and is located in the directory build. Please set POPARCH to bull-gpu-intel, POPDIR to eSalsa-POP, POPEXEDIR to your run directory.




