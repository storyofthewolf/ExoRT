import os
import sys
import subprocess
import argparse

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~ Populate 3d models directory of ExoRT ~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

parser = argparse.ArgumentParser()
parser.add_argument('--all',          action='store_true', help='for all 3d models')
parser.add_argument('--n28archean',   action='store_true', help='for n28archean 3d model')
parser.add_argument('--n42h2o',       action='store_true', help='for n42h2o 3d model')
parser.add_argument('--n68h2o',       action='store_true', help='for n68h2o 3d model')
parser.add_argument('--n68equiv',     action='store_true', help='for n68equiv 3d model')
parser.add_argument('--n84equiv',     action='store_true', help='for n84equiv 3d model')
parser.add_argument('--populate',     action='store_true', help='populate files')
parser.add_argument('--check',        action='store_true', help='check files for differences')
args = parser.parse_args()

# files needed for each 3d build
#src.main/
# exo_init_ref.F90
# exo_radiation_mod.F90  
# mcica_random_numbers.F90  
# planck_mod.F90         
# exo_radiation_cam_intr.F90  
# mcica.F90  
# rayleigh_data.F90  
# sys_rootdir.F90
#
#src.$model
# calc_opd_mod.F90  
# kabs.F90            
# radgrid.F90
# cloud.F90   
# initialize_rad_mod_cam.F90 
# model_specific.F90  
# rad_interp_mod.F90  
# spectral_output_cam.F90
#

                                                                    
models         = ['n28archean','n42h2o','n68h2o','n68equiv','n84equiv']
do_model       = [False, False, False, False, False]
if args.all:          do_model = [True, True,  True,  True,  True]
if args.n28archean:   do_model[0] = True
if args.n42h2o:       do_model[1] = True
if args.n68h2o:       do_model[2] = True
if args.n68equiv:     do_model[3] = True
if args.n84equiv:     do_model[4] = True

srcmodels    = ['','','','','']
threedmodels = ['','','','','']

#haze model not currently included here
#co2 ice nots included here

# master source code directories
srcmain          = '../source/src.main/'
i=0
for m in models:           
    srcmodels[i]        = '../source/src.'+ m + '/' 
    threedmodels[i]      = '../3dmodels/src.cam.' + m + '/'
    i=i+1


###############################################################3
# copy SourceMods files

#print("srcmain", srcmain)
#print("srcmodels", srcmodels)
#print("3dmodels", threedmodels)

model_specified = False
length=len(models)

#check files for differences
if args.check:
    print("checking files...")
    for i in range(length):
        if (do_model[i] == True):              
            print(" ")
            print("<<< <<< <<< <<< <<< <<< <<< ", models[i], " >>> >>> >>> >>> >>> >>> >>>")
            model_specified = True
            print("================ exo_init_ref.F90 ===================")
            f = ['diff', srcmain + 'exo_init_ref.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ exo_radiation_mod.F90 ==============")
            f = ['diff', srcmain + 'exo_radiation_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ mcica_random_numbers.F90 ===========")
            f = ['diff', srcmain + 'mcica_random_numbers.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ planck_mod.F90 =====================")
            f = ['diff', srcmain + 'planck_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ exo_radiation_cam_intr.F90 =========")
            f = ['diff', srcmain + 'exo_radiation_cam_intr.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ mcica.F90 ==========================")
            f = ['diff', srcmain + 'mcica.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ rayleigh_data.F90 ==================")
            f = ['diff', srcmain + 'rayleigh_data.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ sys_rootdir.F90 ===================+")
            f = ['diff', srcmain + 'sys_rootdir.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ calc_opd_mod.F90 ===================")
            f = ['diff', srcmodels[i] + 'calc_opd_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            print("=====================================================")
            print("================ kabs.F90 ===========================")
            f = ['diff', srcmodels[i] + 'kabs.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ radgrid.F90 ========================")
            f = ['diff', srcmodels[i] + 'radgrid.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ cloud.F90 ==========================")
            f = ['diff', srcmodels[i] + 'cloud.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ initialize_rad_mod_cam.F90 =========")
            f = ['diff', srcmodels[i] + 'initialize_rad_mod_cam.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ model_specific.F90 =================")
            f = ['diff', srcmodels[i] + 'model_specific.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ rad_interp_mod.F90 =================")
            f = ['diff', srcmodels[i] + 'rad_interp_mod.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")
            print("================ spectral_output_cam.F90 ============")
            f = ['diff', srcmodels[i] + 'spectral_output_cam.F90', threedmodels[i]]
            subprocess.run(f)             
            print("=====================================================")

# populate 3dmodels directories
if args.populate:
    print("copying files...")
    for i in range(length):
        if (do_model[i] == True):             
            print(" ")
            print("<<< <<< <<< <<< <<< <<< <<< ", models[i], " >>> >>> >>> >>> >>> >>> >>>")
            model_specified = True
            f = ['cp', srcmain + 'exo_init_ref.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'exo_radiation_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'mcica_random_numbers.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'planck_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'exo_radiation_cam_intr.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'mcica.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'rayleigh_data.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmain + 'sys_rootdir.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'calc_opd_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'kabs.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'radgrid.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'cloud.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'initialize_rad_mod_cam.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'model_specific.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'rad_interp_mod.F90', threedmodels[i]]
            subprocess.run(f) 
            f = ['cp', srcmodels[i] + 'spectral_output_cam.F90', threedmodels[i]]
            subprocess.run(f) 
            print("finished populating ",threedmodels[i])

if model_specified == False:
    print("no models were specified")
    print(models, "--all")
if args.check == False and args.populate == False:
    print("nothing done")
    print("--check, --populate")


