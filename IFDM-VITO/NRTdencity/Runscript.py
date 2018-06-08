## Run the NRT chain of ATMO-street for Antwerp.

import shutil,subprocess

## input
ifdm_command='ifdmmodel_1290M'
Runfolder='/projects/dencity/nrt/ifdm/Runfolder/'


## new meteo acquirement
########################
# Replace InputMeteoType1.txt

## new background acquirement
#############################
# Replace Aurora_conc

## coupling step
################
shutil.copy(Runfolder+'ifdm_coupling.conf',Runfolder+'ifdm.conf')
proc=subprocess.Popen('./'+ifdm_command,cwd=Runfolder,shell=True)
proc.wait()


## additive step
###############
shutil.copy(Runfolder+'ifdm_add.conf',Runfolder+'ifdm.conf')
proc=subprocess.Popen('./'+ifdm_command,cwd=Runfolder,shell=True)
proc.wait()

## ospm
#######
shutil.copy(Runfolder+'ifdm_ospm.conf',Runfolder+'ifdm.conf')
proc=subprocess.Popen('./'+ifdm_command,cwd=Runfolder,shell=True)
proc.wait()


## Postprocessing
#################
pol_list=['NO2','PM10','PM25','BC']
for pollutant in pol_list:
    proc=subprocess.Popen('python Postprocessing.py '+pollutant,shell=True)
    proc.wait()
