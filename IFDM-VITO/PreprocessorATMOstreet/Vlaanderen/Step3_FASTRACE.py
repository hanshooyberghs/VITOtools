import os

##################
## RUN FASTRACE ##
##################

cmd='java -Xmx5g -jar "fastrace-cli-1.1.0-20170516.135144-12.jar" "cf_RIO-IFDM-OSPM_2015_netwerk.txt"'
os.system(cmd)



####################
## POSTPROCESSING ##
####################

# Verwerking naar IFDM formaat + toevoegen tunnels / bruggen via standaardscript
cd 'VerwerkingOutput'
run Fastrace2IFDMOSPM.py '../Output/netwerk_all.shp' '../VerwerkteEmissies'