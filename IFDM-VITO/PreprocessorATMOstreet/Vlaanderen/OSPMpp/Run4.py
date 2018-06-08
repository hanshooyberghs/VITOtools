import os

nummers=range(32)

for nummer in nummers:
	# replace values in Lauchscript
	script='Launchscripts/script'+str(nummer+1)	
	cmd='cp Launchscripts/Script_general '+script
	os.system(cmd)	
	cmd="sed -i 's/nummer/"+str(nummer+1)+"/g' "+script
	os.system(cmd)
	
    # Launch
	cmd='qsub -l nodes=1:ppn=5 '+script
	os.system(cmd)
