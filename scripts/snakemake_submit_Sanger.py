#!/software/python-3.5.0/bin/python3
import os
import sys

from snakemake.utils import read_job_properties

#Extraxt job properties
jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)
print(job_properties)

# do something useful with the threads
threads = job_properties["threads"]
job_name = job_properties["rule"]

#Extract memory and construct memory string
mem = str(job_properties.get('resources',{}).get('mem'))
memory_string = ' -R"span[hosts=1] select[mem>' + mem +'] rusage[mem=' + mem +']" -M '+mem

#Construct output file
output_file = os.path.join("FarmOut", job_name + '.%J.txt')

#Construct a submission script
command = "".join(["bsub -G team170 -o ", output_file, " -q ", "normal", 
	" -n ", str(threads), memory_string, " -J ", job_name])

#Run the script
full_command = "{cmd} {script}".format(cmd = command, script=jobscript)
os.system(full_command)
