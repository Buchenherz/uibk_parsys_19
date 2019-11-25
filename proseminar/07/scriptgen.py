for i in range(1,9):
	job_name = "nbody2D_omp"
	output_filename = job_name + f"{i}threads.out"
	howmany_perwhathost = f"openmp {i}"
	what2make = "seq"
	
	out_txt = f"""#!/bin/bash

# Select queue.
#$ -q std.q

# Use the current wd
#$ -cwd

# Name the job
#$ -N {job_name}

# Redirect output stream to this file. Useful for debugging, simple output, etc.
#$ -o {output_filename}

# Redirect error stream to this file.
#$ -e error.dat

# Join the error stream to the output stream.
#$ -j yes

# Specifies the parallel environment. A list is available with qconf -spl. 
#$ -pe {howmany_perwhathost}

# check whether the syntax of the job is okay (do not submit the job)
# #$ -w v

# Load module because the software can be run on another system and not in the current
# ssh session
# load module

module load gcc
make {what2make}

for ((i=1000;i<=10000;i+=1000)); do
	./2D_n-body_simulation_seq_omp 100 100 "$i" {i}
done
	"""
	with open(f"nbody2D_omp_{i}threads.script", "w") as outfile:
		outfile.write(out_txt)