import subprocess
for i in range(1, 9):
    job_name = "pi_omp"
    output_filename = job_name + f"{i}threads.out"
    howmany_perwhathost = f"openmp {i}"
    what2make = "omp"

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
export OMP_NUM_THREADS={i}

for ((i=1000;i<=1000000000;i*=10)); do
	./{job_name} "$i"
done
	"""
    with open(f"{job_name}_{i}threads.script", "w") as outfile:
        outfile.write(out_txt)
        subprocess.run(["qsub", out_txt])