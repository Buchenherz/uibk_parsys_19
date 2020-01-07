names = ["no","mpi","omp","both"]
for i in range(1, 9):
    for j in range(1, 5):
        for k in range(0,4):
            job_name = "hs2d_hybrid_" + names[k]
            output_filename = job_name + f"{names[k]}_{i}threads_{j}ranks.out"
            howmany_perwhathost = f"openmp {i} openmpi-{j}perhost {j}"
            what2make = names[k]

            out_txt = f"""#!/bin/bash

# Select queue.
#$ -q std.q

# Use the current wd
#$ -cwd

# Name the job
#$ -N {job_name}{i}{j}

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

for ((i=100;i<=500;i+=100)); do
	mpiexec -n {j} ./{job_name} "$i" "$i"
done
	"""
            with open(f"{job_name}_{i}threads_{j}ranks.script", "w") as outfile:
                outfile.write(out_txt)
