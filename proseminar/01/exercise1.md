## Exercise 1
### Description
This exercise consists in familiarizing yourself with SGE job submission.
You can find information about LCC2 at https://www.uibk.ac.at/zid/systeme/hpc-systeme/lcc2/ and information about SGE job submission at https://www.uibk.ac.at/zid/systeme/hpc-systeme/common/tutorials/sge-howto.html.
**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an _interactive job_ as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
### Tasks
#### Study how to submit jobs in SGE, how to check their state and how to cancel them.
SGE (Sun Grid Engine) is a system to submit a job to a computing cluster. You may also cancel jobs or check their current state. 
To start a job, you first need a job script, which is similar to a bash script, but with some SGE specific syntax. You start the job by calling `qsub script_name`. This in turn allocates the resources need to run the program, sets up the environment, executes the application and finally frees the allocation.
To cancel / delete a job, you need the the id (or ids) of the job/s you want to cancel: `qdel job_ids`. You can get the id of a job using `qstat`. 
#### Prepare a submission script that starts an arbitrary executable, e.g. `/bin/hostname`
This example is taken and modified from the example discussed in the first proseminar as well as the SGE tutorial on UIBK.
```bash
#!/bin/bash

# Select queue. There are three queues available. 
	# std.q is the general purpose and default queue with a maximum runtime of 240 hours. 
	# short.q is for small tests jobs, has a limited number of CPU slots and a max          	# runtime of 10 hours.
	# bigmem.q is for Leo3e. Max runtime of 240 hours with 512GB memory nodes.
#$ -q std.q

# The batch system should use the current directory as working directory. Makes it 
# easier to handle files 
#$ -cwd

# Name the job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N a_neat_job

# Redirect output stream to this file. Useful for debugging, simple output, etc.
#$ -o output.dat

# Join the error stream to the output stream.
#$ -j yes

# Specifies the parallel environment. A list is available with qconf -spl. In this case
# 2 slots per node, with 8 slots in total
#$ -pe openmpi-2perhost 8

# Load module because the software can be run on another system and not in the current 
# ssh session
module load openmpi/4.0.1

# Execute the software on 8 slots
mpiexec -n 8 /bin/hostname
```
* [SGE tutorial – Universität Innsbruck](https://www.uibk.ac.at/zid/systeme/hpc-systeme/common/tutorials/sge-howto.html#HDR1_1_2)
#### In your opinion, what are the 5 most important parameters available when submitting a job and why? What are possible settings of these parameters, and what effect do they have?
1. `-pe openmpi-Xperhost Y`.
2. `opath`
3. `epath`
4. Queue selection. This is important to not use more resources than needed.
5. `l`. Virtual memory
#### How do you run your program in parallel? What environment setup is required?
You need to specify how many cores / nodes you need within the job script file. You also need to load the module for each machine e.g. you need to add `module load openmpi/4.0.1` to any script running in parallel. Program execution then needs to be done using `mpiexec -n no_of_slots program_binary`