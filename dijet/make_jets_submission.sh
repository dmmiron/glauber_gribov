Initialdir = /atlas/work/dmmiron/dijet

Universe = vanilla
GetEnv = true
Rank = Mips

Executable = /atlas/u/dmmiron/glauber_gribov/dijet/runjobs.sh 
Output = $(Initialdir)/logs/c.$(Cluster).$(Process).out
Error = $(Initialdir)/logs/c.$(Cluster).$(Process).err
Log = $(Initialdir)/logs/c.$(Cluster).$(Process).log
Priority = -10

#process gets turned into b, and phi values to sweep over
# parameters are process, nphi, nsamples, alpha, path, pairs (Bool t/f)
Arguments = $(Process) 7 2000 0 ~/glauber_gribov/dijet/test 0
Queue 112
