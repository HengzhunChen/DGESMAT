DGESMAT (Discontinuous Galerkin Electronic Structure in MATLAB)
==============

To use this code, cd into the dgesmat directory and type `DGESMAT_startup`
in command line or run DGESMAT_startup.m to set up all paths.

You can then cd into the examples directory, and run the corresponding 
example for some simple models. Templates of input configuration file can 
be found in folder src/ESDF.

Currently, DGESMAT contains two submodules:
- discontinuous Galerkin density functional theory (DGDFT)
- plane wave density functional theory (PWDFT)

At this point (August, 2022), DGESMAT appears to need to MATLAB version 
R2021b or higher. Some earlier versions may not work. 

In current implementation, `parfor` is used to save the computation time in 
submodule DGDFT. You may need to set up the parallel settings according to 
the particular question you are trying to solve. Here is a sample code

``` matlab
c = parcluster;
c.NumWorkers = 9;  % number of parallel workers
c.NumThreads = 4;  % number of threads each worker used

if isempty(gcp('nocreate'))
    parpool(c);  % start parallel pool (parpool)
else
    delete(gcp('nocreate'));  % shut down the last parallel pool
    parpool(c);  % start parallel pool (parpool)
end

% computation code

delete(gcp('nocreate'));  % shut down the parallel pool
```

You can cd into the examples directory and see the graphene 54 (G54) mode 
for more details. See also document of `parfor` from MATLAB. You can also
change `parfor` into `for` to remove the parallel computing tools whenever 
you want.

NOTE: To know more infomation about your machine, running `feature('numcores')` 
without assigning its output displays some general debugging information, 
an example is as follow:

```
MATLAB detected: 32 physical cores.
MATLAB detected: 64 logical cores.
MATLAB was assigned: 36 logical cores by the OS.
MATLAB is using: 32 logical cores.
MATLAB is not using all logical cores because hyper-threading is enabled.
```

By default, MATLAB and Parallel Computing Toolbox consider only real cores, 
not hyperthreaded cores. You can override this choice in Parallel Computing 
Toolbox by modifying your local cluster profile.
As the performance is dependent on the problem, here we would recommend 
measuring the performance of your code using different number of workers to 
see if hyperthreading provides any benefit in your case. 

For more details, see also [Definitive answer for hyperthreading and the Parallel Computing Toolbox (PCT)](https://www.mathworks.com/matlabcentral/answers/80129-definitive-answer-for-hyperthreading-and-the-parallel-computing-toolbox-pct#answer_89845).

Here is a plot of run time with respect to the number of threads for 
packages DGDFT (implemented with MPI) and DGESMAT. The mode used in the 
figure is graphene 54 (G54) with element partition `(1, 3, 3)`, Ecut `40` 
and ONCV type pseudopotential, number of SCF outer iteration is `30`. 
Number of threads in horizontal axis is number of processes used in MPI 
or product of number of workers and number of threads each worker in MATLAB. 
In the following example, number of workers are all `9`, which corresponding 
to the number of elements total. 

<center>
<img src="./doc/figure/performace_comparison.png" width="70%">
</center>

From the figure we can see that with suitable number of threads the run time 
of DGESMAT will be about twice of DGDFT for the G54 mode, which is accpetable 
for prototype development of algorithms.