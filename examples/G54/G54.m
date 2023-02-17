% test example G_54: graphene 54

cd ../../
DGESMAT_startup();
cd ./examples/G54

inputFile = 'G54.in';
outFile = "G54_statfile_Ecut40_9_4";

c = parcluster;
c.NumWorkers = 9;
c.NumThreads = 4;

c  % print some information of the parallel cluster

if isempty(gcp('nocreate'))
    parpool(c);  % start parallel pool (parpool)
else
    delete(gcp('nocreate'));  % shut down the last parallel pool
    parpool(c);  % start parallel pool (parpool)
end

spmd
fprintf('Worker %d says Hello World!\n', labindex)
end

dgdft_main(inputFile, outFile);

delete(gcp('nocreate'));