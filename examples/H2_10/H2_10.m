% test example H2_10, 10 is ecut

inputFile = "pwTestH2.in";
outFile = "statfile_pwH2";

pwdft_main(inputFile, outFile);


inputFile = "dgdft.in";
outFile = "statfile_dgH2";

dgdft_main(inputFile, outFile);
