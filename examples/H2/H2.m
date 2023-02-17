% test example H2 with 8 atoms

inputFile = "H2_pwdft.in";
outFile = "statfile_pwH2";

pwdft_main(inputFile, outFile);


inputFile = "H2_dgdft.in";
outFile = "statfile_dgH2";

dgdft_main(inputFile, outFile);
