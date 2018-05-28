INFILE="PNP_model_boundarylayer_test.mph"
OUTFILE="PNP_model_boundarylayer_test_out.mph"
comsol batch -np 4 -inputfile $INFILE -outputfile $OUTFILE
