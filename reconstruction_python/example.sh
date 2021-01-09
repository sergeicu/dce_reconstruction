MRN=<MRN> # copy and paste correct value from `../input/subjectList_MRU.csv`. Remember to never store PHI data on github.   
name=<Name_Surname> # copy and paste correct value from `../input/subjectList_MRU.csv`. Remember to never store PHI data on github.   
outdir=output/
python -W ignore processSubjectRLTI.py -f -no_bmd -d $outdir -s $MRN -l 0.00001 -o reconReg/ -p 40 -stkBinWidth 0.01 -strtFrame 10 -pn $name --referenceSCAN ../input/reconResultsSCAN/
