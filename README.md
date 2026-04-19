The step 1 concentration code is used to convert raw bbNIRS signals from spectrometer from two channels left and right side of forehead. This is converted into concentrations of 
chromophores of Hbo, HHb, CCO and water. This stores in a file which can be sued in step 2 code.
The step 2 code takes input file of step 1 and splits the concentrations into three frequency bands, endogenic, neurogenic and myogenic. This then calculates
total hemoglobin and water free in each frequency bands. After calculating all paramters it is used to calculate correlation and slope in each freuquency bands between
HbT and water free and CCO and water free which is CSF.
