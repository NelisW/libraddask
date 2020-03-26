__author__ = 'DGriffith'
import librad
import glob
# Find and read each uvspec example case in the examples folder
uvINPfiles = glob.glob('examples/UVSPEC_*.INP')
for uvINPfile in uvINPfiles:
    print('Processing ' + uvINPfile)
    uvTestCase = librad.Case(filename=uvINPfile)  # Read the uvspec input file
    try:
        uvTestCase.readout()  # Read the corresponding uvspec output file
    except:
        print('Reading of output failed for ' + uvINPfile)




