## Sample program of analysis.


This code is intended for beginner, so calculation of error bar is not included in this program. 



In order to compile, do

  $make analysis



When you exeute this program, you have to specify at least 2 options, "-n and -output_dir".

ex.  $./analysis -n  10  -output_dir jet_energy

 directory would be generated under data/.



Here I summarize available options.
-----------------------------------------------------------------------
1) "-n " specifies number of file you want to read if your input files are numbering like "xxx_ver0.dat", "xxx_ver1.dat", "xxx_ver2.dat"....

2) "-output_dir" is used to name directory for output.

3) "--input_path" specifies path to input files.

4) "--ext" specifies extension of the input files. ex. "--ext .txt"

--------------------------------------------------------------------



Warning: This analysis program supposes that an input file contains just one event information. If that is not the case, you have to modify a little.





yuuka
