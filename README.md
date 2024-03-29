# Force_Clamp
Matlab files for processing muscle force clamp data files.

M-File which reads data files (either new .dat or old .fcp
files) from muscle force clamp experiments with three steps in
muscle lengths.

Users pick the end of the muscle clamp and the data 20 ms
before the end to the picked end is averaged.  The muscle lengths
are converted to fractions of the total muscle length.  Muscle
velocities are calculated from the muscle length time histories.

The muscle forces are corrected for the passive muscle
forces using two different methods.  A Hill's equation model is
used to fit the muscle velocities and forces.

Plots of the muscle velocity and muscle power with muscle
force are plotted for review before saving the data and fit
parameters to a MS-Excel spreadsheet.

The user selects whether to create a new spreadsheet file or
append to an user selected existing spreadsheet file.

NOTES:

1.  If the user selects an existing output MS-Excel
 spreadsheet, the spreadsheet cannot be open in another
 program (e.g. MS-Excel, text editor, etc.) while using
 this program.

2.  M-files get_tens.m, rd_dat.m, rd_fcp.m, and val2idx.m
 must be in the current path or directory.

