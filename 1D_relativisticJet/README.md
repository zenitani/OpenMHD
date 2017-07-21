Relativistic jet
================

One dimensional relativistic magnetohydrodynamic (RMHD) code
without the normal magnetic field (B_x=0).
Since the 2010 paper, primitive variables have been changed
from v to u = gamma v.


Usage
-------

1. Edit Makefile.
2. Edit main.f90 and model.f90 for paraters.
3. If not found, create an output directory "data/".
4. Make and run the program.

    $ make
    $ ./a.out


License
-----------

All files are distributed under the term of the GNU General Public License version 3 or later.
See the COPYING file.


References
-----------

 * S. Zenitani, M. Hesse, A. Klimas, ApJ, 712, 951 (2010)
 * A. Mignone, M. Ugliano, G. Bodo, MNRAS, 393, 1141 (2009)
 * A. Mignone, J. C. McKinney, MNRAS, 378, 1118 (2007)
