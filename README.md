OpenMHD
=======

OpenMHD is a collection of one/two-dimensional magnetohydrodynamic (MHD) codes.
The latest version is available at the following URLs.

 * <https://sci.nao.ac.jp/MEMBER/zenitani/openmhd-e.html>
 * <https://sci.nao.ac.jp/MEMBER/zenitani/openmhd-j.html>

Usage
-------

1. Edit the Makefile in a project directory.
2. Configure the main file (usually main.f90), the initial configuration file
   (model.f90), and the boundary condition file (bc.f90) appropriately.
3. Make and run the program.
4. The program will output many files in the "data/" directory.
5. Analyze the data. Sample Python and IDL scripts are provided.
   One can use them in the following way.

Python (matplotlib):

    $ ipython3 --pylab
    In [1]: %run -i plot.py

Python (Jupyter):

    $ jupyter-notebook plot.ipynb

IDL:

    $ idl
    IDL> .r batch

License
---------

All files are distributed under the term of the GNU General Public License version 3 or later.
See the COPYING file.


References
-------------

Numerical schemes

 * C. W. Shu, S. Osher, J. Comput. Phys. 77, 439 (1988)
 * B. van Leer, J. Comput. Phys. 23, 276 (1977)
 * T. Miyoshi, K. Kusano, J. Comput. Phys. 208, 315 (2005)
 * K. F. Gurski, SIAM J. Sci. Comput. 25, 2165 (2004)
 * A. Dedner, F. Kemm, D. Kroner, C.-D. Munz, T. Schnitzer, M. Wesenberg, J. Comput. Phys. 175, 645 (2002)

1D Riemann problem

 * M. Brio, C. C. Wu, J. Comput. Phys. 75, 400 (1988)

Corotating interaction region

 * K. Tsubouchi, J. Geophys. Res. 114, A02101 (2009)

Orszag-Tang vortex

 * S. A. Orszag, C.-M. Tang, J. Fluid Mech. 90, 129 (1979)

Kelvin-Helmholtz instability

 * Y. Matsumoto, M. Hoshino, Geophys. Res. Lett. 31, L02807 (2004)

Magnetic reconnection

 * S. Zenitani, T. Miyoshi, Phys. Plasmas 18, 022105 (2011)
 * S. Zenitani, Phys. Plasmas 22, 032114 (2015)

Relativistic jet and relativistic MHD schemes

 * S. Zenitani, M. Hesse, A. Klimas, ApJ, 712, 951 (2010)
 * A. Mignone, M. Ugliano, G. Bodo, MNRAS, 393, 1141 (2009)
 * A. Mignone, J. C. McKinney, MNRAS, 378, 1118 (2007)
