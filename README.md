OpenMHD
=======

OpenMHD is a collection of magnetohydrodynamic (MHD) simulation codes.
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

    $ ipython3
    In [1]: %run -i plot.py

Python (Jupyter Lab):

    $ jupyter-lab plot.ipynb

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

 * [C. W. Shu, S. Osher, J. Comput. Phys. 77, 439 (1988)](https://doi.org/10.1016/0021-9991(88)90177-5)
 * [B. van Leer, J. Comput. Phys. 23, 276 (1977)](https://doi.org/10.1016/0021-9991(77)90095-X)
 * [T. Miyoshi, K. Kusano, J. Comput. Phys. 208, 315 (2005)](https://doi.org//10.1016/j.jcp.2005.02.017)
 * [K. F. Gurski, SIAM J. Sci. Comput. 25, 2165 (2004)](https://doi.org/10.1137/S1064827502407962)
 * [A. Dedner, F. Kemm, D. Kroner, C.-D. Munz, T. Schnitzer, M. Wesenberg, J. Comput. Phys. 175, 645 (2002)](https://doi.org/10.1006/jcph.2001.6961)

1D Riemann problem

 * [M. Brio, C. C. Wu, J. Comput. Phys. 75, 400 (1988)](https://doi.org/10.1016/0021-9991(88)90120-9)

Orszag-Tang vortex

 * [S. A. Orszag, C.-M. Tang, J. Fluid Mech. 90, 129 (1979)](https://doi.org/10.1017/S002211207900210X)

Kelvin-Helmholtz instability

 * [Y. Matsumoto, M. Hoshino, Geophys. Res. Lett. 31, L02807 (2004)](https://doi.org/10.1029/2003GL018195)

Magnetic reconnection

 * [S. Zenitani, T. Miyoshi, Phys. Plasmas 18, 022105 (2011)](https://doi.org/10.1063/1.3554655)
 * [S. Zenitani, Phys. Plasmas 22, 032114 (2015)](https://doi.org/10.1063/1.4916104)

Relativistic jet problem

 * [S. Zenitani, M. Hesse, A. Klimas, Astrophys. J. 712, 951 (2010)](https://doi.org/10.1088/0004-637X/712/2/951)
 * [A. Mignone, M. Ugliano, G. Bodo, MNRAS, 393, 1141 (2009)](https://doi.org/10.1111/j.1365-2966.2008.14221.x)
 * [A. Mignone, J. C. McKinney, MNRAS, 378, 1118 (2007)](https://doi.org/10.1111/j.1365-2966.2007.11849.x)

Loop advection problem

 * [T. A. Gardiner, J. M. Stone, J. Comput. Phys. 227, 4123 (2008)](https://doi.org/10.1016/j.jcp.2007.12.017)

