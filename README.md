======================================================================
The implementation of Visual Difference Predicator (VDP) and
Visual Difference Predicator for High Dynamic Range images (HDR VDP)
======================================================================

The standard VDP implementation was based on:

Scott Daly,
The Visible Differences Predictor: An algorithm for the
assessment of image fidelity.
Digital Image and Human Vision, Cambridge, MA: MIT Press,
A. Watson, Ed., 179-206. (1993)

and

Scott Daly,
The visible differences predictor: an algorithm for the
assessment of image fidelity,
SPIE Vol. 1666 Human Vision, Visual
Processing, and Digital Display III 1992

and

US patent US5,394,483

and some comments from Scott Daly

Note that not all parameters for the model are known so there may be
some discrepancies between this implementation and the orignal
implementation from that paper.

The HDR version (HDR VDP) was based on:

Rafal Mantiuk, Scott Daly, Karol Myszkowski and Hans-Peter Seidel.
Predicting Visible Differences in High Dynamic Range Images - Model and its Calibration.
In: Proc. of Human Vision and Electronic Imaging X, IS&T/SPIE's 17th Annual Symposium on Electronic Imaging 2005. pp. 204-214

and

Rafal Mantiuk, Karol Myszkowski and Hans-Peter Seidel
Visible Difference Predicator for High Dynamic Range Images
In: Proc. of IEEE International Conference on Systems, Man and Cybernetics, 2004. pp. 2763-2769

The sources should compile under all variants of Linux and
cygwin. Compilation under OSX may require some small changes.

Since the original version of VDP has not been maintained for quite
some time, I would strongly recommend using HDR VDP. HDR VDP has the
same capabilities as VDP, but should give more consistent predictions
for much larger dynamic range. It also includes aspect that are
missing in the original VDP, such as the optics of the eye (OTF).

License:
--------

The source code is available under General Public Licence (GPL) (see
COPYING file). If you find this software useful and you have used it
for your project or paper, please state which version you are using
and cite the following paper:

Rafal Mantiuk,Scott Daly, Karol Myszkowski, Hans-Peter Seidel
Predicting Visible Differences in High Dynamic Range Images -
Model and its Calibration
In: Proc. of Human Vision and Electronic Imaging X, IS&T/SPIE's 17th Annual Symposium on Electronic Imaging 2005. pp. 204-214

@inproceedings{mantiuk:2005:PredVisDiff,
  author = {Mantiuk, Rafa{\l} and Daly, Scott and Myszkowski, Karol and Seidel, Hans-Peter},
  editor = {Rogowitz, Bernice E. and Pappas, Thrasyvoulos N. and Daly, Scott J.},
  title = {Predicting Visible Differences in High Dynamic Range Images - Model and its Calibration},
  booktitle = {Human Vision and Electronic Imaging X, IS\&T\/SPIE's 17th Annual Symposium on Electronic Imaging (2005)},
  year = {2005},
  volume = {5666},
  pages = {204--214},
  isbn = {0277-786X},
}

Dependencies:
------------

The code requires two libraries to compile:

* For Fast Fourier Transformation.
  fftw 3.0 and later (fftw3f - 32bit float version)
  http://www.fftw.org/
Note: add --enable-float switch when running ./configure script. Otherwise
  32bit float version of the library won't be generated.

* To read and write HDR images
  pfstools
  http://pfstools.sourceforge.net/

Besides the libraries, netpbm package is needed to write .png files,
which are the output of the 'vdp' frontend.

http://netpbm.sourceforge.net/

Mac:

    brew install fftw

Ubuntu:

   apt-get install 

Compilation:
------------

This software uses standard automake and autoconf scripts to configure,
compile, and install. To configure issue a command:

./configure

The default installation location is /usr/local. If you want to
specify different location, give a parameter:

./configure --prefix=$HOME/local

To compile:

make

To install (for /usr/local must be root user):

make install

If are compiling sources from the CVS, run ./bootstrap script to
generate ./configure script.

Documentation:
--------------

Detailed description how to use vdp can be found in the manual
pages. To view documentation, check 'man vdp', man 'vdpcmp', and 'man
vis'.


Low Dynamic Range Images
------------------------

Starting from the version 1.3 the vdp script can automatically
recognize low dynamic range images and convert them from sRGB to
XYZ. See the man pages of vdp and vdpcmp for more information.

Note
-----------

This implementation of the Visual Difference Predicator (VDP) was
initially based on Karol Myszkowski's SGI code. However, it was
progressively ported, redesigned and reimplemented.

