/**********************************************************************

  FFT.h

  Dominic Mazzoni

  September 2000

  This file contains a few FFT routines, including a real-FFT
  routine that is almost twice as fast as a normal complex FFT,
  and a power spectrum routine which is more convenient when
  you know you don't care about phase information.  It now also
  contains a few basic windowing functions.

  Some of this code was based on a free implementation of an FFT
  by Don Cross, available on the web at:

    http://www.intersrv.com/~dcross/fft.html

  The basic algorithm for his code was based on Numerical Recipes
  in Fortran.  I optimized his code further by reducing array
  accesses, caching the bit reversal table, and eliminating
  double-to-double conversions, and I added the routines to
  calculate a real FFT and a real power spectrum.

  Note: all of these routines use single-precision float.
  I have found that in practice, float work well until you
  get above 8192 samples.  If you need to do a larger FFT,
  you need to use doubles.

**********************************************************************/
/*
  Salvo Ventura - November 2006
  Added more window functions:
    * 4: Blackman
    * 5: Blackman-Harris
    * 6: Welch
    * 7: Gaussian(a=2.5)
    * 8: Gaussian(a=3.5)
    * 9: Gaussian(a=4.5)
*/

/*
  Lars Palo - December 2011 Modifications to original file from
  Audacity projects sources

  Changed all float to double, removed unneccessary
  experimental functions to fit better with LoopAuditioneers
  fft needs.
*/

/*
  Brooklyn Rose Ludlow - August 2023
  Removed remove things unused by LoopFinder
*/

#ifndef M_PI
#define	M_PI		3.14159265358979323846  /* pi */
#endif

/*
 * This is the function you will use the most often.
 * Given an array of doubles, this will compute the power
 * spectrum by doing a Real FFT and then computing the
 * sum of the squares of the real and imaginary parts.
 * Note that the output array is half the length of the
 * input array, and that NumSamples must be a power of two.
 */

void PowerSpectrum(int NumSamples, double *In, double *Out);

/*
 * Computes a FFT of complex input and returns complex output.
 * Currently this is the only function here that supports the
 * inverse transform as well.
 */

void FFT(int NumSamples,
         bool InverseTransform,
         double *RealIn, double *ImagIn, double *RealOut, double *ImagOut);

/*
 * Applies a windowing function to the data in place
 *
 * 0: Rectangular (no window)
 * 1: Bartlett    (triangular)
 * 2: Hamming
 * 3: Hanning
 * 4: Blackman
 * 5: Blackman-Harris
 * 6: Welch
 * 7: Gaussian(a=2.5)
 * 8: Gaussian(a=3.5)
 * 9: Gaussian(a=4.5)
 */

void WindowFunc(int whichFunction, int NumSamples, double *data);

void DeinitFFT();
