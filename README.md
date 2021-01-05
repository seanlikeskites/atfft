# All The FFT!

The purpose of this project is to provide a permissively licenced
interface to various C libraries which calculate Discrete Fourier
Transforms.

In an ideal world one would choose their favorite FFT library
and use it for every project they ever work on. Unfortunately 
we live in a work of conflicting software licences and
unsupported platforms, so at some point most of us will have 
to stray from our preferred implementations. Providing a uniform
interface to as many libraries as possible should help to ease 
this pain.

## Implementations
Here is some useful information about the libraries usable via
the atfft interface:

* [FFTW3](http://www.fftw.org)
  - Swift as the wind!
  - GPL Licenced.

* [FFTS](http://www.github.com/anthonix/ffts)
  - Nicely speedy.
  - Signal lengths of 2^n only.
  - Single precision only.
  - Only builds on 64bit machines.
  - BSD Licenced.

* [Ooura fft4g](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
  - Nice and simple (one source file).
  - Signal lengths of 2^n only.
  - Double precision only.
  - Very permissive licence.

* [GSL](http://www.gnu.org/software/gsl)
  - Part of a much larger library of numerical goodies.
  - Double precision only.
  - GPL Licenced.

* [Kiss FFT](http://kissfft.sourceforge.net) 
  - Simple, as the name might suggest.
  - BSD Licenced.
  - Real transforms only support even lengths.

* [FFmpeg](http://www.ffmpeg.org)
  - Part of libavcodec.
  - Signal lengths of 2^n only.
  - Single precision only.
  - LGPL or GPL Licenced depending on how it's built.

* [vDSP](http://developer.apple.com/library/mac/documentation/Accelerate/Reference/vDSPRef)
  - Apple only.
  - Single and double precision only (no long doubles).
  - Signal lengths of 2^n (and some multiples thereof) only.

* [PFFFT](http://www.bitbucket.org/jpommier/pffft)
  - Single precision only.
  - Transform length must be multiple of 32 (16 for complex transforms.)
  - BDS-like licence.

## Links
* [Documentation](https://seanlikeskites.github.io/atfft/doc/html)
* [Benchmarks](https://seanlikeskites.github.io/atfft/benchmarks/benchmark.html)

## Licence
Copyright (c) 2021 Sean Enderby <sean.enderby@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
