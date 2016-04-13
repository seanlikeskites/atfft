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

* [FFTW3](www.fftw.org)
  - Swift as the wind!
  - GPL Licenced.

* [FFTS](www.github.com/anthonix/ffts)
  - Nicely speedy.
  - Signal lengths of 2^n only.
  - Single precision only.
  - Only builds on 64bit machines.
  - BSD Licenced.

* [Ooura fft4g](www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
  - Nice and simple (one source file).
  - Signal lengths of 2^n only.
  - Double precision only.
  - Very permissive licence.

* [GSL](www.gnu.org/software/gsl)
  - Part of a much larger library of numerical goodies.
  - Double precision only.
  - GPL Licenced.

* [Kiss FFT] (kissfft.sourceforge.net) 
  - Simple, as the name might suggest.
  - BSD Licenced.
  - Only does complex transforms, real transforms 
    are provided by atfft by copying from complex
    transforms.
  - atfft currently only supports Kiss FFT 
    built with -Dkiss_fft_scalar=double.

* [FFmpeg](www.ffmpeg.org)
  - Part of libavcodec.
  - Signal lengths of 2^n only.
  - Single precision only.
  - LGPL or GPL Licenced depending on how it's built.

* [vDSP](developer.apple.com/library/mac/documentation/Accelerate/Reference/vDSPRef)
  - Apple only.
  - Signal lengths of 2^n (and some multiples thereof) only.

## Licence
Copyright (C) 2016 Sean Enderby <sean.enderby@gmail.com>

This program is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it 
and/or modify it under the terms of the Do What The Fuck You Want 
To Public License, Version 2, as published by Sam Hocevar. See 
the COPYING file for more details.

