I'll write some fully fledged documentation here at some point. For now just
know that you should create an atfft structure with atfft_create() and pass that
around to the transform functions (atfft_complex_transform(),
atfft_real_forward_transform() and atfft_real_backward_transform()). The layout of the input and
ouput of the transforms is identical to FFTW's.

Enjoy!
