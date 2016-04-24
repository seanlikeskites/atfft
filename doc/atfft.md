I'll wirte some fully fledges documentation here at some point. For now just
know that you hould create an atfft structure with atfft_create() and pass than
around the the transform functions (atfft_complex_transform(),
atfft_real_forward_transform() and atfft_real_backward_transform()). The layout of the input and
ouput of the transforms is identical to FFTW's.

Enjoy!
