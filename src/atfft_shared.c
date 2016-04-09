#include "../include/atfft.h"

int isPowerOf2 (unsigned int x)
{
    return x && !(x & (x - 1));
}
