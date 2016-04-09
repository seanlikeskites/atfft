#include "../include/atfft.h"

int isEven (unsigned int x)
{
    return !(x % 2);
}

int isPowerOf2 (unsigned int x)
{
    return x && !(x & (x - 1));
}
