#include <math.h>
#include <stdlib.h>
#include <atfft/atfft_dft.h>

int atfft_mod (int a, int n)
{
    int r = a % n;
    return r < 0 ? r + n : r;
}

int atfft_is_prime (unsigned int x)
{
    int sqrtX = (int) sqrt ((double) x);
    int i = 2;

    if (x == 1 || (atfft_is_even (x) && x > 2))
        return 0;

    for (i = 2; i <= sqrtX; ++i)
    {
        if ((x % i) == 0)
            return 0;
    }

    return 1;
}

int atfft_primitive_root_mod_n (unsigned int n)
{
    int g = 2;
    int isRoot = 1;

    /* this algorithm will only work for prime numbers */
    if (!atfft_is_prime (n))
        return -1;

    /* loop over potential primitive roots */
    for (g = 2; g < n; ++g)
    {
        int i = n - 2;
        int m = 1;
        isRoot = 1;

        while (i--)
        {
            m = (m * g) % n;

            if (m == 1)
            {
                isRoot = 0;
                break;
            }
        }

        if (isRoot)
            return g;
    }

    return -1;
}

void atfft_gcd (int a, int b, int *gcd, int *x, int *y)
{
    int absA = abs (a);
    int absB = abs (b);
    int r0 = ATFFT_MAX (absA, absB);
    int r1 = ATFFT_MIN (absA, absB);
    int s0 = absA > absB ? 1 : 0;
    int s1 = 1 - s0;
    int t0 = s1;
    int t1 = s0;

    while (r1)
    {
        int q = r0 / r1;

        int r2 = r0 - q * r1;
        int s2 = s0 - q * s1;
        int t2 = t0 - q * t1;

        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }

    *gcd = r0;
    *x = a < 0 ? -s0 : s0;
    *y = b < 0 ? -t0 : t0;
}

int atfft_mult_inverse_mod_n (int a, int n)
{
    int gcd, x, y;

    atfft_gcd (a % n, n, &gcd, &x, &y);

    if (gcd == 1)
        return atfft_mod (x, n);
    else
        return -1;
}
