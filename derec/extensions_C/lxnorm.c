#include <math.h>
#include <stdio.h>

double lxnorm_m(double u[], double v[], double norm, int len)
    {
        int i;
        double msum=0.;
        for (i=0; i<len; i++ )
        {
            msum = msum + fabs(pow(v[i] - u[i], norm));
        }
        return pow(msum, 1./norm);
    }

double lxnorm_n(double v[], double norm, int len)
    {
        int i;
        double nsum=0.;
        for (i=0; i<len; i++ )
        {
            nsum = nsum + fabs(pow(v[i], norm));
        }
        return pow(nsum, 1/norm);
    }

