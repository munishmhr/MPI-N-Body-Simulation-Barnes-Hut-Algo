#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    /*
        Reading from file and storing mass x and y axis in arrays
     */
    FILE *fp;
    FILE *fe;
    int i = 0,j;
    fp = fopen("tab65536", "r");
    double mass[200000];
    double x[200000];
    double y[200000];
    double z[200000];
    double vX[200000];
    double vY[200000];
    double vZ[200000];

    if (fp != NULL)
    {
        while (!feof(fp))
        {
            fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &mass[i], &x[i], &y[i], &z[i], &vX[i], &vY[i], &vZ[i]);
            //printf(" Reading %lf %lf %lf %lf %lf %lf %lf \n", mass[i], x[i], y[i], z[i], vX[i], vY[i], vZ[i]);
            i++;
        }
    }

    i = 0;
    for(i=0; i<65500; i++)
    {
        mass[i] = mass[i] +1;
        x[i] = x[i] + 1.9;
        y[i] = z[i] + 1.9;
        vX[i] = vX[i] + 1.9;
        vY[i] = vY[i] + 1.9;
        vZ[i] = vZ[i] + 1.9;
        z[i] = z[i] + 1.9;
    }

   // printf("\n%d", i);

    fe = fopen("tab65501","a+");
i=0;
    for(i=0; i<65500; i++)
    {
        //printf("%lf %lf %lf %lf %lf %lf %lf \n", mass[i], x[i], y[i], z[i], vX[i], vY[i], vZ[i]);
        fprintf(fe, "%lf %lf %lf %lf %lf %lf %lf \n", mass[i], x[i], y[i], z[i], vX[i], vY[i], vZ[i]);

    }



    return 0;
}
