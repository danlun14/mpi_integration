#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define RUNGE_BASE4 1.379768
#define MONTE_CARLO_BASE4_1 0.164387
#define MONTE_CARLO_BASE4_2 1.646971

double function(double x)
{
    return sqrt(x * (3 - x)) / (x + 1);
}

void Runge(int rank, int commsize)
{
    double a = 1;
    double b = 1.2;
    int n0 = 100000000;
    double eps = 1E-6;

    if (rank == 0)
    {
        printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
    }

    double sq[2];

    int n = n0;
    int k;
    double delta = 1;

    for (k = 0; delta > eps; n *= 2, k ^= 1)
    {
        int points_per_proc = n / commsize;
        int lb = rank * points_per_proc;
        int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);
        double h = (b - a) / n;

        double s = 0.0;

        sq[k] = 0.0;
        for (int i = lb; i <= ub; i++)
        {
            s += function(a + h * (i + 0.5));
        }

        MPI_Allreduce(&s, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sq[k] *= h;

        if (n > n0)
        {
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }

    if (rank == 0)
    {
        printf("Result: %.12f; Runge rule: EPS %e, n %d\n", sq[k], eps, n / 2);
    }
}

double getrand(unsigned int *seed)
{
    return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
    return exp(x - y);
}

void Monte_Carlo(int rank, int commsize, int n)
{
    srand(rank);
    int in = 0;
    double s = 0;

    double s_loc = 0.0;
    int in_loc = 0;
    unsigned int seed = rank;

    int points_per_proc = n / commsize;
    int lb = rank * points_per_proc;
    int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);

    for (int i = lb; i < ub; i++)
    {
        double x = getrand(&seed) * (-1);
        double y = getrand(&seed);
        in_loc++;
        s_loc += func(x, y);
    }

    MPI_Reduce(&in_loc, &in, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s_loc, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        double v = (1.0 * in) / n;
        double res = (v * s) / in;
        printf("Result: %.12f, n %d \n", res, n);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int commsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    double t = 0;

    t -= MPI_Wtime();
    Runge(rank, commsize);
    t += MPI_Wtime();

    double Tmax;
    MPI_Reduce(&t, &Tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Elapsed time Runge(sec.): %.6f\n", Tmax);
        printf("%d %.6f\n", commsize, RUNGE_BASE4 / Tmax);
    }

    for (int i = 0; i < 2; i++)
    {
        int n = pow(10, 7 + i);

        t = 0;

        t -= MPI_Wtime();
        Monte_Carlo(rank, commsize, n);
        t += MPI_Wtime();

        MPI_Reduce(&t, &Tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            printf("Elapsed time Monte_Carlo n = %d(sec.): %.6f\n", n, Tmax);
            if (i == 0)
            {
                printf("%d %.6f\n", commsize, MONTE_CARLO_BASE4_1 / Tmax);
            }
            else
            {
                printf("%d %.6f\n", commsize, MONTE_CARLO_BASE4_2 / Tmax);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
