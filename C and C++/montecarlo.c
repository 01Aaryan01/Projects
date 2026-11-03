#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* comparison function for qsort */
int compare(const void *a, const void *b) {
    double x = *(double *)a;
    double y = *(double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

int main() {
    int N;
    printf("Enter number of simulations: ");
    scanf("%d", &N);

    double *sim = malloc(N * sizeof(double));
    if (sim == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    srand(time(NULL));  // seed RNG

    for (int i = 0; i < N; i++) {
        double rand_val = (double)rand() / RAND_MAX;
        double cost_impact;

        if (rand_val < 0.15)
            cost_impact = 5;
        else if (rand_val < 0.35)
            cost_impact = 10;
        else if (rand_val < 0.45)
            cost_impact = 30;
        else if (rand_val < 0.70)
            cost_impact = 15;
        else if (rand_val < 0.85)
            cost_impact = 40;
        else
            cost_impact = 0;

        sim[i] = 100 + cost_impact;
    }

    /* sort results */
    qsort(sim, N, sizeof(double), compare);

    int p90_index = (int)(0.9 * N) - 1;
    int p80_index = (int)(0.8 * N) - 1;

    double p90_cost = sim[p90_index];
    double p80_cost = sim[p80_index];

    /* mean */
    double mean = 0.0;
    for (int i = 0; i < N; i++)
        mean += sim[i];
    mean /= N;

    /* variance and std deviation */
    double variance = 0.0;
    for (int i = 0; i < N; i++)
        variance += (sim[i] - mean) * (sim[i] - mean);
    variance /= N;

    double std_dev = sqrt(variance);

    printf("P90 Cost: %.2f\n", p90_cost);
    printf("P80 Cost: %.2f\n", p80_cost);
    printf("Mean Cost: %.2f\n", mean);
    printf("Standard Deviation: %.2f\n", std_dev);

    free(sim);
    return 0;
}
