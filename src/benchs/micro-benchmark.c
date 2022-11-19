#include<stdlib.h>
#include<time.h>
#include<stdio.h>
void main(){
    int i = 0;
    float* A = (float*) malloc(sizeof(float) *100000000);
    float acc;
    struct timespec start, curr;
    for(int n = 100; n <= 100000000; n *= 10){
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
        acc = 0;
        while(i < n)
            A[i++] = (float)rand();
        i = 0;
        while(i < n)
            acc += A[i++];
        i = 0;
        clock_gettime(CLOCK_MONOTONIC_RAW, &curr);
        printf("%lu\n", (curr.tv_sec - start.tv_sec) * 1000000 + (curr.tv_nsec - start.tv_nsec) / 1000);
    }
}
