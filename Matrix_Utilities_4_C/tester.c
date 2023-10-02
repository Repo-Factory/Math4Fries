#include <stdio.h>
#include <stdlib.h>

#include "Matrix_Accel.h"


int main(int argc, char **argv){
    printf("Matrix Accelerator Test\n");

    unsigned int sizeX = 5;
    unsigned int sizeY = 5;

    if(argc > 2){
        sizeX = strtoul(argv[1],NULL,10);
        sizeY = strtoul(argv[2],NULL,10);
    }
    printf("Matrix of size %4u x %4u requested\n", sizeX, sizeY);

    struct Matrix_t *mtx;

    mtx = Create_Matrix(sizeX, sizeY);

    for(int y = 0; y < mtx->size_y; ++y){
        for(int x = 0; x < mtx->size_x; ++x){
            mtx->data[y][x] = y * x;
        }
    }

    
    Print_Matrix_w_Header(mtx);

    printf("Padding with 1s on left\n");

    struct Matrix_t *padd_mtx = Pad_Matrix(mtx,1,1,0,0,0);

    Print_Matrix_w_Header(padd_mtx);

    Destroy_Matrix(mtx);
    Destroy_Matrix(padd_mtx);

    struct Matrix_t *M0;
    struct Matrix_t *M1;

    M0 = Create_Matrix(3,2);
    M1 = Create_Matrix(2,3);

    for(int y = 0; y < M0->size_y; ++y){
        for(int x = 0; x < M0->size_x; ++x){
            M0->data[y][x] = (y + 2) * (x + 1);
            M1->data[x][y] = (y + 3) * (x + 2);
        }
    }

    printf("M0\n");
    Print_Matrix_w_Header(M0);

    printf("M1\n");
    Print_Matrix_w_Header(M1);

    struct Matrix_t *M2;
    M2 = n_Transpose_Matrix(M1);

    printf("M2 is size %u by %u\n", M2->size_x, M2->size_y);
    Print_Matrix_w_Header(M2);

    struct Matrix_t *M3 = Mul_Matrix(M0,M1);
    printf("M3\n");
    Print_Matrix_w_Header(M3);

    Destroy_Matrix(M0);
    Destroy_Matrix(M1);
    Destroy_Matrix(M2);
    Destroy_Matrix(M3);

    return 0;
}