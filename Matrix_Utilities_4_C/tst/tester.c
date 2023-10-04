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
    printf("Matrix of size: %4u x %4u requested\n", sizeX, sizeY);

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

    printf("mtx * mtx\n");
    struct Matrix_t *mtt = Mul_Large_Matrix(mtx,mtx);

    if(mtt){
        printf("mtt = \n");
        Print_Matrix_w_Header(mtt);
    } else {
        printf("Cannot Multiply!\n");
    }

    Destroy_Matrix(mtx);
    Destroy_Matrix(mtt);
    Destroy_Matrix(padd_mtx);

    return 0;
}