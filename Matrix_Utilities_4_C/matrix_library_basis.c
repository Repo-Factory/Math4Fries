/*
    Cubic polynomial regression, or something like that

    Joseph A. De Vico
    23:24   28/07/2023
*/

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>     // Hell yeah


struct Matrix_t {
    uint16_t    dimension;
    uint32_t    presence_bitfield;
    float       **contents;
};


struct Matrix_t *Create_Matrix(uint16_t dimension);
int             Destroy_Matrix(struct Matrix_t *cool_matrix);
void            Print_Matrix(struct Matrix_t *cool_matrix);

#define DIM         4

void main(int argc, char **argv){

    struct Matrix_t *wow_a_matrix;

    if(argc > 1){
        uint16_t idk = (argv[1])[0] - '0';
        if(idk){
            wow_a_matrix = Create_Matrix(idk);
        } else {
            wow_a_matrix = Create_Matrix(DIM);
        }
    } else {
        wow_a_matrix = Create_Matrix(DIM);
    }
    

    for(uint16_t n = 0; n < wow_a_matrix->dimension; n++){
        for(uint16_t m = 0; m < wow_a_matrix->dimension; m++){
            wow_a_matrix->contents[n][m] = (float)(n * m);
        }
    }

    Print_Matrix(wow_a_matrix);
    

    Destroy_Matrix(wow_a_matrix);

}


struct Matrix_t *Create_Matrix(uint16_t dimension){
    float** cool_array_space = malloc(dimension * sizeof(float *));
    for (int n = 0; n < dimension; ++n) {
        cool_array_space[n] = malloc(dimension * sizeof(float));
    }
    
    struct Matrix_t *cool_matrix = (struct Matrix_t *)malloc(sizeof(struct Matrix_t));

    cool_matrix->dimension = dimension;
    cool_matrix->presence_bitfield = 0;
    cool_matrix->contents = cool_array_space;
    
    return cool_matrix;
}

int Destroy_Matrix(struct Matrix_t *cool_matrix){
    int retval = 0;

    for (int n = 0; n < cool_matrix->dimension; ++n) {
        free(cool_matrix->contents[n]);
    }

    free(cool_matrix->contents);
    free(cool_matrix);

    return retval;
}

void Print_Matrix(struct Matrix_t *cool_matrix){
    printf("Contents of Array, size: %2u\n", cool_matrix->dimension);

    for(uint16_t n = 0; n < cool_matrix->dimension; n++){
        printf("%2u -> ", n);
        for(uint16_t m = 0; m < cool_matrix->dimension; m++){
            printf("%.3f\t", cool_matrix->contents[n][m]);
        }
        printf("\n");
    }
}


