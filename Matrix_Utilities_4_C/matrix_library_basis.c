/*
    Cubic polynomial regression, or something like that

    Joseph A. De Vico
    23:24   28/07/2023
*/

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>     // Hell yeah


struct Matrix_t {
    uint8_t     dimension_x;
    uint8_t     dimension_y;
    uint32_t    presence_bitfield;
    float       **contents;
};


struct Matrix_t *Create_Matrix(uint8_t dimension_x, uint8_t dimension_y);
int             Destroy_Matrix(struct Matrix_t *cool_matrix);
void            Print_Matrix(struct Matrix_t *cool_matrix);

#define DIM         4

void main(int argc, char **argv){

    struct Matrix_t *wow_a_matrix;

    if(argc == 2){
        uint16_t idk = (argv[1])[0] - '0';
        if(idk){
            wow_a_matrix = Create_Matrix(idk, idk);
        } else {
            wow_a_matrix = Create_Matrix(DIM, DIM);
        }
    } else if(argc == 3){
        uint16_t idk = (argv[1])[0] - '0';
        uint16_t idk2 = (argv[2])[0] - '0';
        if(idk && idk2){
            wow_a_matrix = Create_Matrix(idk, idk2);
        } else {
            wow_a_matrix = Create_Matrix(DIM, DIM);
        }
    } else {
        wow_a_matrix = Create_Matrix(DIM, DIM);
    }
    

    for(uint16_t n = 0; n < wow_a_matrix->dimension_y; n++){
        for(uint16_t m = 0; m < wow_a_matrix->dimension_x; m++){
            wow_a_matrix->contents[n][m] = (float)(n * m);
        }
    }

    Print_Matrix(wow_a_matrix);
    

    Destroy_Matrix(wow_a_matrix);

}


struct Matrix_t *Create_Matrix(uint8_t dimension_x, uint8_t dimension_y){
    float** cool_array_space = malloc(dimension_y * sizeof(float *));
    for (int n = 0; n < dimension_y; ++n) {
        cool_array_space[n] = malloc(dimension_x * sizeof(float));
    }
    
    struct Matrix_t *cool_matrix = (struct Matrix_t *)malloc(sizeof(struct Matrix_t));

    cool_matrix->dimension_x = dimension_x;
    cool_matrix->dimension_y = dimension_y;
    cool_matrix->presence_bitfield = 0;
    cool_matrix->contents = cool_array_space;
    
    return cool_matrix;
}

int Destroy_Matrix(struct Matrix_t *cool_matrix){
    int retval = 0;

    for (int n = 0; n < cool_matrix->dimension_y; ++n) {
        free(cool_matrix->contents[n]);
    }

    free(cool_matrix->contents);
    free(cool_matrix);

    return retval;
}

void Print_Matrix(struct Matrix_t *cool_matrix){
    printf("Contents of Array, size: %u x %u\n", cool_matrix->dimension_x, cool_matrix->dimension_y);

    for(uint16_t n = 0; n < cool_matrix->dimension_y; n++){
        printf("%2u -> ", n);
        for(uint16_t m = 0; m < cool_matrix->dimension_x; m++){
            printf("%.3f\t", cool_matrix->contents[n][m]);
        }
        printf("\n");
    }
}


