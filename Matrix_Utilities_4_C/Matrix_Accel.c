#include "Matrix_Accel.h"
#include <malloc.h>




// Core Functions
struct  Matrix_t *Create_Matrix(unsigned int size_x, unsigned int size_y){
    struct Matrix_t *A = NULL;
    if(size_x && size_y){
        A = malloc(sizeof(struct Matrix_t));
        if(A){
            A->data = malloc(size_y * sizeof(MATRIX_D_TYPE *));
            if(A->data){
                unsigned int n;
                for(n = 0; n < size_y; ++n){
                    A->data[n] = calloc(sizeof(MATRIX_D_TYPE),size_x);
                    if(A->data[n]) continue;
                    else if(n){
                        // Allocation Fail
                        for(int m = 0; m < n - 1; ++m) free(A->data[m]);
                        free(A);
                        A = NULL;
                        break;
                    }
                }
                if(A){
                    A->size_x = size_x;
                    A->size_y = size_y;
                }
            } else {
                // Allocation Fail
                free(A);
                A = NULL;
            }
        }
    }
    // Return NULL if error, else everything is kosher
    return A;
}


// Destroy Matrix of any size
int Destroy_Matrix(struct Matrix_t *cool_matrix){
    int retval = MTX_FUNC_FAIL;
    
    if(cool_matrix){
        for(int n = 0; n < cool_matrix->size_y; ++n){
            free(cool_matrix->data[n]);
        }
        free(cool_matrix);
        retval = MTX_FUNC_OK;
    }

    return retval;
}

// Print Matrix
void Print_Matrix(struct Matrix_t *cool_matrix){
    if(cool_matrix){
        for(unsigned int y = 0; y < cool_matrix->size_y; ++y){
            printf("%4u -> ", y);
            for(unsigned int x = 0; x < cool_matrix->size_x; ++x){
                MATRIX_D_TYPE pval = cool_matrix->data[y][x];
                #if defined(MTX_TYPE_INT)
                    printf("%4d\t",pval);
                #elif defined(MTX_TYPE_LONG_INT)
                    printf("%4d\t",pval);
                #elif defined(MTX_TYPE_DOUBLE)
                    printf("%.3lf\t",pval);
                #else
                    printf("%.3f\t",pval);
                #endif
            }
            printf("\n");
        }
    }
}

// Fill matrix with a value
int Fill_Matrix(struct Matrix_t *cool_matrix, MATRIX_D_TYPE fill_val){
    int retval = MTX_FUNC_FAIL;
    if(cool_matrix){
        for(unsigned int y = 0; y < cool_matrix->size_y; ++y){
            for(unsigned int x = 0; x < cool_matrix->size_x; ++x){
                cool_matrix->data[y][x] = fill_val;
            }
        }
        retval = MTX_FUNC_OK;
    }
    return retval;
}


// Core Function Wrappers
// Create NxN matrix
struct Matrix_t *Create_Square_Matrix(unsigned int order){
    return Create_Matrix(order, order);
}

// Print matrix with size as a header
void Print_Matrix_w_Header(struct Matrix_t *cool_matrix){
    if(cool_matrix){
        printf("Matrix of size: %4u x %4u\n", cool_matrix->size_x, cool_matrix->size_y);
        printf("----------------------------------------------------------------------------\n");
        Print_Matrix(cool_matrix);
        printf("\n");
    }
}


// Sub-Core Functions
// Swap row pointers
int Swap_Matrix_Rows(struct Matrix_t *matrix, unsigned int row_0, unsigned int row_1){
    int retval = MTX_FUNC_FAIL;
    if(matrix){
        if(row_0 < matrix->size_y && row_1 < matrix->size_y){
            MATRIX_D_TYPE *tmp = matrix->data[row_0];
            matrix->data[row_0] = matrix->data[row_1];
            matrix->data[row_1] = tmp;
            retval = MTX_FUNC_OK;
        }
    }
    return retval;
}

// Create new matrix with requested padding
struct Matrix_t *Pad_Matrix(struct Matrix_t *A, MATRIX_D_TYPE val_2_pad, int num_left, int num_top, int num_bot, int num_right){
    struct Matrix_t *B = NULL;

    if(A){
        B = Create_Matrix(A->size_x + num_left + num_right, A->size_y + num_top + num_bot);
        if(B){
            Fill_Matrix(B, val_2_pad);
            for(unsigned int y = 0; y < A->size_y; ++y){
                for(unsigned int x = 0; x < A->size_x; ++x){
                    B->data[y + num_top][x + num_left] = A->data[y][x];
                }
            }
        }
    }

    return B;
}