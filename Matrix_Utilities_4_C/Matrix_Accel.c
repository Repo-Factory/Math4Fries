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
        free(cool_matrix->data);
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

// Create X by Y matrix of 1s
struct Matrix_t *Create_Ones_Matrix(unsigned int size_x, unsigned int size_y){
    struct Matrix_t *A = Create_Matrix(size_x, size_y);
    Fill_Matrix(A,1);
    return A;
}

// Create X by Y matrix of 0s
struct Matrix_t *Create_Zeros_Matrix(unsigned int size_x, unsigned int size_y){
    struct Matrix_t *A = Create_Matrix(size_x, size_y);
    Fill_Matrix(A,0);
    return A;
}

// Create NxN ID matrix
struct Matrix_t *Create_ID_Matrix(unsigned int order){
    struct Matrix_t *A = Create_Zeros_Matrix(order, order);
    if(A){
        for(unsigned int n = 0; n < order; ++n){
            A->data[n][n] = 1;
        }
    }
    return A;
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


// Basic Maths

// A == B,          A.size == B.size
int Are_Matricies_Equal(struct Matrix_t *A, struct Matrix_t *B){
    int retval = FALSE;

    if(__same_dimension_mtx(A,B)){
        retval = TRUE;
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < A->size_x; ++x){
                if(A->data[y][x] != B->data[y][x]){
                    retval = FALSE;
                    goto exit_equality_check;
                }
            }
        }
    }

    exit_equality_check:

    return retval;
}


// A = A + B,       A.size == B.size
int Add_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    int retval = FALSE;

    if(__same_dimension_mtx(A,B)){
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < A->size_x; ++x){
                A->data[y][x] += B->data[y][x];
            }
        }
        
        retval = TRUE;
    }

    return retval;
}  


// C = A + B,       A.size == B.size
struct Matrix_t *n_Add_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    struct Matrix_t *C = NULL;
    
    if(__same_dimension_mtx(A,B)){
        C = Create_Matrix(A->size_x,A->size_y);
        if(C){
            for(unsigned int y = 0; y < C->size_y; ++y){
                for(unsigned int x = 0; x < C->size_x; ++x){
                    C->data[y][x] = A->data[y][x] + B->data[y][x];
                }
            }
        }
    }

    return C;
}


// A = A - B,       A.size == B.size
int Sub_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    int retval = FALSE;

    if(__same_dimension_mtx(A,B)){
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < A->size_x; ++x){
                A->data[y][x] -= B->data[y][x];
            }
        }
        
        retval = TRUE;
    }

    return retval;
}


// C = A - B,       A.size == B.size
struct Matrix_t *n_Sub_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    struct Matrix_t *C = NULL;
    
    if(__same_dimension_mtx(A,B)){
        C = Create_Matrix(A->size_x,A->size_y);
        if(C){
            for(unsigned int y = 0; y < C->size_y; ++y){
                for(unsigned int x = 0; x < C->size_x; ++x){
                    C->data[y][x] = A->data[y][x] - B->data[y][x];
                }
            }
        }
    }

    return C;
}


// C = A ^ T,       ( cat :p ), does not delete A
struct Matrix_t *n_Transpose_Matrix(struct Matrix_t *cool_matrix){
// Remember: (N x M) ^ T = (M x N)
    struct Matrix_t *C = NULL;
    if(cool_matrix){
        C = Create_Matrix(cool_matrix->size_y,cool_matrix->size_x);
        for(unsigned int y = 0; y < cool_matrix->size_y; ++y){
            for(unsigned int x = 0; x < cool_matrix->size_x; ++x){
                C->data[x][y] = cool_matrix->data[y][x];
            }
        }
    }
    return C;
}


// A = A ^ T, 
int Transpose_Matrix(struct Matrix_t *cool_matrix){
    struct Matrix_t *tee = n_Transpose_Matrix(cool_matrix);
    
    Destroy_Matrix(cool_matrix);

    cool_matrix = tee;

    return TRUE;
}

// SUM(V0[n] * V1[n])
LONG_MATRIX_D_TYPE  _Vector_Dot_Product(MATRIX_D_TYPE *v0, MATRIX_D_TYPE *v1, unsigned int len){
    LONG_MATRIX_D_TYPE sum_ = 0;
    for(unsigned int e = 0; e < len; e++){
        sum_ += v0[e] * v1[e];
    }
    return sum_;
}

// C = A * B,       A.x == B.y, A.y == B.x
struct Matrix_t *Mul_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    struct Matrix_t *C = NULL;
    if(__can_be_multiplied_mtx(A,B)){
        C = Create_Matrix(A->size_y, B->size_x);
        
        // This is speed method... because it's easier
        struct Matrix_t *D = n_Transpose_Matrix(B);

        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < D->size_y; ++x){
                C->data[y][x] = _Vector_Dot_Product(A->data[y],D->data[x], A->size_x);
            }
        }

        Destroy_Matrix(D);
    }
    return C;
}


// A = A .* B,      A.size == B.size
int Hadamard_Product(struct Matrix_t *A, struct Matrix_t *B){
    int retval = MTX_FUNC_FAIL;
    if(__same_dimension_mtx(A,B)){
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < A->size_x; ++x){
                A->data[y][x] *= B->data[y][x];
            }
        }
        retval = MTX_FUNC_OK;
    }
    return retval;
}


// C = A .* B^T,   Ensures more cache hits accelerating large array multiplies
struct Matrix_t *Mul_Large_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    struct Matrix_t *C = NULL;
    if(__can_be_multiplied_mtx(A,B)){
        C = Create_Matrix(A->size_y, B->size_x);
        
        // This is speed method... because it's easier
        struct Matrix_t *D = n_Transpose_Matrix(B);

        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < D->size_y; ++x){
                C->data[y][x] = _Vector_Dot_Product(A->data[y],D->data[x], A->size_x);
            }
        }

        Destroy_Matrix(D);
    }
    return C;
}


// C = A .* B^T,   Arg B comes in trasnposed. 
struct Matrix_t *Mul_TPSD_Matrix(struct Matrix_t *A, struct Matrix_t *B){
    struct Matrix_t *C = NULL;
    if(__same_dimension_mtx(A,B)){
        C = Create_Matrix(A->size_y,B->size_x);
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < B->size_y; ++x){
                C->data[y][x] = _Vector_Dot_Product(A->data[y],B->data[x], A->size_x);
            }
        }
    }
    return C;
}


// RET = SUM(A)
LONG_MATRIX_D_TYPE  Sum_of_Matrix_Elements(struct Matrix_t *A){
    LONG_MATRIX_D_TYPE sum_ = 0;
    if(A){
        for(unsigned int y = 0; y < A->size_y; ++y){
            for(unsigned int x = 0; x < A->size_x; ++x){
                sum_ += A->data[y][x];
            }
        }
    }

    return sum_;
}


// Ordering Functions
// Orders A such that diagonals are non-zero, if B is provided perform same row operations on B
int MTX_Order_NonZero_Diagonals(struct Matrix_t *A, struct Matrix_t *B){
    int retval = MTX_FUNC_FAIL;
    if(A){
        if(__is_square_mtx(A)){
                for(unsigned int n = 0; n < A->size_y; ++n){
                    if(!A->data[n][n]){
                        unsigned int m;
                        for(m = 0; m < A->size_y; ++m){
                            if(A->data[m][n] && A->data[n][m]){
                                Swap_Matrix_Rows(A,m,n);
                            }
                        }
                        if(m == A->size_y){
                            retval = MTX_FUNC_FAIL;
                            goto mtx_not_sortable;
                        } 
                    }
                }
                retval = MTX_FUNC_OK;
            }
        }

    mtx_not_sortable:

    return retval;
}