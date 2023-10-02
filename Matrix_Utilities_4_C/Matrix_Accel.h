/*
    Joseph A. De Vico
    sometime/somedate/2023

    Matrix Basics acceleration library for general purpose use.

    TODO:
        - Univariate/1D matrix acceleration


    Please #define one of the following in your main.c to set the datatype:
    "#define MTX_TYPE_INT"      For INT datatypes with LONG LONG INT summs
    "#define MTX_TYPE_LONG_INT" For LONG LONG INT types with LONG LONG INT sums
    "#define MTX_TYPE_DOUBLE"   For DOUBLE types with DOUBLE sums
    "#define MTX_TYPE_FLOAT"    For FLOAT types with DOUBLE sums, this is the default option with no define
*/

#ifndef MATRIX_ACCEL_h
#define MATRIX_ACCEL_h

#include <stdio.h>
#include <inttypes.h>


#if defined(MTX_TYPE_INT)

typedef int             MATRIX_D_TYPE;
typedef long long int   LONG_MATRIX_D_TYPE;

#elif defined(MTX_TYPE_LONG_INT)

typedef long long int   MATRIX_D_TYPE;
typedef long long int   LONG_MATRIX_D_TYPE;

#elif defined(MTX_TYPE_DOUBLE)

typedef double  MATRIX_D_TYPE;
typedef double  LONG_MATRIX_D_TYPE;

#else   // Default Type

typedef float   MATRIX_D_TYPE;
typedef double  LONG_MATRIX_D_TYPE;

#endif


#define DATA_SIZE       sizeof(MATRIX_D_TYPE)
#define LONG_DATA_SIZE  sizeof(LONG_MATRIX_D_TYPE)

#define MTX_FUNC_OK     0
#define MTX_FUNC_FAIL   -1

#ifndef NULL
#define NULL    0ul
#endif

#ifndef FALSE
#define FALSE   0
#endif

#ifndef TRUE
#define TRUE    1
#endif

struct Matrix_t {
    unsigned int size_x;
    unsigned int size_y;
    // Add virtual transpose flag for 1D accel
    MATRIX_D_TYPE **data;       // Row -> Col
};


// Core Functions
struct  Matrix_t    *Create_Matrix(unsigned int size_x, unsigned int size_y);
struct  Matrix_t    *Copy_Matrix(struct Matrix_t *matrix);
int                 Destroy_Matrix(struct Matrix_t *cool_matrix);
void                Print_Matrix(struct Matrix_t *cool_matrix);
int                 Fill_Matrix(struct Matrix_t *cool_matrix, MATRIX_D_TYPE fill_val);

// Core Function Wrappers
struct  Matrix_t    *Create_Square_Matrix(unsigned int order);
void                Print_Matrix_w_Header(struct Matrix_t *cool_matrix);

// Sub-Core Functions
int                 Swap_Matrix_Rows(struct Matrix_t *matrix, unsigned int row_0, unsigned int row_1);
struct Matrix_t     *Pad_Matrix(struct Matrix_t *A, MATRIX_D_TYPE val_2_pad, int num_left, int num_top, int num_bot, int num_right); // Does not destroy or alter input

// Core Verification Functions
static inline unsigned char __same_dimension_mtx(struct Matrix_t *mat_0, struct Matrix_t *mat_1){
    return (mat_0 && mat_1) ? (mat_0->size_x == mat_1->size_x && mat_0->size_y == mat_1->size_y) : 0;
}

static inline unsigned char __can_be_multiplied_mtx(struct Matrix_t *mat_0, struct Matrix_t *mat_1){
    return (mat_0 && mat_1) ? (mat_0->size_x == mat_1->size_y && mat_0->size_y == mat_1->size_x) : 0;
}

static inline unsigned char __is_square_mtx(struct Matrix_t *mat_0){
    return (mat_0) ? (mat_0->size_x == mat_0->size_y) : 0;
}

static inline unsigned char __is_same_order(struct Matrix_t *mat_0, struct Matrix_t *mat_1){
    return (mat_0 && mat_1) ? ((mat_0->size_x == mat_1->size_x) && (__is_square_mtx(mat_0)) && (__is_square_mtx(mat_1))) : 0;
}

static inline unsigned char __non_zero_diag_mtx(struct Matrix_t *A){
    unsigned int retval = 0;
    if(A && __is_square_mtx(A)){
        for(unsigned int n = 0; n < A->size_x; ++n){
            retval += (A->data[n][n]) ? 1 : 0;
        }
        retval = (retval == A->size_x) ? 1 : 0;
    }
    return (unsigned char)retval;
};

// Basic Maths
int                 Are_Matricies_Equal(struct Matrix_t *A, struct Matrix_t *B);        // A == B,          A.size == B.size
int                 Add_Matrix(struct Matrix_t *A, struct Matrix_t *B);                 // A = A + B,       A.size == B.size
struct Matrix_t     *n_Add_Matrix(struct Matrix_t *A, struct Matrix_t *B);              // C = A + B,       A.size == B.size
int                 Sub_Matrix(struct Matrix_t *A, struct Matrix_t *B);                 // A = A - B,       A.size == B.size
struct Matrix_t     *n_Sub_Matrix(struct Matrix_t *A, struct Matrix_t *B);              // C = A - B,       A.size == B.size
struct Matrix_t     *n_Transpose_Matrix(struct Matrix_t *cool_matrix);                  // C = A ^ T,       ( cat :p ), does not delete A
int                 Transpose_Matrix(struct Matrix_t *cool_matrix);                     // A = A ^ T,       t
LONG_MATRIX_D_TYPE  _Vector_Dot_Product(MATRIX_D_TYPE *v0, MATRIX_D_TYPE *v1, unsigned int len);    // SUM(V0[n] * V1[n])
struct Matrix_t     *Mul_Matrix(struct Matrix_t *A, struct Matrix_t *B);                // C = A * B,       A.x == B.y, A.y == B.x
int                 Hadamard_Product(struct Matrix_t *A, struct Matrix_t *B);           // A = A .* B,      A.size == B.size
struct Matrix_t     *n_Hadamard_Product(struct Matrix_t *A, struct Matrix_t *B);        // M2 = A .* B,     A.size == B.size
struct Matrix_t     *Mul_Large_Matrix(struct Matrix_t *A, struct Matrix_t *B);          // M2 = A .* B^T,   Ensures more cache hits accelerating large array multiplies
LONG_MATRIX_D_TYPE  Sum_of_Matrix_Elements(struct Matrix_t *A);                         // RET = SUM(A)


#endif