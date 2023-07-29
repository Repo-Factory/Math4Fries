/*
    Basic matrix operators and creation tools, will be turned into
    a library soon. The Matrix_t datatype will be used for future projects,
    only the cubic spline 0 version will remain grandfathered in.

    Joseph A. De Vico
    23:24   28/07/2023
*/

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>     // Hell yeah


struct Matrix_Row_t {
#ifdef SUPPORT_LARGE_ARRAYS    
    uint64_t    presence_bitfield;      // Represents NULL or FILLED array elements in bitfield format
#else
    uint32_t    presence_bitfield;      // Run SUPPORT_LARGE_ARRAYS if larger than 32x32 matricies are required
#endif
    float       *column;
};

struct Matrix_t {
    uint8_t     dimension_x;            // x dimension, if you need more than 256 then.... idk man
    uint8_t     dimension_y;            // same but for Y
    struct Matrix_Row_t *row;           // Row contents and bitfield
};


// Standard Utilities
struct Matrix_t         *Create_Matrix(uint8_t dimension_x, uint8_t dimension_y);
struct Matrix_t         *Create_Square_Matrix(uint8_t dimension);
int                     Destroy_Matrix(struct Matrix_t *cool_matrix);
void                    Print_Matrix(struct Matrix_t *cool_matrix);

// Static util, do not call these.
static int _are_matricies_the_same_dimension(struct Matrix_t *cool_matrix_0, struct Matrix_t *cool_matrix_1);



// Standard math operations
// Operations led by n_ indicate a new matrix
//  is created and returned, this is for operations
//  that may require the original to remain untouched
struct Matrix_t         *n_Transpose_Matrix(struct Matrix_t *cool_matrix);
void                    Transpose_Matrix(struct Matrix_t *cool_matrix);
void                    Swap_Matrix_Row(struct Matrix_t *cool_matrix, uint8_t row_0, uint8_t row_1);
struct Matrix_t         *n_Hadamard_Product(struct Matrix_t *cool_matrix_0, struct Matrix_t *cool_matrix_1);



#define DIM         4

void main(int argc, char **argv){

    struct Matrix_t *wow_a_matrix;

    if(argc == 2){
        uint16_t idk = (argv[1])[0] - '0';
        if(idk){
            wow_a_matrix = Create_Square_Matrix(idk);
        } else {
            wow_a_matrix = Create_Square_Matrix(DIM);
        }
    } else if(argc == 3){
        uint16_t idk = (argv[1])[0] - '0';
        uint16_t idk2 = (argv[2])[0] - '0';
        if(idk && idk2){
            wow_a_matrix = Create_Matrix(idk, idk2);
        } else {
            wow_a_matrix = Create_Square_Matrix(DIM);
        }
    } else {
        wow_a_matrix = Create_Matrix(DIM, DIM);
    }
    
    printf("Cool Dimensions!\nFilling With Some Numbers.\n");

    for(uint16_t n = 0; n < wow_a_matrix->dimension_y; n++){
        for(uint16_t m = 0; m < wow_a_matrix->dimension_x; m++){
            wow_a_matrix->row[n].column[m] = (float)(m + (n * wow_a_matrix->dimension_x));
        }
    }

    printf("Matrix is all filled in!\n");
    Print_Matrix(wow_a_matrix);

    printf("\nWow we can swap an array row!\n");
    Swap_Matrix_Row(wow_a_matrix, 1, 3);
    Print_Matrix(wow_a_matrix);

    printf("\nTranspose Complete!\n");
    struct Matrix_t *pog = n_Transpose_Matrix(wow_a_matrix);
    Print_Matrix(pog);

    printf("\nTranspose Back Again!\n");
    Transpose_Matrix(pog);
    Print_Matrix(pog);

    printf("\n .* Operation On Itself! (A = A .* A)!\n");
    struct Matrix_t *wowza = n_Hadamard_Product(pog, pog);
    if(wowza){
        Print_Matrix(wowza);
    }

    Destroy_Matrix(wowza);
    Destroy_Matrix(wow_a_matrix);
    Destroy_Matrix(pog);
}


struct Matrix_t *Create_Matrix(uint8_t dimension_x, uint8_t dimension_y){

    struct Matrix_t *cool_matrix = (struct Matrix_t *)malloc(sizeof(struct Matrix_t));

    cool_matrix->row = malloc(dimension_y * sizeof(struct Matrix_Row_t));
    for(int n = 0; n < dimension_y; ++n){
        cool_matrix->row[n].column = (float *)calloc(dimension_x, sizeof(float));
        cool_matrix->row[n].presence_bitfield = 0;
    }

    cool_matrix->dimension_x = dimension_x;
    cool_matrix->dimension_y = dimension_y;

    return cool_matrix;
}

struct Matrix_t *Create_Square_Matrix(uint8_t dimension){
    return Create_Matrix(dimension, dimension);
}

int Destroy_Matrix(struct Matrix_t *cool_matrix){
    int retval = 0;

    if(cool_matrix){    // Don't wanna free random offsets from NULLPTR....
        for(int n = 0; n < cool_matrix->dimension_y; ++n){
            free(cool_matrix->row[n].column);
            cool_matrix->row[n].column = NULL;
        }

        free(cool_matrix->row);
        cool_matrix->row = NULL;


        free(cool_matrix);
        cool_matrix = NULL;
    } else {
        retval = -1;
    }
    
    return retval;
}

void Print_Matrix(struct Matrix_t *cool_matrix){
    printf("Contents of Array, size: %u x %u\n", cool_matrix->dimension_x, cool_matrix->dimension_y);

    for(uint16_t n = 0; n < cool_matrix->dimension_y; n++){
        printf("%2u -> ", n);
        for(uint16_t m = 0; m < cool_matrix->dimension_x; m++){
            printf("%.3f\t", cool_matrix->row[n].column[m]);
        }
        printf("\n");
    }
}

static int _are_matricies_the_same_dimension(struct Matrix_t *cool_matrix_0, struct Matrix_t *cool_matrix_1){
    return (cool_matrix_0->dimension_x == cool_matrix_1->dimension_x) && (cool_matrix_0->dimension_y == cool_matrix_1->dimension_y);
}

//
//  Math Functions
//

/*
    Returns new matrix as many operations need both A and A^T
*/
struct Matrix_t *n_Transpose_Matrix(struct Matrix_t *cool_matrix){
    // Remember: (N x M) ^ T = (M x N)
    struct Matrix_t *cooler_matrix = Create_Matrix(cool_matrix->dimension_y, cool_matrix->dimension_x);

    for(int n = 0; n < cool_matrix->dimension_y; ++n){
        for(int m = 0; m < cool_matrix->dimension_x; ++m){
            cooler_matrix->row[m].column[n] = cool_matrix->row[n].column[m];
        }
    }

    return cooler_matrix;
}

/*
    Quick wrapper to change the original input data to the transpose.
    This is one way and will destroy the original array, though
    I suppose you can run it back thru and get the original out.
*/
void Transpose_Matrix(struct Matrix_t *cool_matrix){
    struct Matrix_t *retainer = cool_matrix;
    struct Matrix_t *tmp = n_Transpose_Matrix(cool_matrix);

    // Destroy original matrix, we'll rewrite what's important later.
    Destroy_Matrix(cool_matrix);
    
    // Basically retainer holds the pointer to our input pass by reference,
    //  and saves it from being freed/destroyed.
    //  Then we take everything we want from tmp (dimensions, pointer to content)
    //  and free tmp (which frees the row pointer, dim x, and dim y space)
    retainer->dimension_x = tmp->dimension_x;
    retainer->dimension_y = tmp->dimension_y;
    retainer->row = tmp->row;
    
    // Nuke everything, noting that we've assigned the row contents pointer
    //  elsewhere (in retainer) and that pointer will persist.
    free(tmp);
}

void Swap_Matrix_Row(struct Matrix_t *cool_matrix, uint8_t row_0, uint8_t row_1){
    struct Matrix_Row_t tmp;
    tmp = cool_matrix->row[row_0];

    cool_matrix->row[row_0] = cool_matrix->row[row_1];
    cool_matrix->row[row_1] = tmp;
}

/*
    This is an element by element multiply, similar to the .* Oxtave/Matlab operator.
    Matricies MUST be the same size, else return 0
*/
struct Matrix_t *n_Hadamard_Product(struct Matrix_t *cool_matrix_0, struct Matrix_t *cool_matrix_1){
    struct Matrix_t *retptr = NULL;

    if( _are_matricies_the_same_dimension(cool_matrix_0, cool_matrix_1) ){
        retptr = Create_Matrix(cool_matrix_0->dimension_x, cool_matrix_0->dimension_y);
        for(int n = 0; n < retptr->dimension_y; ++n){
            for(int m = 0; m < retptr->dimension_x; ++m){
                retptr->row[n].column[m] = cool_matrix_0->row[n].column[m] * cool_matrix_1->row[n].column[m];

                if(retptr->row[n].column[m]) retptr->row[n].presence_bitfield |= (1u << m);
            }
        }
    }

    return retptr;
}


