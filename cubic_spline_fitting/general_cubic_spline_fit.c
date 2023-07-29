/**********************
 *
 * Alec Stobbs
 * Discord: @PotatoeComrade
 * 07/18/2023
 *
 * The purpose of this program
 * is to produce quadratic spline
 * interpolated values to enhance
 * FIR filter speed and accuracy
 * of CV inputs.
 *
***********************/
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <malloc.h>         // oh yeahhhhh

//#define DEBUG
//#define DEBUG_ARRAY_SORT
#define REAL_DATA
//#define OUTPUT_PRINT
//#define ALEC_DIAG         // Comment out for joseph solver to be used
#define DELAY_SAMPLES   550

#define SRES    50      // resolution of spline interpolation
#define DEPTH   9       // number of origin points to interpolate accross

struct RowParams {
    uint32_t    rbfield;    // bitfield
    uint16_t    posvals;    // sum of row position vals for coarse sort
    float       movability; // How movable is this row? (0 -> only one position 1-> any)
};

struct Cartesian {
    float x[DEPTH - 1][SRES];
    float y[DEPTH - 1][SRES];
};

struct  Cartesian* qSpline(float[], float[]);
float*  Sd0(float);
float*  Sd1(float);
float*  Sd2(float);
float*  solveLU(float A[][4 * (DEPTH - 1)], float[]);
void    arrayPrint(float[], uint16_t, uint16_t);
void    array2Print(float A[][4 * (DEPTH - 1)], uint16_t, uint16_t);
void    diagonal(float[], float A[][4*(DEPTH-1)], float[], int);
void    other_diagonal(float a[][4 * (DEPTH - 1)], float *b, int N);
void    reverse_bin_print(struct RowParams arr[], int N);   // debug presence bitfield
void    calc_presence(float a[][4 * (DEPTH - 1)], struct RowParams arr[], int N);
uint8_t find_rows_2_correct(struct RowParams arr[], uint8_t *row2corr, int N);
uint8_t recursive_solver_2(float a[][4 * (DEPTH - 1)], float *b, struct RowParams ROWZ[], uint8_t *rows_2_correct, int N);
uint8_t check_diag_from_presence(struct RowParams arr[], int N);
void    print_diag_contents(float a[][4 * (DEPTH - 1)], struct RowParams arr[], int N);
void    swapRow(float a[][4 * (DEPTH - 1)], float[], uint16_t, uint16_t);

void main() {
    
#ifdef REAL_DATA
    FILE *dataset_fp;
    dataset_fp = fopen("joseph_output.txt", "r");
    float x[DEPTH];
    float y[DEPTH];

    for(int n = 0; n < DELAY_SAMPLES; n++){
        fscanf(dataset_fp, "%f, %f\n", &x[0], &y[0]);
    }

    for(int n = 0; n < DEPTH; n++){
        fscanf(dataset_fp, "%f, %f\n", &x[n], &y[n]);
        x[n] = n * SRES;
        //printf("%f\t%f\n", x[n], y[n]);
    }

#ifdef DEBUG
    printf("X:\n");
    arrayPrint(x, 0, DEPTH);
    printf("Y:\n");
    arrayPrint(y, 0, DEPTH);
#endif
#else
    float x[DEPTH] = { 0,1,2,3,4 };
    float y[DEPTH] = { 1,3,2,4,1 };
#endif
    
    if(4 * (DEPTH - 1) > 32 + 1){
        printf("Depth too large! Cannot create bitfield.\n");
        return;
    }
    
    
    FILE *fp;
    FILE *fd;
    fp = fopen("pogout.txt", "w+");
    fd = fopen("pogin.txt", "w+");

    struct Cartesian* output;
    output = qSpline(x, y);

#ifdef OUTPUT_PRINT
    printf("Output Points Set:\n");
#endif
    for(int m = 0; m < DEPTH - 1; m++){
        for(int n = 0; n < SRES; n++){
#ifdef OUTPUT_PRINT
            printf("Point [%2u] = {%.2f, %.2f}\n", n, output->x[m][n], output->y[m][n]);
#endif
            fprintf(fp, "%.2f, %.2f\n", output->x[m][n], output->y[m][n]);
        }
        fprintf(fd, "%.2f, %.2f\n", x[m], y[m]);
    }
    fprintf(fd, "%.2f, %.2f\n", x[DEPTH-1], y[DEPTH-1]);

#ifdef DEBUG
    printf("From Dataset:\n");
    printf("X:\n"); 
    arrayPrint(x, 0, DEPTH);
    printf("Y:\n");
    arrayPrint(y, 0, DEPTH);
#endif

#ifdef OUTPUT_PRINT
    printf("x Spline:\n");
    for(uint16_t i = 0; i < DEPTH - 1; ++i){
        printf("[");
        for(uint16_t j = 0; j < SRES; ++j){
            printf("%5.2f, ", output->x[i][j]);
        }
        printf("]\n");
    }
    printf("y Spline:\n");
        for(uint16_t i = 0; i < DEPTH - 1; ++i){
        printf("[");
        for(uint16_t j = 0; j < SRES; ++j){
            printf("%5.2f, ", output->y[i][j]);
        }
        printf("]\n");
    }
#endif

    fclose(fp);
    fclose(fd);
}

/**\brief Quadratic Spline Interpolater
 * \ingroup Equation
 *
 * \param x array of x-axis coordinates
 * \param y array of y-axis coordinates
 *
 * \return Pointer to an array of interpolated values and their respective x coordinates
 * wrapped in a Cartesian Struct
*/
struct Cartesian* qSpline(float x[], float y[]) {

    static const uint16_t n = DEPTH - 1;
    float A[4 * (DEPTH - 1)][4 * (DEPTH - 1)] = { {0}, {0} };
    float b[4 * (DEPTH - 1)] = { 0 };
    uint16_t j = 0;

    // 1.
    for (uint16_t i = 0; i < n; ++i) {
        uint8_t k = 0;
        for (uint16_t h = i * 4; h < (i + 1) * 4; ++h) {
            A[j][h] = Sd0(x[i])[k];
            ++k;
        }
        b[j] = y[i];
        arrayPrint(A[j], i * 4, (i + 1) * 4);
        ++j;
    }

    // 2.
    for (uint16_t i = 0; i < n; ++i) {
        uint8_t k = 0;
        for (uint16_t h = i * 4; h < (i + 1) * 4; ++h) {
            A[j][h] = Sd0(x[i + 1])[k];
            ++k;
        }
        b[j] = y[i + 1];
        arrayPrint(A[j], i * 4, (i + 1) * 4);
        ++j;
    }

    // 3.
    for (uint16_t i = 0; i < n - 1; ++i) {
        uint8_t k = 0;
        for (uint16_t h = i * 4; h < (i + 2) * 4; ++h) {
            A[j][h] = Sd1(x[i + 1])[k];
            ++k;
        }
        b[j] = 0;
        arrayPrint(A[j], i * 4, (i + 2) * 4);
        ++j;
    }

    // 4.
    for (uint16_t i = 0; i < n - 1; ++i) {
        uint8_t k = 0;
        for (uint16_t h = i * 4; h < (i + 2) * 4; ++h) {
            A[j][h] = Sd2(x[i + 1])[k];
            ++k;
        }
        b[j] = 0;
        arrayPrint(A[j], i * 4, (i + 2) * 4);
        ++j;
    }

    // 5.
    //  5.1
    for (uint16_t h = 0; h < 4; ++h) {
        A[j][h] = Sd2(x[0])[h];
    }
    b[j] = 0;
    arrayPrint(A[j], 0, 4);
    ++j;

    //  5.2
    { // Scope bound for temp k variable
        uint8_t k = 0;
        for (uint16_t h = (n - 1) * 4; h < n * 4; ++h) {
            A[j][h] = Sd2(x[n])[k];
            ++k;
        }
    }
    b[j] = 0;

#ifdef DEBUG_ARRAY_SORT
    arrayPrint(A[j], (n - 1) * 4, n * 4);

    printf("b:\n");
    arrayPrint(b, 0, n * 4);

    printf("A:\n");
    for (uint16_t i = 0; i < 4 * n; ++i) {
        arrayPrint(A[i], 0, 4 * n);
    }
#endif
    float* xx = solveLU(A, b);
    //printf("xx:\n");
    //arrayPrint(xx,0,4*n);
    static struct Cartesian sx = { {{0},{0}},{{0},{0}} };

/*
S_ = np.zeros((n, 50))
    for i in range(n):
        S_[i] = xx[i*4] * x_[i] ** 3 + xx[i*4+1] * x_[i] ** 2 + xx[i*4+2] * x_[i] + xx[i*4+3]
*/

    for (uint16_t i = 0; i < n; ++i) {
        float step = (x[i + 1] - x[i]) / SRES;
        for (uint16_t h = 0; h < SRES; ++h) {
            sx.x[i][h] = x[i] + step * h;

            float val, s1, s2, s3 = 0;
            val = sx.x[i][h];
            s1 = val;
            s2 = val*val;
            s3 = val*val*val;
            sx.y[i][h] =    xx[i*4]*s3 + \
                        xx[(i*4)+1]*s2 + \
                        xx[(i*4)+2]*s1 + \
                        xx[(i*4)+3];
        }
    }

    return &sx;
}

/**\brief Fundamental Quadratic Equation S(x)
 * \ingroup Quads
 *
 * a * x^3 + b * x^2 + c_1 * x + d = y
 *
 * \param
 *
 * \return
*/
float* Sd0(float k) {
    static float arr[4];
    arr[0] = k * k * k;
    arr[1] = k * k;
    arr[2] = k;
    arr[3] = 1;

    return arr;
}

/**\brief First Derivative Quadratic Equation S'(x)
 * \ingroup Quads
 *
 *
 *
 * \param
 *
 * \return
*/
float* Sd1(float k) {
    static float arr[8];
    arr[0] = 3 * (k * k);
    arr[1] = 2 * k;
    arr[2] = 1;
    arr[3] = 0;
    arr[4] = -3 * (k * k);
    arr[5] = -2 * k;
    arr[6] = -1;
    arr[7] = 0;

    return arr;
}

/**\brief Second Derivative Quadratic Equation S''(x)
 * \ingroup Quads
 *
 *
 *
 * \param
 *
 * \return
*/
float* Sd2(float k) {
    static float arr[8];
    arr[0] = 6 * k;
    arr[1] = 2;
    arr[2] = 0;
    arr[3] = 0;
    arr[4] = -6 * k;
    arr[5] = -2;
    arr[6] = 0;
    arr[7] = 0;

    return arr;
}

float* solveLU(float A[][4 * (DEPTH - 1)], float b[]) {
    int n = 4 * (DEPTH - 1);
    float l[4 * (DEPTH - 1) + 1][4 * (DEPTH - 1) + 1] = { 0 };
    float u[4 * (DEPTH - 1) + 1][4 * (DEPTH - 1) + 1] = { 0 };
    float AA[4 * (DEPTH - 1) + 1][4 * (DEPTH - 1) + 1] = { 0 };
    float bb[4 * (DEPTH - 1) + 1] = { 0 };
    
    float sum = 0.0;
    static float x[4 * (DEPTH - 1)] = { 0 };

    // copy A and b into a padded matrix
    for (uint16_t i = 1; i < 4 * (DEPTH - 1) + 1; ++i) {
        bb[i] = b[i - 1];
        for (uint16_t j = 1; j < 4 * (DEPTH - 1) + 1; ++j) {
            AA[i][j] = A[i - 1][j - 1];
        }
    }

    //********** LU decomposition *****//
#ifdef ALEC_DIAG
    diagonal(x, A, b, n);
    //diagonal(x, A, b, n);
#else
    other_diagonal(A, b, n);
#ifdef DEBUG_ARRAY_SORT
    printf("WE ARE NOW OUT OF THE DIAGONAL FUNCTION!!!\n");
    printf("a:\n");
    array2Print(A, 0, n);
    printf("b:\n");
    arrayPrint(b, 0, n);
#endif
#endif
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++){
            //printf("DEBUG: a[%2u][%2u] = %f\n", k, k, A[k][k]);
            if (A[k][k] == 0) {
                printf("\nSolution is not exist.\n");
                printf("Broken A:\n");
                array2Print(A, 0, n);
                printf("Broken b:\n");
                arrayPrint(b, 0, n);
                return x;
            }
            float M = A[i][k] / A[k][k];
            //printf("M = %f / %f = %f\n", A[i][k] , A[k][k], M);
            for (int j = k; j < n; j++) {
                //float dbgf = A[i][j];
                A[i][j] -= M * A[k][j]; 
                //printf("A[%2u][%2u] = %f - (%f * %f) = %f - %f = %f\n", i, j, dbgf, M, A[k][j], dbgf, M * A[k][j], A[i][j]);
            }
            
            b[i] -= M * b[k];
        }
        //printf("a:\n");
        //array2Print(A, 0, n);
    }

#ifdef DEBUG_ARRAY_SORT
    printf("Final A:\n");
    array2Print(A, 0, n);
#endif

    for (int i = n - 1; i >= 0; i--) {
        float s = 0;
        for (int j = i; j < n; j++) {
            s = s + A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }

#ifdef DEBUG_ARRAY_SORT
    printf("x:\n");
    arrayPrint(x, 0, n);
#endif
    return x;
}

void arrayPrint(float arr[], uint16_t start, uint16_t stop) {
#ifdef DEBUG
    printf("[");
    for (uint16_t i = start; i < stop - 1; ++i) {
        printf("%5.f,", arr[i]);
    }
    printf("%5.f]\n", arr[stop - 1]);
#endif
}

void array2Print(float arr[][4 * (DEPTH - 1)], uint16_t start, uint16_t stop) {
#ifdef DEBUG
    for (uint16_t i = start; i < stop; ++i) {
        arrayPrint(arr[i], start, stop);
    }
#endif
}


void diagonal(float x[], float a[][4 * (DEPTH - 1)], float b[], int N) {
    // Rearrange rows to make diagonal non-zero
    float temp = 0;
    for (int i = 0; i < N; i++) {
        // Zero on the diag
        if (a[i][i] == 0) {
            int pivot = -1;
            // Go down the column
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                temp = 0.0;
                // Find non-zero value in this column and
                // other diag value that will be swapped must be non-zero
                if (a[j][i] != 0 && a[i][j] != 0){
                    //find largest valid pivot point
                    if (fabs(temp) < fabs(a[j][i])) {
                        temp = a[j][i];
                        pivot = j;
                    }
                }
            }
            // Special Case where simple pivot is not possible
            if (pivot == -1) {
                // if no pivot for row0:
                // 1. find a row1 that fits row0, the other diag value can be zero
                // 2. find a row2 that is interchangable with row1
                // 3. swap row1 and row2
                // 4. swap row1 and row0

                for (int h = 0; h < N; ++h) {
                    // if valid in column, find new row to match
                    if (i == h) continue;
                    if (a[h][i] != 0) {
                        for (int k = 0; k < N; ++k) {
                            // if rows are interchangable, swap
                            if (a[h][k] != 0 && a[k][h] != 0 && h != k) {
                                swapRow(a, b, h, k);
                                swapRow(a, b, k, i);
                                break;
                            }
                        }
                    }
                }
                continue;
            }
            swapRow(a, b, pivot, i);
        }
    }
    printf("A:\n");
    array2Print(a, 0, N);
}

void other_diagonal(float a[][4 * (DEPTH - 1)], float *b, int N){
    //uint16_t presence_bf[N];
    //calc_presence(a, presence_bf, N);   // develop initial presence bitfield

    struct RowParams ROWZ[N];

    calc_presence(a, ROWZ, N);

#ifdef DEBUG_ARRAY_SORT    
    printf("Array Size: %2u\n", N);
    reverse_bin_print(ROWZ, N);  // print presence bitfield
    printf("\n");
#endif

    
    uint16_t diag_ctr = 0;
    // Sort ascending, prob not needed
    while(!check_diag_from_presence(ROWZ, N) && !(diag_ctr > 2 * N)){
        for(int n = 1; n < N; n++){
            if((ROWZ[n].rbfield < ROWZ[n-1].rbfield)){
                swapRow(a, b, n, n - 1);
                calc_presence(a, ROWZ, N);
#ifdef DEBUG_ARRAY_SORT
                printf("-> Swapped Rows: %2u and %2u\n", n, n-1);
                reverse_bin_print(ROWZ, N);  // print presence bitfield
#endif
            }
        }
        diag_ctr += 1;
    }
  

#define SOLVER_3
    uint8_t lowest_movability_fit = 0;
#ifdef SOLVER_1
    for(int n = 0; n < N; n++){
        printf("%2u\n", n);
        //if(n > lowest_movability_fit) lowest_movability_fit = n;
        for(int m = n; m < N; m++){
            if((n != lowest_movability_fit) && (ROWZ[m].rbfield & (1u << n))){
                if(!(ROWZ[lowest_movability_fit].rbfield & (1u << n))){
                    swapRow(a, b, lowest_movability_fit, m);
                    printf("--> Swapped Rows: %2u and %2u\n", lowest_movability_fit, m);
                    lowest_movability_fit = m;
                    calc_presence(a, ROWZ, N);
                    break;
                } else
                if(ROWZ[m].movability <= ROWZ[lowest_movability_fit].movability){
                    swapRow(a, b, lowest_movability_fit, m);
                    printf("---> Swapped Rows: %2u and %2u\n", lowest_movability_fit, m);
                    lowest_movability_fit = m;
                    calc_presence(a, ROWZ, N);
                    break;
                }
            }
            if(check_diag_from_presence(ROWZ, N))break;
        }
        if(check_diag_from_presence(ROWZ, N))break;
    }
#endif
uint8_t rows_2_correct[N];
#ifdef SOLVER_2
    
    uint8_t num_rows = find_rows_2_correct(ROWZ, rows_2_correct, N);
    
    
    for(uint_fast16_t n = 0; n < num_rows ; n++){
        printf("-> Row we're correcting: %2u @ Iter: %2u\n", rows_2_correct[n], n);
        uint_fast16_t m = 0;
        for(; m < N; m++){
            if((ROWZ[rows_2_correct[n]].rbfield & (1lu << m)) && (ROWZ[m].rbfield & (1lu << rows_2_correct[n]))){
                swapRow(a, b, rows_2_correct[n], m);
                printf("--> Swapped Rows: %2u and %2u\n", rows_2_correct[n], m);
                calc_presence(a, ROWZ, N);
                break;
            }
        }
        //if(m == N){     // Failed to find a simple swap match, must do alt strategy
        //    if((ROWZ[rows_2_correct[n]].rbfield & (1lu << m)) && (ROWZ[m].rbfield & (1lu << rows_2_correct[n]))){
        //}
    }
#endif
#ifdef SOLVER_3
    recursive_solver_2(a, b, ROWZ, rows_2_correct, N);

#endif

#ifdef DEBUG_ARRAY_SORT
    printf("\nResulting Presence Array:\n");
    calc_presence(a, ROWZ, N);
    reverse_bin_print(ROWZ, N);  // print presence bitfield
    printf("\n");

    printf("Diagonal Contents:\n");
    print_diag_contents(a, ROWZ, N);

    printf("\nA:\n");
    array2Print(a, 0, N);
    printf("b:\n");
    arrayPrint(b, 0, N);
#endif
#ifdef DEBUG_ARRAY_SORT
    printf("%s\n", (check_diag_from_presence(ROWZ, N)) ? "Ladies and Gentlemen, we got him" : "f u c k o f f l i n e a r a l g (we failed)");
#endif
}

uint8_t find_rows_2_correct(struct RowParams arr[], uint8_t *row2corr, int N){
    uint8_t numrowscounted = 0;
    for(uint_fast16_t n = 0; n < N; n++){
        if(!(arr[n].rbfield & (1lu << n))) {
            row2corr[numrowscounted++] = n;
            //printf("Need to correct Row %2u @ iteration %2u\n", row2corr[numrowscounted - 1], numrowscounted - 1);
        }
    }
    return numrowscounted;
}

uint8_t recursive_solver_2(float a[][4 * (DEPTH - 1)], float *b, struct RowParams ROWZ[], uint8_t *rows_2_correct, int N){
    uint8_t retval = 0;
    uint8_t num_rows = find_rows_2_correct(ROWZ, rows_2_correct, N);
    
    for(uint_fast16_t n = 0; n < num_rows ; n++){
        //num_rows = find_rows_2_correct(ROWZ, rows_2_correct, N);
#ifdef DEBUG_ARRAY_SORT
        printf("-> Row we're correcting: %2u @ Iter: %2u\n", rows_2_correct[n], n);
#endif
        uint_fast16_t m = 0;
        for(; m < N; m++){
            if((ROWZ[rows_2_correct[n]].rbfield & (1lu << m)) && (ROWZ[m].rbfield & (1lu << rows_2_correct[n]))){
                swapRow(a, b, rows_2_correct[n], m);
#ifdef DEBUG_ARRAY_SORT
                printf("--> Swapped Rows: %2u and %2u\n", rows_2_correct[n], m);
#endif                
                calc_presence(a, ROWZ, N);
                break;
            }
        }
        if(m == N){     // Failed to find a simple swap match, must do alt strategy
            swapRow(a, b, rows_2_correct[n], rows_2_correct[n] + 1);
            calc_presence(a, ROWZ, N);
#ifdef DEBUG_ARRAY_SORT
            printf("--> Swapped Rows: %2u and %2u\tFAILOVER\n", rows_2_correct[n], rows_2_correct[n] + 1);
#endif
            retval = recursive_solver_2(a, b, ROWZ, rows_2_correct, N);
        }
    }
#ifdef DEBUG_ARRAY_SORT
    printf("Iteration %3u, solving for %2u rows.\n", retval, num_rows);
#endif
    return retval++;
}


void reverse_bin_print(struct RowParams arr[], int N){
    printf("Bit->");
    for(uint_fast16_t n = 0; n < N; n++) printf("%2u ", n);
    printf("\n");
    for(uint16_t n = 0; n < N; n++){
        printf("%2u -> ", n);
        for(uint16_t m = 0; m < N; m++){
            printf("%1u, ", (arr[n].rbfield & (1 << (m))) ? 1 : 0);
        }
        printf("\t0x%4X\t0x%4X\t%f\n", arr[n].rbfield, arr[n].posvals, arr[n].movability);
    }
}

void print_diag_contents(float a[][4 * (DEPTH - 1)], struct RowParams arr[], int N){
    for(int n = 0; n < N; n++){
        printf("Row %2u -> %sZero! -> a[%2u][%2u] = %f\n", n, (arr[n].rbfield & (1lu << n)) ? "Non-" : "", n, n, a[n][n]);
    }
}

void calc_presence(float a[][4 * (DEPTH - 1)], struct RowParams arr[], int N){
    for(int n = 0; n < N; n++){
        arr[n].rbfield = 0;
        arr[n].movability = 0;
        arr[n].posvals = 0;
        for(int m = 0; m < N; m++){
            if(a[n][m]){
                arr[n].rbfield |= (1u << m);
                arr[n].movability += 1.0 / (float)N ;
                arr[n].posvals += m;
            }
        }
    }
}

uint8_t check_diag_from_presence(struct RowParams arr[], int N){
    uint8_t retval = 0;
    for(int n = 0; n < N; n++){
        if(arr[n].rbfield & (1u << n)) retval += 1;
    }

    if(retval == (N)) printf("Diagonal!!\n");
    return (retval == (N)) ? 1 : 0;
}

void swapRow(float a[][4 * (DEPTH - 1)], float b[], uint16_t r1, uint16_t r2) {
    float temp = 0;
    for (int k = 0; k < 4*(DEPTH-1); k++) {
        temp = a[r1][k];
        a[r1][k] = a[r2][k];
        a[r2][k] = temp;
    }
    temp = b[r1];
    b[r1] = b[r2];
    b[r2] = temp;
}