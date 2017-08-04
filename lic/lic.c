#include <assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <errno.h>

#define XMAX 100
#define YMAX 100
#define LINE_LEN 1024
int xsize, ysize;
float *X, *Y;
float *dx, *dy;
// pix_idx()
// get_vector()
// sl()
// partial_integral()
// kernel
// lic
// image

void coordGrid(float *X, float *Y)
{
    // Make a linspace for X and Y
    for (int i = 0; i < xsize; i++) {
        X[i] = (float)(i)*XMAX / (xsize - 1);
    }
    for (int j = 0; j < ysize; j++) {
        Y[j] = (float)(j)*YMAX / (ysize - 1);
    }
}

void load_vectors(FILE* fd)
{
   char row[LINE_LEN];
       while (fgets (row, LINE_LEN, fd)) puts(row);
   return;
}

/** Given the input coordinates, return the grid index **/
int pixIdx(float *m_coords, int *m_grid_idx)
{
    float tempx = (float)XMAX;
    float tempy = (float)YMAX;
    float diffx, diffy;
    int i, j;
    for (i = 0; i < xsize; i++) {
	diffx = fabs(X[i] - m_coords[0]); 
        // printf("i, Xi, xcoord, diffx = %d, %f, %f, %f, %f\n", i, X[i], m_coords[0], diffx, tempx);
	if (diffx < tempx) {
	    tempx = diffx;
        } else { break; }
    }
    m_grid_idx[0] = i;
    for (j = 0; j < ysize; j++) {
        diffy = fabs(Y[j] - m_coords[1]); 
        // printf("j, Yi, ycoord, diffy = %d, %f, %f, %f, %f\n", j, Y[j], m_coords[1], diffy, tempy);
        if (diffy < tempy) {
	    tempy = diffy;
        } else { break; }
    }
    m_grid_idx[1] = j;
    return 0;
}

int main (int argc, char *argv[])
{
    if(argc != 3) {
        printf("Wrong number of arguments.\n");
	exit(EXIT_FAILURE);
    }
    xsize = atoi(argv[1]);
    ysize = atoi(argv[2]);
    X = (float *)malloc(xsize*sizeof(float));
    Y = (float *)malloc(ysize*sizeof(float));
    dx = (float *)malloc(ysize*xsize*sizeof(float));
    dy = (float *)malloc(ysize*xsize*sizeof(float));
    coordGrid(X, Y);
    
    FILE* dx_fd = fopen("dx.dat", "r");
    if (!dx_fd) {
        perror("Error: "); 
	exit(EXIT_FAILURE);
    } else { load_vectors(dx_fd); }
    fclose(dx_fd);
    FILE* dy_fd = fopen("dy.dat", "r");
    if (!dy_fd) {
        perror("Error: "); 
	exit(EXIT_FAILURE);
    } else {load_vectors(dy_fd); }
    fclose(dy_fd);

    float coords[2] = {30.4, 20.3};
    int grid_idx[2];
    pixIdx(coords, grid_idx);
    printf("xidx = %d, y_idx = %d\n", grid_idx[0], grid_idx[1]);
    return 0;
}
