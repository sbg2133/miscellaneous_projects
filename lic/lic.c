#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<errno.h>

#define XMAX 100
#define YMAX 100
#define LINE_LEN 1024
#define PI 3.14159265

typedef struct {
    int xsize, ysize;
    double *X, *Y;
    double **dx, **dy;
} vector_field_t;
// sl()
// partial_integral()
// kernel
// lic
// image

void coordGrid(vector_field_t *m_field)
{
    m_field->X = (double *)malloc(m_field->xsize*sizeof(double));
    m_field->Y = (double *)malloc(m_field->ysize*sizeof(double));
    // Make a linspace for X and Y
    for (int i = 0; i < m_field->xsize; i++) {
        m_field->X[i] = (double)(i)*XMAX / (m_field->xsize - 1);
    }
    for (int j = 0; j < m_field->ysize; j++) {
        m_field->Y[j] = (double)(j)*YMAX / (m_field->ysize - 1);
    }
}

void load_vectors(vector_field_t *m_field)
{
    FILE* dx_fd = fopen("dx.dat", "r");
    if (!dx_fd) {
        perror("Error: ");
	exit(EXIT_FAILURE);
    }
    FILE* dy_fd = fopen("dy.dat", "r");
    if (!dy_fd) {
        perror("Error: ");
	exit(EXIT_FAILURE);
    }
    m_field->dx = (double**)malloc(m_field->ysize * sizeof(double*));
    m_field->dy = (double**)malloc(m_field->ysize * sizeof(double*));
    // double *vx[ysize], *vy[xsize];
    int i, j;
    for (i = 0; i < m_field->ysize; i++) {
         m_field->dx[i] = (double *)malloc(m_field->xsize * sizeof(double));
         m_field->dy[i] = (double *)malloc(m_field->xsize * sizeof(double));
    }
    for (i = 0; i < m_field->ysize; i++) {
        for (j = 0; j < m_field->xsize; j++) {
            fscanf(dx_fd, "%lf", &m_field->dx[i][j]);
            fscanf(dy_fd, "%lf", &m_field->dy[i][j]);
        }
    }
    fclose(dx_fd);
    fclose(dy_fd);
    return;
}

/** Given the input coordinates, return the grid index **/
int pixIdx(vector_field_t *m_field, double *m_coords, int *grid_idx)
{
    double tempx = (double)XMAX;
    double tempy = (double)YMAX;
    double diffx, diffy;
    int i, j;
    for (i = 0; i < m_field->xsize; i++) {
	diffx = fabs(m_field->X[i] - m_coords[0]);
        // printf("i, Xi, xcoord, diffx = %d, %f, %f, %f, %f\n",
                   // i, m_field->X[i], m_coords[0], diffx, tempx);
	if (diffx < tempx) {
	    tempx = diffx;
        } else { break; }
    }
    grid_idx[0] = i;
    for (j = 0; j < m_field->ysize; j++) {
        diffy = fabs(m_field->Y[j] - m_coords[1]);
        // printf("j, Yi, ycoord, diffy = %d, %f, %f, %f, %f\n",
                       // j, m_field->Y[j], m_coords[1], diffy, tempy);
        if (diffy < tempy) {
	    tempy = diffy;
        } else { break; }
    }
    grid_idx[1] = j;
    return 0;
}

int getVector(vector_field_t *m_field, double *m_coords, double *vector,
                   double angle)
{
    int grid_idx[2];
    pixIdx(m_field, m_coords, grid_idx);
    // printf("idx = %d, %d\n", grid_idx[0], grid_idx[1]);
    vector[0] = m_field->dx[grid_idx[1]][grid_idx[0]];
    vector[1] = m_field->dy[grid_idx[1]][grid_idx[0]];
    angle = atan2(vector[1], vector[0]) * 180. / PI;
    return 0;
}

int advectP(vector_field_t *m_field, int idx, int back)
{
    return 0;
}

int sl(vector_field_t *m_field)
{
    return 0;
}

int kernel(vector_field_t *m_field)
{
    return 0;
}

int partial_integral(vector_field_t *m_field)
{
    return 0;
}

int main (int argc, char *argv[])
{
    if(argc != 3) {
        printf("Wrong number of arguments.\n");
	exit(EXIT_FAILURE);
    }
    vector_field_t field;
    field.xsize = atoi(argv[1]);
    field.ysize = atoi(argv[2]);
    // Create a coordinate grid
    coordGrid(&field);

    // Load in the vector arrays, dx and dy
    load_vectors(&field);
    /* for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < xsize; j++) {
            printf("%f, %f\n", field.dx[i][j], field.dy[i][j]);
        }
    } */

    double coords[2] = {30.4, 20.3};
    /* int grid_idx[2];
    pixIdx(&field, coords, grid_idx);
    printf("xidx = %d, y_idx = %d\n", grid_idx[0], grid_idx[1]); */
    double vector[2];
    double angle;
    getVector(&field, coords, vector, angle);
    printf("vector = %f, %f, angle = %lf\n", vector[0], vector[1], angle);
    return 0;
}
