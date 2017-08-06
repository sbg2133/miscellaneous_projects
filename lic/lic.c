#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<errno.h>

#define XMAX 100 // Number of x coords in image
#define YMAX 100 // Number of y coords in image
#define LINE_LEN 1024 // Buffer length
#define KLEN 31 // Kernel length
#define PSTEPS 25 // Num streamline steps from start, for forward or back
#define PI 3.14159265

typedef struct {
    int xsize, ysize;
    double *X, *Y;
    double **dx, **dy;
    double **texture;
} vector_field_t;

typedef struct {
    double start[2];
    double **forward;
    double **back;
    double *forward_seg;
    double *back_seg;
    double *segs;
    double **streamline;
} streamline_t;

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

void load_texture(vector_field_t *m_field)
{
    FILE* tex_fd = fopen("texture.dat", "r");
    if (!tex_fd) {
        perror("Error: ");
	exit(EXIT_FAILURE);
    }
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

int advectP(vector_field_t *m_field, streamline_t *sl, int m_idx, int back)
{
    double** line;
    double *segs;
    double s;
    double vector[2];
    double angle;
    if (m_idx == 1) {
        getVector(m_field, sl->start, vector, angle);
    } else {
        getVector(m_field, sl->start, vector, angle);
        // if (isnan(vector[0]) || isnan(vector[1])) break;
    }
    double dx = vector[0];
    double dy = vector[1];
    if (back) {
        line = sl->back;
        segs = sl->back_seg;
        dx *= -1;
        dy *= -1;
    } else {
        line = sl->forward;
        segs = sl->forward_seg;
    }
    segs[0] = 0.0;
    int Px = (int)line[m_idx - 1][0];
    int Py = (int)line[m_idx - 1][1];
    double mag = sqrt(dx*dx + dy*dy);
    // arc lengths
    double s_top = ((Py + 1) - line[m_idx - 1][1])*(mag/dy);
    double s_bot = (Py - line[m_idx - 1][0])*(mag/dy);
    double s_right = ((Px + 1) - line[m_idx - 1][0])*(mag/dx);
    double s_left = (Px - line[m_idx - 1][1])*(mag/dx);
    // check for NaN in arc lengths
    if (isnan(s_top) || isnan(s_bot) || isnan(s_right) || isnan(s_left)) {
        return -1;
    }
    // if they're all negative, use min of absolute value
    if (s_top && s_bot && s_right && s_left < 0) {
            s = fabs(s_top);
        if (fabs(s_bot) < s) {
            s = fabs(s_bot);
        } else if (fabs(s_right) < s) {
            s = fabs(s_right);
        } else if (fabs(s_left) < s) {
            s = fabs(s_left);
        }
    // if any are positive, use positive min
    } else {
        s = 100;
        if (0 < s_top < s) {
            s = s_top;
        } else if (0 < s_bot < s) {
            s = s_bot;
        } else if (0 < s_right < s) {
            s = s_right;
        } else if (0 < s_left < s) {
            s = s_left;
        }
    }
    s += 0.08;
    double new_Px = line[m_idx - 1][0] + ((dx/mag)*s);
    double new_Py = line[m_idx - 1][1] + ((dy/mag)*s);
    if (fabs(new_Px - line[m_idx - 1][0]) > 2. ||
              fabs(new_Py - line[m_idx - 1][1]) > 2.) {
        return -1;
    } else {
        line[m_idx][0] = new_Px;
        line[m_idx][1] = new_Py;
        segs[m_idx] = s;
    }
    return 0;
}

/** Creates a streamline centered on start coordinates **/
int streamline(vector_field_t *m_field, streamline_t *m_sl)
{
    m_sl->forward[0] = m_sl->start;
    m_sl->back[0] = m_sl->start;
    // advect forwards froms start
    for (int i = 0; i < PSTEPS; i++) {
        if (advectP(m_field, m_sl, i, 0) < 0) {
            break;
        }
    }
    // advect backwards from start
    for (int j = 0; j < PSTEPS; j++) {
        if (advectP(m_field, m_sl, j, 1) < 0) {
            break;
        }
    }
    return 0;
}

/** Convolution kernel **/
double kern(double k, double s, int hanning)
{
    // boxcar
    if (!hanning) { return k + s; }
    return k + ((cos(s + PI) / KLEN) + 1.)/2.;
}

int partialIntegral(vector_field_t *m_field, streamline_t *m_sl,
                double Fsum, double hsum, int hanning, int back)
{
    double *segs;
    double tex_val;
    if (back) {
        segs = m_sl->back_seg;
    } else { segs = m_sl->forward_seg; }
    double s, k0, k1 = 0.0;
    int L = sizeof(segs)/sizeof(double);
    for (int l = 0; l < L; l++) {
        for (int i = 1; i < L - 2; i++) {
            s += segs[i - 1];
            double s_plus = s + segs[i + 1];
            k1 = kern(k1, s_plus, 1);
            k0 = kern(k0, s, 1);
        }
        double h = k1 - k0;
        hsum += h;
        int m_grid_idx[2];
        if (back) {
            pixIdx(m_field, m_sl->back[l], m_grid_idx);
        } else { pixIdx(m_field, m_sl->forward[l], m_grid_idx); }
        tex_val = m_field->texture[m_grid_idx[0]][m_grid_idx[1]];
        Fsum += tex_val * h;
    }
    return 0;
}

int lic(vector_field_t *m_field, streamline_t *m_sl)
{
    double F_forward, h_forward, F_back, h_back, lic;
    // forward integral
    partialIntegral(m_field, m_sl, F_forward, h_forward, 1, 0);
    // backward integral
    partialIntegral(m_field, m_sl, F_back, h_back, 1, 1);
    if ((F_forward + F_back == 0.) || (h_forward + h_back == 0.)) {
        lic = 0.0;
    }
    lic = (F_forward + F_back) / (h_forward + h_back);
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
