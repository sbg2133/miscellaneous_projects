#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<errno.h>
#include<string.h>

#define XMAX 100 // Number of x coords in image
#define YMAX 100 // Number of y coords in image
#define LINE_LEN 1024 // Buffer length
#define KLEN 31 // Kernel length
#define PSTEPS 10 // Num streamline steps from start, for forward or back
#define PI 3.14159265

typedef struct {
    int xsize, ysize;
    double *X, *Y;
    double **dx, **dy;
    double **texture;
    double **image;
} vector_field_t;

typedef struct {
    double start[2];
    double forward[PSTEPS][2];
    double back[PSTEPS][2];
    double *forward_seg;
    double *back_seg;
    double segs[PSTEPS];
    double streamline[2*PSTEPS];
} streamline_t;

void coordGrid(vector_field_t *m_field)
{
    m_field->X = (double *)malloc(m_field->xsize*sizeof(double));
    m_field->Y = (double *)malloc(m_field->ysize*sizeof(double));
    // Make a linspace for X and Y
    for (int i = 0; i < m_field->xsize; i++) {
        m_field->X[i] = (i*(double)(XMAX)) / (double)m_field->xsize;
    }
    for (int j = 0; j < m_field->ysize; j++) {
        m_field->Y[j] = (j*(double)(YMAX)) / (double)m_field->ysize;
    }
    // Allocate memory for image values
    m_field->image = (double**)malloc(m_field->ysize * sizeof(double*));
    for (int i = 0; i < m_field->ysize; i++) {
         m_field->image[i] = (double *)malloc(m_field->xsize * sizeof(double));
    }
}

void loadVectors(vector_field_t *m_field)
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

void loadTexture(vector_field_t *m_field)
{
    m_field->texture = (double**)malloc(m_field->ysize * sizeof(double*));
    for (int i = 0; i < m_field->ysize; i++) {
         m_field->texture[i] = (double *)malloc(m_field->xsize * sizeof(double));
    }
    FILE* tex_fd = fopen("texture.dat", "r");
    if (!tex_fd) {
        perror("Error: ");
	exit(EXIT_FAILURE);
    }
    for (int i = 0; i < m_field->ysize; i++) {
        for (int j = 0; j < m_field->xsize; j++) {
            fscanf(tex_fd, "%lf", &m_field->texture[i][j]);
        }
    }
    fclose(tex_fd);
}

/** Given the input coordinates, return the grid index **/
int *pixIdx(vector_field_t *m_field, double *m_coords)
{
    int *grid_idx = malloc(2*sizeof(int));
    // printf("Inside pixIdx...\n");
    double tempx = (double)XMAX;
    double tempy = (double)YMAX;
    double diffx, diffy;
    int i, j;
    // printf("xsize = %d\n", m_field->xsize);
    // printf("ysize = %d\n", m_field->ysize);
    // for (int i = 0; i < m_field->xsize; i++) { printf("Xi = %f\n", m_field->X[i]); }
    // for (int j = 0; j < m_field->ysize; j++) { printf("Yi = %f\n", m_field->Y[j]); }
    // printf("Coords = %f, %f\n", m_coords[0], m_coords[1]);
    for (int i = 0; i < m_field->xsize; i++) {
	diffx = fabs(m_field->X[i] - m_coords[0]);
        // printf("Diffx = %f\n", diffx);
	if (diffx < tempx) {
	    tempx = diffx;
        } else {
            grid_idx[0] = i;
            break;
        }
    }
    for (j = 0; j < m_field->ysize; j++) {
        diffy = fabs(m_field->Y[j] - m_coords[1]);
        // printf("Diffy = %f\n", diffy);
        if (diffy < tempy) {
	    tempy = diffy;
        } else {
            grid_idx[1] = j;
            break;
        }
    }
    // printf("returned idx = %d, %d\n", grid_idx[0], grid_idx[1]);
    return grid_idx;
}

double *getVector(vector_field_t *m_field, double *m_coords)
{
    int *idx = malloc(2*sizeof(int));
    double *vector = malloc(3*sizeof(double));
    // printf("Inside getVector\n");
    idx = pixIdx(m_field, m_coords);
    // printf("idx = %d, %d\n", idx[0], idx[1]);
    vector[0] = m_field->dx[idx[1]][idx[0]];
    vector[1] = m_field->dy[idx[1]][idx[0]];
    vector[2] = atan2(vector[1], vector[0]) * 180. / PI;
    return vector;
}

int advectP(vector_field_t *m_field, streamline_t *sl, double line[PSTEPS][2],
           double *segs, int m_idx, int back)
{
    // printf("Inside advectP\n");
    double s;
    double *vector;
    segs[0] = 0.0;
    if (m_idx == 1) {
        vector = getVector(m_field, sl->start);
    } else {
        vector = getVector(m_field, line[m_idx - 1]);
        // if (isnan(vector[0]) || isnan(vector[1])) break;
    }
    double dx = vector[0];
    double dy = vector[1];
    if (back) {
        dx *= -1;
        dy *= -1;
    }
    // printf("forward start = %f\n", sl->forward[m_idx - 1][0]);
    // printf("line start = %f\n", line[m_idx - 1][0]);
    int Px = (int)line[m_idx - 1][0];
    int Py = (int)line[m_idx - 1][1];
    double mag = sqrt(dx*dx + dy*dy);
    // arc lengths
    double s_top = ((Py + 1) - line[m_idx - 1][1])*(mag/dy);
    double s_bot = (Py - line[m_idx - 1][1])*(mag/dy);
    double s_right = ((Px + 1) - line[m_idx - 1][0])*(mag/dx);
    double s_left = (Px - line[m_idx - 1][0])*(mag/dx);
    // check for NaN in arc lengths
    if (isnan(s_top) || isnan(s_bot) || isnan(s_right) || isnan(s_left)) {
        return -1;
    }
    // if they're all negative, use min of absolute value
    if ((s_top < 0) && (s_bot < 0) && (s_right < 0) && (s_left < 0)) {
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
        printf("s = %f\t", s);
    }
    /*if (sl->start[0] == 0.0 && sl->start[1] == 0.0) {
        printf("line = %f, %f \n", line[m_idx][0], line[m_idx][1]);
        // printf("forward = %f, %f \n", sl->forward[m_idx][0], sl->forward[m_idx][1]);
        // for (int k = 0; k < PSTEPS; k++) {
        //    printf("Forward: %f, %f\n", sl->forward[m_idx][k], sl->forward[m_idx][k]);
    }*/
    return 0;
}

/** Creates a streamline centered on start coordinates **/
int streamline(vector_field_t *m_field, streamline_t *m_sl)
{
    // m_sl->forward[0][0] = m_sl->start[0];
    // m_sl->forward[0][1] = m_sl->start[1];
    // m_sl->back[0][0] = m_sl->start[0];
    // m_sl->back[0][1] = m_sl->start[1];
    // advect forwards froms start
    double s_forward[PSTEPS];
    double temp_forward[PSTEPS][2];
    temp_forward[0][0] = m_sl->start[0];
    temp_forward[0][1] = m_sl->start[1];
    int i;
    for (i = 1; i < PSTEPS; i++) {
        if (advectP(m_field, m_sl, temp_forward, s_forward, i, 0) < 0) {
            break;
        }
    }
    for (int k = 0; k < i; k++) {
        // printf("Tempf: %f, %f\n", temp_forward[k][0], temp_forward[k][1]);
        memcpy(m_sl->forward[k], temp_forward[k], sizeof(m_sl->forward)*2);
        // printf("forward: %f, %f\n", m_sl->forward[k][0], m_sl->forward[k][1]);
        // printf("forward seg: %f\n", s_forward[k]);
    }
    m_sl->forward_seg = s_forward;
    // advect backwards from start
    double s_back[PSTEPS];
    double temp_back[PSTEPS][2];
    temp_back[0][0] = m_sl->start[0];
    temp_back[0][1] = m_sl->start[1];
    int j;
    for (j = 1; j < PSTEPS; j++) {
        if (advectP(m_field, m_sl, temp_back, s_back, j, 1) < 0) {
            break;
        }
    }
    for (int k = 0; k < j; k++) {
        // printf("Tempb: %f, %f\n", temp_back[k][0], temp_back[k][1]);
        memcpy(m_sl->back[k], temp_back[k], sizeof(m_sl->back)*2);
        // printf("back: %f, %f\n", m_sl->back[k][0], m_sl->back[k][1]);
        // printf("back seg: %f\n", s_back[k]);
    }
    m_sl->back_seg = s_back;
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
    // printf("Inside partial integral...\n");
    int *m_grid_idx;
    double *segs;
    double tex_val;
    if (back) {
        segs = m_sl->back_seg;
    } else { segs = m_sl->forward_seg; }
    double s, k0, k1 = 0.0;
    double s1;
    for (int l = 0; l < PSTEPS; l++) {
        for (int i = 1; i < PSTEPS - 2; i++) {
            printf("seg = %f\t", segs[i -1]);
            s += segs[i - 1];
            s1 = s + segs[i + 1];
            k1 = kern(k1, s1, 1);
            k0 = kern(k0, s, 1);
        }
        double h = k1 - k0;
        hsum += h;
        // int m_grid_idx[2];
        if (back) {
            // printf("back coords: %f, %f\n", m_sl->back[l][0], m_sl->back[l][1]);
            m_grid_idx = pixIdx(m_field, m_sl->back[l]);
        } else {
            // printf("forward coords: %f, %f\n", m_sl->forward[l][0], m_sl->forward[l][1]);
            m_grid_idx = pixIdx(m_field, m_sl->forward[l]);
        }
        // printf("grid idx = %d %d\n", m_grid_idx[0], m_grid_idx[1]);
        tex_val = m_field->texture[m_grid_idx[0]][m_grid_idx[1]];
        Fsum += tex_val * h;
    }
    return 0;
}

int lic(vector_field_t *m_field)
{
    // printf("xsize = %d\n", m_field->xsize);
    // printf("ysize = %d\n", m_field->ysize);
    double lics[m_field->xsize*m_field->ysize];
    int *m_grid_idx;
    int idx = 0;
    // int m_grid_idx[2];
    for (int i = 0; i < m_field->xsize; i++) {
        for (int j = 0; j < m_field->ysize; j++) {
            streamline_t sl;
            double F_forward, h_forward, F_back, h_back, lic;
            sl.start[0] = m_field->X[i];
            sl.start[1] = m_field->Y[j];
            m_grid_idx = pixIdx(m_field, sl.start);
            streamline(m_field, &sl);
            // forward integral
            partialIntegral(m_field, &sl, F_forward, h_forward, 1, 0);
            // backward integral
            partialIntegral(m_field, &sl, F_back, h_back, 1, 1);
            printf("F_forward, h_forward = %f, %f\n", F_forward, h_forward);
            if ((F_forward + F_back == 0) || (h_forward + h_back == 0)) {
                lic = 0.0;
            } else { lic = (F_forward + F_back) / (h_forward + h_back); }
            lics[idx] = lic;
            if ((idx > 0) && (lic == 0.0)) { lic = lics[idx - 1]; }
            // printf("Idx = %d, %d\n", m_grid_idx[0], m_grid_idx[1]);
            m_field->image[m_grid_idx[0]][m_grid_idx[1]] = lic;
            idx++;
            printf("Start = %f, %f\n", sl.start[0], sl.start[1]);
            printf("LIC = %f\n", lic);
            /*for (int k = 0; k < PSTEPS; i++) {
                printf("Forward: %f, %f\n", sl.forward[0][k], sl.forward[0][k + 1]);
            }*/
        }
    }
    return 0;
}

int main(int argc, char *argv[])
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
    // for (int i = 0; i < field.xsize; i++) { printf("Xi = %f\n", field.X[i]); }
    // for (int j = 0; j < field.ysize; j++) { printf("Yi = %f\n", field.Y[j]); }

    loadTexture(&field);
    // Load in the vector arrays, dx and dy
    loadVectors(&field);
    /* for (int i = 0; i < ysize; i++) {
        for (int j = 0; j < xsize; j++) {
            printf("%f, %f\n", field.dx[i][j], field.dy[i][j]);
        }
    } */

    lic(&field);
    /* for (int i = 0; i < field.ysize; i++) {
        for (int j = 0; j < field.xsize; j++) {
            printf("%f, %f\n", field.image[i][j], field.image[i][j]);
        }
    }*/
    return 0;
}
