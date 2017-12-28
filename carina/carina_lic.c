#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<errno.h>
#include<string.h>

#define PSTEPS 51 // Num streamline steps from start, for forward or back
#define HANNING 1 // Set to 1 to use hanning kernel, 0 for boxcar
#define PI 3.14159265358979323846

typedef struct vector_field {
    int xsize, ysize;
    double *X, *Y;
    double **dx, **dy;
    double **texture;
    double **image;
} vector_field_t;

typedef struct streamline {
    double start[2];
    double forward[PSTEPS][2];
    double back[PSTEPS][2];
    double forward_seg[PSTEPS];
    double back_seg[PSTEPS];
    double forward_len;
    double back_len;
} streamline_t;

void coordGrid(vector_field_t *m_field)
{
    m_field->X = (double *)malloc(m_field->xsize*sizeof(double));
    m_field->Y = (double *)malloc(m_field->ysize*sizeof(double));
    // Make a linspace for X and Y
    for (int i = 0; i < m_field->xsize; i++) {
        m_field->X[i] = (i*(double)(m_field->xsize)) / ((double)m_field->xsize - 1.);
    }
    for (int j = 0; j < m_field->ysize; j++) {
        m_field->Y[j] = (j*(double)(m_field->ysize)) / ((double)m_field->ysize - 1.);
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
    m_field->dx = (double**)malloc(m_field->xsize * sizeof(double*));
    m_field->dy = (double**)malloc(m_field->ysize * sizeof(double*));
    // double *vx[ysize], *vy[xsize];
    int i, j;
    for (i = 0; i < m_field->ysize; i++) {
         m_field->dx[i] = (double *)malloc(m_field->xsize * sizeof(double));
         m_field->dy[i] = (double *)malloc(m_field->ysize * sizeof(double));
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
            // printf("%lf\n", m_field->texture[0][0]);
        }
    }
    fclose(tex_fd);
}

/** Given the input coordinates, return the grid index **/
int *pixIdx(vector_field_t *m_field, double *m_coords)
{
    // printf("ysize = %d\n", m_field->ysize);
    int *grid_idx = malloc(2*sizeof(int));
    // printf("Inside pixIdx...\n");
    double tempx = (double)m_field->xsize;
    double tempy = (double)m_field->ysize;
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
    // segs[0] = 0.0;
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
    // printf("%f, %f, %f, %f\n", s_top, s_bot, s_right, s_left);
    // check for NaN in arc lengths
    if (isnan(s_top) || isnan(s_bot) || isnan(s_right) || isnan(s_left)) {
        return -1;
    }
    // if they're all negative, use min of absolute value
    if ((s_top <= 0.) && (s_bot <= 0.) && (s_right <= 0.) && (s_left <= 0.)) {
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
        if ((s_top >= 0.) && (s_top <= s)) { s = s_top; }
        if ((s_bot >= 0.) && (s_bot <= s)) { s = s_bot; }
        if ((s_right >= 0.) && ( s_right <= s)) { s = s_right; }
        if ((s_left >= 0.) && ( s_left <= s)) { s = s_left; }
    }
    s += 0.08;
    double new_Px = line[m_idx - 1][0] + ((dx/mag)*s);
    double new_Py = line[m_idx - 1][1] + ((dy/mag)*s);
    if ((fabs(new_Px - line[m_idx - 1][0]) > 2.) ||
              (fabs(new_Py - line[m_idx - 1][1]) > 2.)) {
        return -1;
    } else {
        line[m_idx][0] = new_Px;
        line[m_idx][1] = new_Py;
        segs[m_idx - 1] = s;
    }
    // printf("%f\n", segs[m_idx]);
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
    int forward_len;
    int i;
    for (i = 1; i < PSTEPS; i++) {
        if (advectP(m_field, m_sl, temp_forward, s_forward, i, 0) < 0) {
            // printf("i = %d\n", i);
            m_sl->forward_len = i;
            break;
        }
    }
    for (int k = 0; k < PSTEPS; k++) {
        // printf("Tempf: %f, %f\n", temp_forward[k][0], temp_forward[k][1]);
        memcpy(m_sl->forward[k], temp_forward[k], sizeof(m_sl->forward[k]));
        // printf("forward: %f, %f\n", m_sl->forward[k][0], m_sl->forward[k][1]);
        // printf("forward seg: %f\n", s_forward[k]);
    }
    memcpy(m_sl->forward_seg, s_forward, sizeof(m_sl->forward_seg));
    //m_sl->forward_seg = s_forward;
    // advect backwards from start
    double s_back[PSTEPS];
    double temp_back[PSTEPS][2];
    temp_back[0][0] = m_sl->start[0];
    temp_back[0][1] = m_sl->start[1];
    int back_len;
    int j;
    for (j = 1; j < PSTEPS; j++) {
        if (advectP(m_field, m_sl, temp_back, s_back, j, 1) < 0) {
            // printf("j = %d\n", j);
            m_sl->back_len = j;
            break;
        }
    }
    for (int k = 0; k < PSTEPS; k++) {
        // printf("Tempb: %f, %f\n", temp_back[k][0], temp_back[k][1]);
        memcpy(m_sl->back[k], temp_back[k], sizeof(m_sl->back[k]));
        // printf("back: %f, %f\n", m_sl->back[k][0], m_sl->back[k][1]);
    }
    memcpy(m_sl->back_seg, s_back, sizeof(m_sl->back_seg));
    /* for (int k = 0; k < i; k++) {
        printf("temp, sl: %f, %f\n", s_forward[k], m_sl->forward_seg[k]);
    }*/
    /* for (int k = 0; k < i; k++) {
        printf("temp, sl: %f, %f\n", s_back[k], m_sl->back_seg[k]);
    }*/
    return 0;
}

/** Convolution kernel **/
double kern(int i)
{
    // boxcar
    if (HANNING) {
        return 0.5*(1. - cos(i * 2*PI / (PSTEPS - 1)));
    } else {
        return 1.;
    }
}

double *partialIntegral(vector_field_t *m_field, streamline_t *m_sl, int back)
{
    /* for (int k = 0; k < PSTEPS; k++) {
        printf("seg = %f\n", m_sl->forward_seg[k]);
    }*/
    double *sums = malloc(3*sizeof(double));
    int *m_grid_idx;
    double *segs;
    double tex_val;
    double h;
    int N;
    if (back) {
        segs = m_sl->back_seg;
        N = m_sl->back_len;
    } else {
        segs = m_sl->forward_seg;
        N = m_sl->forward_len;
    }
    /* for (int k = 0; k < PSTEPS; k++) {
        printf("k = %d, forward_seg = %f\n", k, m_sl->forward_seg[k]);
        // printf("temp = %f\n", segs[k]);
    }*/
    double s, k0, k1 = 0.0;
    double s1;
    for (int l = 0; l < PSTEPS; l++) {
        /* for (int i = 1; i < PSTEPS - 1; i++) {
            // printf("seg = %f\t", segs[i -1]);
            s += segs[i - 1];
            s1 = s + segs[i + 1];
            // s += segs[i - 1];
            // k1 += s1*kern(k1, i + 1, N);
            // k0 += s*kern(k0, i, N);
            k1 = kern(i + 
        } */
        h = kern(l);
        sums[1] += h;
        if (back) {
            // printf("back coords: %f, %f\n", m_sl->back[l][0], m_sl->back[l][1]);
            m_grid_idx = pixIdx(m_field, m_sl->back[l]);
        } else {
            // printf("forward coords: %f, %f\n", m_sl->forward[l][0], m_sl->forward[l][1]);
            m_grid_idx = pixIdx(m_field, m_sl->forward[l]);
        }
        // printf("grid idx = %d %d\n", m_grid_idx[0], m_grid_idx[1]);
        tex_val = m_field->texture[m_grid_idx[0]][m_grid_idx[1]];
        sums[0] += (tex_val * h);
    }
    // printf("Fsum = %g\t", sums[0]);
    return sums;
}

int lic(vector_field_t *m_field)
{
    // printf("xsize = %d\n", m_field->xsize);
    // printf("ysize = %d\n", m_field->ysize);
    FILE* out = fopen("lic.dat", "wt");
    double lics[m_field->xsize*m_field->ysize];
    int idx = 0;
    // int m_grid_idx[2];
    // printf("ysize = %d\n", m_field->ysize);
    for (int i = 0; i < m_field->xsize; i++) {
        for (int j = 0; j < m_field->ysize; j++) {
            double *forward_sums, *back_sums;
            int *m_grid_idx;
            streamline_t sl;
            double lic;
            sl.start[0] = m_field->X[i];
            sl.start[1] = m_field->Y[j];
            streamline(m_field, &sl);
            /* for (int k = 0; k < PSTEPS; k++) {
                    printf("%f\t", sl.forward_seg[k]);
            }
            printf("\n");*/
            // forward integral
            // printf("i, j, ysize = %d, %d, %d\n", i, j, m_field->ysize);
            forward_sums = partialIntegral(m_field, &sl, 0);
            // backward integral
            back_sums = partialIntegral(m_field, &sl, 1);
            // printf("F_forward, h_forward = %f, %f\n", forward_sums[0], forward_sums[1]);
            if ((forward_sums[0] + back_sums[0] == 0.) ||
                    (forward_sums[1] + back_sums[1] == 0.)) {
                lic = 0.0;
            } else if ((idx > 0) && (lic == 0.0)) {
                lics[idx] = lic;
            } else { lic = (forward_sums[0] + back_sums[0])
                                   /(forward_sums[1] + back_sums[1]);
            }
            m_grid_idx = pixIdx(m_field, sl.start);
            // printf("Idx = %d, %d\n", m_grid_idx[0], m_grid_idx[1]);
            m_field->image[m_grid_idx[0]][m_grid_idx[1]] = lic;
            /* printf("Start = %f, %f\n", sl.start[0], sl.start[1]);
            for (int k = 0; k < PSTEPS; k++) {
                printf("Forward: %f, %f\n", sl.forward[k][0], sl.forward[k][1]);
            }*/
            // printf("LIC = %f\n", lic);
            /*for (int k = 0; k < PSTEPS; k++) {
                    printf("segs: %f\t", sl.forward_seg[k]);
            }*/
            idx++;
            /* if (sl.start[0] == 0.0 && sl.start[1] == 0.0) {
                printf("Start = %f, %f\n", sl.start[0], sl.start[1]);
                printf("dx, dy = %f, %f\n", m_field->dx[0][0], m_field->dy[0][0]);
                printf("LIC = %f\n", lic);
                for (int k = 0; k < PSTEPS; k++) {
                    printf("Forward: %f, %f\t", sl.forward[k][0], sl.forward[k][1]);
                    printf("\n");
                    printf("segs: %f\t", sl.forward_seg[k]);
                }
            }*/
        }
        for (int k = 0; k < m_field->xsize; k++) {
            fprintf(out, "%g ", m_field->image[i][k]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
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
    /* for (int i = 0; i < field.ysize/10; i++) {
        for (int j = 0; j < field.xsize/10; j++) {
            printf("%f, %f\n", field.dx[i][j], field.dy[i][j]);
            // printf("%g\n", field.texture[i][j]);
        }
    }*/

    lic(&field);
    /* for (int i = 0; i < field.ysize; i++) {
        for (int j = 0; j < field.xsize; j++) {
            printf("%f\n", field.image[i][j]);
        }
    }*/
    return 0;
}
