#ifndef MYH_
#define MYH_

#define rows 3
#define FILE_FIELD "fields.txt"
#define FILE_FIELD_1 "fields_1.txt"
#define FILE_FIELD_2 "fields_2.txt"
#define FILE_X "x.txt"
#define FILE_Y "y.txt"
#define FILE_GRID "grid.txt"
#define FILE_DIMS "input.txt"
#define FILE_PARAMS "params.txt"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct
{
    gsl_matrix *r;
    gsl_matrix *x;
    gsl_matrix *y;
    double dx;
    double dy;
    int N;
    int M;

} grid;
typedef struct
{
    double w;
    double r_0;
    double m;
    double sig_max;
} buffer;

typedef struct
{
    gsl_matrix *rho;
    gsl_matrix *u;
    gsl_matrix *v;
    gsl_matrix *p;
    gsl_vector ***Fx;
    gsl_vector ***Fy;
    gsl_vector ***Roex;
    gsl_vector ***Roey;
    gsl_matrix *Qm;
    gsl_matrix *Qu;
    gsl_matrix *Qv;
    gsl_matrix *Ax;
    gsl_matrix *Ay;
    gsl_matrix *RALx;
    gsl_matrix *RALy;
    double time_s;

} state;
typedef struct
{
    double u_0;
    double c_0;
    double rho_0;
    double T_0;
    double f_0;
    double Mach;
    double a;

} params;

typedef struct
{
    double dt;
    double time_tot;
    double time_tot_1;
    double time_tot_2;
    int Ntime;
    int Ntime_1;
    int Ntime_2;

} integr;

void initialization(grid *pgrid, state *pstate, integr *pintegr, params *ppar);
void grid_gen(grid *pgrid, buffer *pbuff);
void alphas(state *pstate, params *ppar);
void roes(state *pstate, grid *pgrid);
void flux(state *pstate, grid *pgrid);
void source(state *pstate, grid *pgrid, params *ppar, buffer *pbuff);
void terminate(state *pstate, grid *pgrid, params *ppar, buffer *pbuff, integr *pintegr);
#endif
