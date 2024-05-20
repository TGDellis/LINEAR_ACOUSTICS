#include "aeroacoustic.h"
const double air_gamma = 1.4;
const double air_R = 287.05;

void initialization(grid *pgrid, state *pstate, integr *pintegr, params *ppar)
{
    int i, j;
    FILE *f_i = fopen(FILE_DIMS, "r");
    fscanf(f_i, "%d %d %lf %lf %lf %lf", &pgrid->N, &pgrid->M, &pintegr->dt, &pintegr->time_tot, &pintegr->time_tot_1, &pintegr->time_tot_2);
    fclose(f_i);

    pintegr->Ntime = pintegr->time_tot / pintegr->dt;
    pintegr->Ntime_1 = pintegr->time_tot_1 / pintegr->dt;
    pintegr->Ntime_2 = pintegr->time_tot_2 / pintegr->dt;

    pstate->u = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->v = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->p = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->rho = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->Qm = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->Qu = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pstate->Qv = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pgrid->x = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pgrid->y = gsl_matrix_calloc(pgrid->N, pgrid->M);
    pgrid->r = gsl_matrix_calloc(pgrid->N, pgrid->M);

    pstate->time_s = 0.;
    pstate->Fx = malloc((pgrid->N - 1) * sizeof(gsl_vector **));
    pstate->Fy = malloc((pgrid->N - 1) * sizeof(gsl_vector **));
    pstate->Roex = malloc((pgrid->N - 1) * sizeof(gsl_vector **));
    pstate->Roey = malloc((pgrid->N - 1) * sizeof(gsl_vector **));
    FILE *f_p = fopen(FILE_PARAMS, "r");
    fscanf(f_p, "%lf %lf %lf %lf", &ppar->Mach, &ppar->T_0, &ppar->f_0, &ppar->a);
    fclose(f_p);

    ppar->c_0 = sqrt(air_gamma * air_R * ppar->T_0);
    ppar->u_0 = ppar->Mach * ppar->c_0;
    ppar->rho_0 = 1.225;
    for (i = 1; i < pgrid->N; i++)
    {
        pstate->Fx[i - 1] = malloc((pgrid->M - 1) * sizeof(gsl_vector *));
        pstate->Fy[i - 1] = malloc((pgrid->M - 1) * sizeof(gsl_vector *));
        pstate->Roex[i - 1] = malloc((pgrid->M - 1) * sizeof(gsl_vector *));
        pstate->Roey[i - 1] = malloc((pgrid->M - 1) * sizeof(gsl_vector *));

        for (j = 1; j < pgrid->M; j++)
        {
            pstate->Fx[i - 1][j - 1] = gsl_vector_calloc(rows);
            pstate->Fy[i - 1][j - 1] = gsl_vector_calloc(rows);
            pstate->Roex[i - 1][j - 1] = gsl_vector_calloc(rows);
            pstate->Roey[i - 1][j - 1] = gsl_vector_calloc(rows);
        }
    }
    pstate->Ax = gsl_matrix_calloc(rows, rows);
    pstate->Ay = gsl_matrix_calloc(rows, rows);
    pstate->RALx = gsl_matrix_calloc(rows, rows);
    pstate->RALy = gsl_matrix_calloc(rows, rows);
}

void grid_gen(grid *pgrid, buffer *pbuff)
{
    int i, j;
    FILE *fg = fopen(FILE_GRID, "r");
    double temp_x;
    double temp_y;
    double temp_r;
    for (i = 0; i < pgrid->N; i++)
    {
        for (j = 0; j < pgrid->N; j++)
        {
            fscanf(fg, "%lf %lf", &temp_x, &temp_y);
            temp_r = sqrt(pow(temp_x, 2) + pow(temp_y, 2));
            gsl_matrix_set(pgrid->x, i, j, temp_x);
            gsl_matrix_set(pgrid->y, i, j, temp_y);
            gsl_matrix_set(pgrid->r, i, j, temp_r);
        }
    }
    fclose(fg);
    pbuff->m = 2.;
    pbuff->r_0 = 100.;
    pbuff->w = gsl_matrix_get(pgrid->x, pgrid->N - 1, 0) - pbuff->r_0;
    pgrid->dx = fabs(gsl_matrix_get(pgrid->x, 1, 0) - gsl_matrix_get(pgrid->x, 0, 0));
    pgrid->dy = fabs(gsl_matrix_get(pgrid->y, 0, 1) - gsl_matrix_get(pgrid->x, 0, 0));
    pbuff->sig_max = 10. / pgrid->dx;
}

void alphas(state *pstate, params *ppar)
{
    gsl_matrix *Rx = gsl_matrix_calloc(rows, rows);
    gsl_matrix *Lx = gsl_matrix_calloc(rows, rows);
    gsl_matrix *Lamx = gsl_matrix_calloc(rows, rows);
    gsl_matrix *temp_x = gsl_matrix_calloc(rows, rows);
    gsl_matrix *Ry = gsl_matrix_calloc(rows, rows);
    gsl_matrix *Ly = gsl_matrix_calloc(rows, rows);
    gsl_matrix *Lamy = gsl_matrix_calloc(rows, rows);
    gsl_matrix *temp_y = gsl_matrix_calloc(rows, rows);

    gsl_matrix_set(pstate->Ax, 0, 0, ppar->u_0);
    gsl_matrix_set(pstate->Ax, 0, 1, ppar->rho_0);
    gsl_matrix_set(pstate->Ax, 1, 0, pow(ppar->c_0, 2) / (ppar->rho_0));
    gsl_matrix_set(pstate->Ax, 1, 1, ppar->u_0);
    gsl_matrix_set(pstate->Ax, 2, 2, ppar->u_0);

    gsl_matrix_set(pstate->Ay, 0, 2, ppar->rho_0);
    gsl_matrix_set(pstate->Ay, 2, 0, pow(ppar->c_0, 2) / (ppar->rho_0));

    gsl_matrix_set(Lamx, 0, 0, fabs(ppar->u_0));
    gsl_matrix_set(Lamx, 1, 1, fabs(ppar->u_0 + ppar->c_0));
    gsl_matrix_set(Lamx, 2, 2, fabs(ppar->u_0 - ppar->c_0));

    gsl_matrix_set(Rx, 0, 1, ppar->rho_0 / ppar->c_0);
    gsl_matrix_set(Rx, 0, 2, -ppar->rho_0 / ppar->c_0);
    gsl_matrix_set(Rx, 1, 1, 1.);
    gsl_matrix_set(Rx, 1, 2, 1.);
    gsl_matrix_set(Rx, 2, 0, 1.);

    gsl_matrix_set(Lx, 0, 2, 1.);
    gsl_matrix_set(Lx, 1, 0, ((ppar->c_0 / (2 * ppar->rho_0))));
    gsl_matrix_set(Lx, 1, 1, 0.5);
    gsl_matrix_set(Lx, 2, 1, 0.5);
    gsl_matrix_set(Lx, 2, 0, -((ppar->c_0 / (2 * ppar->rho_0))));

    gsl_matrix_set(Lamy, 1, 1, fabs(ppar->c_0));
    gsl_matrix_set(Lamy, 2, 2, fabs(-ppar->c_0));

    gsl_matrix_set(Ry, 0, 1, ppar->rho_0 / ppar->c_0);
    gsl_matrix_set(Ry, 0, 2, -ppar->rho_0 / ppar->c_0);
    gsl_matrix_set(Ry, 1, 0, 1.);
    gsl_matrix_set(Ry, 2, 1, 1.);
    gsl_matrix_set(Ry, 2, 2, 1.);

    gsl_matrix_set(Ly, 1, 0, ((ppar->c_0 / (2 * ppar->rho_0))));
    gsl_matrix_set(Ly, 2, 0, -((ppar->c_0 / (2 * ppar->rho_0))));
    gsl_matrix_set(Ly, 0, 1, 1.);
    gsl_matrix_set(Ly, 1, 2, 0.5);
    gsl_matrix_set(Ly, 2, 2, 0.5);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Rx, Lamx, 0.0, temp_x);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ry, Lamy, 0.0, temp_y);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp_x, Lx, 0.0, pstate->RALx);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp_y, Ly, 0.0, pstate->RALy);

    gsl_matrix_free(Lx);
    gsl_matrix_free(Rx);
    gsl_matrix_free(Lamx);
    gsl_matrix_free(temp_x);

    gsl_matrix_free(Ly);
    gsl_matrix_free(Ry);
    gsl_matrix_free(Lamy);
    gsl_matrix_free(temp_y);
}

void roes(state *pstate, grid *pgrid)
{
    int i, j;
    gsl_vector *Ux_diff = gsl_vector_calloc(rows);
    gsl_vector *Uy_diff = gsl_vector_calloc(rows);

    for (i = 1; i < pgrid->N - 2; i++)
    {
        for (j = 1; j < pgrid->M - 2; j++)
        {
            gsl_vector_set(Ux_diff, 0, gsl_matrix_get(pstate->rho, i, j) - gsl_matrix_get(pstate->rho, i - 1, j));
            gsl_vector_set(Uy_diff, 0, gsl_matrix_get(pstate->rho, i, j) - gsl_matrix_get(pstate->rho, i, j - 1));
            gsl_vector_set(Ux_diff, 1, gsl_matrix_get(pstate->u, i, j) - gsl_matrix_get(pstate->u, i - 1, j));
            gsl_vector_set(Uy_diff, 1, gsl_matrix_get(pstate->u, i, j) - gsl_matrix_get(pstate->u, i, j - 1));
            gsl_vector_set(Ux_diff, 2, gsl_matrix_get(pstate->v, i, j) - gsl_matrix_get(pstate->v, i - 1, j));
            gsl_vector_set(Uy_diff, 2, gsl_matrix_get(pstate->v, i, j) - gsl_matrix_get(pstate->v, i, j - 1));
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->RALx, Ux_diff, 0.0, pstate->Roex[i][j]);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->RALy, Uy_diff, 0.0, pstate->Roey[i][j]);
        }
    }
    gsl_vector_free(Ux_diff);
    gsl_vector_free(Uy_diff);
}
void flux(state *pstate, grid *pgrid)
{
    int i, j, k;
    gsl_vector *temp_u = gsl_vector_calloc(rows);
    gsl_vector *temp_uw = gsl_vector_calloc(rows);
    gsl_vector *temp_us = gsl_vector_calloc(rows);

    gsl_vector *temp_fx = gsl_vector_calloc(rows);
    gsl_vector *temp_fy = gsl_vector_calloc(rows);
    gsl_vector *temp_fw = gsl_vector_calloc(rows);
    gsl_vector *temp_fs = gsl_vector_calloc(rows);
    for (i = 1; i < pgrid->N - 2; i++)
    {
        for (j = 1; j < pgrid->M - 2; j++)
        {
            gsl_vector_set(temp_u, 0, gsl_matrix_get(pstate->rho, i, j));
            gsl_vector_set(temp_u, 1, gsl_matrix_get(pstate->u, i, j));
            gsl_vector_set(temp_u, 2, gsl_matrix_get(pstate->v, i, j));
            gsl_vector_set(temp_uw, 0, gsl_matrix_get(pstate->rho, i - 1, j));
            gsl_vector_set(temp_uw, 1, gsl_matrix_get(pstate->u, i - 1, j));
            gsl_vector_set(temp_uw, 2, gsl_matrix_get(pstate->v, i - 1, j));
            gsl_vector_set(temp_us, 0, gsl_matrix_get(pstate->rho, i, j - 1));
            gsl_vector_set(temp_us, 1, gsl_matrix_get(pstate->u, i, j - 1));
            gsl_vector_set(temp_us, 2, gsl_matrix_get(pstate->v, i, j - 1));

            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->Ax, temp_u, 0.0, temp_fx);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->Ay, temp_u, 0.0, temp_fy);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->Ax, temp_uw, 0.0, temp_fw);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->Ay, temp_us, 0.0, temp_fs);
            // printf("\ntemp_fx = %lf\n",gsl_vector_get(temp_fx,0));
            // printf("\ntemp_fx = %lf\n",gsl_vector_get(temp_fx,1));
            // printf("\ntemp_fx = %lf\n",gsl_vector_get(temp_fx,2));

            for (k = 0; k < rows; k++)
            {
                gsl_vector_set(pstate->Fx[i][j], k, 0.5 * (1 * gsl_vector_get(temp_fw, k) + 1 * gsl_vector_get(temp_fx, k) - 1 * gsl_vector_get(pstate->Roex[i][j], k)));
                gsl_vector_set(pstate->Fy[i][j], k, 0.5 * (1 * gsl_vector_get(temp_fs, k) + 1 * gsl_vector_get(temp_fy, k) - 1 * gsl_vector_get(pstate->Roey[i][j], k)));
            }
        }
    }
    gsl_vector_free(temp_fw);
    gsl_vector_free(temp_fx);
    gsl_vector_free(temp_fs);
    gsl_vector_free(temp_fy);
    gsl_vector_free(temp_u);
    gsl_vector_free(temp_uw);
    gsl_vector_free(temp_us);
}

void source(state *pstate, grid *pgrid, params *ppar, buffer *pbuff)
{
    int i, j;
    double temp_Sm;
    double temp_Su;
    double temp_Sv;
    double temp_damp;
    double temp_r;
    double temp_x;
    double temp_y;
    double tiny = 1e-9;

    for (i = 0; i < pgrid->N; i++)
    {
        for (j = 0; j < pgrid->M; j++)
        {
            temp_r = gsl_matrix_get(pgrid->r, i, j);
            temp_x = gsl_matrix_get(pgrid->x, i, j);
            temp_y = gsl_matrix_get(pgrid->y, i, j);

            temp_Sm = 0.;//sin(2. * M_PI * (ppar->f_0) * (pstate->time_s)) * exp(-(ppar->a) * pow(temp_r, 2)) / ppar->c_0;
            temp_Su = 0.; // sin(2. * M_PI * (ppar->f_0) * (pstate->time_s)) * exp(-(ppar->a) * pow(temp_r, 2)) / ppar->c_0;
            temp_Sv = 0.; // sin(2. * M_PI * (ppar->f_0) * (pstate->time_s)) * exp(-(ppar->a) * pow(temp_r, 2)) / ppar->c_0;
            if (fabs(temp_x) < 10 && fabs(temp_y) < 10)
            {
                temp_Su = sin(2. * M_PI * (ppar->f_0) * (pstate->time_s)) * cos(M_PI * temp_x / 20.) * exp(-(ppar->a) * pow(temp_y, 2));
                // temp_Sv = -sin(2. * M_PI * (ppar->f_0) * (pstate->time_s)) * sin(M_PI * temp_y / 20.) * exp(-(ppar->a) * pow(temp_x, 2));
            }
            // temp_Sv = 0;

            if ((temp_r - pbuff->r_0) > tiny)
            {
                temp_damp = pbuff->sig_max * (pow(((temp_r - pbuff->r_0) / pbuff->w), pbuff->m));

                temp_Sm = temp_Sm - temp_damp * (gsl_matrix_get(pstate->rho, i, j));
                temp_Su = temp_Su - temp_damp * (gsl_matrix_get(pstate->u, i, j));
                temp_Sv = temp_Sv - temp_damp * (gsl_matrix_get(pstate->v, i, j));
            }
            gsl_matrix_set(pstate->Qm, i, j, temp_Sm);
            gsl_matrix_set(pstate->Qu, i, j, temp_Su);
            gsl_matrix_set(pstate->Qv, i, j, temp_Sv);
        }
    }
}

void terminate(state *pstate, grid *pgrid, params *ppar, buffer *pbuff, integr *pintegr)
{
    gsl_matrix_free(pstate->u);
    gsl_matrix_free(pstate->v);
    gsl_matrix_free(pstate->rho);
    gsl_matrix_free(pstate->Qm);
    gsl_matrix_free(pstate->Qu);
    gsl_matrix_free(pstate->Qv);
    free(pstate->Fx);
    free(pstate->Fy);
    free(pstate->RALx);
    free(pstate->RALy);
}
