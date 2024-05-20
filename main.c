#include "aeroacoustic.h"

void integration(state *pstate, grid *pgrid, params *ppar, buffer *pbuff, integr *pintegr);

int main()
{
    state state_1;
    state *pstate_1 = &state_1;
    grid grid_1;
    grid *pgrid_1 = &grid_1;
    params params_1;
    params *pparams_1 = &params_1;
    integr integr_1;
    integr *pintegr_1 = &integr_1;
    buffer buff_1;
    buffer *pbuff_1 = &buff_1;

    initialization(pgrid_1, pstate_1, pintegr_1, pparams_1);
    grid_gen(pgrid_1, pbuff_1);
    int i, j;

    printf("\nstarting, Ntime = %d\n", pintegr_1->Ntime);
    printf("\nstarting, Ntime_1 = %d\n", pintegr_1->Ntime_1);
    printf("\nstarting, Ntime_2 = %d\n", pintegr_1->Ntime_2);
    printf("\nc_0 = %lf\n", pparams_1->c_0);
    printf("\nu_0 = %lf\n", pparams_1->u_0);
    printf("\nr_0 = %lf\n", pbuff_1->r_0);
    printf("\nsig_max = %lf\n", pbuff_1->sig_max);
    printf("\nw= %lf\n", pbuff_1->w);
    printf("\ndy = %lf\n", pgrid_1->dy);
    printf("\ndx = %lf\n", pgrid_1->dx);
    printf("\nsin = %lf\n", exp(1. / 2.));

    alphas(pstate_1, pparams_1);
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            printf("\n ax_%d_%d = %lf ", i, j, gsl_matrix_get(pstate_1->Ax, i, j));
        }
    }
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            printf("\n ay_%d_%d = %lf ", i, j, gsl_matrix_get(pstate_1->Ay, i, j));
        }
    }
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            printf("\n ralx_%d_%d = %lf ", i, j, gsl_matrix_get(pstate_1->RALx, i, j));
        }
    }
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            printf("\n raly_%d_%d = %lf ", i, j, gsl_matrix_get(pstate_1->RALy, i, j));
        }
    }

    // source(pstate_1,pgrid_1,pparams_1,pbuff_1);

    integration(pstate_1, pgrid_1, pparams_1, pbuff_1, pintegr_1);
    terminate(pstate_1, pgrid_1, pparams_1, pbuff_1, pintegr_1);
    //
    // printf("\ndx = %lf\n",grid_1.dx);
    // printf("\ndy = %lf\n",grid_1.dy);
    // printf("\nw = %lf\n",buff_1.w);

    return 12;
}

void integration(state *pstate, grid *pgrid, params *ppar, buffer *pbuff, integr *pintegr)
{

    int i, j, l, m, n;
    int t = 0;
    double rk[4] = {0.25, 0.33, 0.5052, 1.};
    double delta_rho;
    double delta_u;
    double delta_v;

    gsl_matrix *rho_0 = gsl_matrix_calloc(pgrid->N, pgrid->M);
    gsl_matrix *u_0 = gsl_matrix_calloc(pgrid->N, pgrid->M);
    gsl_matrix *v_0 = gsl_matrix_calloc(pgrid->N, pgrid->M);

    printf("\n time_tot = %lf\n", pintegr->time_tot);
    printf("\n source Hz= %lf\n", ppar->f_0);
    for (t = 1; t < pintegr->Ntime; t++)
    {

        for (l = 0; l < 4; l++)
        {
            roes(pstate, pgrid);
            flux(pstate, pgrid);
            source(pstate, pgrid, ppar, pbuff);
            for (i = 1; i < pgrid->N - 2; i++)
            {
                for (j = 1; j < pgrid->M - 2; j++)
                {
                    delta_rho = -pintegr->dt * ((pgrid->dy * (gsl_vector_get(pstate->Fx[i + 1][j], 0) - gsl_vector_get(pstate->Fx[i][j], 0)) / pgrid->dx) + (pgrid->dx * (gsl_vector_get(pstate->Fy[i][j + 1], 0) - gsl_vector_get(pstate->Fy[i][j], 0)) / pgrid->dy) - gsl_matrix_get(pstate->Qm, i, j));
                    gsl_matrix_set(pstate->rho, i, j, gsl_matrix_get(rho_0, i, j) + rk[l] * delta_rho);
                    // printf("\n delta_rho_%d_%d = %lf", i, j, delta_rho);

                    delta_u = -pintegr->dt * ((pgrid->dy * (gsl_vector_get(pstate->Fx[i + 1][j], 1) - gsl_vector_get(pstate->Fx[i][j], 1)) / pgrid->dx) + (pgrid->dx * (gsl_vector_get(pstate->Fy[i][j + 1], 1) - gsl_vector_get(pstate->Fy[i][j], 1)) / pgrid->dy) - gsl_matrix_get(pstate->Qu, i, j));
                    gsl_matrix_set(pstate->u, i, j, gsl_matrix_get(u_0, i, j) + rk[l] * delta_u);

                    delta_v = -pintegr->dt * ((pgrid->dy * (gsl_vector_get(pstate->Fx[i + 1][j], 2) - gsl_vector_get(pstate->Fx[i][j], 2)) / pgrid->dx) + (pgrid->dx * (gsl_vector_get(pstate->Fy[i][j + 1], 2) - gsl_vector_get(pstate->Fy[i][j], 2)) / pgrid->dy) - gsl_matrix_get(pstate->Qv, i, j));
                    gsl_matrix_set(pstate->v, i, j, gsl_matrix_get(v_0, i, j) + rk[l] * delta_v);
                }
            }

            // printf("\nset\n");
        }

        pstate->time_s = t * pintegr->dt;
        FILE *ff = fopen(FILE_FIELD, "w");
        fprintf(ff, "\n x y r rho u v Qm Qu zone\n");
        for (i = 0; i < pgrid->N - 1; i++)
        {
            for (j = 0; j < pgrid->M - 1; j++)

            {
                gsl_matrix_set(rho_0, i, j, gsl_matrix_get(pstate->rho, i, j));
                gsl_matrix_set(u_0, i, j, gsl_matrix_get(pstate->u, i, j));
                gsl_matrix_set(v_0, i, j, gsl_matrix_get(pstate->v, i, j));
                if ((gsl_matrix_get(pgrid->r, i, j) - (pbuff->r_0)) < 1e-9)
                {
                    fprintf(ff, "\n %g %g %g %g %g %g %g %g %g 00", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                }
                else
                {

                    fprintf(ff, "\n %g %g %g %g %g %g %g %g %g 10", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                    // printf("\nQm = %g", gsl_matrix_get(pstate->Qu,i,j));
                }
            }
        }
        fclose(ff);

        if (t == pintegr->Ntime_1)
        {
        FILE *f1 = fopen(FILE_FIELD_1, "w");
            fprintf(f1, "\n x y r rho u v Qm Qu zone\n");
            for (i = 0; i < pgrid->N - 1; i++)
            {
                for (j = 0; j < pgrid->M - 1; j++)

                {
                    gsl_matrix_set(rho_0, i, j, gsl_matrix_get(pstate->rho, i, j));
                    gsl_matrix_set(u_0, i, j, gsl_matrix_get(pstate->u, i, j));
                    gsl_matrix_set(v_0, i, j, gsl_matrix_get(pstate->v, i, j));
                    if ((gsl_matrix_get(pgrid->r, i, j) - (pbuff->r_0)) < 1e-9)
                    {
                        fprintf(f1, "\n %g %g %g %g %g %g %g %g %g 00", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                    }
                    else
                    {

                        fprintf(f1, "\n %g %g %g %g %g %g %g %g %g 10", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                        // printf("\nQm = %g", gsl_matrix_get(pstate->Qu,i,j));
                    }
                }
            }
        fclose(f1);
        }
        if (t == pintegr->Ntime_2)
        {
            FILE *f2 = fopen(FILE_FIELD_2, "w");
            fprintf(f2, "\n x y r rho u v Qm Qu zone\n");
            for (i = 0; i < pgrid->N - 1; i++)
            {
                for (j = 0; j < pgrid->M - 1; j++)

                {
                    gsl_matrix_set(rho_0, i, j, gsl_matrix_get(pstate->rho, i, j));
                    gsl_matrix_set(u_0, i, j, gsl_matrix_get(pstate->u, i, j));
                    gsl_matrix_set(v_0, i, j, gsl_matrix_get(pstate->v, i, j));
                    if ((gsl_matrix_get(pgrid->r, i, j) - (pbuff->r_0)) < 1e-9)
                    {
                        fprintf(f2, "\n %g %g %g %g %g %g %g %g %g 00", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                    }
                    else
                    {

                        fprintf(f2, "\n %g %g %g %g %g %g %g %g %g 10", pstate->time_s, gsl_matrix_get(pgrid->x, i, j), gsl_matrix_get(pgrid->y, i, j), gsl_matrix_get(pgrid->r, i, j), gsl_matrix_get(pstate->rho, i, j), gsl_matrix_get(pstate->u, i, j), gsl_matrix_get(pstate->v, i, j), gsl_matrix_get(pstate->Qm, i, j), gsl_matrix_get(pstate->Qu, i, j));
                        // printf("\nQm = %g", gsl_matrix_get(pstate->Qu,i,j));
                    }
                }
            }
        fclose(f2);
        }
        printf("\ntime = %lf", pstate->time_s);
    }
    gsl_matrix_free(u_0);
    gsl_matrix_free(rho_0);
    gsl_matrix_free(v_0);
}