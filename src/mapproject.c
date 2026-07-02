#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifndef PROJ9
// =================================================================================================
// =================================================================================================
// for proj7 or earlier
// =================================================================================================
// =================================================================================================
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H // for newer proj.4 version 6.0.0
#include <proj_api.h>

static  projPJ pj_tmerc, pj_latlong;

void mapproject_initialize_(double *lon_0){
   char s[128];
// pj_tmerc   = pj_init_plus("+proj=tmerc   +lon_0=144.35 +k=1.00010001 +ellps=WGS84 +datum=WGS84 +no_defs");
   sprintf(s, "%s%f%s\n", "+proj=tmerc +lon_0=", *lon_0, " +k=1.00010001 +ellps=WGS84 +datum=WGS84 +no_defs");
   pj_tmerc   = pj_init_plus(s);
   pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs");
}

void mapproject_finalize_(){
   pj_free(pj_tmerc);
   pj_free(pj_latlong);
}

void ll2xy_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
   int i, j, ind;
   for(j=0; j<*ny; j++){
      for(i=0; i<*nx; i++){
         ind = j*(*nx)+i;
         *(x+ind) = (*(lon+ind))*DEG_TO_RAD;
         *(y+ind) = (*(lat+ind))*DEG_TO_RAD;
      }
   }
   pj_transform(pj_latlong, pj_tmerc, (*nx)*(*ny), 1, x, y, NULL);
}

void xy2ll_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
   int i, j, ind;
   for(j=0; j<*ny; j++){
      for(i=0; i<*nx; i++){
         ind = j*(*nx)+i;
         *(lon+ind) = *(x+ind);
         *(lat+ind) = *(y+ind);
      }
   }
   pj_transform(pj_tmerc, pj_latlong, (*nx)*(*ny), 1, lon, lat, NULL);
   for(j=0; j<*ny; j++){
      for(i=0; i<*nx; i++){
         ind = j*(*nx)+i;
         *(lon+ind) /= DEG_TO_RAD;
         *(lat+ind) /= DEG_TO_RAD;
      }
   }
}
#else
// =================================================================================================
// =================================================================================================
// for proj9 or later
// =================================================================================================
// =================================================================================================
#include <math.h>
#include <proj.h>

static PJ_CONTEXT *C;
static PJ *P_tmerc;
static PJ *P_latlong;

#define DEG_TO_RAD (M_PI/180.0)
#define RAD_TO_DEG (180.0/M_PI)

void mapproject_initialize_(double *lon_0){
    char s[128];

    C = proj_context_create();

    sprintf(s, "+proj=tmerc +lon_0=%f +k=1.00010001 +ellps=WGS84 +datum=WGS84 +no_defs", *lon_0);
    P_tmerc   = proj_create(C, s);
    P_latlong = proj_create(C, "+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs");
}

void mapproject_finalize_(){
    proj_destroy(P_tmerc);
    proj_destroy(P_latlong);
    proj_context_destroy(C);
}

void ll2xy_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
    int N = (*nx) * (*ny);

    for(int i = 0; i < N; i++){
        PJ_COORD a = proj_coord(lon[i] * DEG_TO_RAD, lat[i] * DEG_TO_RAD, 0, 0);
        PJ_COORD b = proj_trans(P_tmerc, PJ_FWD, a);
        x[i] = b.xy.x;
        y[i] = b.xy.y;
    }
}

void xy2ll_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
    int N = (*nx) * (*ny);

    for(int i = 0; i < N; i++){
        PJ_COORD a = proj_coord(x[i], y[i], 0, 0);
        PJ_COORD b = proj_trans(P_tmerc, PJ_INV, a);
        lon[i] = b.lp.lam * RAD_TO_DEG;
        lat[i] = b.lp.phi * RAD_TO_DEG;
    }
}
#endif
