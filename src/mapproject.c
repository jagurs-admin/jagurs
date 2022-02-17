#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define ACCEPT_USE_OF_DEPRECATED_PROJ_API_H // for newer proj.4 version 6.0.0
#include <proj_api.h>

static  projPJ pj_tmerc, pj_latlong;

#ifndef _WIN32
void mapproject_initialize_(double *lon_0){
#else
void MAPPROJECT_INITIALIZE(double *lon_0){
#endif
   char s[128];
// pj_tmerc   = pj_init_plus("+proj=tmerc   +lon_0=144.35 +k=1.00010001 +ellps=WGS84 +datum=WGS84 +no_defs");
   sprintf(s, "%s%f%s\n", "+proj=tmerc +lon_0=", *lon_0, " +k=1.00010001 +ellps=WGS84 +datum=WGS84 +no_defs");
   pj_tmerc   = pj_init_plus(s);
   pj_latlong = pj_init_plus("+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs");
}

#ifndef _WIN32
void mapproject_finalize_(){
#else
void MAPPROJECT_FINALIZE(){
#endif
   pj_free(pj_tmerc);
   pj_free(pj_latlong);
}

#ifndef _WIN32
void ll2xy_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
#else
void LL2XY(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
#endif
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

#ifndef _WIN32
void xy2ll_(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
#else
void XY2LL(int *nx, int *ny, double *lon, double *lat, double *x, double *y){
#endif
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
