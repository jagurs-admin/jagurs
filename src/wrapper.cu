#include <stdio.h>
#include <cufft.h>

extern "C" void make_plan_d2z_(cufftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n; onembed[0] = n/2+1;
   cufftPlanMany(plan, 1, N, inembed, 1, n, onembed, 1, n/2+1, CUFFT_D2Z, howmany);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_d2z_(cufftHandle* plan, double *a_in, cuDoubleComplex *a_out){
   cufftExecD2Z(*plan, a_in, a_out);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void make_plan_z2z_(cufftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n; onembed[0] = n;
   cufftPlanMany(plan, 1, N, inembed, 1, n, onembed, 1, n, CUFFT_Z2Z, howmany);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2z_forward_(cufftHandle* plan, cuDoubleComplex *a_in, cuDoubleComplex *a_out){
   cufftExecZ2Z(*plan, a_in, a_out, CUFFT_FORWARD);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2z_backward_(cufftHandle* plan, cuDoubleComplex *a_in, cuDoubleComplex *a_out){
   cufftExecZ2Z(*plan, a_in, a_out, CUFFT_INVERSE);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void make_plan_z2d_(cufftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n/2+1; onembed[0] = n;
   cufftPlanMany(plan, 1, N, inembed, 1, n/2+1, onembed, 1, n, CUFFT_Z2D, howmany);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2d_(cufftHandle* plan, cuDoubleComplex *a_in, double *a_out){
   cufftExecZ2D(*plan, a_in, a_out);
   cudaDeviceSynchronize(); // It must be called!!!
}

extern "C" void destroy_plan_(cufftHandle* plan){
   cufftDestroy(*plan);
}

extern "C" void cudevicesynchronize_(){
   cudaDeviceSynchronize();
}

