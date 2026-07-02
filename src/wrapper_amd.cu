#include <stdio.h>
#include <hipfft/hipfft.h>

extern "C" void make_plan_d2z_(hipfftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n; onembed[0] = n/2+1;
   hipfftPlanMany(plan, 1, N, inembed, 1, n, onembed, 1, n/2+1, HIPFFT_D2Z, howmany);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_d2z_(hipfftHandle* plan, double *a_in, hipDoubleComplex *a_out){
   hipfftExecD2Z(*plan, a_in, a_out);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void make_plan_z2z_(hipfftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n; onembed[0] = n;
   hipfftPlanMany(plan, 1, N, inembed, 1, n, onembed, 1, n, HIPFFT_Z2Z, howmany);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2z_forward_(hipfftHandle* plan, hipDoubleComplex *a_in, hipDoubleComplex *a_out){
   hipfftExecZ2Z(*plan, a_in, a_out, HIPFFT_FORWARD);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2z_backward_(hipfftHandle* plan, hipDoubleComplex *a_in, hipDoubleComplex *a_out){
   hipfftExecZ2Z(*plan, a_in, a_out, HIPFFT_BACKWARD);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void make_plan_z2d_(hipfftHandle* plan, int *n_, int *howmany_){
   int n, howmany;
   int N[1], inembed[1], onembed[1];
   n = *n_;
   howmany = *howmany_;
   N[0] = n; inembed[0] = n/2+1; onembed[0] = n;
   hipfftPlanMany(plan, 1, N, inembed, 1, n/2+1, onembed, 1, n, HIPFFT_Z2D, howmany);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void exec_z2d_(hipfftHandle* plan, hipDoubleComplex *a_in, double *a_out){
   hipfftExecZ2D(*plan, a_in, a_out);
   hipDeviceSynchronize(); // It must be called!!!
}

extern "C" void destroy_plan_(hipfftHandle* plan){
   hipfftDestroy(*plan);
}

extern "C" void cudevicesynchronize_(){
   hipDeviceSynchronize();
}

