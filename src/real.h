#ifdef REAL_DBLE
#define REAL_BYTE 8
#define REAL_MPI  MPI_DOUBLE_PRECISION
#define REAL_FUNC dble
#else
#define REAL_BYTE 4
#define REAL_MPI  MPI_REAL
#define REAL_FUNC real
#endif
