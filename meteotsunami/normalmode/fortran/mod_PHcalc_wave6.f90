module mod_PHcalc_wave6
   implicit none

   real(kind=8) :: srclon, srclat
   integer :: dist_s, dist_e
   real(kind=8) :: tmax, dt

   integer :: maxmode = 13   ! maximum mode number for dispersion curve and synthetic wave calculation
#ifndef JAGURS
   integer :: minmode = 8    ! minimum mode for synthetic wave calculation
#else
   integer :: minmode = 10   ! minimum mode for synthetic wave calculation
#endif

   logical, parameter :: with_Ocean = .true.  ! T: Calculate with ocean layer, F: Without ocean
   logical, parameter :: Only_Ocean = .false. ! T: Calculate only ocean layer, F: Calculate with atmospheric model

   logical :: DCflag = .true.      ! 1: Calculate dispersion curve, 0: Use a priori dispersion curve save as DCfile
   logical, parameter :: DCsaveflag = .true.  ! 1: Save dispersion curve in DCfile, 0: Not save

   logical, parameter :: DCplotflag = .false. ! 1: Plot dispersion curve (phase velocity), 0: Not plot
   logical :: SYNflag = .true.     ! 1: Calculate synthetic waveform, 0: Not calculate

   real(kind=8), parameter :: g0 = 9.8d-3
   real(kind=8), parameter :: Oalpha = 1.5d0
   real(kind=8), parameter :: Orho = 1030.0d0
   real(kind=8) :: OH = 5.0d0

   real(kind=8), parameter :: Re = 6371.0d0
   real(kind=8) :: Samp = 1000d0
   real(kind=8) :: SDt = 2000d0

   real(kind=8), parameter :: omg_searchmax = 0.01d0

   real(kind=8), parameter :: pi = 3.14159265d0
   real(kind=8), parameter :: NaN = transfer(Z'FFFFFFFFFFFFFFFF', 0.d0)
   complex(kind=8), parameter :: im = (0.0d0, 1.0d0)

   integer :: nz
   real(kind=8), allocatable, dimension(:) :: zl, H, alpha, rho, g, ganma, BVf, lamda, beta
   real(kind=8), allocatable, dimension(:) :: t, vp, omega
   real(kind=8) :: Df
   integer :: l_t, l_vp, l_omega
   integer :: sind = 10

   integer, allocatable, dimension(:) :: nkr
   real(kind=8), allocatable, dimension(:) :: kroot, vpkr, omgkr
   integer :: l_out

   character(len=256) :: DCfile, DCtext
#ifndef JAGURS
   character(len=256) :: outfile
#endif
   character(len=256) :: atmfile = 'AtmModel.dat'
contains

   subroutine makearray(a, n, vst, ven, d)
      real(kind=8), allocatable, dimension(:), intent(inout) :: a
      integer, intent(out) :: n
      real(kind=8), intent(in) :: vst, ven, d
      real(kind=8) :: val
      integer :: i
      n = 0
      do
         val = vst + n*d
         if(val > ven) exit
         n = n + 1
      end do
!     write(0,*) '[makearray] n: ', n
      allocate(a(n))
      do i = 1, n
         a(i) = vst + (i - 1)*d
      end do
!     write(0,*) '[makearray] a(1): ', a(1)
!     write(0,*) '[makearray] a(n): ', a(n)
      return
   end subroutine makearray

   subroutine sort(n, val, list)
      real(kind=8), dimension(n), intent(inout) :: val
      integer, dimension(n), intent(inout) :: list
      integer, intent(in) :: n
      real(kind=8), allocatable, dimension(:) :: work1
      integer, allocatable, dimension(:) :: work2
      integer :: i
      allocate(work1(n))
      allocate(work2(n))
      do i = 1, n
         list(i) = i
      end do
      call mergesort(n, val, list, work1, work2, 1, n)
      deallocate(work1)
      deallocate(work2)
      return
   end subroutine sort

   recursive subroutine mergesort(n, val, list, work1, work2, il, ir)
      real(kind=8), dimension(n), intent(inout) :: val, work1
      integer, dimension(n), intent(inout) :: list, work2
      integer, intent(in) :: n, il, ir
      integer :: i, j, k, mid
      if(il >= ir) return
      mid = (il + ir)/2
      call mergesort(n, val, list, work1, work2, il, mid)
      call mergesort(n, val, list, work1, work2, mid+1, ir)
      do i = il, mid
         work1(i) = val(i)
         work2(i) = list(i)
      end do
      do i = mid+1, ir
         j = mid + 1 + ir - i
         work1(i) = val(j)
         work2(i) = list(j)
      end do
      i = il
      j = ir
      do k = il, ir
         if(work1(i) < work1(j)) then
            val(k) = work1(i)
            list(k) = work2(i)
            i = i + 1
         else
            val(k) = work1(j)
            list(k) = work2(j)
            j = j - 1
         end if
      end do
      return
   end subroutine mergesort

   subroutine initialize()
      integer :: nl, k, kst
      real(kind=8) :: dummy0, dummy1, dummy2, dummy3

      open(1, file=trim(atmfile), form='formatted', status='old')

      nl = 0
      do
         read(1,*,end=100) dummy0, dummy1, dummy2, dummy3
          nl = nl + 1
      end do
100   continue

      nz = nl
      kst = 1
      if(with_Ocean) then
         nz = nl + 1
         kst = 2
      end if

      allocate(zl(nz))
      allocate(H(nz))
      allocate(alpha(nz))
      allocate(rho(nz))
      allocate(g(nz))
      allocate(ganma(nz))
      allocate(BVf(nz))
      allocate(lamda(nz))
      allocate(beta(nz))

      rewind(1)

      do k = kst, nz
         read(1,*,end=100) zl(k), H(k), alpha(k), rho(k)
      end do

      close(1)

      if(with_Ocean) then
         zl(1) = 0.0d0
         do k = kst, nz
            zl(k) = zl(k) + OH
         end do
         H(1) =OH
         rho(1) =Orho
         alpha(1) =Oalpha
      end if

      if(only_Ocean) then
         maxmode = 1
         minmode = 1
         nz = 1
      end if

!     write(0,'(a)') 'Making array t'
      call makearray(t, l_t, 0.0d0, tmax, dt)

!     write(0,'(a)') 'Making array vp'
      call makearray(vp, l_vp, 0.1d0, 0.4d0, 0.00001d0)

      Df = 2.0d0*pi/(dble(l_t)*dt)
!     write(0,'(a)') 'Making array omega'
      call makearray(omega, l_omega, Df, 0.05d0, Df)

      do k = 1, nz
         g(k) = g0*(Re/(Re+zl(k)+H(k)/2))**2
         ganma(k) = 1.0d0
      end do

      if(with_Ocean) then
         ganma(2:) = 1.4d0
         do k = 1, nz
            BVF(k) = g(k)*sqrt(ganma(k)-1.0d0)/alpha(k)
         end do
         BVF(1) = 0.003d0
      else
         ganma(:) = 1.4d0
         sind = 1
      end if

      if(only_Ocean) then
         sind = 1
      end if

      do k = 1, nz
         lamda(k) = ganma(k)*g(k)/(2.0d0*alpha(k)**2)
         beta(k) = 2.0d0*alpha(k)**2*BVf(k)/(g(k)*ganma(k))
      end do

      return
   end subroutine initialize

   subroutine finalize()
      deallocate(zl)
      deallocate(H)
      deallocate(alpha)
      deallocate(rho)
      deallocate(g)
      deallocate(ganma)
      deallocate(BVf)
      deallocate(lamda)
      deallocate(beta)

      deallocate(t)
      deallocate(vp)
      deallocate(omega)

      return
   end subroutine finalize

   subroutine calcDC()
      real(kind=8), allocatable, dimension(:,:) :: OMGmat, VPmat, Kmat
      complex(kind=8), allocatable, dimension(:,:) :: A22_
      real(kind=8), allocatable, dimension(:,:,:) :: DC
      real(kind=8), allocatable, dimension(:,:,:) :: DComg
      real(kind=8), allocatable, dimension(:,:,:) :: DCall
      integer :: l_omega_cmp, l_DCall
      integer :: i, j

      write(6,'(a)') 'Calculate dispersion curve'
      allocate(OMGmat(l_vp,l_omega))
      allocate(VPmat(l_vp,l_omega))
      allocate(Kmat(l_vp,l_omega))

      do j = 1, l_omega
         do i = 1, l_vp
            OMGmat(i,j) = omega(j)
            VPmat(i,j) = vp(i)
            Kmat(i,j) = OMGmat(i,j)/VPmat(i,j)
         end do
      end do

      call calcA22()

      call calcDC1()

      call calcDComg()

      call calcDCall()

      call fileout()

      deallocate(OMGmat)
      deallocate(VPmat)
      deallocate(Kmat)
      deallocate(A22_)
      deallocate(DC)
      deallocate(DComg)
      deallocate(DCall)

      return

      contains

      subroutine calcA22()
         complex(kind=8), allocatable, dimension(:,:,:,:,:) :: am 
         complex(kind=8), allocatable, dimension(:,:,:,:) :: A
         real(kind=8) :: delta
         complex(kind=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, a11, a12, a21, a22
         complex(kind=8) :: kr, Pm, Acheck
         integer :: i, j, l

         allocate(am(l_vp,l_omega,2,2,nz))
         allocate(A(l_vp,l_omega,2,2))
         allocate(A22_(l_vp,l_omega))

         do l = 1, nz
            do j = 1, l_omega
               do i = 1, l_vp
                  delta = g(l)**2*Kmat(i,j)**2 - OMGmat(i,j)**4
                  tmp1 = VPmat(i,j)**2/alpha(l)**2
                  tmp2 = VPmat(i,j)**2/beta(l)**2
                  kr = -im*sqrt(Kmat(i,j)**2*(1.0d0 - tmp1) - BVf(l)**2/VPmat(i,j)**2*(1.0d0 - tmp2))
                  Pm = kr*H(l)
                  tmp3 = g(l)/alpha(l)**2
                  tmp4 = alpha(l)**2/VPmat(i,j)**2 - ganma(l)/2.0d0
                  tmp5 = sin(Pm)/kr;
                  am(i,j,1,1,l) = exp(lamda(l)*H(l))*(cos(Pm) + tmp3*tmp4*tmp5);
                  am(i,j,2,2,l) = exp(-lamda(l)*H(l))*(cos(Pm) - tmp3*tmp4*tmp5);
                  tmp6 = g(l)**2*tmp4**2
                  tmp7 = alpha(l)**4*kr**2
                  tmp8 = rho(l)*alpha(l)**4*delta*kr
                  am(i,j,1,2,l) = im*(Kmat(i,j)**3*VPmat(i,j)**3)*(tmp6+tmp7)*sin(Pm)/tmp8
                  am(i,j,2,1,l) = im*rho(l)*delta*sin(Pm)/((Kmat(i,j)**3*VPmat(i,j)**3)*kr)
               end do
            end do
         end do

         A(:,:,:,:) = am(:,:,:,:,1)

         do l = 2, nz
            do j = 1, l_omega
               do i = 1, l_vp
                  a11 = A(i,j,1,1)
                  a12 = A(i,j,1,2)
                  a21 = A(i,j,2,1)
                  a22 = A(i,j,2,2)
                  A(i,j,1,1) = am(i,j,1,1,l)*a11 + am(i,j,1,2,l)*a21
                  A(i,j,2,1) = am(i,j,2,1,l)*a11 + am(i,j,2,2,l)*a21
                  A(i,j,1,2) = am(i,j,1,1,l)*a12 + am(i,j,1,2,l)*a22
                  A(i,j,2,2) = am(i,j,2,1,l)*a12 + am(i,j,2,2,l)*a22
               end do
            end do
         end do

         do j = 1, l_omega
            do i = 1, l_vp
#ifndef NONAN_NOFILT
               a11 = A(i,j,1,1)
               a12 = A(i,j,1,2)
               a21 = A(i,j,2,1)
               a22 = A(i,j,2,2)
               Acheck = a11*a22 - a12*a21
               if(abs(1.0d0 - Acheck) > 0.1d0) then
                  A(i,j,1,1) = NaN
                  A(i,j,1,2) = NaN
                  A(i,j,2,1) = NaN
                  A(i,j,2,2) = NaN
               end if
#endif
               A22_(i,j) = A(i,j,2,2)
            end do
         end do

         deallocate(am)
         deallocate(A)

         return
      end subroutine calcA22

      subroutine calcDC1()
         real(kind=8), allocatable, dimension(:) :: omgtmp, ksort, A22tmp, A22sort, omgsort
         integer, allocatable, dimension(:) :: ksi, zcind, DCind
         integer :: nmode, l_zcind
         real(kind=8) :: tmp, init_sign, old_sign
         integer :: i, j, l, ll

         allocate(DC(l_vp, maxmode, 3))
         allocate(DCind(maxmode))
         DC = NaN
         DCind = 1

         allocate(omgtmp(l_omega))
         allocate(ksort(l_omega))
         allocate(ksi(l_omega))
         allocate(A22tmp(l_omega))
         allocate(A22sort(l_omega))
         allocate(omgsort(l_omega))
         allocate(zcind(l_omega))

         nmode = 0
         do l = 1, l_vp
            do i = 1, l_omega
               ksort(i) = Kmat(l,i)
               omgtmp(i) = OMGmat(l,i)
            end do
            call sort(l_omega, ksort, ksi)
            do i = 1, l_omega
               tmp = dble(A22_(l,i))
               if(tmp == tmp) then
                  A22tmp(i) = sign(1.0d0, dble(A22_(l,i)))
               else
                  A22tmp(i) = NaN
               end if
            end do

            do i = 1, l_omega
               A22sort(i) = A22tmp(ksi(i))
               omgsort(i) = omgtmp(ksi(i))
            end do

            if(l == 1) then
               init_sign = sign(1.0d0, A22sort(1))
               old_sign = init_sign
            else
               init_sign = sign(1.0d0, A22sort(1))
               if(init_sign /= old_sign) then
                  nmode = nmode + 1
                  old_sign = init_sign
               end if
            end if

            l_zcind = 0
            do i = 1, l_omega-1
               if(A22sort(i)*A22sort(i+1) < 0.0d0) then
                  l_zcind = l_zcind + 1
                  zcind(l_zcind) = i
               end if
            end do

            do ll = 1, l_zcind
               if(nmode + ll > maxmode) then
                  exit
               end if
               DC(DCind(nmode+ll), nmode+ll, 1) = vp(l)
               DC(DCind(nmode+ll), nmode+ll, 2) = ksort(zcind(ll))
               DC(DCind(nmode+ll), nmode+ll, 3) = omgsort(zcind(ll))
               DCind(nmode+ll) = DCind(nmode+ll) + 1
            end do
         end do

         deallocate(omgtmp)
         deallocate(ksort)
         deallocate(ksi)
         deallocate(A22tmp)
         deallocate(A22sort)
         deallocate(omgsort)
         deallocate(zcind)

         deallocate(DCind)

         return
      end subroutine calcDC1

      subroutine calcDComg()
         real(kind=8), allocatable, dimension(:) :: vptmp, ksort, A22tmp, A22sort, vpsort
         integer, allocatable, dimension(:) :: ksi, zcind, DCind
         integer :: nmode, l_zcind
         real(kind=8) :: tmp, init_sign, old_sign
         integer :: i, j, l, ll

         l_omega_cmp = 0
         do l = 1, l_omega
            if(omega(l) > omg_searchmax) cycle
            l_omega_cmp = l_omega_cmp + 1
         end do

         allocate(DComg(l_omega_cmp, maxmode, 3))
               allocate(DCind(maxmode))
         DComg = NaN
         DCind = 1

         allocate(vptmp(l_vp))
         allocate(ksort(l_vp))
         allocate(ksi(l_vp))
         allocate(A22tmp(l_vp))
         allocate(A22sort(l_vp))
         allocate(vpsort(l_vp))
         allocate(zcind(l_vp))

         nmode = 0
         do l = 1, l_omega
            if(omega(l) > omg_searchmax) cycle
            do i = 1, l_vp
               vptmp(i) = VPmat(i,l)
               ksort(i) = Kmat(i,l)
            end do
            call sort(l_vp, ksort, ksi)
            do i = 1, l_vp
               tmp = dble(A22_(i,l))
               if(tmp == tmp) then
                  A22tmp(i) = sign(1.0d0, dble(A22_(i,l)))
               else
                  A22tmp(i) = NaN
               end if
            end do

            do i = 1, l_vp
               A22sort(i) = A22tmp(ksi(i))
               vpsort(i) = vptmp(ksi(i))
            end do

            if(l == 1) then
               init_sign = sign(1.0d0, A22sort(1))
               old_sign = init_sign
            else
               init_sign = sign(1.0d0, A22sort(1))
               if(init_sign /= old_sign) then
                  nmode = nmode + 1
                  old_sign = init_sign
               end if
            end if

            l_zcind = 0
            do i = 1, l_vp-1
               if(A22sort(i)*A22sort(i+1) < 0.0d0) then
                  l_zcind = l_zcind + 1
                  zcind(l_zcind) = i
               end if
            end do

            do ll = 1, l_zcind
               if(nmode + ll > maxmode) then
                  exit
               end if
               DComg(DCind(nmode+ll), nmode+ll, 1) = vpsort(zcind(ll))
               DComg(DCind(nmode+ll), nmode+ll, 2) = ksort(zcind(ll))
               DComg(DCind(nmode+ll), nmode+ll, 3) = omega(l)
               DCind(nmode+ll) = DCind(nmode+ll) + 1
            end do
         end do

         deallocate(vptmp)
         deallocate(ksort)
         deallocate(ksi)
         deallocate(A22tmp)
         deallocate(A22sort)
         deallocate(vpsort)
         deallocate(zcind)

         deallocate(DCind)

         return
      end subroutine calcDComg

      subroutine calcDCall()
         integer :: DCallind1, DCallind2
         real(kind=8) :: DCmin
         integer :: i, l

         l_DCall = l_vp
         do l = minmode, maxmode
            DCallind1 = 1
            DCmin = DC(1,l,3)
            do i = 2, l_vp
               if(DCmin > DC(i,l,3)) then
                  DCallind1 = i
                  DCmin = DC(i,l,3)
               end if
            end do

            DCallind2 = 1
            do i = 1, l_omega_cmp
               if(DCmin == DComg(i,maxmode-l+1,3)) DCallind2 = i
            end do

            if(DCallind2 /= 1) then
               l_DCall = max(l_DCall, DCallind1+DCallind2)
            end if
         end do

         allocate(DCall(l_DCall, maxmode, 3))
         DCall = NaN
         DCall(l_vp+1:l_DCall,:,:) = 0.0d0

         do l = minmode, maxmode
            DCallind1 = 1
            DCmin = DC(1,l,3)
            do i = 2, l_vp
               if(DCmin > DC(i,l,3)) then
                  DCallind1 = i
                  DCmin = DC(i,l,3)
               end if
            end do

            DCallind2 = 1
            do i = 1, l_omega_cmp
               if(DCmin == DComg(i,maxmode-l+1,3)) DCallind2 = i
            end do

            if(DCallind2 == 1) then
               DCall(1:DCallind1,l,:) = DC(1:DCallind1,l,:)
            else
               DCall(1:DCallind1,l,:) = DC(1:DCallind1,l,:)
               DCall(DCallind1+1:DCallind1+DCallind2,l,:) = DComg(DCallind2:1:-1,maxmode-l+1,:)
            end if
         end do

         return
      end subroutine calcDCall

      subroutine fileout()
         integer :: ll, l, i, j

         allocate(nkr(maxmode-minmode+1))
         nkr = 0

         ll = 1
         do l = minmode, maxmode
            j = 0
            do i = 1, l_DCall
               if(DCall(i,l,1) == DCall(i,l,1)) j = j + 1
            end do
            nkr(ll) = j
            ll = ll + 1
         end do

         l_out = sum(nkr)

         allocate(kroot(l_out))
         allocate(vpkr(l_out))
         allocate(omgkr(l_out))
         kroot = 0.0d0
         vpkr = 0.0d0
         omgkr = 0.0d0

         j = 0
         do i = 1, l_DCall
            if(DCall(i,minmode,1) == DCall(i,minmode,1)) then
               j = j + 1
               kroot(j) = DCall(i,minmode,2)
               vpkr(j) = DCall(i,minmode,1)
               omgkr(j) = DCall(i,minmode,3)
            end if
         end do

         do l = minmode+1, maxmode
            do i = 1, l_DCall
               if(DCall(i,l,1) == DCall(i,l,1)) then
                  j = j + 1
                  kroot(j) = DCall(i,l,2)
                  vpkr(j) = DCall(i,l,1)
                  omgkr(j) = DCall(i,l,3)
               end if
            end do
         end do

         if(DCsaveflag) then
            write(6,'(a)') 'Save dispersion curve'
            open(100, file=trim(DCfile), form='unformatted')
            write(100) nkr
            write(100) kroot
            write(100) vpkr
            write(100) omgkr
            close(100)

            open(100, file=trim(DCtext), form='formatted')

            write(100,'(a)') '----------------'
            write(100,'(a)') 'mode    nkr     '
            write(100,'(a)') '----------------'
            do l = minmode, maxmode
               write(100,'(2i8)') l, nkr(l-minmode+1)
            end do
!           write(100,'(a)') '----------------'

            write(100,'(a)') '-------------------------------------------------------------'
            write(100,'(a)') 'mode    nkr     kroot          vpkr           omgkr          '
            write(100,'(a)') '-------------------------------------------------------------'
            ll = 1
            do l = minmode, maxmode
               do i = 1, nkr(l-minmode+1)
                  write(100,'(2i8,3f15.8)') l, i, kroot(ll), vpkr(ll), omgkr(ll)
                  ll = ll + 1
               end do
            end do
!           write(100,'(a)') '-------------------------------------------------------------'

            close(100)
         end if

         return
      end subroutine fileout
   end subroutine calcDC

   subroutine readDC()
      allocate(nkr(maxmode-minmode+1))
      write(6,'(a)') 'Load dispersion curve'
      open(100, file=trim(DCfile), form='unformatted')
      read(100) nkr
      l_out = sum(nkr)
      allocate(kroot(l_out))
      allocate(vpkr(l_out))
      allocate(omgkr(l_out))
      read(100) kroot
      read(100) vpkr
      read(100) omgkr
      close(100)
      return
   end subroutine readDC

   subroutine unique(n_in, a_in, n_out, a_out, ind)
      integer, intent(in) :: n_in
      real(kind=8), dimension(n_in), intent(in) :: a_in
      integer, intent(out) :: n_out
      real(kind=8), allocatable, dimension(:), intent(out) :: a_out
      integer, allocatable, dimension(:), intent(out) :: ind
      real(kind=8), dimension(n_in) :: a_tmp
      integer, dimension(n_in) :: ind_tmp
      integer :: i, j
      do i = 1, n_in
         a_tmp(i) = a_in(i)
         ind_tmp(i) = i
      end do
      call sort(n_in, a_tmp, ind_tmp)
      n_out = 1
      do i = 2, n_in
         if(a_tmp(i-1) /= a_tmp(i)) n_out = n_out + 1
      end do
      allocate(a_out(n_out))
      allocate(ind(n_out))
      j = 1
      a_out(1) = a_tmp(1)
      ind(1) = ind_tmp(1)
      do i = 2, n_in
         if(a_tmp(i-1) /= a_tmp(i)) then
            j = j + 1
            a_out(j) = a_tmp(i)
            ind(j) = ind_tmp(i)
         else
            ind(j) = min(ind(j), ind_tmp(i))
         end if
      end do
      return
   end subroutine unique

   subroutine interp1(n_in, x_in, y_in, n_out, x_out, y_out)
      integer, intent(in) :: n_in, n_out
      real(kind=8), dimension(n_in), intent(in) :: x_in, y_in
      real(kind=8), dimension(n_out), intent(in) :: x_out
      real(kind=8), dimension(n_out), intent(out) :: y_out
      real(kind=8) :: x, xmin, xmax, r, x0, x1
      integer :: i, j
      do j = 1, n_out
         x = x_out(j)
         do i = 1, n_in-1
            xmin = x_in(i)
            xmax = x_in(i+1)
            if((x >= xmin) .and. (x < xmax)) then
               r = xmax - xmin
               x0 = (x - xmin)/r
               x1 = 1.0d0 - x0
               y_out(j) = y_in(i)*x1 + y_in(i+1)*x0
               exit
            end if
         end do
      end do
      return
   end subroutine interp1

   subroutine zgemm22(a, b, c)
      complex(kind=8), dimension(2,2), intent(in) :: a, b
      complex(kind=8), dimension(2,2), intent(out) :: c
      c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
      c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
      c(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
      c(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)
      return
   end subroutine zgemm22

   subroutine fft(n_in, a_in, n_out, a_out)
      include 'fftw3.f'
      integer, intent(in) :: n_in
      integer, intent(out) :: n_out
      real(kind=8), dimension(n_in), intent(in) :: a_in
      complex(kind=8), allocatable, dimension(:), intent(out) :: a_out
      integer(kind=8) :: plan
      integer :: i
      n_out = n_in
      allocate(a_out(n_in))
      call dfftw_plan_dft_1d(plan, n_out, a_out, a_out, FFTW_FORWARD, FFTW_MEASURE)
      do i = 1, n_out
         a_out(i) = dcmplx(a_in(i), 0.0d0)
      end do
      call dfftw_execute_dft(plan, a_out, a_out)
      call dfftw_destroy_plan(plan)
      return
   end subroutine fft

   subroutine ifft(n_in, a_in, n_out, a_out)
      include 'fftw3.f'
      integer, intent(in) :: n_in
      integer, intent(out) :: n_out
      complex(kind=8), dimension(n_in), intent(in) :: a_in
      complex(kind=8), allocatable, dimension(:), intent(out) :: a_out
      integer(kind=8) :: plan
      integer :: i
      n_out = n_in
      allocate(a_out(n_in))
      call dfftw_plan_dft_1d(plan, n_out, a_out, a_out, FFTW_BACKWARD, FFTW_MEASURE)
      do i = 1, n_out
         a_out(i) = a_in(i)
      end do
      call dfftw_execute_dft(plan, a_out, a_out)
      a_out = a_out/dble(n_out)
      call dfftw_destroy_plan(plan)
      return
   end subroutine ifft

   subroutine filter(b_data, a_data, n, x, zi, y)
      integer, intent(in) :: n
      real(kind=8), dimension(3), intent(inout) :: b_data, a_data
      real(kind=8), dimension(2), intent(in) :: zi
      real(kind=8), dimension(n), intent(in) :: x
      real(kind=8), dimension(n), intent(out) :: y
      real(kind=8) :: as
      integer :: i, j, a_tmp, naxpy, y_tmp
      if((abs(a_data(1) < huge(a_data(1)))) .and. (a_data(1) == a_data(1)) .and. (a_data(1) /= 0.0d0) .and. (a_data(1) /= 1.0d0)) then
         b_data(1) = b_data(1)/a_data(1)
         b_data(2) = b_data(2)/a_data(1)
         b_data(3) = b_data(3)/a_data(1)
         a_data(2) = a_data(2)/a_data(1)
         a_data(3) = a_data(3)/a_data(1)
         a_data(1) = 1.0d0
      end if
      y(1) = zi(1)
      y(2) = zi(2)
      do i = 3, n
         y(i) = 0.0d0
      end do
      do i = 1, n
         a_tmp = n - i + 1
         if(a_tmp <= 3) then
            naxpy = a_tmp
         else
            naxpy = 3
         end if
         do j = 0, naxpy - 1
            y_tmp = i + j
            y(y_tmp) = y(y_tmp) + x(i)*b_data(j + 1)
         end do
         if(a_tmp - 1 < 2) then
            naxpy = a_tmp - 2
         else
            naxpy = 1
         end if
         as = -y(i)
         do j = 0, naxpy
            y_tmp = i + j + 1
            y(y_tmp) = y(y_tmp) + as*a_data(j + 2)
         end do
      end do
      return
   end subroutine filter

   subroutine filtfilt(b_data, a_data, n, x)
      integer, intent(in) :: n
      real(kind=8), dimension(3), intent(in) :: b_data, a_data
      real(kind=8), dimension(n), intent(inout) :: x
      real(kind=8), dimension(2) :: zi
      real(kind=8), dimension(3) :: a, b
      real(kind=8), dimension(n) :: tmp
      integer :: i
      zi = 0.0d0
      b = b_data
      a = a_data
      call filter(b, a, n, x, zi, tmp)
      do i = 1, n
         x(i) = tmp(n-i+1)
      end do
      b = b_data
      a = a_data
      call filter(b, a, n, x, zi, tmp)
      do i = 1, n
         x(i) = tmp(n-i+1)
      end do
      return
   end subroutine filtfilt

   subroutine tukeywin(L, r, w)
      integer, intent(in) :: L
      real(kind=8), intent(in) :: r
      real(kind=8), dimension(L), intent(out) :: w
      real(kind=8) :: x, v
      integer :: i
      do i = 1, L
         x = dble(i-1)/dble(L-1)
         v = 1.0d0
         if(x < r/2.0d0) then
            v = (1.0d0 + cos(2.0d0*pi*(x - r/2.0d0)/r))/2.0d0
         end if
         if(x >= 1.0d0 - r/2.0d0) then
            v = (1.0d0 + cos(2.0d0*pi*(x - 1.0d0 + r/2.0d0)/r))/2.0d0
         end if
         w(i) = v
      end do
      return
   end subroutine tukeywin

   subroutine butter2h(Wn, a, b)
      real(kind=8), intent(in) :: wn
      real(kind=8), dimension(3), intent(out) :: a, b
      real(kind=8), parameter :: r2 = sqrt(2.0d0)
      real(kind=8) :: t
      t = tan(pi/2.0d0*Wn)
      a(1) = 1.0d0
      a(2) = (2.0d0*t*t - 2.0d0)/(t*t + r2*t + 1.0d0)
      a(3) = (t*t - r2*t + 1.0d0)/(t*t + r2*t + 1.0d0)
      b(1) = 1.0d0/(t*t + r2*t + 1.0d0)
      b(2) = -2.0d0/(t*t + r2*t + 1.0d0)
      b(3) = 1.0d0/(t*t + r2*t + 1.0d0)
      return
   end subroutine butter2h

   subroutine calcSW()
      real(kind=8), allocatable, dimension(:,:) :: baro
      real(kind=8), allocatable, dimension(:) :: FFTf
      complex(kind=8), allocatable, dimension(:) :: FFTbaro
      integer :: l_FFTf
      real(kind=8) :: Nf
      real(kind=8), allocatable, dimension(:) :: omgkrUni
      integer, allocatable, dimension(:) :: indUni
      integer :: jj, ll, dist, nmode, n_in, n_out, ind1, ind2, Find1, Find2, i, l_interp, l
      real(kind=8), allocatable, dimension(:) :: interp_tmp, kint, vpint, omgint
      complex(kind=8), allocatable, dimension(:,:,:,:) :: dam, damo, amkr
      complex(kind=8), allocatable, dimension(:,:,:) :: dA, dAo, Akr
      complex(kind=8), allocatable, dimension(:,:,:) :: As1
      complex(kind=8), allocatable, dimension(:) :: syntmp1, syntmp2
      complex(kind=8), allocatable, dimension(:) :: U
      real(kind=8), allocatable, dimension(:) :: Lt
      complex(kind=8), allocatable, dimension(:) :: P2W, P2Ppm, P2Pm
      complex(kind=8), allocatable, dimension(:) :: Lf, L_
      real(kind=8), allocatable, dimension(:) :: ttmp, taper
      complex(kind=8), allocatable, dimension(:) :: ffttmp
#ifndef JAGURS
      real(kind=8), allocatable, dimension(:,:) :: synwave1D, synwave1D2, synwave1D3, synwave1D4, synwave1D5
#else
      real(kind=8), allocatable, dimension(:,:) :: synwave1D
#endif
      integer :: l_Lf
      real(kind=8) :: kroot2, kroot3, vpkr2, vpkr3, omgkr2, omgkr3, omgkr4, kvpkr3, delta, theta, tmp
      complex(kind=8) :: ddelk, ddelo, dtmp1, kr, dkrk, dkro, Pm, dPmk, dPmo, dPkrk, dPkro
      complex(kind=8) :: dtmp2, dtmp3, dtmp4o, dtmp5, dtmp6, dtmp7k, dtmp8k, dtmp7o, dtmp8o, dtmp9o
      complex(kind=8) :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
      complex(kind=8), dimension(2,2) :: mt0, mt1
      real(kind=8) :: zr, zi, cyr, cyi, umax, umin, dmax, dmin, dtmp
      integer :: outnz, ierr, in, i1, i2, tind1, tind2, j
      real(kind=8), dimension(3) :: a_filt, b_filt
      real(kind=8) :: fs, fc, Wn
      logical :: flag_zero

      write(6,'(a)') 'Calcualte synthetic waveforms'

      nmode = maxmode - minmode + 1
      allocate(baro(l_t, nmode))

      Nf = 2.0d0*pi/(2.0d0*dt)
!     write(0,'(a)') 'Making array FFTf'
      call makearray(FFTf, l_FFTf, 0.0d0, Nf*2-Df, Df)
      allocate(FFTbaro(l_FFTf))

      Find1 = -1
      do i = 1, l_FFTf
         if(FFTf(i) - omega(1) >= 0.0d0) then
            Find1 = i
            exit
         end if
      end do

      Find2 = -1
      do i = 1, l_FFTf
         if(0.0124d0 - FFTf(i) >= 0.0d0) then
            Find2 = i
         end if
      end do

      allocate(synwave1D(dist_e-dist_s+1, l_t))
#ifndef JAGURS
      allocate(synwave1D2(dist_e-dist_s+1, l_t))
      allocate(synwave1D3(dist_e-dist_s+1, l_t))
      allocate(synwave1D4(dist_e-dist_s+1, l_t))
      allocate(synwave1D5(dist_e-dist_s+1, l_t))
#endif

      flag_zero = .false.
      jj = 1
      do dist = dist_s, dist_e
         if(dist == 0) then
            flag_zero = .true.
            jj = jj + 1
            cycle
         end if

         baro = 0.0d0
         do ll = 1, nmode
            if(nkr(ll) == 0) cycle
            FFTbaro = (0.0d0, 0.0d0)
            if(ll == 1) then
               n_in = nkr(1)
               call unique(n_in, omgkr, n_out, omgkrUni, indUni)
            else
               ind1 = sum(nkr(1:ll-1)) + 1
               ind2 = sum(nkr(1:ll-1)) + nkr(ll)
               n_in = ind2 - ind1 + 1
               call unique(n_in, omgkr(ind1), n_out, omgkrUni, indUni)
               indUni = indUni + ind1 - 1
            end if

            l_interp = Find2 - Find1 + 1
            allocate(kint(l_interp))
            allocate(vpint(l_interp))
            allocate(omgint(l_interp))

            allocate(interp_tmp(n_out))

            do i = 1, n_out
               interp_tmp(i) = kroot(indUni(i))
            end do
            call interp1(n_out, omgkrUni, interp_tmp, l_interp, FFTf(Find1), kint)

            do i = 1, n_out
               interp_tmp(i) = vpkr(indUni(i))
            end do
            call interp1(n_out, omgkrUni, interp_tmp, l_interp, FFTf(Find1), vpint)

            deallocate(interp_tmp)

            omgint = FFTf(Find1:Find2)

            allocate(dam(2,2,nz,l_interp))
            allocate(damo(2,2,nz,l_interp))
            allocate(amkr(2,2,nz,l_interp))
            dam = 0.0d0
            damo = 0.0d0
            amkr = 0.0d0

            do l = 1, nz
               do i = 1, l_interp
                  kroot2 = kint(i)**2
                  kroot3 = kint(i)**3
                  vpkr2 = vpint(i)**2
                  vpkr3 = vpint(i)**3
                  omgkr2 = omgint(i)**2
                  omgkr3 = omgint(i)**3
                  omgkr4 = omgint(i)**4
                  kvpkr3 = kroot3*vpkr3

                  delta = g(l)**2*kroot2 - omgkr4
                  ddelk = 2.0d0*g(l)**2*kint(i)
                  ddelo = -4.0d0*omgkr3

                  dtmp1 = kroot2 - omgkr2/alpha(l)**2 - BVf(l)**2*kroot2/omgkr2 + BVf(l)**2/beta(l)**2
                  kr = -im*sqrt(dtmp1)
                  dkrk = -im*(1.0d0 - BVf(l)**2/omgkr2)*kint(i)/sqrt(dtmp1)
                  dkro = -im*(BVf(l)**2*kroot2/omgkr3 - omgint(i)/alpha(l)**2)/sqrt(dtmp1)

                  Pm = kr*H(l)
                  dPmk = dkrk*H(l)
                  dPmo = dkro*H(l)

                  dPkrk = (dPmk*cos(Pm)*kr - dkrk*sin(Pm))/kr**2
                  dPkro = (dPmo*cos(Pm)*kr - dkro*sin(Pm))/kr**2

                  dtmp2 = 2*g(l)*sin(Pm)/kr
                  dtmp3 = g(l)/alpha(l)**2*(alpha(l)**2*kroot2/omgkr2 - ganma(l)/2.0d0)
                  dam(1,1,l,i) = exp(lamda(l)*H(l))*(-dPmk*sin(Pm) + dtmp2*kint(i)/omgkr2 + dtmp3*dPkrk)
                  dam(2,2,l,i) = exp(-lamda(l)*H(l))*(-dPmk*sin(Pm) - dtmp2*kint(i)/omgkr2 - dtmp3*dPkrk)
                  damo(1,1,l,i) = exp(lamda(l)*H(l))*(-dPmo*sin(Pm) - dtmp2*kroot2/omgkr3 + dtmp3*dPkro)
                  damo(2,2,l,i) = exp(-lamda(l)*H(l))*(-dPmo*sin(Pm) + dtmp2*kroot2/omgkr3 - dtmp3*dPkro)

                  dam(2,1,l,i) = im*rho(l)/omgkr3*(ddelk*sin(Pm)/kr + delta*dPkrk)

                  dtmp4o = (ddelo*omgkr3 - 3.0d0*delta*omgkr2)*sin(Pm)/(omgint(i)**6*kr)
                  damo(2,1,l,i) = im*rho(l)*(dtmp4o + delta*dPkro/omgkr3)

                  dtmp5 = alpha(l)**2*kroot2/omgkr2 - ganma(l)/2.0d0
                  dtmp6 = (g(l)**2*dtmp5**2 + alpha(l)**4*kr**2)/delta
                  dtmp7k = 4.0d0*alpha(l)**2*g(l)**2*dtmp5*kint(i)/omgkr2 + 2*alpha(l)**4*kr*dkrk
                  dtmp8k = (-ddelk*dtmp6 + dtmp7k)/delta
                  dam(1,2,l,i) = im*omgkr3/(rho(l)*alpha(l)**4)*(dtmp8k*sin(Pm)/kr +dtmp6*dPkrk)

                  dtmp7o = -4.0d0*alpha(l)**2*g(l)**2*dtmp5*kroot2/omgkr3 + 2.0d0*alpha(l)**4*kr*dkro
                  dtmp8o = (-ddelo*dtmp6 + dtmp7o)/delta
                  dtmp9o = dtmp8o*sin(Pm)/kr + dtmp6*dPkro
                  damo(1,2,l,i) = im/(rho(l)*alpha(l)**4)*(3.0d0*omgkr2*dtmp6*sin(Pm)/kr + omgkr3*dtmp9o)

                  tmp3 = g(l)/alpha(l)**2
                  tmp4 = alpha(l)**2/vpkr2 - ganma(l)/2.0d0
                  tmp5 = sin(Pm)/kr
                  amkr(1,1,l,i) = exp(lamda(l)*H(l))*(cos(Pm) + tmp3*tmp4*tmp5)
                  amkr(2,2,l,i) = exp(-lamda(l)*H(l))*(cos(Pm) - tmp3*tmp4*tmp5)

                  tmp6 = g(l)**2*tmp4**2
                  tmp7 = alpha(l)**4*kr**2
                  tmp8 = rho(l)*alpha(l)**4*delta*kr
                  amkr(1,2,l,i) = im*kvpkr3*(tmp6 + tmp7)*sin(Pm)/tmp8
                  amkr(2,1,l,i) = im*rho(l)*delta*sin(Pm)/(kvpkr3*kr)
               end do
            end do

            allocate(dA(2,2,l_interp))
            allocate(dAo(2,2,l_interp))
            allocate(Akr(2,2,l_interp))

            dA = dam(:,:,1,:) + amkr(:,:,1,:)
            dAo = damo(:,:,1,:) + amkr(:,:,1,:)
            Akr = amkr(:,:,1,:)

            do l = 2, nz
               do i = 1, l_interp
                  call zgemm22(dam (1,1,l,i), Akr(1,1,i), mt0)
                  call zgemm22(amkr(1,1,l,i), dA (1,1,i), mt1)
                  dA(:,:,i) = mt0 + mt1
               end do

               do i = 1, l_interp
                  call zgemm22(damo(1,1,l,i), Akr(1,1,i), mt0)
                  call zgemm22(amkr(1,1,l,i), dAo(1,1,i), mt1)
                  dAo(:,:,i) = mt0 + mt1
               end do

               do i = 1, l_interp
                  call zgemm22(amkr(1,1,l,i), Akr(1,1,i), mt0)
                  Akr(:,:,i) = mt0
               end do
            end do

            allocate(U(l_interp))

            do i = 1, l_interp
               U(i) = -dA(2,2,i)/dAo(2,2,i)
            end do

            allocate(Lt(l_t))
            do i = 1, l_t
               Lt(i) = Samp*sin(2*pi*t(i)/SDt)
            end do

            l = nint(SDt/dt)+1
            do i = l, l_t
               Lt(i) = -0.0d0
            end do

            call fft(l_t, Lt, l_Lf, Lf)

            allocate(L_(l_interp))

            do i = 1, l_interp
               L_(i) = Lf(i+Find1-1)
            end do

            allocate(As1(2,2,l_interp))
            As1(:,:,:) = amkr(:,:,1,:)
            if(sind > 2) then
               do l = 2, sind-1
                  do i = 1, l_interp
                     call zgemm22(amkr(1,1,l,i), As1(1,1,i), mt0)
                     As1(:,:,i) = mt0
                  end do
               end do
            end if
            allocate(syntmp1(l_interp))
            do i = 1, l_interp
               syntmp1(i) = -dconjg(Akr(2,1,i))*(dconjg(As1(2,2,i)) - im*g0*rho(sind)*dconjg(As1(1,2,i))/omgint(i))/dconjg(dA(2,2,i))
            end do
            theta = dist/Re

            allocate(syntmp2(l_interp))
            do i = 1, l_interp
               zr = dist*kint(i)
               zi = 0.0d0
               call zbesh(zr, zi, 0.0d0, 1, 2, 1, cyr, cyi, outnz, ierr)
               syntmp2(i) = sqrt(dble(dist/(Re*sin(theta))))*dcmplx(cyr, cyi)
            end do

            allocate(P2W(l_interp))
            allocate(P2Ppm(l_interp))
            allocate(P2Pm(l_interp))

            do i = 1, l_interp
               P2W(i) = -im*amkr(1,2,1,i)
               P2Ppm(i) = amkr(2,2,1,i)
               P2Pm(i) = P2Ppm(i) + rho(1)*g0*P2W(i)/omgint(i)
            end do

            do i = 1, l_interp
               FFTbaro(i+Find1-1) = P2Ppm(i)*L_(i)*syntmp1(i)*syntmp2(i)
            end do

            in = l_FFTf - ((l_FFTf-1)/2+2) + 1
            do i = 1, in
               i1 = i + ((l_FFTf-1)/2+2) -1
               i2 = -i + ((l_FFTf-1)/2+1) + 1
               FFTbaro(i1) = dconjg(FFTbaro(i2))
            end do

            umax = maxval(dble(U))
            allocate(ttmp(l_t))
            dmin = abs(t(1) - dble(dist)/(umax*1.2d0))
            ttmp(1) = dmin
            do i = 2, l_t
               dtmp = abs(t(i) - dble(dist)/(umax*1.2d0))
               ttmp(i) = dtmp
               dmin = min(dmin, dtmp)
            end do
            tind1 = 1
            do i = 1, l_t
               if(dmin == ttmp(i)) then
                  tind1 = i
                  exit
               end if
            end do
            deallocate(ttmp)

            umin = minval(dble(U))
            allocate(ttmp(l_t))
            dmin = abs(t(1) - dble(dist)/(umin*0.8d0))
            ttmp(1) = dmin
            do i = 2, l_t
               dtmp = abs(t(i) - dble(dist)/(umin*0.8d0))
               ttmp(i) = dtmp
               dmin = min(dmin, dtmp)
            end do
            tind2 = 1
            do i = 1, l_t
               if(dmin == ttmp(i)) then
                  tind2 = i
                  exit
               end if
            end do
            deallocate(ttmp)

            allocate(taper(l_t))
            taper = 0.0d0
            if(tind2 < l_t) then
               j = tind2 - tind1 + 1
               allocate(ttmp(j))
               call tukeywin(j, 0.1d0, ttmp)
               do i = 1, j
                  taper(i+tind1-1) = ttmp(i)
               end do
               deallocate(ttmp)
            else
               j = l_t - tind1 + 1
               allocate(ttmp(j))
               call tukeywin(j, 0.1d0, ttmp)
               do i = 1, j
                  taper(i+tind1-1) = ttmp(i)
               end do
               deallocate(ttmp)
            end if
            call ifft(l_FFTf, FFTbaro, i, ffttmp)
            do i = 1, l_t
!              baro(i,ll) = taper(i)*dble(ffttmp(i))
               baro(i,ll) = dble(ffttmp(i))
            end do

#ifndef NONAN_NOFILT
            fs = 1.0d0/dt
            fc = 1.0d0/5000.0d0
            Wn = fc/(fs/2.0d0)
            call butter2h(Wn, a_filt, b_filt)
            call filtfilt(b_filt, a_filt, l_t, baro(1,ll))
#endif

            deallocate(ffttmp)

            deallocate(taper)

            deallocate(P2W)
            deallocate(P2Ppm)
            deallocate(P2Pm)

            deallocate(syntmp1)
            deallocate(syntmp2)
            deallocate(As1)

            deallocate(Lt)
            deallocate(Lf)
            deallocate(L_)

            deallocate(U)

            deallocate(dA)
            deallocate(dAo)
            deallocate(Akr)

            deallocate(dam)
            deallocate(damo)
            deallocate(amkr)

            deallocate(kint)
            deallocate(vpint)
            deallocate(omgint)
            deallocate(omgkrUni)
            deallocate(indUni)
         end do

         do i = 1, l_t
            tmp = 0.0d0
#ifndef JAGURS
            do j = 3, nmode
               tmp = tmp + baro(i,j)
            end do
            synwave1D(jj,i) = tmp
            synwave1D2(jj,i) = baro(i,3)
            synwave1D3(jj,i) = baro(i,4)
            synwave1D4(jj,i) = baro(i,5)
            synwave1D5(jj,i) = baro(i,6)
#else
            do j = 1, nmode
               tmp = tmp + baro(i,j)
            end do
            synwave1D(jj,i) = tmp
#endif
         end do
         jj = jj +1;
      end do

      if(flag_zero) then
         do i = 1, l_t
#ifndef JAGURS
            synwave1D(1,i) = synwave1D(2,i)
            synwave1D2(1,i) = synwave1D2(2,i)
            synwave1D3(1,i) = synwave1D3(2,i)
            synwave1D4(1,i) = synwave1D4(2,i)
            synwave1D5(1,i) = synwave1D5(2,i)
#else
            synwave1D(1,i) = synwave1D(2,i)
#endif
         end do
      end if

#ifndef JAGURS
      open(100, file=trim(outfile)//'all.dat', form='unformatted')
      write(100) 100.0d0*synwave1D
      close(100)

      open(100, file=trim(outfile)//'0.dat', form='unformatted')
      write(100) 100.0d0*synwave1D5
      close(100)

      open(100, file=trim(outfile)//'1.dat', form='unformatted')
      write(100) 100.0d0*synwave1D4
      close(100)

      open(100, file=trim(outfile)//'2.dat', form='unformatted')
      write(100) 100.0d0*synwave1D3
      close(100)

      open(100, file=trim(outfile)//'3.dat', form='unformatted')
      write(100) 100.0d0*synwave1D2
      close(100)
#else
      open(100, file='normalmode.dat', form='unformatted')
      write(100) 100.0d0*synwave1D
      close(100)
#endif

      deallocate(synwave1D)
#ifndef JAGURS
      deallocate(synwave1D2)
      deallocate(synwave1D3)
      deallocate(synwave1D4)
      deallocate(synwave1D5)
#endif

      deallocate(baro)
      deallocate(FFTf)
      deallocate(FFTbaro)

      return
   endsubroutine calcSW
end module mod_PHcalc_wave6

program test
use mod_PHcalc_wave6
implicit none
#ifndef JAGURS
namelist /normalmode/ srclon, srclat, dist_s, dist_e, tmax, dt, sind, Samp, SDt, DCflag, SYNflag, OH, atmfile
#else
namelist /normalmode/ srclon, srclat, dist_s, dist_e, tmax, dt, sind, minmode, maxmode, Samp, SDt, DCflag, SYNflag, OH, atmfile
#endif

open(1,file='normalmode.namelist',action='read',status='old',form='formatted')
read(1,normalmode)
close(1)

write(DCfile, '(a,i0,a,i0,a)') 'DC_PH5_', minmode, '-', maxmode, '.dat'
write(DCtext, '(a,i0,a,i0,a)') 'DC_PH5_', minmode, '-', maxmode, '.txt'
#ifndef JAGURS
write(outfile, '(a,i0,a,i0,a,i0,a)') 'atm_NM_US_', dist_s, '-', dist_e, 'km_', int(tmax), 'sec_PH5_GR'
#endif

call initialize()
if(nz < sind) then
   write(0,*) 'ERROR! sind must be equal or less than nz!'
   write(0,*) 'nz:   ', nz
   write(0,*) 'sind: ', sind
   stop
end if
if(sind < 1) then
   write(0,*) 'ERROR! sind must be greater than 0!'
   write(0,*) 'sind: ', sind
   stop
end if

if(DCflag) then
   call calcDC()
else
   call readDC()
end if

if(SYNflag) call calcSW()

call finalize()

stop
end program test
