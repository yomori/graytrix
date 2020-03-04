ccc graytrix_v2.2.f   Time-stamp: <2017-05-18 11:17:03 hamana>
!
! This is the fortran program of the full-sky gravitational lensing ray-tracing
! in which the light ray path and magnification matrix are evaluated by the 
! multiple lens-plane algorithm following the CALCLENS (Becker 2012, see also Teyssier et al 2009).
! I thank Becker very much for making the source codes of CALCLENS publicity available.
! The lensing magnification matrix is evaluated by using the recurrence relation (e.g., Hilbert et al 2009).
!
!* required library:
!  Healpix (graytrix_v2.2.f was developed and tested with Healpix_3.11)
!  cfitsio
!
!* compile with gfortran: 
! gfortran graytrix_v2.2.f -o graytrix.exe -O2 -I/usr/local/Healpix_3.11/include -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -L/usr/local/Healpix_3.11/lib -L/usr/local/lib -lhealpix -lhpxgif -lcfitsio -Wl,-R/usr/local/lib
!
!!! usage:
! ./graytrix.exe [options] < input_parameter_file
! [options] are:
! -nres : Healpix's parameter "nres"
! -d    : path to output dir (detauld = 'DATA')
! -f    : output file name (default = 'test')
! -om   : Omega_M (default = 0.279)
! -wde  : dark energy EoS parameter "w_DE" (default = -1.0)
! -th   : thickness of shell [Mpc/h] (default = 150.0)
! -ip   : inverse distance weight power for interpolation (default = 1.0)
! -nt   : number of openmp threads (default = 16)
!
! contents of input_parameter_file should be 
! 4 - number of source planes [integer > 0]
! 0.6 - source plane redshift for 1st source plane 
! 1.0 - source plane redshift for 2nd source plane
! 1.4 - source plane redshift for 3rd source plane
! 1.8 - source plane redshift for 4th source plane (more lines follow for more source plane outputs)
! /work2/cfca2048/PSIS/nres12/b450r000p001 - file name for input potential data (b450r000p001.psis) for nearest plane
! /work2/cfca2048/PSIS/nres12/b450r000p002 - file name for input potential data (b450r000p002.psis) for 2nd plane
! /work2/cfca2048/PSIS/nres12/b450r000p003 
! /work2/cfca2048/PSIS/nres12/b900r000p004
! /work2/cfca2048/PSIS/nres12/b900r000p005
! /work2/cfca2048/PSIS/nres12/b900r000p006
! /work2/cfca2048/PSIS/nres12/b1350r000p007
! /work2/cfca2048/PSIS/nres12/b1350r000p008
! /work2/cfca2048/PSIS/nres12/b1350r000p009
! /work2/cfca2048/PSIS/nres12/b1800r000p010
! /work2/cfca2048/PSIS/nres12/b1800r000p011
! /work2/cfca2048/PSIS/nres12/b1800r000p012
! /work2/cfca2048/PSIS/nres12/b2250r000p013
! /work2/cfca2048/PSIS/nres12/b2250r000p014
! /work2/cfca2048/PSIS/nres12/b2250r000p015
! /work2/cfca2048/PSIS/nres12/b2700r000p016
! /work2/cfca2048/PSIS/nres12/b2700r000p017
! /work2/cfca2048/PSIS/nres12/b2700r000p018
! /work2/cfca2048/PSIS/nres12/b3150r000p019
! /work2/cfca2048/PSIS/nres12/b3150r000p020
! /work2/cfca2048/PSIS/nres12/b3150r000p021
! /work2/cfca2048/PSIS/nres12/b3600r000p022
! /work2/cfca2048/PSIS/nres12/b3600r000p023
! /work2/cfca2048/PSIS/nres12/b3600r000p024 - file name for input potential data (b3600r000p024.psis) for farthest plane
!
!!! example:
!
! ./graytrix.exe -nres 12 -f output_file_name < graytrix.in
!
!!! copyright --- Takashi Hamana, Ryuichi Takahashi (2016,2017)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE cosmoparam
      real omegam
      real wde
      END MODULE cosmoparam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program graytrix
      use alm_tools
      use healpix_types
      use pix_tools
      use udgrade_nr
      use fitstools, only : output_map
      use head_fits 
      use cosmoparam

      implicit real (a-h,o-z)
      implicit integer (i-n)
      
      external dldz
      
!!! variables for fitsio      
      integer(4) :: ft_status,ft_unit,ft_blocksize,ft_bitpix
      integer(4) :: ft_naxis,ft_naxes(2)
      integer(4) :: ft_group,ft_fpixel,ft_nelements
      integer(4) :: ft_pcount,ft_gcount
      character(len=128) :: ft_comment
      logical :: ft_simple,ft_extend,ft_anyf

!!! valiables for Healpix
      integer(i4b) :: nres,nside,lmax,hporder,hpspin
      integer(i4b) :: nsidein,lmaxin
      integer(i8b) :: npix,ipring,ipnest,ipix,kpix,iray
      integer(i8b) :: npixin
      real(sp),allocatable :: hpmap(:,:) !hpmap(0:npix-1,3) output
      real(sp),allocatable :: hpmap_alp(:,:)!(1:2,0:npix-1) input potentials
      real(sp),allocatable :: hpmap_mag(:,:)!(1:3,0:npix-1) input potentials
      character(len=80), dimension(1:64) :: mapheader

      integer(i8b),allocatable :: neighlist(:,:) ! (8,0:npix-1)
      integer(i8b) :: nlist(8)
      integer(i4b) :: nneigh

!!! valiables for Ray-tracing
      real(sp),allocatable :: ampmat(:,:,:,:) !ampmat(2,2,3,0:npix-1) magnification matrix for three planes
      real(sp),allocatable :: umat(:,:) !(3,0:npix-1) phi,ij at ray positions
      real(sp) :: ampmat1(2,2) ! tmp data

      real(sp),allocatable :: conv(:),shear(:,:),omeg(:) !(0:npix-1,1/2) output
      real(sp),allocatable :: dflang(:,:) !(0:npix-1,1:2) output

      real(dp),allocatable :: theta_i(:),phi_i(:) !(0:npix-1) angular ray position on image-plane
      real(dp),allocatable :: theta_s(:),phi_s(:) !(0:npix-1) output - angular ray position on source-plane
      real(dp) :: theta_ray,phi_ray ! angular position of the ray 
      real(dp),allocatable :: xyzray(:,:) ! 3dim coordinate of ray (x,y,z) = xyzray(3,0:npix-1)
      real(dp) :: dvec(3)
      real(sp) :: xbase(3)
      real(dp),allocatable :: dvecpix(:,:) ! dvec of pixel dvecpix(3,0:npix-1)
      real(sp),allocatable :: betaray(:,:) ! beta of ray betaray(3,0:npix-1)
      real(sp) :: betatmp(3)

      real(sp) :: alphanorm,xxyy,tpnorm
      real(sp),allocatable :: phihat(:,:),thetahat(:,:) ! (3,0:npix-1)
      real(sp),allocatable :: alpxyz(:,:) !(3,0:npix-1))
      real(sp),allocatable :: alpxyzpix(:,:) ! alpha of pixels in xyz alpxyzpix(3,0:npix-1)
      real(sp) :: rotvec(3),rotvecnorm,rotmat(3,3)
      
      real(sp) :: dlambda,b4lambda,c4lambda,dxyzray

      real(dp) :: deltaphi,cosalp,sinalp,alpha,cosdelta,sindelta
      real(dp) :: cosoversin
      real(dp) :: cossep(8),cosray,sinray
      real(dp) :: cosmax(4),wgt(4),wgtnorm,dray
      integer(i4b) :: i4neigh(4)
      
!!! source redshifts
      real,allocatable :: dists(:)
      character(len=128) :: snum

      character(len=128) :: fname,fnamei
      character(len=128) :: outfile,outdir,fin,filein,outlog

      character*128 opt,arg

      include 'omp_lib.h'

!
!!!!!! definitions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      chubble0=2.9979e3 ! c/H0 in Mpc/h (H0=100h km/sec/Mpc)

!!! Healpix related
      hporder=1 ! =1 for ring, =2 for nested

!
!!!!!!! set up parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      narg=iargc()
c      write(*,*)narg
!!! cosmological parameters *** FLAT COSMOLOGY ONLY ***
      omegam=0.279
      wde=-1.0

!!! Healpix related
      nres=13

!!! output file names
      outdir='DATA'
      fname='test'

!!! inverse distance weight power 
!!! f(x) = sum_i w_i f_i(x_i) / sum_i w_i
!!! with w_i = 1/d^pidw
!!! where d is the distance between x and x_i
      pidw=4.0
      
!!! RayTrace related: thickness of shells
      thick=150.0 !Mpc/h

!!! number of threds
      nopenmp=32

      do i=1,narg
         call getarg(i,opt)
         call getarg(i+1,arg)
         select case (opt)
         case ('-nres')
            read(arg,*)nres  ! Healpix parameter
         case ('-d')
            outdir=arg
         case ('-f')
            fname=arg
         case ('-om')
            read(arg,*)omegam
         case ('-wde')
            read(arg,*)wde
         case ('-th')
            read(arg,*)thick
         case ('-ip')
            read(arg,*)pidw
         case ('-nt')
            read(arg,*)nopenmp
         end select
      enddo
      pidw=0.5*pidw ! a numerical trick: d^-pidw = (d*d)^(-0.5*pidw)

!!! set parameters for openmp
!$    call omp_set_num_threads(nopenmp)

!!! Healpix related
      nside=2**nres
      npix=12_i8b*nside*nside
      lmax=3*nside
      write(*,*)'##### Healpix parameters:'
      write(*,*)'# nres  =',nres
      write(*,*)'# nside =',nside
      write(*,*)'# npix  =',npix
      write(*,*)'# lmax  =',lmax

!!! setup neighbours list
      allocate(neighlist(8,0:npix-1))
!$omp parallel private(ipnest,nlist,nneigh,j)
!$omp do
      do ipring=0,npix-1
         do j=1,8
            neighlist(j,ipring)=-1
         enddo
         call ring2nest(nside,ipring,ipnest)
         call neighbours_nest(nside,ipnest,nlist,nneigh)
         do j=1,nneigh
            call nest2ring(nside,nlist(j),neighlist(j,ipring))
         enddo
      enddo
!$omp end do
!$omp end parallel

!!! RayTrace related ::: set initial values
      allocate(ampmat(2,2,3,0:npix-1))
      allocate(umat(3,0:npix-1))
      allocate(theta_i(0:npix-1))
      allocate(phi_i(0:npix-1))
      allocate(betaray(3,0:npix-1))
      allocate(dvecpix(3,0:npix-1))
      allocate(xyzray(3,0:npix-1))
      allocate(phihat(3,0:npix-1))
      allocate(thetahat(3,0:npix-1))
      allocate(alpxyzpix(3,0:npix-1))
      allocate(alpxyz(3,0:npix-1))

!!! log
      outlog=ADJUSTL(ADJUSTR(fname)//'.log')
      outlog=ADJUSTL(ADJUSTR(outdir)//'/'//outlog)
      open(20,file=outlog,status='unknown')
      write(20,*)'outdir = ',outdir(1:len_trim(outdir))
      write(20,*)'fname  = ',fname(1:len_trim(fname))
      write(20,*)'nres  =',nres
      write(20,*)'nside =',nside
      write(20,*)'npix  =',npix      
      write(20,*)'omegam =',omegam
      write(20,*)'wde    =',wde
      write(20,*)'thick =',thick
      write(20,*)'pidw =',pidw
      write(20,*)'source redshifts:'      
! source redshifts
      read(*,*)nsplanes
      write(20,*)'n_zs =',nsplanes      
      allocate(dists(nsplanes))
      do izs=1,nsplanes
         read(*,*)zs
         write(20,*)'zs [',izs,'] =',zs
         azs=1.0/(zs+1.0)
         call qromb(dldz,azs,1.0,ss)
         dists(izs)=ss*chubble0 ! distance to source plane
         write(20,*)'ds [',izs,'] =',dists(izs)
      enddo
      write(20,*)'input files:'

!!! Healpix 
      allocate(hpmap_alp(1:2,0:npix-1))
      allocate(hpmap_mag(1:3,0:npix-1))
      
!
!!! start ray-tracing 
!
      koutflag=0
      isplanes=1
      kplane=0
!
!!!! initial setting
!
c      distl=0.5*thick ! distance to the 1st lens-plane
      distlmin=0.0
      distlmax=thick
      distl=0.75*(distlmax**4-distlmin**4)/(distlmax**3-distlmin**3) ! light-cone weighted mean distance 
                                                                     ! change made in v1.1
!$omp parallel private(dvec,xxyy,tpnorm,j)
!$omp do      
      do ipix=0,npix-1
         call pix2vec_ring(nside,ipix,dvec)
         do j=1,3
            xyzray(j,ipix)=distl*dvec(j)
            betaray(j,ipix)=dvec(j)
            dvecpix(j,ipix)=dvec(j)
         enddo
ccc   theta-phi unit vector
         xxyy=dvec(1)*dvec(1)+dvec(2)*dvec(2)
         tpnorm=sqrt(xxyy)
         thetahat(3,ipix)=-tpnorm
         tpnorm=1.0/tpnorm
         phihat(1,ipix)=-dvec(2)*tpnorm
         phihat(2,ipix)=dvec(1)*tpnorm
         phihat(3,ipix)=0.0
         tpnorm=dvec(3)*tpnorm
         thetahat(1,ipix)=dvec(1)*tpnorm
         thetahat(2,ipix)=dvec(2)*tpnorm
         call vec2ang(dvec,theta_i(ipix),phi_i(ipix)) ! image position
      enddo
!$omp end do
!$omp end parallel

!$omp parallel private(j)
!$omp do 
      do ipix=0,npix-1
         do j=1,2
            ampmat(1,1,j,ipix)=1.0 ! A(1,1)
            ampmat(2,1,j,ipix)=0.0 ! A(2,1)
            ampmat(1,2,j,ipix)=0.0 ! A(1,2)
            ampmat(2,2,j,ipix)=1.0 ! A(2,2)
         enddo
      enddo
!$omp end do
!$omp end parallel

!
!!!!! loop for lens planes
!     
 111  read(*,'(a128)',end=999)fin ! loop for lens planes
      kplane=kplane+1
      write(*,*)'#**********',kplane,'lens-plane **********'      
c      distl=(kplane-0.5)*thick ! distance to this lens-plane
      distlmin=(kplane-1)*thick
      distlmax=kplane*thick
      distl=0.75*(distlmax**4-distlmin**4)/(distlmax**3-distlmin**3) ! change made in v1.1
      distlinv=1.0/distl
      write(*,*)'# distance to this lens-plane =',distl
c      distllast=distl-thick ! distance to the last lens-plane      
      distlmin=(kplane-2)*thick
      distlmax=(kplane-1)*thick
      distllast=0.75*(distlmax**4-distlmin**4)
     $     /(distlmax**3-distlmin**3) ! change made in v1.1
c      distlnext=(kplane+0.5)*thick ! distance to the next lens-plane
      distlmin=kplane*thick
      distlmax=(kplane+1)*thick
      distlnext=0.75*(distlmax**4-distlmin**4)
     $     /(distlmax**3-distlmin**3) ! change made in v1.1
      if (distlnext.ge.dists(isplanes)) then
         koutflag=1
c         distlnext=dists
      endif

!
!!!!!! reading psis file
!
      filein=ADJUSTL(ADJUSTR(fin))
      write(*,*)'##### reading psis data: ',filein(1:len_trim(filein))
      !filein=ADJUSTL(ADJUSTR(fin)//'.psis')
      write(20,*)kplane,filein(1:len_trim(filein))
      open(10,file=filein,status='old',form='unformatted')
      read(10)nsidein,npixin
      write(*,*)'# nsidein = ',nsidein
      write(*,*)'# npixin = ',npixin
      read(10)hpmap_alp
      read(10)hpmap_mag
      close(10)

!
!!! covariant derivatives 
!!! see Appendic A of Becker (2013) MNRAS, 435, 115
!
!$omp parallel private(cosoversin)
!$omp do
      do ipix=0,npix-1
         cosoversin=dcos(theta_i(ipix))/dsin(theta_i(ipix))
         hpmap_mag(2,ipix)=hpmap_mag(2,ipix)
     $        -cosoversin*hpmap_alp(2,ipix)
         hpmap_mag(3,ipix)=hpmap_mag(3,ipix)
     $        +cosoversin*hpmap_alp(1,ipix)
      enddo
!$omp end do
!$omp end parallel

ccc alpha in 3-vector
!$omp parallel private(j)
!$omp do
      do ipix=0,npix-1
         do j=1,3
            alpxyzpix(j,ipix)=hpmap_alp(1,ipix)*thetahat(j,ipix)
     $           +hpmap_alp(2,ipix)*phihat(j,ipix)
         enddo
      enddo
!$omp end do
!$omp end parallel
!
!!!!!!! for each ray, searching for 4 nearest pixels from the ray 
!!!!!!! and interpolating the alpha and magnification matrix
!
      write(*,*)'##### interpolating psi,i & psi,ij at ray positions'
!
      if (kplane.eq.1) then ! for the first lens plane (image plane) no interpolation is needed 
!$omp parallel private(j)
!$omp do
         do ipix=0,npix-1 
            umat(1,ipix)=-hpmap_mag(1,ipix) ! umat(1)=phi_11=kappa+gamma1
            umat(2,ipix)=-hpmap_mag(2,ipix) ! umat(2)=phi_12=phi_12=-gamma2
            umat(3,ipix)=-hpmap_mag(3,ipix) ! umat(3)=phi_22=kappa-gamma1
            do j=1,3
               alpxyz(j,ipix)=alpxyzpix(j,ipix)
            enddo
         enddo
!$omp end do
!$omp end parallel

      else ! for i>1 plane, serching for neighbor pixels and doing interpolation 

c clear array
!$omp parallel
!$omp do
         do ipix=0,npix-1
            do j=1,3
               umat(j,ipix)=0.0
               alpxyz(j,ipix)=0.0
            enddo
         enddo
!$omp end do
!$omp end parallel
!
!$omp parallel
!$omp+private(dvec,ipring,i4neigh,theta_ray,phi_ray,cosray,sinray,dray,
!$omp+cosmax,wgt,kpix,nlist,cossep,nneigh,wgtnorm,iray,ampmat1,
!$omp+j,k,l,m,n)
!$omp do
         do ipix=0,npix-1
            do j=1,3 
               dvec(j)=xyzray(j,ipix)*distlinv
            enddo
            call vec2pix_ring(nside,dvec,ipring) !!! pixel nearest to the ray
            i4neigh(1)=ipring   ! nearest pixel
            call vec2ang(dvec,theta_ray,phi_ray) !!! (theta,phi) of ray
            cosray=dcos(theta_ray)
            sinray=dsin(theta_ray)
            cosmax(1)=dcos(theta_i(ipring))*cosray
     $           +dsin(theta_i(ipring))*sinray
     $           *dcos(phi_i(ipring)-phi_ray)
            if (cosmax(1).gt.1.0_dp) then
               dray=0.0_dp
            else
               dray=dacos(cosmax(1))
            endif
            wgt(1)=(dray*dray)**pidw

!!! searching for next nearest 3 pixels            
            cosmax(2)=-1.1_dp  ! searching for 2nd nearest pixel 
            do k=1,7 ! loop for surrounding 8 pixels (1-7)
               kpix=neighlist(k,ipring)
               nlist(k)=kpix
               cossep(k)=dcos(theta_i(kpix))*cosray
     $              +dsin(theta_i(kpix))*sinray
     $              *dcos(phi_i(kpix)-phi_ray)
               if (cossep(k).gt.cosmax(2)) then
                  cosmax(2)=cossep(k)
                  i4neigh(2)=kpix ! 2nd nearest pixel
               endif
            enddo
            nneigh=7
            if (neighlist(8,ipring).ge.0) then ! for the surrounding 8th pixel, if exits
               nneigh=8
               kpix=neighlist(8,ipring)
               nlist(8)=kpix
               cossep(8)=dcos(theta_i(kpix))*cosray
     $              +dsin(theta_i(kpix))*sinray
     $              *dcos(phi_i(kpix)-phi_ray)
               if (cossep(8).gt.cosmax(2)) then
                  cosmax(2)=cossep(8)
                  i4neigh(2)=kpix ! 2nd nearest pixel
               endif
            endif
            if (cosmax(2).gt.1.0_dp) then
               dray=0.0_dp
            else
               dray=dacos(cosmax(2))
            endif
            wgt(2)=(dray*dray)**pidw
c     
            cosmax(3)=-1.1_dp  ! searching for 3rd nearest pixel
            do k=1,nneigh
               if (nlist(k).ne.i4neigh(2).and.
     $              cossep(k).gt.cosmax(3)) then
                  cosmax(3)=cossep(k)
                  i4neigh(3)=nlist(k) ! 3rd nearest pixel
               endif
            enddo
            if (cosmax(3).gt.1.0_dp) then
               dray=0.0_dp
            else
               dray=dacos(cosmax(3))
            endif
            wgt(3)=(dray*dray)**pidw
c     
            cosmax(4)=-1.1_dp ! searching for 4th nearest pixel
            do k=1,nneigh
               if (nlist(k).ne.i4neigh(2).and.nlist(k).ne.i4neigh(3).
     $              and.cossep(k).gt.cosmax(4)) then
                  cosmax(4)=cossep(k)
                  i4neigh(4)=nlist(k) ! 4th nearest pixel
               endif
            enddo
            if (cosmax(4).gt.1.0_dp) then
               dray=0.0_dp
            else
               dray=dacos(cosmax(4))
            endif
            wgt(4)=(dray*dray)**pidw
c     
            wgtnorm=1.0_dp
            do l=2,4
               wgt(l)=wgt(1)/wgt(l)
               wgtnorm=wgtnorm+wgt(l)
            enddo
            wgt(1)=1.0_dp
            do l=1,4
               wgt(l)=wgt(l)/wgtnorm
            enddo
ccc transport the magnification matrix from ray's xbase() to the image position's xbase() 
ccc and computing alpxyz at ray-position from 4 neighbors with the weight wgt
            do l=1,4
               iray=i4neigh(l)
c for magnification matrix
               ampmat1(1,1)=-hpmap_mag(1,iray)
               ampmat1(2,1)=-hpmap_mag(2,iray)
               ampmat1(1,2)=-hpmap_mag(2,iray)
               ampmat1(2,2)=-hpmap_mag(3,iray)
               if (wgt(l).lt.1.0_dp) 
     $              call transampmat(dvecpix(1:3,iray),
     $              dvecpix(1:3,ipix),ampmat1)
               umat(1,ipix)=umat(1,ipix)+wgt(l)*ampmat1(1,1)
               umat(2,ipix)=umat(2,ipix)+wgt(l)*
     $              0.5*(ampmat1(1,2)+ampmat1(2,1))
               umat(3,ipix)=umat(3,ipix)+wgt(l)*ampmat1(2,2)
c for deflection angle in xyz coordinate (thus no base transformation is needed) 
               do j=1,3
                  alpxyz(j,ipix)=alpxyz(j,ipix)+wgt(l)*alpxyzpix(j,iray)
               enddo
            enddo
         enddo
!$omp end do
!$omp end parallel
      endif

      if (koutflag.eq.1) goto 511
!     
!!!!!!!!! computing kappa & gamma on the next place
!!!!!!!!! with multiple-lens plane recurrence relation !!!!
!
 115  distl2distlnext=distlnext/distl
      write(*,*)'##### computing kappa and gamma for the plane =',kplane
      dratinrec2=distl*(distlnext-distllast)
     $     /(distlnext*(distl-distllast))
      dratinrec1=1.0-dratinrec2
      dratinrec3=-(distlnext-distl)/distlnext
!$omp parallel private(m,n)
!$omp do
      do ipix=0,npix-1
c computing A(i,j) on this source plane
         ampmat(1,1,3,ipix)=dratinrec1*ampmat(1,1,1,ipix) ! A11
     $        +dratinrec2*ampmat(1,1,2,ipix)
     $        +dratinrec3*(umat(1,ipix)*ampmat(1,1,2,ipix)
     $        +umat(2,ipix)*ampmat(2,1,2,ipix))
         ampmat(2,1,3,ipix)=dratinrec1*ampmat(2,1,1,ipix) ! A21
     $        +dratinrec2*ampmat(2,1,2,ipix)
     $        +dratinrec3*(umat(2,ipix)*ampmat(1,1,2,ipix)
     $        +umat(3,ipix)*ampmat(2,1,2,ipix))
         ampmat(1,2,3,ipix)=dratinrec1*ampmat(1,2,1,ipix) ! A12
     $        +dratinrec2*ampmat(1,2,2,ipix)
     $        +dratinrec3*(umat(1,ipix)*ampmat(1,2,2,ipix)+
     $        umat(2,ipix)*ampmat(2,2,2,ipix))
         ampmat(2,2,3,ipix)=dratinrec1*ampmat(2,2,1,ipix) ! A22
     $        +dratinrec2*ampmat(2,2,2,ipix)
     $        +dratinrec3*(umat(2,ipix)*ampmat(1,2,2,ipix)+
     $        umat(3,ipix)*ampmat(2,2,2,ipix))
c updating
         do n=1,2
            do m=1,2
               ampmat(m,n,1,ipix)=ampmat(m,n,2,ipix)
               ampmat(m,n,2,ipix)=ampmat(m,n,3,ipix)
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel
!     
!!!!!!!!! computing ray positions on the next lens plane 
!!!!!!!!! using the ray position and deflection angle at this plane  
!
      write(*,*)'##### ray-position at the next plane'
      c4lambda=distl*distl-distlnext*distlnext
!$omp parallel 
!$omp+private(alphanorm,xbase,rotvec,rotvecnorm,rotmat,betatmp,
!$omp+b4lambda,dlambda,dxyzray,dvec,j,k)
!$omp do
      do ipix=0,npix-1
         alphanorm=sqrt(alpxyz(1,ipix)*alpxyz(1,ipix)
     $        +alpxyz(2,ipix)*alpxyz(2,ipix)
     $        +alpxyz(3,ipix)*alpxyz(3,ipix))
         if (alphanorm.gt.0.0) then
            do j=1,3
               xbase(j)=xyzray(j,ipix)/distl
            enddo
ccc   compute xbase x alpha
            rotvec(1)=xbase(2)*alpxyz(3,ipix)-xbase(3)*alpxyz(2,ipix)
            rotvec(2)=xbase(3)*alpxyz(1,ipix)-xbase(1)*alpxyz(3,ipix)
            rotvec(3)=xbase(1)*alpxyz(2,ipix)-xbase(2)*alpxyz(1,ipix)
            rotvecnorm=1.0/sqrt(rotvec(1)*rotvec(1)
     $           +rotvec(2)*rotvec(2)+rotvec(3)*rotvec(3))
            do j=1,3
               rotvec(j)=rotvec(j)*rotvecnorm
            enddo
ccc   compute rotation matrix
            call gen_rotmat(rotvec,alphanorm,rotmat)
ccc   compute beta for the next plane
            do j=1,3
               betatmp(j)=0.0
               do k=1,3
                  betatmp(j)=betatmp(j)+rotmat(j,k)*betaray(k,ipix)
               enddo
            enddo
         else
            do j=1,3
               betatmp(j)=betaray(j,ipix)  ! RT
            enddo
         endif
ccc   compute lambda
         b4lambda=0.0
         do j=1,3 
            b4lambda=b4lambda+xyzray(j,ipix)*betatmp(j)
         enddo
         dlambda=sqrt(b4lambda*b4lambda-c4lambda)-b4lambda
         if (dlambda.lt.0.0) dlambda=c4lambda/dlambda
         dxyzray=0.0
         do j=1,3
            xyzray(j,ipix)=xyzray(j,ipix)+dlambda*betatmp(j)
            dxyzray=dxyzray+xyzray(j,ipix)*xyzray(j,ipix)
            betaray(j,ipix)=betatmp(j)
         enddo
         dxyzray=1.0/sqrt(dxyzray)
         do j=1,3
            dvec(j)=xyzray(j,ipix)*dxyzray ! normalized position
            xyzray(j,ipix)=dvec(j)*distlnext
         enddo
      enddo
!$omp end do
!$omp end parallel
      goto 111 ! loop for lens planes
      close(20)

!     
!!!!!!!!! computing kappa & gamma on the next place
!!!!!!!!! with multiple-lens plane recurrence relation !!!!
!
 511  write(*,*)'##### output for the source plane =',isplanes
      write(*,*)'# redshift of this source plane = ',zs
      write(*,*)'# distance to this source plane = ',dists(isplanes)
      write(snum,*)isplanes
      fnamei=ADJUSTL(ADJUSTR(fname)//'.zs')
      fnamei=ADJUSTL(ADJUSTR(fnamei)//ADJUSTL(snum))
      write(*,*)'# file name = ',fnamei(1:len_trim(fnamei))
      write(*,*)'# computing kappa & gamma'
      distl2distlnext=dists(isplanes)/distl      
      dratinrec2=distl*(dists(isplanes)-distllast)
     $     /(dists(isplanes)*(distl-distllast))
      dratinrec1=1.0-dratinrec2
      dratinrec3=-(dists(isplanes)-distl)/dists(isplanes)
!$omp parallel private(m,n)
!$omp do
      do ipix=0,npix-1
c computing A(i,j) on this source plane
         ampmat(1,1,3,ipix)=dratinrec1*ampmat(1,1,1,ipix) ! A11
     $        +dratinrec2*ampmat(1,1,2,ipix)
     $        +dratinrec3*(umat(1,ipix)*ampmat(1,1,2,ipix)
     $        +umat(2,ipix)*ampmat(2,1,2,ipix))
         ampmat(2,1,3,ipix)=dratinrec1*ampmat(2,1,1,ipix) ! A21
     $        +dratinrec2*ampmat(2,1,2,ipix)
     $        +dratinrec3*(umat(2,ipix)*ampmat(1,1,2,ipix)
     $        +umat(3,ipix)*ampmat(2,1,2,ipix))
         ampmat(1,2,3,ipix)=dratinrec1*ampmat(1,2,1,ipix) ! A12
     $        +dratinrec2*ampmat(1,2,2,ipix)
     $        +dratinrec3*(umat(1,ipix)*ampmat(1,2,2,ipix)+
     $        umat(2,ipix)*ampmat(2,2,2,ipix))
         ampmat(2,2,3,ipix)=dratinrec1*ampmat(2,2,1,ipix) ! A22
     $        +dratinrec2*ampmat(2,2,2,ipix)
     $        +dratinrec3*(umat(2,ipix)*ampmat(1,2,2,ipix)+
     $        umat(3,ipix)*ampmat(2,2,2,ipix))
      enddo
!$omp end do
!$omp end parallel

cccc output
      allocate(conv(0:npix-1))
      allocate(shear(0:npix-1,1:2))
      allocate(omeg(0:npix-1))
!$omp parallel
!$omp do
      do ipix=0,npix-1
         conv(ipix)=1.0-0.5*(ampmat(1,1,3,ipix)+ampmat(2,2,3,ipix))
         shear(ipix,1)=0.5*(ampmat(2,2,3,ipix)-ampmat(1,1,3,ipix))
         shear(ipix,2)=-0.5*(ampmat(1,2,3,ipix)+ampmat(2,1,3,ipix))
         omeg(ipix)=0.5*(ampmat(2,1,3,ipix)-ampmat(1,2,3,ipix))
      enddo
!$omp end do
!$omp end parallel
       
!!!! output kappa & gamma map
      outfile=ADJUSTL(ADJUSTR(fnamei)//'.mag.dat')
      outfile=ADJUSTL(ADJUSTR(outdir)//'/'//outfile)
      open(10,file=outfile,status='unknown',form='unformatted')
      write(10)nside,npix
      write(10)conv
      write(10)(shear(ipix,1),ipix=0,npix-1)
      write(10)(shear(ipix,2),ipix=0,npix-1)
      write(10)omeg
      close(10)
      deallocate(omeg)
      
!!!! output fits 
      allocate(hpmap(0:npix-1,3))
!$omp parallel
!$omp do
      do ipix=0,npix-1 
         hpmap(ipix,1)=conv(ipix)
         hpmap(ipix,2)=shear(ipix,1)
         hpmap(ipix,3)=shear(ipix,2)
      enddo
!$omp end do
!$omp end parallel
      deallocate(conv)
      deallocate(shear)      

!!! check if fitsfile exists, and if it does, delete it 
      outfile=ADJUSTL(ADJUSTR(fnamei)//'.mag.fits')
      outfile=ADJUSTL(ADJUSTR(outdir)//'/'//outfile)
      ft_status=0
      ft_blocksize=1
      call ftgiou(ft_unit,ft_status)
      call ftopen(ft_unit,outfile,1,ft_blocksize,ft_status)
      if (ft_status.eq.0) then 
         call ftdelt(ft_unit,ft_status)
      elseif (ft_status.eq.104)then 
         ft_status=0
         call ftcmsg
      else ! there was some other error opening the file; delete the file anyway
         ft_status=0
         call ftcmsg
         call ftdelt(ft_unit,ft_status)
      endif
      call ftclos(ft_unit,ft_status)
      call ftfiou(ft_unit,ft_status)

!!! output fits
      call write_minimal_header(mapheader, 'MAP', nside=nside,
     $     order=hporder)
      call add_card(mapheader,'TTYPE1','kappa','label for field 1')
      call add_card(mapheader,'TTYPE2','shear1','label for field 2')
      call add_card(mapheader,'TTYPE3','shear2','label for field 3')

      call output_map(hpmap,mapheader,outfile)
      deallocate(hpmap)

!     
!!!!!!!!! computing ray positions on the next lens plane 
!!!!!!!!! using the ray position and deflection angle at this plane  
!
      write(*,*)'##### ray-position at the source plane'
      allocate(theta_s(0:npix-1))
      allocate(phi_s(0:npix-1))
      c4lambda=distl*distl-dists(isplanes)*dists(isplanes)
!$omp parallel 
!$omp+private(alphanorm,xbase,rotvec,rotvecnorm,rotmat,betatmp,
!$omp+b4lambda,dlambda,dxyzray,dvec,j,k)
!$omp do
      do ipix=0,npix-1
         alphanorm=sqrt(alpxyz(1,ipix)*alpxyz(1,ipix)
     $        +alpxyz(2,ipix)*alpxyz(2,ipix)
     $        +alpxyz(3,ipix)*alpxyz(3,ipix))
         if (alphanorm.gt.0.0) then
            do j=1,3
               xbase(j)=xyzray(j,ipix)/distl
            enddo
ccc   compute xbase x alpha
            rotvec(1)=xbase(2)*alpxyz(3,ipix)-xbase(3)*alpxyz(2,ipix)
            rotvec(2)=xbase(3)*alpxyz(1,ipix)-xbase(1)*alpxyz(3,ipix)
            rotvec(3)=xbase(1)*alpxyz(2,ipix)-xbase(2)*alpxyz(1,ipix)
            rotvecnorm=1.0/sqrt(rotvec(1)*rotvec(1)
     $           +rotvec(2)*rotvec(2)+rotvec(3)*rotvec(3))
            do j=1,3
               rotvec(j)=rotvec(j)*rotvecnorm
            enddo
ccc   compute rotation matrix
            call gen_rotmat(rotvec,alphanorm,rotmat)
ccc   compute beta for the next plane
            do j=1,3
               betatmp(j)=0.0
               do k=1,3
                  betatmp(j)=betatmp(j)+rotmat(j,k)*betaray(k,ipix)
               enddo
            enddo
         else
            do j=1,3
               betatmp(j)=betaray(j,ipix)  ! RT
            enddo
         endif
ccc   compute lambda
         b4lambda=0.0
         do j=1,3 
            b4lambda=b4lambda+xyzray(j,ipix)*betatmp(j)
         enddo
         dlambda=sqrt(b4lambda*b4lambda-c4lambda)-b4lambda
         if (dlambda.lt.0.0) dlambda=c4lambda/dlambda
         dxyzray=0.0
         do j=1,3
            dvec(j)=xyzray(j,ipix)+dlambda*betatmp(j)
            dxyzray=dxyzray+dvec(j)*dvec(j)
         enddo
         dxyzray=1.0/sqrt(dxyzray)
         do j=1,3
            dvec(j)=dvec(j)*dxyzray ! normalized position
         enddo
         call vec2ang(dvec,theta_s(ipix),phi_s(ipix))
      enddo
!$omp end do
!$omp end parallel

!!!! output source position
      outfile=ADJUSTL(ADJUSTR(fnamei)//'.src.dat')
      outfile=ADJUSTL(ADJUSTR(outdir)//'/'//outfile)
      open(10,file=outfile,status='unknown',form='unformatted')
      write(10)nside,npix
      write(10)(theta_s(ipix),ipix=0,npix-1)   ! RT
      write(10)(phi_s(ipix),ipix=0,npix-1)
      close(10)
      deallocate(theta_s)
      deallocate(phi_s)
      
      koutflag=0
      isplanes=isplanes+1
      if (isplanes.gt.nsplanes) goto 999
      goto 115 ! return to the loop for the lens planes

c!!!! compute deflection angle
c      allocate(dflang(0:npix-1,1:2))
c!$omp parallel private(deltaphi,cosalp,alpha,sinalp,cosdelta,sindelta)
c!$omp do
c      do ipix=0,npix-1
c         deltaphi=phi_s(ipix)-phi_i(ipix)
c         cosalp=dcos(theta_i(ipix))*dcos(theta_s(ipix))
c     $        +dsin(theta_i(ipix))*dsin(theta_s(ipix))*dcos(deltaphi)
c         if (cosalp.gt.1.0_dp) then
c            alpha=0.0_dp
c         else
c            alpha=dacos(cosalp)
c         endif
c         sinalp=dsin(alpha)
c         if (sinalp.eq.0.0_dp) then
c            dflang(ipix,1)=0.0
c            dflang(ipix,2)=0.0
c         else
c            cosdelta=(dcos(theta_s(ipix))-dcos(theta_i(ipix))*cosalp)
c     $           /(dsin(theta_i(ipix))*sinalp)
c            sindelta=dsin(deltaphi)*dsin(theta_s(ipix))/sinalp
c            dflang(ipix,1)=-alpha*cosdelta
c            dflang(ipix,2)=alpha*sindelta
c         endif
c      enddo
c!$omp end do
c!$omp end parallel
c!!!! output deflection angle
c      outfile=ADJUSTL(ADJUSTR(fname)//'.alp.dat')
c      outfile=ADJUSTL(ADJUSTR(outdir)//'/'//outfile)
c      open(10,file=outfile,status='unknown',form='unformatted')
c      write(10)nside,npix
c      write(10)(dflang(ipix,1),ipix=0,npix-1)
c      write(10)(dflang(ipix,2),ipix=0,npix-1)
c      close(10)


!!! ray-tracing done, closing arrays
 999  deallocate(neighlist)
      deallocate(betaray)
      deallocate(dvecpix)
      deallocate(phihat)
      deallocate(thetahat)
      deallocate(alpxyzpix)
      deallocate(alpxyz)
      deallocate(umat)
      deallocate(hpmap_alp)
      deallocate(hpmap_mag)
      deallocate(xyzray)
      deallocate(theta_i)
      deallocate(phi_i)

      end program graytrix

ccc
      subroutine transampmat(xbase,dvec,ampmat)
ccc This is the FORTRAN version of the paratrans_tangtensor of CALCLENS (Becker 2012).
ccc I thank Becker very much for making the source codes of CALCLENS publicity available.
      use healpix_types
      implicit real (a-h,o-z)
      implicit integer (i-n)
      real(dp) :: xbase(3),dvec(3) !!!! should be normalized !!!
      real(sp) :: ampmat(2,2),t1(2,2)
      real(dp) :: axis(3),cosangle,sinangle,p(3)
      real(dp) :: ephi_dvec(3),etheta_dvec(3),rephi_xbase(3)
      real(dp) :: norm,sinpsi,cospsi,r(2,2),rt(2,2)
      
      axis(1)=xbase(2)*dvec(3)-xbase(3)*dvec(2)
      axis(2)=xbase(3)*dvec(1)-xbase(1)*dvec(3)
      axis(3)=xbase(1)*dvec(2)-xbase(2)*dvec(1)
      cosangle=xbase(1)*dvec(1)+xbase(2)*dvec(2)+xbase(3)*dvec(3)
      sinangle=dsqrt(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
      if (sinangle.ne.0.0_dp) then
         forall (j=1:3) axis(j)=axis(j)/sinangle
      else
         axis(1)=1.0_dp
         axis(2)=0.0_dp
         axis(3)=0.0_dp
      endif

      p(1)=-xbase(2)
      p(2)=xbase(1)
      p(3)=0.0_dp
      call gen_crotmax(p,rephi_xbase,axis,cosangle,sinangle)
      
      ephi_dvec(1)=-dvec(2)
      ephi_dvec(2)=dvec(1)
      ephi_dvec(3)=0.0_dp
      
      etheta_dvec(1)=dvec(3)*dvec(1)
      etheta_dvec(2)=dvec(3)*dvec(2)
      etheta_dvec(3)=-1.0_dp*(dvec(1)*dvec(1)+dvec(2)*dvec(2))

      if (dvec(3).gt.1.0_dp) dvec(3)=1.0_dp
      if (xbase(3).gt.1.0_dp) xbase(3)=1.0_dp
      norm=dsqrt((1.0_dp-dvec(3)*dvec(3))*(1.0_dp-xbase(3)*xbase(3)))

      sinpsi=(rephi_xbase(1)*etheta_dvec(1)
     $     +rephi_xbase(2)*etheta_dvec(2)
     $     +rephi_xbase(3)*etheta_dvec(3))/norm
      cospsi=(rephi_xbase(1)*ephi_dvec(1)
     $     +rephi_xbase(2)*ephi_dvec(2)
     $     +rephi_xbase(3)*ephi_dvec(3))/norm
      
      r(1,1)=cospsi
      r(1,2)=-1.0_dp*sinpsi
      r(2,1)=sinpsi
      r(2,2)=cospsi
      rt(1,1)=r(1,1)
      rt(1,2)=r(2,1)
      rt(2,1)=r(1,2)
      rt(2,2)=r(2,2)
      
      do i=1,2
         do j=1,2
            t1(i,j)=ampmat(i,1)*r(1,j)+ampmat(i,2)*r(2,j)
         enddo
      enddo
      do i=1,2
         do j=1,2
            ampmat(i,j)=rt(i,1)*t1(1,j)+rt(i,2)*t1(2,j)
         enddo
      enddo

      return
      end

      subroutine gen_crotmax(vec,rvec,axis,cosangle,sinangle)
      use healpix_types
      implicit real (a-h,o-z)
      implicit integer (i-n)
      real(dp) :: vec(3),rvec(3),axis(3),cosangle,sinangle
      real(dp) :: axisdotvec,axiscrossvec(3)
      axisdotvec=axis(1)*vec(1)+axis(2)*vec(2)+axis(3)*vec(3)
      axiscrossvec(1)=axis(2)*vec(3)-axis(3)*vec(2)
      axiscrossvec(2)=axis(3)*vec(1)-axis(1)*vec(3)
      axiscrossvec(3)=axis(1)*vec(2)-axis(2)*vec(1)
      rvec(1)=vec(1)*cosangle+axis(1)*axisdotvec*(1.0_dp-cosangle)
     $     +axiscrossvec(1)*sinangle
      rvec(2)=vec(2)*cosangle+axis(2)*axisdotvec*(1.0_dp-cosangle)
     $     +axiscrossvec(2)*sinangle
      rvec(3)=vec(3)*cosangle+axis(3)*axisdotvec*(1.0_dp-cosangle)
     $     +axiscrossvec(3)*sinangle
      return
      end

ccc   
      SUBROUTINE gen_rotmat(v,t,r)
      use healpix_types
      implicit real (a-h,o-z)
      implicit integer (i-n)
      real(sp) :: v(3),t,r(3,3)
      real(sp) :: c,s,c1
      real(sp) :: v1,v2,v3,s1,s2,s3
      c=cos(t)
      s=sin(t)
      c1=1.0-c
      v1=v(1)*c1
      v2=v(2)*c1
      v3=v(3)*c1
      s1=v(1)*s
      s2=v(2)*s
      s3=v(3)*s
      r(1,1)=v(1)*v1+c
      r(2,2)=v(2)*v2+c
      r(3,3)=v(3)*v3+c
      r(1,2)=v(1)*v2-s3
      r(2,1)=v(2)*v1+s3
      r(1,3)=v(1)*v3+s2
      r(3,1)=v(3)*v1-s2
      r(2,3)=v(2)*v3-s1
      r(3,2)=v(3)*v2+s1
      return
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
ccc  Newton-Raphson method
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newlap(alm,aout,ain)
      implicit real (a-h,o-z)
      implicit integer (i-n)
      parameter (eps=1.0e-5)
      external dldz
      rhs=alm
      a=ain
 356  call qromb(dldz,a,1.0,solv)
      err=abs(solv/rhs-1.0)
      if (err.gt.eps) then
         dela=-(rhs-solv)/dldz(a)
         a=a+dela
         goto 356
      endif
      aout=a
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
ccc  package of solving the integral equation of 
ccc  affine parameter - redshit (scale factor) relation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dldz(a)
      use cosmoparam
      implicit real (a-h,o-z)
      implicit integer (i-n)      
      dldz=1.0/sqrt(a*omegam+(1.0-omegam)*a**(1.0-3.0*wde))
      return
      end
ccc
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      write(*,*)'too many steps in qromb'
      stop
      END
ccc
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
ccc
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*)'failure in polint'
             stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

