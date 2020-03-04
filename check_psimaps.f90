!
! This is the fortran program of the full-sky gravitational lensing
! ray-tracing
! in which the light ray path and magnification matrix are evaluated by
! the 
! multiple lens-plane algorithm following the CALCLENS (Becker 2012, see
! also Teyssier et al 2009).
! I thank Becker very much for making the source codes of CALCLENS
! publicity available.
! The lensing magnification matrix is evaluated by using the recurrence
! relation (e.g., Hilbert et al 2009).
!
!* required library:
!  Healpix (graytrix_v2.2.f was developed and tested with Healpix_3.11)
!  cfitsio
!
!* compile with gfortran: 
! gfortran graytrix_v2.2.f -o graytrix.exe -O2
! -I/usr/local/Healpix_3.11/include -DGFORTRAN -fno-second-underscore
! -fopenmp -fPIC -L/usr/local/Healpix_3.11/lib -L/usr/local/lib
! -lhealpix -lhpxgif -lcfitsio -Wl,-R/usr/local/lib
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
! -ip   : inverse distance weight power for interpolation (default =
! 1.0)
! -nt   : number of openmp threads (default = 16)
!
! contents of input_parameter_file should be 
! 4 - number of source planes [integer > 0]
! 0.6 - source plane redshift for 1st source plane 
! 1.0 - source plane redshift for 2nd source plane
! 1.4 - source plane redshift for 3rd source plane
! 1.8 - source plane redshift for 4th source plane (more lines follow
! for more source plane outputs)
! /work2/cfca2048/PSIS/nres12/b450r000p001 - file name for input
! potential data (b450r000p001.psis) for nearest plane
! /work2/cfca2048/PSIS/nres12/b450r000p002 - file name for input
! potential data (b450r000p002.psis) for 2nd plane
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
! /work2/cfca2048/PSIS/nres12/b3600r000p024 - file name for input
! potential data (b3600r000p024.psis) for farthest plane
!
!!! example:
!
! ./graytrix.exe -nres 12 -f output_file_name < graytrix.in
!
!!! copyright --- Takashi Hamana, Ryuichi Takahashi (2016,2017)
! 

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
      character(len=1024) :: path,dir_out,infile,outfile
!!! valiables for Healpix
      integer(i4b) :: nres,nside,lmax,hporder,hpspin,narg
      integer(i4b) :: nsidein,lmaxin
      integer(i8b) :: npix,ipring,ipnest,ipix,kpix,iray
      integer(i8b) :: npixin
      real(sp),allocatable :: hpmap(:,:) !hpmap(0:npix-1,3) output
      real(sp),allocatable :: hpmap_alp(:,:)!(1:2,0:npix-1) input
      real(sp),allocatable :: hpmap_mag(:,:)!(1:3,0:npix-1) input
      character(len=80), dimension(1:64) :: mapheader

      integer(i8b),allocatable :: neighlist(:,:) ! (8,0:npix-1)
      integer(i8b) :: nlist(8)
      integer(i4b) :: nneigh

!!! valiables for Ray-tracing
      real(sp),allocatable :: ampmat(:,:,:,:) !ampmat(2,2,3,0:npix-1)
      real(sp),allocatable :: umat(:,:) !(3,0:npix-1) phi,ij at ray
      real(sp) :: ampmat1(2,2) ! tmp data

      real(sp),allocatable :: conv(:),shear(:,:),omeg(:) !(0:npix-1,1/2)
      real(sp),allocatable :: dflang(:,:) !(0:npix-1,1:2) output

      real(dp),allocatable :: theta_i(:),phi_i(:) !(0:npix-1) angular
      real(dp),allocatable :: theta_s(:),phi_s(:) !(0:npix-1) output -
      real(dp) :: theta_ray,phi_ray ! angular position of the ray 
      real(dp),allocatable :: xyzray(:,:) ! 3dim coordinate of ray
      real(dp) :: dvec(3)
      real(sp) :: xbase(3)
      real(dp),allocatable :: dvecpix(:,:) ! dvec of pixel
      real(sp),allocatable :: betaray(:,:) ! beta of ray
      real(sp) :: betatmp(3)

      real(sp) :: alphanorm,xxyy,tpnorm
      real(sp),allocatable :: phihat(:,:),thetahat(:,:) ! (3,0:npix-1)
      real(sp),allocatable :: alpxyz(:,:) !(3,0:npix-1))
      real(sp),allocatable :: alpxyzpix(:,:) ! alpha of pixels in xyz
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
      character(len=128) :: outdir,fin,filein,outlog

      character*128 opt,arg


      nside=2**13
      npix=12_i8b*nside*nside

narg=iargc()
write(*,*)narg

do i=1,narg
    call getarg(i,opt)
    call getarg(i+1,arg)
    select case (opt)
    case ('-infile')    ! input prefix
        read (arg,'(A)') infile
    case ('-outfile')
        read (arg,'(A)') outfile      ! output prefix
    end select
enddo


      allocate(hpmap_alp(1:2,0:npix-1))
      allocate(hpmap_mag(1:3,0:npix-1))
      
      open(10,file=infile,status='old',form='unformatted')
      read(10)nsidein,npixin
      read(10)hpmap_alp
      read(10)hpmap_mag
      close(10)


      allocate(hpmap(0:npix-1,1:1))
!$omp parallel private(j)
!$omp do
      do ipix=0,npix-1
         do j=1,1
            hpmap(ipix,j)=hpmap_mag(j,ipix)
         enddo
      enddo
!$omp end do
!$omp end parallel
      call write_minimal_header(mapheader,'MAP',nside=8192,ordering='RING')
      call add_card(mapheader,'TTYPE1','convergence','label for field 1')
      call output_map(hpmap,mapheader,outfile)
      deallocate(hpmap)

!      allocate(hpmap(0:npix-1,1:1))
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,1
!            hpmap(ipix,j)=hpmap_mag(2,ipix)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
!      call write_minimal_header(mapheader,'MAP',nside=8192,ordering='RING')
!      call add_card(mapheader,'TTYPE1','convergence','label for field1')
!      call output_map(hpmap,mapheader,'hpmap_mag2.fits')
!      deallocate(hpmap)


!      allocate(hpmap(0:npix-1,1:1))
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,1
!            hpmap(ipix,j)=hpmap_mag(3,ipix)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
!      call write_minimal_header(mapheader,'MAP',nside=8192,ordering='RING')
!      call add_card(mapheader,'TTYPE1','convergence','label for field1')
!      call output_map(hpmap,mapheader,'hpmap_mag3.fits')
!      deallocate(hpmap)

!      allocate(hpmap(0:npix-1,1:1))
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,1
!            hpmap(ipix,j)=hpmap_alp(1,ipix)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
!      call write_minimal_header(mapheader,'MAP',nside=8192,ordering='RING')
!      call add_card(mapheader,'TTYPE1','convergence','label for field1')
!      call output_map(hpmap,mapheader,'hpmap_alp1.fits')
!      deallocate(hpmap)

!      allocate(hpmap(0:npix-1,1:1))
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,1
!            hpmap(ipix,j)=hpmap_alp(2,ipix)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
!      call write_minimal_header(mapheader,'MAP',nside=8192,ordering='RING')
!      call add_card(mapheader,'TTYPE1','convergence','label for field1')
!      call output_map(hpmap,mapheader,'hpmap_alp2.fits')
!      deallocate(hpmap)

      end program graytrix
