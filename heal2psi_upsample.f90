
program heal2psi
use alm_tools
USE healpix_types
use pix_tools
use fitstools
use omp_lib

IMPLICIT NONE

integer(i4b) :: nres,nside,lmax,hporder,hpspin,nmaps,iter,narg,nsideout
integer(i8b) :: npix,ipring,ipnest,ipix,np,npixtot,l,m,npixout
character(len=80) :: header
complex(spc), allocatable, dimension(:,:,:) :: hpalm
real(sp), allocatable :: hpmap(:,:)
real(dp), allocatable, dimension(:,:) :: dw8 
real(dp), dimension(2) :: hpz 
real(sp),allocatable :: hpmap_phi(:),hpmap_alp(:,:)
real(sp),allocatable :: hpmap_mag(:,:)
integer(4) :: j,i
real(dp), dimension(2) :: z
character(len=1024) :: fname,indir,outdir,fin,filein,outlog,inprefix,outprefix
character(len=1024) :: path,dir_out,infile,outfile
character*256 opt,arg
!implicit real (a-h,o-z)
!implicit integer (i-n)  

narg=iargc()
write(*,*)narg

do i=1,narg
    call getarg(i,opt)
    call getarg(i+1,arg)
    select case (opt)
    case ('-infile')    ! input prefix
!    WRITE (*,*) arg
        read (arg,'(A)') infile 
    case ('-nside')
        read (arg,'(I8)') nside
!        read(arg,*)nside  ! HEALPix nside input maps 
    case ('-lmax')
        read (arg,'(I8)') lmax   ! HEALPix lmax to do SHTs
    case ('-outfile')
        read (arg,'(A)') outfile      ! output prefix
    end select
enddo
!write(*,*)'*** main configuration' 
write(*,"(A)",advance="no") "-nside of input maps:"
write(*,'(I8)') nside
write(*,"(A)",advance="no") "-lmax of SHTS:"
write(*,'(I8)') lmax
write(*,"(A)",advance="no") "-input prefix:"
write(*,'(A)') infile
write(*,"(A)",advance="no") "-output prefix:"
write(*,'(A)') outfile

!inprefix='/home/users/yomori/oak/mdpl2/dens/cmbkappa/mdpl2_cmbkappa_11_z0.098.fits'
!outprefix='aaa.psis'

!nside=8192
!npix=805306368
!lmax=1024
nmaps=1
!iter=3
!allocate(hpalm(1:1,0:lmax,0:lmax))
!allocate(hpmap_phi(0:npix-1))
!allocate(hpmap_alp(0:npix-1,1:2))
!allocate(hpmap_mag(0:npix-1,1:3))

!!! Healpix related
hporder=1 ! integer =1 for ring, =2 for nested
hpz(1)=-1.0_dp
hpz(2)=1.0_dp


!!! read phi alm !!!
!call fits2alms('alms.fits', 65*66/2, hpalm, 1, header, 80, 3)
!use pix_tools, only: nside2npix
!use fitstools, only: getsize_fits, input_map

!npixtot = getsize_fits('//project2/chihway/sims/MDPL2/cmbkappa/mdpl2_cmbkappa_4_229.fits',nmaps=nmaps, nside=nside)
nsideout=16384
npix  = nside2npix(nside)
npixout = nside2npix(nsideout)

allocate(hpmap(0:npix-1,1:1))
allocate(hpalm(1:1,0:lmax,0:lmax))
allocate(dw8(1:2*nside,1:1))
dw8=1.0_dp
write(*,*)'*** reading fits file'
call input_map(infile, hpmap, npix, nmaps)
!call input_map('//project2/chihway/sims/MDPL2/cmbkappa/mdpl2_cmbkappa_4_229.fits', hpmap, npix, nmaps)
!print *, hpmap(npix,1)
!do, i=1:npix
!        write *, hpmap(i)
!enddo
write(*,*)'*** doing a map2alm transform'
!z = sin(10.0_dp * DEG2RAD)
hpz = (-1.d0,1.d0)
call map2alm(nside,lmax,lmax,hpmap(:,1),hpalm,hpz,dw8)
deallocate(dw8)
deallocate(hpmap)

write(*,*)'*** converting klm -> plm'
!$omp parallel private(l)
!$omp do
      do m=0,lmax
         hpalm(1,0,m)=0.0 !psi_lm=0 for l=0
         do l=1,lmax
            hpalm(1,l,m)=hpalm(1,l,m)/(0.5*l*(l+1))
         enddo
      enddo
!$omp end do
!$omp end parallel


allocate(hpmap_phi(0:npixout-1))
allocate(hpmap_alp(0:npixout-1,1:2))
allocate(hpmap_mag(0:npixout-1,1:3))

write(*,*)'*** computing first and second derivatives of pot'
call alm2map_der(nsideout,lmax,lmax,hpalm,hpmap_phi,hpmap_alp,hpmap_mag)
deallocate(hpmap_phi)
deallocate(hpalm)

write(*,*)'*** hpmap_alp(0:npix-1,1:2) => hpmap(1:2,0:npix-1)'
allocate(hpmap(1:2,0:npixout-1))
hpmap=transpose(hpmap_alp)
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,2
!            hpmap(j,ipix)=hpmap_alp(ipix,j)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
deallocate(hpmap_alp)

!outfile=ADJUSTL(ADJUSTR(pfname)//'.psis')
!outfile=ADJUSTL(ADJUSTR(outdir)//'/'//outfile)

write(*,"(A,A)",advance="no") "*** writing to psis file:",outfile
open(10,file=outfile,status='unknown',form='unformatted')
write(10)nsideout,npixout
write(10)hpmap
deallocate(hpmap)


write(*,*)'*** changing hpmap_mag(0:npix-1,1:3) => hpmap(1:3,0:npix-1)'
allocate(hpmap(1:3,0:npixout-1))
hpmap=transpose(hpmap_mag)
!!$omp parallel private(j)
!!$omp do
!      do ipix=0,npix-1
!         do j=1,3
!            hpmap(j,ipix)=hpmap_mag(ipix,j)
!         enddo
!      enddo
!!$omp end do
!!$omp end parallel
deallocate(hpmap_mag)
write(10)hpmap
close(10)
write(*,*)'*** finished everything'
end program heal2psi
