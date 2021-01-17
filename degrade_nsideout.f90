program degrade
use alm_tools
use healpix_types
use pix_tools
use udgrade_nr
use fitstools, only: getsize_fits, input_map, output_map
use head_fits
use omp_lib

IMPLICIT NONE

integer(i4b) :: nres,nside,lmax,hporder,hpspin,nmaps,narg
integer(i4b) :: nsidein,lmaxin,nsideout
integer(i8b) :: npix,ipring,ipnest,ipix,kpix,iray
integer(i8b) :: npixin,npixout
integer(4) :: j,i

real(dp), dimension(2) :: hpz
real(dp), allocatable, dimension(:,:) :: dw8

real(sp),allocatable :: hpmap(:,:) !hpmap(0:npix-1,3) output
real(sp),allocatable :: hpmapout(:,:) !hpmap(0:npix-1,3) output
complex(spc), allocatable, dimension(:,:,:) :: hpalm

character(len=128) :: fname,fnamei
character(len=128) :: outfile,outdir,fin,filein,outlog,fileouta
character(len=128) :: fileoutb
character(len=80), dimension(1:64) :: mapheader

character*128 opt,arg

narg=iargc()
do i=1,narg
    call getarg(i,opt)
    call getarg(i+1,arg)
    select case (opt)
    case ('-infile')
        fin=arg
    end select
enddo

filein=ADJUSTL(ADJUSTR(fin)//'.mag.dat')
fileouta=ADJUSTL(ADJUSTR(fin)//'.kg1g2.fits')
fileoutb=ADJUSTL(ADJUSTR(fin)//'.omega.fits')
write(*,*)'##### reading data: ',filein
write(*,*)'##### writing: ',fileouta
write(*,*)'##### writing: ',fileoutb

hporder=1 ! integer =1 for ring, =2 for nested
!hpz(1)=-1.0_dp
!hpz(2)=1.0_dp
hpz = (-1.d0,1.d0)

nmaps   = 4
lmax    = 24576
nsidein = 16384
nsideout= 8192
npixin  = nside2npix(nsidein)
npixout = nside2npix(nsideout)


open(12,file=filein,status="old",form='unformatted')

read(12) nside,npix

allocate(hpmap(0:npixin-1,1:4))

read(12) hpmap(0:npixin-1,1)  ! convergence
read(12) hpmap(0:npixin-1,2)  ! shear1
read(12) hpmap(0:npixin-1,3)  ! shear2
read(12) hpmap(0:npixin-1,4)  ! rotation
close(12)

!allocate(hpmap(0:npixin-1,1:nmaps))
!call input_map(filein, hpmap, npixin, nmaps)

!!!!!!! kappa/gamma1/gamma2 !!!!!!!!!!!!
allocate(hpmapout(0:npixout-1,1:3))
allocate(hpalm(1:3,0:lmax,0:lmax))
allocate(dw8(1:2*nsidein,1:3))
dw8=1.0_dp
call map2alm(nsidein, lmax,lmax,hpmap(:,1:3),hpalm,hpz,dw8)
call alm2map(nsideout,lmax,lmax,hpalm,hpmapout)

call write_minimal_header(mapheader, 'MAP', nside=nsideout,order=hporder)
call add_card(mapheader,'TTYPE1','kappa','label for field 1')
call add_card(mapheader,'TTYPE2','shear1','label for field 2')
call add_card(mapheader,'TTYPE3','shear2','label for field 3')
call output_map(hpmapout,mapheader,fileouta)
deallocate(hpmapout)
deallocate(hpalm)
deallocate(dw8)

!!!!!!!! omgega file !!!!!!!!!!
!allocate(hpmapout(0:npixout-1,1:1))
!allocate(hpalm(1:1,0:lmax,0:lmax))
!allocate(dw8(1:2*nsidein,1:1))
!dw8=1.0_dp
!call map2alm(nsidein, lmax,lmax,hpmap(:,4:4),hpalm,hpz,dw8)
!call alm2map(nsideout,lmax,lmax,hpalm, hpmapout)

!call write_minimal_header(mapheader, 'MAP', nside=nsideout,order=hporder)
!call add_card(mapheader,'TTYPE1','omega','label for field 1')
!call output_map(hpmapout,mapheader,fileoutb)
!deallocate(hpmap)
!deallocate(hpmapout)
!deallocate(hpalm)
!deallocate(dw8)

end program degrade
