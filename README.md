# graytrix
(I am not the original code developer, I just wrote the helper functions)

Graytrix is written by T. Hamana and the original page can be found here http://th.nao.ac.jp/MEMBER/hamanatk/GRayTrix/

If you use Graytrix for scientific work, please reference the following papers in which Graytrix algorithm is described:

1. Cosmological constraints from Subaru weak lensing cluster counts
Takashi Hamana, Junya Sakurai, Michitaro Koike, Lance Miller
PASJ, 67, 34 (2015)

2. Probing cosmology with weak lensing selected clusters I: Halo approach and all-sky simulations
Masato Shirasaki, Takashi Hamana, Naoki Yoshida
MNRAS, 453, 3043


Steps
--------
1. Project all the particles (x/y/z -> vec2pix -> HEALPix maps)
2. Run make_kappaslices.py to convert shells into local convergence<br>
   ```python
   python make_kappaslices.py ${dir_in} ${dir_out} ${cambpar}
   ```
3. Run heal2psis.x to compute the first and second derivatives<br>
   ```python
   ./heal2psi.x -nside ${nside} -lmax ${lmax} -infile ${inputfilename} -outfile ${outputfile}```
4. Run graytrix.x
   ```python
   ./graytrix.x -nres ${nres} -d ${outdir} -om ${omegamatter} -th ${shellwidth} -nt ${OMP_NUM_THREADS} -f ${prefix} < infile
   ```
   where<BR>
   ```nres```: HEALPix resolution 4096->12, 8192->13<BR>
   ```outdir```: output directory<BR>
   ```omegamatter```: \Omega_{\rm m}<BR>
   ```shellwidth```: width of the shells in Mpc/h<BR>
   ```OMP_NUM_THREAD```: number of threads to use<BR>
   ```prefix```: prefix of the output files<BR>

WARNING: I have been finding that ifort gives unknown spurious results resulting in random corrupted healpix files. This can be avoided by compiling gfortran.
