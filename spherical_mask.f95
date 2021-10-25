module params
  ! Variables to be shared by all program units
  real :: pad_finest
  real, dimension(3) :: x0, lx, xc, shift
  
  integer :: levelmin, levelmax, nlevel, npad
  ! integer :: ix, iy, iz
  ! integer :: nx, ny, nz

  real :: omega_m, omega_v, boxlen, H0, astart, xi, yi, zi, r, r2
end module params

program ref_mask

  ! NOTES
  !
  ! off_abs is the absolute offset of the grids at ilevel from the
  ! origin, measured in ilevel grids
  !
  ! off_rel is the relative offset between the grids at ilevel and the
  ! grids at ilevel-1, measured in ilevel-1 grids

  ! TODO
  !
  ! perhaps bulding the refinement hierarchy is necessary all the
  ! way down past levelmin, in case the zoom region isn't consistent
  ! below there?

  use params

  implicit none
  
  type :: masks
     integer, pointer :: data(:, :, :)
  end type masks

  type(masks), allocatable :: rm(:)

  real, allocatable :: slab(:, :)
  real, allocatable :: data_w(:, :, :), data_coarse(:, :, :)

  character(len=800):: nml_file
  character(len=100) :: path
  character(len=4):: lnum
  character(len=6):: ldim
  character(len=3):: xyz = 'xyz'
  
  integer, allocatable :: n(:, :)
  integer, allocatable :: off_abs(:, :), off_rel(:, :)

  integer :: mv
  real :: dx
  real, dimension(3) :: x = 0.0
  integer, dimension(3) :: i_fine = 0

  logical :: ff
  logical :: check_point, check
  logical :: is_refined, is_in_mask
  
  integer :: ilevel
  integer :: ii, jj, kk

  ! Convenience and header variables for writing
  real :: dx_w
  real :: pvar_val
  real, dimension(3) :: off_w
  integer, dimension(3) :: n_w, n_coarse
  ! g9p header variables
  integer :: n1_g9p
  real :: dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p

  real :: ainv, adot, tmp1, growth, tmp, fupper, flower, growthVel, vFact
  integer :: ixc, iyc, izc, ixv, iyv, izv, k, kv, iax, nf
  integer :: ix, iy, iz, nx, ny, nz

  real, external :: func_dtdacube, rombint

  namelist /mask_params/ levelmin, levelmax, xi, yi, zi, r, npad, path

  ! Read the namelist
  if(command_argument_count() .lt. 1) then
     print*,'Usage: ./cubic_mask namelist'
     stop 1
  endif
  
  call get_command_argument(1, nml_file)

  open(11,file=trim(nml_file),status='old')
  read(11,nml=mask_params)
  close(11)
  print *,'Read the namelist:'
  print *,'---- levelmin', levelmin, 'levelmax', levelmax
  ! print *,'---- ix', ix, 'iy', iy, 'iz', iz
  ! print *,'---- nx', nx, 'ny', ny, 'nz', nz
  print *,'---- x', xi, 'y',  yi, 'z',  zi, 'r', r
  print *,'---- npad', npad
  print *,'---- path ', trim(path)
  print *,''

  nlevel = levelmax - levelmin
  pvar_val = 0.5
  
  ! Read in some parameters from the g9p file
  write(ldim, '(I5)') 2**levelmin
  open(13,file=trim(path)//'g9p'//trim(adjustl(ldim))//'_delta',form='unformatted',status='old')
  read(13) n1_g9p, n1_g9p, n1_g9p, dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, astart, omega_m, omega_v, H0
  print *, 'Read in parameters from ', 'g9p'//trim(adjustl(ldim))//'_delta'
  print *, '---- astart', astart, 'zstart', (1./astart) - 1.
  print *, '---- H0', H0, 'omega_m', omega_m, 'omega_v', omega_v
  print *, ''
  ! Calculate box length from the cell spacing on the coarsest level
  boxlen = (2**levelmin) * dx_g9p

  ! Calculate the region params
  r2 = r * r
  nf = 2 ** levelmax
  dx = 1.0 / real(nf)
  pad_finest = (real(npad) + 1.0) * dx

  ix = floor((xi - r) * nf)
  iy = floor((yi - r) * nf)
  iz = floor((zi - r) * nf)

  nx = ceiling((xi + r) * nf) - ix
  ny = ceiling((yi + r) * nf) - iy
  nz = ceiling((zi + r) * nf) - iz

  x0(1) = real(ix) * dx
  x0(2) = real(iy) * dx
  x0(3) = real(iz) * dx
  lx(1) = real(nx) * dx
  lx(2) = real(ny) * dx
  lx(3) = real(nz) * dx

  xc(1) = xi
  xc(2) = yi
  xc(3) = zi
  
  ! do ii = 1, 3
  !    xc(ii) = x0(ii) + 0.5*lx(ii)
  !    shift(ii) = xc(ii) - 0.5
  ! end do
  
  ! TESTING
  ! print *, 'x0', x0, 'lx', lx !, 'shift', shift
  
  allocate(n(nlevel+1, 3))
  allocate(off_abs(nlevel+1, 3))
  allocate(off_rel(nlevel+1, 3))

  ! Calculate extents and offsets
  call hierarchy(n, off_abs, off_rel)

  ! Build the refmask hierarchy
  allocate(rm(nlevel+1))
  do ilevel = 1, nlevel+1
     allocate(rm(ilevel)%data(0:n(ilevel, 1)-1, 0:n(ilevel, 2)-1, 0:n(ilevel, 3)-1))
     ! Initialise to zero
     rm(ilevel)%data = 0
  end do
  
  ! TESTING
  ! do ilevel = 1, nlevel+1
  !    print *, 'level ', levelmin + ilevel - 1
  !    print *, size(rm(ilevel)%data), shape(rm(ilevel)%data)
  ! end do

  ! Now we can start to build the refinement mask
  do ilevel = nlevel+1, 1, -1
     dx = 1.0 / (2**(ilevel+levelmin-1))

     ! TESTING
     ! print *, 'level ', ilevel + levelmin -1, 'dx ', dx
     
     ! Check whether each index is within the refinement
     ! region. Stepping in twos because if icell is in the mask then
     ! its immediate neighbour must be

     ! Calculate position of cell centre
     do ii = 0, n(ilevel, 1)-1, 2
        x(1) = (real(off_abs(ilevel, 1)) + real(ii) + 0.5)*dx ! + shift(1)
        do jj = 0, n(ilevel, 2)-1, 2
           x(2) = (real(off_abs(ilevel, 2)) + real(jj) + 0.5)*dx! + shift(2)
           do kk = 0, n(ilevel, 3)-1, 2
              x(3) = (real(off_abs(ilevel, 3)) + real(kk) + 0.5)*dx! + shift(3)

              ! Initialise the mask value to be outside of the mask
              mv = -1

              ! Check whether the point is within the mask
              if (check_point(x, ilevel + levelmin -1) .or. (ilevel .eq. 1)) then
                 mv = 1
              end if

              ! Update the value of that point and its neighbours
              rm(ilevel)%data(ii, jj, kk) = mv
              rm(ilevel)%data(ii+1, jj, kk) = mv
              rm(ilevel)%data(ii, jj+1, kk) = mv
              rm(ilevel)%data(ii, jj, kk+1) = mv
              rm(ilevel)%data(ii+1, jj+1, kk) = mv
              rm(ilevel)%data(ii+1, jj, kk+1) = mv
              rm(ilevel)%data(ii, jj+1, kk+1) = mv
              rm(ilevel)%data(ii+1, jj+1, kk+1) = mv           
           end do
        end do
     end do
  end do

  ! Now check whether the fine cells are flagged
  do ilevel = 1, nlevel, 1
     do ii = 1, n(ilevel, 1)
        do jj = 1, n(ilevel, 2)
           do kk = 1, n(ilevel, 3)
              ! Initialise fine flagged check
              ff = .false.

              ! Calculate corresponding fine coordinates for current
              ! grid cell
              
              i_fine(1) = 2*ii - 2*off_rel(ilevel+1, 1)
              i_fine(2) = 2*jj - 2*off_rel(ilevel+1, 2)
              i_fine(3) = 2*kk - 2*off_rel(ilevel+1, 3)

              ! Determine whether the fine descendants of this coarse
              ! cell are flagged

              ! if (((i_fine(1) .ge. 1) .and. (i_fine(1) .le. n(ilevel+1, 1))) &
              !      .and. ((i_fine(2) .ge. 1) .and. (i_fine(2) .le. n(ilevel+1, 2))) &
              !      .and. ((i_fine(3) .ge. 1) .and. (i_fine(3) .le. n(ilevel+1, 3)))) then

              if (((i_fine(1) .ge. 0) .and. (i_fine(1) .lt. n(ilevel+1, 1))) &
                   .and. ((i_fine(2) .ge. 0) .and. (i_fine(2) .lt. n(ilevel+1, 2))) &
                   .and. ((i_fine(3) .ge. 0) .and. (i_fine(3) .lt. n(ilevel+1, 3)))) then

                 ff = ff .or. rm(ilevel+1)%data(i_fine(1), i_fine(2), i_fine(3)) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1)+1, i_fine(2), i_fine(3)) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1), i_fine(2)+1, i_fine(3)) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1), i_fine(2), i_fine(3)+1) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1)+1, i_fine(2)+1, i_fine(3)) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1), i_fine(2)+1, i_fine(3)+1) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1)+1, i_fine(2), i_fine(3)+1) > 0
                 ff = ff .or. rm(ilevel+1)%data(i_fine(1)+1, i_fine(2)+1, i_fine(3)+1) > 0
                 
                 ! Now update refmask value
                 if (ff) then
                    rm(ilevel)%data(ii, jj, kk) = 2  ! Then this coarse cell is refined

                    ! Then these fine cells are in the mask
                    rm(ilevel+1)%data(i_fine(1), i_fine(2), i_fine(3)) = 1
                    rm(ilevel+1)%data(i_fine(1)+1, i_fine(2), i_fine(3)) = 1
                    rm(ilevel+1)%data(i_fine(1), i_fine(2)+1, i_fine(3)) = 1
                    rm(ilevel+1)%data(i_fine(1), i_fine(2), i_fine(3)+1) = 1
                    rm(ilevel+1)%data(i_fine(1)+1, i_fine(2)+1, i_fine(3)) = 1
                    rm(ilevel+1)%data(i_fine(1), i_fine(2)+1, i_fine(3)+1) = 1
                    rm(ilevel+1)%data(i_fine(1)+1, i_fine(2), i_fine(3)+1) = 1
                    rm(ilevel+1)%data(i_fine(1)+1, i_fine(2)+1, i_fine(3)+1) = 1
                 end if
              end if
           end do
        end do
     end do
  end do

  
  ! TESTING
  print *, 'For level', ilevel + levelmin - 1
  print *, '----', real(sum(rm(ilevel)%data, mask=rm(ilevel)%data .gt. 0))/real(size(rm(ilevel)%data)), 'of cells in mask'

 
  ! With the refinement mask in hand, we can write out the IC files,
  ! starting with the maximum level
  write(lnum, '(I4)') nlevel + levelmin + 1000
  open(12,file='level_'//lnum(2:4)//'/ic_refmap',form='unformatted')
  open(15,file='level_'//lnum(2:4)//'/ic_pvar_00001',form='unformatted')
  
  ! Write header
  dx_w = boxlen / (2**levelmin + nlevel -1)
  do ii = 1, 3
     n_w(ii) = n(nlevel+1, ii)
     off_w(ii) = off_abs(nlevel+1, ii) * dx_w
  end do
  write(12) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0
  write(15) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0

  ! Convert the finest level to a grafic refmap
  allocate(data_w(0:n_w(1)-1, 0:n_w(2)-1, 0:n_w(3)-1))
  do ii = 0, n_w(1)-1
     do jj = 0, n_w(2)-1
        do kk = 0, n_w(3)-1
           if ((is_in_mask(rm(nlevel+1)%data(ii, jj, kk))) .and. &
                (.not. is_refined(rm(nlevel+1)%data(ii, jj, kk)))) then
              data_w(ii, jj, kk) = 1.0
           else
              data_w(ii, jj, kk) = 0.0
           end if
        end do
     end do
  end do

  ! Write the data
  do kk = 0, n_w(3)-1
     write(12) (((data_w(ii, jj, kk)), ii=0, n_w(1)-1), jj=0,n_w(2)-1)
     write(15) (((data_w(ii, jj, kk) * pvar_val), ii=0, n_w(1)-1), jj=0,n_w(2)-1)
  end do

  close(12)
  close(15)
  
  ! Now need to do the refmap for coarser levels
  do ilevel = nlevel, 1, -1
     ! Start by calculating the header variables
     dx_w = boxlen / (2**(ilevel + levelmin -1))
     do ii = 1, 3
        n_coarse(ii) = n(ilevel, ii)
        ! Technically off_w should be off_coarse, but we don't use
        ! off_coarse in restrict_mask so we can skip defining it as
        ! coarse first
        off_w(ii) = off_abs(ilevel, ii) * dx_w
     end do

     ! TESTING
     print *, 'Restricting mask for level ', ilevel+levelmin-1
     
     allocate(data_coarse(0:n_coarse(1)-1, 0:n_coarse(2)-1, 0:n_coarse(3)-1))
     call restrict_mask(n_coarse, n_w, off_rel(ilevel+1, :), data_coarse, data_w)

     ! Now reassign variables
     n_w = n_coarse
     deallocate(data_w)
     allocate(data_w(0:n_w(1)-1, 0:n_w(2)-1, 0:n_w(3)-1))
     data_w = data_coarse

     ! Since restrict_mask modifies n_w to n_coarse and data_w to
     ! data_coarse, we can go ahead and write them out
     write(lnum, '(I4)') ilevel + levelmin + 999
     open(12,file='level_'//lnum(2:4)//'/ic_refmap',form='unformatted')
     open(15,file='level_'//lnum(2:4)//'/ic_pvar_00001',form='unformatted')

     ! Write header
     write(12) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0
     write(15) n_w(1), n_w(2), n_w(3), dx_w, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0

     ! Write data
     do kk = 0, n_w(3)-1
        write(12) (((data_w(ii, jj, kk)), ii=0, n_w(1)-1), jj=0,n_w(2)-1)
        write(15) (((data_w(ii, jj, kk) * pvar_val), ii=0, n_w(1)-1), jj=0,n_w(2)-1)
     end do

     ! Now we can deallocate data_coarse ready for the next loop
     deallocate(data_coarse)
     close(12)
     close(15)
  end do

  deallocate(data_w)

  ! Now we are done with the refmap we can deallocate to free up memory
  deallocate(rm)

  ! Now write the rest of the IC files
  do ilevel = nlevel+1, 1, -1
     print *, 'Writing IC files for level ', ilevel + levelmin -1
     ! Utility variables for loading files
     ! Maximum possible number of cells in this level along each dimension
     write(lnum, '(I4)') ilevel + levelmin - 1 + 1000
     write(ldim, '(I5)') 2**(ilevel + levelmin - 1)
     
     ! Bits of the following do loop is (almost) directly pinched from
     ! Sergey Pilipenko's cubic_mask3.f90
     do iax = 0,3
        if(iax.eq.0) then
           open(13,file=trim(path)//'g9p'//trim(adjustl(ldim))//'_delta',form='unformatted',status='old')
           open(12,file='level_'//lnum(2:4)//'/ic_deltab',form='unformatted')
        else 
           open(13,file=trim(path)//'g9p'//trim(adjustl(ldim))//'_vel'//xyz(iax:iax),form='unformatted',status='old')
           open(12,file='level_'//lnum(2:4)//'/ic_velc'//xyz(iax:iax),form='unformatted')
           open(14,file='level_'//lnum(2:4)//'/ic_posc'//xyz(iax:iax),form='unformatted')
        endif
        read(13) n1_g9p, n1_g9p, n1_g9p, dx_g9p, x1o_g9p, x2o_g9p, x3o_g9p, astart, omega_m, omega_v, H0

        ! Calculate grid offsets in comoving Mpc and extent in grid cells
        do ii = 1, 3
           n_w(ii) = n(ilevel, ii)
           off_w(ii) = off_abs(ilevel, ii) * dx_g9p
        end do

        write(12) n_w(1), n_w(2), n_w(3), dx_g9p, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0
        if(iax.gt.0) write(14) n_w(1), n_w(2), n_w(3), dx_g9p, off_w(1), off_w(2), off_w(3), astart, omega_m, omega_v, H0

        allocate(slab(n1_g9p, n1_g9p))

        ainv = 1. / astart
        adot = sqrt(omega_m*(ainv-1) + omega_v*(astart*astart-1) + 1)
        tmp1 = rombint(func_dtdacube, 1e-6, astart, 1e-5)
        growth = 2.5 * omega_m * adot * tmp1 / astart
        tmp  = 1. - omega_m - omega_v
        fupper = 2.5 * omega_m / growth - 1.5 * omega_m * ainv - tmp
        flower = omega_m * ainv + omega_v * astart * astart + tmp
        growthVel = fupper/flower

        vFact = 1.  / ( (adot) * 100. * (growthVel) )

        ixv = nint(x1o_g9p/dx_g9p)
        iyv = nint(x2o_g9p/dx_g9p)
        izv = nint(x3o_g9p/dx_g9p)

        ! Calculate i*c, since on levelmin we start from 1
        if (ilevel .ne. 1) then
           ixc = off_abs(ilevel, 1)
           iyc = off_abs(ilevel, 2)
           izc = off_abs(ilevel, 3)
        else
           ixc = 1
           iyc = 1
           izc = 1
        end if

        do k = 1, n1_g9p
           read(13) slab
           kv = k + izv
           ! if (ilevel .eq. 1) kv
           ! LC - swapped .ge. for .gt. and .lt. for .le. to account for 1-indexing and added +1
           ! LC 7/5 - swapped back
           if((kv .ge. izc) .and. (kv .lt. izc + n_w(3))) then
              write(12) slab(ixc-ixv:ixc+n_w(1)-1-ixv, &
                   iyc-iyv:iyc+n_w(2)-1-iyv)
              if(iax.gt.0) then
                 write(14) slab(ixc-ixv:ixc+n_w(1)-1-ixv, &
                      iyc-iyv:iyc+n_w(2)-1-iyv)  * vFact
              endif
           endif
        enddo
        deallocate(slab)

        close(13)
        close(12)
        if(iax.gt.0) close(14)
     enddo
  end do

  ! Finally, write the namelist
  call write_namelist()
  
! Using internal subroutines since the size of n, off_abs, and off_rel
! will be dynamically allocated and hence we need expicit interfaces
contains
  
  subroutine hierarchy(n, off_abs, off_rel)
    ! Determines the refinement hierarchy by calculating the extent
    ! and offsets of the grids in each level.
    !
    ! Notes
    !
    ! Padding is taken care of in here
    !
    ! Arguments
    !
    !    n (int) [nlevel, 3] - extent of each level in units of that
    !    level's grid cells
    !
    !    off_abs (int) [nlevel, 3] - offset of each level from the
    !    origin, in units of that level's grid cells
    !
    !    off_rel (int) [nlevel, 3] - offset of each level from the
    !    next coarser level, in units of the next coarser level's grid
    !    cells

    use params
    
    integer, dimension(:, :), intent(out) :: n
    integer, dimension(:, :), intent(out) :: off_abs
    integer, dimension(:, :), intent(out) :: off_rel
    integer :: ilevel
    integer, dimension(3) :: il = 0 ! LHS of grid region
    integer, dimension(3) :: ir = 0 ! RHS of grid region
  
    ! Initialise the relevant variables to 0
    n = 0
    off_abs = 0
    off_rel = 0

    ! We can immediately set the extent of the coarsest level
    n(1, :) = 2**levelmin

    ! Start the calculations with the finest level
    print *, 'Calculating extent of level ', levelmax
    il(1) = ix
    il(2) = iy
    il(3) = iz
    ir(1) = il(1) + nx
    ir(2) = il(2) + ny
    ir(3) = il(3) + nz
  
    ! Align with coarser grid, +1 in the mod is because Fortran is
    ! one-indexed
    do ii = 1, 3
       il(ii) = il(ii) - mod(il(ii), 2)
       ir(ii) = ir(ii) + mod(ir(ii), 2)
       n(nlevel+1, ii) = ir(ii) - il(ii)
       off_abs(nlevel+1, ii) = il(ii)
    end do

    ! TESTING
    ! print *, 'FINEST LEVEL'
    ! print *, 'ilx ', il(1), 'ixmax ', ix, 'nx ', n(nlevel+1, 1), 'nx ', nx
    ! print *, 'ily ', il(2), 'iymax ', iy, 'ny ', n(nlevel+1, 2), 'ny ', ny
    ! print *, 'ilz ', il(3), 'izmax ', iz, 'nz ', n(nlevel+1, 3), 'nz ', nz

    ! Find the region and extent of coarser levels
    do ilevel = nlevel, 2, -1
       print *, 'Calculating extent of level ', ilevel + levelmin - 1

       ! Calculate the bounding region of the current level
       do ii = 1, 3
          il(ii) = int(real(il(ii))/2.0 - real(npad))
          ir(ii) = int(real(ir(ii))/2.0 + real(npad))
        
          ! Check the alignment with coarser grids
          il(ii) = il(ii) - mod(il(ii), 2)
          ! ir(ii) = ir(ii) - mod(ir(ii) + 1, 2) 
          ir(ii) = ir(ii) + mod(ir(ii), 2) ! LC 2/8 -- swapped to + mod

          ! Store
          n(ilevel, ii) = ir(ii) - il(ii)
          off_abs(ilevel, ii) = il(ii)
       end do
    end do
     
    ! Find the relative offsets between levels
    do ilevel = nlevel+1,  2, -1
       print *, 'Calculating relative offset for level ', ilevel + levelmin - 1

       do ii = 1, 3
          off_rel(ilevel, ii) = (off_abs(ilevel, ii) / 2) - off_abs(ilevel-1, ii)
       end do
    end do



    ! Now check the absolute offsets are correct (doing the integer
    ! division before might have messed some things up)

    ! LC - testing
    
    do ilevel = 2, nlevel+1, 1
       print *, 'Checking absolute offset for level ', ilevel + levelmin - 1

       do ii = 1, 3
          print*, 'before'
          print*, 'ilevel', ilevel, 'ii', ii, 'off_abs', off_abs(ilevel, ii)
          print*, 'off_rel', off_rel(ilevel, ii)
          off_abs(ilevel, ii) = 2*off_rel(ilevel, ii) + 2*off_abs(ilevel-1, ii)

          ! Re-check alignment with coarse grids
          ! off_abs(ilevel, ii) = off_abs(ilevel, ii) - mod(off_abs(ilevel, ii) + 1, 2)

          print*, 'after'
          print*, 'ilevel', ilevel, 'ii', ii, 'off_abs', off_abs(ilevel, ii)

       end do
    end do
    
  end subroutine hierarchy

  subroutine restrict_mask(n_coarse, n_w, off_rel, data_coarse, data_w)
    ! Restricts the coarse mask to be consistent with the fine
    ! mask. Also produces a grafic refinment map that is ready to be
    ! written. As a final step, n_w is set equal to n_coarse and
    ! data_w is set equal to data_coarse.
    !
    ! Arguments
    !
    !     n_coarse (int) [3] - the extent of the ``coarser'' level,
    !     which is the level that is about to be written
    !
    !     n_w (int) [3] - the extent of the ``finer'' level, which is
    !     the level that has just been written, and at the end is set
    !     equal to n_coarse
    !
    !     off_rel (int) [3] - the relative offset between the
    !     ``coarser'' and ``finer'' level, which is accessed by
    !     off_rel(i_level_coarser, :)
    !
    !     data_coarse (real) [n_coarse(1), n_coarse(2), n_coarse(3)] -
    !     empty allocated array into which the restricted ``coarser''
    !     mask will be written, this is only a temporary variable
    !     since, in the end, data_w is set equal to it
    !
    !     data_w (real) [n_w(1), n_w(2), n_w(3)] - the mask from the
    !     ``finer'' level, which should have already been written, and
    !     at the end is set equal to data_coarse
    
    integer, dimension(3), intent(inout) :: n_coarse, n_w
    integer, dimension(3), intent(in) :: off_rel
    ! real, dimension(:, :, :), intent(inout) :: data_coarse
    real, allocatable, dimension(:, :, :), intent(inout) :: data_w, data_coarse

    integer :: ii, jj, kk, ii_c, jj_c, kk_c
    ! Initialise coarse mask
    data_coarse = 0.0

    ! TESTING
    ! print*, 'off_rel', off_rel
    ! print*, 'size(data_coarse)', size(data_coarse), 'shape(data_coarse)', shape(data_coarse)
    
    ! Restrict, loop is over the fine (write) extent
    do ii = 0, n_w(1) - 1
       ii_c = ii/2 + off_rel(1)

       ! TESTING
       ! print *, 'ii_c', ii_c
       
       do jj = 0, n_w(2) - 1
          jj_c = jj/2 + off_rel(2)

          ! TESTING
          ! print *, 'jj_c', jj_c
          
          do kk = 0, n_w(3) - 1
             kk_c = kk/2 + off_rel(3)

             ! TESTING
             ! print *, 'kk_c', kk_c


             if (data_w(ii, jj, kk) .gt. 0.0) then
                data_coarse(ii_c, jj_c, kk_c) = data_coarse(ii_c, jj_c, kk_c) + 1.0
             end if
          end do
       end do
    end do

    ! Now build the grafic map for the coarse level
    do ii = 0, n_coarse(1)-1
       do jj = 0, n_coarse(2)-1
          do kk = 0, n_coarse(3)-1
             if (data_coarse(ii, jj, kk) .gt. 0.0) then
                data_coarse(ii, jj, kk) = 1.0
             end if
          end do
       end do
    end do
    
    ! ! Now reassign the coarse variables to the write variables
    ! n_w = n_coarse

    ! deallocate(data_w)
    ! allocate(data_w(n_w(1), n_w(2), n_w(3)))
    ! data_w = data_coarse
             
  end subroutine restrict_mask

  
  subroutine write_namelist
    use params

    integer :: ilev
    character(len=64) :: nexpand, str_levelmin
    character(len=3), dimension(10) :: i_nexpand

    open(20, file='nml.nml')!, status='new')

    ! RUN_PARAMS
    write(20, '(a)') '&RUN_PARAMS'
    write(20, '(a)') 'cosmo=.true.'
    write(20, '(a)') 'pic=.true.'
    write(20, '(a)') 'poisson=.true.'
    write(20, '(a)') 'hydro=.true.'
    write(20, '(a)') 'nrestart=0'
    write(20, '(a)') 'nremap=10'
    write(20, '(a)') 'ncontrol=1'
    write(20, '(a)') 'verbose=.false.'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! INIT_PARAMS
    write(20, '(a)') '&INIT_PARAMS'
    write(20, '(a)') "filetype='grafic'"
    do ilev = levelmin, levelmax
       write(20, fmt='(a, i1, a, i3.3, a)') "initfile(",ilev-levelmin+1,")='level_",ilev,"'"
    end do
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! AMR_PARAMS
    write(20, '(a)') '&AMR_PARAMS'
    write(str_levelmin, '(i2)') levelmin
    str_levelmin = trim(adjustl(str_levelmin))
    write(20, fmt='(a, a)') 'levelmin=', str_levelmin
    write(20, fmt='(a, i2)') 'levelmax=', levelmax + 9
    write(20, '(a)') 'ngridmax=100000'
    write(20, '(a)') 'npartmax=200000'
    if (levelmin .eq. levelmax) then
       write(20, '(a)') 'nexpand=1'
    else
       nexpand='nexpand='
       do ilev=levelmin, levelmax-2
          write(i_nexpand(ilev-levelmin+1), '(i2, a)') (npad-1)/2, ','
          nexpand = trim(adjustl(nexpand))//trim(adjustl(i_nexpand(ilev-levelmin+1)))
       end do
       nexpand = trim(nexpand)//'1,1'
       write(20, '(a)') trim(nexpand)
    end if

    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! REFINE_PARAMS
    write(20, '(a)') '&REFINE_PARAMS'
    write(20, fmt='(a, i2, a)') 'm_refine=', levelmax - levelmin + 10,'*8.,'
    write(20, fmt='(a, i1)') 'ivar_refine=',6  ! ic_pvar_00001
    write(20, fmt='(a, f5.3)') 'var_cut_refine=', 0.01
    write(20, fmt='(a, e11.5)') 'mass_cut_refine=', 2.0/(2.0 ** (levelmax * 3.0))
    write(20, '(a)') 'interpol_var=1'
    write(20, '(a)') 'interpol_type=0'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! OUTPUT_PARAMS
    write(20, '(a)') '&OUTPUT_PARAMS'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    ! HYDRO_PARAMS
    write(20, '(a)') '&HYDRO_PARAMS'
    write(20, '(a)') '/'
    write(20, '(a)') ''

    close(20)

  end subroutine write_namelist
 
end program ref_mask


logical function check_point(xx, ilevel)
  use params
  integer, intent(in) :: ilevel
  real, dimension(3), intent(in) :: xx
  ! real, intent(in) :: r2  ! r^2
  
  real :: rr = 0.
  integer :: ii, jj, kk
  logical :: check

  ! Initialise
  ! check = .true.
  
  do ii = 1, 3
     rr = rr + (xx(ii) - xc(ii)) ** 2.
  end do
  
     ! Do extra padding check
     ! if (xx .lt. -0.5) thena
     !    xx = xx + 1.0
     ! else if (xx .gt. 0.5) then
     !    xx = xx - 1.0
     ! end if

     ! TESTING
     ! print *, 'dx after', xx
     ! print *, 'lx_ii', lx_ii
     ! print *, 'x_ii', x(ii)
     
  if (rr .lt. r2) then
     check_point = .true.
  else
     check_point = .false.
  end if
     !check = check .and. ((xx .ge. pad_finest) .and. (xx .le. lx(ii) - pad_finest))
     

     ! TESTING
     ! if ((check .eqv. .true.) .and. (ii .eq. 3)) then
     !    print *, 'true', count
     ! end if

  ! TESTING
  ! if (check) then
  !    print *, 'check', check
  ! end if 

  ! check_point = check
end function check_point


logical function is_refined(x)
  ! Function to check whether a point in the mask is in the mask and
  ! refined (corresponds to a mask value equal to 2)
  integer, intent(in) :: x

  is_refined = (x .eq. 2)
end function is_refined


logical function is_in_mask(x)
  ! Function to check whether a point in the mask is in the mask and
  ! refined (corresponds to a mask value equal to 2)
  integer, intent(in) :: x

  is_in_mask = (x .ge. 0)
end function is_in_mask

! The following two functions are also directly pinched from Sergey
! Pilipenko's cubic_mask3.f90
real function func_dtdacube(a)
  !  real,intent(in)::a
  use params
  
  func_dtdacube = 1. / (omega_m * (1./a - 1.0) + omega_v * (a * a - 1.0) + 1.0)**1.5
end function func_dtdacube

real function rombint(f,a,b,tol)
        parameter (MAXITER=30,MAXJ=5)
        implicit real (a-h,o-z)
        dimension g(MAXJ+1)
        real a, b, tol
        external f
!        interface
!          real function f(x)
!            real,intent(in)::x
!          end function f
!        end interface
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol) write(*,*) 'Rombint failed to converge; integral, error=',rombint,error
        return
end function rombint
