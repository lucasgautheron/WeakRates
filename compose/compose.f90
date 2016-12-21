!***********************************************************************
!***********************************************************************
! Fortran90 program from the CompOSE website http://compose.obspm.fr
!    for reading, checking, interpolating and writing EoS tables.
! Required standard input files
!    from website: eos.t, eos.nb, eos.yq, eos.thermo;
!    supplied by user: eos.parameters, eos_quantities
! Optional additional input files
!    from website: eos.compo, eos.micro.
! See manual for details and alternative options.
! CompOSE core team, version 0.25, 2016/06/15
!***********************************************************************
!***********************************************************************
PROGRAM compose
!***********************************************************************
!***********************************************************************
! Stefan Typel for the CompOSE core team, version 0.14, 2012/12/20
! working with CompOSE EoS tables
implicit none
integer :: iwr

! 1: write information and report to terminal, else not
iwr = 1

! read and analyze eos data tables
! to be called only once
!+++++++++++++++++++++++
call init_eos_table(iwr)
!+++++++++++++++++++++++

! define quantities in eos to be interpolated
!+++++++++++++++++++++++++
call define_eos_table(iwr)
!+++++++++++++++++++++++++

! generate eos table
!++++++++++++++++++++++
call get_eos_table(iwr)
!++++++++++++++++++++++

! output quantities
! eos_thermo(1): p [MeV fm^-3]
! eos_thermo(2): S [dimensionless]
! eos_thermo(3): mu_b - mu_n [MeV]
! eos_thermo(4): mu_q [MeV]
! eos_thermo(5): mu_l [MeV]
! eos_thermo(6): F/m_n-1 [dimensionless]
! eos_thermo(7): E/m_n-1 [dimensionless]
! eos_thermo(8): H/m_n-1 [dimensionless]
! eos_thermo(9): G/m_n-1 [dimensionless]
! eos_thermo(10): dp/dnb|E [MeV]
! eos_thermo(11): dp/dE|nb [fm^-3]
! eos_thermo(12): c_s^2 [dimensionless]
! eos_thermo(13): c_V [dimensionless]
! eos_thermo(14): c_p [dimensionless]
! eos_thermo(15): Gamma [dimensionless]
! eos_thermo(16): alpha_p [MeV^-1]
! eos_thermo(17): beta_V [fm^-3]
! eos_thermo(18): kappa_T [MeV^-1 fm^3]
! eos_thermo(19): kappa_S [MeV^-1 fm^3]
! eos_thermo_add(1), ..., eos_thermo(n_add)
! eos_compo_p(1), ..., eos_compo_p(n_p)
! eos_compo_q(1,1), ..., eos_compo_q(n_q,1): A^av
! eos_compo_q(1,2), ..., eos_compo_q(n_q,2): Z^av
! eos_compo_q(1,3), ..., eos_compo_q(n_q,3): Y^av
! eos_micro(1), ..., eos_micro(n_m)
! eos_thermo_err(1): absolute error in f/n [MeV]
! eos_thermo_err(2): relative error in f/n [dimensionless]
! eos_thermo_err(3): absolute error in e/n [MeV]
! eos_thermo_err(4): relative error in e/n [dimensionless]
! eos_thermo_err(5): absolute error in p/n [MeV]
! eos_thermo_err(6): relative error in p/n [dimensionless]
! eos_thermo_err(7): absolute error in s/n [dimensionless]
! eos_thermo_err(8): relative error in s/n [dimensionless]

stop
!***********************************************************************
end PROGRAM compose
!***********************************************************************
!***********************************************************************
! subroutines
!***********************************************************************
!***********************************************************************
SUBROUTINE init_eos_table(iwr)
! Stefan Typel for the CompOSE core team, version 0.09, 2012/11/02
implicit none
integer :: nbl,iwr

! maximum number of subtables
nbl = 10

call read_eos_tables_tny(iwr)

call read_eos_table_thermo(iwr,nbl)

call read_eos_table_compo(iwr,nbl)

call read_eos_table_micro(iwr,nbl)

call get_diff_rules()

call init_ipl_rule()

call get_eos_report(iwr)

return
end SUBROUTINE init_eos_table
!***********************************************************************
SUBROUTINE define_eos_table(iwr)
! Stefan Typel for the CompOSE core team, version 0.06, 2012/11/01
USE compose_internal
implicit none
integer :: i,j,ierr,iunit,ierror,iwr

ierr = 0
iunit = 20
open(unit=iunit,file='eos.quantities',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   if (iwr == 1) then
      write(*,*)
      write(*,*) 'reading selection of quantities'
   end if

   read(iunit,*)
   read(iunit,*) n_qty,n_add
   if ((n_qty < 0).or.(n_qty > dim_qtyt)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 61
   end if
   if ((n_add < 0).or.(n_qty > dim_qtyt)) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 62
   end if
   if (ierr == 0) then
      read(iunit,*)
      read(iunit,*) (idx_qty(i),i=1,n_qty,1),(idx_add(j),j=1,n_add,1)
      do i=1,n_qty,1
         if ((idx_qty(i) < 0).or.(idx_qty(i) > dim_qtyt)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 63
         end if
      end do
      do i=1,n_add,1
         if ((idx_add(i) < 0).or.(idx_add(i) > nadd_max)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 64
         end if
      end do
      
      read(iunit,*)
      read(iunit,*) n_p,n_q
      if ((n_p < 0).or.(n_p > dim_qtyp)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 65
      end if
      if ((n_q < 0).or.(n_q > dim_qtyq)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 66
      end if
      if (ierr == 0) then
         read(iunit,*)
         read(iunit,*) (idx_p(i),i=1,n_p,1),(idx_q(j),j=1,n_q,1)

         read(iunit,*)
         read(iunit,*) n_m
         if ((n_m < 0).or.(n_m > dim_qtym)) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 67
         end if
         if (ierr == 0) then
            read(iunit,*)
            read(iunit,*) (idx_m(i),i=1,n_m,1)

            read(iunit,*)
            read(iunit,*) n_err
            if ((n_err < 0).or.(n_err > 6)) then
               ierr = ierr+1
               if (ierr < dim_err) error_msg(ierr) = 68
            end if
            if (ierr == 0) then
               read(iunit,*)
               read(iunit,*) (idx_err(i),i=1,n_err,1)
               if ((idx_err(i) < 1).or.(idx_err(i) > 8)) idx_err(i) = i

               read(iunit,*)
               read(iunit,*) iout
#ifndef hdf5
               if ((iout < 1).or.(iout > 2)) then
                  ierr = ierr+1
                  if (ierr < dim_err) error_msg(ierr) = 69
               end if
#endif
            end if
         end if
         
      end if
      if (n_add > nadd_max) n_add = nadd_max
      if (n_p > np_max) n_p = np_max
      if (n_q > nq_max) n_q = nq_max
      if (n_m > nm_max) n_m = nm_max
      
   end if
else
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 60
endif
close(unit=iunit)

if (idx_ex(2) /= 1) then
   n_p = 0
   n_q = 0
end if
if (idx_ex(3) /= 1) n_m = 0

if (iwr == 1) then
   write(*,*) 'number of extracted regular quantities    =',n_qty
   write(*,*) 'indices of regular quantities             :',&
        (idx_qty(j),j=1,n_qty,1)
   write(*,*)
   write(*,*) 'number of extracted additional quantities =',n_add
   write(*,*) 'indices of additional quantities          :',&
        (idx_add(j),j=1,n_add,1)
   write(*,*)
   write(*,*) 'number of pairs for composition           =',n_p
   write(*,*) 'indices of pairs                          :',&
        (idx_p(j),j=1,n_p,1)
   write(*,*)
   write(*,*) 'number of quadruples for composition      =',n_q
   write(*,*) 'indices of quadruples                     :',&
        (idx_q(j),j=1,n_q,1)
   write(*,*)
   write(*,*) 'number of microscopic quantities          =',n_m
   write(*,*) 'indices of microscopic quantities         :',&
        (idx_m(j),j=1,n_m,1)
   write(*,*)
   write(*,*) 'number of error quantities                =',n_err
   write(*,*) 'indices of error quantities               :',&
        (idx_err(j),j=1,n_err,1)
   write(*,*)
   if (iout == 1) then
      write(*,*) 'format of output table                    =',&
           iout,' (ASCII)'
   else
      write(*,*) 'format of output table                    =',&
           iout,' (HDF5)'
   end if
end if

call write_errors(ierr)

return
end SUBROUTINE define_eos_table
!***********************************************************************
SUBROUTINE get_eos_table(iwr)
! Stefan Typel for the CompOSE core team, version 0.16, 2012/12/19
USE compose_internal
#if defined (hdf5)
USE hdfparameters
#endif

implicit none
integer :: i,k,i_tny,i_beta,i_tab,ipl_t,ipl_n,ipl_y,itest,iwr,&
     i1,i2,i3,i4,i5,i6,j_t,j_nb,j_yq,n_t,n_nb,n_yq,i_t,i_nb,i_yq,&
     ierr,iunit,ierror,iunit2,n_tny
double precision :: t_min,t_max,d_t,nb_min,nb_max,d_nb,&
     yq_min,yq_max,d_yq,t,n,y,eos_q(0:dim_qtyq2)

if (iwr == 1) then
   write(*,*)
   write(*,*) 'begin generating eos table'
end if

itest = 0

eos_compo_p(0) = 0.d00
eos_q(0) = 0.d00
eos_micro(0) = 0.d00

ierr = 0
   
if (itest == 0) then
   
   iunit = 20
   open(unit=iunit,file='eos.parameters',&
        status='old',action='read',iostat=ierror)
   if (ierror == 0) then
      if (iwr ==1) then
         write(*,*)
         write(*,*) 'reading parameters'
      end if
      
      read(iunit,*)
      read(iunit,*) ipl_t,ipl_n,ipl_y
      if ((ipl_t < 1).or.(ipl_t > 3)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 70
      end if
      if ((ipl_n < 1).or.(ipl_n > 3)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 71
      end if
      if ((ipl_y < 1).or.(ipl_y > 3)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 72
      end if
      read(iunit,*)
      read(iunit,*) i_beta
      read(iunit,*)
      read(iunit,*) i_tab
      read(iunit,*)
      if (i_tab == 0) then
         read(iunit,*) n_tny
         if (n_tny < 1) then
            ierr = ierr +1
            if (ierr < dim_err) error_msg(ierr) = 80
         end if
      else
         read(iunit,*) t_min,t_max,n_t,i_t
         if (i_t == 0) then
            if (t_min < 0.d00) then
               ierr = ierr +1
               if (ierr < dim_err) error_msg(ierr) = 81
            end if
         else
            if (t_min <= 0.d00) then
               ierr = ierr +1
               if (ierr < dim_err) error_msg(ierr) = 82
            end if
         end if
         if (t_max < t_min) then
            ierr = ierr +1
            if (ierr < dim_err) error_msg(ierr) = 83
         end if
         if (n_t < 1) n_t = 1
         if (i_t == 0) then
            if (n_t == 1) then
               d_t = 0.d00
            else
               d_t = (t_max-t_min)/dble(n_t-1)
            end if
         else
            if (n_t == 1) then
               d_t = 1.d00
            else
               d_t = exp(log(t_max/t_min)/dble(n_t-1))
            end if
         end if
         read(iunit,*) nb_min,nb_max,n_nb,i_nb
         if (nb_min <= 0.d00) then
            ierr = ierr +1
            if (ierr < dim_err) error_msg(ierr) = 84
         end if
         if (nb_max < nb_min) then
            ierr = ierr +1
            if (ierr < dim_err) error_msg(ierr) = 85
         end if
         if (n_nb < 1) n_nb = 1
         if (i_nb == 0) then
            if (n_nb == 1) then
               d_nb = 0.d00
            else
               d_nb = (nb_max-nb_min)/dble(n_nb-1)
            end if
         else
            if (n_nb == 1) then
               d_nb = 1.d00
            else
               d_nb = exp(log(nb_max/nb_min)/dble(n_nb-1))
            end if
         end if
         if (i_beta == 1) then
            yq_min = 0.d00
            yq_max = 0.d00
            d_yq = 0.d00
            n_yq = 1
         else
            read(iunit,*) yq_min,yq_max,n_yq,i_yq
            if (yq_min < 0.d00) then
               ierr = ierr +1
               if (ierr < dim_err) error_msg(ierr) = 86
            end if
            if (yq_max < yq_min) then
               ierr = ierr +1
               if (ierr < dim_err) error_msg(ierr) = 87
            end if
            if (yq_max > 1.d00) then
               ierr = ierr +1
               if (ierr < dim_err) error_msg(ierr) = 88
            end if
            if (n_yq < 1) n_yq = 1
            if (i_yq == 0) then
               if (n_yq == 1) then
                  d_yq = 0.d00
               else
                  d_yq = (yq_max-yq_min)/dble(n_yq-1)
               end if
            else
               if (n_yq == 1) then
                  d_yq = 1.d00
               else
                  d_yq = exp(log(yq_max/yq_min)/dble(n_yq-1))
               end if
            end if
         end if
      end if

      ! read eos.parameters and output in ASCII or HDF5 format
      iunit2 = 21
      if (iout == 1) then
         ! ASCII
         open(unit=iunit2,file='eos.table',&
              status='unknown',action='write',iostat=ierror)
      else
         ! HDF5
#if defined (hdf5)
         IF(i_tab == 0) then
            call initialise_hdf5(n_tny,1,1,i_tab)
         else
            call initialise_hdf5(n_nb,n_t,n_yq,i_tab)
         end IF
#endif 
      end if
      if (ierror == 0) then
         if (ierr == 0) then
            if (i_tab == 0) then
               ! use list of parameters from file eos.parameters
               do i_tny=1,n_tny,1
                  read(iunit,*) t,n,y
                  !+++++++++++++++++++++++++++++++++++++++++++
                  call get_eos(t,n,y,ipl_t,ipl_n,ipl_y,i_beta)
                  !+++++++++++++++++++++++++++++++++++++++++++
                  if (y > -1.d00) then
                     ! no output if no solution of beta equilibrium
                     if (iout == 1) then
                        ! ASCII
                        do i=1,n_q,1
                           eos_q(3*(i-1)+1) = eos_compo_q(i,1)
                           eos_q(3*(i-1)+2) = eos_compo_q(i,2)
                           eos_q(3*(i-1)+3) = eos_compo_q(i,3)
                        end do
                        write(iunit2,*) t,n,y,&
                             (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                             (eos_thermo_add(i2),i2=1,n_add,1),&
                             (eos_compo_p(i3),i3=1,n_p,1),&
                             (eos_q(i4),i4=1,3*n_q,1),&
                             (eos_micro(i5),i5=1,n_m,1),&
                             (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)
                     else
                        ! HDF5
#if defined (hdf5)
                        nb_hdf5(i_tny) = n
                        t_hdf5(i_tny) = t
                        y_q_hdf5(i_tny) = y
                        IF(n_qty.ne.0) then
                           thermo_hdf5(i_tny,1,1,1:n_qty) = eos_thermo(idx_qty(1:n_qty))
                           index_thermo(1:n_qty) = idx_qty(1:n_qty)
                        end IF
                        IF(n_add.ne.0) then
                           thermo_hdf5_add(i_tny,1,1,1:n_add) = eos_thermo_add(1:n_add)
                           index_thermo_add(1:n_add) = idx_add(1:n_add)
                        end IF
                        IF(n_p.ne.0) then
                           yi_hdf5(i_tny,1,1,1:n_p) = eos_compo_p(1:n_p)
                           index_yi(1:n_p) = idx_p(1:n_p)
                        end IF
                        IF(n_q.ne.0) then
                           aav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,1)
                           zav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,2)
                           yav_hdf5(i_tny,1,1,1:n_q) = eos_compo_q(1:n_q,3)
                           index_av(1:n_q) = idx_q(1:n_q)
                        end IF
                        IF(n_m.ne.0) then
                           micro_hdf5(i_tny,1,1,1:n_m) = eos_micro(1:n_m)
                           index_micro(1:n_m) = idx_m(1:n_m)
                        end IF
                        IF(n_err.ne.0) then
                           err_hdf5(i_tny,1,1,1:n_err) = eos_thermo_err(idx_err(1:n_err))
                           index_err(1:n_err) = idx_err(1:n_err)
                        end IF
#endif                        
                     end if
                  end if
               end do
            else
               ! use cycle form of parameters from file eos.parameters
               do j_t=1,n_t,1

                  if (i_t == 0) then
                     t = t_min+d_t*dble(j_t-1)
                  else
                     t = t_min*(d_t**(j_t-1))
                  end if
#if defined (hdf5)
                  IF(iout .ne.1) t_hdf5(j_t) = t
#endif
                  do j_nb=1,n_nb,1
                     if (i_nb == 0) then
                        n = nb_min+d_nb*dble(j_nb-1)
                     else
                        n = nb_min*(d_nb**(j_nb-1))
                     end if
#if defined (hdf5)
                     IF(iout.ne.1) nb_hdf5(j_nb) = n
#endif
                     do j_yq=1,n_yq,1
                        if (i_yq == 0) then
                           y = yq_min+d_yq*dble(j_yq-1)
                        else
                           y = yq_min*(d_yq**(j_yq-1))
                        end if
#if defined (hdf5)
                        IF(iout.ne.1) y_q_hdf5(j_yq) = y
#endif

                        !+++++++++++++++++++++++++++++++++++++++++++
                        call get_eos(t,n,y,ipl_t,ipl_n,ipl_y,i_beta)
                        !+++++++++++++++++++++++++++++++++++++++++++
                        if (y > -1.d00) then
                           ! no output if no solution of beta equilibrium
                           if (iout == 1) then
                              !ASCII
                              do i=1,n_q,1
                                 eos_q(3*(i-1)+1) = eos_compo_q(i,1)
                                 eos_q(3*(i-1)+2) = eos_compo_q(i,2)
                                 eos_q(3*(i-1)+3) = eos_compo_q(i,3)
                              end do
                              write(iunit2,*) t,n,y,&
                                   (eos_thermo(idx_qty(i1)),i1=1,n_qty,1),&
                                   (eos_thermo_add(i2),i2=1,n_add,1),&
                                   (eos_compo_p(i3),i3=1,n_p,1),&
                                   (eos_q(i4),i4=1,3*n_q,1),&
                                   (eos_micro(i5),i5=1,n_m,1),&
                                   (eos_thermo_err(idx_err(i6)),i6=1,n_err,1)
!hallo test
!                              write(iunit2,*) t,n,y,&
!                                   eos_thermo(idx_qty(1)),&
!                                   eos_thermo(idx_qty(2))*m_n+0.5d00*(m_n-m_p),&
!                                   eos_thermo(idx_qty(3))*m_n+0.5d00*(m_n-m_p)
!                                  eos_thermo(idx_qty(1)),&
!                                  eos_thermo(idx_qty(2)),&
!                                  eos_thermo(idx_qty(3))
                           else
                              ! HDF5
#if defined (hdf5)
                              nb_hdf5(j_nb) = n

                              y_q_hdf5(j_yq) = y
                              IF(n_qty.ne.0) then
                                 thermo_hdf5(j_nb,j_t,j_yq,1:n_qty) = eos_thermo(idx_qty(1:n_qty))
                                 index_thermo(1:n_qty) = idx_qty(1:n_qty)
                              end IF
                              IF(n_add.ne.0) then
                                 thermo_hdf5_add(j_nb,j_t,j_yq,1:n_add) = eos_thermo_add(1:n_add)
                                 index_thermo_add(1:n_add) = idx_add(1:n_add)
                              end IF
                              IF(n_p.ne.0) then
                                 yi_hdf5(j_nb,j_t,j_yq,1:n_p) = eos_compo_p(1:n_p)
                                 index_yi(1:n_p) = idx_p(1:n_p)
                              end IF
                              IF(n_q.ne.0) then
                                 aav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,1)
                                 zav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,2)
                                 yav_hdf5(j_nb,j_t,j_yq,1:n_q) = eos_compo_q(1:n_q,3)
                                 index_av(1:n_q) = idx_q(1:n_q)
                              end IF
                              IF(n_m.ne.0) then
                                 micro_hdf5(j_nb,j_t,j_yq,1:n_m) = eos_micro(1:n_m)
                                 index_micro(1:n_m) = idx_m(1:n_m)
                              end IF
                              IF(n_err.ne.0) then
                                 err_hdf5(j_nb,j_t,j_yq,1:n_err) = eos_thermo_err(idx_err(1:n_err))
                                 index_err(1:n_err) = idx_err(1:n_err)
                              end IF
#endif                        
                           end if
                        end if
                     end do
                  end do
               end do
            end if
         end if
      end if
   end if
   if (iout == 1) then
      ! ASCII
      close(unit=iunit2)
   else
      ! HDF5
#if defined (hdf5)
      call write_hdf5(i_tab)
      call close_hdf5
! testing purposes
      call read_hdf5(i_tab)
      call close_hdf5
!****************

#endif
   end if
   close(unit=iunit)   

else
   ! only for testing
   ! input:
   ! temperature T in MeV
   t = 1.d00
   ! baryon number density n_b in fm^-3
   n = 0.01d00
   ! baryonic charge fraction Y_q (dimensionless)
   y = 0.45d00
   ! interpolation rules
   ! 1: first order with continuity of function at grid points
   ! 2: second order with continuity of function and first derivative
   ! 3: third order with continuity of function, first and second derivative
   ! interpolation rule for temperature T
   ipl_t = 3
   ! interpolation rule for baryon number density n_b
   ipl_n = 3
   ! interpolation rule for baryonic charge fraction y_q
   ipl_y = 3
   ! calculation of eos
   do i=1,1000,1
      ! i = 0
      ! t = 0.1d00+0.005d00*i
      t = 0.5+0.001*i
      ! n = 0.0001d00*i
      ! n = 0.01d00+0.00001d00*i
      ! n = 0.01001d00
      ! y = 0.01d00+0.0005d00*i
      ! y = 0.35d00+0.0001d00*i
      !++++++++++++++++++++++++++++++++++++++
      call get_eos(t,n,y,ipl_t,ipl_n,ipl_y,0)
      !++++++++++++++++++++++++++++++++++++++
      ! write(30,*) t,n,y,eos_thermo_add(1)
      write(30,*) t,n,y,(eos_thermo(k),k=1,7,1)
      ! write(30,*) t,n,y,(eos_thermo(k),k=1,7,1)
      ! write(31,*) t,n,y,(idx_add(k),eos_thermo_add(k),k=1,n_add,1)
      ! write(32,*) t,n,y,(eos_thermo_err(2*k-1),k=1,3,1)
      ! write(33,*) t,n,y,(eos_thermo_err(2*k),k=1,3,1)
      ! if (n_p > 0)&
      ! write(34,*) t,n,y,(idx_p(k),eos_compo_p(k),k=1,n_p,1)
      ! if (n_q > 0)&
      ! write(35,*) t,n,y,(idx_q(k),&
      ! eos_compo_q(k,1),eos_compo_q(k,2),eos_compo_q(k,3),k=1,n_q,1)
      ! if (n_m > 0)&
      ! write(36,*) t,n,y,(idx_m(k),eos_micro(k),k=1,n_m,1)
   end do

end if

call write_errors(ierr)

if (iwr == 1) then
   write(*,*)
   write(*,*) 'end generating eos table'
   write(*,*)
end if

return
end SUBROUTINE get_eos_table
!***********************************************************************
SUBROUTINE get_eos(t,n,y,ipl_t,ipl_n,ipl_y,i_beta)
! Stefan Typel for the CompOSE core team, version 0.09, 2012/12/17
USE compose_internal
implicit none
integer :: ierr,ipl_t,ipl_n,ipl_y,i_beta
double precision :: t,n,y

! error counter
ierr = 0

! determination of beta equilibrium conditions not supported yet
if (i_beta == 1) then
! determination of beta equilibrium
   if ((eos_type == 1).or.(eos_type == 2)) then
      if (eos_type == 2) t = 0.d00
      call get_eos_beta(t,n,y,ipl_t,ipl_n,ipl_y,ierr)
   else
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 90
   end if
else
   if ((eos_type == 2).or.(eos_type > 5)) t = 0.d00
   if ((eos_type == 3).or.(eos_type == 6)) y = 0.5d00
   if ((eos_type == 4).or.(eos_type > 7)) y = 0.d00
   call get_eos_sub(t,n,y,ipl_t,ipl_n,ipl_y,ierr,0)
end if

call write_errors(ierr)

return
end SUBROUTINE get_eos
!***********************************************************************
SUBROUTINE get_eos_beta(t,n,y,ipl_t,ipl_n,ipl_y,ierr)
! Stefan Typel for the CompOSE core team, version 0.08, 2012/11/02
USE compose_internal
implicit none
integer :: ipl_t,ipl_n,ipl_y,ierr
double precision :: t,n,y,y_min,y_max,f_min,f_max,f,d_f,d_y

! error counter
ierr = 0

! boundaries
y_min = tab_para(idx_min(3),3)*1.00000001d00
call get_eos_sub(t,n,y_min,ipl_t,ipl_n,ipl_y,ierr,0)
f_min = eos_thermo(5)
!write(*,*) y_min,f_min

y_max = tab_para(idx_max(3),3)*0.99999999d00
call get_eos_sub(t,n,y_max,ipl_t,ipl_n,ipl_y,ierr,0)
f_max = eos_thermo(5)
!write(*,*) y_max,f_max

if ((f_min*f_max > 0.d00).or.(y_min > y_max)) then
! no beta equilibrium
   y = -2.d00
   return
else
! determination of zero
   y = 0.d00
   d_y = 1.d00
   do while (dabs(d_y) > 1.d-12)
      d_y = y
      d_f = f_max-f_min
      if (dabs(d_f) > 1.d-08) then
         y = (f_max*y_min-f_min*y_max)/d_f
      else
         y = 0.5d00*(y_max+y_min)
      end if
      call get_eos_sub(t,n,y,ipl_t,ipl_n,ipl_y,ierr,0)
      f = eos_thermo(5)
      d_y = d_y-y
      !      write(*,*) y,f,d_y
      if (f_min*f > 0.d00) then
         y_min = y
         f_min = f
      else
         y_max = y
         f_max = f
      end if
   end do
endif

call get_eos_sub(t,n,y,ipl_t,ipl_n,ipl_y,ierr,0)

return
end SUBROUTINE get_eos_beta
!***********************************************************************
SUBROUTINE get_eos_sub(t,n,y,ipl_t,ipl_n,ipl_y,ierr,inmp)
! Stefan Typel for the CompOSE core team, version 0.16, 2012/12/18
USE compose_internal
implicit none
integer :: ierr,inmp,ipl_t,ipl_n,ipl_y,m(3),ip,ic,i1,i2,i3,ipl(3)
double precision :: t,n,y,x(3),q(3)

! parameters
x(1) = t
x(2) = n
x(3) = y

! interpolation rules
ipl(1) = ipl_t
ipl(2) = ipl_n
ipl(3) = ipl_y
do ip=1,3,1
   if ((ipl(ip) < 1).or.(ipl(ip) > 3)) ipl(ip) = 3
end do

! grid parameters
do ip = 1,3,1
   if (idx_ipl(ip) == 1) call get_eos_grid_para(ip,x,m,q,ierr)
end do
if ((eos_type == 2).or.(eos_type > 5)) then
! T = 0 EoS
   m(1) = 0
   q(1) = 0.d00
end if
if ((eos_type == 3).or.(eos_type == 6)) then
! symmetric nuclear matter EoS
   m(3) = idx_min(3)
   q(3) = 0.d00
end if
if ((eos_type == 4).or.(eos_type == 7)) then
! neutron matter matter EoS
   m(3) = idx_min(3)
   q(3) = 0.d00
end if
if ((eos_type == 5).or.(eos_type == 8)) then
! beta equilibrium EoS
   m(3) = 0
   q(3) = 0.d00
end if

!do ip=1,3,1
!   write(*,*) ip,x(ip),m(ip),q(ip)
!end do

! test for existence of corner points
ic = 1
! 3 dim
if ((idx_ipl(1) == 1).and.(idx_ipl(2) == 1).and.(idx_ipl(3) == 1)) then
   do i1=m(1),m(1)+1,1
      do i2=m(2),m(2)+1,1
         do i3=m(3),m(3)+1,1
            ic = ic*idx_arg(i1,i2,i3,0)
         end do
      end do
   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 53
   end if
end if
! 2 dim
if ((idx_ipl(1) == 0).and.(idx_ipl(2) == 1).and.(idx_ipl(3) == 1)) then
   !   do i1=m(1),m(1)+1,1
   i1 = m(1)
   do i2=m(2),m(2)+1,1
      do i3=m(3),m(3)+1,1
         ic = ic*idx_arg(i1,i2,i3,0)
      end do
   end do
   !   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 52
   end if
end if
if ((idx_ipl(1) == 1).and.(idx_ipl(2) == 0).and.(idx_ipl(3) == 1)) then
   do i1=m(1),m(1)+1,1
      !      do i2=m(2),m(2)+1,1
      i2 = m(2)
      do i3=m(3),m(3)+1,1
         ic = ic*idx_arg(i1,i2,i3,0)
      end do
      !      end do
   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 52
   end if
end if
if ((idx_ipl(1) == 1).and.(idx_ipl(2) == 1).and.(idx_ipl(3) == 0)) then
   do i1=m(1),m(1)+1,1
      do i2=m(2),m(2)+1,1
         !         do i3=m(3),m(3)+1,1
         i3 = m(3)
         ic = ic*idx_arg(i1,i2,i3,0)
         !         end do
      end do
   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 52
   end if
end if
! 1 dim
if ((idx_ipl(1) == 0).and.(idx_ipl(2) == 0).and.(idx_ipl(3) == 1)) then
   !   do i1=m(1),m(1)+1,1
   i1 = m(1)
   !      do i2=m(2),m(2)+1,1
   i2 = m(2)
   do i3=m(3),m(3)+1,1
      ic = ic*idx_arg(i1,i2,i3,0)
   end do
   !      end do
   !   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 51
   end if
end if
if ((idx_ipl(1) == 0).and.(idx_ipl(2) == 1).and.(idx_ipl(3) == 0)) then
   !   do i1=m(1),m(1)+1,1
   i1 = m(1)
   do i2=m(2),m(2)+1,1
      !         do i3=m(3),m(3)+1,1
      i3 = m(3)
      ic = ic*idx_arg(i1,i2,i3,0)
      !         end do
   end do
   !   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 51
   end if
end if
if ((idx_ipl(1) == 1).and.(idx_ipl(2) == 0).and.(idx_ipl(3) == 0)) then
   do i1=m(1),m(1)+1,1
      !      do i2=m(2),m(2)+1,1
      i2 = m(2)
      !         do i3=m(3),m(3)+1,1
      i3 = m(3)
      ic = ic*idx_arg(i1,i2,i3,0)
      !         end do
      !      end do
   end do
   if (ic /= 1) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 51
   end if
end if

if (ierr == 0) call eos_interpol(x,m,q,ipl,inmp)

return
end SUBROUTINE get_eos_sub
!***********************************************************************
SUBROUTINE get_eos_grid_para(ip,x,m,q,ierr)
! Stefan Typel for the CompOSE core team, version 0.05, 2012/11/05
USE compose_internal
implicit none
integer :: ierr,ip,ia,im,m(3)
double precision x(3),q(3)

im = idx_min(ip)
do ia=idx_min(ip),idx_max(ip),1
   if (x(ip) >= (tab_para(ia,ip))) im = ia
end do
if (im > (idx_max(ip)-1)) im = idx_max(ip)-1
m(ip) = im
q(ip) = (x(ip)-tab_para(im,ip))/(tab_para(im+1,ip)-tab_para(im,ip))
if ((q(ip) < -1.d-06).or.((q(ip)-1.d00) > 1.d-06)) then
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 40+ip
end if

return
end SUBROUTINE get_eos_grid_para
!***********************************************************************
SUBROUTINE eos_interpol(x,m,q,ipl,inmp)
! Stefan Typel for the CompOSE core team, version 0.19, 2016/06/15
USE compose_internal
implicit none
integer :: m(3),ipl(3),inmp,i1,i2,i3,i1p,i2p,i3p,iq,iq2,is,&
     igp,ngp,igp3,ngp3,igp_r,igp3_r,&
     idx1(100),idx2(100),idx3(10),ir(-4:5),irp,&
     idx_thermo2(dim_qty2),&
     get_ipl_rule
double precision :: x(3),q(3),vec(-4:5),vec3(-4:5,3),&
     mat(dim_ipl,-4:5,-4:5),mat3(dim_ipl,-4:5,-4:5,3),&
     qx,qy,dg(0:2,0:2),v_thermo(dim_qty2,0:3),tmp,dh(0:2),&
     vp_compo(dim_qtyp),vq_compo(dim_qtyq,3),v_micro(dim_qtym)

! list of grid points and reference point in T and n_b
igp = 0
igp_r = 0
do i1=-4,5,1 ! T
   i1p = i1+m(1)
   if ((i1p >= idx_min(1)).and.(i1p <= idx_max(1))) then
      do i2=-4,5,1 ! n_b
         i2p = i2+m(2)
         if ((i2p >= idx_min(2)).and.(i2p <= idx_max(2))) then
            igp = igp+1
            idx1(igp) = i1
            idx2(igp) = i2
            if (q(1) <= 0.5d00) then
               if (q(2) <= 0.5d00) then
                  if ((i1 == 0).and.(i2 == 0)) igp_r = igp
               else
                  if ((i1 == 0).and.(i2 == 1)) igp_r = igp
               end if
            else
               if (q(2) <= 0.5d00) then
                  if ((i1 == 1).and.(i2 == 0)) igp_r = igp
               else
                  if ((i1 == 1).and.(i2 == 1)) igp_r = igp
               end if
            end if
         end if
      end do
   end if
end do
ngp = igp
! list of grid points and reference point in Y_q
if (idx_ipl(3) == 1) then
   igp3 = 0
   igp3_r = 0
   do i3=-4,5,1 ! Y_q
      i3p = i3+m(3)
      if ((i3p >= idx_min(3)).and.(i3p <= idx_max(3))) then
         igp3 = igp3+1
         idx3(igp3) = i3
         if (q(3) <= 0.5d00) then
            if (i3 == 0) igp3_r = igp3
         else
            if (i3 == 1) igp3_r = igp3
         end if
      end if
   end do
   ngp3 = igp3
end if

do iq=1,(dim_reg+nadd_max),1
   idx_thermo2(iq) = iq
end do
do iq=(dim_reg+1+nadd_max),dim_qty2,1
   idx_thermo2(iq) = 0
end do

if ((eos_type == 2).or.(eos_type > 5)) then
   ! for interpolation in n_b
   call get_idx_arg1(m)
else
   ! for interpolation in T and n_b
   call get_idx_arg2(m)
   call get_diff_rules2(m,ipl)
end if

! eos_thermo
! direct interpolation in Y_q if needed
do iq=1,dim_qty2,1
   do i1=-4,5,1
      do i2=-4,5,1
         mat(iq,i1,i2) = 0.d00
      end do
   end do
   if (idx_thermo2(iq) > 0) then
      if (idx_ipl(3) == 1) then
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            do i3=-4,5,1
               vec(i3) = 0.d00
            end do

            do igp3=1,ngp3,1
               i3p = idx3(igp3)+m(3)
               vec(idx3(igp3)) = tab_thermo(i1p,i2p,i3p,iq)
               irp = idx_arg(i1p,i2p,i3p,3)
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(3))
            end do
            call get_interpol_yq(m,q,ir,vec,dh,ipl)
            mat(iq,idx1(igp),idx2(igp)) = dh(0)
         end do
      else
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            i3p = m(3)
            mat(iq,idx1(igp),idx2(igp)) = tab_thermo(i1p,i2p,i3p,iq)
         end do
      end if
   end if
end do

if ((eos_type == 2).or.(eos_type > 5)) then
   ! direct interpolation in n_b for T = 0
   do iq=1,dim_qty2,1
      if (idx_thermo2(iq) > 0) then
         do i3=-4,5,1
            vec(i3) = 0.d00
         end do
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            if (i1p /= 0) then
               write(*,*) 'error'
               stop
            end if
            vec(idx2(igp)) = mat(iq,idx1(igp),idx2(igp))
            irp = idx_arg1(idx2(igp),2)
            ir(idx2(igp)) = get_ipl_rule(irp,idx2(igp),ipl(2))
         end do
         call get_interpol_nb(m,q,ir,vec,dh,ipl)
         v_thermo(iq,0) = dh(0)
         v_thermo(iq,1) = dh(1)
         v_thermo(iq,2) = 0.d00
         v_thermo(iq,3) = dh(2)
      else 
         v_thermo(iq,0) = 0.d00
         v_thermo(iq,1) = 0.d00
         v_thermo(iq,2) = 0.d00
         v_thermo(iq,3) = 0.d00
      end if
   end do
else
   ! direct interpolation in T and n_b
   qx = q(1)
   qy = q(2)
   do iq=1,dim_qty2,1
      if (idx_thermo2(iq) > 0) then
         do i1=-4,5,1
            do i2=-4,5,1
               df(0,0,i1,i2) = mat(iq,i1,i2)
            end do
         end do
         call get_derivatives(ipl)
         call get_coefficients()
         call get_interpol_xy(qx,qy,dg,1)
         v_thermo(iq,0) = dg(0,0)
         v_thermo(iq,1) = dg(0,1)
         v_thermo(iq,2) = dg(1,0)
         v_thermo(iq,3) = dg(0,2)
      else
         v_thermo(iq,0) = 0.d00
         v_thermo(iq,1) = 0.d00
         v_thermo(iq,2) = 0.d00
         v_thermo(iq,3) = 0.d00
      end if
   end do
end if

! p pressure
eos_thermo(1) = v_thermo(1,0)*x(2)
! S entropy ber baryon
eos_thermo(2) = v_thermo(2,0)
! mu_b-m_n baryon number chemical potential
eos_thermo(3) = v_thermo(3,0)*m_n
! mu_q charge chemical potential
eos_thermo(4) = v_thermo(4,0)*m_n
! mu_l lepton number chemical potential
eos_thermo(5) = v_thermo(5,0)*m_n
! F/m_n-1 free energy per baryon
eos_thermo(6) = v_thermo(6,0)
! E/m_n-1 internal energy per baryon
eos_thermo(7) = v_thermo(7,0)
!eos_thermo(7) = v_thermo(6,0)+x(1)*eos_thermo(2)/m_n

if (inmp == 1) then
   eos_thermo_add(1) = v_thermo(1,1)
   eos_thermo_add(2) = v_thermo(1,3)
   return
endif

! delta (f/n) = f/n+p/n-(mu_b+y_q*(mu_q+mu_e), mu_e = mu_l-mu_q

if (incl_l == 1) then
   ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_l)
   eos_thermo_err(1) = eos_thermo(6)*m_n+v_thermo(1,0)&
        -(eos_thermo(3)+x(3)*eos_thermo(5))
   !   write(*,*) eos_thermo(6)*m_n,v_thermo(1,0)&
   !        -(eos_thermo(3)+x(3)*eos_thermo(5))
else 
   ! delta (f/n) = f/n+p/n-(mu_b+y_q*mu_q)
   eos_thermo_err(1) = eos_thermo(6)*m_n+v_thermo(1,0)&
        -(eos_thermo(3)+x(3)*eos_thermo(4))
   !write(*,*) eos_thermo(6)*m_n,v_thermo(1,0)&
   !     -(eos_thermo(3)+x(3)*eos_thermo(4))
end if

! (delta f/n)/(f/n)
if (eos_thermo(6) /= -1.d00) then
   eos_thermo_err(2) = eos_thermo_err(1)/((eos_thermo(6)+1.d00)*m_n)
else
   eos_thermo_err(2) = 0.d00
end if

if (incl_l == 0) then
! delta (e/n) = f/n+p/n-(mu_b+y_q*mu_q)-Ts
   eos_thermo_err(3) = eos_thermo(7)*m_n+v_thermo(1,0)&
        -(eos_thermo(3)+x(3)*eos_thermo(4))&
        -x(1)*eos_thermo(2)
else
! delta (e/n) = e/n+p/n-(mu_b+y_q*mu_l)-Ts
   eos_thermo_err(3) = eos_thermo(7)*m_n+v_thermo(1,0)&
        -(eos_thermo(3)+x(3)*eos_thermo(5))&
        -x(1)*eos_thermo(2)
end if

! (delta e/n)/(e/n)
if (eos_thermo(7) /= -1.d00) then
   eos_thermo_err(4) = eos_thermo_err(3)/((eos_thermo(7)+1.d00)*m_n)
else
   eos_thermo_err(4) = 0.d00
end if

! delta (p/n) = p/n-n*d(f/n)/dn
eos_thermo_err(5) = v_thermo(1,0)-x(2)*v_thermo(6,1)*m_n
if (v_thermo(1,0) /= 0.d00) then
   eos_thermo_err(6) = eos_thermo_err(5)/v_thermo(1,0)
else
   eos_thermo_err(6) = 0.d00
end if

! delta (s/n) = s/n+d(f/n)/dt
eos_thermo_err(7) = eos_thermo(2)+v_thermo(6,2)*m_n
if (eos_thermo(2) /= 0.d00) then
   eos_thermo_err(8) = eos_thermo_err(7)/eos_thermo(2)
else
   eos_thermo_err(8) = 0.d00
end if

tmp = v_thermo(1,0)/m_n
! H/m_n-1  enthalpy per baryon
eos_thermo(8) = eos_thermo(7)+tmp
! G/m_n-1 free enthalpy per baryon
eos_thermo(9) = eos_thermo(6)+tmp
!do iq = 10,dim_qty2,1
do iq = 10,dim_qtyt,1
  eos_thermo(iq) = 0.d00
end do

! beta_v
eos_thermo(17) = v_thermo(1,2)*x(2)
! kappa_t
eos_thermo(18) = x(2)*(x(2)*v_thermo(1,1)+v_thermo(1,0))
if (eos_thermo(18) > 0.d00) then
   eos_thermo(18) = 1.d00/eos_thermo(18)
else
   eos_thermo(18) = 0.d00
end if
! alpha_p
eos_thermo(16) = eos_thermo(17)*eos_thermo(18)
! c_v
eos_thermo(13) = v_thermo(2,2)*x(1)
! c_p
eos_thermo(14) = eos_thermo(13)+eos_thermo(16)*eos_thermo(17)*x(1)/x(2)
! gamma
if (eos_thermo(13) > 0.d00) then
   eos_thermo(15) = eos_thermo(14)/eos_thermo(13)
else
   eos_thermo(16) = 1.d00
end if
! kappa_s
eos_thermo(19) = eos_thermo(18)/eos_thermo(15)
! dp/de|n
eos_thermo(11) = x(1)*v_thermo(2,2)
if (eos_thermo(11) > 0.d00) then
   eos_thermo(11) = x(2)*v_thermo(1,2)/eos_thermo(11)
else
   eos_thermo(11) = 0.d00
end if
! cs^2 and dp/dn|e
eos_thermo(12) = eos_thermo(19)*(eos_thermo(8)+1.d00)*m_n*x(2)
if (eos_thermo(12) > 0.d00) then
   eos_thermo(12) = 1.d00/eos_thermo(12)
   eos_thermo(10) = (eos_thermo(8)+1.d00)*m_n*eos_thermo(12)&
        -v_thermo(1,0)*eos_thermo(11)/x(2)
else
   eos_thermo(12) = 0.d00
   eos_thermo(10) = 0.d00
end if

! additional quantities
do iq=1,nadd_max,1
   eos_thermo_add(iq) = v_thermo(iq+dim_reg,0)
end do
!endif

! eos_compo
do iq=1,dim_qtyp,1
   idx_compo_p(iq) = 0
   eos_compo_p(iq) = 0.d00
end do
do iq=1,dim_qtyq,1
   idx_compo_q(iq) = 0
   eos_compo_q(iq,1) = 0.d00
   eos_compo_q(iq,2) = 0.d00
   eos_compo_q(iq,3) = 0.d00
end do
!if ((idx_ex(2) == 1).and.(inmp /= 1)) then
if (idx_ex(2) == 1) then
   ! pairs
   ! direct interpolation in Y_q if needed
   do iq=1,np_max,1
      do i1=-4,5,1
         do i2=-4,5,1
            mat(iq,i1,i2) = 0.d00
         end do
      end do
   end do
   
   do iq=1,np_max,1
      if (idx_ipl(3) == 1) then
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            do i3=-4,5,1
               vec(i3) = 0.d00
            end do
            do igp3=1,ngp3,1
               i3p = idx3(igp3)+m(3)
               vec(idx3(igp3)) = 0.d00
               do iq2=1,np_max,1
                  if (idxp_compo(i1p,i2p,i3p,iq2) == idx_p(iq)) then
                     vec(idx3(igp3)) = tabp_compo(i1p,i2p,i3p,iq2)
                  end if
               end do
               irp = idx_arg(i1p,i2p,i3p,3)
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(3))
            end do
            call get_interpol_yq(m,q,ir,vec,dh,ipl)
            mat(iq,idx1(igp),idx2(igp)) = dh(0)
         end do
      else
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            i3p = m(3)
            mat(iq,idx1(igp),idx2(igp)) = 0.d00
            do iq2=1,np_max,1
               if (idxp_compo(i1p,i2p,i3p,iq2) == idx_p(iq)) then
                  mat(iq,idx1(igp),idx2(igp)) = tabp_compo(i1p,i2p,i3p,iq2)
               end if
            end do
         end do
      end if
   end do

   if ((eos_type == 2).or.(eos_type > 5)) then
      ! direct interpolation in n_b for T = 0
      do iq=1,np_max,1
         do i3=-4,5,1
            vec(i3) = 0.d00
         end do
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            if (i1p /= 0) then
               write(*,*) 'error'
               stop
            end if
            vec(idx2(igp)) = mat(iq,idx1(igp),idx2(igp))
            irp = idx_arg1(idx2(igp),2)
            ir(idx2(igp)) = get_ipl_rule(irp,idx2(igp),ipl(2))
         end do
         call get_interpol_nb(m,q,ir,vec,dh,ipl)
         vp_compo(iq) = dh(0)
      end do
   else
      ! direct interpolation in T and n_b
      qx = q(1)
      qy = q(2)
      do iq=1,np_max,1
         do i1=-4,5,1
            do i2=-4,5,1
               df(0,0,i1,i2) = mat(iq,i1,i2)
            end do
         end do
         call get_derivatives(ipl)
         call get_coefficients()
         call get_interpol_xy(qx,qy,dg,0)
         vp_compo(iq) = dg(0,0)
      end do
   end if
   do iq=1,np_max,1
      eos_compo_p(iq) = vp_compo(iq)
      idx_compo_p(iq) = idx_p(iq)
   end do

   ! quadruples
   ! direct interpolation in Y_q if needed
   do iq=1,nq_max,1
      do i1=-4,5,1
         do i2=-4,5,1
            mat(iq,i1,i2) = 0.d00
            do is=1,3,1
               mat3(iq,i1,i2,is) = 0.d00
            end do
         end do
      end do
   end do
   
   do iq=1,nq_max,1
      if (idx_ipl(3) == 1) then
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            do i3=-4,5,1
               vec(i3) = 0.d00
               vec3(i3,1) = 0.d00
               vec3(i3,2) = 0.d00
               vec3(i3,3) = 0.d00
            end do
            do igp3=1,ngp3,1
               i3p = idx3(igp3)+m(3)
               do is=1,3,1
                  vec3(idx3(igp3),is) = 0.d00
               end do
               do iq2=1,nq_max,1
                  if (idxq_compo(i1p,i2p,i3p,iq2) == idx_q(iq)) then
                     do is=1,3,1
                        vec3(idx3(igp3),is) =&
                             tabq_compo(i1p,i2p,i3p,3*(iq2-1)+is)
                     end do
                  end if
               end do
               irp = idx_arg(i1p,i2p,i3p,3)
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(3))
            end do
            do is=1,3,1
               do i3=-4,5,1
                  vec(i3) = vec3(i3,is)
               end do
               call get_interpol_yq(m,q,ir,vec,dh,ipl)
               mat3(iq,idx1(igp),idx2(igp),is) = dh(0)
            end do
         end do
      else
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            i3p = m(3)
            do is=1,3,1
               mat3(iq,idx1(igp),idx2(igp),is) = 0.d00
            end do
            do iq2=1,nq_max,1
               if (idxp_compo(i1p,i2p,i3p,iq2) == idx_q(iq)) then
                  do is=1,3,1
                     mat3(iq,idx1(igp),idx2(igp),is) =&
                          tabq_compo(i1p,i2p,i3p,3*(iq2-1)+is)
                  end do
               end if
            end do
         end do
      end if
   end do

   if ((eos_type == 2).or.(eos_type > 5)) then
      ! direct interpolation in n_b for T = 0
      do iq=1,nq_max,1
         do i3=-4,5,1
            vec(i3) = 0.d00
            do is=1,3,1
               vec3(i3,is) = 0.d00
            end do
         end do
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            if (i1p /= 0) then
               write(*,*) 'error'
               stop
            end if
            do is=1,3,1
               vec3(idx2(igp),is) = mat3(iq,idx1(igp),idx2(igp),is)
            end do
            irp = idx_arg1(idx2(igp),2)
            ir(idx2(igp)) = get_ipl_rule(irp,idx2(igp),ipl(2))
         end do
         do is=1,3,1
            do i3=-4,5,1
               vec(i3) = vec3(i3,is)
            end do
            call get_interpol_nb(m,q,ir,vec,dh,ipl)
            vq_compo(iq,is) = dh(0)
         end do
      end do
   else
      ! direct interpolation in T and n_b
      qx = q(1)
      qy = q(2)
      do iq=1,nq_max,1
         do is=1,3,1
            do i1=-4,5,1
               do i2=-4,5,1
                  df(0,0,i1,i2) = mat3(iq,i1,i2,is)
               end do
            end do
            call get_derivatives(ipl)
            call get_coefficients()
            call get_interpol_xy(qx,qy,dg,0)
            vq_compo(iq,is) = dg(0,0)
         end do
      end do
   end if
   do iq=1,nq_max,1
      do is=1,3,1
         eos_compo_q(iq,is) = vq_compo(iq,is)
      end do
      idx_compo_q(iq) = idx_q(iq)
   end do
end if

! eos_micro
do iq=1,dim_qtym,1
   idx_micro(iq) = 0
   eos_micro(iq) = 0.d00
end do
!if ((idx_ex(3) == 1).and.(inmp /= 1)) then
if (idx_ex(3) == 1) then
   ! direct interpolation in Y_q if needed
   do iq=1,nm_max,1
      do i1=-4,5,1
         do i2=-4,5,1
            mat(iq,i1,i2) = 0.d00
         end do
      end do
   end do
   
   do iq=1,nm_max,1
      if (idx_ipl(3) == 1) then
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            do i3=-4,5,1
               vec(i3) = 0.d00
            end do
            do igp3=1,ngp3,1
               i3p = idx3(igp3)+m(3)
               vec(idx3(igp3)) = 0.d00
               do iq2=1,nm_max,1
                  if (idx_mic(i1p,i2p,i3p,iq2) == idx_m(iq)) then
                     vec(idx3(igp3)) = tab_mic(i1p,i2p,i3p,iq2)
                  end if
               end do
               irp = idx_arg(i1p,i2p,i3p,3)
               ir(idx3(igp3)) = get_ipl_rule(irp,idx3(igp3),ipl(3))
            end do
            call get_interpol_yq(m,q,ir,vec,dh,ipl)
            mat(iq,idx1(igp),idx2(igp)) = dh(0)
         end do
      else
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            i3p = m(3)
            mat(iq,idx1(igp),idx2(igp)) = 0.d00
            do iq2=1,nm_max,1
               if (idx_mic(i1p,i2p,i3p,iq2) == idx_m(iq)) then
                  mat(iq,idx1(igp),idx2(igp)) = tab_mic(i1p,i2p,i3p,iq2)
               end if
            end do
         end do
      end if
   end do

   if ((eos_type == 2).or.(eos_type > 5)) then
      ! direct interpolation in n_b for T = 0
      do iq=1,nm_max,1
         do i3=-4,5,1
            vec(i3) = 0.d00
         end do
         do igp=1,ngp,1
            i1p = idx1(igp)+m(1)
            i2p = idx2(igp)+m(2)
            if (i1p /= 0) then
               write(*,*) 'error'
               stop
            end if
            vec(idx2(igp)) = mat(iq,idx1(igp),idx2(igp))
            irp = idx_arg1(idx2(igp),2)
            ir(idx2(igp)) = get_ipl_rule(irp,idx2(igp),ipl(2))
         end do
         call get_interpol_nb(m,q,ir,vec,dh,ipl)
         v_micro(iq) = dh(0)
      end do
   else
      ! direct interpolation in T and n_b
      qx = q(1)
      qy = q(2)
      do iq=1,nm_max,1
         do i1=-4,5,1
            do i2=-4,5,1
               df(0,0,i1,i2) = mat(iq,i1,i2)
            end do
         end do
         call get_derivatives(ipl)
         call get_coefficients()
         call get_interpol_xy(qx,qy,dg,0)
         v_micro(iq) = dg(0,0)
      end do
   end if
   do iq=1,nm_max,1
      eos_micro(iq) = v_micro(iq)
      idx_micro(iq) = idx_m(iq)
   end do
end if

return
end SUBROUTINE eos_interpol
!***********************************************************************
SUBROUTINE init_ipl_rule()
! Stefan Typel for the CompOSE core team, version 0.07, 2012/10/04
USE compose_internal
implicit none
integer :: ir

! third order
do ir=-8,8,1
   ipl_rule(0,3,ir) = ir
   ipl_rule(1,3,ir) = ir
end do

! second order
do ir=-8,8,1
   ipl_rule(0,2,ir) = 6
end do
ipl_rule(0,2, 3) =  7
ipl_rule(0,2,-3) = -7
ipl_rule(0,2, 5) =  7
ipl_rule(0,2,-5) = -7
ipl_rule(0,2, 7) =  7
ipl_rule(0,2,-7) = -7
ipl_rule(0,2, 8) =  8
ipl_rule(0,2,-8) = -8
do ir=-8,8,1
   ipl_rule(1,2,ir) = ipl_rule(0,2,ir)
end do

! first order
do ir=-8,8,1
   ipl_rule(0,1,ir) = -8
   ipl_rule(1,1,ir) = 8
end do
ipl_rule(0,1,-3) = -8
ipl_rule(0,1,-5) = -8
ipl_rule(0,1,-7) = -8
ipl_rule(0,1,-8) = -8
ipl_rule(1,1, 3) =  8
ipl_rule(1,1, 5) =  8
ipl_rule(1,1, 7) =  8
ipl_rule(1,1, 8) =  8

return
end SUBROUTINE init_ipl_rule
!***********************************************************************
SUBROUTINE get_interpol_xy(qx,qy,dg,order)
! Stefan Typel for the CompOSE core team, version 0.02, 2012/10/05
USE compose_internal
implicit none
integer :: order,ix,iy
double precision :: qx,qy,dg(0:2,0:2),xx(-2:5),yy(-2:5)

do ix=-2,-1,1
   xx(ix) = 0.d00
   yy(ix) = 0.d00
end do
xx(0) = 1.d00
yy(0) = 1.d00
do ix=1,5,1
   xx(ix) = xx(ix-1)*qx
   yy(ix) = yy(ix-1)*qy
end do

do ix=0,2,1
   do iy=0,2,1
      dg(ix,iy) = 0.d00
   end do
end do

if (order == 0) then
   ! function
   do ix=0,5,1
      do iy=0,5,1
         dg(0,0) = dg(0,0)+fc(ix,iy)*xx(ix)*yy(iy)
      end do
   end do
else
   ! function and derivatives
   do ix=0,5,1
      do iy=0,5,1
         dg(0,0) = dg(0,0)+fc(ix,iy)*xx(ix)*yy(iy)
         dg(1,0) = dg(1,0)+fc(ix,iy)*xx(ix-1)*yy(iy)*dble(ix)
         dg(0,1) = dg(0,1)+fc(ix,iy)*xx(ix)*yy(iy-1)*dble(iy)
         dg(2,0) = dg(2,0)+fc(ix,iy)*xx(ix-2)*yy(iy)*dble(ix*(ix-1))
         dg(1,1) = dg(1,1)+fc(ix,iy)*xx(ix-1)*yy(iy-1)*dble(iy*ix)
         dg(0,2) = dg(0,2)+fc(ix,iy)*xx(ix)*yy(iy-2)*dble(iy*(iy-1))
         dg(2,1) = dg(2,1)+fc(ix,iy)*xx(ix-2)*yy(iy-1)*dble(ix*(ix-1)*iy)
         dg(1,2) = dg(1,2)+fc(ix,iy)*xx(ix-1)*yy(iy-2)*dble(ix*iy*(iy-1))
         dg(2,2) = dg(2,2)+fc(ix,iy)*xx(ix-2)*yy(iy-2)*&
              dble(ix*(ix-1)*iy*(iy-1))
      end do
   end do
   
   ! rescaling
   if (dx > 0.d00) then
      dg(1,0) = dg(1,0)/dx
      dg(2,0) = dg(2,0)/dx**2
   end if
   if (dy > 0.d00) then
      dg(0,1) = dg(0,1)/dy
      dg(0,2) = dg(0,2)/dy**2
   end if
   if ((dx > 0.d00).and.(dy > 0.)) then
      dg(1,1) = dg(1,1)/(dx*dy)
      dg(2,1) = dg(2,1)/((dx**2)*dy)
      dg(1,2) = dg(1,2)/(dx*(dy**2))
      dg(2,2) = dg(2,2)/(dx*dy)**2
   end if
end if

return
end SUBROUTINE get_interpol_xy
!***********************************************************************
SUBROUTINE get_coefficients()
! Stefan Typel for the CompOSE core team, version 0.03, 2012/10/05
USE compose_internal
implicit none

fc(0,0) = df(0,0,0,0)
fc(0,1) = df(0,1,0,0)
fc(0,2) = 0.5d00*df(0,2,0,0)
fc(0,3) = -10.d00*(df(0,0,0,0)-df(0,0,0,1))&
     -6.d00*df(0,1,0,0)-4.d00*df(0,1,0,1)&
     -1.5d00*df(0,2,0,0)+0.5d00*df(0,2,0,1)
fc(0,4) =  15.d00*(df(0,0,0,0)-df(0,0,0,1))&
     +8.d00*df(0,1,0,0)+7.d00*df(0,1,0,1)&
     +1.5d00*df(0,2,0,0)-1.d00*df(0,2,0,1)
fc(0,5) =  -6.d00*(df(0,0,0,0)-df(0,0,0,1))&
     -3.d00*(df(0,1,0,0)+df(0,1,0,1))&
     -0.5d00*(df(0,2,0,0)-df(0,2,0,1))
fc(1,0) = df(1,0,0,0)
fc(1,1) = df(1,1,0,0)
fc(1,2) = 0.5d00*df(1,2,0,0)
fc(1,3) = -10.d00*(df(1,0,0,0)-df(1,0,0,1))&
     -6.d00*df(1,1,0,0)-4.d00*df(1,1,0,1)&
     -1.5d00*df(1,2,0,0)+0.5d00*df(1,2,0,1)
fc(1,4) =  15.d00*(df(1,0,0,0)-df(1,0,0,1))&
     +8.d00*df(1,1,0,0)+7.d00*df(1,1,0,1)&
     +1.5d00*df(1,2,0,0)-1.d00*df(1,2,0,1)
fc(1,5) =  -6.d00*(df(1,0,0,0)-df(1,0,0,1))&
     -3.d00*(df(1,1,0,0)+df(1,1,0,1))&
     -0.5d00*(df(1,2,0,0)-df(1,2,0,1))
fc(2,0) = 0.5d00*df(2,0,0,0)
fc(2,1) = 0.5d00*df(2,1,0,0)
fc(2,2) = 0.25d00*df(2,2,0,0)
fc(2,3) = -5.d00*(df(2,0,0,0)-df(2,0,0,1))&
     -3.d00*df(2,1,0,0)-2.d00*df(2,1,0,1)&
     -0.75d00*df(2,2,0,0)+0.25d00*df(2,2,0,1)
fc(2,4) =  7.5d00*(df(2,0,0,0)-df(2,0,0,1))&
     +4.d00*df(2,1,0,0)+3.5d00*df(2,1,0,1)&
     +0.75d00*df(2,2,0,0)-0.5d00*df(2,2,0,1)
fc(2,5) = -3.d00*(df(2,0,0,0)-df(2,0,0,1))&
     -1.5d00*(df(2,1,0,0)+df(2,1,0,1))&
     -0.25d00*(df(2,2,0,0)-df(2,2,0,1))
fc(3,0) =  -10.d00*(df(0,0,0,0)-df(0,0,1,0))&
     -6.d00*df(1,0,0,0)-4.d00*df(1,0,1,0)&
     -1.5d00*df(2,0,0,0)+0.5d00*df(2,0,1,0)
fc(3,1) = -10.d00*(df(0,1,0,0)-df(0,1,1,0))&
     -6.d00*df(1,1,0,0)-4.d00*df(1,1,1,0)&
     -1.5d00*df(2,1,0,0)+0.5d00*df(2,1,1,0)
fc(3,2) = -5.d00*(df(0,2,0,0)-df(0,2,1,0))&
     -3.d00*df(1,2,0,0)-2.d00*df(1,2,1,0)&
     -0.75d00*df(2,2,0,0)+0.25d00*df(2,2,1,0)
fc(3,3) =  100.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     +60.d00*(df(1,0,0,0)-df(1,0,0,1))+40.d00*(df(1,0,1,0)-df(1,0,1,1))&
     +60.d00*(df(0,1,0,0)-df(0,1,1,0))+40.d00*(df(0,1,0,1)-df(0,1,1,1))&
     +15.d00*(df(2,0,0,0)-df(2,0,0,1))-5.d00*(df(2,0,1,0)-df(2,0,1,1))&
     +36.d00*df(1,1,0,0)+24.d00*(df(1,1,1,0)+df(1,1,0,1))+16.d00*df(1,1,1,1)&
     +15.d00*(df(0,2,0,0)-df(0,2,1,0))-5.d00*(df(0,2,0,1)-df(0,2,1,1))&
     +9.d00*(df(2,1,0,0)+df(1,2,0,0))-3.d00*(df(2,1,1,0)+df(1,2,0,1))&
     +6.d00*(df(2,1,0,1)+df(1,2,1,0))-2.d00*(df(2,1,1,1)+df(1,2,1,1))&
     +2.25d00*df(2,2,0,0)-0.75d00*(df(2,2,1,0)+df(2,2,0,1))&
     +0.25d00*df(2,2,1,1)
fc(3,4) = -150.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     -90.d00*(df(1,0,0,0)-df(1,0,0,1))-60.d00*(df(1,0,1,0)-df(1,0,1,1))&
     -80.d00*(df(0,1,0,0)-df(0,1,1,0))-70.d00*(df(0,1,0,1)-df(0,1,1,1))&
     -22.5d00*(df(2,0,0,0)-df(2,0,0,1))+7.5d00*(df(2,0,1,0)-df(2,0,1,1))&
     -48.d00*df(1,1,0,0)-32.d00*df(1,1,1,0)&
     -42.d00*df(1,1,0,1)-28.d00*df(1,1,1,1)&
     -15.d00*(df(0,2,0,0)-df(0,2,1,0))+10.d00*(df(0,2,0,1)-df(0,2,1,1))&
     -12.d00*df(2,1,0,0)+4.d00*df(2,1,1,0)&
     -10.5d00*df(2,1,0,1)+3.5d00*df(2,1,1,1)&
     -9.d00*df(1,2,0,0)-6.d00*(df(1,2,1,0)-df(1,2,0,1))+4.d00*df(1,2,1,1)&
     -2.25d00*df(2,2,0,0)+0.75d00*df(2,2,1,0)&
     +1.5d00*df(2,2,0,1)-0.5d00*df(2,2,1,1)
fc(3,5) =   60.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     +36.d00*(df(1,0,0,0)-df(1,0,0,1))+24.d00*(df(1,0,1,0)-df(1,0,1,1))&
     +30.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
     +9.d00*(df(2,0,0,0)-df(2,0,0,1))-3.d00*(df(2,0,1,0)-df(2,0,1,1))&
     +18.d00*(df(1,1,0,0)+df(1,1,0,1))+12.d00*(df(1,1,1,0)+df(1,1,1,1))&
     +5.d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
     +4.5d00*(df(2,1,0,0)+df(2,1,0,1))-1.5d00*(df(2,1,1,0)+df(2,1,1,1))&
     +3.d00*(df(1,2,0,0)-df(1,2,0,1))+2.d00*(df(1,2,1,0)-df(1,2,1,1))&
     +0.75d00*(df(2,2,0,0)-df(2,2,0,1))&
     -0.25d00*(df(2,2,1,0)-df(2,2,1,1))
fc(4,0) =   15.d00*(df(0,0,0,0)-df(0,0,1,0))&
     +8.d00*df(1,0,0,0)+7.d00*df(1,0,1,0)&
     +1.5d00*df(2,0,0,0)-1.d00*df(2,0,1,0)
fc(4,1) = 15.d00*(df(0,1,0,0)-df(0,1,1,0))&
     +8.d00*df(1,1,0,0)+7.d00*df(1,1,1,0)&
     +1.5d00*df(2,1,0,0)-1.d00*df(2,1,1,0)
fc(4,2) = 7.5d00*(df(0,2,0,0)-df(0,2,1,0))&
     +4.d00*df(1,2,0,0)+3.5d00*df(1,2,1,0)&
     +0.75d00*df(2,2,0,0)-0.5d00*df(2,2,1,0)
fc(4,3) = -150.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     -80.d00*(df(1,0,0,0)-df(1,0,0,1))-70.d00*(df(1,0,1,0)-df(1,0,1,1))&
     -90.d00*(df(0,1,0,0)-df(0,1,1,0))-60.d00*(df(0,1,0,1)-df(0,1,1,1))&
     -15.d00*(df(2,0,0,0)-df(2,0,0,1))+10.d00*(df(2,0,1,0)-df(2,0,1,1))&
     -48.d00*df(1,1,0,0)-42.d00*df(1,1,1,0)&
     -32.d00*df(1,1,0,1)-28.d00*df(1,1,1,1)&
     -22.5d00*(df(0,2,0,0)-df(0,2,1,0))+7.5d00*(df(0,2,0,1)-df(0,2,1,1))&
     -9.d00*df(2,1,0,0)+6.d00*(df(2,1,1,0)-df(2,1,0,1))+4.d00*df(2,1,1,1)&
     -12.d00*df(1,2,0,0)-10.5d00*df(1,2,1,0)&
     +4.d00*df(1,2,0,1)+3.5d00*df(1,2,1,1)&
     -2.25d00*df(2,2,0,0)+1.5d00*df(2,2,1,0)&
     +0.75d00*df(2,2,0,1)-0.5d00*df(2,2,1,1)
fc(4,4) =  225.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     +120.d00*(df(1,0,0,0)-df(1,0,0,1))+105.d00*(df(1,0,1,0)-df(1,0,1,1))&
     +120.d00*(df(0,1,0,0)-df(0,1,1,0))+105.d00*(df(0,1,0,1)-df(0,1,1,1))&
     +22.5d00*(df(2,0,0,0)-df(2,0,0,1))-15.d00*(df(2,0,1,0)-df(2,0,1,1))&
     +64.d00*df(1,1,0,0)+56.d00*(df(1,1,1,0)+df(1,1,0,1))+49.d00*df(1,1,1,1)&
     +22.5d00*(df(0,2,0,0)-df(0,2,1,0))-15.d00*(df(0,2,0,1)-df(0,2,1,1))&
     +12.d00*(df(2,1,0,0)+df(1,2,0,0))-8.d00*(df(2,1,1,0)+df(1,2,0,1))&
     +10.5d00*(df(2,1,0,1)+df(1,2,1,0))-7.d00*(df(2,1,1,1)+df(1,2,1,1))&
     +2.25d00*df(2,2,0,0)-1.5d00*(df(2,2,1,0)+df(2,2,0,1))&
     +1.d00*df(2,2,1,1)
fc(4,5) =  -90.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     -48.d00*(df(1,0,0,0)-df(1,0,0,1))-42.d00*(df(1,0,1,0)-df(1,0,1,1))&
     -45.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
     -9.d00*(df(2,0,0,0)-df(2,0,0,1))+6.d00*(df(2,0,1,0)-df(2,0,1,1))&
     -24.d00*(df(1,1,0,0)+df(1,1,0,1))-21.d00*(df(1,1,1,0)+df(1,1,1,1))&
     -7.5d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
     -4.5d00*(df(2,1,0,0)+df(2,1,0,1))+3.d00*(df(2,1,1,0)+df(2,1,1,1))&
     -4.d00*df(1,2,0,0)-3.5d00*(df(1,2,1,0)-df(1,2,1,1))+4.d00*df(1,2,0,1)&
     -0.75d00*(df(2,2,0,0)-df(2,2,0,1))&
     +0.5d00*(df(2,2,1,0)-df(2,2,1,1))
fc(5,0) =   -6.d00*(df(0,0,0,0)-df(0,0,1,0))&
     -3.d00*(df(1,0,0,0)+df(1,0,1,0))&
     -0.5d00*(df(2,0,0,0)-df(2,0,1,0))
fc(5,1) = -6.d00*(df(0,1,0,0)-df(0,1,1,0))&
     -3.d00*(df(1,1,0,0)+df(1,1,1,0))&
     -0.5d00*(df(2,1,0,0)-df(2,1,1,0))
fc(5,2) = -3.d00*df(0,2,0,0)+3.d00*df(0,2,1,0)&
     -1.5d00*df(1,2,0,0)-1.5d00*df(1,2,1,0)&
     -0.25d00*df(2,2,0,0)+0.25d00*df(2,2,1,0)
fc(5,3) =   60.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     +30.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
     +36.d00*(df(0,1,0,0)-df(0,1,1,0))+24.d00*(df(0,1,0,1)-df(0,1,1,1))&
     +5.d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
     +18.d00*(df(1,1,0,0)+df(1,1,1,0))+12.d00*(df(1,1,0,1)+df(1,1,1,1))&
     +9.d00*(df(0,2,0,0)-df(0,2,1,0))-3.d00*(df(0,2,0,1)-df(0,2,1,1))&
     +3.d00*(df(2,1,0,0)-df(2,1,1,0))+2.d00*(df(2,1,0,1)-df(2,1,1,1))&
     +4.5d00*df(1,2,0,0)+4.5d00*df(1,2,1,0)&
     -1.5d00*(df(1,2,0,1)+df(1,2,1,1))&
     +0.75d00*(df(2,2,0,0)-df(2,2,1,0))&
     -0.25d00*(df(2,2,0,1)-df(2,2,1,1))
fc(5,4) =  -90.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     -45.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
     -48.d00*(df(0,1,0,0)-df(0,1,1,0))-42.d00*(df(0,1,0,1)-df(0,1,1,1))&
     -7.5d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
     -24.d00*(df(1,1,0,0)+df(1,1,1,0))-21.d00*(df(1,1,0,1)+df(1,1,1,1))&
     -9.d00*(df(0,2,0,0)-df(0,2,1,0))+6.d00*(df(0,2,0,1)-df(0,2,1,1))&
     -4.d00*(df(2,1,0,0)-df(2,1,1,0))-3.5d00*(df(2,1,0,1)-df(2,1,1,1))&
     -4.5d00*(df(1,2,0,0)+df(1,2,1,0))+3.d00*(df(1,2,0,1)+df(1,2,1,1))&
     -0.75d00*(df(2,2,0,0)-df(2,2,1,0))&
     +0.5d00*(df(2,2,0,1)-df(2,2,1,1))
fc(5,5) =   36.d00*(df(0,0,0,0)-df(0,0,1,0)-df(0,0,0,1)+df(0,0,1,1))&
     +18.d00*(df(1,0,0,0)+df(1,0,1,0)-df(1,0,0,1)-df(1,0,1,1))&
     +18.d00*(df(0,1,0,0)-df(0,1,1,0)+df(0,1,0,1)-df(0,1,1,1))&
     +3.d00*(df(2,0,0,0)-df(2,0,1,0)-df(2,0,0,1)+df(2,0,1,1))&
     +9.d00*(df(1,1,0,0)+df(1,1,1,0)+df(1,1,0,1)+df(1,1,1,1))&
     +3.d00*(df(0,2,0,0)-df(0,2,1,0)-df(0,2,0,1)+df(0,2,1,1))&
     +1.5d00*(df(2,1,0,0)-df(2,1,1,0)+df(2,1,0,1)-df(2,1,1,1))&
     +1.5d00*(df(1,2,0,0)+df(1,2,1,0)-df(1,2,0,1)-df(1,2,1,1))&
     +0.25d00*(df(2,2,0,0)-df(2,2,1,0)-df(2,2,0,1)+df(2,2,1,1))

return
end SUBROUTINE get_coefficients
!***********************************************************************
SUBROUTINE get_derivatives(ipl)
! Stefan Typel for the CompOSE core team, version 0.02, 2012/10/05
USE compose_internal
implicit none
integer :: i1,i2,ix,iy,i1p,i2p,iq,ipl(3)
double precision :: tmp

do i1=-4,5,1
   do i2=-4,5,1
      do ix=1,2,1
         do iy=0,2,1
            df(ix,iy,i1,i2) = 0.d00
         end do
      end do
      do ix=0,2,1
         do iy=1,2,1
            df(ix,iy,i1,i2) = 0.d00
         end do
      end do
   end do
end do

do i1=-4,5,1
   do i2=-4,5,1
      ! first derivatives
      ! x
      do iq=-4,4,1
         i1p = i1+iq
         if ((i1p >= -4).and.(i1p <= 5)) then
            df(1,0,i1,i2) = df(1,0,i1,i2)+df(0,0,i1p,i2)*d1x(i1,i2,iq)
         end if
      end do
      ! y
      do iq=-4,4,1
         i2p = i2+iq
         if ((i2p >= -4).and.(i2p <= 5)) then
            df(0,1,i1,i2) = df(0,1,i1,i2)+df(0,0,i1,i2p)*d1y(i1,i2,iq)
         end if
      end do
      ! second derivatives
      ! xx
      do iq=-4,4,1
         i1p = i1+iq
         if ((i1p >= -4).and.(i1p <= 5)) then
            df(2,0,i1,i2) = df(2,0,i1,i2)+df(0,0,i1p,i2)*d2x(i1,i2,iq)
         end if
      end do
      ! yy
      do iq=-4,4,1
         i2p = i2+iq
         if ((i2p >= -4).and.(i2p <= 5)) then
            df(0,2,i1,i2) = df(0,2,i1,i2)+df(0,0,i1,i2p)*d2y(i1,i2,iq)
         end if
      end do
   end do
end do
! mixed derivatives
do i1=0,1,1
   do i2=0,1,1
      ! second derivatives
      ! xy
      do iq=-4,4,1
         i1p = i1+iq
         if ((i1p >= -4).and.(i1p <= 5)) then
            df(1,1,i1,i2) = df(1,1,i1,i2)+df(0,1,i1p,i2)*d1x(i1,i2,iq)
         end if
      end do
      do iq=-4,4,1
         i2p = i2+iq
         if ((i2p >= -4).and.(i2p <= 5)) then
            df(1,1,i1,i2) = df(1,1,i1,i2)+df(1,0,i1,i2p)*d1y(i1,i2,iq)
         end if
      end do
      df(1,1,i1,i2) = 0.5d00*df(1,1,i1,i2)
      ! third derivatives
      ! xxy
      do iq=-4,4,1
         i2p = i2+iq
         if ((i2p >= -4).and.(i2p <= 5)) then
            df(2,1,i1,i2) = df(2,1,i1,i2)+df(2,0,i1,i2p)*d1y(i1,i2,iq)
         end if
      end do
      ! xyy
      do iq=-4,4,1
         i1p = i1+iq
         if ((i1p >= -4).and.(i1p <= 5)) then
            df(1,2,i1,i2) = df(1,2,i1,i2)+df(0,2,i1p,i2)*d1x(i1,i2,iq)
         end if
      end do
      ! fourth derivative
      ! xxyy
      do iq=-4,4,1
         i1p = i1+iq
         if ((i1p >= -4).and.(i1p <= 5)) then
            df(2,2,i1,i2) = df(2,2,i1,i2)+df(0,2,i1p,i2)*d2x(i1,i2,iq)
         end if
      end do
      do iq=-4,4,1
         i2p = i2+iq
         if ((i2p >= -4).and.(i2p <= 5)) then
            df(2,2,i1,i2) = df(2,2,i1,i2)+df(2,0,i1,i2p)*d2y(i1,i2,iq)
         end if
      end do
      df(2,2,i1,i2) = 0.5d00*df(2,2,i1,i2)
   end do
end do

! rescaling
do i1=0,1,1
   do i2=0,1,1
      df(1,0,i1,i2) = df(1,0,i1,i2)*dx
      df(0,1,i1,i2) = df(0,1,i1,i2)*dy
      df(2,0,i1,i2) = df(2,0,i1,i2)*dx2
      df(0,2,i1,i2) = df(0,2,i1,i2)*dy2
      df(1,1,i1,i2) = df(1,1,i1,i2)*dx*dy
      df(2,1,i1,i2) = df(2,1,i1,i2)*dx2*dy
      df(1,2,i1,i2) = df(1,2,i1,i2)*dx*dy2
      df(2,2,i1,i2) = df(2,2,i1,i2)*dx2*dy2
   end do
end do

if (ipl(2) == 2) then
   do i2=0,2,1
      do i1=0,1,1
         tmp = 6.d00*(df(i2,0,i1,1)-df(i2,0,i1,0))
         df(i2,2,i1,0) = tmp-2.d00*df(i2,1,i1,1)-4.d00*df(i2,1,i1,0)
         df(i2,2,i1,1) = 4.d00*df(i2,1,i1,1)+2.d00*df(i2,1,i1,0)-tmp
      end do
   end do
end if

if (ipl(1) == 2) then
   do i2=0,2,1
      do i1=0,1,1
         tmp = 6.d00*(df(0,i2,1,i1)-df(0,i2,0,i1))
         df(2,i2,0,i1) = tmp-2.d00*df(1,i2,1,i1)-4.d00*df(1,i2,0,i1)
         df(2,i2,1,i1) = 4.d00*df(1,i2,1,i1)+2.d00*df(1,i2,0,i1)-tmp
      end do
   end do
end if

if (ipl(2) == 1) then
   do i1=0,1,1
      do i2=0,2,1
         tmp = df(i2,0,i1,1)-df(i2,0,i1,0)
         df(i2,1,i1,0) = tmp
         df(i2,1,i1,1) = tmp
         df(i2,2,i1,0) = 0.d00
         df(i2,2,i1,1) = 0.d00
      end do
   end do
end if

if (ipl(1) == 1) then
   do i1=0,1,1
      do i2=0,2,1
         tmp = df(0,i2,1,i1)-df(0,i2,0,i1)
         df(1,i2,0,i1) = tmp
         df(1,i2,1,i1) = tmp
         df(2,i2,0,i1) = 0.d00
         df(2,i2,1,i1) = 0.d00
      end do
   end do
end if

return
end SUBROUTINE get_derivatives
!***********************************************************************
SUBROUTINE get_diff_rules2(m,ipl)
! Stefan Typel for the CompOSE core team, version 0.02, 2012/04/05
USE compose_internal
implicit none
integer :: m(3),ipl(3),i1,i2,i3,ix,iy,ixp,iyp,i1p,i2p,iq,ir,irp,&
     get_ipl_rule

do i1=-4,5,1
   do i2=-4,5,1
      do i3=-4,4,1
         d1x(i1,i2,i3) = 0.d00
         d2x(i1,i2,i3) = 0.d00
         d1y(i1,i2,i3) = 0.d00
         d2y(i1,i2,i3) = 0.d00
      end do
   end do
end do

ix = 1
ixp = ix+3

iy = 2
iyp = iy+3

! x
do i1=-4,5,1
   i1p = i1+m(ix)
   if ((i1p >= idx_min(ix)).and.(i1p <= idx_max(ix))) then
      do i2=-4,5,1
         irp = idx_arg2(i1,i2,ix)
         ir = get_ipl_rule(irp,i1,ipl(ix))
         do iq=-4,4,1
            d1x(i1,i2,iq) = r1d(ix,i1p,iq,ir)
            d2x(i1,i2,iq) = r2d(ix,i1p,iq,ir)
         end do
      end do
   end if
end do
dx = tab_para(m(ix),ixp)
dx2 = dx*dx

! y
do i2=-4,5,1
   i2p = i2+m(iy)
   if ((i2p >= idx_min(iy)).and.(i2p <= idx_max(iy))) then
      do i1=-4,5,1
         irp = idx_arg2(i1,i2,iy)
         ir = get_ipl_rule(irp,i2,ipl(iy))
         do iq=-4,4,1
            d1y(i1,i2,iq) = r1d(iy,i2p,iq,ir)
            d2y(i1,i2,iq) = r2d(iy,i2p,iq,ir)
         end do
      end do
   end if
end do
dy = tab_para(m(iy),iyp)
dy2 = dy*dy

return
end SUBROUTINE get_diff_rules2
!***********************************************************************
SUBROUTINE get_idx_arg1(m)
! Stefan Typel for the CompOSE core team, version 0.01, 2012/03/27
USE compose_internal
implicit none
integer :: m(3),i2,ip,i2p,varg(-4:5)

! initialisation
do i2=-4,5,1 ! n_b
   do ip=0,2,1
      idx_arg1(i2,ip) = 0
   end do
end do

if (idx_ipl(3) == 1) then
   do i2=-4,5,1 ! n_b
      i2p = i2+m(2)
      if ((i2p >= idx_min(2)).and.(i2p <= idx_max(2))) then
         if ((idx_arg(m(1),i2p,m(3),0) == 1).and.&
              (idx_arg(m(1),i2p,m(3)+1,0) == 1)) then
            idx_arg1(i2,0) = 1
         end if
      end if
   end do
else
   do i2=-4,5,1 ! n_b
      i2p = i2+m(2)
      if ((i2p >= idx_min(2)).and.(i2p <= idx_max(2))) then
         if (idx_arg(m(1),i2p,m(3),0) == 1) then
            idx_arg1(i2,0) = 1
         end if
      end if
   end do
end if

! x
do i2=-4,4,1
   if ((idx_arg1(i2,0) == 1).and.&
        (idx_arg1(i2+1,0) == 1)) then
      varg(i2) = 1
   else
      varg(i2) = 0
   end if
   varg(5) = 0
end do
do i2=-4,4,1
   if ((varg(i2) == 1).and.(varg(i2+1) >= 1)) then
      varg(i2) = 2
   end if
end do
do i2=-4,3,1
   if ((varg(i2) == 2).and.(varg(i2+2) >= 1)) then
      varg(i2) = 3
   end if
end do
do i2=-4,2,1
   if ((varg(i2) == 3).and.(varg(i2+3) >= 1)) then
      varg(i2) = 4
   end if
end do
do i2=-4,5,1
   if (varg(i2) == 4) then
      idx_arg1(i2+4,2) = 3
      idx_arg1(i2+3,2) = 2
      idx_arg1(i2+2,2) = 1
      if (idx_arg1(i2+1,2) == 0) then
         idx_arg1(i2+1,2) = -2
      end if
      if (idx_arg1(i2  ,2) == 0) then
         idx_arg1(i2  ,2) = -3
      end if
   else
      if (varg(i2) == 3) then
         if (idx_arg1(i2+3,2) == 0) then
            idx_arg1(i2+3,2) =  5
         end if
         if (idx_arg1(i2+2,2) == 0) then
            idx_arg1(i2+2,2) =  4
         end if
         if (idx_arg1(i2+1,2) == 0) then
            idx_arg1(i2+1,2) = -4
         end if
         if (idx_arg1(i2  ,2) == 0) then
            idx_arg1(i2  ,2) = -5
         end if
      else
         if (varg(i2) == 2) then
            if (idx_arg1(i2+2,2) == 0) then
               idx_arg1(i2+2,2) =  7
            end if
            if (idx_arg1(i2+1,2) == 0) then
               idx_arg1(i2+1,2) =  6
            end if
            if (idx_arg1(i2  ,2) == 0) then
               idx_arg1(i2  ,2) = -7
            end if
         else
            if (varg(i2) == 1) then
               if (idx_arg1(i2+1,2) == 0) then
                  idx_arg1(i2+1,2) =  8
               end if
               if (idx_arg1(i2  ,2) == 0) then
                  idx_arg1(i2+1,2) = -8
               end if
            end if
         end if
      end if
   end if
end do

return
end SUBROUTINE get_idx_arg1
!***********************************************************************
SUBROUTINE get_idx_arg2(m)
! Stefan Typel for the CompOSE core team, version 0.02, 2012/03/27
USE compose_internal
implicit none
integer :: m(3),i1,i2,ip,i1p,i2p,varg(-4:5)

! initialisation
do i1=-4,5,1 ! T
   do i2=-4,5,1 ! n_b
      do ip=0,2,1
         idx_arg2(i1,i2,ip) = 0
      end do
   end do
end do

if (idx_ipl(3) == 1) then
   do i1=-4,5,1 ! T
      i1p = i1+m(1)
      if ((i1p >= idx_min(1)).and.(i1p <= idx_max(1))) then
         do i2=-4,5,1 ! n_b
            i2p = i2+m(2)
            if ((i2p >= idx_min(2)).and.(i2p <= idx_max(2))) then
               if ((idx_arg(i1p,i2p,m(3),0) == 1).and.&
                    (idx_arg(i1p,i2p,m(3)+1,0)== 1)) then
                  idx_arg2(i1,i2,0) = 1
               end if
            end if
         end do
      end if
   end do
else
   do i1=-4,5,1 ! T
      i1p = i1+m(1)
      if ((i1p >= idx_min(1)).and.(i1p <= idx_max(1))) then
         do i2=-4,5,1 ! n_b
            i2p = i2+m(2)
            if ((i2p >= idx_min(2)).and.(i2p <= idx_max(2))) then
               if (idx_arg(i1p,i2p,m(3),0) == 1) then
                  idx_arg2(i1,i2,0) = 1
               end if
            end if
         end do
      end if
   end do
end if

! x
do i2=-4,5,1
   do i1=-4,4,1
      if ((idx_arg2(i1,i2,0) == 1).and.&
           (idx_arg2(i1+1,i2,0) == 1)) then
         varg(i1) = 1
      else
         varg(i1) = 0
      end if
      varg(5) = 0
   end do
   do i1=-4,4,1
      if ((varg(i1) == 1).and.(varg(i1+1) >= 1)) then
         varg(i1) = 2
      end if
   end do
   do i1=-4,3,1
      if ((varg(i1) == 2).and.(varg(i1+2) >= 1)) then
         varg(i1) = 3
      end if
   end do
   do i1=-4,2,1
      if ((varg(i1) == 3).and.(varg(i1+3) >= 1)) then
         varg(i1) = 4
      end if
   end do
   do i1=-4,5,1
      if (varg(i1) == 4) then
         idx_arg2(i1+4,i2,1) = 3
         idx_arg2(i1+3,i2,1) = 2
         idx_arg2(i1+2,i2,1) = 1
         if (idx_arg2(i1+1,i2,1) == 0) then
            idx_arg2(i1+1,i2,1) = -2
         end if
         if (idx_arg2(i1  ,i2,1) == 0) then
            idx_arg2(i1  ,i2,1) = -3
         end if
      else
         if (varg(i1) == 3) then
            if (idx_arg2(i1+3,i2,1) == 0) then
               idx_arg2(i1+3,i2,1) =  5
            end if
            if (idx_arg2(i1+2,i2,1) == 0) then
               idx_arg2(i1+2,i2,1) =  4
            end if
            if (idx_arg2(i1+1,i2,1) == 0) then
               idx_arg2(i1+1,i2,1) = -4
            end if
            if (idx_arg2(i1  ,i2,1) == 0) then
               idx_arg2(i1  ,i2,1) = -5
            end if
         else
            if (varg(i1) == 2) then
               if (idx_arg2(i1+2,i2,1) == 0) then
                  idx_arg2(i1+2,i2,1) =  7
               end if
               if (idx_arg2(i1+1,i2,1) == 0) then
                  idx_arg2(i1+1,i2,1) =  6
               end if
               if (idx_arg2(i1  ,i2,1) == 0) then
                  idx_arg2(i1  ,i2,1) = -7
               end if
            else
               if (varg(i1) == 1) then
                  if (idx_arg2(i1+1,i2,1) == 0) then
                     idx_arg2(i1+1,i2,1) =  8
                  end if
                  if (idx_arg2(i1  ,i2,1) == 0) then
                     idx_arg2(i1  ,i2,1) = -8
                  end if
               end if
            end if
         end if
      end if
   end do
end do

! y
do i1=-4,5,1
   do i2=-4,4,1
      if ((idx_arg2(i1,i2,0) == 1).and.&
           (idx_arg2(i1,i2+1,0) == 1)) then
         varg(i2) = 1
      else
         varg(i2) = 0
      end if
      varg(5) = 0
   end do
   do i2=-4,4,1
      if ((varg(i2) == 1).and.(varg(i2+1) >= 1)) then
         varg(i2) = 2
      end if
   end do
   do i2=-4,3,1
      if ((varg(i2) == 2).and.(varg(i2+2) >= 1)) then
         varg(i2) = 3
      end if
   end do
   do i2=-4,2,1
      if ((varg(i2) == 3).and.(varg(i2+3) >= 1)) then
         varg(i2) = 4
      end if
   end do
   do i2=-4,5,1
      if (varg(i2) == 4) then
         idx_arg2(i1,i2+4,2) = 3
         idx_arg2(i1,i2+3,2) = 2
         idx_arg2(i1,i2+2,2) = 1
         if (idx_arg2(i1,i2+1,2) == 0) then
            idx_arg2(i1,i2+1,2) = -2
         end if
         if (idx_arg2(i1,i2  ,2) == 0) then
            idx_arg2(i1,i2  ,2) = -3
         end if
      else
         if (varg(i2) == 3) then
            if (idx_arg2(i1,i2+3,2) == 0) then
               idx_arg2(i1,i2+3,2) =  5
            end if
            if (idx_arg2(i1,i2+2,2) == 0) then
               idx_arg2(i1,i2+2,2) =  4
            end if
            if (idx_arg2(i1,i2+1,2) == 0) then
               idx_arg2(i1,i2+1,2) = -4
            end if
            if (idx_arg2(i1,i2  ,2) == 0) then
               idx_arg2(i1,i2  ,2) = -5
            end if
         else
            if (varg(i2) == 2) then
               if (idx_arg2(i1,i2+2,2) == 0) then
                  idx_arg2(i1,i2+2,2) =  7
               end if
               if (idx_arg2(i1,i2+1,2) == 0) then
                  idx_arg2(i1,i2+1,2) =  6
               end if
               if (idx_arg2(i1,i2  ,2) == 0) then
                  idx_arg2(i1,i2  ,2) = -7
               end if
            else
               if (varg(i2) == 1) then
                  if (idx_arg2(i1,i2+1,2) == 0) then
                     idx_arg2(i1,i2+1,2) =  8
                  end if
                  if (idx_arg2(i1,i2  ,2) == 0) then
                     idx_arg2(i1,i2+1,2) = -8
                  end if
               end if
            end if
         end if
      end if
   end do
end do

return
end SUBROUTINE get_idx_arg2
!***********************************************************************
SUBROUTINE get_interpol_yq(m,q,ir,f,dh,ipl)
! Stefan Typel for the CompOSE core team, version 0.05, 2012/10/05
USE compose_internal
implicit none
integer :: m(3),i3,i3p,iq,ir(-4:5),ipl(3)
double precision :: q(3),f(-4:5),fx(0:1),fxx(0:1),ffc(0:5),dh(0:2),&
     tmp0,tmp1,tmp2,tmp3

! derivatives
do i3=0,1,1
   fx(i3)  = 0.d00
   fxx(i3) = 0.d00
   i3p = i3+m(3)
   do iq=-4,4,1
      fx(i3)  = fx(i3) +f(i3+iq)*r1d(3,i3p,iq,ir(i3))
      fxx(i3) = fxx(i3)+f(i3+iq)*r2d(3,i3p,iq,ir(i3))
   end do
   ! rescaling
   tmp0 = tab_para(m(3),6)
   fx(i3) = fx(i3)*tmp0
   fxx(i3) = fxx(i3)*tmp0*tmp0
end do
if (ipl(3) == 2) then
   tmp1 = 6.d00*(f(1)-f(0))
   fxx(0) = tmp1-2.d00*fx(1)-4.d00*fx(0)
   fxx(1) = 4.d00*fx(1)+2.d00*fx(0)-tmp1
end if
if (ipl(3) == 1) then
   tmp1 = f(1)-f(0)
   fx(0) = tmp1
   fx(1) = tmp1
   fxx(0) = 0.d00
   fxx(1) = 0.d00
end if

! coefficients
ffc(0) = f(0)
ffc(1) = fx(0)
ffc(2) = 0.5d00*fxx(0)
tmp1 = f(1)-f(0)-fx(0)-0.5d00*fxx(0)
tmp2 = fx(1)-fx(0)-fxx(0)
tmp3 = fxx(1)-fxx(0)
ffc(3) =  10.d00*tmp1-4.d00*tmp2+0.5d00*tmp3
ffc(4) = -15.d00*tmp1+7.d00*tmp2-tmp3
ffc(5) =   6.d00*tmp1-3.d00*tmp2+0.5d00*tmp3

! interpolation
! function
dh(0) = ffc(0)+q(3)*(ffc(1)+q(3)*(ffc(2)&
     +q(3)*(ffc(3)+q(3)*(ffc(4)+q(3)*ffc(5)))))
! first derivative
dh(1) = ffc(1)+q(3)*(2.d00*ffc(2)&
     +q(3)*(3.d00*ffc(3)+q(3)*(4.d00*ffc(4)+q(3)*5.d00*ffc(5))))
! second derivative
dh(2) = 2.d00*ffc(2)&
     +q(3)*(6.d00*ffc(3)+q(3)*(12.d00*ffc(4)+q(3)*20.d00*ffc(5)))
! rescaling
dh(1) = dh(1)/tmp0
dh(2) = dh(2)/(tmp0*tmp0)

return
end SUBROUTINE get_interpol_yq
!***********************************************************************
SUBROUTINE get_interpol_nb(m,q,ir,f,dh,ipl)
! Stefan Typel for the CompOSE core team, version 0.05, 2012/10/05
USE compose_internal
implicit none
integer :: m(3),i2,i2p,iq,ir(-4:5),ipl(3)
double precision :: q(3),f(-4:5),fx(0:1),fxx(0:1),ffc(0:5),dh(0:2),&
     tmp0,tmp1,tmp2,tmp3

! derivatives
do i2=0,1,1
   fx(i2)  = 0.d00
   fxx(i2) = 0.d00
   i2p = i2+m(2)
   do iq=-4,4,1
      fx(i2)  = fx(i2) +f(i2+iq)*r1d(2,i2p,iq,ir(i2))
      fxx(i2) = fxx(i2)+f(i2+iq)*r2d(2,i2p,iq,ir(i2))
   end do
   ! rescaling
   tmp0 = tab_para(m(2),5)
   fx(i2) = fx(i2)*tmp0
   fxx(i2) = fxx(i2)*tmp0*tmp0
end do
if (ipl(2) == 2) then
   tmp1 = 6.d00*(f(1)-f(0))
   fxx(0) = tmp1-2.d00*fx(1)-4.d00*fx(0)
   fxx(1) = 4.d00*fx(1)+2.d00*fx(0)-tmp1
end if
if (ipl(2) == 1) then
   tmp1 = f(1)-f(0)
   fx(0) = tmp1
   fx(1) = tmp1
   fxx(0) = 0.d00
   fxx(1) = 0.d00
end if

! coefficients
ffc(0) = f(0)
ffc(1) = fx(0)
ffc(2) = 0.5d00*fxx(0)
tmp1 = f(1)-f(0)-fx(0)-0.5d00*fxx(0)
tmp2 = fx(1)-fx(0)-fxx(0)
tmp3 = fxx(1)-fxx(0)
ffc(3) =  10.d00*tmp1-4.d00*tmp2+0.5d00*tmp3
ffc(4) = -15.d00*tmp1+7.d00*tmp2-tmp3
ffc(5) =   6.d00*tmp1-3.d00*tmp2+0.5d00*tmp3

! interpolation
! functions
dh(0) = ffc(0)+q(2)*(ffc(1)+q(2)*(ffc(2)&
     +q(2)*(ffc(3)+q(2)*(ffc(4)+q(2)*ffc(5)))))
! fisrt derivative
dh(1) = ffc(1)+q(2)*(2.d00*ffc(2)&
     +q(2)*(3.d00*ffc(3)+q(2)*(4.d00*ffc(4)+q(2)*5.d00*ffc(5))))
! second derivative
dh(2) = 2.d00*ffc(2)&
     +q(2)*(6.d00*ffc(3)+q(2)*(12.d00*ffc(4)+q(2)*20.d00*ffc(5)))
! rescaling
dh(1) = dh(1)/tmp0
dh(2) = dh(2)/(tmp0*tmp0)

return
end SUBROUTINE get_interpol_nb
!***********************************************************************
SUBROUTINE get_diff_rules()
! Stefan Typel for the CompOSE core team, version 0.08, 2012/12/19
USE compose_internal
implicit none
integer :: ip,in,it,iy,ia,ia_min,ia_max,ii,ir,iap,iq,iarg(0:dim_a)
double precision :: z0,zm,zm2,zp,zp2

! choose maximum differentiation rule

! used grid points  -4 -3 -2 -1  0  1  2  3  4
! case        1            *  *  *  *  *
! case        2         *  *  *  *  *
! case       -2               *  *  *  *  *
! case        3      *  *  *  *  *
! case       -3                  *  *  *  *  *
! case        4            *  *  *  *
! case       -4               *  *  *  *
! case        5         *  *  *  *
! case       -5                  *  *  *  *
! case        6               *  *  *
! case        7            *  *  *
! case       -7                  *  *  *
! case        8               *  *
! case       -8                  *  *

do it=idx_min(1),idx_max(1),1
   do in=idx_min(2),idx_max(2),1
      do iy=idx_min(3),idx_max(3),1
         do ip=1,3,1
            idx_argx(it,in,iy,ip) = idx_arg(it,in,iy,ip)
         end do
      end do
   end do
end do

do it=idx_min(1),idx_max(1),1
   do in=idx_min(2),idx_max(2),1
      do iy=idx_min(3),idx_max(3),1
         if ((idx_argx(it,in,iy,1) == idx_ex(1)).and.&
              (idx_argx(it,in,iy,2) == idx_ex(2)).and.&
              (idx_argx(it,in,iy,3) == idx_ex(3))) then
            idx_arg(it,in,iy,0) = 1
         end if
         do ip=1,3,1
            idx_arg(it,in,iy,ip)= 0
         end do
      end do
   end do
end do

! temperature
ip = 1
do in=1,dim_n,1
   do iy=0,dim_y,1
      iarg(0) = 0
!      do it=1,(dim_t-1),1
      do it=0,(dim_t-1),1
         if ((idx_arg(it,in,iy,0) == 1).and.&
              (idx_arg(it+1,in,iy,0) == 1)) then
            iarg(it) = 1
         else
            iarg(it) = 0
         end if
      end do
      iarg(dim_t) = 0
!      do it=1,(dim_t-1),1
      do it=0,(dim_t-1),1
         if ((iarg(it) == 1).and.(iarg(it+1) == 1)) then
            iarg(it) = 2
         end if
      end do
!      do it=1,(dim_t-2),1
      do it=0,(dim_t-2),1
         if ((iarg(it) == 2).and.(iarg(it+2) >= 1)) then
            iarg(it) = 3
         end if
      end do
!      do it=1,(dim_t-3),1
      do it=0,(dim_t-3),1
         if ((iarg(it) == 3).and.(iarg(it+3) >= 1)) then
            iarg(it) = 4
         end if
      end do
!      do it=1,dim_t,1
      do it=0,dim_t,1
         if (iarg(it) == 4) then
            idx_arg(it+4,in,iy,ip) = 3
            idx_arg(it+3,in,iy,ip) = 2
            idx_arg(it+2,in,iy,ip) = 1
            if (idx_arg(it+1,in,iy,ip) == 0) then
               idx_arg(it+1,in,iy,ip) = -2
            end if
            if (idx_arg(it,in,iy,ip) == 0) then
               idx_arg(it,in,iy,ip) = -3
            end if
         else
            if (iarg(it) == 3) then
               if (idx_arg(it+3,in,iy,ip) == 0) then
                  idx_arg(it+3,in,iy,ip) = 5
               end if
               if (idx_arg(it+2,in,iy,ip) == 0) then
                  idx_arg(it+2,in,iy,ip) = 4
               end if
               if (idx_arg(it+1,in,iy,ip) == 0) then
                  idx_arg(it+1,in,iy,ip) = -4
               end if
               if (idx_arg(it,in,iy,ip) == 0) then
                  idx_arg(it,in,iy,ip) = -5
               end if
            else
               if (iarg(it) == 2) then
                  if (idx_arg(it+2,in,iy,ip) == 0) then
                     idx_arg(it+2,in,iy,ip) = 7
                  end if
                  if (idx_arg(it+1,in,iy,ip) == 0) then
                     idx_arg(it+1,in,iy,ip) = 6
                  end if
                  if (idx_arg(it,in,iy,ip) == 0) then
                     idx_arg(it,in,iy,ip) = -7
                  end if
               else
                  if (iarg(it) == 1) then
                     if (idx_arg(it+1,in,iy,ip) == 0) then
                        idx_arg(it+1,in,iy,ip) = 8
                     end if
                     if (idx_arg(it,in,iy,ip) == 0) then
                        idx_arg(it,in,iy,ip) = -8
                     end if
                  end if
               end if
            end if
         end if
      end do
   end do
end do

! baryon number density
ip = 2
do it=0,dim_t,1
   do iy=0,dim_y,1
      do in=1,(dim_n-1),1
         if ((idx_arg(it,in,iy,0) == 1).and.&
              (idx_arg(it,in+1,iy,0) == 1)) then
            iarg(in) = 1
         else
            iarg(in) = 0
         end if
      end do
      iarg(dim_n) = 0
      do in=1,(dim_n-1),1
         if ((iarg(in) == 1).and.(iarg(in+1) == 1)) then
            iarg(in) = 2
         end if
      end do
      do in=1,(dim_n-2),1
         if ((iarg(in) == 2).and.(iarg(in+2) >= 1)) then
            iarg(in) = 3
         end if
      end do
      do in=1,(dim_n-3),1
         if ((iarg(in) == 3).and.(iarg(in+3) >= 1)) then
            iarg(in) = 4
         end if
      end do
      do in=1,dim_n,1
         if (iarg(in) == 4) then
            idx_arg(it,in+4,iy,ip) = 3
            idx_arg(it,in+3,iy,ip) = 2
            idx_arg(it,in+2,iy,ip) = 1
            if (idx_arg(it,in+1,iy,ip) == 0) then
               idx_arg(it,in+1,iy,ip) = -2
            end if
            if (idx_arg(it,in,iy,ip) == 0) then
               idx_arg(it,in,iy,ip) = -3
            end if
         else
            if (iarg(in) == 3) then
               if (idx_arg(it,in+3,iy,ip) == 0) then
                  idx_arg(it,in+3,iy,ip) = 5
               end if
               if (idx_arg(it,in+2,iy,ip) == 0) then
                  idx_arg(it,in+2,iy,ip) = 4
               end if
               if (idx_arg(it,in+1,iy,ip) == 0) then
                  idx_arg(it,in+1,iy,ip) = -4
               end if
               if (idx_arg(it,in,iy,ip) == 0) then
                  idx_arg(it,in,iy,ip) = -5
               end if
            else
               if (iarg(in) == 2) then
                  if (idx_arg(it,in+2,iy,ip) == 0) then
                     idx_arg(it,in+2,iy,ip) = 7
                  end if
                  if (idx_arg(it,in+1,iy,ip) == 0) then
                     idx_arg(it,in+1,iy,ip) = 6
                  end if
                  if (idx_arg(it,in,iy,ip) == 0) then
                     idx_arg(it,in,iy,ip) = -7
                  end if
               else
                  if (iarg(in) == 1) then
                     if (idx_arg(it,in+1,iy,ip) == 0) then
                        idx_arg(it,in+1,iy,ip) = 8
                     end if
                     if (idx_arg(it,in,iy,ip) == 0) then
                        idx_arg(it,in,iy,ip) = -8
                     end if
                  end if
               end if
            end if
         end if
      end do
   end do
end do

! hadronic charge fraction
ip = 3
!do it=1,dim_t,1
do it=0,dim_t,1
   do in=1,dim_n,1
      iarg(0) = 0
!      do iy=1,(dim_y-2),1
      do iy=0,(dim_y-2),1
         if ((idx_arg(it,in,iy,0) == 1).and.&
              (idx_arg(it,in,iy+1,0) == 1)) then
            iarg(iy) = 1
         else
            iarg(iy) = 0
         end if
      end do
      iarg(dim_y-1) = 0
      iarg(dim_y) = 0
!      do iy=1,(dim_y-2),1
      do iy=0,(dim_y-2),1
         if ((iarg(iy) == 1).and.(iarg(iy+1) == 1)) then
            iarg(iy) = 2
         end if
      end do
!      do iy=1,(dim_y-3),1
      do iy=0,(dim_y-3),1
         if ((iarg(iy) == 2).and.(iarg(iy+2) >= 1)) then
            iarg(iy) = 3
         end if
      end do
!      do iy=1,(dim_y-4),1
      do iy=0,(dim_y-4),1
         if ((iarg(iy) == 3).and.(iarg(iy+3) >= 1)) then
            iarg(iy) = 4
         end if
      end do
!      do iy=1,dim_y,1
      do iy=0,(dim_y-4),1
         if (iarg(iy) == 4) then
            idx_arg(it,in,iy+4,ip) = 3
            idx_arg(it,in,iy+3,ip) = 2
            idx_arg(it,in,iy+2,ip) = 1
            if (idx_arg(it,in,iy+1,ip) == 0) then
               idx_arg(it,in,iy+1,ip) = -2
            end if
            if (idx_arg(it,in,iy,ip) == 0) then
               idx_arg(it,in,iy,ip) = -3
            end if
         else
            if (iarg(iy) == 3) then
               if (idx_arg(it,in,iy+3,ip) == 0) then
                  idx_arg(it,in,iy+3,ip) = 5
               end if
               if (idx_arg(it,in,iy+2,ip) == 0) then
                  idx_arg(it,in,iy+2,ip) = 4
               end if
               if (idx_arg(it,in,iy+1,ip) == 0) then
                  idx_arg(it,in,iy+1,ip) = -4
               end if
               if (idx_arg(it,in,iy,ip) == 0) then
                  idx_arg(it,in,iy,ip) = -5
               end if
            else
               if (iarg(iy) == 2) then
                  if (idx_arg(it,in,iy+2,ip) == 0) then
                     idx_arg(it,in,iy+2,ip) = 7
                  end if
                  if (idx_arg(it,in,iy+1,ip) == 0) then
                     idx_arg(it,in,iy+1,ip) = 6
                  end if
                  if (idx_arg(it,in,iy,ip) == 0) then
                     idx_arg(it,in,iy,ip) = -7
                  end if
               else
                  if (iarg(iy) == 1) then
                     if (idx_arg(it,in,iy+1,ip) == 0) then
                        idx_arg(it,in,iy+1,ip) = 8
                     end if
                     if (idx_arg(it,in,iy,ip) == 0) then
                        idx_arg(it,in,iy,ip) = -8
                     end if
                  end if
               end if
            end if
         end if
      end do
   end do
end do
if (1 == 0) then
it = 0
do in=1,dim_n,1
   iarg(0) = 0
   do iy=0,(dim_y-1),1
      if ((idx_arg(it,in,iy,0) == 1).and.&
           (idx_arg(it,in,iy+1,0) == 1)) then
         iarg(iy) = 1
      else
         iarg(iy) = 0
      end if
   end do
   iarg(dim_y) = 0
   do iy=0,(dim_y-1),1
      if ((iarg(iy) == 1).and.(iarg(iy+1) == 1)) then
         iarg(iy) = 2
      end if
   end do
   do iy=0,(dim_y-2),1
      if ((iarg(iy) == 2).and.(iarg(iy+2) >= 1)) then
         iarg(iy) = 3
      end if
   end do
   do iy=0,(dim_y-3),1
      if ((iarg(iy) == 3).and.(iarg(iy+3) >= 1)) then
         iarg(iy) = 4
      end if
   end do
   do iy=0,dim_y,1
      if (iarg(iy) == 4) then
         idx_arg(it,in,iy+4,ip) = 3
         idx_arg(it,in,iy+3,ip) = 2
         idx_arg(it,in,iy+2,ip) = 1
         if (idx_arg(it,in,iy+1,ip) == 0) then
            idx_arg(it,in,iy+1,ip) = -2
         end if
         if (idx_arg(it,in,iy,ip) == 0) then
            idx_arg(it,in,iy,ip) = -3
         end if
      else
         if (iarg(iy) == 3) then
            if (idx_arg(it,in,iy+3,ip) == 0) then
               idx_arg(it,in,iy+3,ip) = 5
            end if
            if (idx_arg(it,in,iy+2,ip) == 0) then
               idx_arg(it,in,iy+2,ip) = 4
            end if
            if (idx_arg(it,in,iy+1,ip) == 0) then
               idx_arg(it,in,iy+1,ip) = -4
            end if
            if (idx_arg(it,in,iy,ip) == 0) then
               idx_arg(it,in,iy,ip) = -5
            end if
         else
            if (iarg(iy) == 2) then
               if (idx_arg(it,in,iy+2,ip) == 0) then
                  idx_arg(it,in,iy+2,ip) = 7
               end if
               if (idx_arg(it,in,iy+1,ip) == 0) then
                  idx_arg(it,in,iy+1,ip) = 6
               end if
               if (idx_arg(it,in,iy,ip) == 0) then
                  idx_arg(it,in,iy,ip) = -7
               end if
            else
               if (iarg(iy) == 1) then
                  if (idx_arg(it,in,iy+1,ip) == 0) then
                     idx_arg(it,in,iy+1,ip) = 8
                  end if
                  if (idx_arg(it,in,iy,ip) == 0) then
                     idx_arg(it,in,iy,ip) = -8
                  end if
               end if
            end if
         end if
      end if
   end do
end do
end if

! define differentiation rules
do ip=1,3,1
   choose_parameter: select case (ip)
   case (1) ! T
      ia_min = 0
      ia_max = dim_t
   case (2) ! n_b
      ia_min = 1
      ia_max = dim_n
   case (3) ! Y_q
      ia_min = 0
      ia_max = dim_y
   end select choose_parameter
! initialisation with zero
   do ia=ia_min,ia_max,1
      do ii=-4,4,1
         do ir=-8,8,1
            r1d(ip,ia,ii,ir) = 0.d00
            r2d(ip,ia,ii,ir) = 0.d00
         end do
      end do

! five-point formulas
      do iq=-2,2,1
         choose_rule5: select case (iq)
         case (2)
            ir = 3
         case (1)
            ir = 2
         case (0)
            ir = 1
         case (-1)
            ir = -2
         case (-2)
            ir = -3
         end select choose_rule5
         iap = ia-iq
         if ((iap >= (ia_min+2)).and.(iap <= (ia_max-2))) then
            zp2 = tab_para(iap+2,ip)-tab_para(ia,ip)
            zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
            z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
            zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
            zm2 = tab_para(iap-2,ip)-tab_para(ia,ip)
! first derivative
            r1d(ip,ia, 2-iq,ir) = &
                 -((zp+z0)*zm*zm2+zp*z0*(zm+zm2))/&
                 ((zp2-zp)*(zp2-z0)*(zp2-zm)*(zp2-zm2))
            r1d(ip,ia, 1-iq,ir) = &
                 -((zp2+z0)*zm*zm2+zp2*z0*(zm+zm2))/&
                 ((zp-zp2)*(zp-z0)*(zp-zm)*(zp-zm2))
            r1d(ip,ia,  -iq,ir) = &
                 -((zp2+zp)*zm*zm2+zp2*zp*(zm+zm2))/&
                 ((z0-zp2)*(z0-zp)*(z0-zm)*(z0-zm2))
            r1d(ip,ia,-1-iq,ir) = &
                 -((zp2+zp)*z0*zm2+zp2*zp*(z0+zm2))/&
                 ((zm-zp2)*(zm-zp)*(zm-z0)*(zm-zm2))
            r1d(ip,ia,-2-iq,ir) = &
                 -((zp2+zp)*z0*zm+zp2*zp*(z0+zm))/&
                 ((zm2-zp2)*(zm2-zp)*(zm2-z0)*(zm2-zm))
! second derivative
            r2d(ip,ia, 2-iq,ir) = 2.d00*&
                 (zp*z0+(zp+z0)*(zm+zm2)+zm*zm2)/&
                 ((zp2-zp)*(zp2-z0)*(zp2-zm)*(zp2-zm2))
            r2d(ip,ia, 1-iq,ir) = 2.d00*&
                 (zp2*z0+(zp2+z0)*(zm+zm2)+zm*zm2)/&
                 ((zp-zp2)*(zp-z0)*(zp-zm)*(zp-zm2))
            r2d(ip,ia,  -iq,ir) = 2.d00*&
                 (zp2*zp+(zp2+zp)*(zm+zm2)+zm*zm2)/&
                 ((z0-zp2)*(z0-zp)*(z0-zm)*(z0-zm2))
            r2d(ip,ia,-1-iq,ir) = 2.d00*&
                 (zp2*zp+(zp2+zp)*(z0+zm2)+z0*zm2)/&
                 ((zm-zp2)*(zm-zp)*(zm-z0)*(zm-zm2))
            r2d(ip,ia,-2-iq,ir) = 2.d00*&
                 (zp2*zp+(zp2+zp)*(z0+zm)+z0*zm)/&
                 ((zm2-zp2)*(zm2-zp)*(zm2-z0)*(zm2-zm))
         end if
      end do

! four-point formulas
      do iq=-2,2,1
         choose_rule4: select case (iq)
         case (2)
            ir = 5
         case (1)
            ir = 4
         case (-1)
            ir = -4
         case (-2)
            ir = -5
         end select choose_rule4
         iap = ia-iq
         if (iq > 0) then
            if ((iap >= (ia_min+1)).and.(iap <= (ia_max-2))) then
               zp2 = tab_para(iap+2,ip)-tab_para(ia,ip)
               zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
               z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
               zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
               r1d(ip,ia, 2-iq,ir) =&
                    (zp*z0+zp*zm+z0*zm)/&
                    ((zp2-zp)*(zp2-z0)*(zp2-zm))
               r1d(ip,ia, 1-iq,ir) =&
                    (zp2*z0+zp2*zm+z0*zm)/&
                    ((zp-zp2)*(zp-z0)*(zp-zm))
               r1d(ip,ia,  -iq,ir) =&
                    (zp2*zp+zp2*zm+zp*zm)/&
                    ((z0-zp2)*(z0-zp)*(z0-zm))
               r1d(ip,ia,-1-iq,ir) =&
                    (zp2*zp+zp2*z0+zp*z0)/&
                    ((zm-zp2)*(zm-zp)*(zm-z0))
               r2d(ip,ia, 2-iq,ir) = 2.d00*&
                    (-zp-z0-zm)/&
                    ((zp2-zp)*(zp2-z0)*(zp2-zm))
               r2d(ip,ia, 1-iq,ir) = 2.d00*&
                    (-zp2-z0-zm)/&
                    ((zp-zp2)*(zp-z0)*(zp-zm))
               r2d(ip,ia,  -iq,ir) = 2.d00*&
                    (-zp2-zp-zm)/&
                    ((z0-zp2)*(z0-zp)*(z0-zm))
               r2d(ip,ia,-1-iq,ir) = 2.d00*&
                    (-zp2-zp-z0)/&
                    ((zm-zp2)*(zm-zp)*(zm-z0))
            end if
         else if (iq < 0) then
            if ((iap >= (ia_min+2)).and.(iap <= (ia_max-1))) then
               zp  = tab_para(iap+1,ip)-tab_para(ia,ip)
               z0  = tab_para(iap  ,ip)-tab_para(ia,ip)
               zm  = tab_para(iap-1,ip)-tab_para(ia,ip)
               zm2 = tab_para(iap-2,ip)-tab_para(ia,ip)
               r1d(ip,ia, 1-iq,ir) =&
                    (z0*zm+z0*zm2+zm*zm2)/&
                    ((zp-z0)*(zp-zm)*(zp-zm2))
               r1d(ip,ia,  -iq,ir) =&
                    (zp*zm+zp*zm2+zm*zm2)/&
                    ((z0-zp)*(z0-zm)*(z0-zm2))
               r1d(ip,ia,-1-iq,ir) =&
                    (zp*z0+zp*zm2+z0*zm2)/&
                    ((zm-zp)*(zm-z0)*(zm-zm2))
               r1d(ip,ia,-2-iq,ir) =&
                    (zp*z0+zp*zm+z0*zm)/&
                    ((zm2-zp)*(zm2-z0)*(zm2-zm))
               r2d(ip,ia, 1-iq,ir) = 2.d00*&
                    (-z0-zm-zm2)/&
                    ((zp-z0)*(zp-zm)*(zp-zm2))
               r2d(ip,ia,  -iq,ir) = 2.d00*&
                    (-zp-zm-zm2)/&
                    ((z0-zp)*(z0-zm)*(z0-zm2))
               r2d(ip,ia,-1-iq,ir) = 2.d00*&
                    (-zp-z0-zm2)/&
                    ((zm-zp)*(zm-z0)*(zm-zm2))
               r2d(ip,ia,-2-iq,ir) = 2.d00*&
                    (-zp-z0-zm)/&
                    ((zm2-zp)*(zm2-z0)*(zm2-zm))
            end if
         end if
      end do

! three-point formulas
      do iq=-1,1,1
         choose_rule3: select case (iq)
         case (1)
            ir = 7
         case (0)
            ir = 6
         case (-1)
            ir = -7
         end select choose_rule3
         iap = ia-iq
         if ((iap >= (ia_min+1)).and.(iap <= (ia_max-1))) then
            zp = tab_para(iap+1,ip)-tab_para(ia,ip)
            z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
            zm = tab_para(iap-1,ip)-tab_para(ia,ip)
            r1d(ip,ia, 1-iq,ir) =&
                 -(z0+zm)/&
                 ((zp-z0)*(zp-zm))
            r1d(ip,ia,  -iq,ir) =&
                 -(zp+zm)/&
                 ((z0-zp)*(z0-zm))
            r1d(ip,ia,-1-iq,ir) =&
                 -(zp+z0)/&
                 ((zm-zp)*(zm-z0))
            r2d(ip,ia, 1-iq,ir) = 2.d00/&
                 ((zp-z0)*(zp-zm))
            r2d(ip,ia,  -iq,ir) = 2.d00/&
                 ((z0-zp)*(z0-zm))
            r2d(ip,ia,-1-iq,ir) = 2.d00/&
                 ((zm-zp)*(zm-z0))
         end if
      end do

! two-point formulas
      do iq=-1,1,1
         choose_rule2: select case (iq)
         case (1)
            ir = 8
         case (-1)
            ir = -8
         end select choose_rule2
         iap = ia-iq
         if (iq > 0) then
            if ((iap >= (ia_min)).and.(iap <= (ia_max-1))) then
               zp = tab_para(iap+1,ip)-tab_para(ia,ip)
               z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
               r1d(ip,ia, 1-iq,ir) =&
                    1.d00/(zp-z0)
               r1d(ip,ia,  -iq,ir) =&
                    1.00/(z0-zp)
            end if
         else if (iq < 0) then
            if ((iap >= (ia_min+1)).and.(iap <= (ia_max))) then
               z0 = tab_para(iap  ,ip)-tab_para(ia,ip)
               zm = tab_para(iap-1,ip)-tab_para(ia,ip)
               r1d(ip,ia,  -iq,ir) =&
                    1.d00/(z0-zm)
               r1d(ip,ia,-1-iq,ir) =&
                    1.d00/(zm-z0)
            end if
         end if

      end do
   end do
end do

return
end SUBROUTINE get_diff_rules
!***********************************************************************
SUBROUTINE get_eos_report(iwr)
! Stefan Typel for the CompOSE core team, version 0.19, 2012/12/19
USE compose_internal
implicit none
integer iwr,icc
integer, parameter :: dim_q=20,dim_d=16
integer :: ierror,ierr,ip,iq,it,in,iy,iflag,icnt(0:3),ifa,itmp,&
     nqty,nphase,iphase(dim_q),nphase_miss,&
     nqtyp,iqtyp(dim_q),np_miss,&
     nqtyq,iqtyq(dim_q),nq_miss,&
     nqtym,iqtym(dim_q),nm_miss,&
     distr(2,-1:dim_d)
double precision :: nsat,bsat,k,kp,j,l,ksym,tmp,dfa_max(2),&
     dlnt,f1,f2,fn,en,nb,yq,t

! error counter
ierr = 0

! dimension of eos table and interpolation in parameters
eos_dim = 0
do ip=1,3,1
   idx_ipl(ip) = 1
   if (idx_max(ip) > idx_min(ip)) eos_dim = eos_dim+1
end do

! type of eos table
if (eos_dim == 3) eos_type = 1 ! general purpose eos
if (eos_dim == 2) then
   if (idx_min(1) == idx_max(1)) then
      if (dabs(tab_para(idx_min(1),1)) < 1.d-08) then
         eos_type = 2 ! zero temperature
         tab_para(idx_min(1),1) = 0.d00
         idx_ipl(1) = 0
      else
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 30
      end if
   end if
   if (idx_min(3) == idx_max(3)) then
      if (dabs(tab_para(idx_min(3),3)-0.5d00) < 1.d-08) then
         eos_type = 3 ! symmetric nuclear matter
         tab_para(idx_min(3),3) = 0.5d00
         idx_ipl(3) = 0
      else
         if (dabs(tab_para(idx_min(3),3)) < 1.d-08) then
            eos_type = 4 ! neutron matter
            tab_para(idx_min(3),3) = 0.d00
            idx_ipl(3) = 0
         else
            if ((dabs(tab_para(idx_min(3),3)-0.5d00)-0.5d00) > 1.d-08) then
               eos_type = 5 ! beta-equilibrated matter
               tab_para(idx_min(3),3) = 0.d00
               idx_min(3) = 0
               idx_max(3) = 0
               idx_ipl(3) = 0
            else
               ierr = ierr+1
               if (ierr < dim_err) error_msg(ierr) = 31
            end if
         end if
      end if
   end if
end if
if (eos_dim == 1) then
   if (idx_min(1) == idx_max(1)) then
      if (dabs(tab_para(idx_min(1),1)) > 1.d-08) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 30
      else
         if (idx_min(3) == idx_max(3)) then
            if (dabs(tab_para(idx_min(3),3)-0.5d00) < 1.d-08) then
               eos_type = 6 ! zero temperature, symmetric nuclear matter
               tab_para(idx_min(3),3) = 0.5d00
               idx_ipl(1) = 0
               idx_ipl(3) = 0
            else
               if (dabs(tab_para(idx_min(3),3)) < 1.d-08) then
                  eos_type = 7 ! zero temperature, neutron matter
                  tab_para(idx_min(3),3) = 0.d00
                  idx_ipl(1) = 0
                  idx_ipl(3) = 0
               else
                  if ((dabs(tab_para(idx_min(3),3)-0.5d00)-0.5d00) > 1.d-08) then
                     eos_type = 8 ! zero temperature, beta-equilibrated matter
                     tab_para(idx_min(3),3) = 0.d00
                     idx_min(3) = 0
                     idx_max(3) = 0
                     idx_ipl(1) = 0
                     idx_ipl(3) = 0
                  else
                     ierr = ierr+1
                     if (ierr < dim_err) error_msg(ierr) = 31
                  end if
               end if
            end if
         end if
      end if
   end if
end if

! completeness and consistency check
do ip=0,3,1
   icnt(ip) = 0
end do
do it=idx_min(1),idx_max(1),1
   do in=idx_min(2),idx_max(2),1
      do iy=idx_min(3),idx_max(3),1
         if ((idx_argx(it,in,iy,1) == idx_ex(1)).and.&
              (idx_argx(it,in,iy,2) == idx_ex(2)).and.&
              (idx_argx(it,in,iy,3) == idx_ex(3))) then
            idx_arg(it,in,iy,0) = 1
            icnt(0) = icnt(0)+1
         else
            do ip=1,3,1
               if ((idx_ex(ip) == 1).and.(idx_argx(it,in,iy,ip) /= 1)) then
                  icnt(ip) = icnt(ip)+1
               end if
            end do
         end if
      end do
   end do
end do



if (iwr == 1) then
   write(*,*)
   if (idx_ex(1) == 1) then
      write(*,*) icnt(1),' entries in eos.thermo missing'
   else
      write(*,*) 'no table eos.thermo'
      stop
   end if
   if (idx_ex(2) == 1) then
      write(*,*) icnt(2),' entries in eos.compo missing'
   else
      write(*,*) 'no table eos.compo'
   end if
   if (idx_ex(3) == 1) then
      write(*,*) icnt(3),' entries in eos.micro missing'
   else
      write(*,*) 'no table eos.micro'
   end if
   write(*,*)
   write(*,*) icnt(0),' complete entries in eos tables'
end if
   
! additional quantities in eos.thermo
nqty = nadd_max

! quantities in eos.compo
! phases
do ip=1,dim_q,1
   iphase(ip)= -100
end do
nphase = 0
nphase_miss = 0
if (idx_ex(2) == 1) then
   do it=idx_min(1),idx_max(1),1
      do in=idx_min(2),idx_max(2),1
         do iy=idx_min(3),idx_max(3),1
            iflag = 0
            do ip=1,nphase,1
               if (idx_arg(it,in,iy,4) == iphase(ip)) iflag = 1
            end do
            if (iflag == 0) then
               if (nphase < dim_q) then
                  nphase = nphase+1
                  iphase(nphase) = idx_arg(it,in,iy,4)
               else
                  nphase_miss = nphase_miss+1
               end if
            end if
         end do
      end do
   end do
end if

! particles
do ip=1,dim_q,1
   iqtyp(ip) = -1
end do
nqtyp = 0
np_miss = 0
if (idx_ex(2) == 1) then
   do it=idx_min(1),idx_max(1),1
      do in=idx_min(2),idx_max(2),1
         do iy=idx_min(3),idx_max(3),1
            do ip=1,np_max,1
               iflag = 0
               do iq=1,nqtyp,1
                  if (idxp_compo(it,in,iy,ip) == iqtyp(iq)) iflag = 1
               end do
               if (iflag == 0) then
                  if (idxp_compo(it,in,iy,ip) > -1) then
                     if (nqtyp < dim_q) then
                        nqtyp = nqtyp+1
                        iqtyp(nqtyp) = idxp_compo(it,in,iy,ip)
                     else
                        np_miss = np_miss+1
                     end if
                  end if
               end if
            end do
         end do
      end do
   end do
end if

! particle sets
do ip=1,dim_q,1
   iqtyq(ip) = -1
end do
nqtyq = 0
nq_miss = 0
if (idx_ex(2) == 1) then
   do it=idx_min(1),idx_max(1),1
      do in=idx_min(2),idx_max(2),1
         do iy=idx_min(3),idx_max(3),1
            do ip=1,nq_max,1
               iflag = 0
               do iq=1,nqtyq,1
                  if (idxq_compo(it,in,iy,ip) == iqtyq(iq)) iflag = 1
               end do
               if (iflag == 0) then
                  if (idxq_compo(it,in,iy,ip) > -1) then
                     if (nqtyq < dim_q) then
                        nqtyq = nqtyq+1
                        iqtyq(nqtyq) = idxq_compo(it,in,iy,ip)
                     else
                        nq_miss = nq_miss+1
                     end if
                  end if
               end if
            end do
         end do
      end do
   end do
end if

! quantities in eos.micro
do ip=1,dim_q,1
   iqtym(ip) = -1
end do
nqtym = 0
nm_miss = 0
if (idx_ex(3) == 1) then
   do it=idx_min(1),idx_max(1),1
      do in=idx_min(2),idx_max(2),1
         do iy=idx_min(3),idx_max(3),1
            do ip=1,nm_max,1
               iflag = 0
               do iq=1,nqtym,1
                  if (idx_mic(it,in,iy,ip) == iqtym(iq)) iflag = 1
               end do
               if (iflag == 0) then
                  if (idx_mic(it,in,iy,ip) > -1) then
                     if (nqtym < dim_q) then
                        nqtym = nqtym+1
                        iqtym(nqtym) = idx_mic(it,in,iy,ip)
                     else
                        nm_miss = nm_miss+1
                     end if
                  end if
               end if
            end do
         end do
      end do
   end do
end if

nsat = 0.d00
bsat = 0.d00
k    = 0.d00
kp   = 0.d00
j    = 0.d00
l    = 0.d00
ksym = 0.d00
open(unit=23,file='eos.beta',&
     status='unknown',action='write',iostat=ierror)
if (incl_l /= 1) then
   if ((eos_type == 1).or.(eos_type == 2)) then
      ! nuclear matter parameters
      call get_eos_nmp(nsat,bsat,k,kp,j,l,ksym,ierr)
   end if
else
! EoS of beta-equilibrated matter at lowest temperature
   if ((eos_type == 1).or.(eos_type == 2)) then
      t = tab_para(idx_min(1),1)
      yq = 0.5d00
      do in=idx_min(2),idx_max(2),1
         nb = tab_para(in,2)
         !++++++++++++++++++++++++++++
         call get_eos(t,nb,yq,3,3,3,1)
         !++++++++++++++++++++++++++++
         if (yq > 0.d00) then
            write(23,*) nb,yq,(eos_thermo(6)+1.d00)*m_n*nb,eos_thermo(1)
         end if
      end do
   end if
end if
close (unit=23)

! thermodynamic consistency check and distribution of errors
! delta(f/n)/(f/n) and delta(e/n)/(e/n)
icc = 1
if (icc == 1) then
   dlnt = dlog(10.d00)
   ifa = 0
   do ip=1,2,1
      dfa_max(ip) = 0.d00
      do iq=0,dim_d,1
         distr(ip,iq) = 0
      end do
   end do
   if (incl_l == 1) then
      f1 = 1.d00
      f2 = 0.d00
   else
      f1 = 0.d00
      f2 = 1.d00
   end if
   if (idx_ex(1) == 1) then
      do it=idx_min(1),idx_max(1),1
         do in=idx_min(2),idx_max(2),1
            do iy=idx_min(3),idx_max(3),1
               if (idx_arg(it,in,iy,0) == 1) then
                  ifa = ifa+1
                  fn = m_n*(tab_thermo(it,in,iy,3)&
                       +tab_para(iy,3)*(f1*tab_thermo(it,in,iy,5)&
                       +f2*tab_thermo(it,in,iy,4)))&
                       -tab_thermo(it,in,iy,1)
                  en = fn+tab_para(it,1)*tab_thermo(it,in,iy,2)
                  fn = m_n*tab_thermo(it,in,iy,6)-fn
                  fn = fn/(m_n*(tab_thermo(it,in,iy,6)+1.d00))
                  en = m_n*tab_thermo(it,in,iy,7)-en
                  en = en/(m_n*(tab_thermo(it,in,iy,7)+1.d00))
                  do ip=1,2,1
                     if (ip == 1) then
                        tmp = fn
                     else
                        tmp = en
                     end if
                     if (dabs(tmp) > dabs(dfa_max(ip))) then
                        dfa_max(ip) = tmp
                     end if
                     tmp = dabs(tmp)
                     if (tmp < 1.d00) then
                        if (tmp < 1.d-18) then
                           itmp = dim_d
                        else
                           itmp = int(-dlog(tmp)/dlnt)
                           if (itmp > dim_d) itmp = dim_d
                        end if
                     else
                        itmp = -1
                     end if
                     distr(ip,itmp) = distr(ip,itmp)+1
                  end do
               end if
            end do
         end do
      end do
   end if
end if



open(unit=21,file='eos.report',&
     status='unknown',action='write',iostat=ierror)

write(21,*) '# dimension and type of eos table'
write(21,*) eos_dim,eos_type
do ip=1,3,1
   if (ip == 1)& 
        write(21,*) '# number of grid points in temperature T,',&
        ' minimum temperature, maximum temperature [MeV]'
   if (ip == 2)&
        write(21,*) '# number of grid points in baryon number density n_b,',&
        ' minimum baryon number density, maximum baryon number density [fm^-3]'
   if (ip == 3)&
        write(21,*) '# number of grid points in charge fraction Y_q,',&
        ' minimum charge fraction, maximum charge fraction [dimensionless]'
   write(21,'(3x,i4,3x,e16.8,3x,e16.8)') (idx_max(ip)-idx_min(ip)+1),&
        tab_para(idx_min(ip),ip),tab_para(idx_max(ip),ip)
end do
do ip=1,3,1
   if (ip == 1)& 
        write(21,*) '# thermodynamical data exist (1) or not (0)'
   if (ip == 2)&
        write(21,*) '# composition data exist (1) or not (0)'
   if (ip == 3)&
        write(21,*) '# microscopic data exist (1) or not (0)'
   write(21,*) idx_ex(ip)
end do
write(21,*) '# total number of table entries'
write(21,*) icnt(0)
write(21,*) '# number of additional quantities in eos.thermo'
write(21,*) nqty
write(21,*) '# number of phases, number of not listed phases'
write(21,*) nphase,nphase_miss
write(21,*) '# indices of phases (see data sheet)'
if (nphase > 0) then
   write(21,*) (iphase(ip),ip=1,nphase,1)
else
   write(21,*)
end if
write(21,*) '# number of particles, number of not listed particles'
write(21,*) nqtyp,np_miss
write(21,*) '# indices of particles (see manual)'
if (nqtyp > 0) then
   write(21,*) (iqtyp(ip),ip=1,nqtyp,1)
else
   write(21,*)
end if
write(21,*) '# number of particle sets, number of not listed particle sets'
write(21,*) nqtyq,nq_miss
write(21,*) '# indices of particle sets (see manual and data sheet)'
if (nqtyq > 0) then
   write(21,*) (iqtyq(ip),ip=1,nqtyq,1)
else
   write(21,*)
end if
write(21,*) '# number of quantities in eos.micro, number of not listed quantities'
write(21,*) nqtym,nm_miss
write(21,*) '# indices of quantities (see manual)'
if (nqtym > 0) then
   write(21,*) (iqtym(ip),ip=1,nqtym,1)
else
   write(21,*)
end if
write(21,*) '# saturation density of symmetric nuclear matter [fm^-3]'
write(21,'(f10.5)') nsat
write(21,*) '# binding energy B_sat per baryon at saturation [MeV]'
write(21,'(f10.3)') bsat
write(21,*) '# incompressibility K [MeV]'
write(21,'(f10.3)') k
write(21,*) '# skewness K^prime [MeV]'
write(21,'(f10.3)') kp
write(21,*) '# symmetry energy J [MeV]'
write(21,'(f10.3)') j
write(21,*) '# symmetry energy slope parameter L [MeV]'
write(21,'(f10.3)') l
write(21,*) '# symmetry incompressibility K_sym [MeV]'
write(21,'(f10.3)') ksym

if (icc == 1) then
   tmp = dble(ifa)
   open(unit=22,file='eos.errdistr',&
        status='unknown',action='write',iostat=ierror)
   if (ifa > 0) then
      ! for xmgrace file
      write(22,*) '@target G0.S0'
      write(22,*) '@type xy'
      do iq=0,dim_d,1
         fn = dble(distr(1,iq))/tmp
         if (fn < 1.d-30) fn = 1.d-30
         write(22,*) 10.d00**(-iq),fn
      end do
      write(22,*) 10.d00**(-dim_d-1),fn
      write(22,*) '&'
      write(22,*) '@target G0.S1'
      write(22,*) '@type xy'
      do iq=0,dim_d,1
         en = dble(distr(2,iq))/tmp
         if (en < 1.d-30) en = 1.d-30
         write(22,*) 10.d00**(-iq),en
      end do
      write(22,*) 10.d00**(-dim_d-1),en
      write(22,*) '&'
   else
      write(*,*) 'error in error distribution'
      stop
   end if
   close(22)
   write(21,*) '# maximum relative error in free energy per baryon'
   write(21,'(e16.8)') dfa_max(1)
   write(21,*) '# maximum relative error in internal energy per baryon'
   write(21,'(e16.8)') dfa_max(2)
end if

close(unit=21)


call write_errors(ierr)

return
end SUBROUTINE get_eos_report
!***********************************************************************
SUBROUTINE get_eos_nmp(nsat,bsat,k,kp,j,l,ksym,ierr)
! Stefan Typel for the CompOSE core team, version 0.05, 2012/12/17
USE compose_internal
implicit none
integer :: i,idx,ierr,ipl,irule
double precision :: t,nb,y,f,fp,nsat,bsat,k,kp,j,l,ksym,&
     nb_min,nb_max,f_min,f_max,d_nb,d_f,y_min,y_max,&
     g1(-2:2),g2(-2:2),g3(-2:2),d_b,d_b2

ipl = 3

irule = 5

!if ((incl_l /= 1).and.((eos_type == 1).or.(eos_type == 2))) then
! determination of saturation density
! find intervall around minimum in free energy
   t = tab_para(idx_min(1),1)
   y = 0.5d00
   idx = -1
   fp = 1.d30
   do i=idx_max(2),idx_min(2),-1
      nb = tab_para(i,2)
      call get_eos_sub(t,nb,y,ipl,ipl,ipl,ierr,1)
      f = eos_thermo(6)*m_n
      if (f < fp) then
         fp = f
         idx = i
      end if
   end do
   if ((idx_min(2)<idx).and.(idx_max(2)>idx)) then
      nb_min = tab_para((idx-1),2)
      call get_eos_sub(t,nb_min,y,ipl,ipl,ipl,ierr,1)
      f_min = eos_thermo(1)
      nb_max = tab_para((idx+1),2)
      call get_eos_sub(t,nb_max,y,ipl,ipl,ipl,ierr,1)
      f_max = eos_thermo(1)
      ! determination of zero in pressure
      nb = 0.d00
      d_nb = 1.d00
      do while (dabs(d_nb) > 1.d-12)
         d_nb = nb
         d_f = f_max-f_min
!         if (dabs(d_f) > 1.d-08) then
!            nb = (f_max*nb_min-f_min*nb_max)/d_f
!         else
            nb = 0.5d00*(nb_max+nb_min)
!         end if
         call get_eos_sub(t,nb,y,ipl,ipl,ipl,ierr,1)
         f = eos_thermo(1)
         d_nb = d_nb-nb
         if (f_min*f > 0.d00) then
            nb_min = nb
            f_min = f
         else
            nb_max = nb
            f_max = f
         end if
      end do
      ! parameters for symmetric nuclear matter
      nsat = nb
      bsat = -(eos_thermo(6)*m_n+0.5d00*(m_n-m_p))
      k = 9.d00*nsat*eos_thermo_add(1)
      kp = 27.d00*nsat*nsat*eos_thermo_add(2)-6.d00*k
      ! isospin related parameters
      y_min = tab_para(idx_min(3),3)
      y_max = tab_para(idx_max(3),3)
      if ((y_min < 0.48d00).and.(y_max > 0.52d00)) then
         g1(0) = eos_thermo(6)*m_n
         g2(0) = eos_thermo(1)
         g3(0) = eos_thermo_add(1)
         d_b = 0.01d00
         d_b2 = d_b*d_b
         if (irule == 3) then
            y = 0.5d00-d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(-1) = eos_thermo(6)*m_n
            g2(-1) = eos_thermo(1)
            g3(-1) = eos_thermo_add(1)
            y = 0.5d00+d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(1) = eos_thermo(6)*m_n
            g2(1) = eos_thermo(1)
            g3(1) = eos_thermo_add(1)
            j = 0.125d00*(g1(-1)-2.d00*g1(0)+g1(1))/d_b2
            l = 0.375d00*(g2(-1)-2.d00*g2(0)+g2(1))/(d_b2*nsat)
            ksym = 1.125d00*nsat*(g3(-1)-2.d00*g3(0)+g3(1))/d_b2-3.d00*l
         else
            y = 0.5d00-2.d00*d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(-2) = eos_thermo(6)*m_n
            g2(-2) = eos_thermo(1)
            g3(-2) = eos_thermo_add(1)
            y = 0.5d00-d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(-1) = eos_thermo(6)*m_n
            g2(-1) = eos_thermo(1)
            g3(-1) = eos_thermo_add(1)
            y = 0.5d00+d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(1) = eos_thermo(6)*m_n
            g2(1) = eos_thermo(1)
            g3(1) = eos_thermo_add(1)
            y = 0.5d00+2.d00*d_b
            call get_eos_sub(t,nsat,y,ipl,ipl,ipl,ierr,1)
            g1(2) = eos_thermo(6)*m_n
            g2(2) = eos_thermo(1)
            g3(2) = eos_thermo_add(1)
            j = (-g1(-2)+16.d00*g1(-1)-30.d00*g1(0)+16.d00*g1(1)-g1(2))&
                 /(96.d00*d_b2)
            l = (-g2(-2)+16.d00*g2(-1)-30.d00*g2(0)+16.d00*g2(1)-g2(2))&
                 /(32.d00*d_b2*nsat)
            ksym = 3.d00*nsat*&
                 (-g3(-2)+16.d00*g3(-1)-30.d00*g3(0)+16.d00*g3(1)-g3(2))&
                 /(32.d00*d_b2)-3.d00*l
         endif
      end if
   end if

!end if

return
end SUBROUTINE get_eos_nmp
!***********************************************************************
SUBROUTINE read_eos_tables_tny(iwr)
! Stefan Typel for the CompOSE core team, version 0.10, 2012/11/05
USE compose_internal
implicit none
integer :: iwr,it,in,iy,ip,ipp,ia,iunit,ierror,itmp,itmpp,ierr
integer :: icnt(3)
double precision :: tmp

! counter for total number of read table entries
icnt(1) = 0
icnt(2) = 0
icnt(3) = 0

! error counter
ierr = 0

! initialization
do ia=0,dim_a,1
   do ip=1,6,1
      tab_para(ia,ip) = -1.d00
   end do
end do
do it=0,dim_t,1
   do in=1,dim_n,1
      do iy=0,dim_y,1
         do ip=0,3,1
            idx_arg(it,in,iy,ip) = 0
         end do
         idx_arg(it,in,iy,4) = -1
      end do
   end do
end do

! temperature
ip = 1
iunit = 20+ip
open(unit=iunit,file='eos.t',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   if (iwr == 1) then
      write(*,*)
      write(*,*) 'reading parameter table for temperature'
   end if
   read(iunit,*) itmp
   if (itmp < 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 2
   end if
   read(iunit,*) itmpp
   if (itmpp > dim_t) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 3
   end if
   if (itmpp < itmp) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 4
   end if
else
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 1
end if
if (ierr == 0) then
   do it=itmp,itmpp,1
      read(iunit,*) tmp
      tab_para(it,ip) = tmp
      icnt(ip) = icnt(ip)+1
   end do
   do it=(itmp+1),itmpp,1
      if ((ierr == 0).and.(tab_para(it,ip) < tab_para(it-1,ip))) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 5
      end if
   end do
   if (ierr == 0) then
      idx_min(ip) = itmp
      idx_max(ip) = itmpp
   end if
end if
close(unit=iunit)
if (iwr == 1) then
   write(*,*) icnt(ip),&
        ' entries of parameter table for temperatures read'
end if
ipp = ip+3
do it=0,(dim_t-1),1
   tab_para(it,ipp) = tab_para(it+1,ip)-tab_para(it,ip)
end do

! baryon number density
ip = 2
iunit = 20+ip
open(unit=iunit,file='eos.nb',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   if (iwr == 1) then
      write(*,*)
      write(*,*) 'reading parameter table for baryon number density'
   end if
   read(iunit,*) itmp
   if (itmp <= 0) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 7
   end if
   read(iunit,*) itmpp
   if (itmpp > dim_n) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 8
   end if
   if (itmpp < itmp) then
      ierr = ierr+1
      if (ierr < dim_err) error_msg(ierr) = 9
   end if
else
   ierr = ierr+1
   error_msg(ierr) = 6
end if
if (ierr == 0) then
   do in=itmp,itmpp,1
      read(iunit,*) tmp
      tab_para(in,ip) = tmp
      icnt(ip) = icnt(ip)+1
   end do
   do it=(itmp+1),itmpp,1
      if ((ierr == 0).and.(tab_para(it,ip) < tab_para(it-1,ip))) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 10
      end if
   end do
   if (ierr == 0) then
      idx_min(ip) = itmp
      idx_max(ip) = itmpp
   end if
end if
close(unit=iunit)
if (iwr == 1) then
   write(*,*) icnt(ip),&
        ' entries of parameter table for baryon number densities read'
end if
ipp = ip+3
do in=1,(dim_n-1),1
   tab_para(in,ipp) = tab_para(in+1,ip)-tab_para(in,ip)
end do

! hadronic charge fraction
ip = 3
iunit = 20+ip
open(unit=iunit,file='eos.yq',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   if (iwr == 1) then
      write(*,*)
      write(*,*) 'reading parameter table for charge fraction'
   end if
   read(iunit,*) itmp
   read(iunit,*) itmpp
   if ((itmp < 0).or.(itmpp < 0)) then
      if ((itmp < -1).or.(itmpp < -1)) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 16
      end if
   else
      if (itmp < 0) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 12
      end if
      if (itmpp > dim_y) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 13
      end if
      if (itmpp < itmp) then
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 14
      end if
   end if
else
   ierr = ierr+1
   if (ierr < dim_err) error_msg(ierr) = 11
end if
if (ierr == 0) then
! beta equilibrium
   if (itmp == -1) then
      tab_para(0,ip) = -1.d00
      icnt(ip) = icnt(ip)+1
      idx_min(ip) = -1
      idx_max(ip) = -1
   else 
      do iy=itmp,itmpp,1
         read(iunit,*) tmp
         tab_para(iy,ip) = tmp
         icnt(ip) = icnt(ip)+1
      end do
      do it=(itmp+1),itmpp,1
         if ((ierr == 0).and.(tab_para(it,ip) < tab_para(it-1,ip))) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 15
         end if
      end do
      if (ierr == 0) then
         idx_min(ip) = itmp
         idx_max(ip) = itmpp
      end if
   end if
end if
close(unit=iunit)
if (iwr == 1) then
   write(*,*) icnt(ip),&
        ' entries of parameter table for charge fraction read'
end if
ipp = ip+3
do iy=0,(dim_y-1),1
   tab_para(iy,ipp) = tab_para(iy+1,ip)-tab_para(iy,ip)
end do

call write_errors(ierr)

return
end SUBROUTINE read_eos_tables_tny
!***********************************************************************
SUBROUTINE read_eos_table_thermo(iwr,nbl)
! Stefan Typel for the CompOSE core team, version 0.18, 2012/12/18
USE compose_internal
implicit none
integer :: iwr,nbl
integer :: ibl,ierror,icnt,it,in,iy,iv,iunit,iflag,itmp,nqty,iw,ierr
!integer :: iqty(dim_qty)
double precision :: tmp(dim_reg),qty(dim_qty)

! indices and stored quantities in eos_thermo
! 1: p/n_b
! 2: s/n_b
! 3: mu_b/m_n-1
! 4: mu_q/m_n
! 5: mu_e/m_n
! 6: f/(n_b*m_n)-1
! 7: e/(n_b_m_n)-1

! counter for total number of read table entries
icnt = 0

! error counter
ierr = 0

! initialization
do it=0,dim_t,1
   do in=1,dim_n,1
      do iy=0,dim_y,1
         do iv=0,4,1
            idx_arg(it,in,iy,iv) = 0
         end do
         do iv=1,dim_qty2,1
            tab_thermo(it,in,iy,iv) = 0.d00
         end do
      end do
   end do
end do

idx_ex(1) = 0
open(unit=20,file='eos.thermo',&
     status='old',action='read',iostat=ierror)
if (ierror == 0) then
   itmp = 0
else
   itmp = nbl-1
end if

do ibl=0,itmp,1
   ! flag for stop reading file
   iflag = 0
   iunit = 20+ibl
   
   choose_file: select case (ibl)
   case (0)
      if (itmp > 0) then
         open(unit=iunit,file='eos.thermo0',&
              status='old',action='read',iostat=ierror)
      end if
   case (1)
      open(unit=iunit,file='eos.thermo1',&
              status='old',action='read',iostat=ierror)
   case (2)
      open(unit=iunit,file='eos.thermo2',&
           status='old',action='read',iostat=ierror)
   case (3)
      open(unit=iunit,file='eos.thermo3',&
           status='old',action='read',iostat=ierror)
   case (4)
      open(unit=iunit,file='eos.thermo4',&
           status='old',action='read',iostat=ierror)
   case (5)
      open(unit=iunit,file='eos.thermo5',&
           status='old',action='read',iostat=ierror)
   case (6)
      open(unit=iunit,file='eos.thermo6',&
           status='old',action='read',iostat=ierror)
   case (7)
      open(unit=iunit,file='eos.thermo7',&
           status='old',action='read',iostat=ierror)
   case (8)
      open(unit=iunit,file='eos.thermo8',&
           status='old',action='read',iostat=ierror)
   case (9)
      open(unit=iunit,file='eos.thermo9',&
           status='old',action='read',iostat=ierror)
   end select choose_file
   
   if (ierror == 0) then
      idx_ex(1) = 1
      nadd_max = 0
      if (iwr == 1) then
         write(*,*)
         write(*,*) 'reading eos table(s) with thermodynamical properties'
      end if
      read(iunit,*,iostat=ierror) m_n,m_p,incl_l
      if ((ierror /= 0).or.(m_n < 0.d00).or.(m_p < 0.d00)) then 
         iflag = 1
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 27
      end if
      do while (iflag == 0)
         read(iunit,*,iostat=ierror) it,in,iy,(tmp(iv),iv=1,dim_reg,1),&
              nqty,(qty(iw),iw=1,nqty,1)         
         if (ierror /= 0) then 
            iflag = 1
         else 
            if ((nqty < 0).or.(nqty > dim_qty)) then
               iflag = 1
               ierr = ierr+1
               if (ierr < dim_err) error_msg(ierr) = 21
            end if
            idx_arg(it,in,iy,1) = 1
            do iv=1,dim_reg,1
               tab_thermo(it,in,iy,iv) = tmp(iv)
            end do
            if (nqty > nadd_max) nadd_max = nqty
            do iw=1,nqty,1
               if ((nqty > -1).or.(nqty < (dim_qty+1))) then
                  tab_thermo(it,in,iy,dim_reg+iw) = qty(iw)
               end if
            end do
            icnt = icnt+1
         end if
      end do
   end if
   
   close(unit=iunit)
end do
if (iwr == 1) then
   write(*,*) icnt,' entries of thermodynamical table read'
end if

call write_errors(ierr)
   
return
end SUBROUTINE read_eos_table_thermo
!***********************************************************************
SUBROUTINE read_eos_table_compo(iwr,nbl)
! Stefan Typel for the CompOSE core team, version 0.12, 2012/10/02
USE compose_internal
implicit none
integer :: iwr,nbl
integer :: ibl,ierror,icnt,it,in,iy,iphase,iunit,iflag,itmp,ierr,ic,&
     iv,iw,np,nq,alloc_status
integer :: iqtyp(dim_qtyp),iqtyq(dim_qtyq)
double precision :: qtyp(dim_qtyp),qtyq(dim_qtyq,3)

! counter for total number of read table entries
icnt = 0

! error counter
ierr = 0

! reading cycles
idx_ex(2) = 0
do ic=1,2,1
   if (ic == 1) then
      np_max = 0
      nq_max = 0
   end if
   open(unit=20,file='eos.compo',&
        status='old',action='read',iostat=ierror)
   if (ierror == 0) then
      itmp = 0
   else
      itmp = nbl-1
   end if
   
   do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = 20+ibl
      
      choose_file: select case (ibl)
      case (0)
         if (itmp > 0) then
            open(unit=iunit,file='eos.compo0',&
                 status='old',action='read',iostat=ierror)
         end if
      case (1)
         open(unit=iunit,file='eos.compo1',&
              status='old',action='read',iostat=ierror)
      case (2)
         open(unit=iunit,file='eos.compo2',&
              status='old',action='read',iostat=ierror)
      case (3)
         open(unit=iunit,file='eos.compo3',&
              status='old',action='read',iostat=ierror)
      case (4)
         open(unit=iunit,file='eos.compo4',&
              status='old',action='read',iostat=ierror)
      case (5)
         open(unit=iunit,file='eos.compo5',&
              status='old',action='read',iostat=ierror)
      case (6)
         open(unit=iunit,file='eos.compo6',&
              status='old',action='read',iostat=ierror)
      case (7)
         open(unit=iunit,file='eos.compo7',&
              status='old',action='read',iostat=ierror)
      case (8)
         open(unit=iunit,file='eos.compo8',&
              status='old',action='read',iostat=ierror)
      case (9)
         open(unit=iunit,file='eos.compo9',&
              status='old',action='read',iostat=ierror)
      end select choose_file
      
      if (ierror == 0) then
         idx_ex(2) = 1
         if (ic == 1) then
            if (iwr == 1) then
               write(*,*)
               write(*,*) 'reading eos table(s) with composition'
            end if
         end if
         do while (iflag == 0)
            read(iunit,*,iostat=ierror) it,in,iy,iphase,&
                 np,(iqtyp(iv),qtyp(iv),iv=1,np,1),&
                 nq,(iqtyq(iw),qtyq(iw,1),qtyq(iw,2),qtyq(iw,3),iw=1,nq,1)
            
            if (ic == 1) then
               if (ierror /= 0) then 
                  iflag = 1
               else 
                  if (np > np_max) np_max = np
                  if (nq > nq_max) nq_max = nq
                  icnt = icnt+1
               end if
            else
               if (ierror /= 0) then 
                  iflag = 1
               else 
                  idx_arg(it,in,iy,2) = 1
                  idx_arg(it,in,iy,4) = iphase
                  do iv=1,np,1
                     idxp_compo(it,in,iy,iv) = iqtyp(iv)
                     tabp_compo(it,in,iy,iv) = qtyp(iv)
                  end do
                  do iv=(np+1),np_max,1
                     idxp_compo(it,in,iy,iv) = -1
                     tabp_compo(it,in,iy,iv) = 0.d00
                  end do
                  do iw=1,nq,1
                     idxq_compo(it,in,iy,iw) = iqtyq(iw)
                     tabq_compo(it,in,iy,3*(iw-1)+1) = qtyq(iw,1)
                     tabq_compo(it,in,iy,3*(iw-1)+2) = qtyq(iw,2)
                     tabq_compo(it,in,iy,3*(iw-1)+3) = qtyq(iw,3)
                  end do
                  do iw=(nq+1),nq_max,1
                     idxq_compo(it,in,iy,iw) = -1
                  end do
                  !                     icnt = icnt+1
               end if
            end if
         end do
      end if
      close(unit=iunit)
   end do
   if (ic == 1) then
      if ((idx_ex(2) == 1).and.(iwr == 1)) then
         write(*,*) icnt,' entries of composition table read'
      end if
      if (np_max > dim_qtyp) then
         np_max = dim_qtyp
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 22
      end if
      if (nq_max > dim_qtyq) then
         nq_max = dim_qtyq
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 23
      end if
      if ((idx_ex(2) == 1).and.(iwr == 1)) then
         write(*,*) 'maximum number of pairs in eos.compo      =',np_max
         write(*,*) 'maximum number of quadruples in eos.compo =',nq_max
      end if
      if (ierr == 0) then
         allocate(idxp_compo(0:dim_t,1:dim_n,0:dim_y,np_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 24
         end if
         allocate(idxq_compo(0:dim_t,1:dim_n,0:dim_y,nq_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 24
         end if
         allocate(tabp_compo(0:dim_t,1:dim_n,0:dim_y,np_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 24
         end if
         allocate(tabq_compo(0:dim_t,1:dim_n,0:dim_y,3*nq_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 24
         end if
      end if
      if (ierr == 0) then
         ! initialization
         do it=0,dim_t,1
            do in=1,dim_n,1
               do iy=0,dim_y,1
                  idx_arg(it,in,iy,2) = 0
                  do iv=1,np_max,1
                     idxp_compo(it,in,iy,iv) = -1
                  end do
                  do iw=1,nq_max,1
                     idxq_compo(it,in,iy,iw) = -1
                  end do
                  do iv=1,np_max,1
                     tabp_compo(it,in,iy,iv) = 0.d00
                  end do
                  do iw=1,3*nq_max,1
                     tabq_compo(it,in,iy,iw) = 0.d00
                  end do
               end do
            end do
         end do
      end if
   end if
end do

call write_errors(ierr)

return
end SUBROUTINE read_eos_table_compo
!***********************************************************************
SUBROUTINE read_eos_table_micro(iwr,nbl)
! Stefan Typel for the CompOSE core team, version 0.09, 2012/10/02
USE compose_internal
implicit none
integer :: iwr,nbl
integer :: ibl,ierror,icnt,it,in,iy,iunit,iflag,itmp,ierr,ic,&
     iv,nm,alloc_status
integer :: iqtym(dim_qtym)
double precision :: qtym(dim_qtym)

! counter for total number of read table entries
icnt = 0

! error counter
ierr = 0

! reading cycles
idx_ex(3) = 0
do ic=1,2,1
   if (ic == 1) nm_max = 0
   open(unit=20,file='eos.micro',&
        status='old',action='read',iostat=ierror)
   if (ierror == 0) then
      itmp = 0
   else
      itmp = nbl-1
   end if
   
   do ibl=0,itmp,1
      ! flag for stop reading file
      iflag = 0
      iunit = 20+ibl
      
      choose_file: select case (ibl)
      case (0)
         if (itmp > 0) then
            open(unit=iunit,file='eos.micro0',&
                 status='old',action='read',iostat=ierror)
         end if
      case (1)
         open(unit=iunit,file='eos.micro1',&
              status='old',action='read',iostat=ierror)
      case (2)
         open(unit=iunit,file='eos.micro2',&
              status='old',action='read',iostat=ierror)
      case (3)
         open(unit=iunit,file='eos.micro3',&
              status='old',action='read',iostat=ierror)
      case (4)
         open(unit=iunit,file='eos.micro4',&
              status='old',action='read',iostat=ierror)
      case (5)
         open(unit=iunit,file='eos.micro5',&
              status='old',action='read',iostat=ierror)
      case (6)
         open(unit=iunit,file='eos.micro6',&
              status='old',action='read',iostat=ierror)
      case (7)
         open(unit=iunit,file='eos.micro7',&
              status='old',action='read',iostat=ierror)
      case (8)
         open(unit=iunit,file='eos.micro8',&
              status='old',action='read',iostat=ierror)
      case (9)
         open(unit=iunit,file='eos.micro9',&
              status='old',action='read',iostat=ierror)
      end select choose_file
      
      if (ierror == 0) then
         idx_ex(3) = 1
         if ((ic == 1).and.(iwr == 1)) then
            write(*,*)
            write(*,*) 'reading eos table(s) with microscopic information'
         end if
         do while (iflag == 0)
            read(iunit,*,iostat=ierror) it,in,iy,&
                 nm,(iqtym(iv),qtym(iv),iv=1,nm,1)
            if (ic == 1) then
               if (ierror /= 0) then 
                  iflag = 1
               else 
                  if (nm > nm_max) nm_max = nm
                  icnt = icnt+1
               end if
            else
               if (ierror /= 0) then 
                  iflag = 1
               else 
                  idx_arg(it,in,iy,3) = 1
                  do iv=1,nm,1
                     idx_mic(it,in,iy,iv) = iqtym(iv)
                     tab_mic(it,in,iy,iv) = qtym(iv)
                  end do
                  do iv=(nm+1),nm_max,1
                     idx_mic(it,in,iy,iv) = -1
                     tab_mic(it,in,iy,iv) = 0.d00
                  end do
                  !                     icnt = icnt+1
               end if
            end if
         end do
      end if
      
      close(unit=iunit)
   end do
   if (ic == 1) then
      if ((idx_ex(3) == 1).and.(iwr == 1)) then
         write(*,*) icnt,' entries of microscopic table read'
      end if
      if (nm_max < 1) nm_max = 1
      if (nm_max > dim_qtym) then
         nm_max = dim_qtym
         ierr = ierr+1
         if (ierr < dim_err) error_msg(ierr) = 25
      end if
      if ((idx_ex(3) == 1).and.(iwr == 1)) then
         write(*,*) 'maximum number of pairs in eos.micro      =',nm_max
      end if
      if (ierr == 0) then
         allocate(idx_mic(0:dim_t,1:dim_n,0:dim_y,nm_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 26
         end if
         allocate(tab_mic(0:dim_t,1:dim_n,0:dim_y,nm_max),stat=alloc_status)
         if (alloc_status /= 0) then
            ierr = ierr+1
            if (ierr < dim_err) error_msg(ierr) = 26
         end if
      end if
      if (ierr == 0) then
         ! initialization
         do it=0,dim_t,1
            do in=1,dim_n,1
               do iy=0,dim_y,1
                  idx_arg(it,in,iy,3) = 0
                  do iv=1,nm_max,1
                     idx_mic(it,in,iy,iv) = -1
                     tab_mic(it,in,iy,iv) = 0.d00
                  end do
               end do
            end do
         end do
      end if
   end if
end do

call write_errors(ierr)

return
end SUBROUTINE read_eos_table_micro
!***********************************************************************
SUBROUTINE write_errors(ierr)
! Stefan Typel for the CompOSE core team, version 0.10, 2012/11/01
USE compose_internal
implicit none
integer :: ierr,i,imax

error_msg(0) = ierr
if (ierr > 0) then
   write(*,*) '!!! warning !!!'
   write(*,*) ierr,' error(s) detected'
   imax = ierr
   if (imax > dim_err) imax = dim_err
   do i=1,imax,1
      write(*,*) 'error(',i,') =',error_msg(i)

      choose_error: select case (error_msg(i))
      case (1)
         write(*,*) 'input file eos.t does not exist'
      case (2)
         write(*,*) 'number in first line of file eos.t too small'
      case (3)
         write(*,*) 'number in second line of file eos.t too large'
      case (4)
         write(*,*) 'maximum index smaller than minimum index in eos.t'
      case (5)
         write(*,*) 'temperature in eos.t not increasing'
      case (6)
         write(*,*) 'input file eos.nb does not exist'
      case (7)
         write(*,*) 'number in first line of file eos.nb too small'
      case (8)
         write(*,*) 'number in second line of file eos.nb too large'
      case (9)
         write(*,*) 'maximum index smaller than minimum index in eos.nb'
      case (10)
         write(*,*) 'baryon number density in eos.nb not increasing'
      case (11)
         write(*,*) 'input file eos.yq does not exist'
      case (12)
         write(*,*) 'number in first line of file eos.yq too small'
      case (13)
         write(*,*) 'number in second line of file eos.yq too large'
      case (14)
         write(*,*) 'maximum index smaller than minimum index in eos.yq'
      case (15)
         write(*,*) 'baryon charge fraction in eos.yq not increasing'
      case (16)
         write(*,*) 'baryon charge fraction indices out of range'
      case (20)
         write(*,*) 'quantity index out of range in eos.thermo'
      case (21)
         write(*,*) 'number of additional quantities in eos.thermo out of range'
      case (22)
         write(*,*) 'number of pairs in eos.compo out of range'
      case (23)
         write(*,*) 'number of quadruples in eos.compo out of range'
      case (24)
         write(*,*) 'error in allocation of memory for compo tables'
      case (25)
         write(*,*) 'number of pairs in eos.micro out of range'
      case (26)
         write(*,*) 'error in allocation of memory for micro tables'
      case (27)
         write(*,*) 'error in reading masses from file eos.thermo'
      case (30)
         write(*,*) 'wrong temperature in parameter table for T=0'
      case (31)
         write(*,*) 'wrong asymmetry in parameter table'
      case (41)
         write(*,*) 'temperature out of range'
      case (42)
         write(*,*) 'baryon number density out of range'
      case (43)
         write(*,*) 'baryonic charge fraction out of range'
      case (51)
         write(*,*) 'not all points of interpolation line exist'
      case (52)
         write(*,*) 'not all points of interpolation square exist'
      case (53)
         write(*,*) 'not all points of interpolation cube exist'
      case (60)
         write(*,*) 'input file eos.quantities does not exist'
      case (61)
         write(*,*) 'number of quantities out of range'
      case (62)
         write(*,*) 'number of additional quantities out of range'
      case (63)
         write(*,*) 'index of quantities out of range'
      case (64)
         write(*,*) 'index of additional quantities out of range'
      case (65)
         write(*,*) 'number of pairs out of range'
      case (66)
         write(*,*) 'number of quadruples out of range'
      case (67)
         write(*,*) 'number of microscopic quantities out of range'
      case (68)
         write(*,*) 'number of error quantities out of range'
      case (69)
         write(*,*) 'index for output format out of range'
      case (70)
         write(*,*) 'index for interpolation in T out of range'
      case (71)
         write(*,*) 'index for interpolation in n_b out of range'
      case (72)
         write(*,*) 'index for interpolation in Y_q out of range'
      case (80)
         write(*,*) 'number of table entries too small'
      case (81)
         write(*,*) 'minimum of parameter T < 0'
      case (82)
         write(*,*) 'minimum of parameter T <= 0'
      case (83)
         write(*,*) 'maximum of parameter T < minimum of parameter T'
      case (84)
         write(*,*) 'minimum of parameter n_b <= 0'
      case (85)
         write(*,*) 'maximum of parameter T < minimum of parameter n_b'
      case (86)
         write(*,*) 'minimum of parameter Y_q < 0'
      case (87)
         write(*,*) 'maximum of parameter Y_q < minimum of parameter Y_q'
      case (88)
         write(*,*) 'maximum of parameter Y_q > 1'
      case (90)
         write(*,*) 'option for beta equilibrium not available'
      end select choose_error

   end do
   write(*,*) 'program terminated'
   write(*,*)
   stop
end if

return
end SUBROUTINE write_errors
!***********************************************************************
INTEGER FUNCTION get_ipl_rule(ir,idx,ipl)
! Stefan Typel for the CompOSE core team, version 0.01, 2012/04/05
USE compose_internal
implicit none
integer :: ir,idx,ipl

if (idx > 0) then
   get_ipl_rule = ipl_rule(1,ipl,ir)
else
   get_ipl_rule = ipl_rule(0,ipl,ir)
end if

return
end FUNCTION get_ipl_rule
!***********************************************************************
