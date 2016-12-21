MODULE precision
      
  integer,parameter :: dp = selected_real_kind(15)

END MODULE precision


MODULE hdfparameters

  USE precision

  real(dp), allocatable :: nb_hdf5(:),t_hdf5(:),y_q_hdf5(:)
  integer :: m_nb,n_temp,o_y_q

  real(dp), allocatable :: thermo_hdf5(:,:,:,:), thermo_hdf5_add(:,:,:,:)
  real(dp), allocatable :: yi_hdf5(:,:,:,:),aav_hdf5(:,:,:,:),zav_hdf5(:,:,:,:),yav_hdf5(:,:,:,:)
  real(dp), allocatable :: micro_hdf5(:,:,:,:),err_hdf5(:,:,:,:)
  integer,allocatable :: index_thermo(:),index_thermo_add(:),index_err(:)
  integer, allocatable :: index_yi(:),index_av(:),index_micro(:)

contains 

  subroutine initialise_hdf5(n_nb,n_t,n_yq,i_tab)

    USE compose_internal 
    IMPLICIT NONE
    
    integer n_nb,n_t,n_yq,i_tab
    
    m_nb = n_nb
    n_temp = n_t
    o_y_q = n_yq

! memory allocation for the data arrays, thermo

    IF(n_qty.ne.0) then
       allocate(thermo_hdf5(m_nb,n_temp,o_y_q,n_qty))
       thermo_hdf5 = 0._dp
       allocate(index_thermo(n_qty))
       index_thermo = 0
    else
       write(*,*) 'n_qty equal to zero, no thermodynamical quantities'
    end IF
    IF(n_add.ne.0) then
       allocate(thermo_hdf5_add(m_nb,n_temp,o_y_q,n_add))
       thermo_hdf5_add = 0._dp
       allocate(index_thermo_add(n_add))
       index_thermo_add = 0
    else
       write(*,*) 'n_add equal to zero, no additional thermodynamic quantities'
    end IF

! memory allocation for the data arrays, composition

    IF(n_p.ne.0) then
       allocate(yi_hdf5(m_nb,n_temp,o_y_q,n_p))
       allocate(index_yi(n_p))
       yi_hdf5 = 0._dp
       index_yi = 0
    else
       write(*,*) 'n_p equal to zero, no pairs'
    end if
    IF(n_q.ne.0) Then
       allocate(yav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(aav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(zav_hdf5(m_nb,n_temp,o_y_q,n_q))
       allocate(index_av(n_q))
       yav_hdf5 = 0._dp
       aav_hdf5 = 0._dp
       zav_hdf5 = 0._dp
       index_av = 0
    else
       write(*,*) 'n_q equal to zero, no quadruples'
    end if

! memory allocation for the data arrays, microscopic
    IF(n_m.ne.0) THEN
       allocate(micro_hdf5(m_nb,n_temp,o_y_q,n_m))
       allocate(index_micro(n_m))
       micro_hdf5 = 0._dp
       index_micro = 0
    else
       write(*,*) 'n_m equal to zero, no microscopic quantities'
    end if

    IF(n_err.ne.0) THEN
       allocate(err_hdf5(m_nb,n_temp,o_y_q,n_err))
       allocate(index_err(n_err))
       err_hdf5 = 0._dp
       index_err = 0
    else
       write(*,*) 'n_m equal to zero, no error quantities'
    end if


    IF(i_tab.eq.0) then
       allocate(nb_hdf5(m_nb))
       allocate(t_hdf5(m_nb))
       allocate(y_q_hdf5(m_nb))
    else
       allocate(nb_hdf5(m_nb))
       allocate(t_hdf5(n_temp))
       allocate(y_q_hdf5(o_y_q))
    end IF
    nb_hdf5 = 0._dp
    t_hdf5 = 0._dp
    y_q_hdf5 = 0._dp
  end subroutine initialise_hdf5

  subroutine write_hdf5(itab)

    USE compose_internal
    IMPLICIT NONE
    
    integer itab

    integer :: openfile = 1 ! =1 if file for writing hdf5 has not been opened
                            ! =2 if file is already open
                            ! the temperature, density and y_q arrays are 
                            ! only written the first time (i.e. openfile = 1)
    
    write(*,*) 'writing ',n_qty,' thermodynamic quantities into file '

    call hdf5_write_thermo(m_nb,n_temp,o_y_q, n_qty, nb_hdf5, t_hdf5, y_q_hdf5, &
                             thermo_hdf5,index_thermo,openfile,itab)
    openfile = 2

    write(*,*) 'writing ',n_add,' additional thermodynamic quantities into file'
    call hdf5_write_thermo_add(m_nb,n_temp,o_y_q, n_add, nb_hdf5, t_hdf5, y_q_hdf5,&
                             thermo_hdf5_add,index_thermo_add,openfile,itab)

    write(*,*) 'writing ',n_p,' pairs into file'
      
    call hdf5_write_compo_p(m_nb,n_temp,o_y_q, n_p, nb_hdf5, t_hdf5, y_q_hdf5, &
                     yi_hdf5,index_yi,openfile,itab)
    
    write(*,*) 'writing ',n_q,' quadruples into file'

    call hdf5_write_compo_q(m_nb,n_temp,o_y_q, n_q, nb_hdf5, t_hdf5, y_q_hdf5, &
                              yav_hdf5,aav_hdf5,zav_hdf5,index_av,openfile,itab)
    write(*,*) 'writing ',n_m,' microscopic quantities into file'
    call hdf5_write_micro(m_nb,n_temp,o_y_q, n_m, nb_hdf5, t_hdf5, y_q_hdf5, &
                              micro_hdf5,index_micro,openfile,itab)


    write(*,*) 'writing ',n_err,' error quantities into file'
    call hdf5_write_err(m_nb,n_temp,o_y_q, n_err, nb_hdf5, t_hdf5, y_q_hdf5, &
                     err_hdf5,index_err,openfile,itab)


  end subroutine write_hdf5


  subroutine read_hdf5(itab)

    USE compose_internal
    IMPLICIT NONE

    integer itab
    
    integer :: openfile = 1 ! =1 if file for reading hdf5 has not been opened
                            ! =2 if file is already open
                            ! the temperature, density and y_q arrays are 
                            ! only read the first time (i.e. openfile = 1)
    character(30) n_qty_name,qty_name,index_name

    integer m_nb_rd,n_temp_rd,o_y_q_rd,n_qty_rd,n_add_rd,n_p_rd,n_q_rd,n_m_rd,n_err_rd,i,j,k,i_tab_rd

    ! First read the number of different quantities stored

    call hdf5_read_dimensions(m_nb_rd,n_temp_rd,o_y_q_rd,n_qty_rd,n_add_rd,&
                              n_p_rd,n_q_rd,n_m_rd,n_err_rd,i_tab_rd,openfile)

    IF(itab.ne.i_tab_rd) then
       write(*,*) 'error reading file, i_tab is',i_tab_rd,' not ',itab
       stop
    end IF

    ! If necessary you can check here if the numbers read are in agreement
    ! with your expectation. If not, just initialise the dimensions
    ! For the moment we use it for testing purposes, thus we check that
    ! everything has the correct dimensions
    IF(m_nb.ne.m_nb_rd) then
       write(*,*) 'error reading the number of points in nb',m_nb,m_nb_rd
       stop
    end IF
    IF(itab.eq.0) then
       IF(m_nb.ne.n_temp_rd) then
          write(*,*) 'error reading the number of points in T',m_nb,n_temp_rd
          stop
       end IF
    else
       IF(n_temp.ne.n_temp_rd) then
          write(*,*) 'error reading the number of points in T',n_temp,n_temp_rd
          stop
       end IF
    end IF
    IF(itab.eq.0) then
       IF(m_nb.ne.o_y_q_rd) then
          write(*,*) 'error reading the number of points in yq',m_nb,o_y_q_rd
          stop
       end IF
    else
       IF(o_y_q.ne.o_y_q_rd) then
          write(*,*) 'error reading the number of points in yq',o_y_q,o_y_q_rd
          stop
       end IF
    end IF
    IF(n_qty.ne.n_qty_rd) then
       write(*,*) 'error reading the number of thermodynamical quantities',n_qty,n_qty_rd
       stop
    end IF
    IF(n_add.ne.n_add_rd) then
       write(*,*) 'error reading the number of additional thermo quantities',n_add,n_add_rd
       stop
    end IF
    IF(n_p.ne.n_p_rd) then
       write(*,*) 'error reading the number of pairs',n_p,n_p_rd
       stop
    end IF
    IF(n_q.ne.n_q_rd) then
       write(*,*) 'error reading the number of quadruples',n_q,n_q_rd
       stop
    end IF
    IF(n_m.ne.n_m_rd) then
       write(*,*) 'error reading the number of microscopic quantities',n_m,n_m_rd
       stop
    end IF
    IF(n_err.ne.n_err_rd) then
       write(*,*) 'error reading the number of error quantities',n_err,n_err_rd
       stop
    end IF

    ! Part not used for the moment
!!$    m_nb = m_nb_rd
!!$    n_temp = n_temp_rd
!!$    o_y_q = o_y_q_rd
!!$    n_qty = n_qty_rd
!!$    n_add = n_add_rd
!!$    n_p = n_p_rd
!!$    n_q = n_q_rd
!!$    n_m = n_m_rd
!!$    n_err = n_err_rd

    IF(itab == 0) then
       call initialise_hdf5(m_nb_rd,1,1,itab)
    else
       call initialise_hdf5(m_nb_rd,n_temp_rd,o_y_q_rd,itab)
    end IF


    IF(n_qty.ne.0) then
       write(*,*) 'reading ',n_qty,' thermodynamic quantities from file'
       qty_name = 'thermo'
       qty_name = trim(qty_name)//CHAR(0)
       index_name = 'index_thermo'
       index_name = trim(index_name)//CHAR(0)
       call hdf5_read(nb_hdf5, t_hdf5, y_q_hdf5, &
                             thermo_hdf5,index_thermo,openfile,&
                            qty_name,index_name)

       openfile = 2
      
    end IF

    IF(n_add.ne.0) then
       write(*,*) 'reading ',n_add,' additional thermodynamic quantities from file'
       qty_name = 'thermo_add'
       qty_name = trim(qty_name)//CHAR(0)
       index_name = 'index_thermo_add'
       index_name = trim(index_name)//CHAR(0)
       call hdf5_read(nb_hdf5, t_hdf5, y_q_hdf5, &
                             thermo_hdf5_add,index_thermo_add,openfile,&
                            qty_name,index_name)
       openfile = 2
    end IF

    IF(n_p.ne.0) then
       write(*,*) 'reading ',n_p,' pairs from file'
       qty_name = 'yi'
       qty_name = trim(qty_name)//CHAR(0)
       index_name = 'index_yi'
       index_name = trim(index_name)//CHAR(0)
       call hdf5_read(nb_hdf5, t_hdf5, y_q_hdf5, &
                             yi_hdf5,index_yi,openfile,&
                            qty_name,index_name)
       openfile = 2
    end IF

    IF(n_q.ne.0) then
       write(*,*) 'reading ',n_q,' quadruples from file'
       call hdf5_read_av(nb_hdf5, t_hdf5, y_q_hdf5, &
                              yav_hdf5,aav_hdf5,zav_hdf5,index_av,openfile)
       openfile = 2
    end IF

    IF(n_m.ne.0) then
       write(*,*) 'reading ',n_m,' microscopic quantities from file'
       qty_name = 'micro'
       qty_name = trim(qty_name)//CHAR(0)
       index_name = 'index_micro'
       index_name = trim(index_name)//CHAR(0)
       call hdf5_read(nb_hdf5, t_hdf5, y_q_hdf5, &
                             micro_hdf5,index_micro,openfile,&
                            qty_name,index_name)
       openfile = 2
    end IF


    IF(n_err.ne.0) then
       write(*,*) 'reading ',n_err,' error quantities from file'
       qty_name = 'error'
       qty_name = trim(qty_name)//CHAR(0)
       index_name = 'index_err'
       index_name = trim(index_name)//CHAR(0)
       call hdf5_read(nb_hdf5, t_hdf5, y_q_hdf5, &
                             err_hdf5,index_err,openfile,&
                            qty_name,index_name)
       openfile = 2
    end IF

    open(33,file='readtest.d',status='unknown')
    IF(itab.eq.0) then
       DO i = 1,m_nb
          write(33,*) nb_hdf5(i),t_hdf5(i),y_q_hdf5(i),thermo_hdf5(i,1,1,1:n_qty),yi_hdf5(i,1,1,1:n_p),&
            aav_hdf5(i,1,1,1:n_q),zav_hdf5(i,1,1,1:n_q),yav_hdf5(i,1,1,1:n_q)&
            ,micro_hdf5(i,1,1,1:n_m),err_hdf5(i,1,1,1:n_err)
       end DO
    else
       DO i = 1,n_temp
          DO j = 1,m_nb
             DO k = 1,o_y_q
                write(33,*) t_hdf5(i),nb_hdf5(j),y_q_hdf5(k),thermo_hdf5(j,i,k,1:n_qty),yi_hdf5(j,i,k,1:n_p),&
            aav_hdf5(j,i,k,1:n_q),zav_hdf5(j,i,k,1:n_q),yav_hdf5(j,i,k,1:n_q)&
            ,micro_hdf5(j,i,k,1:n_m),err_hdf5(j,i,k,1:n_err)
             end DO
          end DO
       end DO
    end IF

  end subroutine read_hdf5

  subroutine close_hdf5

    IMPLICIT NONE

    IF(allocated(thermo_hdf5)) deallocate(thermo_hdf5)
    IF(allocated(thermo_hdf5_add)) deallocate(thermo_hdf5_add)
    IF(allocated(yi_hdf5)) deallocate(yi_hdf5)
    IF(allocated(yav_hdf5)) deallocate(yav_hdf5)
    IF(allocated(zav_hdf5)) deallocate(zav_hdf5)
    IF(allocated(aav_hdf5)) deallocate(aav_hdf5)
    IF(allocated(micro_hdf5)) deallocate(micro_hdf5)
    IF(allocated(err_hdf5)) deallocate(err_hdf5)
    IF(allocated(index_thermo)) deallocate(index_thermo)
    IF(allocated(index_thermo_add)) deallocate(index_thermo_add)
    IF(allocated(index_err)) deallocate(index_err)
    IF(allocated(index_yi)) deallocate(index_yi)
    IF(allocated(index_av)) deallocate(index_av)
    IF(allocated(index_micro)) deallocate(index_micro)
    IF(allocated(nb_hdf5)) deallocate(nb_hdf5)
    IF(allocated(t_hdf5)) deallocate(t_hdf5)
    IF(allocated(y_q_hdf5)) deallocate(y_q_hdf5)

  end subroutine close_hdf5

end MODULE hdfparameters
