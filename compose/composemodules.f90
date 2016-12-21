!***********************************************************************
!***********************************************************************
! modules
!***********************************************************************
!***********************************************************************
MODULE eos_tables
! Stefan Typel for the CompOSE core team, version 0.18, 2013/02/06

! parameters
! dim_t:     maximum index number in temperature
! dim_n:     maximum index number in baryon number density
! dim_y:     maximum index number in baryonic charge fraction
! dim_a      = Max(dim_t,dim_n,dim_y)
! dim_reg:   number of standard thermodynamic quantities in eos.thermo (7)
! dim_qty:   maximum number of additional quantities in eos.thermo
! dim_qtyp:  maximum number of different pairs in eos.compo
! dim_qtyq:  maximum number of different quadruples in eos.compo
! dim_qtym:  maximum number of different quantities in eos.micro
! dim_qtye:  maximum number of error quantities
! dim_qty2   = dim_qty+dim_reg
! dim_qtyt   = 19 = regular number of quantities in eos_thermo
! dim_qtyq2  = 3*dim_qtyq
! dim_ipl    = Max(dim_qty2,dim_qtyp,dim_qtyq,dim_qtym)
!            = maximum number of interpolated quantities

implicit none

integer, parameter :: dim_t=201,dim_n=401,dim_y=100,dim_a=401,&
     dim_reg=7,dim_qty=15,dim_qtyp=15,dim_qtyq=2,&
     dim_qty2=21,dim_qtyq2=6,dim_qtyt=19,&
     dim_qtym=16,dim_qtye=8,dim_ipl=21

integer :: idx_min(3),idx_max(3),&
     nadd_max,np_max,nq_max,nm_max,&
     idx_thermo(0:dim_t,1:dim_n,0:dim_y,dim_qty),&
     idx_compo_p(dim_qtyp),idx_compo_q(dim_qtyq),&
     idx_micro(dim_qtym)

integer, allocatable :: idxp_compo(:,:,:,:),idxq_compo(:,:,:,:),&
     idx_mic(:,:,:,:)

double precision :: tab_para(0:dim_a,6),&
     tab_thermo(0:dim_t,1:dim_n,0:dim_y,dim_qty2),&
     eos_thermo(dim_qtyt),eos_thermo_add(dim_qty),&
     eos_compo_p(0:dim_qtyp),eos_compo_q(dim_qtyq,3),&
     eos_micro(0:dim_qtym),eos_thermo_err(dim_qtye)

double precision, allocatable :: tabp_compo(:,:,:,:),&
     tabq_compo(:,:,:,:),tab_mic(:,:,:,:)

! neutron mass in MeV/c^{2}, standard value
double precision :: m_n =  939.565379d00

! proton mass in MeV/c^{2}, standard value
double precision :: m_p =  938.272046d00

END MODULE eos_tables
!***********************************************************************
MODULE compose_internal
! Stefan Typel for the CompOSE core team, version 0.09, 2012/12/18
USE eos_tables
implicit none

integer, parameter :: dim_err=100 

integer :: eos_dim,eos_type,idx_ex(0:3),idx_ipl(3),&
     idx_arg2(-4:5,-4:5,0:2),idx_arg1(-4:5,0:2),&
     idx_arg(0:dim_t,1:dim_n,0:dim_y,0:4),&
     idx_argx(0:dim_t,1:dim_n,0:dim_y,3),&
     ipl_rule(0:1,3,-8:8),&
     error_msg(0:dim_err),&
     n_qty,n_add,idx_qty(dim_qtyt),idx_add(dim_qty),&
     n_p,n_q,idx_p(dim_qtyp),idx_q(dim_qtyq),&
     n_m,idx_m(dim_qtym),&
     n_err,idx_err(dim_qtye),iout,&
     incl_l

double precision :: d1x(-4:5,-4:5,-4:4),d2x(-4:5,-4:5,-4:4),&
     d1y(-4:5,-4:5,-4:4),d2y(-4:5,-4:5,-4:4),dx,dy,dx2,dy2,&
     r1d(3,-1:dim_a,-4:4,-8:8),r2d(3,-1:dim_a,-4:4,-8:8),&
     df(0:3,0:3,-4:5,-4:5),fc(0:5,0:5)

END MODULE compose_internal
!***********************************************************************
