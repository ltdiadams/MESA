! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none

      ! these routines are called by the standard run_star check_model
      contains

      !include 'standard_run_star_extras.inc'
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

!-------------------------------------------------------------------------------
         s% use_other_mlt = .true.
         !s% use_other_eos = .true.
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_mlt => my_other_mlt

         !s% other_eosDT_get => my_eosDT_get
         !s% other_eosDT_get_T => my_eosDT_get_T
         !s% other_eosDT_get_Rho => my_eosDT_get_Rho
         !s% other_eosPT_get => my_eosPT_get
         !s% other_eosPT_get_T => my_eosPT_get_T
         !s% other_eosPT_get_Pgas => my_eosPT_get_Pgas
         !s% other_eosPT_get_Pgas_for_Rho => my_eosPT_get_Pgas_for_Rho

!-------------------------------------------------------------------------------
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.

      end subroutine extras_controls

      ! None of the following functions are called unless you set their
      ! function point in extras_control.


      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup



! extra subroutine from other_mlt
!------------------------------------------------------------------------------------------

      subroutine my_other_mlt(  &
            id, k, cgrav, m, mstar, r, L, X, &
            T_face, rho_face, P_face, &
            chiRho_face, chiT_face, &
            Cp_face, opacity_face, grada_face, &
            alfa, beta, d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, & ! f_face = alfa*f_00 + beta*f_m1
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &
            gradr_factor, gradL_composition_term, &
            alpha_semiconvection, semiconvection_option, &
            thermohaline_coeff, thermohaline_option, &
            dominant_iso_for_thermohaline, &
            mixing_length_alpha, alt_scale_height, remove_small_D_limit, &
            MLT_option, Henyey_y_param, Henyey_nu_param, &
            gradT_smooth_low, gradT_smooth_high, smooth_gradT, &
            normal_mlt_gradT_factor, &
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, just_gradr, &
            mixing_type, mlt_basics, mlt_partials1, ierr)
         use star_lib, only: star_mlt_eval
         integer, intent(in) :: id ! id for star

         integer, intent(in) :: k ! cell number or 0 if not for a particular cell
         !integer :: k

         real(dp), intent(in) :: &
            cgrav, m, mstar, r, L, X, &
            T_face, rho_face, P_face, &
            chiRho_face, chiT_face, &
            Cp_face, opacity_face, grada_face, &
            alfa, beta, d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, &
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &
            gradr_factor, gradL_composition_term, &
            alpha_semiconvection, thermohaline_coeff, mixing_length_alpha, &
            Henyey_y_param, Henyey_nu_param, &
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, remove_small_D_limit, &
            gradT_smooth_low, gradT_smooth_high, normal_mlt_gradT_factor
         logical, intent(in) :: alt_scale_height
         character (len=*), intent(in) :: thermohaline_option, MLT_option, semiconvection_option
         integer, intent(in) :: dominant_iso_for_thermohaline
         logical, intent(in) :: just_gradr, smooth_gradT
         integer, intent(out) :: mixing_type
         real(dp), intent(inout) :: mlt_basics(:) ! (num_mlt_results)
         real(dp), intent(inout), pointer :: mlt_partials1(:) ! =(num_mlt_partials, num_mlt_results)
         integer, intent(out) :: ierr
         !ierr = 0
         !real(8), parameter :: PI_8 = 4 * atan (1.0_8)
         real :: B0, mu_not, P_mag
         !real, dimension(:), allocatable :: msv
         integer :: factor, i

      type (star_info), pointer :: s

      ierr = 0
      !mixing_type = 0

      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

        B0 = s% x_ctrl(1)
        !mu_not = 4*PI_8*1.e-7
        !mu_not = 4*pi*1.e-7*(1.11259408*1.e-15) ! converting from H/m to cgs units
        mu_not = (4*pi*1.e-7)*((1.e7)/(4*pi))

        P_mag = (B0**2)/(2*mu_not)
        !write(*,*) 'Pinitial', s% P(k)

        !s% P(k) = s% P(k) + P_mag kinda works...


        !if (P_mag /= 0) then
          !do i=1, s% nz
            !s% P(k) = P_face + P_mag

            s% P(k) = s% Prad(k) + s% Pgas(k) + P_mag
            !s% P(k) = s% P(k) + P_mag !same as above

            !write(*,*) 'Ploop', s% P(i)
          !enddo
        !endif
        !write(*,*) 'Pafter', s% P(k)
        !write(*,*) 'P_face', P_face
        !write(*,*) 'P_mag', P_mag
        !write(*,*) 'P_face + P_mag', P_face + P_mag


      !write(*,*) 'Pressure', s% P(k)

         call star_mlt_eval(  &
            id, k, cgrav, m, mstar, r, L, X, &
            T_face, rho_face, P_face, &
            chiRho_face, chiT_face, &
            Cp_face, opacity_face, grada_face, &
            alfa, beta, d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, & ! f_face = alfa*f_00 + beta*f_m1
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &
            gradr_factor, gradL_composition_term, &
            alpha_semiconvection, semiconvection_option, &
            thermohaline_coeff, thermohaline_option, &
            dominant_iso_for_thermohaline, &
            mixing_length_alpha, alt_scale_height, remove_small_D_limit, &
            MLT_option, Henyey_y_param, Henyey_nu_param, &
            gradT_smooth_low, gradT_smooth_high, smooth_gradT, &
            normal_mlt_gradT_factor, &
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, just_gradr, &
            mixing_type, mlt_basics, mlt_partials1, ierr)
      end subroutine my_other_mlt


      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 1 !changed from 0 to 1
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)

         !--------------------------------------
         ! need to use the star mass defined from inlist? to find the mean value
         ! of the magnetic field at the magnetic equator given by: ((mu_not * mass)/(4*pi*(radius^3)))


         !real :: radius, m_star, mu_not, big_constant, r_star, mean_mag
         real(8), parameter :: PI_8 = 4 * atan (1.0_8)
         real :: B0, mu_not, P_mag
         !real, dimension(:), allocatable :: msv
         integer :: factor, i
         !--------------------------------------

         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         B0 = s% x_ctrl(1)
         !mu_not = 4*PI_8*1.e-7*(1.11259408*1.e-15) ! converting from H/m to cgs units
         mu_not = (4*pi*1.e-7)*((1.e7)/(4*pi))

         P_mag = (B0**2)/(2*mu_not)

         names(1) = 'P_mag'
         vals(1) = P_mag


         !mu_not = s% center_mu

         !r_star = 10.0**s % log_surface_radius !units of msun
         !m_star = s% mstar !in g

         !big_constant = mu_not/4*PI_8
         !factor = m_star/radius**3

         ! B0 = (mu_not*mass)/(4*pi*radius^3)

         !mean_mag = big_constant*factor

         !names(1) = 'r_star'
         !names(2) = 'm_star'
         !names(3) = 'B0_outer_edge'
         !vals(1) = r_star
         !vals(2) = m_star
         !vals(3) = mean_mag


         !vals(1) = value/target_rad  ! compactness = 2.5 msun / radius (1000 km)
         !vals(1) = factor*big_constant
         !vals(1) = factor*big_constant

        !note: do NOT add the extras names to history_columns.list
        ! the history_columns.list is only for the built-in log column options.
        ! it must not include the new column names you are adding here.

        ! deallocate(msv)

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: B0, mu_not, P_mag, P_div, P_with, P_without
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = "P"
         names(2) = "P_without"
         !names(2) = "P_after"
         names(3) = "P_div"

         call star_ptr(id, s, ierr)
         if (ierr /=0) return

         B0 = s% x_ctrl(1)
         !mu_not = (4*pi*1.e-7)!*(1.11259408*1.e-15) ! converting from H/m to cgs units
         mu_not = (4*pi*1.e-7)*((1.e7)/(4*pi))

         P_mag = (B0**2)/(2*mu_not)

         do k=1, s% nz

           P_without = s% P(k) - P_mag
           P_with = s% P(k)
           !s% P(k) = s% P(k) + P_mag

           !P_div = (s% P(k) - P_mag) / (s% P(k))
           P_div = P_without / P_with

           !vals(k,1) = s% P(k)
           vals(k,1) = P_with

           !vals(k,2) = s% P(k) - P_mag
           vals(k,2) = P_without

           vals(k,3) = P_div
           !vals(k,2) = s% P(k) + P_mag
           !vals(k,3) = s% P(k) / (s% P(k) + P_mag)

         enddo

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols=0
      end subroutine how_many_extra_history_header_items

      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra history header item
      !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items

      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         !-----------------------------------------------------------------------
         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
