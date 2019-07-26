! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************


      module mlt_info

      use const_def
      use num_lib
      use utils_lib
      use star_private_def

      implicit none

      private
      public :: set_mlt_vars, do1_mlt, do1_mlt_eval, set_grads, Get_results, &
         set_gradT_excess_alpha, adjust_gradT_fraction, adjust_gradT_excess, &
         switch_to_no_mixing, switch_to_radiative, check_for_redo_MLT

      logical, parameter :: dbg = .false.
      integer, parameter :: kdbg = -1
      
      integer, parameter :: nvbs = num_mlt_partials

      contains


      subroutine set_mlt_vars(s, nzlo, nzhi, ierr)
         use star_utils, only: start_time, update_time
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, op_err, time0
         real(dp) :: total, opacity, gamma1, Cv, chiRho, chiT, Cp, &
            grada, P, xh, gradL_composition_term
         include 'formats'
         ierr = 0
         if (dbg) write(*, *) 'doing set_mlt_vars'
         gradL_composition_term = -1d0
         opacity = -1d0
         chiRho = -1d0
         chiT = -1d0
         Cp = -1d0
         grada = -1d0
         P = -1d0
         xh = -1d0 
         if (s% doing_timing) call start_time(s, time0, total)
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(guided)
         do k = nzlo, nzhi
            op_err = 0
            call do1_mlt(s, k, s% alpha_mlt(k), gradL_composition_term, &
               opacity, chiRho, chiT, Cp, grada, P, xh, &
               op_err)
            if (op_err /= 0) then
               ierr = op_err
               if (s% report_ierr) write(*,2) 'set_mlt_vars failed', k
            end if
            if (s% make_gradr_sticky_in_newton_iters .and. s% newton_iter > 3) then
               if (.not. s% fixed_gradr_for_rest_of_newton_iters(k)) &
                  s% fixed_gradr_for_rest_of_newton_iters(k) = &
                     (s% mlt_mixing_type(k) == no_mixing)
            end if            
         end do
!$OMP END PARALLEL DO
         if (s% doing_timing) call update_time(s, time0, total, s% time_mlt)
      end subroutine set_mlt_vars


      subroutine do1_mlt(s, k, mixing_length_alpha, gradL_composition_term_in, &
            opacity_face_in, chiRho_face_in, &
            chiT_face_in, Cp_face_in, grada_face_in, P_face_in, xh_face_in, &
            ierr)
         ! get convection info for point k
         !use mlt_lib
         use eos_def
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha, &
            gradL_composition_term_in, opacity_face_in, &
            chiRho_face_in, chiT_face_in, &
            Cp_face_in, grada_face_in, P_face_in, xh_face_in
         integer, intent(out) :: ierr

         real(dp) :: m, mstar, L, r, dlnm_dlnq, v0, thc, asc, Q_face, &
            a, b, Pgas_div_P_limit, da_dlnd, da_dlnT, db_dlnd, db_dlnT, &
            max_q_for_Pgas_div_P_limit, min_q_for_Pgas_div_P_limit, &
            mlt_basics(num_mlt_results), prev_conv_vel, max_conv_vel, dt, &
            alfa, beta, d_alfa_dq00, d_alfa_dqm1, d_alfa_dqp1, &
            T_face, rho_face, P_face, Cv_face, gamma1_face, &
            chiRho_face, chiT_face, Cp_face, opacity_face, grada_face, v, &
            cv_old, gradr_factor, f, xh_face, tau_face, &
            d_grada_face_dlnd00, d_grada_face_dlnT00, &
            d_grada_face_dlndm1, d_grada_face_dlnTm1, &
            gradL_composition_term, dlnT, dlnP, &
            abs_du_div_cs, cs, porosity_F, &
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
            normal_mlt_gradT_factor
         real(dp), target :: mlt_partials1_ary(num_mlt_partials*num_mlt_results)
         real(dp), pointer :: mlt_partials1(:), mlt_partials(:,:), vel(:)
         integer :: i, mixing_type, h1, nz, k_T_max
         real(dp), parameter :: conv_vel_mach_limit = 0.9d0
         logical :: Schwarzschild_stable, Ledoux_stable
         character (len=32) :: MLT_option
         logical, parameter :: just_gradr = .false.

         include 'formats'

         ierr = 0
         if (s% gamma_law_hydro > 0d0) then  
            call set_no_MLT_results
            return
         end if

         nz = s% nz
         
         if (k < 1 .or. k > nz) then
            write(*,3) 'bad k for do1_mlt', k, nz
            ierr = -1
            return
            stop 'mlt'
         end if
         
         MLT_option = s% MLT_option

         mlt_partials1 => mlt_partials1_ary
         mlt_partials(1:num_mlt_partials,1:num_mlt_results) => &
            mlt_partials1(1:num_mlt_partials*num_mlt_results)

         ! TODO: for mass_flag might need some adjustement when using mass_corrections
         m = s% m_grav(k)
         mstar = s% m_grav(1)

         if (m < 0) then
            write(*,2) 'mlt nz', s% nz
            write(*,2) 's% q(k)', k, s% q(k)
            write(*,2) 's% m(k)', k, s% m(k)
            write(*,2) 's% m_grav(k)', k, s% m_grav(k)
         end if

         dlnm_dlnq = 1
         r = s% r(k)
         L = s% L(k)

         if (is_bad_num(L)) then
            write(*,2) 'do1_mlt L', k, L
            call mesa_error(__FILE__,__LINE__)
         end if

         if (s% rotation_flag .and. s% mlt_use_rotation_correction) then
            gradr_factor = s% ft_rot(k)/s% fp_rot(k)
         else
            gradr_factor = 1d0
         end if
         if (is_bad_num(gradr_factor)) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) 'do1_mlt_eval gradr_factor', k, gradr_factor
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'gradr_factor', k, gradr_factor
               stop 'do1_mlt_eval'
            end if
            return
         end if

         ! alfa is the fraction coming from k; (1-alfa) from k-1.
         if (k == 1) then
            alfa = 1d0
            d_alfa_dq00 = 0d0
            d_alfa_dqm1 = 0d0
            d_alfa_dqp1 = 0d0
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            if (s% mass_flag) then
               if (k < s% nz) then
                  ! remember that dq(k) = q(k) - q(k+1)
                  ! so here alfa = (q(k-1) - q(k))/(q(k-1) - q(k+1))
                  d_alfa_dq00 = -1d0/(s% q(k-1) - s% q(k+1))
                  d_alfa_dqm1 = 1d0/(s% q(k-1) - s% q(k+1)) - (s% q(k-1) - s% q(k))/pow2(s% q(k-1) - s% q(k+1))
                  d_alfa_dqp1 = (s% q(k-1) - s% q(k))/pow2(s% q(k-1) - s% q(k+1))
               else
                  ! remember that dq(nz) = q(nz)
                  ! so here alfa = (q(k-1) - q(k))/q(k-1)
                  d_alfa_dq00 = -1d0/s% q(k-1)
                  d_alfa_dqm1 = 1d0/s% q(k-1) - (s% q(k-1) - s% q(k))/pow2(s% q(k-1))
                  d_alfa_dqp1 = 0d0
               end if
            else
               d_alfa_dq00 = 0d0
               d_alfa_dqm1 = 0d0
               d_alfa_dqp1 = 0d0
            end if
         end if
         beta = 1d0 - alfa
         h1 = s% net_iso(ih1)
         
         if (k > 0) then
            d_chiRho_00_dlnd = s% d_eos_dlnd(i_chiRho, k)
            d_chiRho_00_dlnT = s% d_eos_dlnT(i_chiRho, k)
            d_chiT_00_dlnd = s% d_eos_dlnd(i_chiT, k)
            d_chiT_00_dlnT = s% d_eos_dlnT(i_chiT, k)
            d_Cp_00_dlnd = s% d_eos_dlnd(i_Cp, k)
            d_Cp_00_dlnT = s% d_eos_dlnT(i_Cp, k)
            d_grada_00_dlnd = s% d_eos_dlnd(i_grad_ad, k)
            d_grada_00_dlnT = s% d_eos_dlnT(i_grad_ad, k)
            d_opacity_00_dlnd = s% d_opacity_dlnd(k)
            d_opacity_00_dlnT = s% d_opacity_dlnT(k)
         else
            d_chiRho_00_dlnd = 0d0
            d_chiRho_00_dlnT = 0d0
            d_chiT_00_dlnd = 0d0
            d_chiT_00_dlnT = 0d0
            d_Cp_00_dlnd = 0d0
            d_Cp_00_dlnT = 0d0
            d_grada_00_dlnd = 0d0
            d_grada_00_dlnT = 0d0
            d_opacity_00_dlnd = 0d0
            d_opacity_00_dlnT = 0d0
         end if
         
         if (k > 1) then
            d_chiRho_m1_dlnd = s% d_eos_dlnd(i_chiRho, k-1)
            d_chiRho_m1_dlnT = s% d_eos_dlnT(i_chiRho, k-1)
            d_chiT_m1_dlnd = s% d_eos_dlnd(i_chiT, k-1)
            d_chiT_m1_dlnT = s% d_eos_dlnT(i_chiT, k-1)
            d_Cp_m1_dlnd = s% d_eos_dlnd(i_Cp, k-1)
            d_Cp_m1_dlnT = s% d_eos_dlnT(i_Cp, k-1)
            d_grada_m1_dlnd = s% d_eos_dlnd(i_grad_ad, k-1)
            d_grada_m1_dlnT = s% d_eos_dlnT(i_grad_ad, k-1)
            d_opacity_m1_dlnd = s% d_opacity_dlnd(k-1)
            d_opacity_m1_dlnT = s% d_opacity_dlnT(k-1)
         else
            d_chiRho_m1_dlnd = 0d0
            d_chiRho_m1_dlnT = 0d0
            d_chiT_m1_dlnd = 0d0
            d_chiT_m1_dlnT = 0d0
            d_Cp_m1_dlnd = 0d0
            d_Cp_m1_dlnT = 0d0
            d_grada_m1_dlnd = 0d0
            d_grada_m1_dlnT = 0d0
            d_opacity_m1_dlnd = 0d0
            d_opacity_m1_dlnT = 0d0
         end if

         opacity_face = opacity_face_in
         chiRho_face = chiRho_face_in
         chiT_face = chiT_face_in
         Cp_face = Cp_face_in
         grada_face = grada_face_in
         P_face = P_face_in
         xh_face = xh_face_in
         gradL_composition_term = gradL_composition_term_in      
         rho_00 = s% rho(k)
         T_00 = s% T(k)
         P_00 = s% P(k)
         chiRho_00 = s% chiRho(k)
         chiT_00 = s% chiT(k)
         Cp_00 = s% Cp(k)
         opacity_00 = s% opacity(k)
         grada_00 = s% grada(k)

         if (alfa == 1d0) then
            
            rho_m1 = 0d0
            T_m1 = 0d0
            P_m1 = 0d0
            chiRho_m1 = 0d0
            chiT_m1 = 0d0
            Cp_m1 = 0d0
            opacity_m1 = 0d0
            grada_m1 = 0d0

            T_face = T_00
            rho_face = rho_00
            tau_face = s% tau_start(k)
            if (P_face < 0) P_face = P_00
            if (chiRho_face < 0) chiRho_face = chiRho_00
            if (chiT_face < 0) chiT_face = chiT_00
            if (Cp_face < 0) Cp_face = Cp_00
            if (opacity_face < 0) opacity_face = opacity_00
            if (grada_face < 0) grada_face = grada_00
            if (h1 /= 0 .and. xh_face < 0) xh_face = s% xa(h1, k)
            s% actual_gradT(k) = 0

         else
         
            rho_m1 = s% rho(k-1)
            T_m1 = s% T(k-1)
            P_m1 = s% P(k-1)
            chiRho_m1 = s% chiRho(k-1)
            chiT_m1 = s% chiT(k-1)
            Cp_m1 = s% Cp(k-1)
            opacity_m1 = s% opacity(k-1)
            grada_m1 = s% grada(k-1)

            tau_face = alfa*s% tau_start(k) + beta*s% tau_start(k-1)
            T_face = alfa*T_00 + beta*T_m1
            rho_face = alfa*rho_00 + beta*rho_m1
            if (P_face < 0) P_face = alfa*P_00 + beta*P_m1
            if (chiRho_face < 0) chiRho_face = alfa*chiRho_00 + beta*chiRho_m1
            if (chiT_face < 0) chiT_face = alfa*chiT_00 + beta*chiT_m1
            if (Cp_face < 0) Cp_face = alfa*Cp_00 + beta*Cp_m1
            if (opacity_face < 0) opacity_face = alfa*opacity_00 + beta*opacity_m1
            if (grada_face < 0) grada_face = alfa*grada_00 + beta*grada_m1
            if (h1 /= 0 .and. xh_face < 0) xh_face = alfa*s% xa(h1,k) + beta*s% xa(h1,k-1)
            dlnT = (T_m1 - T_00)/T_face
            dlnP = (P_m1 - P_00)/P_face
            if (abs(dlnP) > 1d-20) then
               s% actual_gradT(k) = dlnT/dlnP
            else
               s% actual_gradT(k) = 0
            end if
         end if

         s% grada_face(k) = alfa*grada_00 + beta*grada_m1
         !s% d_grada_face_dlnd00(k) = alfa*d_grada_00_dlnd
         !s% d_grada_face_dlnT00(k) = alfa*d_grada_00_dlnT
         !if (k == 1) then
         !   s% d_grada_face_dlndm1(k) = 0d0
         !   s% d_grada_face_dlnTm1(k) = 0d0
         !else
         !   s% d_grada_face_dlndm1(k) = beta*d_grada_m1_dlnd
         !   s% d_grada_face_dlnTm1(k) = beta*d_grada_m1_dlnT
         !end if

         if (Cp_face <= 0d0) then
            ierr = -1
            if (.not. s% report_ierr) return
            write(*,2) 'T_face', k, T_face
            write(*,2) 'rho_face', k, rho_face
            write(*,2) 'P_face', k, P_face
            write(*,2) 'chiRho_face', k, chiRho_face
            write(*,2) 'chiT_face', k, chiT_face
            write(*,2) 'opacity_face', k, opacity_face
            write(*,2) 'grada_face', k, grada_face
            write(*,2) 'tau_face', k, tau_face
            write(*,2) 'Cp_face', k, Cp_face
            write(*,2) 'Cp_00', k, Cp_00
            write(*,2) 'Cp_m1', k, Cp_m1
            write(*,2) 'alfa', k, alfa
            write(*,2) 'beta', k, beta
            write(*,*) 'MLT_option ', trim(MLT_option)
            stop 'mlt'
         end if
         
         tau_face = 1d0 ! disable this for now.

         if (s% use_porosity_with_dPrad_dm_form) then
            porosity_F = max(s% min_porosity_factor, s% porosity_factor(k))
            
            opacity_face = opacity_face/porosity_F
            opacity_00 = opacity_00/porosity_F
            opacity_m1 = opacity_m1/porosity_F
            d_opacity_00_dlnd = d_opacity_00_dlnd/porosity_F
            d_opacity_00_dlnT = d_opacity_00_dlnT/porosity_F
            d_opacity_m1_dlnd = d_opacity_m1_dlnd/porosity_F
            d_opacity_m1_dlnT = d_opacity_m1_dlnT/porosity_F
         else
            porosity_F = 1d0
         end if

         if (s% RTI_flag .and. tau_face > 0.667d0) then
            if (s% alpha_RTI(k) > 0d0) then ! RTI takes priority over MLT
               call set_no_mixing
               return
            end if
         end if
         
         if (s% lnT_start(k)/ln10 > s% max_logT_for_mlt) then
            call set_no_mixing
            return
         end if
         
         if (s% no_MLT_below_shock) then ! check for outward shock above k
            if (s% u_flag) then
               vel => s% u
            else
               vel => s% v
            end if
            do i=k-1,1,-1
               cs = s% csound(i)
               if (vel(i+1) >= cs .and. vel(i) < cs) then
                  call set_no_mixing
                  return
               end if
            end do
         end if
         
         if (s% no_MLT_below_T_max) then
            k_T_max = maxloc(s% T_start(1:nz),dim=1)
            if (k > k_T_max) then
               call set_no_mixing
               return
            end if
         end if
         
         if (s% make_gradr_sticky_in_newton_iters) then
            if (s% fixed_gradr_for_rest_of_newton_iters(k)) then
               call set_no_mixing
               return
            end if
         end if

         thc = s% thermohaline_coeff
         asc = s% alpha_semiconvection
         if (s% center_h1 > s% semiconvection_upper_limit_center_h1) asc = 0

         cv_old = -1
         if (s% have_previous_conv_vel .and. .not. s% conv_vel_flag) then
            if (s% generations >= 2) then
               if (.not. is_bad(s% conv_vel_old(k))) &
                  cv_old = s% conv_vel_old(k)
            else if (s% use_previous_conv_vel_from_file) then
               if (.not. is_bad(s% prev_conv_vel_from_file(k))) &
                  cv_old = s% prev_conv_vel_from_file(k)
            end if
         end if
         
         if (T_face >= s% min_T_for_acceleration_limited_conv_velocity &
               .and. T_face <= s% max_T_for_acceleration_limited_conv_velocity &
               .and. cv_old >= 0 .and. .not. s% conv_vel_flag) then
            prev_conv_vel = cv_old
            dt = s% dt
         else
            prev_conv_vel = -1
            dt = -1
         end if

         max_conv_vel = s% csound_face(k)*s% max_conv_vel_div_csound
         if (prev_conv_vel >= 0) then
            if (must_limit_conv_vel(s,k)) then
               max_conv_vel = prev_conv_vel
            end if
         else if ( &
               s% dt < s% min_dt_for_increases_in_convection_velocity) then
            prev_conv_vel = 0d0
            max_conv_vel = 1d-2*s% dt*s% cgrav(k)
         end if

         if (s% csound_start(k) > 0d0) then
            if (s% u_flag) then        
               abs_du_div_cs = 0d0
               if (s% u_start(k)/1e5 > s% max_v_for_convection) then
                  max_conv_vel = 0d0              
               else if (s% q(k) > s% max_q_for_convection_with_hydro_on) then
                  max_conv_vel = 0d0
               else if ((abs(s% u_start(k))) >= &
                     s% csound_start(k)*s% max_v_div_cs_for_convection) then
                  max_conv_vel = 0d0              
               else
                  if (k == 1) then
                     abs_du_div_cs = 1d99
                  else if (k < nz) then
                     abs_du_div_cs = max(abs(s% u_start(k) - s% u_start(k+1)), &
                         abs(s% u_start(k) - s% u_start(k-1))) / s% csound_start(k)
                  end if
                  if (abs_du_div_cs > s% max_abs_du_div_cs_for_convection) then
                     ! main purpose is to force radiative in shock face
                     max_conv_vel = 0d0
                  end if
               end if
            else if (s% v_flag) then
               if (s% v_start(k)/1e5 > s% max_v_for_convection) then
                  max_conv_vel = 0d0              
               else if (s% q(k) > s% max_q_for_convection_with_hydro_on) then
                  max_conv_vel = 0d0
               else if ((abs(s% v_start(k))) >= &
                     s% csound_start(k)*s% max_v_div_cs_for_convection) then
                  max_conv_vel = 0d0
               end if
            end if
            if (max_conv_vel == 0d0) then
               MLT_option = 'none'
            end if
         end if

         if (s% use_Ledoux_criterion .and. gradL_composition_term < 0) then
            gradL_composition_term = s% gradL_composition_term(k)
         else
            gradL_composition_term = 0d0
         end if

         if (.not. s% conv_vel_flag) then
            normal_mlt_gradT_factor = 1d0
         else if (abs(s% max_q_for_normal_mlt_gradT_full_on - &
                  s% min_q_for_normal_mlt_gradT_full_off) < 1d-10) then
            if (s% q(k) > s% min_q_for_normal_mlt_gradT_full_off) then
               normal_mlt_gradT_factor = 1d0
            else
               normal_mlt_gradT_factor = 0d0
            end if
         else
            normal_mlt_gradT_factor = &
               (s% q(k) - s% min_q_for_normal_mlt_gradT_full_off)/&
               (s% max_q_for_normal_mlt_gradT_full_on - s% min_q_for_normal_mlt_gradT_full_off)
            normal_mlt_gradT_factor = min(1d0, normal_mlt_gradT_factor)
            normal_mlt_gradT_factor = max(0d0, normal_mlt_gradT_factor)
         end if

         call do1_mlt_eval(s, k, &
            s% cgrav(k), m, mstar, r, L, xh_face, &            
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
            asc, s% semiconvection_option, thc, s% thermohaline_option, &
            s% dominant_iso_for_thermohaline(k), &
            mixing_length_alpha, s% alt_scale_height_flag, s% remove_small_D_limit, &
            MLT_option, s% Henyey_MLT_y_param, s% Henyey_MLT_nu_param, &
            s% gradT_smooth_low, s% gradT_smooth_high, s% smooth_gradT, &
            normal_mlt_gradT_factor, &
            prev_conv_vel, max_conv_vel, s% mlt_accel_g_theta, dt, tau_face, just_gradr, &
            mixing_type, mlt_basics, mlt_partials1, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) then
!$OMP critical
               write(*,*) 'ierr in do1_mlt_eval for k', k
               call show_stuff(.true.)
               stop 
!$OMP end critical
            end if
            return
         end if

         s% mlt_mixing_type(k) = mixing_type
         s% mlt_mixing_length(k) = mlt_basics(mlt_Lambda)
         s% mlt_Gamma(k) = mlt_basics(mlt_Gamma)
         s% mlt_D(k) = mlt_basics(mlt_D)
         s% mlt_D_semi(k) = mlt_basics(mlt_D_semi)
         s% mlt_D_thrm(k) = mlt_basics(mlt_D_thrm)
         s% gradr(k) = mlt_basics(mlt_gradr)
         s% scale_height(k) = mlt_basics(mlt_scale_height)
         s% gradL(k) = mlt_basics(mlt_gradL)
         s% conv_dP_term(k) = mlt_basics(mlt_conv_dP_term)
         s% mlt_cdc(k) = s% mlt_D(k)*pow2(pi4*r*r*rho_face)

         if ((s% scale_height(k) <= 0d0 .and. MLT_option /= 'none') &
                .or. is_bad_num(s% scale_height(k))) then
            ierr = -1
            if (.not. s% report_ierr) return
!$OMP critical
            write(*,2) 's% scale_height(k)', k, s% scale_height(k)
            call show_stuff(.true.)
            write(*,*)
            write(*,2) 's% scale_height(k)', k, s% scale_height(k)
            write(*,2) 's% grada(k)', k, s% grada(k)
            write(*,2) 's% gradr(k)', k, s% gradr(k)
            write(*,2) 's% gradT(k)', k, s% gradT(k)
            write(*,2) 's% gradL(k)', k, s% gradL(k)
            write(*,2) 'max_conv_vel', k, max_conv_vel
            stop 'mlt info'
!$OMP end critical
         end if

         call store_gradr_partials

         s% d_conv_dP_term_dlnR(k) = mlt_partials(mlt_dlnR, mlt_conv_dP_term)
         s% d_conv_dP_term_dL(k) = mlt_partials(mlt_dL, mlt_conv_dP_term)
         s% d_conv_dP_term_dlnd00(k) = mlt_partials(mlt_dlnd00, mlt_conv_dP_term)
         s% d_conv_dP_term_dlnT00(k) = mlt_partials(mlt_dlnT00, mlt_conv_dP_term)
         if (k == 1) then
            s% d_conv_dP_term_dlndm1(k) = 0d0
            s% d_conv_dP_term_dlnTm1(k) = 0d0
         else
            s% d_conv_dP_term_dlndm1(k) = mlt_partials(mlt_dlndm1, mlt_conv_dP_term)
            s% d_conv_dP_term_dlnTm1(k) = mlt_partials(mlt_dlnTm1, mlt_conv_dP_term)
         end if

         s% gradT(k) = mlt_basics(mlt_gradT)         
         s% d_gradT_dlnR(k) = mlt_partials(mlt_dlnR, mlt_gradT)
         s% d_gradT_dL(k) = mlt_partials(mlt_dL, mlt_gradT)
         s% d_gradT_dln_cvpv0(k) = mlt_partials(mlt_cv_var, mlt_gradT)
         s% d_gradT_dq00(k) = mlt_partials(mlt_dq00, mlt_gradT)
         s% d_gradT_dqm1(k) = mlt_partials(mlt_dqm1, mlt_gradT)
         s% d_gradT_dqp1(k) = mlt_partials(mlt_dqp1, mlt_gradT)
         if (is_bad(s% d_gradT_dln_cvpv0(k))) then
            write(*,2) 's% d_gradT_dln_cvpv0(k)', k, s% d_gradT_dln_cvpv0(k)
            stop 'mlt_info'
         end if
         s% d_gradT_dlnd00(k) = mlt_partials(mlt_dlnd00, mlt_gradT)  
         s% d_gradT_dlnT00(k) = mlt_partials(mlt_dlnT00, mlt_gradT) 
         if (k == 1) then
            s% d_gradT_dlndm1(k) = 0d0
            s% d_gradT_dlnTm1(k) = 0d0
         else
            s% d_gradT_dlndm1(k) = mlt_partials(mlt_dlndm1, mlt_gradT)
            s% d_gradT_dlnTm1(k) = mlt_partials(mlt_dlnTm1, mlt_gradT)
         end if
         if (s% mass_flag) then
            s% d_gradT_dq00(k) = mlt_partials(mlt_dq00, mlt_gradT)
            s% d_gradT_dqm1(k) = mlt_partials(mlt_dqm1, mlt_gradT)
            s% d_gradT_dqp1(k) = mlt_partials(mlt_dqp1, mlt_gradT)
         end if

         s% mlt_vc(k) = mlt_basics(mlt_convection_velocity)     
         s% d_mlt_vc_dlnR(k) = mlt_partials(mlt_dlnR, mlt_convection_velocity)
         s% d_mlt_vc_dL(k) = mlt_partials(mlt_dL, mlt_convection_velocity)
         s% d_mlt_vc_dlnd00(k) = mlt_partials(mlt_dlnd00, mlt_convection_velocity)  
         s% d_mlt_vc_dlnT00(k) = mlt_partials(mlt_dlnT00, mlt_convection_velocity) 
         if (k == 1) then
            s% d_mlt_vc_dlndm1(k) = 0d0
            s% d_mlt_vc_dlnTm1(k) = 0d0
         else
            s% d_mlt_vc_dlndm1(k) = mlt_partials(mlt_dlndm1, mlt_convection_velocity)
            s% d_mlt_vc_dlnTm1(k) = mlt_partials(mlt_dlnTm1, mlt_convection_velocity)
         end if

         if (mixing_type == 0 .and. s% mlt_vc(k) /= 0d0) then
            write(*,2) 'mixing_type mlt_vc', mixing_type, s% mlt_vc(k)
            stop 'MLT'
         end if

         Schwarzschild_stable = (s% gradr(k) < grada_face)
         Ledoux_stable = (s% gradr(k) < s% gradL(k))

         if (s% q(k) <= s% qmax_zero_non_radiative_luminosity .or. &
             s% m(k) <= s% m_force_radiative_core) then
            f = 0d0
         else if (s% mlt_gradT_fraction >= 0d0 .and. s% mlt_gradT_fraction <= 1d0) then
            f = s% mlt_gradT_fraction
         else
            f = s% adjust_mlt_gradT_fraction(k)
         end if
         call adjust_gradT_fraction(s, k, f)

         ! Note that L_conv can be negative when conv_vel /= 0 in a convectively
         ! stable region. This can happen when using conv_vel as a variable
         if (mlt_option /= 'none') then
            s% L_conv(k) = s% L(k) * (1d0 - s% gradT(k)/s% gradr(k)) ! C&G 14.109
         else
            s% L_conv(k) = 0d0
         end if
         
         if (k == 1 .or. mlt_option == 'none') then
            s% grad_superad(k) = 0d0
         else
            Q_face = chiT_face/(T_face*chiRho_face)
            s% grad_superad(k) = &
               4*pi*r*r*s% scale_height(k)*rho_face* &
                  (Q_face/Cp_face*(s% P(k-1)-s% P(k)) - (s% lnT(k-1)-s% lnT(k)))/s% dm_bar(k)
            ! grad_superad = area*Hp_face*rho_face*(Q_face/Cp_face*dP - dlogT)/dmbar
            ! Q_face = chiT_face/(T_face*chiRho_face)
            if (abs(s% lnP(k-1)-s% lnP(k)) < 1d-10) then
               s% grad_superad_actual(k) = 0
            else
               s% grad_superad_actual(k) = &
                  (s% lnT(k-1)-s% lnT(k))/(s% lnP(k-1)-s% lnP(k)) - grada_face
            end if
         end if

         if (is_bad_num(s% d_gradT_dlnT00(k))) then
            if (s% report_ierr) then
               write(*,2) 's% d_gradT_dlnT00(k)', k, s% d_gradT_dlnT00(k)
               return
            end if
            if (s% stop_for_bad_nums) then
!$OMP critical
               write(*,2) 's% d_gradT_dlnT00(k)', k, s% d_gradT_dlnT00(k)
               call show_stuff(.true.)
               stop 'mlt info'
!$OMP end critical
            end if
            ierr = -1
            return
         end if
         
         !if (k == 100) then
         !   call show_stuff(.true.)
         !   stop 'mlt'
         !end if
         
         !if (k == 738) write(*,4) 'gradT', k, s% newton_iter, mixing_type, s% gradT(k)
         
         !if (.false. .and. k == s% hydro_test_partials_k .and. &
         !      s% newton_iter == s% hydro_dump_iter_number) then
         !   call show_stuff(.true.)
         !   stop 'mlt'
         !end if

         contains
         
         subroutine get_mlt_eval_gradr_info(ierr)
            !use mlt_lib, only: mlt_eval_gradr_info
            integer, intent(out) :: ierr
            logical, parameter :: just_get_gradr = .true.
            include 'formats'

            call do1_mlt_eval(s, k, &
               s% cgrav(k), m, mstar, r, L, xh_face, &            
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
               asc, s% semiconvection_option, thc, s% thermohaline_option, &
               s% dominant_iso_for_thermohaline(k), &
               mixing_length_alpha, s% alt_scale_height_flag, s% remove_small_D_limit, &
               MLT_option, s% Henyey_MLT_y_param, s% Henyey_MLT_nu_param, &
               s% gradT_smooth_low, s% gradT_smooth_high, s% smooth_gradT, &
               normal_mlt_gradT_factor, &
               prev_conv_vel, max_conv_vel, s% mlt_accel_g_theta, dt, tau_face, just_get_gradr, &
               mixing_type, mlt_basics, mlt_partials1, ierr)
            if (ierr /= 0) return
            s% gradr(k) = mlt_basics(mlt_gradr)
            call store_gradr_partials
            if (is_bad(s% gradr(k))) then
               ierr = -1
               if (.not. s% report_ierr) return
!$OMP critical
               write(*,2) 's% gradr(k)', k, s% gradr(k)
               write(*,2) 'P_face', k, P_face
               write(*,2) 'opacity_face', k, opacity_face
               write(*,2) 'L', k, L
               write(*,2) 'm', k, m
               write(*,2) 's% cgrav(k)', k, s% cgrav(k)
               write(*,2) 'tau_face', k, tau_face
               write(*,2) 'T_face', k, T_face
               write(*,2) 'r', k, r
               write(*,2) 'rho_face', k, rho_face
               stop 'get_mlt_eval_gradr_info'
!$OMP end critical
            end if 
         end subroutine get_mlt_eval_gradr_info

         
         subroutine store_gradr_partials
            s% d_gradr_dlnR(k) = mlt_partials(mlt_dlnR, mlt_gradr)
            s% d_gradr_dL(k) = mlt_partials(mlt_dL, mlt_gradr)
            s% d_gradr_dlnd00(k) = mlt_partials(mlt_dlnd00, mlt_gradr)
            s% d_gradr_dlnT00(k) = mlt_partials(mlt_dlnT00, mlt_gradr)
            if (k == 1) then
               s% d_gradr_dlndm1(k) = 0d0
               s% d_gradr_dlnTm1(k) = 0d0
            else
               s% d_gradr_dlndm1(k) = mlt_partials(mlt_dlndm1, mlt_gradr)
               s% d_gradr_dlnTm1(k) = mlt_partials(mlt_dlnTm1, mlt_gradr)
            end if
         end subroutine store_gradr_partials


         subroutine show_stuff(with_results)
            logical, intent(in) :: with_results
            real(dp) :: vsem, Lambda, D, radiative_conductivity
            include 'formats'
            write(*,*)
            write(*,*) 'do1_mlt info for k, nz', k, s% nz
            write(*,2) 's% model_number', s% model_number
            write(*,2) 's% newton_iter', s% newton_iter
            write(*,*)
            write(*,1) 'cgrav =', s% cgrav(k)
            write(*,1) 'm =', m
            write(*,1) 'r =', r
            write(*,1) 'T =', T_face
            write(*,1) 'rho =', rho_face
            write(*,1) 'L =', L
            write(*,1) 'P =', P_face
            write(*,1) 'chiRho =', chiRho_face
            write(*,1) 'chiT =', chiT_face
            write(*,1) 'Cp =', Cp_face
            write(*,1) 'X =', xh_face
            write(*,1) 'opacity =', opacity_face
            write(*,1) 'grada =', grada_face
            write(*,*)
            write(*,1) 'gradr_factor =', gradr_factor
            write(*,1) 'gradL_composition_term =', s% gradL_composition_term(k)
            write(*,*)
            write(*,1) 'alpha_semiconvection =', asc
            write(*,'(a)') 'semiconvection_option = "' // trim(s% semiconvection_option) // '" '
            write(*,*)
            write(*,1) 'thermohaline_coeff =', thc
            write(*,'(a)') 'thermohaline_option = "' // trim(s% thermohaline_option) // '"'
            write(*,2) 'dominant_iso_for_thermohaline =', s% dominant_iso_for_thermohaline(k)
            write(*,*)
            write(*,1) 'mixing_length_alpha =', mixing_length_alpha
            if (s% alt_scale_height_flag) then
               write(*,'(a50)') '         alt_scale_height = .true.'
            else
               write(*,'(a50)') '         alt_scale_height = .false.'
            end if
            write(*,*)
            write(*,1) 'Henyey_y_param =', s% Henyey_MLT_y_param
            write(*,1) 'Henyey_nu_param =', s% Henyey_MLT_nu_param
            write(*,*)
            write(*,'(a)') "MLT_option = '" // trim(s% MLT_option) // "'"
            write(*,*)
            write(*,1) 'prev_conv_vel =', prev_conv_vel
            write(*,1) 'max_conv_vel =', min(1d99,max_conv_vel)
            write(*,*)
            write(*,1) 'gradT_smooth_low =', s% gradT_smooth_low
            write(*,1) 'gradT_smooth_high =', s% gradT_smooth_high
            write(*,*) 'smooth_gradT =', s% smooth_gradT
            write(*,*)
            write(*,1) 'dt =', dt
            write(*,1) 'tau =', tau_face
            write(*,*)
            write(*,*) '--------------------------------------'
            write(*,*)
            write(*,*)
            write(*,*)

            write(*,1) 'logRho =', s% lnd(k)/ln10
            write(*,1) 'logT =', s% lnT(k)/ln10
            write(*,1) 'x =', s% x(k)
            write(*,1) 'z =', 1d0 - (s% x(k) + s% y(k))
            write(*,1) 'abar =', s% abar(k)
            write(*,1) 'zbar =', s% zbar(k)
            write(*,*)
            write(*,*)
            write(*,3) 'k, nz', k, s% nz
            write(*,*)
            write(*,*)
            if (k > 1) then
               write(*,2) 's% opacity(k)', k, s% opacity(k)
               write(*,2) 's% opacity(k-1)', k-1, s% opacity(k-1)
               write(*,1) 'alfa', alfa
               write(*,1) 'beta', beta
               write(*,1) 'alfa', alfa*s% opacity(k)
               write(*,1) 'beta', beta*s% opacity(k-1)
               write(*,1) 'opacity_face', opacity_face
            end if

            if (ierr /= 0 .or. .not. with_results) return
            write(*,1) 's% gradr(k)', s% gradr(k)
            write(*,1) 's% gradT(k)', s% gradT(k)
            write(*,1) 's% gradL(k)', s% gradL(k)
            write(*,1) 's% gradL(k) - grada_face', s% gradL(k) - grada_face
            write(*,*)
            write(*,1) 's% mlt_D(k)', s% mlt_D(k)
            write(*,1) 's% mlt_D_semi(k)', s% mlt_D_semi(k)
            write(*,1) 's% mlt_D_thrm(k)', s% mlt_D_thrm(k)
            write(*,1) 's% mlt_vc(k)', s% mlt_vc(k)
            write(*,2) 's% mlt_mixing_type(k)', s% mlt_mixing_type(k)
            write(*,1) 's% mlt_mixing_length(k)', s% mlt_mixing_length(k)
            write(*,1) 's% d_gradT_dlnd00(k)', s% d_gradT_dlnd00(k)
            write(*,1) 's% d_gradT_dlnT00(k)', s% d_gradT_dlnT00(k)
            write(*,1) 's% d_gradT_dlndm1(k)', s% d_gradT_dlndm1(k)
            write(*,1) 's% d_gradT_dlnTm1(k)', s% d_gradT_dlnTm1(k)
            write(*,1) 's% d_gradT_dlnR(k)', s% d_gradT_dlnR(k)
            write(*,1) 's% d_gradT_dL(k)', s% d_gradT_dL(k)
            write(*,*)
            write(*,1) 's% d_gradr_dlnd00(k)', s% d_gradr_dlnd00(k)
            write(*,1) 's% d_gradr_dlnT00(k)', s% d_gradr_dlnT00(k)
            write(*,1) 's% d_gradr_dlndm1(k)', s% d_gradr_dlndm1(k)
            write(*,1) 's% d_gradr_dlnTm1(k)', s% d_gradr_dlnTm1(k)
            write(*,1) 's% d_gradr_dlnR(k)', s% d_gradr_dlnR(k)
            write(*,1) 's% d_gradr_dL(k)', s% d_gradr_dL(k)
            write(*,*)
            write(*,1) 's% conv_dP_term(k)', s% conv_dP_term(k)
            write(*,1) 's% d_conv_dP_term_dlnd00(k)', s% d_conv_dP_term_dlnd00(k)
            write(*,1) 's% d_conv_dP_term_dlnT00(k)', s% d_conv_dP_term_dlnT00(k)

            write(*,*) 'Schwarzschild_stable', Schwarzschild_stable
            write(*,*) 'Ledoux_stable', Ledoux_stable
            write(*,*)

         end subroutine show_stuff

         subroutine set_no_MLT_results()
            s% mlt_mixing_type(k) = no_mixing
            s% mlt_mixing_length(k) = 0d0
            s% mlt_vc(k) = 0d0
            s% mlt_Gamma(k) = 0d0
            s% grada_face(k) = 0d0
            if (s% zero_gravity) then
               s% scale_height(k) = 0
            else
               s% scale_height(k) = P_face*r*r/(s% cgrav(k)*m*rho_face)
            end if
            s% gradL(k) = 0d0
            s% L_conv(k) = 0d0
            s% conv_dP_term(k) = 0d0
            s% d_conv_dP_term_dlnR(k) = 0d0
            s% d_conv_dP_term_dL(k) = 0d0
            s% d_conv_dP_term_dlnd00(k) = 0d0
            s% d_conv_dP_term_dlnT00(k) = 0d0
            s% d_conv_dP_term_dlndm1(k) = 0d0
            s% d_conv_dP_term_dlnTm1(k) = 0d0
            s% gradT(k) = 0d0
            s% d_gradT_dlnR(k) = 0d0
            s% d_gradT_dL(k) = 0d0
            if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
            if (s% mass_flag) then
               s% d_gradT_dq00(k) = 0d0
               s% d_gradT_dqm1(k) = 0d0
               s% d_gradT_dqp1(k) = 0d0
            end if
            s% d_gradT_dlnd00(k) = 0d0
            s% d_gradT_dlnT00(k) = 0d0
            s% d_gradT_dlndm1(k) = 0d0
            s% d_gradT_dlnTm1(k) = 0d0
            s% mlt_D(k) = 0d0
            s% mlt_D_semi(k) = 0d0
            s% mlt_D_thrm(k) = 0d0
            s% mlt_cdc(k) = 0d0
            s% actual_gradT(k) = 0
            s% grad_superad(k) = 0
         end subroutine set_no_MLT_results

         subroutine set_no_mixing()
            call get_mlt_eval_gradr_info(ierr)
            if (ierr /= 0) return
            s% mlt_mixing_type(k) = no_mixing
            s% mlt_mixing_length(k) = 0d0
            s% mlt_vc(k) = 0d0
            s% mlt_Gamma(k) = 0d0
            s% grada_face(k) = 0d0
            s% scale_height(k) = P_face*r*r/(s% cgrav(k)*m*rho_face)
            s% gradL(k) = 0d0
            s% L_conv(k) = 0d0
            s% conv_dP_term(k) = 0d0
            s% d_conv_dP_term_dlnR(k) = 0d0
            s% d_conv_dP_term_dL(k) = 0d0
            s% d_conv_dP_term_dlnd00(k) = 0d0
            s% d_conv_dP_term_dlnT00(k) = 0d0
            s% d_conv_dP_term_dlndm1(k) = 0d0
            s% d_conv_dP_term_dlnTm1(k) = 0d0
            s% gradT(k) = s% gradr(k)
            s% d_gradT_dlnR(k) = s% d_gradr_dlnR(k)
            s% d_gradT_dL(k) = s% d_gradr_dL(k)
            if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
            if (s% mass_flag) then
               s% d_gradT_dq00(k) = 0d0
               s% d_gradT_dqm1(k) = 0d0
               s% d_gradT_dqp1(k) = 0d0
            end if
            s% d_gradT_dlnd00(k) = s% d_gradr_dlnd00(k)
            s% d_gradT_dlnT00(k) = s% d_gradr_dlnT00(k)
            s% d_gradT_dlndm1(k) = s% d_gradr_dlndm1(k)
            s% d_gradT_dlnTm1(k) = s% d_gradr_dlnTm1(k)
            s% mlt_D(k) = 0d0
            s% mlt_D_semi(k) = 0d0
            s% mlt_D_thrm(k) = 0d0
            s% mlt_cdc(k) = 0d0
            s% actual_gradT(k) = 0
            s% grad_superad(k) = 0
         end subroutine set_no_mixing

      end subroutine do1_mlt


      logical function must_limit_conv_vel(s,k0)
         type (star_info), pointer :: s
         integer, intent(in) :: k0
         real(dp) :: alfa, beta, one_m_f
         integer :: k, wdth
         include 'formats'
         must_limit_conv_vel = .false.
         if (s% q(k0) <= s% max_q_for_limit_conv_vel .or. &
             s% m(k0) <= s% max_mass_in_gm_for_limit_conv_vel .or. &
             s% r(k0) <= s% max_r_in_cm_for_limit_conv_vel) then
            must_limit_conv_vel = .true.
            return
         end if
         if (.not. s% v_flag) return
         wdth = s% width_for_limit_conv_vel
         if (wdth < 0) return
         do k = max(2,k0-wdth),min(s% nz,k0+wdth)
            if (s% csound(k) < s% v(k) .and. s% v(k) <= s% csound(k-1)) then
               must_limit_conv_vel = .true.
               return
            end if
         end do
      end function must_limit_conv_vel


      subroutine adjust_gradT_fraction(s,k,f)
         ! replace gradT by combo of grada_face and gradr
         ! then check excess
         use eos_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: f
         integer, intent(in) :: k

         real(dp) :: alfa, beta, one_m_f

         include 'formats'

         if (f >= 0 .and. f <= 1) then

            ! alfa is the fraction coming from k; (1-alfa) from k-1.
            if (k == 1) then
               alfa = 1
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            end if
            beta = 1 - alfa

            if (f == 0d0) then
               s% gradT(k) = s% gradr(k)
               s% d_gradT_dlnR(k) = s% d_gradr_dlnR(k)
               s% d_gradT_dL(k) = s% d_gradr_dL(k)
               s% d_gradT_dlnd00(k) = s% d_gradr_dlnd00(k)
               s% d_gradT_dlnT00(k) = s% d_gradr_dlnT00(k)
               if (k > 1) then
                  s% d_gradT_dlndm1(k) = s% d_gradr_dlndm1(k)
                  s% d_gradT_dlnTm1(k) = s% d_gradr_dlnTm1(k)
               end if
            else ! mix
               one_m_f = 1d0 - f
               s% gradT(k) = f*s% grada_face(k) + one_m_f*s% gradr(k)
               s% d_gradT_dlnR(k) = one_m_f*s% d_gradr_dlnR(k)
               s% d_gradT_dL(k) = one_m_f*s% d_gradr_dL(k)
               s% d_gradT_dlnd00(k) = &
                  f*alfa*s% d_eos_dlnd(i_grad_ad, k) + one_m_f*s% d_gradr_dlnd00(k)
               s% d_gradT_dlnT00(k) = &
                  f*alfa*s% d_eos_dlnT(i_grad_ad, k) + one_m_f*s% d_gradr_dlnT00(k)
               if (k > 1) then
                  s% d_gradT_dlndm1(k) = &
                     f*beta*s% d_eos_dlnd(i_grad_ad, k-1) + one_m_f*s% d_gradr_dlndm1(k)
                  s% d_gradT_dlnTm1(k) = &
                     f*beta*s% d_eos_dlnT(i_grad_ad, k-1) + one_m_f*s% d_gradr_dlnTm1(k)
               end if
            end if
            if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
            ! TODO: must check q derivatives when using mass_flag

         end if

         s% gradT_sub_grada(k) = s% gradT(k) - s% grada_face(k) ! gradT_excess

         call adjust_gradT_excess(s, k)

      end subroutine adjust_gradT_fraction


      subroutine adjust_gradT_excess(s, k)
         use eos_def
         type (star_info), pointer :: s
         integer, intent(in) :: k

         real(dp) :: alfa, beta, log_tau, gradT_excess_alpha, &
            d_grada_face_dlnd00, d_grada_face_dlnT00, &
            d_grada_face_dlndm1, d_grada_face_dlnTm1

         include 'formats'

         !s% gradT_excess_alpha is calculated at start of step and held constant during iterations
         ! gradT_excess_alpha = 0 means no efficiency boost; = 1 means full efficiency boost

         gradT_excess_alpha = s% gradT_excess_alpha
         s% gradT_excess_effect(k) = 0d0

         if (gradT_excess_alpha <= 0 .or. &
             s% gradT_sub_grada(k) <= s% gradT_excess_f1) return

         if (s% lnT(k)/ln10 > s% gradT_excess_max_logT) return

         log_tau = log10_cr(s% tau(k))
         if (log_tau < s% gradT_excess_max_log_tau_full_off) return

         ! boost efficiency of energy transport

         ! alfa is the fraction coming from k; (1-alfa) from k-1.
         if (k == 1) then
            alfa = 1
         else
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         end if
         beta = 1 - alfa

         ! grada_face = alfa*s% grada(k) + beta*s% grada(k-1)
         d_grada_face_dlnd00 = alfa*s% d_eos_dlnd(i_grad_ad, k)
         d_grada_face_dlnT00 = alfa*s% d_eos_dlnT(i_grad_ad, k)
         if (k > 1) then
            d_grada_face_dlndm1 = beta*s% d_eos_dlnd(i_grad_ad, k-1)
            d_grada_face_dlnTm1 = beta*s% d_eos_dlnT(i_grad_ad, k-1)
         else
            d_grada_face_dlndm1 = 0
            d_grada_face_dlnTm1 = 0
         end if

         if (log_tau < s% gradT_excess_min_log_tau_full_on) &
            gradT_excess_alpha = gradT_excess_alpha* &
               (log_tau - s% gradT_excess_max_log_tau_full_off)/ &
               (s% gradT_excess_min_log_tau_full_on - s% gradT_excess_max_log_tau_full_off)

         alfa = s% gradT_excess_f2 ! for full boost, use this fraction of gradT
         if (gradT_excess_alpha < 1) then ! only partial boost, so increase alfa
            ! alfa goes to 1 as gradT_excess_alpha goes to 0
            ! alfa unchanged as gradT_excess_alpha goes to 1
            alfa = alfa + (1d0 - alfa)*(1d0 - gradT_excess_alpha)
         end if
         beta = 1d0 - alfa

         s% gradT(k) = alfa*s% gradT(k) + beta*s% grada_face(k)
         s% gradT_excess_effect(k) = beta

         s% d_gradT_dlnR(k) = alfa*s% d_gradT_dlnR(k)
         s% d_gradT_dL(k) = alfa*s% d_gradT_dL(k)
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = alfa*s% d_gradT_dln_cvpv0(k)
         !TODO: much check q derivatives when using mass_flag

         s% d_gradT_dlnd00(k) = &
            alfa*s% d_gradT_dlnd00(k) + beta*d_grada_face_dlnd00
         s% d_gradT_dlnT00(k) = &
            alfa*s% d_gradT_dlnT00(k) + beta*d_grada_face_dlnT00

         if (k > 1) then
            s% d_gradT_dlndm1(k) = &
               alfa*s% d_gradT_dlndm1(k) + beta*d_grada_face_dlndm1
            s% d_gradT_dlnTm1(k) = &
               alfa*s% d_gradT_dlnTm1(k) + beta*d_grada_face_dlnTm1
         end if

      end subroutine adjust_gradT_excess


      subroutine set_grads(s, ierr)
         use chem_def, only: chem_isos
         use star_utils, only: smooth
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, nz, j, cid, max_cid
         real(dp) :: val, max_val, A, Z
         real(dp), pointer, dimension(:) :: dlnP, dlnd, dlnT

         include 'formats'

         ierr = 0
         nz = s% nz
         call do_alloc(ierr)
         if (ierr /= 0) return

         do k = 2, nz
            dlnP(k) = s% lnP(k-1) - s% lnP(k)
            dlnd(k) = s% lnd(k-1) - s% lnd(k)
            dlnT(k) = s% lnT(k-1) - s% lnT(k)
         end do
         dlnP(1) = dlnP(2)
         dlnd(1) = dlnd(2)
         dlnT(1) = dlnT(2)

         call smooth(dlnP,nz)
         call smooth(dlnd,nz)
         call smooth(dlnT,nz)

         s% grad_density(1) = 0
         s% grad_temperature(1) = 0
         do k = 2, nz
            if (dlnP(k) >= 0) then
               s% grad_density(k) = 0
               s% grad_temperature(k) = 0
            else
               s% grad_density(k) = dlnd(k)/dlnP(k)
               s% grad_temperature(k) = dlnT(k)/dlnP(k)
            end if
         end do

         call smooth(s% grad_density,nz)
         call smooth(s% grad_temperature,nz)

         do k=1,nz
            s% gradL_composition_term(k) = s% unsmoothed_brunt_B(k)
         end do
         call smooth_gradL_composition_term

         call dealloc

         do k=3,nz-2
            max_cid = 0
            max_val = -1d99
            do j=1,s% species
               cid = s% chem_id(j)
               A = dble(chem_isos% Z_plus_N(cid))
               Z = dble(chem_isos% Z(cid))
               val = (s% xa(j,k-2) + s% xa(j,k-1) - s% xa(j,k) - s% xa(j,k+1))*(1d0 + Z)/A
               if (val > max_val) then
                  max_val = val
                  max_cid = cid
               end if
            end do
            s% dominant_iso_for_thermohaline(k) = max_cid
         end do
         s% dominant_iso_for_thermohaline(1:2) = &
            s% dominant_iso_for_thermohaline(3)
         s% dominant_iso_for_thermohaline(nz-1:nz) = &
            s% dominant_iso_for_thermohaline(nz-2)


         contains

         subroutine smooth_gradL_composition_term
            use star_utils, only: weighed_smoothing, threshold_smoothing
            logical, parameter :: preserve_sign = .false.
            real(dp), pointer, dimension(:) :: work
            integer :: k
            include 'formats'
            ierr = 0
            work => dlnd
            if (s% num_cells_for_smooth_gradL_composition_term <= 0) return
            call threshold_smoothing( &
               s% gradL_composition_term, s% threshold_for_smooth_gradL_composition_term, s% nz, &
               s% num_cells_for_smooth_gradL_composition_term, preserve_sign, work)
         end subroutine smooth_gradL_composition_term

         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               dlnP, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnd, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dlnT, nz, nz_alloc_extra, 'mlt', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays

      end subroutine set_grads


      subroutine switch_to_no_mixing(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% mlt_mixing_type(k) = no_mixing
         s% mlt_mixing_length(k) = 0
         s% mlt_D(k) = 0
         s% mlt_D_semi(k) = 0
         s% mlt_D_thrm(k) = 0
         s% mlt_cdc(k) = 0d0
         s% mlt_vc(k) = 0
      end subroutine switch_to_no_mixing


      subroutine switch_to_radiative(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% gradT(k) = s% gradr(k)
         s% d_gradT_dlnd00(k) = s% d_gradr_dlnd00(k)
         s% d_gradT_dlnT00(k) = s% d_gradr_dlnT00(k)
         s% d_gradT_dlndm1(k) = s% d_gradr_dlndm1(k)
         s% d_gradT_dlnTm1(k) = s% d_gradr_dlnTm1(k)
         s% d_gradT_dlnR(k) = s% d_gradr_dlnR(k)
         s% d_gradT_dL(k) = s% d_gradr_dL(k)
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
         !TODO: must check mass derivatives when mass_flag = true
      end subroutine switch_to_radiative


      subroutine switch_to_adiabatic(s,k)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         integer, intent(in) :: k
         s% gradT(k) = s% grada(k)
         s% d_gradT_dlnd00(k) = s% d_eos_dlnd(i_grad_ad, k)
         s% d_gradT_dlnT00(k) = s% d_eos_dlnT(i_grad_ad, k)
         s% d_gradT_dlndm1(k) = 0
         s% d_gradT_dlnTm1(k) = 0
         s% d_gradT_dlnR(k) = 0
         s% d_gradT_dL(k) = 0
         if (s% conv_vel_flag) s% d_gradT_dln_cvpv0(k) = 0d0
         !TODO: must check mass derivatives when mass_flag = true
      end subroutine switch_to_adiabatic


      subroutine set_gradT_excess_alpha(s, ierr)
         use alloc
         use star_utils, only: get_Lrad_div_Ledd, after_C_burn
         use chem_def, only: ih1, ihe4

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp) :: beta, lambda, phi, tmp, alpha, alpha2
         real(dp) :: &
            beta_limit, &
            lambda1, &
            beta1, &
            lambda2, &
            beta2, &
            dlambda, &
            dbeta

         integer :: k, k_beta, k_lambda, nz, h1, he4

         logical, parameter :: dbg = .false.

         include 'formats'

         if (dbg) write(*,*) 'enter set_gradT_excess_alpha'

         ierr = 0
         if (.not. s% okay_to_reduce_gradT_excess) then
            s% gradT_excess_alpha = 0
            if (dbg) write(*,1) 'okay_to_reduce_gradT_excess'
            return
         end if
         nz = s% nz

         h1 = s% net_iso(ih1)
         if (h1 /= 0) then
            if (s% xa(h1,nz) > s% gradT_excess_max_center_h1) then
               s% gradT_excess_alpha = 0
               if (dbg) write(*,1) 'gradT_excess_max_center_h1'
               return
            end if
         end if

         he4 = s% net_iso(ihe4)
         if (he4 /= 0) then
            if (s% xa(he4,nz) < s% gradT_excess_min_center_he4) then
               s% gradT_excess_alpha = 0
               if (dbg) write(*,1) 'gradT_excess_min_center_he4'
               return
            end if
         end if

         beta = 1d0 ! beta = min over k of Pgas(k)/P(k)
         k_beta = 0
         do k=1,nz
            tmp = s% Pgas(k)/s% P(k)
            if (tmp < beta) then
               k_beta = k
               beta = tmp
            end if
         end do

         beta = beta*(1d0 + s% xa(1,nz))

         s% gradT_excess_min_beta = beta
         if (dbg) write(*,2) 'gradT_excess_min_beta', k_beta, beta

         lambda = 0d0 ! lambda = max over k of Lrad(k)/Ledd(k)
         do k=2,k_beta
            tmp = get_Lrad_div_Ledd(s,k)
            if (tmp > lambda) then
               k_lambda = k
               lambda = tmp
            end if
         end do
         lambda = min(1d0,lambda)
         s% gradT_excess_max_lambda = lambda
         if (dbg) write(*,2) 'gradT_excess_max_lambda', k_lambda, lambda

         lambda1 = s% gradT_excess_lambda1
         beta1 = s% gradT_excess_beta1
         lambda2 = s% gradT_excess_lambda2
         beta2 = s% gradT_excess_beta2
         dlambda = s% gradT_excess_dlambda
         dbeta = s% gradT_excess_dbeta

         if (dbg) then
            write(*,1) 'lambda1', lambda1
            write(*,1) 'lambda2', lambda2
            write(*,1) 'lambda', lambda
            write(*,*)
            write(*,1) 'beta1', beta1
            write(*,1) 'beta2', beta2
            write(*,1) 'beta', beta
            write(*,*)
         end if

         ! alpha is fraction of full boost to apply
         ! depends on location in (beta,lambda) plane

         if (lambda1 < 0) then
            alpha = 1
         else if (lambda >= lambda1) then
            if (beta <= beta1) then
               alpha = 1
            else if (beta < beta1 + dbeta) then
               alpha = (beta1 + dbeta - beta)/dbeta
            else ! beta >= beta1 + dbeta
               alpha = 0
            end if
         else if (lambda >= lambda2) then
            beta_limit = beta2 + &
               (lambda - lambda2)*(beta1 - beta2)/(lambda1 - lambda2)
            if (beta <= beta_limit) then
               alpha = 1
            else if (beta < beta_limit + dbeta) then
               alpha = (beta_limit + dbeta - beta)/dbeta
            else
               alpha = 0
            end if
         else if (lambda > lambda2 - dlambda) then
            if (beta <= beta2) then
               alpha = 1
            else if (beta < beta2 + dbeta) then
               alpha = (lambda - (lambda2 - dlambda))/dlambda
            else ! beta >= beta2 + dbeta
               alpha = 0
            end if
         else ! lambda <= lambda2 - dlambda
            alpha = 0
         end if

         if (s% generations > 1 .and. lambda1 >= 0) then ! time smoothing
            s% gradT_excess_alpha = &
               (1d0 - s% gradT_excess_age_fraction)*alpha + &
               s% gradT_excess_age_fraction*s% gradT_excess_alpha_old
            if (s% gradT_excess_max_change > 0d0) then
               if (s% gradT_excess_alpha > s% gradT_excess_alpha_old) then
                  s% gradT_excess_alpha = min(s% gradT_excess_alpha, s% gradT_excess_alpha_old + &
                     s% gradT_excess_max_change)
               else
                  s% gradT_excess_alpha = max(s% gradT_excess_alpha, s% gradT_excess_alpha_old - &
                     s% gradT_excess_max_change)
               end if
            end if
         else
            s% gradT_excess_alpha = alpha
         end if

         if (s% gradT_excess_alpha < 1d-4) s% gradT_excess_alpha = 0d0
         if (s% gradT_excess_alpha > 0.9999d0) s% gradT_excess_alpha = 1d0

         if (dbg) then
            write(*,1) 'gradT excess new', alpha
            write(*,1) 's% gradT_excess_alpha_old', s% gradT_excess_alpha_old
            write(*,1) 's% gradT_excess_alpha', s% gradT_excess_alpha
            write(*,*)
         end if

      end subroutine set_gradT_excess_alpha


      subroutine check_for_redo_MLT(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, k_bot
         real(dp) :: bot_Hp, bot_r, top_Hp, top_r, dr
         logical :: dbg

         include 'formats'

         ! check_for_redo_MLT assumes that nzlo = 1, nzhi = nz
         ! that is presently true; make sure that assumption doesn't change
         if (.not. ((nzlo.eq.1).and.(nzhi.eq.s%nz))) then
            write(*,*) 'nzlo != 1 or nzhi != nz'
            call mesa_error(__FILE__,__LINE__)
         endif

         ierr = 0
         dbg = .false.

         bot_Hp = 0; bot_r = 0; top_Hp = 0; top_r = 0; dr = 0

         in_convective_region = (s% mlt_mixing_type(nzhi) == convective_mixing)
         k_bot = nzhi
         bot_r = s% r(k_bot)
         bot_Hp = s% scale_height(k_bot)

         do k=nzhi-1, nzlo+1, -1
            if (in_convective_region) then
               if (s% mlt_mixing_type(k) /= convective_mixing) then
                  call end_of_convective_region
               end if
            else ! in non-convective region
               if (s% mlt_mixing_type(k) == convective_mixing) then ! start of a convective region
                  k_bot = k+1
                  in_convective_region = .true.
                  bot_r = s% r(k_bot)
                  bot_Hp = s% scale_height(k_bot)
               end if
            end if
         end do

         if (in_convective_region) then
            k = 1 ! end at top
            call end_of_convective_region
         end if


         contains


         subroutine end_of_convective_region()
            integer :: kk, op_err, mix_type
            real(dp) :: Hp
            logical :: end_dbg

            9 format(a40, 3i7, 99(1pd26.16))
            include 'formats'

            in_convective_region = .false.

            end_dbg = .false.

            top_r = s% r(k)
            top_Hp = s% scale_height(k)
            dr = top_r - bot_r
            Hp = (bot_Hp + top_Hp)/2

            if (dr < s% alpha_mlt(k)*min(top_Hp, bot_Hp) .and. &
                  s% redo_conv_for_dr_lt_mixing_length) then
!$OMP PARALLEL DO PRIVATE(kk,op_err)
               do kk = k, k_bot
                  op_err = 0
                  call redo1_mlt(s,kk,dr,op_err)
                  if (op_err /= 0) ierr = op_err
               end do
!$OMP END PARALLEL DO
            else if (s% limit_mixing_length_by_dist_to_bdy > 0) then
!$OMP PARALLEL DO PRIVATE(kk,op_err)
               do kk = k, k_bot
                  op_err = 0
                  call redo1_mlt(s, kk, s% limit_mixing_length_by_dist_to_bdy* &
                     min(top_r - s% r(kk), s% r(kk) - bot_r), op_err)
                  if (op_err /= 0) ierr = op_err
               end do
!$OMP END PARALLEL DO
            end if


         end subroutine end_of_convective_region


         subroutine redo1_mlt(s,k,dr,ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: k
            real(dp), intent(in) :: dr
            integer, intent(out) :: ierr
            real(dp) :: Hp, reduced_alpha
            include 'formats'
            ierr = 0
            if (dr >= s% mlt_mixing_length(k)) return
            ! if convection zone is smaller than mixing length
            ! redo MLT with reduced alpha so mixing_length = dr
            Hp = s% scale_height(k)
            reduced_alpha = dr/Hp
            call do1_mlt(s, k, reduced_alpha, -1d0, &
               -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, -1d0, &
               ierr)
         end subroutine redo1_mlt


      end subroutine check_for_redo_MLT


      subroutine do1_mlt_eval(s, k, &
            cgrav, m, mstar, r, L, X, &            
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

         type (star_info), pointer :: s
         integer, intent(in) :: k
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
            gradT_smooth_low, gradT_smooth_high, &
            normal_mlt_gradT_factor
         logical, intent(in) :: alt_scale_height
         character (len=*), intent(in) :: thermohaline_option, MLT_option, semiconvection_option
         integer, intent(in) :: dominant_iso_for_thermohaline
         logical, intent(in) :: just_gradr, smooth_gradT
         integer, intent(out) :: mixing_type
         real(dp), intent(inout) :: mlt_basics(:) ! (num_mlt_results)
         real(dp), intent(inout), pointer :: mlt_partials1(:) ! =(num_mlt_partials, num_mlt_results)
         integer, intent(out) :: ierr
         
         real(dp), pointer :: mlt_partials(:,:)
         include 'formats'

         if (s% use_other_mlt) then
            call s% other_mlt(  &
               s% id, k, cgrav, m, mstar, r, L, X, &            
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
            return
         end if
         
         mlt_partials(1:num_mlt_partials,1:num_mlt_results) => &
            mlt_partials1(1:num_mlt_partials*num_mlt_results)
         
         ierr = 0
         mlt_basics = 0
         mlt_partials = 0
         call Get_results(s, k, &
            cgrav, m, mstar, r, L, X, &            
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
            mixing_type, &
            mlt_basics(mlt_gradT), mlt_partials(:,mlt_gradT), &
            mlt_basics(mlt_gradr), mlt_partials(:,mlt_gradr), &
            mlt_basics(mlt_gradL), mlt_partials(:,mlt_gradL), &
            mlt_basics(mlt_scale_height), mlt_partials(:,mlt_scale_height), &
            mlt_basics(mlt_Lambda), mlt_partials(:,mlt_Lambda), &
            mlt_basics(mlt_convection_velocity), mlt_partials(:,mlt_convection_velocity), &
            mlt_basics(mlt_D), mlt_partials(:,mlt_D), &
            mlt_basics(mlt_D_semi), mlt_partials(:,mlt_D_semi), &
            mlt_basics(mlt_D_thrm), mlt_partials(:,mlt_D_thrm), &
            mlt_basics(mlt_gamma), mlt_partials(:,mlt_gamma), &
            mlt_basics(mlt_conv_dP_term), mlt_partials(:,mlt_conv_dP_term), &
            mlt_basics(mlt_unsmoothed_gradT), mlt_partials(:,mlt_unsmoothed_gradT), &
            ierr)
         if (mlt_basics(mlt_D) < 0) then
            !stop 'do1_mlt_eval: mlt_basics(mlt_D) < 0'
         end if
         
         if (is_bad(mlt_partials(mlt_cv_var,mlt_gradT))) then
            write(*,2) 'mlt_partials(mlt_cv_var,mlt_gradT)', k, mlt_partials(mlt_cv_var,mlt_gradT)
            ierr = -1
            !stop 'do1_mlt_eval'
         end if
         
      end subroutine do1_mlt_eval
      

      subroutine Get_results(ss, kz, &
            cgrav, m, mstar, r, L, xh, &            
            T, rho, P, chiRho, chiT, Cp, opacity, grada, &            
            a_00, a_m1, d_a_00_dq00, d_a_00_dqm1, d_a_00_dqp1, &
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
            prev_conv_vel, max_conv_vel, g_theta, dt, tau, just_gradr, mixing_type, &
            gradT, d_gradT_dvb, &
            gradr, d_gradr_dvb, &
            gradL, d_gradL_dvb, &
            scale_height, d_scale_height_dvb, &
            Lambda, d_Lambda_dvb, &
            conv_vel, d_conv_vel_dvb, & ! convection velocity
            D, d_D_dvb, &
            D_semi, d_D_semi_dvb, &
            D_thrm, d_D_thrm_dvb, &
            Gamma, d_Gamma_dvb, &
            conv_P, d_conv_P_dvb, &
            unsmoothed_gradT, d_unsmoothed_gradT_dvb, &
            ierr)

         use utils_lib, only: is_bad
         use chem_def, only: chem_isos
         
         type (star_info), pointer :: ss
         integer, intent(in) :: kz
         real(dp), intent(in) :: &
            cgrav, m, mstar, r, L, xh, &            
            T, rho, P, chiRho, chiT, Cp, opacity, grada, &            
            a_00, a_m1, d_a_00_dq00, d_a_00_dqm1, d_a_00_dqp1, &
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
         logical, intent(in) :: smooth_gradT, just_gradr
         
         integer, intent(out) :: mixing_type
         real(dp), intent(inout) :: gradT, d_gradT_dvb(:)
         real(dp), intent(inout) :: gradr, d_gradr_dvb(:)
         real(dp), intent(inout) :: gradL, d_gradL_dvb(:)
         real(dp), intent(inout) :: scale_height, d_scale_height_dvb(:)
         real(dp), intent(inout) :: Lambda, d_Lambda_dvb(:)
         real(dp), intent(inout) :: conv_vel, d_conv_vel_dvb(:)
         real(dp), intent(inout) :: D, d_D_dvb(:), D_semi, d_D_semi_dvb(:), D_thrm, d_D_thrm_dvb(:)
         real(dp), intent(inout) :: Gamma, d_Gamma_dvb(:) ! convective efficiency
         real(dp), intent(inout) :: conv_P, d_conv_P_dvb(:)
         real(dp), intent(inout) :: unsmoothed_gradT, d_unsmoothed_gradT_dvb(:)
         
         integer, intent(out) :: ierr

         real(dp) :: scale_height1, scale_height2
         real(dp) :: Pg, Pr, dP_dvb(nvbs), dPg_dvb(nvbs), dPr_dvb(nvbs), dRho_dvb(nvbs)
         real(dp) :: dT_dvb(nvbs), dL_dvb(nvbs), alpha, phi, dgrad, denom, tmp
         
         real(dp) :: grav, d_grav_dvb(nvbs)
         real(dp) :: diff_grads, d_diff_grads_dvb(nvbs)
         real(dp) :: convective_conductivity, d_cc_dvb(nvbs)
         real(dp) :: radiative_conductivity, d_rc_dvb(nvbs)
         real(dp) :: surf, dsurf_dvb(nvbs)
         real(dp) :: beta, d_beta_dvb(nvbs)
         real(dp) :: chi, d_chi_dvb(nvbs)
         real(dp) :: D_div_B, d_D_div_B_dvb(nvbs)
         real(dp) :: Q, dQ_dvb(nvbs)
         real(dp) :: A, dA_dvb(nvbs)
         real(dp) :: Bcubed, d_Bcubed_dvb(nvbs)
         real(dp) :: Zeta, d_Zeta_dvb(nvbs)
         real(dp) :: d_Cp_dvb(nvbs)
         real(dp) :: dR_dvb(nvbs)
         real(dp) :: d_opacity_dvb(nvbs)
         real(dp) :: d_grada_dvb(nvbs)
         real(dp) :: Dconv, d_Dconv_dvb(nvbs)
         real(dp) :: delta, d_delta_dvb(nvbs)
         real(dp) :: f, f0, d_f0_dvb(nvbs)
         real(dp) :: f1, d_f1_dvb(nvbs)
         real(dp) :: f2, d_f2_dvb(nvbs)
         real(dp) :: x, d_x_dvb(nvbs)

         real(dp) :: d_chiT_dvb(nvbs), d_chiRho_dvb(nvbs)

         real(dp) :: a0, omega, theta, s_gradr, f_gradr, dilute_factor
         real(dp) :: d_omega_dvb(nvbs), d_a0_dvb(nvbs), d_theta_dvb(nvbs)
         
         integer :: i
         real(dp), parameter :: tiny = 1d-30, min_D_th = 1d-3
         character (len=256) :: message        
         logical ::  quit
         real(dp) :: diff_grad, K, gamma0, L_ratio, frac, s, &
            dilution_factor, conv_tau, init_conv_vel
         real(dp) :: K_T, K_mu, nu_rad, nu_mol, nu, grad_mu, R0, r_th, H_P

         real(dp) :: scale_factor, interp_factor, dinterp_factor
         real(dp) :: gradT_temp, d_gradT_temp_dvb(nvbs)
         
         logical :: test_partials, debug

         include 'formats.dek'

         !test_partials = (kz == ss% hydro_test_partials_k)
         test_partials = .false.
      
         ierr = 0
         gradT_temp = 0
         debug = .false.
         
         dL_dvb(:) = 0d0
         dL_dvb(mlt_dL) = 1d0
         
         dR_dvb(:) = 0d0
         dR_dvb(mlt_dlnR) = r
         
         call set1_dvb(dT_dvb, &
            0d0, a_00*T_00, 0d0, a_m1*T_m1, &
            d_a_00_dq00*(T_00-T_m1), &
            d_a_00_dqm1*(T_00-T_m1), &
            d_a_00_dqp1*(T_00-T_m1))
            
         call set1_dvb(dRho_dvb, &
            a_00*rho_00, 0d0, a_m1*rho_m1, 0d0, &
            d_a_00_dq00*(rho_00-rho_m1), &
            d_a_00_dqm1*(rho_00-rho_m1), &
            d_a_00_dqp1*(rho_00-rho_m1))
            
         call set1_dvb(dP_dvb, &
            a_00*P_00*chiRho_00, a_00*P_00*chiT_00, &
            a_m1*P_m1*chiRho_m1, a_m1*P_m1*chiT_m1, &
            d_a_00_dq00*(P_00-P_m1), &
            d_a_00_dqm1*(P_00-P_m1), &
            d_a_00_dqp1*(P_00-P_m1))
            
         call set1_dvb(d_chiT_dvb, &
            a_00*d_chiT_00_dlnd, a_00*d_chiT_00_dlnT, &
            a_m1*d_chiT_m1_dlnd, a_m1*d_chiT_m1_dlnT, &
            d_a_00_dq00*(chiT_00-chiT_m1), &
            d_a_00_dqm1*(chiT_00-chiT_m1), &
            d_a_00_dqp1*(chiT_00-chiT_m1))
            
         call set1_dvb(d_chiRho_dvb, &
            a_00*d_chiRho_00_dlnd, a_00*d_chiRho_00_dlnT, &
            a_m1*d_chiRho_m1_dlnd, a_m1*d_chiRho_m1_dlnT, &
            d_a_00_dq00*(chiRho_00-chiRho_m1), &
            d_a_00_dqm1*(chiRho_00-chiRho_m1), &
            d_a_00_dqp1*(chiRho_00-chiRho_m1))
            
         call set1_dvb(d_Cp_dvb, &
            a_00*d_Cp_00_dlnd, a_00*d_Cp_00_dlnT, &
            a_m1*d_Cp_m1_dlnd, a_m1*d_Cp_m1_dlnT, &
            d_a_00_dq00*(Cp_00-Cp_m1), &
            d_a_00_dqm1*(Cp_00-Cp_m1), &
            d_a_00_dqp1*(Cp_00-Cp_m1))
            
         call set1_dvb(d_opacity_dvb, &
            a_00*d_opacity_00_dlnd, a_00*d_opacity_00_dlnT, &
            a_m1*d_opacity_m1_dlnd, a_m1*d_opacity_m1_dlnT, &
            d_a_00_dq00*(opacity_00-opacity_m1), &
            d_a_00_dqm1*(opacity_00-opacity_m1), &
            d_a_00_dqp1*(opacity_00-opacity_m1))
            
         call set1_dvb(d_grada_dvb, &
            a_00*d_grada_00_dlnd, a_00*d_grada_00_dlnT, &
            a_m1*d_grada_m1_dlnd, a_m1*d_grada_m1_dlnT, &
            d_a_00_dq00*(grada_00-grada_m1), &
            d_a_00_dqm1*(grada_00-grada_m1), &
            d_a_00_dqp1*(grada_00-grada_m1))

         if (is_bad(d_grada_dvb(mlt_dlnd00))) then
            ierr = -1
            if (.not. ss% report_ierr) return
!$OMP critical
            write(*,2) 'd_grada_dvb(mlt_dlnd00)', kz, d_grada_dvb(mlt_dlnd00)
            write(*,2) 'a_00', kz, a_00
            write(*,2) 'a_m1', kz, a_m1
            write(*,2) 'grada', kz, grada
            write(*,2) 'd_grada_00_dlnd', kz, d_grada_00_dlnd
            write(*,2) 'd_grada_00_dlnT', kz, d_grada_00_dlnT
            write(*,2) 'd_grada_m1_dlnd', kz, d_grada_m1_dlnd
            write(*,2) 'd_grada_m1_dlnT', kz, d_grada_m1_dlnT
            stop 'MLT'
!$OMP end critical
            return
         end if

         Pr = one_third*crad*T*T*T*T
         if (debug) write(*,1) 'Pr', Pr
         call set1_dvb(dPr_dvb, &
            0d0, 4d0*Pr*a_00*T_00/T, &
            0d0, 4d0*Pr*a_m1*T_m1/T, &
            4d0*Pr/T*d_a_00_dq00*(T_00-T_m1), &
            4d0*Pr/T*d_a_00_dqm1*(T_00-T_m1), &
            4d0*Pr/T*d_a_00_dqp1*(T_00-T_m1))
         
         !gradr = eval_Paczynski_gradr(P,opacity,L,m,cgrav,Pr,tau,T,r,rho)
         gradr = P*opacity*L / (16*pi*clight*m*cgrav*Pr)
         if (tau < 2d0/3d0) then ! B. Paczynski, 1969, Acta Astr., vol. 19, 1., eqn 14.
            s_gradr = (2*crad*T*T*T*sqrt(r))/(3*cgrav*m*rho)*pow_cr(L/(8*pi*boltz_sigma), 0.25d0) ! eqn 15
            f_gradr = 1 - 1.5d0*tau ! Paczynski, 1969, eqn 8
            dilute_factor = (1 + f_gradr*s_gradr*(4*pi*cgrav*clight*m)/(opacity*L))/(1 + f_gradr*s_gradr)
            gradr = gradr*dilute_factor
         end if
         gradr = gradr*gradr_factor
         d_gradr_dvb = gradr*(dP_dvb/P + d_opacity_dvb/opacity - dPr_dvb/Pr)
         d_gradr_dvb(mlt_dL) = gradr_factor*P*opacity / (16*pi*clight*m*cgrav*Pr)
         
         if (is_bad(gradr)) then
            ierr = -1
            if (.not. ss% report_ierr) return
!$OMP critical
            write(*,2) 'gradr', kz, gradr
            write(*,2) 'P', kz, P
            write(*,2) 'L', kz, L
            write(*,2) 'opacity', kz, opacity
            write(*,2) 'm', kz, m
            write(*,2) 'Pr', kz, Pr
            write(*,2) '16*pi*clight*m*cgrav*Pr', kz, 16*pi*clight*m*cgrav*Pr
            write(*,2) 'P*opacity*L', kz, P*opacity*L
            write(*,2) 'gradr_factor', kz, gradr_factor
            write(*,2) 'tau', kz, tau
            !write(*,2) '', kz, 
            stop 'mlt'
!$OMP end critical
         end if
         
         if (just_gradr) return

         grav = cgrav*m / (r*r)
         d_grav_dvb = 0d0
         d_grav_dvb(mlt_dlnR) = -2*grav
         d_grav_dvb(mlt_dq00) = grav*mstar/m
         
         if (grav < 0) then
            write(*,1) 'grav', grav
            write(*,1) 'cgrav', cgrav
            write(*,1) 'm', m
            write(*,1) 'r', r
            stop 'mlt'
         end if

         scale_height = P / (grav*rho)
         d_scale_height_dvb = scale_height*(dP_dvb/P - d_grav_dvb/grav - dRho_dvb/Rho)
         if (alt_scale_height) then
            ! consider sound speed*hydro time scale as an alternative scale height
            ! (this comes from Eggleton's code.)
            scale_height2 = sqrt(P/cgrav)/rho
            if (scale_height2 < scale_height) then
               scale_height = scale_height2
               d_scale_height_dvb = scale_height*(0.5d0*dP_dvb/P - dRho_dvb/Rho)
            end if
         end if
         H_P = scale_height

         if (scale_height <= 0d0 .or. is_bad(scale_height)) then
            ierr = -1
            return
!$OMP critical
            write(*,1) 'scale_height', scale_height
            stop 'set_convective_mixing'
!$OMP end critical
         end if

         if (is_bad(d_scale_height_dvb(mlt_dlnd00))) then
            ierr = -1
            return
!$OMP critical
            write(*,1) 'd_scale_height_dvb(mlt_dlnd00)', d_scale_height_dvb(mlt_dlnd00)
            stop 'set_convective_mixing'
!$OMP end critical
            return
         end if
         
         surf = pi4*r*r
         if (debug) write(*,1) 'surf', surf
         dsurf_dvb = 8*pi*r*dR_dvb
         
         ! Ledoux temperature gradient (same as Schwarzschild if composition term = 0)
         gradL = grada + gradL_composition_term
         d_gradL_dvb = d_grada_dvb ! ignore partials of composition term
         
         diff_grads = gradr - gradL ! convective if this is > 0
         d_diff_grads_dvb = d_gradr_dvb - d_gradL_dvb
         if (is_bad(d_diff_grads_dvb(mlt_dlnT00))) then
            ierr = -1
            return
!$omp critical
            write(*,1) 'd_grada_dvb(mlt_dlnT00)', d_grada_dvb(mlt_dlnT00)
            write(*,1) 'd_gradr_dvb(mlt_dlnT00)', d_gradr_dvb(mlt_dlnT00)
            write(*,1) 'd_gradL_dvb(mlt_dlnT00)', d_gradL_dvb(mlt_dlnT00)
            write(*,1) 'd_diff_grads_dvb(mlt_dlnT00)', d_diff_grads_dvb(mlt_dlnT00)
            stop 'MLT'
!$omp end critical
         end if

         Pg = P - Pr
         if (debug) write(*,1) 'Pg', Pg
         if (Pg < tiny) then
            call set_no_mixing
            return
         end if
         
         dPg_dvb = dP_dvb - dPr_dvb

         beta = Pg / P
         if (debug) write(*,1) 'beta', beta
         d_beta_dvb = beta*(dPg_dvb/Pg - dP_dvb/P)
         
         if (debug) write(*,1) 'scale_height', scale_height

         ! mixing length, Lambda
         Lambda = mixing_length_alpha*scale_height
         if (debug) write(*,1) 'Lambda', Lambda
         d_Lambda_dvb = mixing_length_alpha*d_scale_height_dvb
                  
         if (mixing_length_alpha <= 0) then
            call set_no_mixing
            return
         end if
         
         if (MLT_option == 'none') then
            call set_no_mixing
            return
         end if
         
         if (opacity < 1d-10 .or. P < 1d-20 .or. T < 1d-10 .or. Rho < 1d-20 &
               .or. m < 1d-10 .or. r < 1d-10 .or. cgrav < 1d-10 .or. &
               max_conv_vel == 0d0) then
            if (.false.) then
               write(*,2) 'special set no mixing', kz
               write(*,*) 'opacity < 1d-10', opacity < 1d-10
               write(*,*) 'P < 1d-20', P < 1d-20
               write(*,*) 'T < 1d-10', T < 1d-10
               write(*,*) 'Rho < 1d-20', Rho < 1d-20
               write(*,*) 'm < 1d-10', m < 1d-10
               write(*,*) 'r < 1d-10', r < 1d-10
               write(*,*) 'cgrav < 1d-10', cgrav < 1d-10
               write(*,*) 'max_conv_vel == 0d0', max_conv_vel == 0d0       
               write(*,*) "MLT_option == 'none' ", MLT_option == 'none'      
               stop 'mlt'
            end if
            call set_no_mixing
            return
         end if

         ! 'Q' param  C&G 14.24
         Q = chiT/chiRho
         dQ_dvb = Q*( d_chiT_dvb/chiT - d_chiRho_dvb/chiRho )
         if (Q <= 0) then
            call set_no_mixing
            return
         end if
                     
         radiative_conductivity = (4*crad*clight/3)*T*T*T / (opacity*rho) ! erg / (K cm sec)
         if (debug) write(*,1) 'radiative_conductivity', radiative_conductivity
         d_rc_dvb = radiative_conductivity*(3d0*dT_dvb/T - dRho_dvb/rho - d_opacity_dvb/opacity)

         if (smooth_gradT .and. diff_grads > gradT_smooth_low .and. diff_grads < gradT_smooth_high) then
            ! adjust diff_grads so its positive in the blending region,
            ! to compute a gradT from MLT in the Ledoux stable region.
            diff_grads = diff_grads - gradT_smooth_low
            call set_convective_mixing
            gradT_temp = gradT
            d_gradT_temp_dvb = d_gradT_dvb
            diff_grads = diff_grads + gradT_smooth_low
         end if
         
         if (diff_grads <= 0d0) then ! not convective (Ledoux stable)    
            call set_no_mixing ! also sets gradT = gradr    
            if (gradL_composition_term < 0) then ! composition unstable
               call set_thermohaline
               D = D_thrm
               d_D_dvb = d_D_thrm_dvb
               if (debug) write(*,1) 'after set_thermohaline D_thrm', D_thrm
               if (ss% conv_vel_flag .and. ss% conv_vel_ignore_thermohaline) then
                  conv_vel = 0d0
                  d_conv_vel_dvb = 0d0
                  !mixing_type = no_mixing
               end if
            else if (gradr > grada) then ! Schw unstable
               call set_semiconvection
               D = D_semi
               d_D_dvb = d_D_semi_dvb
               if (debug) write(*,1) 'after set_semiconvection D_semi', D_semi
               if (ss% conv_vel_flag .and. ss% conv_vel_ignore_semiconvection) then
                  conv_vel = 0d0
                  d_conv_vel_dvb = 0d0
                  !mixing_type = no_mixing
               end if
            else
               call time_limit_conv_vel
            end if
            if (debug) write(*,1) 'remove_small_D_limit', remove_small_D_limit
            if (D < remove_small_D_limit .or. is_bad(D)) then
               call set_no_mixing
            end if
            if (debug) write(*,1) 'final D', D
            unsmoothed_gradT = gradT
            d_unsmoothed_gradT_dvb = d_gradT_dvb
            ! blend at the boundary of convection to get a smooth gradT
            if (smooth_gradT .and. diff_grads >= gradT_smooth_low) then
               scale_factor = (diff_grads-gradT_smooth_low)/(-gradT_smooth_low)
               interp_factor = 10d0*pow_cr(scale_factor,3d0) &
                  -15*pow_cr(scale_factor,4d0)+6*pow_cr(scale_factor,5d0)
               dinterp_factor = (30d0*pow_cr(scale_factor,2d0) &
                  -60*pow_cr(scale_factor,3d0)&
                  +30*pow_cr(scale_factor,4d0)) &
                  /(-gradT_smooth_low)
               d_gradT_dvb = (1-interp_factor)*d_gradr_dvb - dinterp_factor * d_diff_grads_dvb * gradr &
                  + interp_factor*d_gradT_temp_dvb + dinterp_factor * d_diff_grads_dvb * gradT_temp
               gradT = (1-interp_factor)*gradr + interp_factor*gradT_temp
            end if
            if (conv_vel > 0d0) then
               D = conv_vel*Lambda/3     ! diffusion coefficient [cm^2/sec]
               if (debug) write(*,1) 'D', D
               d_D_dvb = (d_conv_vel_dvb*Lambda + conv_vel*d_Lambda_dvb)/3
               call set_conv_P
               if (mixing_type == no_mixing) then
                  mixing_type = softened_convective_mixing
               end if
            end if
            if (ss% conv_vel_flag) then
               if (normal_mlt_gradT_factor > 0d0) then
                  gradT_temp = gradT
                  d_gradT_temp_dvb = d_gradT_dvb
               end if
               if (normal_mlt_gradT_factor < 1d0) then
                  call revise_using_cv_var_variable
                  if (ierr /= 0) return
               end if
               if (normal_mlt_gradT_factor > 0d0 .and. normal_mlt_gradT_factor < 1d0) then
                  scale_factor = normal_mlt_gradT_factor
                  interp_factor = 10d0*pow_cr(scale_factor,3d0) &
                     -15*pow_cr(scale_factor,4d0)+6*pow_cr(scale_factor,5d0)
                  gradT = (1-interp_factor)*gradT + interp_factor*gradT_temp
                  d_gradT_dvb = (1-interp_factor)*d_gradT_dvb + interp_factor*d_gradT_temp_dvb
               end if
               if (kz > 0) then
                  if (ss% conv_vel(kz) > 0d0 .and. mixing_type == no_mixing) &
                     mixing_type = softened_convective_mixing
               end if
            end if
            return            
         end if
         
         call set_convective_mixing
         call set_conv_P
         if (quit) return
         
         D = conv_vel*Lambda/3     ! diffusion coefficient [cm^2/sec]
         if (debug) write(*,1) 'D', D
         d_D_dvb = (d_conv_vel_dvb*Lambda + conv_vel*d_Lambda_dvb)/3

         mixing_type = convective_mixing

         ! blend at the boundary of convection to get a smooth gradT
         unsmoothed_gradT = gradT
         d_unsmoothed_gradT_dvb = d_gradT_dvb
         if (smooth_gradT) then
            if (diff_grads <= gradT_smooth_high) then
               scale_factor = diff_grads/gradT_smooth_high
               interp_factor = -scale_factor*scale_factor*scale_factor*(-10d0 + scale_factor*(15d0 - 6d0*scale_factor))
               dinterp_factor = 30d0*(scale_factor - 1d0)*(scale_factor - 1d0)*scale_factor*scale_factor &
                  /gradT_smooth_high
               d_gradT_dvb = (1-interp_factor)*d_gradT_temp_dvb &
                  - dinterp_factor * d_diff_grads_dvb * gradT_temp &
                  + interp_factor*d_gradT_dvb + dinterp_factor * d_diff_grads_dvb * gradT
               gradT = (1-interp_factor)*gradT_temp + interp_factor*gradT

               !d_conv_vel_dvb = interp_factor*d_conv_vel_dvb + dinterp_factor * d_diff_grads_dvb * conv_vel
               !conv_vel = interp_factor * conv_vel
               !d_D_dvb = interp_factor*d_D_dvb + dinterp_factor * d_diff_grads_dvb * D
               !D = interp_factor * D
               !d_conv_P_dvb = interp_factor*d_conv_P_dvb + dinterp_factor * d_diff_grads_dvb * conv_P
               !conv_P = interp_factor * conv_P
            end if
         end if
         
         if (debug .or. D < 0) then
            write(*,*) 'get_gradT: convective_mixing'
            write(*,1) 'D', D
            write(*,1) 'conv_vel', conv_vel
            write(*,1) 'Pg/P', Pg/P
            write(*,1) 'H_P', H_P
            write(*,1) 'scale_height', scale_height
            write(*,1) 'scale_height/H_P', scale_height/H_P
            write(*,1) 'm/Msun', m/Msun
            write(*,1) 'r/Rsun', r/Rsun
            write(*,1) 'T', T
            write(*,1) 'rho', rho
            write(*,1) 'grada', grada
            write(*,1) 'chiT', chiT
            write(*,1) 'chiRho', chiRho
            write(*,2) 'mixing_type', mixing_type
            write(*,*)
            if (.not. debug) stop 'MLT: get_gradT'
         end if

         if (D < remove_small_D_limit .or. is_bad(D)) then
            call set_no_mixing
         end if

         if (ss% conv_vel_flag) then
            if (normal_mlt_gradT_factor > 0d0) then
               gradT_temp = gradT
               d_gradT_temp_dvb = d_gradT_dvb
            end if
            if (normal_mlt_gradT_factor < 1d0) then
               call revise_using_cv_var_variable
               if (ierr /= 0) return
            end if
            if (normal_mlt_gradT_factor > 0d0 .and. normal_mlt_gradT_factor < 1d0) then
               scale_factor = normal_mlt_gradT_factor
               interp_factor = 10d0*pow_cr(scale_factor,3d0) &
                  -15*pow_cr(scale_factor,4d0)+6*pow_cr(scale_factor,5d0)
               gradT = (1-interp_factor)*gradT + interp_factor*gradT_temp
               d_gradT_dvb = (1-interp_factor)*d_gradT_dvb + interp_factor*d_gradT_temp_dvb
            end if
            if (kz > 0) then ! e.g. pre_ms_model calls with kz = 0
               if (ss% conv_vel(kz) > 0d0 .and. mixing_type == no_mixing) then
                  mixing_type = softened_convective_mixing
               end if
            end if
         end if
         
         contains

         
         subroutine set1_dvb(dvb, dlnd00, dlnT00, dlndm1, dlndTm1, dq00, dqm1, dqp1)
            real(dp), intent(inout) :: dvb(nvbs)
            real(dp), intent(in) :: dlnd00, dlnT00, dlndm1, dlndTm1, dq00, dqm1, dqp1
            dvb(mlt_dlnd00) = dlnd00
            dvb(mlt_dlnT00) = dlnT00
            dvb(mlt_dlndm1) = dlndm1
            dvb(mlt_dlnTm1) = dlndTm1
            dvb(mlt_dlnR) = 0d0
            dvb(mlt_dL) = 0d0
            dvb(mlt_cv_var) = 0d0
            dvb(mlt_dq00) = dq00
            dvb(mlt_dqm1) = dqm1
            dvb(mlt_dqp1) = dqp1
         end subroutine set1_dvb
         
         subroutine set_conv_P
            if (conv_vel == 0d0) then
               conv_P = 0d0
               d_conv_P_dvb = 0d0
            else
               conv_P = rho*conv_vel*conv_vel/(3*P) ! see C&G (14.69) for Pturb/P
               d_conv_P_dvb = rho*2*d_conv_vel_dvb*conv_vel/(3*P) + &
                  dRho_dvb*conv_vel*conv_vel/(3*P) - conv_P*dP_dvb/P
            end if
         end subroutine set_conv_P
         
         subroutine set_thermohaline
            real(dp) :: kipp_D, tau, Pr, BGS_C, Nu_mu, l2, lmda, phi, r_guess, &
            sqrt_1_plus_phi, sqrt_Pr
   
            logical, parameter :: dbg = .false.
            include 'formats'
                        
            if (dbg) write(*,*) 'set_thermohaline ' // trim(thermohaline_option)
            
            diff_grad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable               
            K = 4*crad*clight*T*T*T/(3*opacity*rho) ! thermal conductivity
            
            if (thermohaline_option == 'Kippenhahn') then
            
               ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
               D_thrm = -thermohaline_coeff*3*K/(2*rho*cp)*gradL_composition_term/diff_grad
            
            else if (thermohaline_option == 'Brown_Garaud_Stellmach_13' .or. &
                     thermohaline_option == 'Traxler_Garaud_Stellmach_11') then
                     
               call get_diff_coeffs(K_T,K_mu,nu)

               R0 = (gradr - grada)/gradL_composition_term
               Pr = nu/K_T
               tau = K_mu/K_T
               r_th = (R0 - 1d0)/(1d0/tau - 1d0)

               if (r_th >= 1d0) then ! stable if R0 >= 1/tau
                  D_thrm = 0d0
               else if (thermohaline_option == 'Traxler_Garaud_Stellmach_11') then 
                  ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
                  ! also see Denissenkov. ApJ 723:563–579, 2010.
                  D_thrm = 101d0*sqrt(K_mu*nu)*exp_cr(-3.6d0*r_th)*pow_cr(1d0 - r_th,1.1d0) ! eqn 24
               else             
                  ! if (thermohaline_option == 'Brown_Garaud_Stellmach_13') then
                  D_thrm = K_mu*(Numu(R0,r_th,pr,tau) - 1d0)
                  ! evbauer 07/18: changed from K_mu*Numu(R0,r_th,pr,tau), Pascale signed off
               endif
               D_thrm = thermohaline_coeff*D_thrm
               
            else
                 
               D_thrm = 0
               ierr = -1
               write(*,*) 'unknown value for MLT thermohaline_option' // trim(thermohaline_option)
               return   
               
            end if
            
            conv_vel = 3*D_thrm/Lambda
            d_conv_vel_dvb = 0
            d_D_thrm_dvb = 0
            
            call time_limit_conv_vel
            if (conv_vel > max_conv_vel) conv_vel = max_conv_vel
            if (init_conv_vel /= conv_vel) D_thrm = conv_vel*Lambda/3
            
            if (D_thrm < min_D_th .or. D_thrm <= 0) then
               call set_no_mixing
               return
            end if
            
            mixing_type = thermohaline_mixing 
            
         end subroutine set_thermohaline
         
         subroutine time_limit_conv_vel
            real(dp), dimension(nvbs) :: d_vconv_accel_dvb, d_brunt_timescale_dvb
            real(dp), dimension(nvbs) :: d_diff_grads_s_dvb
            real(dp) :: new_conv_vel, vconv_accel, l
            real(dp) :: brunt_timescale, diff_grads_s
            include 'formats'

            init_conv_vel = conv_vel
            if (dt <= 0 .or. conv_vel <= prev_conv_vel .or. conv_vel == 0d0) return
            !if (dt <= 0) return

            if (g_theta > 0) then ! max accel is grav*g_theta
               if (conv_vel > prev_conv_vel) then
                  !increase convective velocity as needed
                  new_conv_vel = prev_conv_vel + dt*grav*g_theta
                  if (new_conv_vel < conv_vel) then
                     conv_vel = new_conv_vel
                     d_conv_vel_dvb = dt*d_grav_dvb*g_theta
                  end if
               else
                  !reduce convective velocity as needed
                  new_conv_vel = prev_conv_vel - dt*grav*g_theta
                  if (new_conv_vel > conv_vel) then 
                     conv_vel = new_conv_vel
                     d_conv_vel_dvb = -dt*d_grav_dvb*g_theta
                  end if
               end if
            else
               ! Arnett, W.D., 1969, Ap. and Space Sci, 5, 180.
               l = Lambda
               if (conv_vel > 0d0) then
                  vconv_accel = 2d0*(conv_vel*conv_vel-prev_conv_vel*prev_conv_vel)/l
                  d_vconv_accel_dvb = 4d0*conv_vel*d_conv_vel_dvb/l &
                     -2d0*(conv_vel*conv_vel-prev_conv_vel*prev_conv_vel)*d_Lambda_dvb/l*l
                  if (conv_vel > prev_conv_vel) then
                     !increase convective velocity as needed
                     new_conv_vel = prev_conv_vel + dt*vconv_accel
                     if (new_conv_vel < conv_vel) then
                        conv_vel = new_conv_vel
                        d_conv_vel_dvb = dt*d_vconv_accel_dvb
                     end if
                  else
                     !reduce convective velocity as needed
                     new_conv_vel = prev_conv_vel + dt*vconv_accel
                     if (new_conv_vel > conv_vel) then
                        conv_vel = new_conv_vel
                        d_conv_vel_dvb = dt*d_vconv_accel_dvb
                     end if
                  end if
               else ! conv_vel == 0d0
                  if (prev_conv_vel == 0d0) return
                  if (smooth_gradT) then
                     diff_grads_s = diff_grads - gradT_smooth_low
                  else
                     diff_grads_s = diff_grads
                  end if
                  d_diff_grads_s_dvb = d_diff_grads_dvb
                  if (diff_grads_s < 0d0) then
                     brunt_timescale = 1/sqrt(-chiT/chiRho*diff_grads_s*grav/scale_height)
                  else
                     brunt_timescale = 1d99
                  end if
                  if (l/prev_conv_vel > brunt_timescale) then
                     ! reduce conv_vel on the brunt_timescale
                     new_conv_vel = prev_conv_vel*exp_cr(-dt/brunt_timescale)
                     d_brunt_timescale_dvb = 0.5d0*brunt_timescale*(&
                        -d_chiT_dvb/chiT + d_chiRho_dvb/chiRho - d_diff_grads_s_dvb/diff_grads_s &
                        -d_grav_dvb/grav + d_scale_height_dvb/scale_height)
                     d_conv_vel_dvb = d_brunt_timescale_dvb*&
                        (dt*exp_cr(-dt/brunt_timescale)*prev_conv_vel)/powi_cr(brunt_timescale,2)
                  else
                     new_conv_vel = prev_conv_vel*l/(2*prev_conv_vel*dt+l)
                     d_conv_vel_dvb = (2*dt*prev_conv_vel*prev_conv_vel)/&
                        (4*dt*dt*prev_conv_vel*prev_conv_vel+4*dt*l*prev_conv_vel+l*l) &
                        *d_Lambda_dvb
                  end if
                  conv_vel = new_conv_vel
               end if
            end if

            if (conv_vel > max_conv_vel) then
               conv_vel = max_conv_vel
               d_conv_vel_dvb = 0d0
            end if

         end subroutine time_limit_conv_vel

         subroutine get_diff_coeffs(kt,kmu,vis)

         use chem_def, only: chem_isos
         use crlibm_lib

         real(dp) :: kt,kmu,vis,qe4
         real(dp) :: loglambdah,loglambdacx,loglambdacy,ccx,ccy
         real(dp) :: Bcoeff
         real(dp) :: chemA,chemZ,acx,acy       
         real(dp), parameter :: sqrt5 = sqrt(5d0)
                 
         kt = K/(Cp*rho)       ! thermal diffusivity (assumes radiatively dominated)
         qe4=qe*qe*qe*qe
         
         ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
         loglambdah = -19.26 - 0.5*log_cr(rho) + 1.5*log_cr(T) - 0.5*log_cr(1. + 0.5*(1+xh)) 
         nu_rad = 4*crad*T*T*T*T/(15*clight*opacity*rho*rho) ! radiative viscosity
         nu_mol = 0.406*sqrt(amu)*pow_cr(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
         ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
         ! Assumes pure H. Still trying to work out what it would be for a mixture. 
         vis = nu_mol + nu_rad   ! total viscosity

         ! The following is from Proffitt & Michaud, 1993.
         ! Their constant B (equation 15)
         Bcoeff = (15./16.)*sqrt(2.*amu/(5*pi))*pow_cr(boltzm,2.5d0)/qe4
         ! Extract what species drives the thermohaline concvection
         chemA = chem_isos%Z_plus_N(dominant_iso_for_thermohaline)
         chemZ = chem_isos%Z(dominant_iso_for_thermohaline)

         if(chemZ.gt.2) then
         ! This is if the driving chemical is NOT He.
            ! Log Lambda for H-dominant chem mixture (equation 10)
            loglambdacx = loglambdah - log_cr(chemz)  
            ! Log Lambda for He-dominant chem mixture (equation 10)
            loglambdacy = loglambdah - log_cr(2.*chemz)
            ! Calculation of C_ij coeffs (equation 12)
            ccx = log_cr(exp_cr(1.2*loglambdacx)+1.)/1.2
            ccy = log_cr(exp_cr(1.2*loglambdacy)+1.)/1.2
            ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He
            acx = (1.*chemA)/(1.+chemA)
            acy = 4*chemA/(4.+chemA)
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = 2*Bcoeff*pow_cr(T,2.5d0)/(sqrt5*rho*chemZ*chemZ)/ &
               (xh*sqrt(acx)*ccx + (1-xh)*sqrt(acy)*ccy)

         else
            ! Log Lambda for H-He mixture (equation 10)
            loglambdah = -19.26 - log_cr(2.d0) - 0.5*log_cr(rho) + &
               1.5*log_cr(T) - 0.5*log_cr(1. + 0.5*(1+xh)) 
            ! Calculation of C_ij coeffs (equation 12)
            ccy = log_cr(exp_cr(1.2*loglambdah)+1.)/1.2
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = (Bcoeff*pow_cr(T,2.5d0)/(rho*ccy))*(3+xh)/((1+xh)*(3+5*xh)*(0.7+0.3*xh))
            
         endif
         ! write(57,*) kt,kmu,vis,chemZ

         end subroutine get_diff_coeffs

         double precision function numu(R0,r_th,prandtl,diffratio)
         !Function calculates Nu_mu from input parameters, following Brown et al. 2013.
         !Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

            real(dp), intent(in) :: R0,r_th,prandtl,diffratio
            real(dp) :: maxl2,maxl,lambdamax
            real(dp) :: myvars(2)
            integer :: ierr, iter

            ! Initialize guess using estimates from Brown et al. 2013
            call analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
                  
            myvars(1) = maxl
            myvars(2) = lambdamax

           !Call Newton relaxation algorithm
           call NR(myvars,prandtl,diffratio,R0,ierr)
          
           !If the growth rate is negative, then try another set of parameters as first guess.  
           !Repeat as many times as necessary until convergence is obtained.
           iter = 1
           do while((myvars(2)<0).or.(ierr /= 0)) 
              !write(*,*) 'Alternative', r_th,prandtl,diffratio,iter
           !Reset guess values
              myvars(1) = maxl
              myvars(2) = lambdamax
           !Call relaxation for slightly different Pr, tau, R0.
              call NR(myvars,prandtl*(1.+iter*1.d-2),diffratio,R0/(1.+iter*1.d-2),ierr)
           !If it converged this time, call NR for the real parameters.
              if(ierr.eq.0) call NR(myvars,prandtl,diffratio,R0,ierr)
              !write(*,*) prandtl,diffratio,R0,myvars(1),myvars(2),ierr
              !Otherwise, increase counter and try again.
              iter = iter + 1            
           enddo

           !Plug solution into "l^2" and lambda.
           maxl2 = myvars(1)*myvars(1)
           lambdamax = myvars(2) 
           !write(*,*) prandtl,diffratio,r_th,maxl2,lambdamax

           !Calculate Nu_mu using Formula (33) from Brown et al, with C = 7.
           numu = 1. + 49.*lambdamax*lambdamax/(diffratio*maxl2*(lambdamax+diffratio*maxl2))

         return
         end function numu 

         subroutine thermohaline_rhs(myx,myf,myj,prandtl,diffratio,R0)
         ! This routine is needed for the NR solver.
         ! Inputs the two following equations for lambda and maxl2:
         ! lambda^3 + a_2 lambda^2 + a_1 lambda + a_0 = 0 (eq. 19 of Brown et al.)
         ! b_2 lambda^2 + b_1 lambda + b_0 = 0 (eq. 20 of Brown et al.)
         ! Inputs f, the equations, and j, their jacobian.
         ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

         real(dp), intent(in) :: myx(2),  prandtl, diffratio, R0
         real(dp), intent(out) :: myf(2), myj(2,2)
         real(dp) :: a_2,a_1,a_0,b_2,b_1,b_0,myterm,myx1_2,myx1_3,myx1_4
 
         !This inputs the coefficients.
         b_2 = 1.+prandtl+diffratio
         myx1_2 = myx(1)*myx(1)
         myx1_3 = myx1_2*myx(1)
         myx1_4 = myx1_3*myx(1)
         a_2 = myx1_2*b_2
         myterm = diffratio*prandtl+prandtl+diffratio
         b_1 = 2*myx1_2*myterm
         a_1 = myx1_4*myterm + prandtl*(1. - (1./R0))
         b_0 = 3.*myx1_4*diffratio*prandtl + prandtl*(diffratio - (1./R0))
         a_0 = myx1_4*myx1_2*diffratio*prandtl + myx1_2*prandtl*(diffratio - (1./R0))

!         write(*,*) a_2,a_1,a_0,b_2,b_1,b_0

         !These are equations 19 and 20
         myf(1) = ((myx(2) + a_2)*myx(2) + a_1)*myx(2) + a_0
         myf(2) = b_2*myx(2)*myx(2) + b_1*myx(2) + b_0
         
         !These are their Jacobians for the NR relaxation.
         myj(1,1) = 2*myx(1)*b_2*myx(2)*myx(2) + &
            4*myx1_3*myterm*myx(2) + 6*myx1_4*myx(1)*diffratio*prandtl   &
              + 2*myx(1)*prandtl*(diffratio - (1./R0))
         myj(1,2) = 3*myx(2)*myx(2) + 2*a_2*myx(2) + a_1
         myj(2,1) = 4*myx(1)*myterm*myx(2) + 12.*myx1_3*diffratio*prandtl
         myj(2,2) = 2*b_2*myx(2) + b_1
 
         return
         end subroutine thermohaline_rhs               

         subroutine analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
         !Inputs analytical estimates for l and lambda from Brown et al. 2013.

         real(dp) :: prandtl, diffratio, maxl, lambdamax, r_th, phi, maxl4, maxl6
         
         phi = diffratio/prandtl

         if(r_th .lt. 0.5) then
            if(r_th .gt. prandtl) then
               maxl = pow_cr((1./(1.+phi)) - 2.*dsqrt(r_th*phi)/pow_cr(1.+phi,2.5d0),0.25d0)   
                  ! Equation (B14)
               maxl4 = maxl*maxl*maxl*maxl
               maxl6 = maxl4*maxl*maxl
               lambdamax = 2*prandtl*phi*maxl6/(1.-(1.+phi)*maxl4)    ! Equation (B11)
            else
               maxl = dsqrt(dsqrt(1./(1.+phi)) - dsqrt(prandtl)*(1.+phi/((1.+phi)*(1.+phi))))  
                  ! Equation (B5)
               lambdamax = dsqrt(prandtl) - prandtl*dsqrt(1.+phi)   !Equation (B5)
            endif
         else
            maxl = pow_cr((1d0/3d0)*(1d0-r_th) + (1d0-r_th)*(1d0-r_th)*(5d0-4d0*phi)/27d0,0.25d0)
               ! Equation (B19) carried to next order (doesn't work well otherwise)
            maxl4 = maxl*maxl*maxl*maxl
            maxl6 = maxl4*maxl*maxl
            lambdamax = 2d0*prandtl*phi*maxl6/(1d0-(1d0+phi)*maxl4) ! Equation (B11)
         endif
         if(lambdamax<0) then   ! shouldn't be needed, but just as precaution
            maxl = 0.5d0
            lambdamax = 0.5d0
         endif
         
         return
         end subroutine analytical_estimate_th

         subroutine NR(xrk,prandtl,diffratio,R0,ierr)
         ! Newton Relaxation routine used to solve cubic & quadratic in thermohaline case.
         ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

         real(dp), parameter :: acy = 1.d-13 ! accuracy of NR solution.
         integer, parameter :: niter = 20  ! max number of iterations allowed before giving up.
         integer, parameter :: &  !array dimension input parameters for dgesvx
               n = 2, &
               nrhs = 1, &
               lda = n, &
               ldaf = n, &
               ldb = n, &
               ldx = n

         integer :: iter,ierr
         real(dp) :: xrk(2), f(2) ! Functions f 
         real(dp) :: j(2,2) ! Jacobian
         real(dp) :: err,errold ! Error at each iteration
         real(dp) :: x1_sav,x2_sav
         real(dp) :: prandtl, diffratio, R0
         real(dp) :: A(lda,n), AF(ldaf,n), R(n), C(n), B(ldb,nrhs), X(ldx,nrhs), &
               rcond, ferr(nrhs), berr(nrhs), work(4*n)
         character :: fact, trans, equed
         integer :: ipiv(n), iwork(n)

         include 'formats'

         !Initialize flags and other counters.
         ierr = 0
         iter = 0
         err = 0d0
         errold = 0d0
         !Save input guess (probably not necessary here but useful in other routines)
         x1_sav = xrk(1)
         x2_sav = xrk(2)

         !While error is too large .and. decreasing, iterate.
         do while ((err.gt.acy).and.(ierr.eq.0).and.(iter.lt.niter))
            call thermohaline_rhs(xrk,f,j,prandtl,diffratio,R0)    
            
            fact = 'E'
            trans = 'N'
            equed = ''
               
            A  = j
            B(1,1) = f(1)
            B(2,1) = f(2)

            call dgesvx( fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, &
               equed, r, c, B, ldb, x, ldx, rcond, ferr, berr, &
               work, iwork, ierr )
 
            if (ierr /= 0) then
               !write(*,*) 'dgesvx failed in thermohaline routine', iter
               !write(*,2) j(1,1),j(1,2)
               !write(*,2) j(2,1),j(2,2)
            else
               iter = iter + 1
               f(1) = X(1,1)
               f(2) = X(2,1)
               err = dsqrt(f(1)*f(1)+f(2)*f(2)) ! Calculate the new error
               ! If, after a while, the error is still not decreasing, give up and exit NR.
               ! Otherwise, continue.
               if((iter.gt.5).and.(err.gt.errold)) then              
                  ! Write(*,2) 'Error not decreasing at iter', iter, err, errold
                  ierr = 1
                  ! Reset xs and exit loop.
                  xrk(1) = x1_sav
                  xrk(2) = x2_sav                   
               else
                  xrk = xrk - f ! The solution is now in f, so update x 
                  errold = err
               endif
            endif
         enddo
            
         !if(iter.eq.niter) write(*,2) 'Failed to converge'
         return
         end subroutine NR
         
         subroutine set_convective_mixing
            ! need to set gradT, d_gradT_dvb, conv_vel, d_conv_vel_dvb
            include 'formats.dek'
            real(dp) ff1, ff2, ff3, ff4, ff5, aa, bb, y0, xres, a1, a2, sqrt_x
            real(dp) :: A_0, A_1, A_2, A_numerator, A_denom, inv_sqrt_x
            real(dp), dimension(nvbs) :: &
               dA_0_dvb, dA_1_dvb, dA_2_dvb, dA_numerator_dvb, dA_denom_dvb, &
               d_inv_sqrt_x_dvb
            
            real(qp) :: q1, q2, q3
            real(qp), dimension(nvbs) :: qd_dvb1, qd_dvb2, qd_dvb3
            
            real(dp), parameter :: two_13 = 1.2599210498948730d0 ! = pow_cr(2d0,1d0/3d0)
            real(dp), parameter :: four_13 = 1.5874010519681994d0 ! = pow_cr(4d0,1d0/3d0)

! options for MLT_option are:
!    'ML1'        Bohm-Vitense 1958 MLT
!    'ML2'        Bohm and Cassinelli 1971 MLT
!    'Mihalas'    Mihalas 1978, Kurucz 1979 MLT
!    'Henyey'     Henyey, Rardya, and Bodenheimer 1965 MLT
! Values of the f1..f4 coefficients are taken from Table 1 of Ludwig et al. 1999, A&A, 346, 111
! with the following exception: their value of f3 for Henyey convection is f4/8 when it should be
! 8*f4, i.e., f3=32*pi**2/3 and f4=4*pi**2/3. f3 and f4 are related to the henyey y parameter, so
! for the 'Henyey' case they are set based on the value of Henyey_y_param. The f1..f4 parameters
! have been renamed with a double ff, i.e., ff1..ff4, to avoid a collision of variable names with
! f1 in the cubic root solver.
            
            quit = .false.

            x = Q*Rho / (2*P)
            d_x_dvb = x*(drho_dvb/rho + dQ_dvb/Q - dP_dvb/P)
         
            convective_conductivity = Cp*grav*Lambda*Lambda*Rho*(sqrt(x)) / 9 ! erg / (K cm sec)
            
            if (convective_conductivity < 0) then
               ierr = -1
               return
               write(*,1) 'MLT error: convective_conductivity', convective_conductivity
               write(*,1) 'Cp', Cp
               write(*,1) 'grav', grav
               write(*,1) 'Lambda', Lambda
               write(*,1) 'Rho', Rho
               write(*,1) 'x', x
               stop 'mlt'
            end if
            
            if (debug) write(*,1) 'convective_conductivity', convective_conductivity
            d_cc_dvb = convective_conductivity* &
                 (d_Cp_dvb/Cp + d_grav_dvb/grav + &
                     2*d_Lambda_dvb/Lambda + dRho_dvb/rho + d_x_dvb/(2*x))

            if (MLT_option == 'Cox') then ! this assumes optically thick

               a0 = 9d0/4d0
               d_a0_dvb = 0d0

               ! 'A' param is ratio of convective to radiative conductivities   C&G 14.98
               A = convective_conductivity / radiative_conductivity !  unitless.

               if (debug) write(*,1) 'A', A
               dA_dvb = (d_cc_dvb - d_rc_dvb*A) / radiative_conductivity
               
               if (A < 0 .or. is_bad(A)) then
                  write(*,*) "MLT_option == 'Cox'", MLT_option == 'Cox'
                  write(*,1) 'A', A
                  write(*,1) 'convective_conductivity', convective_conductivity
                  write(*,1) 'radiative_conductivity', radiative_conductivity
                  stop 'mlt'
               end if

            else
            
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1./Henyey_nu_param
                  ff2=1./2.
                  ! popular values for y are 1/3 or 3/(4*pi**2)
                  ff3=8./Henyey_y_param
                  ff4=1./Henyey_y_param
               case ('ML1')
                  ff1=1./8.
                  ff2=1./2.
                  ff3=24.0
                  ff4=0.0
               case ('ML2')
                  ff1=1.
                  ff2=2.
                  ff3=16.0
                  ff4=0.0
               case ('Mihalas')
                  ff1=1./8.
                  ff2=1./2.
                  ff3=16.0
                  ff4=2.0
               case default
                  write(*,'(3a)') 'Error: ',trim(MLT_option), &
                     ' is not an allowed MLT version for convection'
                  write(*,*)
                  return
               end select
            
               omega = Lambda*Rho*opacity !dimensionless
               d_omega_dvb = omega*( d_Lambda_dvb/Lambda + dRho_dvb/Rho + d_opacity_dvb/opacity)

               ! the variable theta in no longer needed
               ! theta = omega / ( 1d0 + Henyey_y_param*omega**2 )
               ! d_theta_dvb = d_omega_dvb*(1d0 - Henyey_y_param*omega**2 ) /
               ! ( ( 1d0 + Henyey_y_param*omega**2 )**2 )

               ! a0 = 0.75d0*omega*theta
               !d_a0_dvb = a0*( d_omega_dvb/omega + d_theta_dvb/theta )
               a0 = (3./16.)*ff2*ff3/(1.+ff4/(omega*omega))
               d_a0_dvb = a0*2*ff4*d_omega_dvb/(ff4 + omega*omega)/omega
               !ignore d_a0_dvb
               d_a0_dvb = 0d0

               ! A = sqrt(P*Q*rho/Henyey_nu_param)*(Cp*mixing_length_alpha)/
               !        (2*crad*clight*T**3*theta)               
               A_0 = sqrt(ff1*P*Q*rho)
               dA_0_dvb = ff1*(dP_dvb*Q*rho + P*dQ_dvb*rho + P*Q*drho_dvb)/(2*A_0)
               
               A_1 = 4*A_0*Cp
               dA_1_dvb = 4*(dA_0_dvb*Cp + A_0*d_Cp_dvb)
               
               A_2 = mixing_length_alpha*omega*(1.+ff4/(omega*omega))
               dA_2_dvb = mixing_length_alpha*(1-ff4/(omega*omega))*d_omega_dvb

               if (is_bad(dA_2_dvb(mlt_dlnT00))) then
                  ierr = -1
                  write(*,1) 'dA_2_dvb(mlt_dlnT00)', dA_2_dvb(mlt_dlnT00)
                  stop 'MLT'
                  return
               end if
               
               A_numerator = A_1*A_2
               dA_numerator_dvb = dA_1_dvb*A_2 + A_1*dA_2_dvb

               if (is_bad(dA_numerator_dvb(mlt_dlnT00))) then
                  ierr = -1
                  return
                  write(*,1) 'dA_numerator_dvb(mlt_dlnT00)', dA_numerator_dvb(mlt_dlnT00)
                  stop 'MLT'
               end if
                     
               A_denom = ff3*crad*clight*T*T*T
               dA_denom_dvb = A_denom*3*dT_dvb/T
               
               A = A_numerator/A_denom                     
               dA_dvb = dA_numerator_dvb/A_denom - A_numerator*dA_denom_dvb/(A_denom*A_denom)

               if (is_bad(dA_dvb(mlt_dlnT00))) then
                  ierr = -1
                  return
                  write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
                  stop 'MLT'
               end if
               
               if (A < 0) then
                  ierr = -1
                  return
                  write(*,1) 'A', A
                  write(*,1) 'A_numerator', A_numerator
                  write(*,1) 'A_denom', A_denom
                  write(*,1) 'A_1', A_1
                  write(*,1) 'ff3', ff3
                  write(*,1) 'A_0', A_0
                  stop 'mlt'
               end if
            
            end if

            ! 'B' param  C&G 14.81
            Bcubed = (A*A / a0)*diff_grads         
            d_Bcubed_dvb = (A*A / a0)*d_diff_grads_dvb + (2*A*dA_dvb / a0)*diff_grads &
               - Bcubed*d_a0_dvb/a0

            if (is_bad(d_Bcubed_dvb(mlt_dlnT00))) then
               ierr = -1
               return
!$omp critical
               write(*,1) 'd_diff_grads_dvb(mlt_dlnT00)', d_diff_grads_dvb(mlt_dlnT00)
               write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
               write(*,1) 'd_Bcubed_dvb(mlt_dlnT00)', d_Bcubed_dvb(mlt_dlnT00)
               write(*,1) 'a0', a0
               write(*,1) 'diff_grads', diff_grads
               write(*,1) 'A', A
               write(*,1) 'Bcubed', Bcubed
               stop 'MLT'
!$omp end critical
            end if
         
            if (debug) write(*,1) 'Bcubed', Bcubed

            ! now solve cubic equation for convective efficiency, Gamma
            ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82, 
            ! rewritten in terms of Gamma
            ! leave it to Mathematica to find an expression for the root we want (with a0 = 9/4)
         
            delta = a0*Bcubed
            d_delta_dvb = a0*d_Bcubed_dvb + Bcubed*d_a0_dvb
         
            if (debug) write(*,1) 'a0', a0
            if (debug) write(*,1) 'delta', delta
      
            f = -2 + 9*a0 + 27*a0*a0*delta
            if (debug) write(*,1) 'f', f
            if (f > 1d100) then
               f0 = f
               d_f0_dvb = 27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb
            else
               f0 = f*f + 4*(-1 + 3*a0)*(-1 + 3*a0)*(-1 + 3*a0)
               if (f0 < 0d0) then
                  f0 = f
                  d_f0_dvb = 27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb
               else
                  f0 = sqrt(f0)         
                  !d_f0_dvb = (f*(27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb) &
                  !   + 18*(-1 + 3*a0)*(-1 + 3*a0)*d_a0_dvb)/f0!27*a0*a0*f*d_delta_dvb / f0
                  d_f0_dvb = 27*a0*a0*f*d_delta_dvb / f0 &
                     + (f*(9+54*a0*delta)+ 18*(-1 + 3*a0)*(-1 + 3*a0))*d_a0_dvb/f0
               end if
            end if
         
            if (debug) write(*,1) 'f0', f0

            f1 = -2 + 9*a0 + 27*a0*a0*delta + f0  

            if (is_bad(f1)) then
               ierr = -1
               return
!$omp critical
               write(*,1) 'f1', f1
               write(*,1) 'a0', a0
               write(*,1) 'delta', delta
               write(*,1) 'f0', f0
               stop 'MLT: bad f1'
!$omp end critical
            end if   

            if (f1 < 0) then
               call set_no_mixing
               call time_limit_conv_vel
               return
            end if   
            f1 = pow_cr(f1,one_third)     
            d_f1_dvb = (27*a0*a0*d_delta_dvb + (9+54*a0*delta)*d_a0_dvb + d_f0_dvb) / (3*f1*f1)
            f2 = 2*two_13*(1 - 3*a0) / f1       
            d_f2_dvb = -f2*d_f1_dvb / f1 - 6*two_13*d_a0_dvb / f1

            Gamma = (four_13*f1 + f2 - 2) / (6*a0)
            d_Gamma_dvb = (four_13*d_f1_dvb + d_f2_dvb) / (6*a0) - Gamma*d_a0_dvb/a0

            if (is_bad(Gamma)) then
               ierr = -1
               return
!$omp critical
               write(*,1) 'Gamma', Gamma
               write(*,1) 'f1', f1
               write(*,1) 'f2', f2
               write(*,1) 'a0', a0
               write(*,1) 'd_f1_dvb', d_f1_dvb
               write(*,1) 'd_f2_dvb', d_f2_dvb
               stop 'MLT: bad f1'
!$omp end critical
!               call set_no_mixing
!               quit = .true.
!               return
            end if

            if (Gamma < 0) then
               call set_no_mixing
               call time_limit_conv_vel
               return
            end if
            
            ! average convection velocity, vbar   C&G 14.86b
            ! vbar = vsound*Sqrt(Q)*alpha*Gamma / (2*Sqrt(2*Gamma1)*A)
            ! vsound = Sqrt(Gamma1*P / rho), so
            ! vbar = Sqrt(Q*P / (8*rho))*alpha*Gamma / A

            x = Q*P / (8*rho)
            sqrt_x = sqrt(x)
            conv_vel = mixing_length_alpha*sqrt_x*Gamma / A
            if (conv_vel > max_conv_vel) then
               conv_vel = max_conv_vel
               d_conv_vel_dvb = 0
            else
               d_conv_vel_dvb = 0.5d0*conv_vel* &
                 (-2*dA_dvb/A + 2*d_Gamma_dvb/Gamma + &
                 dP_dvb/P + dQ_dvb/Q - drho_dvb/rho)
            end if
            
            

            if (test_partials) then
               ss% hydro_test_partials_val = conv_vel
               ss% hydro_test_partials_var = ss% i_lum
               ss% hydro_test_partials_dval_dx = d_conv_vel_dvb(mlt_dL)
            end if
            
            
            

            conv_tau = 1d99
            call time_limit_conv_vel
            
            if (init_conv_vel /= conv_vel .or. conv_vel == max_conv_vel) then
               ! need to recalculate Gamma to match modified conv_vel
               if (A <= 1d-99 .or. sqrt_x <= 1d-99) then
                  Gamma = 1d25
                  d_Gamma_dvb = 0
                  if (dbg) write(*,1) 'A or sqrt_x too small', A, sqrt_x
               else      
                  if (dbg) write(*,*) 'recalculate Gamma to match modified conv_vel'
                  inv_sqrt_x = 1d0/sqrt_x
                  d_inv_sqrt_x_dvb = &
                     (Q*P*drho_dvb/rho - Q*dP_dvb - P*dQ_dvb)/(16d0*x*sqrt_x*rho)
                  Gamma = conv_vel*A*inv_sqrt_x/mixing_length_alpha  
                  !d_Gamma_dvb = ( &
                  !   d_conv_vel_dvb*A*inv_sqrt_x + &
                  !   conv_vel*dA_dvb*inv_sqrt_x + &
                  !   conv_vel*A*d_inv_sqrt_x_dvb)/mixing_length_alpha  
                  ! the "correct" form breaks edep in example_ccsn_IIp.
                  d_Gamma_dvb = 0d0
               end if
            end if

            if (dbg) &
               write(*,1) 'prev/init init/final conv_vel', &
                  prev_conv_vel/init_conv_vel, &
                  prev_conv_vel, init_conv_vel, conv_vel
            
            if (debug) write(*,1) 'conv_vel', conv_vel
            if (conv_vel < 0) then
               ierr = -1
               return
!$omp critical
               write(*,1) 'conv_vel', conv_vel
               write(*,1) 'mixing_length_alpha', mixing_length_alpha
               write(*,1) 'x', x
               write(*,1) 'A', A
               write(*,1) 'Gamma', Gamma
               stop 'MLT: set_convective_mixing'
!$omp end critical
            end if
            
            Zeta = Gamma*Gamma*Gamma/Bcubed  ! C&G 14.80
            d_Zeta_dvb = (3d0*Gamma*Gamma*d_Gamma_dvb - d_Bcubed_dvb*Zeta)/Bcubed

            !Zeta = a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
            !d_Zeta_dvb = d_Gamma_dvb*(((Gamma**2+2*Gamma)*a0)/(pow4(Gamma)*a0**2&
            !   +(2*pow3(Gamma)+2*Gamma**2)*a0+Gamma**2+2*Gamma+1)) + &
            !   d_a0_dvb*(pow2(Gamma)/(Gamma*(Gamma*a0+1)+1)&
            !             -(pow4(Gamma)*a0)/pow2(Gamma*(Gamma*a0+1)+1))

            !write(*,*) "Compare Zeta", kz, log10_cr(T), Zeta, &
            !   a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
            
            ! Zeta must be >= 0 and < 1
            if (Zeta < 0d0) then
               Zeta = 0
               d_Zeta_dvb = 0
            else if (Zeta >= 1d0) then
               Zeta = 1d0
               d_Zeta_dvb = 0
            end if
            
            !gradT = (1d0 - Zeta)*gradr + Zeta*gradL ! C&G 14.79 with gradL for grada
            !d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_gradL_dvb + &
            !            (gradL - gradr)*d_Zeta_dvb
            gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79
            d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_grada_dvb + &
                        (grada - gradr)*d_Zeta_dvb

            if (is_bad(gradT)) then
               call set_no_mixing
               quit = .true.
               return
            end if
         
         end subroutine set_convective_mixing   

         subroutine set_semiconvection ! Langer 1983 & 1985
            real(dp) :: alpha, bc, LG, &
               a0, a1, a2, a3, a4, a5, a6, a, &
               b1, b2, b3, b4, b5, b6, b7, b, div, bsq
            real(dp), dimension(nvbs) :: &
               d_bc_dvb, d_LG_dvb, d_a0_dvb, d_a1_dvb, d_a2_dvb, d_a3_dvb, d_a4_dvb, &
               d_a5_dvb, d_a6_dvb, d_a_dvb, d_b1_dvb, d_b2_dvb, d_b3_dvb, d_b4_dvb, &
               d_b5_dvb, d_b6_dvb, d_b7_dvb, d_b_dvb, d_div_dvb
            
            include 'formats.dek'
            if (dbg) write(*,*) 'check for semiconvection'
            call set_no_mixing ! sets gradT = gradr
            D_semi = alpha_semiconvection*radiative_conductivity/(6*Cp*rho) &
                  *(gradr - grada)/(gradL - gradr)
            if (D_semi <= 0) then
               if (dbg) then
                  write(*,1) 'set_no_mixing D_semi', D_semi
                  write(*,1) 'alpha_semiconvection', alpha_semiconvection
                  write(*,1) 'radiative_conductivity', radiative_conductivity
                  write(*,1) 'gradr - grada', gradr - grada
                  write(*,1) 'gradL - gradr', gradL - gradr
                  stop
               end if
               call set_no_mixing
               return
            end if
            d_D_semi_dvb = 0 ! not used, so skip for now.
            conv_vel = 3*D_semi/Lambda 
            d_conv_vel_dvb = 0
            call time_limit_conv_vel
            if (conv_vel > max_conv_vel) conv_vel = max_conv_vel
            if (init_conv_vel /= conv_vel) D_semi = conv_vel*Lambda/3
            if (D_semi <= 0) then
               call set_no_mixing
               return
            end if
            
            mixing_type = semiconvective_mixing
            if (dbg) write(*,2) 'mixing_type', mixing_type
            
            if (semiconvection_option == 'Langer_85 mixing; gradT = gradr') return
            if (semiconvection_option /= 'Langer_85') then
               write(*,*) 'MLT: unknown values for semiconvection_option ' // &
                  trim(semiconvection_option)
               ierr = -1
               return
            end if
            
            
!            Solve[{
!                  L/Lrad - Lsc/Lrad - 1 == 0, 
!                  Lrad == grad LG, 
!                  gradMu == (4 - 3*beta)/beta*gradL_composition_term,
!                  Lsc/Lrad == alpha (grad - gradA)/(2 grad (gradL - grad))
!                              (grad - gradA - (beta (8 - 3 beta))/bc gradMu)}, 
!                  grad, {Lsc, Lrad, gradMu}] // Simplify
                  
            alpha = min(1d0, alpha_semiconvection)

            bc = 32 - 24*beta - beta*beta
            d_bc_dvb = - 24*d_beta_dvb - 2*d_beta_dvb*beta
            
            LG = (16d0/3d0*pi*clight*m*cgrav*crad*T*T*T*T)/(P*opacity)
            d_LG_dvb = LG*(4d0*dT_dvb/T - dP_dvb/P - d_opacity_dvb/opacity)
            
            a0 = alpha*gradL_composition_term*LG
            d_a0_dvb = alpha*gradL_composition_term*d_LG_dvb
            
            a1 = -2*bc*L
            d_a1_dvb = -2*L*d_bc_dvb
            d_a1_dvb(mlt_dL) = d_a1_dvb(mlt_dL) - 2*bc
            
            a2 = 2*alpha*bc*grada*LG
            d_a2_dvb = 2*alpha*(d_bc_dvb*grada*LG + bc*d_grada_dvb*LG + bc*grada*d_LG_dvb)
            
            a3 = -2*bc*gradL*LG
            d_a3_dvb = -2*(d_bc_dvb*gradL*LG + bc*d_gradL_dvb*LG + bc*gradL*d_LG_dvb)
            
            a4 = 32*a0
            d_a4_dvb = 32*d_a0_dvb
            
            a5 = -36*beta*a0
            d_a5_dvb = -36*(d_beta_dvb*a0 + beta*d_a0_dvb)
            
            a6 = 9*beta*beta*a0
            d_a6_dvb = 9*(2*beta*d_beta_dvb*a0 + beta*beta*d_a0_dvb)
            
            a = a1 + a2 + a3 + a4 + a5 + a6
            d_a_dvb = d_a1_dvb + d_a2_dvb + d_a3_dvb + d_a4_dvb + d_a5_dvb + d_a6_dvb 
                           
            b1 = 32 - 36*beta + 9*beta*beta
            d_b1_dvb = - 36*d_beta_dvb + 18*beta*d_beta_dvb
            
            b2 = b1*a0
            d_b2_dvb = d_b1_dvb*a0 + b1*d_a0_dvb
            
            b3 = -2*gradL*L + alpha*grada*grada*LG
            d_b3_dvb = -2*d_gradL_dvb*L + alpha*(2*grada*d_grada_dvb*LG + grada*grada*d_LG_dvb)
            d_b3_dvb(mlt_dL) = d_b3_dvb(mlt_dL) - 2*gradL
            
            b4 = (-alpha*gradA + gradL)*LG
            d_b4_dvb = (-alpha*d_grada_dvb + d_gradL_dvb)*LG + (-alpha*gradA + gradL)*d_LG_dvb
            
            b5 = -b2 + 2*bc*(L + b4)
            d_b5_dvb = -d_b2_dvb + 2*d_bc_dvb*(L + b4) + 2*bc*d_b4_dvb
            d_b5_dvb(mlt_dL) = d_b5_dvb(mlt_dL) + 2*bc
            
            b6 = b2*grada + bc*b3
            d_b6_dvb = d_b2_dvb*grada + b2*d_grada_dvb + d_bc_dvb*b3 + bc*d_b3_dvb
            
            b7 = -4*(-2 + alpha)*bc*LG*b6
            d_b7_dvb = -4*(-2 + alpha)*(d_bc_dvb*LG*b6 + bc*d_LG_dvb*b6 + bc*LG*d_b6_dvb)
            
            b = b7 + b5*b5
            d_b_dvb = d_b7_dvb + 2*b5*d_b5_dvb
            
            div = 2*(-2 + alpha)*bc*LG
            d_div_dvb = 2*(-2 + alpha)*(d_bc_dvb*LG + bc*d_LG_dvb)

            bsq = sqrt(b)
            gradT = (a + bsq)/div
            d_gradT_dvb = -gradT*d_div_dvb/div + d_a_dvb/div + 0.5d0*d_b_dvb/(div*bsq)
            
         end subroutine set_semiconvection
                  
         subroutine set_no_mixing
            ! assumes have set gradr, scale_height, gradL, and Lambda.
            mixing_type = no_mixing
            gradT = gradr
            d_gradT_dvb = d_gradr_dvb
            conv_vel = 0
            d_conv_vel_dvb = 0
            D = 0
            d_D_dvb = 0
            D_semi = 0
            d_D_semi_dvb = 0
            D_thrm = 0
            d_D_thrm_dvb = 0
            conv_P = 0
            d_conv_P_dvb = 0
            Gamma = 0
            d_Gamma_dvb = 0
            unsmoothed_gradT = 0
            d_unsmoothed_gradT_dvb = 0
         end subroutine set_no_mixing         
         
         subroutine show_args
 1          format(a30,1pe26.16)
            
            write(*,1) 'cgrav = ', cgrav
            write(*,1) 'm = ', m
            write(*,1) 'r = ', r 
            write(*,1) 'T = ', T 
            write(*,1) 'Rho = ', Rho 
            write(*,1) 'L  = ', L 
            write(*,1) 'P = ', P
            write(*,1) 'chiRho = ', chiRho 
            write(*,1) 'chiT = ', chiT
            write(*,1) 'Cp = ', Cp 
            write(*,1) 'xh = ', xh
            write(*,1) 'opacity = ', opacity 
            write(*,1) 'grada = ', grada
            write(*,1) 'mixing_length_alpha = ', mixing_length_alpha
            
         end subroutine show_args


         subroutine revise_using_cv_var_variable
            ! change D, d_D_dvb, gradT, d_gradT_dvb
            include 'formats.dek'
            real(dp) ff1, ff2, ff3, ff4, sqrt_x, tmp
            real(dp) :: cv_var, A_0, A_1, A_2, A_numerator, A_denom, &
               inv_sqrt_x, save_gradT
            real(dp), dimension(nvbs) :: &
               dA_0_dvb, dA_1_dvb, dA_2_dvb, dA_numerator_dvb, dA_denom_dvb, &
               d_inv_sqrt_x_dvb, d_cv_var_dvb, d_save_dvb
            real(qp) :: q1, q2, q3
            real(qp), dimension(nvbs) :: dq1_dvb, dq2_dvb, dq3_dvb
            integer :: j

            quit = .false.
            if (kz == 0) return
            if (ss% use_321x_dconv_vel_dt_equation .and. &
                ss% lnT_start(kz)/ln10 < ss% min_logT_for_321x_dconv_vel_dt) then
               return
            end if
            
            if (ss% use_321x_dconv_vel_dt_equation) then
               cv_var = ss% conv_vel(kz)
               d_cv_var_dvb = 0
               d_cv_var_dvb(mlt_cv_var) = 1d0
            else
               !!Pablo: TODO, not sure if this helps, use velocity from middle of the step
               !!NOTE, needed to change hydro_vars as well to make conv_vel_start available
               cv_var = 0.5d0*(ss% conv_vel(kz)+ss% conv_vel_start(kz))
               d_cv_var_dvb = 0
               d_cv_var_dvb(mlt_cv_var) = 0.5d0*(ss% conv_vel(kz)+ss% conv_vel_v0)
               !cv_var = ss% conv_vel(kz)
               !d_cv_var_dvb = 0
               !d_cv_var_dvb(mlt_cv_var) = ss% conv_vel(kz)+ss% conv_vel_v0
            end if
            
            d_D_dvb(mlt_cv_var) = 0d0
            d_gradT_dvb(mlt_cv_var) = 0d0

            D = cv_var*Lambda/3     ! diffusion coefficient [cm^2/sec]
            if (debug) write(*,1) 'D', D
            d_D_dvb = (d_cv_var_dvb*Lambda + cv_var*d_Lambda_dvb)/3

            x = Q*Rho / (2d0*P) ! using x as a temporary variable here
            d_x_dvb = x*(drho_dvb/rho + dQ_dvb/Q - dP_dvb/P)
         
            convective_conductivity = Cp*grav*Lambda*Lambda*Rho*(sqrt(x)) / 9 ! erg / (K cm sec)
            
            if (convective_conductivity < 0) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,1) 'MLT error: convective_conductivity', convective_conductivity
                  write(*,1) 'Cp', Cp
                  write(*,1) 'grav', grav
                  write(*,1) 'Lambda', Lambda
                  write(*,1) 'Rho', Rho
                  write(*,1) 'x', x
                  stop 'mlt'
               else
                  return
               end if
            end if
            
            if (debug) write(*,1) 'convective_conductivity', convective_conductivity
            d_cc_dvb = convective_conductivity* &
                 (d_Cp_dvb/Cp + d_grav_dvb/grav + &
                     2*d_Lambda_dvb/Lambda + dRho_dvb/rho + d_x_dvb/(2*x))

            if (MLT_option == 'Cox') then ! this assumes optically thick

               a0 = 9d0/4d0
               d_a0_dvb = 0d0

               ! 'A' param is ratio of convective to radiative conductivities   C&G 14.98
               A = convective_conductivity / radiative_conductivity !  unitless.

               if (debug) write(*,1) 'A', A
               dA_dvb = (d_cc_dvb - d_rc_dvb*A) / radiative_conductivity
               
               if (A < 0 .or. is_bad(A)) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,*) "MLT_option == 'Cox'", MLT_option == 'Cox'
                     write(*,1) 'A', A
                     write(*,1) 'convective_conductivity', convective_conductivity
                     write(*,1) 'radiative_conductivity', radiative_conductivity
                     stop 'mlt'
                  else
                     return
                  end if
               end if

            else
            
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1d0/Henyey_nu_param
                  ff2=0.5d0
                  ! popular values for y are 1/3 or 3/(4*pi**2)
                  ff3=8d0/Henyey_y_param
                  ff4=1d0/Henyey_y_param
               case ('ML1')
                  ff1=1d0/8d0
                  ff2=1d0/2d0
                  ff3=24d0
                  ff4=0d0
               case ('ML2')
                  ff1=1d0
                  ff2=2d0
                  ff3=16d0
                  ff4=0d0
               case ('Mihalas')
                  ff1=1d0/8d0
                  ff2=1d0/2d0
                  ff3=16d0
                  ff4=2d0
               case default
                  write(*,'(3a)') 'Error: ',trim(MLT_option), &
                     ' is not an allowed MLT version for convection'
                  write(*,*)
                  return
               end select
            
               omega = Lambda*Rho*opacity !dimensionless
               d_omega_dvb = omega*( d_Lambda_dvb/Lambda + dRho_dvb/Rho + d_opacity_dvb/opacity)

               ! the variable theta in no longer needed
               ! theta = omega / ( 1d0 + Henyey_y_param*omega**2 )
               ! d_theta_dvb = d_omega_dvb*(1d0 - Henyey_y_param*omega**2 ) /
               ! ( ( 1d0 + Henyey_y_param*omega**2 )**2 )

               ! a0 = 0.75d0*omega*theta
               !d_a0_dvb = a0*( d_omega_dvb/omega + d_theta_dvb/theta )
               a0 = (3d0/16d0)*ff2*ff3/(1d0+ff4/(omega*omega))
               d_a0_dvb = a0*2d0*ff4*d_omega_dvb/(ff4 + omega*omega)/omega

               ! A = sqrt(P*Q*rho/Henyey_nu_param)*(Cp*mixing_length_alpha)/
               !        (2*crad*clight*T**3*theta)               
               A_0 = sqrt(ff1*P*Q*rho)
               dA_0_dvb = ff1*(dP_dvb*Q*rho + P*dQ_dvb*rho + P*Q*drho_dvb)/(2*A_0)
               
               A_1 = 4d0*A_0*Cp
               dA_1_dvb = 4d0*(dA_0_dvb*Cp + A_0*d_Cp_dvb)
               
               A_2 = mixing_length_alpha*omega*(1d0+ff4/(omega*omega))
               dA_2_dvb = mixing_length_alpha*(1d0-ff4/(omega*omega))*d_omega_dvb

               if (is_bad(dA_2_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_2_dvb(mlt_dlnT00)', dA_2_dvb(mlt_dlnT00)
                     stop 'MLT'
                  else
                     return
                  end if
               end if
               
               A_numerator = A_1*A_2
               dA_numerator_dvb = dA_1_dvb*A_2 + A_1*dA_2_dvb

               if (is_bad(dA_numerator_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_numerator_dvb(mlt_dlnT00)', dA_numerator_dvb(mlt_dlnT00)
                     stop 'MLT'
                  else
                     return
                  end if
               end if
                     
               A_denom = ff3*crad*clight*T*T*T
               dA_denom_dvb = A_denom*3d0*dT_dvb/T
               
               A = A_numerator/A_denom                     
               dA_dvb = dA_numerator_dvb/A_denom - A_numerator*dA_denom_dvb/(A_denom*A_denom)

               if (is_bad(dA_dvb(mlt_dlnT00))) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'dA_dvb(mlt_dlnT00)', dA_dvb(mlt_dlnT00)
                     stop 'MLT'
                  else
                     return
                  end if
               end if
               
               if (A < 0) then
                  ierr = -1
                  if (ss% report_ierr) then
                     write(*,1) 'A', A
                     write(*,1) 'A_numerator', A_numerator
                     write(*,1) 'A_denom', A_denom
                     write(*,1) 'A_1', A_1
                     write(*,1) 'ff3', ff3
                     write(*,1) 'A_0', A_0
                     stop 'mlt'
                  else
                     return
                  end if
               end if
            
            end if

            ! average convection velocity, vbar   C&G 14.86b
            ! cv_var = vsound*Sqrt(Q)*alpha*Gamma / (2*Sqrt(2*Gamma1)*A)
            ! vsound = Sqrt(Gamma1*P / rho), so
            ! cv_var = Sqrt(Q*P / (8*rho))*alpha*Gamma / A
            ! cv_var = sqrt_x*alpha*Gamma / A
            ! Gamma = cv_var*A/(sqrt_x*alpha)

            x = Q*P / (8d0*rho) ! using x as a temporary variable here
            sqrt_x = sqrt(x)
            
            if (A <= 1d-99 .or. sqrt_x <= 1d-99) then
               Gamma = 1d25
               d_Gamma_dvb = 0
               !if (dbg) &
                  write(*,2) 'mlt: A or sqrt_x too small', kz, A, sqrt_x
            else      
               if (dbg) write(*,*) 'calculate Gamma using cv_var'
               inv_sqrt_x = 1d0/sqrt_x
               d_inv_sqrt_x_dvb = &
                  (Q*P*drho_dvb/rho - Q*dP_dvb - P*dQ_dvb)/(16d0*x*sqrt_x*rho)
               if (.false.) then
                  write(*,2) 'rel_diff new old Gamma', kz, &
                     (cv_var*A*inv_sqrt_x/mixing_length_alpha - Gamma)/Gamma, &
                     cv_var*A*inv_sqrt_x/mixing_length_alpha, Gamma
                  write(*,2) 'old d_Gamma_dvb(mlt_dL)', kz, d_Gamma_dvb(mlt_dL)
                  write(*,2) 'old d_Gamma_dvb(mlt_dlnR)', kz, d_Gamma_dvb(mlt_dlnR)
               end if
               Gamma = cv_var*A*inv_sqrt_x/mixing_length_alpha  
               d_Gamma_dvb = ( &
                  d_cv_var_dvb*A*inv_sqrt_x + &
                  cv_var*dA_dvb*inv_sqrt_x + &
                  cv_var*A*d_inv_sqrt_x_dvb)/mixing_length_alpha  
               if (.false.) then
                  write(*,2) 'new d_Gamma_dvb(mlt_dL)', kz, d_Gamma_dvb(mlt_dL)
                  write(*,2) 'new d_Gamma_dvb(mlt_dlnR)', kz, d_Gamma_dvb(mlt_dlnR)
               end if
            end if

            
            if (.true.) then
            
               ! C&G, eq. (14.78), but rewritten in a way that prevents
               ! the multiplication of terms that go as ~Gamma. This is because
               ! Gamma can be very large, and Gamma^2 can actually overflow and
               ! produce a NaN.
               tmp = 1d0/(1d0 + a0*Gamma - a0*Gamma/(Gamma+1d0))
               Zeta = 1d0 - tmp
               d_Zeta_dvb =d_Gamma_dvb*(((a0*Gamma/(Gamma+1d0))*tmp)&
                  *((Gamma+2d0)/(Gamma+1d0)*tmp)) + &
                  d_a0_dvb*((Gamma/(Gamma+1))*(Gamma*tmp)*tmp)
            
            else if (.false.) then ! quad precision
            
               q1 = Gamma
               dq1_dvb = d_Gamma_dvb
               q2 = a0
               q3 = q2*q1*q1/(1+q1*(1+q2*q1)) ! C&G, eq. (14.78)
               dq3_dvb = dq1_dvb*(((q1*q1+2*q1)*q2)/(q1*q1*q1*q1*q2*q2&
                  +(2*q1*q1*q1+2*q1*q1)*q2+q1*q1+2*q1+1))
               Zeta = q3
               d_Zeta_dvb = dq3_dvb
               
            else

               Zeta = a0*Gamma*Gamma/(1+Gamma*(1+a0*Gamma)) ! C&G, eq. (14.78)
               d_Zeta_dvb = d_Gamma_dvb*(((Gamma*Gamma+2*Gamma)*a0)/(pow4(Gamma)*a0*a0&
                  +(2*pow3(Gamma)+2*Gamma*Gamma)*a0+Gamma*Gamma+2*Gamma+1))

            end if

            ! TESTING
            ! 'B' param  C&G 14.81
            !Bcubed = (A*A / a0)*diff_grads         
            !d_Bcubed_dvb = (A*A / a0)*d_diff_grads_dvb + (2*A*dA_dvb / a0)*diff_grads
            !Zeta = Gamma*Gamma*Gamma/Bcubed  ! C&G 14.80
            !d_Zeta_dvb = (3d0*Gamma*Gamma*d_Gamma_dvb - d_Bcubed_dvb*Zeta)/Bcubed
            
            !if (kz == 33) then
            if (.false.) then
               write(*,2) 'cv_var', kz, cv_var
               write(*,2) 'ss% mlt_vc(kz)', kz, ss% mlt_vc(kz)
               write(*,2) 'standard mlt gradT', kz, gradT
               write(*,2) 'new gradT', kz, (1d0 - Zeta)*gradr + Zeta*grada
               write(*,2) 'new Zeta', kz, Zeta
               write(*,2) 'ss% gradr(kz)', kz, ss% gradr(kz)
               write(*,2) 'gradr', kz, gradr
               write(*,2) 'ss% grada(kz)', kz, ss% grada(kz)
               write(*,2) 'grada', kz, grada
               write(*,2) 'ss% L(kz)', kz, ss% L(kz)
               write(*,2) 'ss% r(kz)', kz, ss% r(kz)
               write(*,*)
            end if
            
            save_gradT = gradT
            d_save_dvb = d_gradT_dvb        

            !gradT = (1d0 - Zeta)*gradr + Zeta*gradL ! C&G 14.79 with gradL for grada
            !d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_gradL_dvb + &
            !            (gradL - gradr)*d_Zeta_dvb
            gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79
            d_gradT_dvb = (1d0 - Zeta)*d_gradr_dvb + Zeta*d_grada_dvb + &
                        (grada - gradr)*d_Zeta_dvb
            
            if (.false. .and. Zeta > 0.45 .and. Zeta < 0.55) then
!$OMP critical
               write(*,2) 'Zeta', kz, Zeta
               write(*,2) 'rel_diff new old gradT', kz, (gradT - save_gradT)/gradT, gradT, save_gradT
               do j=1,num_mlt_partials
                  tmp = d_gradT_dvb(mlt_cv_var)*d_conv_vel_dvb(j) + d_gradT_dvb(j)
                  write(*,3) 'rel_diff new old gradT partial ' // trim(mlt_partial_str(j)), j, kz, &
                     (tmp-d_save_dvb(j))/tmp, tmp, d_save_dvb(j)
               end do
               stop 'revise_using_cv_var_variable'
!$OMP end critical
            end if
                        
            if (is_bad(d_gradT_dvb(mlt_cv_var))) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'd_gradT_dvb(mlt_cv_var)', kz, d_gradT_dvb(mlt_cv_var)
                  stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if

            if (is_bad(gradT)) then
               ierr = -1
               if (ss% report_ierr) then
                  write(*,2) 'gradT', kz, gradT
                  stop 'revise_using_cv_var_variable'
               else
                  return
               end if
            end if
         
         end subroutine revise_using_cv_var_variable


      end subroutine Get_results


      end module mlt_info
