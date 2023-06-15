! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

! Note: this file has been modified by Kayla Martin (incl. comments throughout), based on the work of Glenn-Michael Oomen (incl. comments and code).


      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

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


         ! Glenn edit:
         ! Use evolve_to_PAGB to just get the stellar wind when evolving towards
         ! the post-AGB phase and to remove the envelope once the core mass is reached.
         ! Use do_PAGB_evolution to evolve the post-AGB star with accretion with a specific
         ! chemical composition, while also including a stellar wind.

         s% other_adjust_mdot => evolve_to_PAGB         ! use to evolve to start of post-AGB phase and after 50,000 K to end of WD phase
         !s% other_adjust_mdot => do_PAGB_evolution      ! use to evolve WITH accretion applied


      end subroutine extras_controls



      subroutine evolve_to_PAGB(id, ierr)
         ! This routine evovles the star to the post-AGB with only a
         ! stellar wind applied.
         use star_def
         use star_lib, only: star_ptr
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: L1, M1, R1, T1, SCwind, AGBwind, CSPNwind, wind, &
                     x, g, z, logP
         real(dp), dimension(40) :: Tcond, abun
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! OOPS
            return
         end if      

         ! Compute the stellar wind for the (post-)AGB star.
         ! This follows the wind prescription from Miller Bertolami (2016):
         L1 = s% L_phot*Lsun
         M1 = s% mstar
         T1 = s% Teff
         R1 = sqrt(L1/(pi*crad*clight*T1*T1*T1*T1))
         g = 10**(s% photosphere_logg)    !g = 10**(s% log_surface_gravity)
         x = (log10(T1)-3.7)/0.2
         z = 1 - s% surface_h1 - s% surface_he4

         ! RGB wind from Schroder & Cuntz (2005)
         SCwind = 8d-14 * (T1/4000)**3.5 * (1 + 27478.2333/(g*4300.0))&
                   * ((L1/Lsun)*(R1/Rsun)/(M1/Msun))

         ! AGB wind based on pulsation period
         logP = -1.92 - 0.73 * log10(M1/Msun) + 1.86 * log10(R1/Rsun)
         if (logP > 2.0) then
            if (s% surface_c12 * 4 < s% surface_o16 * 3) then
                AGBwind = 10**(-9 + 0.0032*(10**logP))
                AGBwind = max(AGBwind, SCwind)
            else
                AGBwind = 10**(-16.54 + 4.08*logP)
                AGBwind = max(AGBwind, SCwind)
            end if
         else
            AGBwind = SCwind
            write(*,*) 'Using AGBwind = SCwind', AGBwind
         end if

         ! Hot central star wind
         CSPNwind = 9.778d-15 * (L1/Lsun)**1.674 * (z/0.02)**(-2.0/3.0)

         ! Transition from RGB wind to CSPN wind between logT of 3.7 and 3.9
         if (x <= 0.0) then
            wind = AGBwind
         else if (0.0 < x .and. x < 1.0) then
            wind = 10**((1 - x) * log10(AGBwind) + x * log10(CSPNwind))
         else if (1.0 <= x) then
            wind = CSPNwind
            write(*,*) 'Using wind = CSPNwind', wind
         end if

         ! Final mass change is the stellar wind.
         ! To remove the envelope of the AGB star,
         ! change this to -1d-4 in order to get rid of
         ! the envelope in a timescale of 10000 years.
         s% mass_change = -wind    ! use when evolving to desired He core mass, and to end of WD phase
         !s% mass_change = -1d-4    ! use to remove envelope rapidly down to 0.02 Msun (for Post-AGB star) in 10,000 years


      end subroutine evolve_to_PAGB


 
      subroutine do_PAGB_evolution(id, ierr)
         ! This subroutine is to evolve the star through the post-AGB
         ! with accretion. Here is where we set the initial accretion
         ! and disk parameters (i.e., disk mass, accretion rate,
         ! temperature to start accretion).
         use star_def
         use star_lib, only: star_ptr
         use chem_lib, only: chem_get_iso_id
         use chem_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: L1, M1, R1, T1, SCwind, AGBwind, CSPNwind, wind, x, g, logP, z, sum_frac, &
                    time_since_start_accr, Mdot0, Tstart, Md
         real(dp), dimension(40) :: Tcond, abun
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! OOPS
            return
         end if

         ! Give the parameters for the model -- convert solar mass to grams (1 solar mass = 1.989e33 grams).
         
         Mdot0 = 1d-7 * 1.989d33   ! initial accretion rate
         !Mdot0 = 5d-7 * 1.989d33
         !Mdot0 = 1d-6 * 1.989d33
         
         Md = 1d-2 * 1.989d33      ! initial disc mass
         !Md = 3d-2 * 1.989d33
         !Md = 7d-3 * 1.989d33
         
         Tstart = 4000d0           ! temperature at which accretion starts (change to imply the star leaves AGB interaction with varying radii)


         time_since_start_accr = 0
         ! Here, I will assign the time coordinate for accretion.
         ! Logical operator s% PAGB is set to .true. in star_data.inc and becomes .false. once accretion has started.
         ! s% accretion_time is defined as the starting time for accretion.
         if (.not. s% PAGB) then
            time_since_start_accr = (s% star_age - s% accretion_time)
            write(*,*) 'Star age: ', s% star_age, 'years'
            write(*,*) 'accretion_time: ', s% accretion_time, 'years'
            write(*,*) 'Time since start post-AGB phase:', time_since_start_accr, 'years'
         end if

         ! Compute the accretion rate (Eq. 5 in Oomen et al. 2019)
         ! Note that we define Mdot0 here as 2x that of Oomen et al. (2019)
         ! Logical operator s% deplete starts as .false., but once Teff of the post-AGB star
         ! becomes larger than the starting temperature, accretion will keep on going
         ! even if Teff drops below the starting temperature again as the result of accretion.
         ! Also, s% deplete will initiate the definition of the accretion abundances below.
         if (s% Teff > Tstart .or. s% deplete) then
            s% deplete = .true.
            s% accrete = 0.5 * Mdot0 * (1 + (2*Mdot0/Md)*time_since_start_accr)**(-1.5)
            
            write(*,*) 'Initial accretion rate Mdot0 = ', Mdot0/1.989d33, 'Msun/yr'
            write(*,*) 'Initial disk mass = ', Md/1.989d33, 'Msun'
            write(*,*) 'Starting temperature for accretion = ', Tstart, 'K'
            write(*,*) 'Accretion rate (s% accrete, from run_star_extras.f90) = ', s% accrete, 'grams/yr = ', s% accrete/1.989d33, 'Msun/yr'
         else
            s% deplete = .false.
            s% accrete = 0d0 * 1.989d33
            write(*,*) 'Using s% deplete = false (accretion rate = ', s% accrete/1.989d33, 'Msun/yr)'
         end if
         
         ! This is to remove superadiabaticity in the more massive post-AGB stars
         ! in order to improve convergence. We defined s% gradT_excess_min_logT because
         ! we want to only reduce the gradT excess at the bottom of the convective envelope,
         ! which is the problematic region for massive post-AGB stars (see Miller Bertolami 2016).
         ! It appears this function has been implemented in a more intelligent way in
         ! the newest MESA version by softening the T gradient at the boundaries of convective zones.
         s% gradT_excess_min_logT = 4.9

         ! Compute the stellar wind for the (post-)AGB star.
         ! This follows the wind prescription from Miller Bertolami (2016):
         L1 = s% L_phot*Lsun
         M1 = s% mstar
         T1 = s% Teff
         R1 = sqrt(L1/(pi*crad*clight*T1*T1*T1*T1))
         g = 10**(s% photosphere_logg)   !g = 10**(s% log_surface_gravity)
         x = (log10(T1)-3.7)/0.2
         z = 1 - s% surface_h1 - s% surface_he4

         ! RGB wind from Schroder & Cuntz (2005)
         SCwind = 8d-14 * (T1/4000)**3.5 * (1 + 27478.2333/(g*4300.0))&
                   * ((L1/Lsun)*(R1/Rsun)/(M1/Msun))

         ! AGB wind based on pulsation period
         logP = -1.92 - 0.73 * log10(M1/Msun) + 1.86 * log10(R1/Rsun)
         if (logP > 2.0) then
            if (s% surface_c12 * 4 < s% surface_o16 * 3) then
                AGBwind = 10**(-9 + 0.0032*(10**logP))
                AGBwind = max(AGBwind, SCwind)
            else
                AGBwind = 10**(-16.54 + 4.08*logP)
                AGBwind = max(AGBwind, SCwind)
            end if
         else
            AGBwind = SCwind
            write(*,*) 'Using AGBwind = SCwind', AGBwind
         end if

         ! Hot central star wind
         CSPNwind = 9.778d-15 * (L1/Lsun)**1.674 * (z/0.02)**(-2.0/3.0)

         ! Transition from RGB wind to CSPN wind between logT of 3.7 and 3.9
         if (x <= 0.0) then
            wind = AGBwind
         else if (0.0 < x .and. x < 1.0) then
            wind = 10**((1 - x) * log10(AGBwind) + x * log10(CSPNwind))
         else if (1.0 <= x) then
            wind = CSPNwind
            write(*,*) 'Using wind = CSPNwind', wind
         end if

         ! Final mass change is the stellar wind + accretion from disc
         s% mass_change = -wind + (s% accrete / Msun)
         !write(*,*) 'Wind = -', wind, 'Msun/yr'
         write(*,*) 'Final mass change (stellar wind + accretion) = ', s% mass_change, 'grams/yr'

         ! Here we define the chemical composition of the accreted material. This block
         ! should only be called once.
         ! I added some lines to assign the chemical composition of the accreted material
         ! s% xa(j,k) where j is species and k is coordinate.
         if (s% deplete .and. s% PAGB) then
            s% accrete_same_as_surface = .false.
            s% accrete_given_mass_fractions = .true.
            s% num_accretion_species = s% species

            ! Make the composition of the accreted gas in Fig. 2 in Oomen et al. (2019).
            ! Define the condensation temperatures:
            do i=1,20 ! elements until Na
                Tcond(i) = 0
            end do
            do i = 21, 23 ! Na
                Tcond(i) = 958
            end do
            do i = 24, 26 ! Mg
                Tcond(i) = 1336
            end do
            do i = 27, 29 ! Al
                Tcond(i) = 1653
            end do
            Tcond(30) = 1310 ! Si
            Tcond(31) = 664  ! S
            Tcond(32) = 1517 ! Ca
            Tcond(33) = 1659 ! Sc
            Tcond(34) = 1582 ! Ti
            Tcond(35) = 1296 ! Cr
            Tcond(36) = 1158 ! Mn
            Tcond(37) = 1334 ! Fe
            Tcond(38) = 1353 ! Ni
            Tcond(39) = 1037 ! Cu
            Tcond(40) = 726  ! Zn


            ! Get the accretion abundances (linear in logspace):
            do i = 1,s% species
                if (Tcond(i) <= 726) then      ! change 726 value (Tturnoff temperature) here and below, for different object depletion patterns
                    abun(i) = 1d0
                else
                    abun(i) = -4*(Tcond(i) - 726d0)/(1582d0-726d0)  ! scale wrt [Zn/Ti] = 4 (i.e., maximum depletion/scaling factor); change for different observed object [Zn/Ti] ratios
                    abun(i) = 10**abun(i)
                end if
            end do


            do i=1,s% species
                s% accretion_species_id(i) = chem_isos% name(s% chem_id(i))
                s% accretion_species_xa(i) = s% xa(i,1) * abun(i) ! takes the abundance at the surface and multiply to deplete
                write(*,*) s% accretion_species_id(i), 'has acrretion abundance', &
                                 s% accretion_species_xa(i), 'which was multiplied by', &
                                 log10(abun(i))
            end do

            ! These accretion fractions are more or less representative for Maas et al. (2005).
            sum_frac = sum(s% accretion_species_xa(1:s% species))
            write(*,*) 'The sum of the acrretion abundances is', sum_frac

            s% PAGB = .false. ! Set to false such that this section is not called anymore
            s% accretion_time = s% star_age ! Start of the post-AGB accretion timescale
        end if

        ! timestep control:
        ! We want to avoid that a too large fraction of the outer convective region becomes accreted in one timestep,
        ! since this would cause the outer layers to become entirely composed of accreted material. We impose that
        ! at most 10% of the outer convective zone becomes accreted in one timestep.
        if (s% deplete) then
            do i=1,s% species
                if (s% accretion_species_id(i) == 'ti48') then
                    if (s% xa(i,1) > 2 * s% accretion_species_xa(i)) then
                        if (s% conv_mx1_top == 1.0) then
                            s% max_years_for_timestep = 0.1 * ((s% conv_mx1_top - s% conv_mx1_bot) * s% star_mass) / s% mass_change
                        else
                            s% max_years_for_timestep = 0.1 * ((s% conv_mx2_top - s% conv_mx2_bot) * s% star_mass) / s% mass_change
                        end if
                    else
                        s% max_years_for_timestep = 0.0
                    end if
                end if
            end do
        end if



      end subroutine do_PAGB_evolution



! ----------------------------- Here starts the 'standard_run_star_extras.inc' stuff -----------------

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
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


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
      
