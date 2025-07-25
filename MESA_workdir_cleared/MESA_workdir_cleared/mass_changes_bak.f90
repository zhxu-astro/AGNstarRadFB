module mass_changes

   use star_lib
   use star_def
   use const_def
   use math_lib
   use agn_stars
   use auto_diff
   use rotation

   implicit none

   public
   contains




    real(dp) function super_eddington_mdot(Ledd, M, L, R, Eddington_factor) result(loss)
       use star_def
       ! Parameters for new Ledd_factor
       real(dp), intent(in) :: Ledd, M, L, R, Eddington_factor
       !real(dp), parameter :: kM = 10.0_dp, kN = 0.5_dp, lambda_0 = 0.5_dp
       real(dp) :: vesc2

       if (R <= 0.0d0 .or. M <= 0.0d0) then
!           write(*,*) "Error: Non-positive radius or mass in super_eddington_mdot."
           loss = 0.0d0
           return
       end if

       vesc2 = 2d0 * standard_cgrav * M / R

       if (vesc2 <= 0.0d0) then
!           write(*,*) "Error: Non-positive vesc2 in super_eddington_mdot."
           loss = 0.0d0
           return
       end if

       !loss = -(L / vesc2) * (1d0 + tanh((L - Ledd) / (Eddington_factor * Ledd)))
        loss = -(L / vesc2) * (1- pow2(1d0 - sqrt(Eddington_factor)) / (1d0 + pow(Eddington_factor,10d0) ) ) * 1
     
    end function super_eddington_mdot


   type(auto_diff_real_4var_order1) function super_eddington_mdot_autodiff(Ledd, M,L,R,Eddington_factor) result(loss)
     use star_def
     type(auto_diff_real_4var_order1), intent(in) :: M , L , R, Eddington_factor, Ledd
     type(auto_diff_real_4var_order1) :: vesc2
     vesc2 = 2d0*standard_cgrav*M / R
     loss = -(L / vesc2) * (1- pow2(1d0 - sqrt(Eddington_factor)) / (1d0 + pow(Eddington_factor,10d0) ) ) * 1
     
     !loss = -(L / vesc2) * (1d0 + tanh((L - Ledd) / (Eddington_factor * Ledd)))
   end function super_eddington_mdot_autodiff

   subroutine eval_super_eddington_wind(s, L, M, R, ierr)
      type (star_info), pointer :: s
      real(dp), intent(in) :: L, M, R
      integer, intent(out) :: ierr

      real(dp) :: Ledd, Leff, vesc2
      include 'formats'

      ierr = 0
      s% super_eddington_wind_mdot = 0
      if (s% super_eddington_scaling_factor <= 0) return

      Ledd = s% prev_Ledd
      if (s% x_logical_ctrl(7)) then
         write(*,*) '----------------------------------------'
         write(*,*) 'J change Leddington: Ledd, j/j_c', Ledd, s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit)
         write(*,*) '----------------------------------------'
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if

      Leff = L/s% super_eddington_wind_Ledd_factor
      if (Leff <= Ledd) return
      vesc2 = s% cgrav(1)*M/R  ! GM/R vs. 2GM/R ?
      s% super_eddington_wind_mdot = s% super_eddington_scaling_factor*(Leff - Ledd)/vesc2
!      write(*,'(a60,i12,1p2e12.4)') 'super eddington wind: lg_Mdot, L/Ledd', s% model_number, log10(s% super_eddington_wind_mdot/(Msun/secyer)), L/Ledd
       
   end subroutine eval_super_eddington_wind

   ! call as: set_xa(id, ierr, min(gain, abs(loss)) * s% dt, abs(loss)*s%dt)
   subroutine set_xa(id, ierr, replenish, loss)
      use chem_lib, only: chem_get_iso_id
      integer, intent(in) :: id
      real(dp), intent(in) :: replenish ! Assume this amount is what gets replenished
      real(dp), intent(in) :: loss ! Assume this amount is what gets lost
      type (star_info), pointer :: s
      real(dp), allocatable, dimension(:) :: xa
      real(dp) :: pre_burn_frac, frac_h, frac_he
      integer, intent(out) :: ierr
      integer :: j, i, species, cid
      include 'formats'

      ierr = 0
      call star_ptr(id, s, ierr)
      species = s% species ! Number or isotopes in nuclear network

      ! Loop to determine and initialize names of isotopes of accreted material
      allocate(xa(species))
      xa(1:species) = 0
      do j=1,s% num_accretion_species
         if (len_trim(s% accretion_species_id(j)) == 0) cycle
         cid = chem_get_iso_id(s% accretion_species_id(j))
         if (cid <= 0) cycle
         i = s% net_iso(cid)
         if (i == 0) cycle
         xa(i) = s% accretion_species_xa(j)
      end do

      ! Loss in this routine is abs(loss)*dt
      lost_mass = 0d0
      do i=1,s%nz
         if (s%q(i) > 1d0 - loss / s%m(1)) then ! Store composition of lost material
            do j=1,species
               lost_mass(j) =  lost_mass(j) + s%xa(j,i) * s%dm(i)
            end do
         end if

         if (s%q(i) > 1d0 - replenish / s%m(1)) then !
            do j=1,species
               s%xa(j,i) = xa(j) ! Turns mass fraction of material down to 1- replenish to reflect omposition of accreted material
            end do
         end if
      end do
   end subroutine set_xa   




   subroutine other_adjust_mdot(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr

      real(dp) :: Ledd, Ledd_over_M, Leff, vesc2, loss

      real(dp) :: R_bondi, H_disk, erf_arg, erf_f1, erf_f2, vortpar, tidepar, R_AGN, R_Hill 
      real(dp) :: gain ! g/s
      real(dp) :: previous_radius, prev_luminosity
      ! for chen's prescription
      real(dp) :: mdot_rad_old, mdot_rad, mdot_bh, mdot_k, ratio_kb, chen_x, Leff_Ledd_ratio, Ledd_factor, facA

      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)

      ! if (s%model_number < 100) return

      prev_luminosity = s% job% extras_rpar(3)
      previous_radius = s% job% extras_rpar(4)

      ! Calculate super-eddington mass loss
      Ledd = eval_Ledd(s, ierr)
      Leff = prev_luminosity
      vesc2 = 2d0*standard_cgrav*s% m(1) / previous_radius
      loss = super_eddington_mdot(Ledd, s% m(1),Leff,previous_radius,0.1d0)

      if (s% x_logical_ctrl(7)) then
         write(*,*) '----------------------------------------'
         write(*,*) 'J change Leddington: Ledd, j/j_c', Ledd, s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit)
         write(*,*) '----------------------------------------'
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if


      ! Calculate accretion

      ! Bondi radius
      R_bondi = 2d0 * standard_cgrav * s%m(1) / pow2(const_csb)

      ! Mdot
      gain = pi * pow2(R_bondi) * const_rho * const_csb
      

      if (s% x_logical_ctrl(1)) then
         if (Leff < Ledd) then
            gain = gain * (1d0 - Leff/Ledd)**2d0 ! Suggested by Alexander Dittmann
            ! "As long as the radiative flux obeys the same scaling law with distance as the force of gravity
            ! (a decent approximation for the optically thin case), the net effect can be interpreted as reducing
            ! the effective mass by a factor of (1-L/L_edd). This then reduces the Bondi radius, and then the Bondi accretion rate, accordingly.
         else
            gain = 0.0
        end if
      else
        gain = gain * (1d0 - tanh(abs(Leff / Ledd)))
      end if

      if (use_mass_cutoff) then
          ! gain = gain * exp(-pow(s%m(1) / Msun / cutoff_mass, 2d0/3d0))
!          ! write(*,*) 'cutoff factor', exp(-pow(s%m(1) / Msun / cutoff_mass, 2d0/3d0))

          ! real(dp) :: cutoff_mass

          ! Get the cutoff mass from control parameters
          ! cutoff_mass = s% x_ctrl(10)

          ! Update gain using the tanh function for a smooth transition
          ! Transition width is set to 5% of the cutoff mass
          gain = gain * 0.5d0 * (1.0d0 - tanh((s%m(1) / Msun - cutoff_mass) / (0.05d0 * cutoff_mass)))

          ! Output the cutoff factor for debugging purposes
!          write(*,*) 'cutoff factor', 0.5d0 * (1.0d0 - tanh((s%m(1) / Msun - cutoff_mass) / (0.05d0 * cutoff_mass)))
      end if


      if (s% x_logical_ctrl(3)) then
            H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
            erf_arg = -R_bondi / H_disk
            erf_f1 = (-2d0/pow(pi, 0.5d0) )*pow( 1d0 - exp(-erf_arg*erf_arg), 0.5d0)
            erf_f2 = pow(pi, 0.5d0)*0.5d0 + (3.1d1/2d2)*exp(-erf_arg*erf_arg) - (3.41d2/8d3)*exp(-2d0*erf_arg*erf_arg)
            gain = 0.5d0 * gain * pow( pi, 0.5d0 ) * erf_f1*erf_f2 / erf_arg
            ! Adjust accretion rate to take into account changes in disk density vertically. 
            ! Comes from averaging exp(-(z/H)^2) over the surface of a sphere with r = r_Bondi
            ! ERF is approximated by a Burmann series using c1 = 31/200 and c2 = -341/8000
      end if

      if (s% x_logical_ctrl(4)) then
            H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
            R_AGN = pow( standard_cgrav * Mass_AGN * Msun / (Omega_AGN * Omega_AGN) , 1d0/3d0 )
            vortpar = (R_bondi * R_bondi * Omega_AGN)/(2d0 * const_csb * R_AGN) !avg vorticity parameter at R_bondi
            vortpar = (2d0 / (pi * vortpar))*asinh( pow(2d0*vortpar, 1d0/3d0))  !Krumholz, McKee, Klein 2005 eqn. A7
            gain = gain * pow(1d0 + pow(vortpar, -10d0), -0.1d0) !pick minimum factor
            ! Adjust accretion rate to take into account shear
      end if

      if (s% x_logical_ctrl(5)) then
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0)
            tidepar = pow(R_Hill / R_bondi, 2d0)
            gain = gain * min(1d0 , tidepar) !pick minimum factor
            ! Use min(Bondi, Hill) for accretion rate. Similar to Rosenthal et al. 2020
      end if

      if (s% x_logical_ctrl(6)) then 
            H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN 
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0) 
            tidepar = min(1d0, R_Hill/R_bondi) 
            tidepar = min(tidepar, H_disk/R_bondi) 
            tidepar = tidepar * min(1d0, R_Hill/R_bondi) 
            tidepar = tidepar * min(1d0, exp(2d0*(1d0-pow(R_Hill/H_disk, 2d0))))
            !Dobbs-Dixon, Li, Lin 2007, approx eqn. 28
            gain = gain * tidepar
            ! Adjust accretion rate to take into account the tidal barrier
      end if

      if (s% x_logical_ctrl(10)) then 
         ! gain = mdot_bondi
         ! Cspeed = 2.99792e10
         ! cs = sqrt(T-inf*Rgas)
         ! G = standard_cgrav
         ! rho = const_rho
         ! lambda_star = abs(Leff / Ledd)
         ! kappacrit = 0.34_dp
         ! arad = 7.567601e-15
         ! Rgas = 8.314
         ! Prad = 4.0_dp * arad * (cs**2 / Rgas)**4 / 3.0_dp
         ! mdot_rad := pi * const_csb * (1.0_dp - lambda_star) * G * M%val / (kappacrit * Prad) * s%job%x_ctrl(1)
         ! P_rad = crad * pow4(T_surf) / 3d0
         
         ! mdot_bondi = 4 * pi * R_bondi**2 * const_csb * const_rho
         
        Leff_Ledd_ratio = abs(Leff / Ledd)
        
        if (Leff_Ledd_ratio >= 1d0) then
              chen_x = 0
         ! Leff_Ledd_ratio = min(Leff_Ledd_ratio, 1.0d0 - 1e-6_dp)
        else
            ! Compute A factor, assume power=8
            facA = pow(Leff_Ledd_ratio / 0.5_dp, 8)
            facA = ( facA ) / (1.0_dp + facA) 

            ! Compute new Ledd_factor
            Ledd_factor = abs(1.0_dp - facA * Leff_Ledd_ratio**0.5_dp)**(1.0_dp / 0.5_dp) / (1.0_dp + Leff_Ledd_ratio**10_dp)


            write(*,*) "Debug: Leff_Ledd_ratio =", Leff_Ledd_ratio
            write(*,*) "Debug: Ledd_factor =", Ledd_factor


             mdot_rad = Ledd_factor * Ledd * (sqrt(2.0_dp) / (12.0_dp * const_csb**2))
             
             ! mdot_bh = (1.0_dp - lambda_star) * Ledd * R / (G * M) 
             !mdot_bh = (1.0_dp - Leff_Ledd_ratio) * Ledd * s%r(1) / (1e-16_dp + standard_cgrav * s%m(1)) ! 1e-16 is avoid NaN for zero-denominator
             mdot_bh = Ledd_factor * Ledd * s%r(1) / (1e-16_dp + standard_cgrav * s%m(1)) ! 1e-16 is avoid NaN for zero-denominator
             mdot_k = 1.0_dp / (1.0_dp / mdot_bh + 1.0_dp / mdot_rad)
             ratio_kb = mdot_k / (1e-6_dp + 4 * pi * R_bondi**2 * const_csb * const_rho)
             chen_x = ratio_kb * (sqrt(4.0_dp + ratio_kb) - sqrt(ratio_kb)) / (sqrt(4.0_dp + ratio_kb ) + sqrt(ratio_kb) ) 
         end if 
         
        !write(*,*) '----------------------------------------'
        !write(*,*) '--- Debugging Chen Factor Calculation ---'
        !write(*,*) 'mdot_rad_old:', mdot_rad_old
        !write(*,*) 'mdot_rad:', mdot_rad
        !write(*,*) 'mdot_rad_diff:', mdot_rad_old -  mdot_rad
        !write(*,*) 'mdot_bh:', mdot_bh
        !write(*,*) 'mdot_k:', mdot_k
        write(*,*) 'r(1):', s%r(1)
        write(*,*) 'm(1):', s%m(1)
        !write(*,*) 'ratio_kb:', ratio_kb
        
        !write(*,*) 'chen_x:', chen_x
        !write(*,*) '----------------------------------------'
        !write(*,*) '----------mass_change.f90----------------'
        write(*,*) '----------------------------------------'
        write(*,*) "mdot_rad:", mdot_rad, "mdot_bh:", mdot_bh, 'mdot_k:', mdot_k, "chen_x:", chen_x
        write(*,*) '----------------------------------------'
        
         gain = gain * chen_x
      end if
      
      if (s% x_logical_ctrl(12)) then 
          gain = gain * (1d0 - max(1e-2, pow(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit), 2d0)))
      end if

      gain = gain * mdot_blend
      loss = loss * (1-Ledd_factor)
      s% mstar_dot = gain + loss

      last_gain = gain * s%dt

     ! Improved scheme for handling joint loss and gain
     if (loss < 0d0) then
     ! Lose mass. Then gain an equal amount of fresh material.
          call set_xa(id, ierr, min(gain, abs(loss)) * s% dt, abs(loss)*s%dt)
     end if


     call update_J_dist(id, gain, loss, R_bondi)

    s% job% extras_rpar(i_gain) = gain  ! We don't seem to use this, remove?
    s% job% extras_rpar(i_loss) = loss  ! We don't seem to use this, remove?
     ! Save variables for history

     s% xtra1_array(1) = gain
     s% xtra1_array(2) = loss
     s% xtra1_array(9) = Leff / Ledd
    

    ! In other_adjust_mdot
!    write(*,*) "Gain:", gain, "Loss:", loss, "mstar_dot:", s% mstar_dot

   end subroutine other_adjust_mdot



end module mass_changes
