module atmosphere

   use auto_diff
   use star_lib
   use star_def
   use const_def
   use math_lib
   use agn_stars
   use utils_lib
   use mass_changes

   implicit none

   public
   contains




   subroutine other_surface_PT(id, skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

         integer, intent(in) :: id
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

      logical, parameter :: dbg = .false.
      type (star_info), pointer :: s
      real(dp) :: minT, blend_minT, Ledd_rdp
      ! real(dp), parameter :: kM = 10.0_dp, kN = 0.5_dp, lambda_0 = 0.5_dp
      type(auto_diff_real_4var_order1) :: M, R, L, kappa
      type(auto_diff_real_4var_order1) :: R_ph, R_bondi, R_acc, vesc2, Mdot, loss, super_eddington_factor, L_incident
      type(auto_diff_real_4var_order1) :: H_disk, erf_arg, erf_f1, erf_f2, vortpar, R_AGN, R_Hill, tidepar
      type(auto_diff_real_4var_order1) :: L_stream, L_shock, Ledd, T_surf, P_surf, P_rad, P_acc, P_outflow, P_simple
      type(auto_diff_real_4var_order1) :: simple_Teff, simple_lnT, simple_lnP, agn_Teff, agn_lnT, agn_lnP
      type(auto_diff_real_4var_order1) :: lnT, lnP, overall_Teff
               ! for chen's prescription
      type(auto_diff_real_4var_order1) :: Leff_Ledd_ratio, facA, Lambda_factor, mdot_rad, mdot_bh, mdot_k, ratio_kb, chen_x
      ! Basic setup
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Setup autodiff variables
      M = s%m(1)
      M%d1val1 = 1d0

      R = s%r(1)
      R%d1val2 = 1d0

      L = s%L(1)
      L%d1val3 = 1d0

      !kappa = const_kappa
      kappa = 0.2d0*(1d0+s%accretion_species_xa(1))
      kappa%d1val4 = 0d0

      ! Basic quantities
      R_bondi = 2d0*standard_cgrav * M / pow2(const_csb) ! Bondi radius
      R_acc = R_bondi
      Mdot = pi * pow2(R_bondi) * const_rho * const_csb ! Mdot gain

      R_ph = pow2(kappa * const_rho) * pow3(R_bondi) ! Photosphere radius
      R_ph = R_ph + R
      R_ph = min(R_ph, R_bondi)

      ! Adjust gain down due to radiation pressure
      Ledd = eval_Ledd(s, ierr)*(M/M%val)

      if (s% x_logical_ctrl(7)) then
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if

      if (s% x_logical_ctrl(1)) then
         if (L < Ledd) then
            Mdot = Mdot * pow2(1d0 - L/Ledd) ! Suggested by Alexander Dittmann
         else
            Mdot = 0d0
         end if
      else
         Mdot = Mdot * (1d0 - tanh(abs(L / Ledd)))
      end if

      ! ! Adjust gain down due to gap opening in the disk
      ! if (s% x_logical_ctrl(2)) then
      !    Mdot = Mdot * exp(-pow(M / Msun / cutoff_mass, 2d0/3d0)) ! Form suggested by Doug Lin
      ! end if

      ! Adjust gain down due to gap opening in the disk with smooth transition
      if (s% x_logical_ctrl(2)) then
         ! real(dp) :: cutoff_mass, cutoff_mass_low, cutoff_mass_high, transition_width
         ! cutoff_mass = s% x_ctrl(10)
         ! cutoff_mass_low = 0.95 * cutoff_mass  
         ! cutoff_mass_high = 1.05 * cutoff_mass
         ! transition_width = (cutoff_mass_high - cutoff_mass_low) / 2.0

         ! Smooth transition using tanh
         Mdot = Mdot * 0.5d0 * (1.0d0 - tanh((M / Msun - cutoff_mass) / (0.05d0 * cutoff_mass)))
      end if



      if (s% x_logical_ctrl(3)) then
         H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
         erf_arg = -R_bondi / H_disk
         erf_f1 = (-2d0/pow(pi, 0.5d0) )*pow( 1d0 - exp(-erf_arg*erf_arg), 0.5d0)
         erf_f2 = pow(pi, 0.5d0)*0.5d0 + (3.1d1/2d2)*exp(-erf_arg*erf_arg) - (3.41d2/8d3)*exp(-2d0*erf_arg*erf_arg)
         Mdot = 0.5d0 * Mdot * pow( pi, 0.5d0 ) * erf_f1 * erf_f2 / erf_arg
         ! Adjust accretion rate to take into account changes in disk density vertically.
         ! Comes from averaging exp(-(z/H)^2) over the surface of a sphere with r = r_Bondi
         ! ERF is approximated by a Burmann series using c1 = 31/200 and c2 = -341/8000
      end if

      if (s% x_logical_ctrl(4)) then
            R_AGN = pow( standard_cgrav * Mass_AGN * Msun / (Omega_AGN * Omega_AGN) , 1d0/3d0 )
            vortpar = (R_bondi * R_bondi * Omega_AGN)/(2d0 * const_csb * R_AGN) !avg vorticity parameter at R_bondi
            vortpar = (2d0 / (pi * vortpar))*asinh( pow(2d0*vortpar, 1d0/3d0))  !Krumholz, McKee, Klein 2005 eqn. A7
            Mdot = Mdot * pow(1d0 + pow(vortpar, -10d0), -0.1d0) !pick minimum factor
            ! Adjust accretion rate to take into account shear
      end if

      if (s% x_logical_ctrl(5)) then
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0)
            tidepar = pow(R_Hill / R_bondi, 2d0)
            Mdot = Mdot * min(1d0, tidepar) !pick minimum factor
            R_acc = min(R_Hill, R_Bondi)
            ! Use min(Bondi, Hill) for accretion rate. Similar to Rosenthal et al. 2020
      end if

      if (s% x_logical_ctrl(6)) then
            H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0)
            R_acc = min(R_Hill, R_Bondi)
            tidepar = min(1d0, R_Hill/R_bondi)
            tidepar = min(tidepar, H_disk/R_bondi)
            tidepar = tidepar * min(1d0, R_Hill/R_bondi)
            tidepar = tidepar * min(1d0, exp(2d0*(1d0-pow(R_Hill/H_disk, 2d0))))
            !Dobbs-Dixon, Li, Lin 2007, approx eqn. 28
            Mdot = Mdot * tidepar
            ! Adjust accretion rate to take into account the tidal barrier
      end if

      ! Calculate super-eddington mass loss
      vesc2 = 2d0*standard_cgrav*M / R
      super_eddington_factor = 0.1d0
      loss = super_eddington_mdot_autodiff(Ledd,M,L,R,super_eddington_factor)

      ! Set a minimum on Teff.
      minT = const_Tb

      ! Various luminosities
      L_shock = Mdot * M * standard_cgrav / R
      L_stream = L + L_shock
      L_incident = 4 * pi * pow2(R_ph) * boltz_sigma * pow4(minT)
      L_stream = L_stream + L_incident
      L_stream = max(L_stream, L_incident) ! Avoid NaN's due to bad near-surface grad(T)

      ! Evaluate simple photosphere
      simple_Teff = pow(L_stream / (4d0 * pi * pow2(R) * boltz_sigma), 0.25d0)
      simple_lnP = log(standard_cgrav * M / (pow2(R) * kappa))
      P_simple = exp(simple_lnP)
      simple_lnT = log(simple_Teff) ! First cell is tau=2/3 in this simple solution.

      ! Evaluate complicated photosphere
      agn_Teff = pow(L_stream / (4 * pi * boltz_sigma * pow2(R_ph)), 0.25d0)
      T_surf = agn_Teff * pow(R_ph / R, 5d0 / 8d0)
      agn_lnT = log(T_surf)

      ! P_surf
      P_acc = const_rho * pow2(const_csb) * pow(R_acc / R, 2.5d0)
      P_rad = crad * pow4(T_surf) / 3d0
      P_outflow = - loss * pow(vesc2, 0.5d0) / (4 * pi * pow2(R))
      P_surf = P_acc + P_rad + P_outflow + P_simple
      agn_lnP = log(P_surf)

      if (s% x_logical_ctrl(10)) then
         ! Mdot = pi * pow2(R_bondi) * const_rho * const_csb ! Mdot gain
         ! tidepar = * P_rad%val 
         ! mdot_bondi = pi * pow2(R_bondi) * const_rho * const_csb
         !mdot_rad = pi * const_csb * standard_cgrav * const_rho* s%m(1)/P_rad
         !mdot_rad = mdot_rad * (1d0 - abs(L / Ledd))/0.34d0 
         !mdot_bh = pi * s%r(1) * const_csb * (1d0 - abs(L / Ledd))/0.34d0 
         ! mdot_bh = mdot_bh 

        Leff_Ledd_ratio = abs(L / Ledd)
        if (Leff_Ledd_ratio >= 1.0d0) then
            chen_x = 0
        else  
            ! New Lambda_factor computation
            ! --- 2025‑05‑23 patch start ---
            facA = pow8(Leff_Ledd_ratio / 0.75_dp)            ! 阈值由 0.5 → 0.75
            ! 旧公式保留注释以便回溯
            ! Lambda_factor = pow2( abs(1.0_dp - (facA/(1.0_dp+facA)) * sqrt(Leff_Ledd_ratio)) ) / (1.0_dp + pow(Leff_Ledd_ratio,10))
            Lambda_factor = pow2(1.0_dp - (facA/(1.0_dp+facA)))   ! 简化公式，去除 / (1+(L/L_Edd)^10)
            ! --- 2025‑05‑23 patch end ---
             ! Leff_Ledd_ratio = min(Leff_Ledd_ratio, 1.0d0 - 1e-6_dp)
             !mdot_rad = (1.0_dp - Leff_Ledd_ratio) * Ledd * const_rho / (4.0_dp * 7.567601e-15_dp * 6.0e3_dp**4) !4_dp**4)
             ! Lambda_factor = pow2(1.0_dp - sqrt(Leff_Ledd_ratio)) / (1.0_dp + pow(Leff_Ledd_ratio,10) ) !((1.0_dp - Leff_Ledd_ratio**2d0)**0.5) / (1.0_dp + pow(Leff_Ledd_ratio,10d0))
             !mdot_rad = (1.0_dp - Leff_Ledd_ratio) * Ledd * (sqrt(2.0_dp) / (12.0_dp * const_csb**2))
             !mdot_bh = (1.0_dp -  Leff_Ledd_ratio) * Ledd * s%r(1) / (1e-6_dp + standard_cgrav * s%m(1))
             mdot_rad = Lambda_factor * Ledd * (sqrt(2.0_dp) / (12.0_dp * const_csb**2))
             mdot_bh = Lambda_factor * Ledd * s%r(1) / (1e-6_dp + standard_cgrav * s%m(1))
             mdot_k = 1d0 / (1d0 / mdot_bh + 1d0 / mdot_rad)
             ratio_kb = mdot_k / (1e-6_dp + pi * pow2(R_bondi) * const_rho * const_csb)
             chen_x = ratio_kb * (sqrt(4d0+ratio_kb) - sqrt(ratio_kb)) / (sqrt(4d0+ratio_kb) + sqrt(ratio_kb))
        end if 
         
                 ! Print variables influencing chen_x
!        ! write(*,*) '----------------------------------------'
!        ! write(*,*) '--- Debugging Chen Factor Calculation ---'
!        ! write(*,*) 'lambda', Leff_Ledd_ratio
!        ! write(*,*) 'mdot_rad:', mdot_rad
!        ! write(*,*) 'mdot_bh:', mdot_bh
!        ! write(*,*) 'mdot_k:', mdot_k
!        ! write(*,*) 'R_bondi:', R_bondi
!        ! write(*,*) 'const_rho:', const_rho
!        ! write(*,*) 'const_csb:', const_csb
!        ! write(*,*) 'ratio_kb', ratio_kb
!        ! write(*,*) 'chen_x:', chen_x
!        ! write(*,*) '----------------------------------------'
!        ! write(*,*) '----------atmosphere.f90----------------'
        
        loss = loss * (1.0_dp - Lambda_factor)
         Mdot = Mdot * chen_x
      end if
      
      
      if (s% x_logical_ctrl(12)) then 
          Mdot = Mdot * (1d0 - max(1e-2, pow(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit), 2d0)))
      end if
    

      ! blend gets set in extras_start_step
      lnT = agn_lnT
      lnP = agn_lnP
      overall_Teff = agn_Teff

      ! Format for output
      lnT_surf = lnT%val
      dlnT_dlnM = lnT%d1val1 * M%val
      dlnT_dlnR = lnT%d1val2 * R%val
      dlnT_dL = lnT%d1val3
      dlnT_dlnkap = lnT%d1val4 * kappa%val

      lnP_surf = lnP%val
      dlnP_dlnM = lnP%d1val1 * M%val
      dlnP_dlnR = lnP%d1val2 * R%val
      dlnP_dL = lnP%d1val3
      dlnP_dlnkap = lnP%d1val4 * kappa%val

      !Teff = overall_Teff%val

      !if (is_bad(lnT_surf) .or. is_bad(lnP_surf) .or. is_bad(Teff) .or. dbg) then
      if (is_bad(lnT_surf) .or. is_bad(lnP_surf) .or. dbg) then
!         write(*,*) 'Bad ATM'
!         write(*,*) '----------- Inputs -----------'
!         write(*,*) 'M',M
!         write(*,*) 'R',R
!         write(*,*) 'L',L
!         write(*,*) 'kap',kappa
!         write(*,*) 'Mdot', Mdot
!         write(*,*) '----------- Derived -----------'
!         write(*,*) 'Rph', R_ph
!         write(*,*) 'R_Bondi', R_bondi
!         write(*,*) 'L_stream', L_stream
!         write(*,*) 'Blend', atm_blend
!         write(*,*) '----------- Outputs -----------'
!         write(*,*) 'lnT',lnT
!         write(*,*) 'lnP',lnP
!         write(*,*) 'Teff',overall_Teff
      end if

      s% job% extras_rpar(3) = L%val
      s% job% extras_rpar(4) = R%val

      ! Save variables for history

      s% xtra1_array(4) = L_shock%val    ! Accretion Shock Luminosity
      s% xtra1_array(5) = L_stream%val   ! Stream Luminosity
      s% xtra1_array(6) = R_ph%val       ! Photospheric Radius
      s% xtra1_array(7) = R_bondi%val    ! R Bondi
      s% xtra1_array(10) = R%val         ! R Shock (R 'Surface' of MESA model)
      s% xtra1_array(11) = lnT%val       ! log2 T Shock (T 'Surface' of MESA model)

   end subroutine other_surface_PT



end module atmosphere
