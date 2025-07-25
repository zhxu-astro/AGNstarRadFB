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

   !--------------------------------------------------------------------
   ! Function: super_eddington_mdot
   ! Purpose : 计算超爱丁顿条件下的质量损失率
   ! Inputs  : Ledd  - 爱丁顿极限亮度
   !           M, L, R - 质量、亮度和半径
   !           Eddington_factor - 爱丁顿因子（控制参数）
   ! Output  : loss  - 质量损失率
   !--------------------------------------------------------------------
   real(dp) function super_eddington_mdot(Ledd, M, L, R, Eddington_factor) result(loss)
       use star_def
       real(dp), intent(in) :: Ledd, M, L, R, Eddington_factor
       real(dp) :: vesc2

       if (R <= 0.0d0 .or. M <= 0.0d0) then
           loss = 0.0d0
           return
       end if

       vesc2 = 2d0 * standard_cgrav * M / R
       if (vesc2 <= 0.0d0) then
           loss = 0.0d0
           return
       end if

       ! 采用一种修正形式计算质量损失率
       loss = -(L / vesc2) * ( 1d0 - pow2(1d0 - sqrt(Eddington_factor)) / (1d0 + pow(Eddington_factor,10d0)) )
   end function super_eddington_mdot

   !--------------------------------------------------------------------
   ! Auto-differentiation 版本的超爱丁顿质量损失函数
   !--------------------------------------------------------------------
   type(auto_diff_real_4var_order1) function super_eddington_mdot_autodiff(Ledd, M, L, R, Eddington_factor) result(loss)
       use star_def
       type(auto_diff_real_4var_order1), intent(in) :: Ledd, M, L, R, Eddington_factor
       type(auto_diff_real_4var_order1) :: vesc2

       vesc2 = 2d0 * standard_cgrav * M / R
       loss = -(L / vesc2) * ( 1d0 - pow2(1d0 - sqrt(Eddington_factor)) / (1d0 + pow(Eddington_factor,10d0)) )
   end function super_eddington_mdot_autodiff

   !--------------------------------------------------------------------
   ! Subroutine: eval_super_eddington_wind
   ! Purpose   : 计算星体超爱丁顿风的质量损失率
   ! 注意：如果开启控制标志7 (use_modified_ledd) 则对 Ledd 进行修正
   !--------------------------------------------------------------------
   subroutine eval_super_eddington_wind(s, L, M, R, ierr)
      type(star_info), pointer :: s
      real(dp), intent(in) :: L, M, R
      integer, intent(out) :: ierr
      real(dp) :: Ledd, Leff, vesc2
      include 'formats'

      ierr = 0
      s%super_eddington_wind_mdot = 0.0d0
      if (s%super_eddington_scaling_factor <= 0.0d0) return

      ! 初始 Ledd 取自先前记录
      Ledd = s%prev_Ledd

      ! 【模式1】若开启控制标志7，则根据角动量比例修正 Ledd
      if (s%x_logical_ctrl(7)) then
         write(*,*) '----------------------------------------'
         write(*,*) 'J change Leddington: Ledd, j/j_c', Ledd, s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit)
         write(*,*) '----------------------------------------'
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if

      ! 有效亮度（考虑风的系数）
      Leff = L / s%super_eddington_wind_Ledd_factor
      if (Leff <= Ledd) return

      ! 注意：此处采用 s%cgrav(1) 计算逃逸速度平方（GM 或 2GM，根据具体定义）
      vesc2 = s%cgrav(1) * M / R
      s%super_eddington_wind_mdot = s%super_eddington_scaling_factor * (Leff - Ledd) / vesc2
   end subroutine eval_super_eddington_wind

   !--------------------------------------------------------------------
   ! Subroutine: set_xa
   ! Purpose   : 调整质量变化时的化学成分（补充与损失）
   !--------------------------------------------------------------------
   subroutine set_xa(id, ierr, replenish, loss)
      use chem_lib, only: chem_get_iso_id
      integer, intent(in) :: id
      real(dp), intent(in) :: replenish   ! 补充的质量量
      real(dp), intent(in) :: loss         ! 损失的质量量
      type(star_info), pointer :: s
      real(dp), allocatable, dimension(:) :: xa
      integer, intent(out) :: ierr
      integer :: j, i, species, cid
      include 'formats'

      ierr = 0
      call star_ptr(id, s, ierr)
      species = s%species

      allocate(xa(species))
      xa = 0.0d0
      do j = 1, s%num_accretion_species
         if (len_trim(s%accretion_species_id(j)) == 0) cycle
         cid = chem_get_iso_id(s%accretion_species_id(j))
         if (cid <= 0) cycle
         i = s%net_iso(cid)
         if (i == 0) cycle
         xa(i) = s%accretion_species_xa(j)
      end do

      ! 对每个 zone 处理质量损失与补充
      lost_mass = 0.0d0
      do i = 1, s%nz
         if (s%q(i) > 1d0 - loss / s%m(1)) then
            do j = 1, species
               lost_mass(j) = lost_mass(j) + s%xa(j,i) * s%dm(i)
            end do
         end if
         if (s%q(i) > 1d0 - replenish / s%m(1)) then
            do j = 1, species
               s%xa(j,i) = xa(j)
            end do
         end if
      end do
   end subroutine set_xa

   !--------------------------------------------------------------------
   ! Subroutine: other_adjust_mdot
   ! Purpose   : 根据不同方案调整星体的增减质量率
   !
   ! 这里主要区分两种模式：
   !  A. 【模式1+模式10】如果开启控制标志7，则先用角动量比例修正 Ledd，
   !     随后如果控制标志10开启，采用 Chen 方案计算 mdot（利用修正后的 Ledd）
   !
   !  B. 默认直接采用 (1 - Leff/Ledd) 或 tanh 修正计算 mdot（类似 1 - J/Jc 的效果）
   !--------------------------------------------------------------------
   subroutine other_adjust_mdot(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr

      real(dp) :: Ledd, Leff, vesc2, loss
      real(dp) :: R_bondi, H_disk, erf_arg, erf_f1, erf_f2, vortpar, tidepar, R_AGN, R_Hill
      real(dp) :: gain  ! 以 g/s 为单位的增质量率
      real(dp) :: previous_radius, prev_luminosity
      ! 用于 Chen 方案的变量
      real(dp) :: mdot_rad, mdot_bh, mdot_k, ratio_kb, chen_x, Leff_Ledd_ratio, Ledd_factor, facA

      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)

      prev_luminosity = s%job%extras_rpar(3)
      previous_radius = s%job%extras_rpar(4)

      ! 计算超爱丁顿条件下的质量损失（loss）
      Ledd = eval_Ledd(s, ierr)
      ! 初始化 Ledd_factor 为 1（默认不修正）
      Ledd_factor = 1.0d0
      Leff = prev_luminosity
      vesc2 = 2d0 * standard_cgrav * s%m(1) / previous_radius
      loss = super_eddington_mdot(Ledd, s%m(1), Leff, previous_radius, 0.1d0)

      !【模式1】如果开启控制标志7，则用角动量比修正 Ledd
      if (s%x_logical_ctrl(7)) then
         write(*,*) '----------------------------------------'
         write(*,*) 'J change Leddington: Ledd, j/j_c', Ledd, s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit)
         write(*,*) '----------------------------------------'
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if

      ! 计算 Bondi 半径
      R_bondi = 2d0 * standard_cgrav * s%m(1) / pow2(const_csb)
      ! 基本 Bondi 增质量率
      gain = pi * pow2(R_bondi) * const_rho * const_csb

      !【默认模式】根据 Leff 与 Ledd 的关系调整增质量率
      if (s%x_logical_ctrl(1)) then
         if (Leff < Ledd) then
            gain = gain * (1d0 - Leff / Ledd)**2d0  ! 参考 Alexander Dittmann 的建议
         else
            gain = 0.0d0
         end if
      else
         gain = gain * (1d0 - tanh(abs(Leff / Ledd)))
      end if

      ! 若启用了质量截断，进一步修正增质量率（此处用平滑 tanh 过渡）
      if (use_mass_cutoff) then
         gain = gain * 0.5d0 * (1d0 - tanh((s%m(1) / Msun - cutoff_mass) / (0.05d0 * cutoff_mass)))
      end if

      ! 以下块为针对 AGN 环境中盘结构的各种修正
      if (s%x_logical_ctrl(3)) then
         H_disk = sqrt(2d0) * const_csb / Omega_AGN
         erf_arg = -R_bondi / H_disk
         erf_f1 = (-2d0 / sqrt(pi)) * sqrt(1d0 - exp(-erf_arg**2))
         erf_f2 = sqrt(pi)*0.5d0 + (31d0/200d0)*exp(-erf_arg**2) - (341d0/8000d0)*exp(-2d0*erf_arg**2)
         gain = 0.5d0 * gain * sqrt(pi) * erf_f1 * erf_f2 / erf_arg
      end if

      if (s%x_logical_ctrl(4)) then
         H_disk = sqrt(2d0) * const_csb / Omega_AGN
         R_AGN = (standard_cgrav * Mass_AGN * Msun / (Omega_AGN**2))**(1d0/3d0)
         vortpar = (R_bondi**2 * Omega_AGN) / (2d0 * const_csb * R_AGN)
         vortpar = (2d0 / (pi * vortpar)) * asinh((2d0 * vortpar)**(1d0/3d0))
         gain = gain * (1d0 + vortpar**(-10d0))**(-0.1d0)
      end if

      if (s%x_logical_ctrl(5)) then
         R_Hill = (standard_cgrav * s%m(1) / (3d0 * Omega_AGN**2))**(1d0/3d0)
         tidepar = (R_Hill / R_bondi)**2d0
         gain = gain * min(1d0, tidepar)
      end if

      if (s%x_logical_ctrl(6)) then 
         H_disk = sqrt(2d0) * const_csb / Omega_AGN 
         R_Hill = (standard_cgrav * s%m(1) / (3d0 * Omega_AGN**2))**(1d0/3d0) 
         tidepar = min(1d0, R_Hill / R_bondi)
         tidepar = min(tidepar, H_disk / R_bondi)
         tidepar = tidepar * min(1d0, R_Hill / R_bondi)
         tidepar = tidepar * min(1d0, exp(2d0*(1d0 - (R_Hill / H_disk)**2d0)))
         gain = gain * tidepar
      end if

      !【模式2：Chen 方案】如果开启控制标志10，则采用 Chen 方案修正增质量率
      if (s%x_logical_ctrl(10)) then 
         Leff_Ledd_ratio = abs(Leff / Ledd)
         if (Leff_Ledd_ratio >= 1d0) then
              chen_x = 0.0d0
         else
            ! === 2025‑05‑23 Chen‑factor patch start ===
            ! 新计算：使用与大气层一致的 facA 定义，便于统一 Eddington 修正
            facA = (Leff_Ledd_ratio / 0.75_dp)**8
            facA = facA / (1.0_dp + facA)
            Ledd_factor = abs(1.0_dp - facA)**2.0_dp             ! 新公式

            ! ---- 旧公式备份（保留以便比较） ----
            ! facA = (Leff_Ledd_ratio / 0.5_dp)**8
            ! facA = facA / (1.0_dp + facA)
            ! Ledd_factor = abs(1.0_dp - facA * sqrt(Leff_Ledd_ratio)) &
            !               ** (1.0_dp / 0.5_dp) / (1.0_dp + Leff_Ledd_ratio**10)
            ! -----------------------------------------
            ! === 2025‑05‑23 Chen‑factor patch end ===
            
              write(*,*) "Debug: Leff_Ledd_ratio =", Leff_Ledd_ratio
              write(*,*) "Debug: Ledd_factor =", Ledd_factor
              
              mdot_rad = Ledd_factor * Ledd * (sqrt(2.0_dp) / (12.0_dp * const_csb**2))
              mdot_bh = Ledd_factor * Ledd * s%r(1) / (1e-16_dp + standard_cgrav * s%m(1))
              mdot_k = 1d0 / (1d0 / mdot_bh + 1d0 / mdot_rad)
              ratio_kb = mdot_k / (1e-6_dp + 4d0 * pi * R_bondi**2 * const_csb * const_rho)
              chen_x = ratio_kb * (sqrt(4d0 + ratio_kb) - sqrt(ratio_kb)) / (sqrt(4d0 + ratio_kb) + sqrt(ratio_kb))
         end if
         
         write(*,*) 'r(1):', s%r(1)
         write(*,*) 'm(1):', s%m(1)
         write(*,*) "mdot_rad:", mdot_rad, "mdot_bh:", mdot_bh, "mdot_k:", mdot_k, "chen_x:", chen_x
         gain = gain * chen_x
      end if
      
      !【额外修正】如果开启控制标志12，则根据角动量比 (J/Jc)² 进一步修正增质量率
      if (s%x_logical_ctrl(12)) then 
          gain = gain * (1d0 - max(1e-2, (s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))**2))
      end if

      gain = gain * mdot_blend
      loss = loss * (1d0 - Ledd_factor)
      s%mstar_dot = gain + loss

      last_gain = gain * s%dt

      ! 当有质量损失时，同时更新补充材料的化学成分
      if (loss < 0d0) then
          call set_xa(id, ierr, min(gain, abs(loss)) * s%dt, abs(loss) * s%dt)
      end if

      call update_J_dist(id, gain, loss, R_bondi)

      s%job%extras_rpar(i_gain) = gain
      s%job%extras_rpar(i_loss) = loss
      s%xtra1_array(1) = gain
      s%xtra1_array(2) = loss
      s%xtra1_array(9) = Leff / Ledd

   end subroutine other_adjust_mdot

end module mass_changes
