module rotation

   use star_lib
   use star_def
   use const_def
   use math_lib
   use agn_stars
   use utils_lib

   implicit none

   logical :: first_call = .true.
   real(dp) :: time_since_last_update = 0.0_dp
   real(dp), parameter :: update_interval = 1d10! Define the update interval

   public
   contains

   subroutine update_J_dist(id, gain, loss, R_bondi)
      real(dp) :: mu_J, var_J, dt, R_bondi, M, R, gain, loss
      real(dp) :: R_acc, h, R_Hill, j_gain_std, j_gain_avg, R_ph
      real(dp) :: J_crit, sigma
      real(dp) :: I_star, lambda, tau, vt
      real(dp) :: jdot_gain, jdot_loss, net_delta_j
      real(dp) :: random_factor  ! Gaussian random variable
      real(dp) :: am_loss_factor  ! 角动量损失系数
      real(dp) :: original_sign, clipped_mu_J
      
      ! Declare these additional variables at the module or subroutine level
        real(dp) :: random_factor_old = 0.0_dp
        real(dp) :: random_factor_new = 0.0_dp
        real(dp) :: smoothing_factor
        real(dp) :: update_interval_local
        logical :: initialized_random = .false.

      integer, intent(in) :: id
      integer :: ierr, i
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      

    ! Assume s%x_ctrl(13) stores your update interval from inlist
    update_interval_local = s%x_ctrl(13)

      if (ierr /= 0) return

!      write(*,*) 'id', id

      M = s%m(1)
      R = s%r(1)
      dt = s%dt
      time_since_last_update = time_since_last_update + dt

      I_star = 0d0
      do i=1,s%nz
         I_star = I_star + (2d0 / 3d0) * pow2(s%r(i)) * s%dm(i)
      end do
      lambda = M * pow2(R) / I_star

      h = sqrt(2d0) * const_csb / Omega_AGN

      ! Rh = a * (Mstar / 3 MBH)^{1/3}
      ! = a * [(cs^2 R_bondi / 2 G) (1/3) (G / a^3 Omega^2)]^{1/3}
      ! = [(cs^2 R_bondi / Omega^2) (1/6)]^{1/3}
      ! = [(R_bondi h^2) (1/12)]^{1/3}
      R_Hill = pow((1d0 / 12d0) * pow2(h) * R_bondi, 1d0/3d0)
!      write(*,*) 'R_Hill', R_Hill
!      write(*,*) 'R_Bondi', R_Bondi
!      ! write(*,*) 'R_Star', R
      R_acc = min(R_Hill, R_Bondi)
!      write(*,*) 'R_acc', R_acc
!      write(*,*) 'j_R_acc', pow2(R_acc) * Omega_AGN

      vt = const_csb * min(1d0, pow(R_acc / h, 1d0/3d0))
      tau = min(h, R_acc) / vt


    ! Initialize random factors on the very first call
    if (.not. initialized_random) then
       call random_number_uniform(random_factor_old)
       call random_number_uniform(random_factor_new)
       time_since_last_update = 0.0_dp
       initialized_random = .true.
       write(*,*) 'Initial random_factor set:', random_factor_old
    end if

    ! Check if it's time to update the random number
    if (time_since_last_update >= update_interval_local) then
       random_factor_old = random_factor_new
       call random_number_uniform(random_factor_new)
       time_since_last_update = 0.0_dp
       write(*,*) 'Updated random_factor_old:', random_factor_old, &
                  'random_factor_new:', random_factor_new
    end if

    ! Compute smoothing factor (linear interpolation)
    smoothing_factor = time_since_last_update / update_interval_local
    smoothing_factor = min(max(smoothing_factor, 0.0_dp), 1.0_dp)  ! Clamp between 0 and 1

    ! Smooth transition between old and new random factor
    random_factor = (1.0_dp - smoothing_factor) * random_factor_old &
                    + smoothing_factor * random_factor_new

    write(*,*) 'Smoothly transitioned random_factor:', random_factor
      ! Time step quantities
      j_gain_avg = random_factor * pow(standard_cgrav * M * R, 0.5d0)
      j_gain_std = abs(random_factor) * vt * vt * tau  ! Adjusted for random factor
      
!      write(*,*) 'j_gain_avg', j_gain_avg
!      ! write(*,*) 'j_gain_std', j_gain_std

      ! Previous step quantities
      mu_J = s%job%extras_rpar(i_mu_J)
      var_J = s%job%extras_rpar(i_var_J)
      jdot_gain = s%job%extras_rpar(i_jdot_gain)
      jdot_loss = s%job%extras_rpar(i_jdot_loss)

      ! Compute J_crit
      ! J_crit = I Omega_crit
      ! Omega_crit = Sqrt[G M / R^3]
      ! So
      ! J_crit =  I Sqrt[G M / R^3]
      J_crit = (1d0 / lambda) * pow(standard_cgrav * pow3(M) * R, 0.5d0)

      ! Initialize mu_J to J_crit if it is the first call
      if (first_call) then
          mu_J = 0.001d0 * J_crit
          first_call = .false.
!          write(*,*) 'Initialize mu_J ~ J_crit if it is the first call'
      end if


       am_loss_factor = 1.0d0  ! 默认值
       !if (loss > gain) am_loss_factor = 0.1d0
       
      !  Compute next step
       jdot_loss = loss * (abs(mu_J) / M) * lambda * am_loss_factor 
       jdot_loss = jdot_loss * (sign(1.0d0, mu_J))
       jdot_gain = gain * j_gain_avg

       net_delta_j = dt * ( jdot_gain + jdot_loss )

        ! 1. 保存原始符号（1或-1）
        original_sign = sign(1.0d0, mu_J)

        ! 2. 计算理论新值
        mu_J = mu_J + dt * (jdot_gain + jdot_loss)

        ! 3. 智能截断逻辑（数学方式实现）
        clipped_mu_J = max(0.0d0, mu_J * original_sign)  ! 当符号反转时结果为0
        mu_J = clipped_mu_J * original_sign              ! 恢复原始符号方向

       !mu_J = mu_J + net_delta_j
       var_J = var_J + dt * (pow2(gain * j_gain_std) / tau + lambda * loss * (var_J / M))
        
        if (mod(s%model_number, 100) == 0) then
            write(*,*) 'AM loss factor:', am_loss_factor, ' Loss/Gain:', loss/gain
        end if
       write(*,*) 'loss_jdot', jdot_loss
       write(*,*) 'gain_jdot', jdot_gain
       write(*,*) 'net_delta_j', net_delta_j

       ! If the mean J exceeds critical, we limit it to critical.
       if (mu_J > J_crit) then
          ! In practice this is what the truncated Gaussian turns into.
          ! We approach critical so quickly that we overshoot dramatically, leading
          ! to numerical problems with a fancier approach.
          mu_J = J_crit
          var_J = 0d0
       end if

     write(*,*) 'mu_J', mu_J
!  !     write(*,*) 'var_J', var_J
!  !     write(*,*) 'J_crit', J_crit
     write(*,*) 'J/J_crit', mu_J / J_crit
!      !write(*,*) 'j_gain_std', j_gain_std
!      !write(*,*) 'j_gain_avg', j_gain_avg
!  !     write(*,*) 'gain', gain
!  !     write(*,*) 'loss', loss
!!      write(*,*) 'tau', tau
!     !  write(*,*) 'lambda', lambda
!      ! write(*,*) 'sigma', sigma

      if (is_bad(mu_J)) then
          write(*,*) 'Bad mu_J, halting.'
          stop
      end if

      s%job%extras_rpar(i_mu_J) = mu_J
      s%job%extras_rpar(i_var_J) = var_J
      s%job%extras_rpar(i_J_crit) = J_crit

     ! Save jdot_gain and jdot_loss to star_info structure
     !s%job%extras_rpar(i_jdot_gain) = 0.0_dp  ! or keep as is if necessary
     !s%job%extras_rpar(i_jdot_loss) = 0.0_dp  ! or keep as is if necessary
      ! Save jdot_gain and jdot_loss to star_info structure
      s%job%extras_rpar(i_jdot_gain) = jdot_gain
      s%job%extras_rpar(i_jdot_loss) = jdot_loss
      
   end subroutine update_J_dist

   ! ! Function to generate a Gaussian random variable using Box-Muller transform
   ! subroutine random_number_gaussian(random_value)
   !    real(dp), intent(out) :: random_value
   !    real(dp) :: u1, u2
   !    call random_number(u1)
   !    call random_number(u2)
   !    random_value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   ! end subroutine random_number_gaussian

   ! Function to generate a uniform random variable in range [-1, 1]
   subroutine random_number_uniform(random_value)
      real(dp), intent(out) :: random_value
      call random_number(random_value)
      random_value = 2.0_dp * random_value - 1.0_dp  ! Transform to range [-1, 1]
   end subroutine random_number_uniform

end module rotation
