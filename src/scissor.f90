!> Scissor correction: every conduction band of the window is shifted rigidly, E_c -> E_c + Delta (keyword Scissor,
!! eV). It is applied at the response stage to the single-particle data read from the .omesp, which always holds
!! the unshifted model, so one .omesp serves any Delta. Only conduction-valence (cv) pairs change: same-side
!! energy differences, intraband velocities and same-side position elements are untouched by a rigid shift.
!!
!! Keyword Scissor_mode chooses what is kept fixed (w_cv = E_c - E_v before the shift, g = w_cv/(w_cv + Delta)):
!!   r_unchanged (default): r_cv kept, v_cv -> v_cv/g. The scissor as a nonlocal operator Delta P_c commuting with
!!     the position operator in the length gauge: Levine & Allan, PRL 63, 1719 (1989); second order, Nastos et al.,
!!     PRB 72, 045223 (2005).
!!   v_unchanged: v_cv kept, r_cv -> g r_cv, so the generalised derivative picks up r_cv d_a g with
!!     d_a g = Delta (v^a_cc - v^a_vv)/(w_cv + Delta)^2. This is what Xatu's own Kubo conductivity does (the
!!     BSE diagonal and the transition energies are shifted, the velocity matrix elements are not).
!! Band index convention of the window arrays: bands 1..nv_ex valence, nv_ex+1..nb conduction.
module scissor
  use parser_input_file,      only: scissor_ev, scissor_r_unchanged
  use parser_optics_xatu_dim, only: nv_ex
  implicit none
  private
  public :: scissor_on, scissor_ha, scissor_report
  public :: scissor_mesh_ev, scissor_mesh_ex, scissor_mesh_shift_cov, scissor_mesh_stencil

  real(8), parameter :: ha_ev = 27.211385d0
  real(8), parameter :: wmin = 1.0d-8          ! Hartree: a cv pair closer than this is a metal at this filling

contains

  logical function scissor_on()
    scissor_on = scissor_ev /= 0.0d0
  end function scissor_on

  real(8) function scissor_ha()
    scissor_ha = scissor_ev/ha_ev
  end function scissor_ha

  !> One log line naming the shift and the regime, printed where the shift is applied.
  subroutine scissor_report(where)
    character(len=*), intent(in) :: where
    if (.not. scissor_on()) return
    if (scissor_r_unchanged) then
      write(*,'(a,f10.5,a)') '    Scissor: conduction bands shifted by ', scissor_ev, &
           ' eV, r_unchanged (v_cv scaled by (w_cv + Delta)/w_cv) -- '//trim(where)
    else
      write(*,'(a,f10.5,a)') '    Scissor: conduction bands shifted by ', scissor_ev, &
           ' eV, v_unchanged (r_cv scaled by w_cv/(w_cv + Delta)) -- '//trim(where)
    end if
  end subroutine scissor_report

  !> Stops if a conduction band of the window is not above a valence band at some k (a scissor needs a gap: the
  !! rescalings divide by w_cv). A negative Scissor that closes the gap is not checked.
  subroutine require_gap(w, where)
    real(8), intent(in) :: w
    character(len=*), intent(in) :: where
    if (w > wmin) return
    write(*,'(a,es10.3,a)') ' ERROR (scissor): a conduction band lies at or below a valence band (E_c - E_v = ', &
         w*ha_ev, ' eV) in '//trim(where)//'.'
    write(*,*) '       A scissor shifts the conduction bands of an insulator; check Nfermi and Bandlist.'
    stop 1
  end subroutine require_gap

  !> Energies ek(nk, nb) (Hartree) and velocities v(nk, 3, nb, nb) of the window, as the linear Kubo kernel and the
  !! shift kernel use them: conduction energies + Delta; r_unchanged also scales v_cv (both orders) by 1/g.
  subroutine scissor_mesh_ev(nk, nb, ek, v)
    integer,    intent(in)    :: nk, nb
    real(8),    intent(inout) :: ek(nk, nb)
    complex(8), intent(inout) :: v(nk, 3, nb, nb)
    real(8) :: d, w
    integer :: k, n, m
    if (.not. scissor_on()) return
    d = scissor_ha()
    do k = 1, nk
      do m = 1, nv_ex
        do n = nv_ex + 1, nb
          w = ek(k,n) - ek(k,m)
          call require_gap(w, 'the single-particle energies')
          if (scissor_r_unchanged) then
            v(k,:,n,m) = v(k,:,n,m)*(w + d)/w
            v(k,:,m,n) = v(k,:,m,n)*(w + d)/w
          end if
        end do
      end do
      ek(k, nv_ex+1:nb) = ek(k, nv_ex+1:nb) + d
    end do
  end subroutine scissor_mesh_ev

  !> Excitonic route: the window velocities vme(nk,3,nb,nb) and positions xme(nk,3,nb,nb) from which P_n and X_n
  !! are built. ek is NOT shifted here: the exciton energies come from Xatu (which applies its own scissor) and the
  !! single-particle energies are only used for degeneracy tests on the same side of the gap.
  !! r_unchanged: v_cv -> v_cv/g (P_n changes, X_n does not); v_unchanged: r_cv -> g r_cv (X_n changes).
  subroutine scissor_mesh_ex(nk, nb, ek, vme, xme)
    integer,    intent(in)    :: nk, nb
    real(8),    intent(in)    :: ek(nk, nb)
    complex(8), intent(inout) :: vme(nk, 3, nb, nb), xme(nk, 3, nb, nb)
    real(8) :: d, w, s
    integer :: k, n, m
    if (.not. scissor_on()) return
    d = scissor_ha()
    do k = 1, nk
      do m = 1, nv_ex
        do n = nv_ex + 1, nb
          w = ek(k,n) - ek(k,m)
          call require_gap(w, 'the single-particle energies')
          if (scissor_r_unchanged) then
            s = (w + d)/w
            vme(k,:,n,m) = vme(k,:,n,m)*s; vme(k,:,m,n) = vme(k,:,m,n)*s
          else
            s = w/(w + d)
            xme(k,:,n,m) = xme(k,:,n,m)*s; xme(k,:,m,n) = xme(k,:,m,n)*s
          end if
        end do
      end do
    end do
  end subroutine scissor_mesh_ex

  !> Response = shift (covariant): energies ek, raw velocities vraw and the block-covariant generalised derivative
  !! rgen(nk, a, c, n, m) = r^{c;a}_nm. The kernel builds r = -i v/(E_n - E_m) itself. r_unchanged: v_cv scaled so
  !! that r_cv and rgen stay. v_unchanged: v kept, so the kernel's r becomes g r, and rgen -> g rgen + r d_a g
  !! (d_a g = Delta (v^a_cc - v^a_vv)/(w + Delta)^2; zero along a non-periodic direction, where v^a_nn = 0).
  !! g is real, so the r d_a g term drops out of the b <-> c symmetrised shift current; it is added so that rgen
  !! stays the derivative of the r the kernel builds.
  subroutine scissor_mesh_shift_cov(nk, nb, ek, vraw, rgen)
    integer,    intent(in)    :: nk, nb
    real(8),    intent(inout) :: ek(nk, nb)
    complex(8), intent(inout) :: vraw(nk, 3, nb, nb), rgen(nk, 3, 3, nb, nb)
    complex(8), parameter :: ci = (0.0d0, 1.0d0)
    real(8) :: d, w, g, dg(3)
    complex(8) :: rnm(3), rmn(3)
    integer :: k, n, m, a
    if (.not. scissor_on()) return
    d = scissor_ha()
    if (.not. scissor_r_unchanged) then
      do k = 1, nk
        do m = 1, nv_ex
          do n = nv_ex + 1, nb
            w = ek(k,n) - ek(k,m)
            call require_gap(w, 'the single-particle energies')
            g = w/(w + d)
            do a = 1, 3
              dg(a) = d*dble(vraw(k,a,n,n) - vraw(k,a,m,m))/(w + d)**2
            end do
            rnm = -ci*vraw(k,:,n,m)/w
            rmn = ci*vraw(k,:,m,n)/w
            do a = 1, 3
              rgen(k,a,:,n,m) = g*rgen(k,a,:,n,m) + rnm*dg(a)
              rgen(k,a,:,m,n) = g*rgen(k,a,:,m,n) + rmn*dg(a)
            end do
          end do
        end do
      end do
    end if
    call scissor_mesh_ev(nk, nb, ek, vraw)
  end subroutine scissor_mesh_shift_cov

  !> Covariant method B (shg, electrooptic, rectification, general): energies nbe(nk, p, n) and positions
  !! nbr(nk, p, c, n, m) at the centre (p = 7) and the six mesh neighbours. Everything is built from r and E at
  !! each point, so the shift is applied point by point and the k-derivative of the transition fields follows:
  !! conduction energies + Delta at every point; v_unchanged also scales r_cv by g at every point. Neighbours of a
  !! non-periodic direction are stored as zeros and are left as they are.
  subroutine scissor_mesh_stencil(nk, nb, nbe, nbr)
    integer,    intent(in)    :: nk, nb
    real(8),    intent(inout) :: nbe(nk, 7, nb)
    complex(8), intent(inout) :: nbr(nk, 7, 3, nb, nb)
    real(8) :: d, w, g
    integer :: k, p, n, m
    if (.not. scissor_on()) return
    d = scissor_ha()
    do k = 1, nk
      do p = 1, 7
        if (all(nbe(k,p,:) == 0.0d0)) cycle
        if (.not. scissor_r_unchanged) then
          do m = 1, nv_ex
            do n = nv_ex + 1, nb
              w = nbe(k,p,n) - nbe(k,p,m)
              call require_gap(w, 'the single-particle energies')
              g = w/(w + d)
              nbr(k,p,:,n,m) = g*nbr(k,p,:,n,m); nbr(k,p,:,m,n) = g*nbr(k,p,:,m,n)
            end do
          end do
        end if
        nbe(k,p,nv_ex+1:nb) = nbe(k,p,nv_ex+1:nb) + d
      end do
    end do
  end subroutine scissor_mesh_stencil

end module scissor
