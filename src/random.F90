module random
  implicit none
  private
  public :: h

contains

  ! 主函数：h(k1, k2, r1, r2) ∈ (-1, 1)，满足奇对称性
  real(8) function h(k1, k2, r1, r2) result(val)
    real(8), intent(in) :: k1, k2, r1, r2
    real(8) :: g_pos, g_neg

    g_pos = hash_to_uniform(k1, k2, r1, r2)
    g_neg = hash_to_uniform(-k1, -k2, r1, r2)
    val = g_pos - g_neg   ! ∈ (-1, 1)
  end function h

  ! 辅助函数：将 (a,b,c,d) 映射为 [0,1) 的确定性伪随机值
  real(8) function hash_to_uniform(k1, k2, r1, r2) result(u)
    real(8), intent(in) :: k1, k2, r1, r2
    integer(8) :: i1, i2, i3, i4, h

    i1 = transfer(k1, i1)
    i2 = transfer(k2, i2)
    i3 = transfer(r1*614615.d0, i3)
    i4 = transfer(r2*406064.d0, i4)

    ! 混合所有输入
    h = i1
    h = ieor(h, shiftl(i2, 13))
    h = ieor(h, shiftl(i3, 23))
    h = ieor(h, shiftl(i4, 7))

    h = h * 7919_8
    h = ieor(h, ishft(h, -11))
    h = h * 9187_8
    h = ieor(h, ishft(h, 17))

    ! 转为 real 并缩放到 (-0.5, 0.5)
    u = dble(h) * (1.0d0 / 9223372036854775808.0d0)  ! 1 / 2^63
  end function hash_to_uniform

end module random