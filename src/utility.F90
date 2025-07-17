!+---------------------------------------------------------------------+
!| This module contains utility subroutines                            |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 10-08-2022  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module utility
  !
  use stlaio,  only: get_unit
  !
  implicit none
  !
  contains
  !+-------------------------------------------------------------------+
  !| This subroutine is to init a text file, either create a new file  |
  !| to resume an old file.                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2023: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine listinit(filename,handle,firstline)
    !
    use commvar,   only: nstep
    use strings,   only: split
    !
    character(len=*),intent(in) :: filename
    integer,intent(out) :: handle
    character(len=*),intent(in),optional :: firstline
    !
    character(len=16),allocatable :: args(:)
    logical :: fex
    logical :: lrestart = .false.
    integer :: nargs,ns,ferr,n
    character(len=120) :: txtformat
    !
    inquire(file=filename,exist=fex)
    handle=get_unit()
    !
    open(handle,file=filename)
    !
    if(lrestart .and. fex) then
      ! resume a file
      ns=0
      read(handle,*)
      ! first line is alway a text 
      !
      do while(ns<nstep)
        !
        read(handle,*,iostat=ferr)ns
        !
        if(ferr< 0) then
          print*,' ** ns,nstep=',ns,nstep
          print*,' ** end of file is reached.'
          exit
        endif
        !
      enddo
      !
      backspace(handle)
      write(*,'(A,I0)')'   ** resume'//filename//'at step: ',ns
      !
    else
      ! create a file
      call split(firstline,args,' ')
      !
      nargs=size(args)
      !
      write(txtformat,'(A,I0,A)')'(',nargs,'(1X,A20))'
      !
      write(handle,txtformat)(trim(args(n)),n=1,nargs)
      !
    endif
    !
  end subroutine listinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listinit.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to write listing data.                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2023: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine listwrite(handle,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,&
                        var11,var12,var13,var14,var15,var16,var17)
    !
    use commvar, only: nstep,time
    !
    integer,intent(in) :: handle
    real(8),intent(in),optional :: var1,var2,var3,var4,var5,var6,      &
                                   var7,var8,var9,var10,var11,var12,   &
                                   var13,var14,var15,var16,var17
    real(8) :: ref_tim = 1.d0
    !
    if(present(var17)) then
      write(handle,"(1X,I20,18(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12,  &
                                          var13,var14,var15,var16,var17
    elseif(present(var16)) then
      write(handle,"(1X,I20,17(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12,  &
                                          var13,var14,var15,var16
    elseif(present(var15)) then
      write(handle,"(1X,I20,16(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12,  &
                                          var13,var14,var15
    elseif(present(var14)) then
      write(handle,"(1X,I20,15(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12,  &
                                          var13,var14
    elseif(present(var13)) then
      write(handle,"(1X,I20,14(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12,  &
                                          var13
    elseif(present(var12)) then
      write(handle,"(1X,I20,13(1X,E20.13E2))")nstep,time/ref_tim,var1,var2,     &
                                          var3,var4,var5,var6,var7,     &
                                          var8,var9,var10,var11,var12
    elseif(present(var11)) then
      write(handle,"(1X,I20,12(1X,E20.13E2))")nstep,time/ref_tim,       &
                 var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11
    elseif(present(var10)) then
      write(handle,"(1X,I20,11(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8,var9,var10
    elseif(present(var9)) then
      write(handle,"(1X,I20,10(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8,var9
    elseif(present(var8)) then
      write(handle,"(1X,I20,9(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8
    elseif(present(var7)) then
      write(handle,"(1X,I20,8(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7
    elseif(present(var6)) then
      write(handle,"(1X,I20,7(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6
    elseif(present(var5)) then
      write(handle,"(1X,I20,6(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5
    elseif(present(var4)) then
      write(handle,"(1X,I20,5(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4
    elseif(present(var3)) then
      write(handle,"(1X,I20,4(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3
    elseif(present(var2)) then
      write(handle,"(1X,I20,3(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2
    elseif(present(var1)) then
      write(handle,"(1X,I20,2(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1
    else
      stop ' !! error @ listwrite'
    endif
    !
  end subroutine listwrite
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listwrite.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function Verifies that a character string represents a       |
  !|  numerical value                                                  |
  !+-------------------------------------------------------------------+
  ! ref: http://fcode.cn/code_gen-115-1.html
  !+-------------------------------------------------------------------+
  Integer Function IsNum(zval)
    ! 确定字符是否是数值类型：
    ! 0-非数值的字符串
    ! 1-整数(integer)
    ! 2-小数(fixed point real)
    ! 3-指数类型实数(exponent type real)
    ! 4-双精度实数指数形式(exponent type double)
    Character (Len=*), Intent (In) :: zval
    !
    Integer :: num, nmts, nexp, kmts, ifexp, ichr
    !
    Integer, Parameter :: kint = 1 ! integer
    Integer, Parameter :: kfix = 2 ! fixed point real
    Integer, Parameter :: kexp = 3 ! exponent type real
    Integer, Parameter :: kdbl = 4 ! exponent type double
    !
    ! initialise
    num = 0  ! 数字的格式，最后传递给ISNUM返回
    nmts = 0 ! 整数或浮点数的数字个数
    nexp = 0 ! 指数形式的数字个数
    kmts = 0 ! 有+-号为1，否则为0
    ifexp = 0! 似乎没用
    ! loop over characters
    ichr = 0
    !
    Do
    
      If (ichr>=len(zval)) Then
    
        ! last check
    
        If (nmts==0) Exit
    
        If (num>=kexp .And. nexp==0) Exit
    
        isnum = num
    
        Return
    
      End If
    
      ichr = ichr + 1
    
      Select Case (zval(ichr:ichr))
    
        ! process blanks
    
      Case (' ')
    
        Continue
    
        ! process digits
    
      Case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
    
        If (num==0) num = kint
    
        If (num<kexp) Then
    
          nmts = nmts + 1
    
          ! 整数或浮点数+1
    
        Else
    
          nexp = nexp + 1
    
          ! 指数形式+1
    
        End If
    
        ! process signs
    
      Case ('+', '-')
    
        If (num==0) Then
    
          If (kmts>0) Exit
    
          ! 出现2个符号，非数字
    
          kmts = 1
    
          num = kint
    
        Else
    
          If (num<kexp) Exit
    
          If (ifexp>0) Exit
    
          ifexp = 1
    
        End If
    
        ! process decimal point
    
      Case ('.')
    
        If (num/=kint .And. ichr/=1) Exit
    
        ! 前面不是整数，小数点也不是第一个字符，则非数字
    
        num = kfix
    
        ! process exponent
    
      Case ('e', 'E')
    
        If (num>=kexp) Exit
    
        If (nmts==0) Exit
    
        num = kexp
      Case ('d', 'D')
    
        If (num>=kexp) Exit
    
        If (nmts==0) Exit
    
        num = kdbl
    
        ! any other character means the string is non-numeric
    
      Case Default
    
        Exit
    
      End Select
    
    End Do
    
    ! if this point is reached, the string is non-numeric
    
    isnum = 0
    
    Return
    
  End Function IsNum
  !+-------------------------------------------------------------------+
  !| The end of the Function IsNum.                                    |
  !+-------------------------------------------------------------------+
  !
  function rnorm_box_muller(mode) result(variates) 
    !
    use constdef
    ! coded formulas from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    !
    character(len=*),intent(in) :: mode
    !
    ! return two uncorrelated standard normal variates
    integer,allocatable :: seed(:)
    integer :: rantime(8)
    !
    integer :: n
    real(8) :: variates(2)
    real(8) :: u(2), factor, arg
    !
    logical,save :: firstcall=.true.
    !
    if(mode=='sync') then
      ! all processor generate same random numbers
      if(firstcall) then
        !
        call random_seed(size = n)
        allocate(seed(n))
        call date_and_time(values=rantime) 
        !  use date and minutes for synthetisation 
        !    1      2    3    4     5      6      7    8
        !-----------------------------------------------
        ! 2023     10    5   60    23      6     28  962
        !-----------------------------------------------
        ! year  month date   ??  hour minute second msec
        !
        seed=0
        seed(1:6)=rantime(1:6)
        !
        call random_seed(put=seed)
        !
        deallocate(seed)
        !
        firstcall=.false.
      endif
      !
    endif
    !
    do
       call random_number(u)
       if (u(1) > 0.d0) exit
    end do
    factor = sqrt(-2 * log(u(1)))
    arg = 2.d0*pi*u(2)
    variates = factor * [cos(arg),sin(arg)]
    !
  end function rnorm_box_muller
  !
end module utility
!+---------------------------------------------------------------------+
!| The end of the module utility                                       |
!+---------------------------------------------------------------------+