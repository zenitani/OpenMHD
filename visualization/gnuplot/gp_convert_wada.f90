!-----------------------------------------------------------------------
!     data converter for gnuplot
! -*- compile-command: "gfortran -O2 gp_convert_wave.f90" -*-
!-----------------------------------------------------------------------
!     2016/03/03  T. Wada      original version
!     2016/03/21  S. Zenitani  modified for the latest data format
!-----------------------------------------------------------------------
program gp_convert_wada
  implicit none
  include 'param.h'
  integer :: argc  
  character(256),dimension(100) :: argvarray
  character*32 :: TYPE
  character*256 :: INPUT
  character*256 :: OUTPUT
  integer :: i, j
  integer :: is=0, ie=0, js=0, je=0
  real(8) :: dx,dy,dz,dvx,dvy,dvz
  real(8) :: dvarx,dvary,dvarz !dro,den用変数
  !読み取るデータ用の変数の定義 output.f90 と合わせる
  real(8) :: t
  integer :: ix, jx !x、y方向のグリッドの数
  real(8),allocatable,dimension(:) :: x,y
  real(8),allocatable,dimension(:,:,:) :: U
  real(8),allocatable,dimension(:,:,:) :: V
  real(8),allocatable,dimension(:,:) :: jx_,jy_,jz_ !電流密度
  real(8),allocatable,dimension(:,:) :: ARY0 !Az,div_vなど汎用
  real(8) :: Az0,vmag,Ez,vpara,vpara2,vperp,tmpvar
  real(8) :: coef
  real(8) :: x0,x1,y0,y1
  integer :: xyswitch
  !gnuplotの都合でixjとjxiのデータを作らないといけない
  !0なら内側のループがi、1なら内側ループはj
  !実行時の引数から入力ファイルと出力ファイルを決める
  argc = iargc()
  !使い方
  if (argc<3) then
     write(*,*) 'USAGE'
     write(*,*) './a.out TYPE INPUT OUTPUT x0 x1 y0 y1 > output.dat'
     write(*,*) './a.out TYPEBIN INPUT OUTPUT x0 x1 y0 y1'
     write(*,*) './a.out TYPEIDL INPUT OUTPUT x0 x1 y0 y1'
!     write(*,*) 'TYPE: crop_profile'
!     write(*,*) './a.out crop_profile INPUT OUTPUT x0 x1 y0 y1 xyswitch'
!     write(*,*) 'xyswitch=0 < for x profile'
!     write(*,*) 'xyswitch=1 < for y profile'
!./data-read-test00 rho data/field-00000.dat output | gzip > tmp.dat.gz
     stop
  end if
  !引数は最大100までとしておく
  if (argc .ge. 100) then
     write(*,*) 'Number of argument is too large. argc=', argc
     stop
  endif
  !引数の代入
  do i=1, argc
     call getarg(i, argvarray(i))
  end do
!  write(*,*) argc
!  stop
  TYPE=argvarray(1)

  write(*,*)"TYPE=",TYPE
!  if (TYPE .eq. "all" ) then
!     write(*,*)'aaa'
!  end if
!  stop


  INPUT=argvarray(2)
  OUTPUT=argvarray(3)

  read (argvarray(4),*) x0
  read (argvarray(5),*) x1
  read (argvarray(6),*) y0
  read (argvarray(7),*) y1
!  if( (argc .eq. 8) .and. (TYPE .eq. "crop_profile" ) )then
!     read (argvarray(8),*) xyswitch
!     write(*,*) '#crop region: x0,x1,y0,y1=',x0,x1,y0,y1
!     write(*,*) '#xyswitch=',xyswitch
!  end if

  !データ読み取り
  open (16, file=INPUT,form='unformatted',access='stream')
  read(16) t
  read(16) ix
  read(16) jx


!  write(*,*) 'ix=',ix,'jx=',jx
!  stop

  !ixとjxから動的に配列を確保
  allocate(x(1:ix))
  allocate(y(1:jx))
  allocate(U(1:ix,1:jx,1:var1))
  allocate(V(1:ix,1:jx,1:var2))
  allocate(jx_(1:ix,1:jx))
  allocate(jy_(1:ix,1:jx))
  allocate(jz_(1:ix,1:jx))
  allocate(ARY0(1:ix,1:jx))
  !配列の読み取り開始、fortranって簡単でよいね
  read(16) x
  read(16) y

!  write(*,*) 'x,y=',x(ix),y(jx)
!  stop

  read(16) U(1:ix,1:jx,mx)
  read(16) U(1:ix,1:jx,my)
  read(16) U(1:ix,1:jx,mz)
  read(16) U(1:ix,1:jx,en)
  read(16) U(1:ix,1:jx,ro)
  read(16) U(1:ix,1:jx,bx)
  read(16) U(1:ix,1:jx,by)
  read(16) U(1:ix,1:jx,bz)
  read(16) U(1:ix,1:jx,ps)
  read(16) V(1:ix,1:jx,vx)
  read(16) V(1:ix,1:jx,vy)
  read(16) V(1:ix,1:jx,vz)
  read(16) V(1:ix,1:jx,pr)
  close(16)

  !中心差分で電流密度の計算
  !$omp parallel do &
  !$omp private(i,dx,dy,dz)
  do j=1+1,jx-1
     do i=1+1,ix-1
        jx_(i,j)=0.d0
        jy_(i,j)=0.d0
        jz_(i,j)=0.d0
        dx=x(i+1)-x(i-1)
        dy=y(j+1)-y(j-1)
        dz=0.d0
        jx_(i,j) = (U(i,j+1,bz)-U(i,j-1,bz))/dy - 0.d0
        jy_(i,j) = 0.d0                         - (U(i+1,j,bz)-U(i-1,j,bz))/dx
        jz_(i,j) = (U(i+1,j,by)-U(i-1,j,by))/dx - (U(i,j+1,bx)-U(i,j-1,bx))/dy
     end do
  end do
  !$omp end  parallel do
  !$omp parallel do
  do i=1,ix !x方向の境界の電流
     jx_(i,1)=jx_(i,1+1)
     jy_(i,1)=jy_(i,1+1)
     jz_(i,1)=jz_(i,1+1)
     jx_(i,jx)=jx_(i,jx-1)
     jy_(i,jx)=jy_(i,jx-1)
     jz_(i,jx)=jz_(i,jx-1)
  end do
  !$omp end  parallel do
  !$omp parallel do
  do j=1,jx !y方向の境界の電流xc
     jx_(1,j)=jx_(1+1,j)
     jy_(1,j)=jy_(1+1,j)
     jz_(1,j)=jz_(1+1,j)
     jx_(ix,j)=jx_(ix-1,j)
     jy_(ix,j)=jy_(ix-1,j)
     jz_(ix,j)=jz_(ix-1,j)
  end do
  !$omp end  parallel do

  dx=x(2)-x(1) !全領域でdx=constだから
  dy=y(2)-y(1) !全領域でdy=constだから
  !Azの計算
  if( (TYPE .eq. "Azbin") .or. (TYPE .eq. "Azidl") .or. (TYPE .eq. "crop_profile") ) then
     !Azd(y)=B0*D*log(cosh(y/D))
     !Azu(y)=B0*D/k*log(cosh(y/D))
     !Az0=log(cosh(y(1)))
     coef=1.d0
     if (abs(y(1))<700) then
        Az0=coef*log(cosh(y(1)))
     else
        !y>700だとcoshが大きな値になってエラーが出るので近似式に変える
        Az0=coef*(abs(y(1))-log(2.0))
     end if
     !$omp parallel do &
     !$omp private(i)
     do j=1,jx
        do i=1,ix
           !        ARY0(i,j)=0.d0
           ARY0(i,j)=Az0
        end do
     end do
     !$omp end  parallel do
     !ここは並列化するならsumだかreductionだかにしないといけない
     do j=1,jx
        do i=1,ix
           if(i+1 < ix) then
              !x方向へByを積分
              ARY0(i+1,j)=ARY0(i,j)-coef*U(i,j,by)*dx
           else
              ARY0(1,j)=ARY0(ix,j)-coef*U(i,j,by)*dx
           end if
        end do
     end do
     do i=1,ix
        do j=1,jx
           if(j+1 < jx) then
              !y方向へBxを積分
              ARY0(i,j+1)=ARY0(i,j)+coef*U(i,j,bx)*dy
           else
              ARY0(i,1)=ARY0(i,jx)+coef*U(i,j,bx)*dy
           end if
        end do
     end do
  end if

  !div vの計算
  if( (TYPE .eq. "divvbin") .or. (TYPE .eq. "divvidl") ) then     
     !$omp parallel do &
     !$omp private(i,dx,dy,dz,dvx,dvy,dvz)
     do j=1+1,jx-1
        do i=1+1,ix-1
           dx=x(i+1)-x(i-1)
           dy=y(j+1)-y(j-1)
           dz=0
           dvx=( V(i+1,j,vx)-V(i-1,j,vx) )
           dvy=( V(i,j+1,vy)-V(i,j-1,vy) )
           dvz=0 !使わない
           ARY0(i,j)=dvx/dx+dvy/dy
        end do
     end do
     !$omp end  parallel do
     !$omp parallel do
     do i=1,ix !x方向の境界のdiv_v
        ARY0(i,1)=ARY0(i,1+1)
        ARY0(i,jx)=ARY0(i,jx-1)
     end do
     !$omp end  parallel do
     !$omp parallel do
     do j=1,jx !y方向の境界のdiv_v
        ARY0(1,j)=ARY0(1+1,j)
        ARY0(ix,j)=ARY0(ix-1,j)
     end do
     !$omp end  parallel do
  end if

  !vnroの計算
  !vnenの計算
  !en=cv*log(pr/ro^gamma)なので厳密にはenに関連した量
  !roよりも勾配が激しくて2Dマップにしたときコントラストがつく
  if( (TYPE .eq. "vnrobin") .or. (TYPE .eq. "vnroidl") .or. (TYPE .eq. "vnenbin") .or. (TYPE .eq. "vnenidl") ) then
     !$omp parallel do &
     !$omp private(i,dx,dy,dz,dvarx,dvary,dvarz,tmpvar)
     do j=1+1,jx-1
        do i=1+1,ix-1
           dx=x(i+1)-x(i-1)
           dy=y(j+1)-y(j-1)
           dz=0
           if( (TYPE .eq. "vnrobin") .or. (TYPE .eq. "vnroidl") ) then
              dvarx=( U(i+1,j,ro)-U(i-1,j,ro) )
              dvary=( U(i,j+1,ro)-U(i,j-1,ro) )
              dvarz=0 !使わない
              tmpvar=sqrt(dvarx/dx*dvarx/dx+dvary/dy*dvary/dy)
              ARY0(i,j)=(dvarx/dx*V(i,j,vx)+dvary/dy*V(i,j,vy))/tmpvar
           end if
           if( (TYPE .eq. "vnenbin") .or. (TYPE .eq. "vnenidl") ) then
              dvarx=( V(i+1,j,pr)/(U(i+1,j,ro)**gamma)-V(i-1,j,pr)/(U(i-1,j,ro)**gamma) )
              dvary=( V(i,j+1,pr)/(U(i,j+1,ro)**gamma)-V(i,j-1,pr)/(U(i,j-1,ro)**gamma) )
              dvarz=0 !使わない
              tmpvar=sqrt(dvarx/dx*dvarx/dx+dvary/dy*dvary/dy)
              ARY0(i,j)=(dvarx/dx*V(i,j,vx)+dvary/dy*V(i,j,vy))/tmpvar
           end if
        end do
     end do
     !$omp end  parallel do
     !$omp parallel do
     do i=1,ix !x方向の境界のvnro
        ARY0(i,1)=ARY0(i,1+1)
        ARY0(i,jx)=ARY0(i,jx-1)
     end do
     !$omp end  parallel do
     !$omp parallel do
     do j=1,jx !y方向の境界のvnro
        ARY0(1,j)=ARY0(1+1,j)
        ARY0(ix,j)=ARY0(ix-1,j)
     end do
     !$omp end  parallel do
  end if

  !極端に大きな値がないか確認
!!$  do i=1,ix,100
!!$     do j=1,jx,100
!!$        if(ARY0(i,j)>10000) then
!!$           write(*,*)"Az>1000",i,j,ARY0(i,j)
!!$        end if
!!$     end do
!!$  end do
!!$  do i=1,ix,1000
!!$    do j=1,jx,100
!!$        if((j>12600).and.(j<12900)) then
!!$           write(*,*)"Az>1000",i,j,ARY0(i,j)
!!$        end if
!!$     end do
!!$  end do
!!$  stop

  !gnuplot用可視化データの出力
!  write(*,*)"#x,y,ro,pr,vx,vy,vz,bx,by,bz,ps,jx,jy,jz",TYPE
  write(*,*)"#TYPE=",TYPE
  write(*,*)"#INPUT=",INPUT
  write(*,*)"#OUTPUT=",OUTPUT
  write(*,*)"#t=",t
  write(*,*)"#ix=",ix
  write(*,*)"#jx=",jx

  !後でグラフを描くときにセル数が必要なのでカウント
  do i=1,ix
     if((x0.le.x(i)).and.(x(i).le.x1)) then
        if(is.eq.0) then
           is=i
        endif
        ie=i
     end if
  end do
  do j=1,jx
     if((y0.le.y(j)).and.(y(j).le.y1)) then
        if(js.eq.0) then
           js=j
        endif
        je=j
     end if
  end do
  write(*,*)'#ixcnt=',ie-is+1,' jxcnt=',je-js+1

  !本当は単精度でいい
  !unformated binaryで
  !open(100,file="tmpro.dat",form='unformatted')
  !write(100) U(:,:,ro)
  !close(100)

  !gnuplot用のアスキーデータ出力
  !大規模データの場合は全部バイナリにする。
  !空行の出力が出てしまうので使わない場合はコメントアウトする。
!$  do i=is,ie
!$     do j=js,je
!$        !write(*,*) x(i),y(j),U(i,j,ro),V(i,j,pr),&
!$        !V(i,j,vx),V(i,j,vy),V(i,j,vz),&
!$        !U(i,j,bx),U(i,j,by),U(i,j,bz),U(i,j,ps),&
!$        !jz_(i,j),jy_(i,j),jz_(i,j),ARY0(i,j)
!$        !二次元シミュレーションなので必要最低限な成分だけにする
!$        if (TYPE .eq. "all") then !全部書き出し
!$           write(*,*) x(i),y(j),U(i,j,ro),V(i,j,pr),&
!$                V(i,j,vx),V(i,j,vy),&
!$                U(i,j,bx),U(i,j,by),U(i,j,ps),&
!$                jz_(i,j),ARY0(i,j)
!$        end if
!$        if (TYPE .eq. "ro") then !密度
!$           write(*,*) x(i),y(j),U(i,j,ro)
!$           !ファイルサイズあまりかわらない気がするasciiでは無意味？
!$           !write(*,*) real(x(i)),real(y(j)),real(U(i,j,ro))
!$        end if
!$        if (TYPE .eq. "pr") then !圧力
!$           write(*,*) x(i),y(j),V(i,j,pr)
!$        end if
!$        if (TYPE .eq. "pt") then !全圧
!$           tmpvar=V(i,j,pr)+(U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz))/2.d0
!$           write(*,*) x(i),y(j),tmpvar
!$        end if
!$        if (TYPE .eq. "pm") then !磁気圧
!$           tmpvar=(U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz))/2.d0
!$           write(*,*) x(i),y(j),tmpvar
!$        end if
!$        if (TYPE .eq. "enem") then !磁気圧と磁気エネルギーは同じ？
!$           tmpvar=(U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz))/2.d0
!$           write(*,*) x(i),y(j),tmpvar
!$        end if
!$        if (TYPE .eq. "enekin") then
!$           tmpvar=U(i,j,ro)/2.0*(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz))
!$           write(*,*) x(i),y(j),tmpvar
!$        end if
!$        if (TYPE .eq. "eneint") then 
!$           tmpvar=1.0/(gamma-1.0)*V(i,j,pr)/U(i,j,ro)
!$           write(*,*) x(i),y(j),tmpvar
!$        end if
!$        if (TYPE .eq. "vx") then !vx
!$           write(*,*) x(i),y(j),V(i,j,vx)
!$        end if
!$        if (TYPE .eq. "vy") then !vy
!$           write(*,*) x(i),y(j),V(i,j,vy)
!$        end if
!$        if (TYPE .eq. "vmag") then !vmag
!$           vmag=sqrt(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz))
!$           write(*,*) x(i),y(j),vmag
!$        end if
!$        if (TYPE .eq. "vecv") then !vvec
!$           vmag=sqrt(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy))
!$           if (vmag .le. 2e-3) then !しきい値の設定は微調整が必要かも
!$              write(*,*) x(i),y(j),0,0
!$           else
!$              write(*,*) x(i),y(j),V(i,j,vx)/vmag,V(i,j,vy)/vmag
!$           end if
!$        end if
!$        if (TYPE .eq. "bx") then !bx
!$           write(*,*) x(i),y(j),U(i,j,bx)
!$        end if
!$        if (TYPE .eq. "by") then !by
!$           write(*,*) x(i),y(j),U(i,j,by)
!$        end if
!$        if (TYPE .eq. "jz") then !jz
!$           write(*,*) x(i),y(j),jz_(i,j)
!$        end if
!$        if (TYPE .eq. "Az") then !Az
!$           write(*,*) x(i),y(j),ARY0(i,j)
!$        end if
!$        if (TYPE .eq. "T") then !温度
!$           write(*,*) x(i),y(j),V(i,j,pr)/(U(i,j,ro)+1e-10)
!$        end if           
!$     end do     
!$     write(*,*) !pm3d、contour用にx方向に1増えるときには空行が必要
!$     !ix*iy行を1ブロックとして同じ行のブロックが空行おきで続くこと
!$  end do

  !gnuplot用のバイナリデータ出力
  !http://stackoverflow.com/questions/15878211/pm3d-in-gnuplot-with-binary-data
  !splot 'tmp.bin' binary record=(642,-1) format='%float'
  !642=ix
  !asciiの出力の時、空のバイナリファイルができてしまうけど気にしないこと
  !にする
  open(unit=100,file=OUTPUT,action='write',status='replace',access='stream',form='unformatted')

  do j=js,je
     do i=is,ie
        if (TYPE .eq. "robin") then
           write(100) real(x(i)),real(y(j)),real(U(i,j,ro))
        elseif (TYPE .eq. "prbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,pr))
        elseif (TYPE .eq. "ptbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,pr)+(U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz))/2.d0)
        elseif ((TYPE .eq. "pmbin") .or. (TYPE .eq. "enembin")) then
           write(100) real(x(i)),real(y(j)),real((U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz))/2.d0)
        elseif (TYPE .eq. "enekinbin") then
           write(100) real(x(i)),real(y(j)),real(U(i,j,ro)/2.0*(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz)))
        elseif (TYPE .eq. "eneintbin") then 
           write(100) real(x(i)),real(y(j)),real(1.0/(gamma-1.0)*V(i,j,pr)/U(i,j,ro))
        elseif (TYPE .eq. "vxbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,vx))
        elseif (TYPE .eq. "vybin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,vy))
        elseif (TYPE .eq. "vzbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,vz))
        elseif (TYPE .eq. "vmagbin") then
           vmag=sqrt(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz))
           write(100) real(x(i)),real(y(j)),real(vmag)
        elseif (TYPE .eq. "divvbin") then
           write(100) real(x(i)),real(y(j)),real(ARY0(i,j))
        elseif (TYPE .eq. "vnrobin") then
           write(100) real(x(i)),real(y(j)),real(ARY0(i,j))
        elseif (TYPE .eq. "vnenbin") then
           write(100) real(x(i)),real(y(j)),real(ARY0(i,j))
        elseif (TYPE .eq. "vecvbin") then
           vmag=sqrt(V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz))
           if (vmag .le. 2e-3) then !しきい値の設定は微調整が必要かも
              write(100) real(x(i)),real(y(j)),real(0.0),real(0.0)
           else
              write(100) real(x(i)),real(y(j)),real(V(i,j,vx)/vmag),real(V(i,j,vy)/vmag)
           end if
        elseif (TYPE .eq. "bxbin") then
           write(100) real(x(i)),real(y(j)),real(U(i,j,bx))
        elseif (TYPE .eq. "bybin") then
           write(100) real(x(i)),real(y(j)),real(U(i,j,by))
        elseif (TYPE .eq. "bzbin") then
           write(100) real(x(i)),real(y(j)),real(U(i,j,bz))
        elseif (TYPE .eq. "jzbin") then
           write(100) real(x(i)),real(y(j)),real(jz_(i,j))
        elseif (TYPE .eq. "Azbin") then              
           write(100) real(x(i)),real(y(j)),real(ARY0(i,j))
        elseif (TYPE .eq. "Tbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,pr)/(U(i,j,ro)+1e-10))
        elseif (TYPE .eq. "enbin") then
           write(100) real(x(i)),real(y(j)),real(V(i,j,pr)/(U(i,j,ro)**(gamma)+1e-10))
        elseif (TYPE .eq. "Ezbin") then
           Ez=-(V(i,j,vx)*U(i,j,by)-U(i,j,bx)*V(i,j,vy))
           write(100) real(x(i)),real(y(j)),real(Ez)
        elseif (TYPE .eq. "vparabin") then !磁場方向に対して
           tmpvar=U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz)
           vpara=(V(i,j,vx)*U(i,j,bx)+V(i,j,vy)*U(i,j,by)+V(i,j,vz)*U(i,j,bz))/sqrt(tmpvar)
           write(100) real(x(i)),real(y(j)),real(vpara)
        elseif (TYPE .eq. "vperpbin") then !磁場方向に対して
           tmpvar=U(i,j,bx)*U(i,j,bx)+U(i,j,by)*U(i,j,by)+U(i,j,bz)*U(i,j,bz)
           vpara2=((V(i,j,vx)*U(i,j,bx)+V(i,j,vy)*U(i,j,by)+V(i,j,vz)*U(i,j,bz))**2)/tmpvar
           tmpvar=V(i,j,vx)*V(i,j,vx)+V(i,j,vy)*V(i,j,vy)+V(i,j,vz)*V(i,j,vz)
           vperp=sqrt(max(0.d0,tmpvar-vpara2))
           write(100) real(x(i)),real(y(j)),real(vperp)
        endif
     end do
  end do
  close(100)

  !IDL用の簡単な出力
  if (TYPE .eq. "roidl") then !IDL用密度
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(U(:,:,ro))
     close(100)
  elseif (TYPE .eq. "pridl") then !IDL用圧力
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,pr))
     close(100)
  elseif (TYPE .eq. "ptidl") then !IDL用全圧
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,pr)+(U(:,:,bx)*U(:,:,bx)+U(:,:,by)*U(:,:,by)+U(:,:,bz)*U(:,:,bz))/2.d0)
     close(100)
  elseif ((TYPE .eq. "pmidl").or.(TYPE .eq. "enemidl")) then !IDL用磁気圧・磁気エネルギー
     open(100,file=OUTPUT,form='unformatted')
     write(100) real((U(:,:,bx)*U(:,:,bx)+U(:,:,by)*U(:,:,by)+U(:,:,bz)*U(:,:,bz))/2.d0)
     close(100)
  elseif (TYPE .eq. "enekinidl") then !IDL用運動エネルギー
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(U(:,:,ro)/2.0*(V(:,:,vx)*V(:,:,vx)+V(:,:,vy)*V(:,:,vy)+V(:,:,vz)*V(:,:,vz)))
     close(100)
  elseif (TYPE .eq. "eneintidl") then !IDL用内部エネルギー
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(1.0/(gamma-1.0)*V(:,:,pr))
     close(100)
  elseif (TYPE .eq. "vxidl") then !IDL用vx
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,vx))
     close(100)
  elseif (TYPE .eq. "vyidl") then !IDL用vy
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,vy))
     close(100)
  elseif (TYPE .eq. "vzidl") then !IDL用vz
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,vz))
     close(100)
  elseif (TYPE .eq. "vmagidl") then !IDL用vmag
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(sqrt(V(:,:,vx)*V(:,:,vx)+V(:,:,vy)*V(:,:,vy)+V(:,:,vz)*V(:,:,vz)))
     close(100)
  elseif (TYPE .eq. "divvidl") then !IDL用divv
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(ARY0(:,:))
     close(100)
  elseif (TYPE .eq. "vnroidl") then !IDL用vnro
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(ARY0(:,:))
     close(100)
  elseif (TYPE .eq. "vnenidl") then !IDL用vnen
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(ARY0(:,:))
     close(100)
  elseif (TYPE .eq. "bxidl") then !IDL用bx
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(U(:,:,bx))
     close(100)
  elseif (TYPE .eq. "byidl") then !IDL用by
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(U(:,:,by))
     close(100)
  elseif (TYPE .eq. "bzidl") then !IDL用bz
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(U(:,:,bz))
     close(100)
  elseif (TYPE .eq. "jzidl") then !IDL用jz
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(jz_(:,:))
     close(100)
  elseif (TYPE .eq. "Azidl") then !IDL用Az
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(ARY0(:,:))
     close(100)
  elseif (TYPE .eq. "Tidl") then !温度
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,pr)/(U(:,:,ro)+1e-10))
     close(100)
  elseif (TYPE .eq. "enidl") then !エントロピー
     open(100,file=OUTPUT,form='unformatted')
     write(100) real(V(:,:,pr)/(U(:,:,ro)**(gamma)+1e-10))
     close(100)
  end if

!  if (argc .eq. 8) then !crop有り
!     !2Dグラフはこの順番にしないとlinesで線にならないっぽい
!     !splotで断面図はグラフの目盛数字と凡例の位置調整がうまくいかないのでやめた
!     !x profileのためにはループをこの順番にしないといけない?
!     if (xyswitch .eq. 0) then
!        do j=js,je
!           do i=js,ie
!              write(*,*) x(i),y(j),U(i,j,mx),U(i,j,my),U(i,j,mz),U(i,j,en),&
!                   U(i,j,bx),U(i,j,by),U(i,j,bz),U(i,j,ps),&
!                   jz_(i,j),ARY0(i,j),V(i,j,pr)/(U(i,j,ro)+1e-10)
!           end do
!        end do
!     end if
!     !y profileのためにはループをこの順番にしないといけない?
!     if (xyswitch .eq. 1) then
!        do i=is,ie        
!           do j=js,je
!              write(*,*) x(i),y(j),U(i,j,mx),U(i,j,my),U(i,j,mz),U(i,j,en),&
!                   U(i,j,bx),U(i,j,by),U(i,j,bz),U(i,j,ps),&
!                   jz_(i,j),ARY0(i,j),V(i,j,pr)/(U(i,j,ro)+1e-10)
!           end do
!        end do
!     end if
!  end if
!  stop

end program gp_convert_wada
