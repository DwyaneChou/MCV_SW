MODULE projection_mod

  USE constants_mod ,only: radius, PI

  implicit none

contains

  subroutine pointProjPlane2Sphere(lambda,theta,x,y,k)

    implicit none

    integer,intent(in ) :: k            ! Patch index
    real   ,intent(in ) :: x,y          ! local coordinate on patch
    real   ,intent(out) :: lambda,theta ! spherical coordinate
    
    ! Local variables
    real :: a,b,r

    select case(k)      
    case(1:4)
      lambda = x + dble(k - 1) * pi / 2.
      theta  = atan2( tan(y) * cos(x), 1. )
    case(5)
      a = tan(x)
      b = tan(y)
      
      lambda = atan2(a , -b           )
      theta  = atan2(1., sqrt(a*a+b*b))
    case(6)
      a = tan(x)
      b = tan(y)
      
      lambda =  atan2(a , b            )
      theta  = -atan2(1., sqrt(a*a+b*b))  
    end select

    return
  end subroutine pointProjPlane2Sphere

  subroutine pointProjSphere2Plane(x,y,lambda,theta,k)

    implicit none

    integer, intent(in ) :: k            ! Patch index
    real   , intent(in ) :: lambda,theta ! spherical coordinate
    real   , intent(out) :: x,y          ! local coordinate on patch

    select case(k)
    case(1:4)
      x = atan(tan(lambda-float(k-1)*pi/2.))
      y = atan(tan(theta)/cos(lambda-float(k-1)*pi/2.))
    case(5:6)
      x = atan((-1.)**(k+1)*sin(lambda)/tan(theta))
      y = atan(-cos(lambda)/tan(theta))
    end select

    return
  end subroutine pointProjSphere2Plane

  subroutine contravProjSphere2Plane(contrav1,contrav2,sv1,sv2,k,lambda,theta)

    implicit none

    integer,intent(in) :: k
    real,intent(out)   :: contrav1,contrav2
    real,intent(in)    :: sv1,sv2,lambda,theta
  ! Local variables
    real    :: ia(2,2)
  ! -------------------------

    call matrixIA(ia,k,lambda,theta)
          
    contrav1=ia(1,1)*sv1+ia(1,2)*sv2
    contrav2=ia(2,1)*sv1+ia(2,2)*sv2

    return
  end subroutine contravProjSphere2Plane

  subroutine contravProjPlane2Sphere(sv1,sv2,contrav1,contrav2,k,lambda,theta)

    implicit none

    integer k
    real,intent(in)       :: contrav1,contrav2,lambda,theta
    real,intent(out)      :: sv1,sv2
  ! Local variables
    real    :: a(2,2)
  ! -------------------------
          
    call matrixA(a,k,lambda,theta)
    
    sv1 = a(1,1) * contrav1 + a(1,2) * contrav2
    sv2 = a(2,1) * contrav1 + a(2,2) * contrav2

    return
  end subroutine contravProjPlane2Sphere

  subroutine matrixIA(ima,k,lambda,theta)

    implicit none

    integer,intent(in) :: k     
    real,intent(in)    :: lambda,theta
    real,intent(out)   :: ima(2,2)
  ! Local variables
    real   :: alambda,atheta,a,b,c,d,temp
  ! -------------------------
    
    ima = 0.

    if (k <= 4) then
      alambda=lambda-dble(k-1)*pi/2.
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=d*d*b*b+c*c
      ima(1,1)=1./d
      ima(1,2)=0.
      ima(2,1)=a*c/temp
      ima(2,2)=b/temp
    else if (k==5) then
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=c+a*a*d*d/c
      ima(1,1)=b/temp
      ima(1,2)=-a/c/temp
      temp=c+b*b*d*d/c
      ima(2,1)=a/temp
      ima(2,2)=b/c/temp
    else
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=c+a*a*d*d/c
      ima(1,1)=-b/temp
      ima(1,2)=a/c/temp
      temp=c+b*b*d*d/c
      ima(2,1)=a/temp
      ima(2,2)=b/c/temp
    endif

    ima=ima/radius
          
    return
  end subroutine matrixIA

  subroutine matrixA(ma,k,lambda,theta)

    implicit none

    integer,intent(in)    :: k
    real,intent(in)       :: lambda,theta
    real,intent(out)      :: ma(2,2)
  ! Local variables
    real     :: alambda,atheta,a,b,c,d,temp,r
  ! -------------------------
    
    ma = 0.d0

    if (k <= 4) then
      alambda=lambda-dble(k-1)*pi/2.
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      ma(1,1)=d
      ma(1,2)=0.d0
      ma(2,1)=-c*d*a/b
      ma(2,2)=b*d*d+c*c/b
    else if (k==5) then
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=1.+a*a*d*d/c/c
      ma(1,1)=b*c*temp
      ma(2,1)=-c*c*a*temp
      temp=1.+b*b*d*d/c/c
      ma(1,2)=a*c*temp
      ma(2,2)=b*c*c*temp
    else
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=1.+a*a*d*d/c/c
      ma(1,1)=-b*c*temp
      ma(2,1)=c*c*a*temp
      temp=1.+b*b*d*d/c/c
      ma(1,2)=a*c*temp
      ma(2,2)=b*c*c*temp
    endif

    ma = ma*radius

    return
  end subroutine matrixA

  subroutine matrixG(mg,x,y)

    implicit none

    real,intent(in)    :: x,y
    real,intent(out)   :: mg(2,2)
  ! Local variables
    real     :: rho
  ! -------------------------
          
    rho=sqrt(1.+tan(x)**2.+tan(y)**2.)

    mg(1,1) = (1.+tan(x)**2.)
    mg(1,2) = -tan(x)*tan(y)
    mg(2,1) = mg(1,2)
    mg(2,2) = (1.+tan(y)**2.)
    
    mg      = mg/(rho**4.*cos(x)**2.*cos(y)**2.)

    mg = radius**2 * mg

    return
  end subroutine matrixG

  subroutine matrixIG(img,x,y)

    implicit none

    real,intent(in)          :: x,y
    real,intent(out)         :: img(2,2)
  ! Local variables
    real                     :: rho
  ! -------------------------
    
    img = 0.

    rho=sqrt(1.+tan(x)**2.+tan(y)**2.)

    img(1,1)=1.+tan(y)**2.
    img(1,2)=tan(x)*tan(y)
    img(2,1)=img(1,2)
    img(2,2)=1.+tan(x)**2.

    img=img/(radius**2)*(rho**2.*cos(x)**2.*cos(y)**2.)

    return
  end subroutine matrixIG

  subroutine computeJacobian(jab,x,y)

    implicit none

    real,intent(in)      :: x,y
    real,intent(out)     :: jab
  ! Local variables
    real                 :: rho
  ! -----------------------

    rho = sqrt(1 + tan(x)**2 + tan(y)**2)
    
    jab = radius**2 /( cos(x)**2 * cos(y)**2 * rho**3 )

    return
  end subroutine computeJacobian

END MODULE projection_mod

