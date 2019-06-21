MODULE projection_mod
  use constants_mod
  use parameters_mod
  implicit none

contains

  subroutine pointProjPlane2Sphere(lambda,theta,x,y,iPatch)

    implicit none

    integer,intent(in ) :: iPatch       ! Patch index
    real   ,intent(in ) :: x,y          ! local coordinate on patch in degree
    real   ,intent(out) :: lambda,theta ! spherical coordinate
    
    ! Local variables
    real :: a,b

    select case(iPatch)      
    case(1:4)
      lambda = x + dble(iPatch - 1) * pi / 2.
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

  subroutine pointProjSphere2Plane(x,y,lambda,theta,iPatch)

    implicit none

    integer, intent(in ) :: iPatch            ! Patch index
    real   , intent(in ) :: lambda,theta ! spherical coordinate
    real   , intent(out) :: x,y          ! local coordinate on patch

    select case(iPatch)
    case(1:4)
      x = atan(tan(lambda-float(iPatch-1)*pi/2.))
      y = atan(tan(theta)/cos(lambda-float(iPatch-1)*pi/2.))
    case(5:6)
      x = atan((-1.)**(iPatch+1)*sin(lambda)/tan(theta))
      y = atan(-cos(lambda)/tan(theta))
    end select

    return
  end subroutine pointProjSphere2Plane
  
  subroutine covProjSphere2Plane(cov1, cov2, sv1, sv2, matrixIA, matrixG)

    implicit none
    
    real, intent(out) :: cov1,cov2
    real, intent(in ) :: sv1,sv2
    real, intent(in ) :: matrixIA(2,2)
    real, intent(in ) :: matrixG (2,2)
    
    real matrix(2,2)
    
    matrix = matmul(matrixG,matrixIA)
          
    cov1 = matrix(1,1)*sv1 + matrix(1,2)*sv2
    cov2 = matrix(2,1)*sv1 + matrix(2,2)*sv2

    return
  end subroutine covProjSphere2Plane
  
  subroutine covProjPlane2Sphere(sv1, sv2, cov1, cov2, matrixA, matrixIG)

    implicit none

    real, intent(out) :: sv1,sv2
    real, intent(in)  :: cov1,cov2
    real, intent(in)  :: matrixA (2,2)
    real, intent(in)  :: matrixIG(2,2)
    
    real matrix(2,2)
    
    matrix = matmul(matrixA,matrixIG)
    
    sv1 = matrix(1,1) * cov1 + matrix(1,2) * cov2
    sv2 = matrix(2,1) * cov1 + matrix(2,2) * cov2

    return
  end subroutine covProjPlane2Sphere
  
  subroutine contravProjSphere2Plane(contrav1, contrav2, sv1, sv2, matrixIA)

    implicit none
    
    real   ,intent(in)  :: sv1,sv2
    real   ,intent(in)  :: matrixIA(2,2)
    real   ,intent(out) :: contrav1,contrav2
          
    contrav1 = matrixIA(1,1)*sv1 + matrixIA(1,2)*sv2
    contrav2 = matrixIA(2,1)*sv1 + matrixIA(2,2)*sv2

    return
  end subroutine contravProjSphere2Plane

  subroutine contravProjPlane2Sphere(sv1, sv2, contrav1, contrav2, matrixA)

    implicit none

    real,intent(in)  :: contrav1,contrav2
    real,intent(in)  :: matrixA(2,2)
    real,intent(out) :: sv1,sv2
    
    sv1 = matrixA(1,1) * contrav1 + matrixA(1,2) * contrav2
    sv2 = matrixA(2,1) * contrav1 + matrixA(2,2) * contrav2

    return
  end subroutine contravProjPlane2Sphere

  subroutine contrav2cov(cov1,cov2,contrav1,contrav2, matrixG)
    real, intent(in ) :: contrav1
    real, intent(in ) :: contrav2
    real, intent(in ) :: matrixG(2,2)
    real, intent(out) :: cov1
    real, intent(out) :: cov2
    
    cov1 = matrixG(1,1) * contrav1 + matrixG(1,2) * contrav2
    cov2 = matrixG(2,1) * contrav1 + matrixG(2,2) * contrav2
    
  end subroutine contrav2cov
  
  subroutine cov2contrav(contrav1,contrav2,cov1,cov2,matrixIG)
    real, intent(in ) :: cov1
    real, intent(in ) :: cov2
    real, intent(in ) :: matrixIG(2,2)
    real, intent(out) :: contrav1
    real, intent(out) :: contrav2
    
    contrav1 = matrixIG(1,1) * cov1 + matrixIG(1,2) * cov2
    contrav2 = matrixIG(2,1) * cov1 + matrixIG(2,2) * cov2
    
  end subroutine cov2contrav
  
  subroutine calc_matrixG(matrixG,x,y)

    implicit none

    real,intent(in)    :: x,y
    real,intent(out)   :: matrixG(2,2)
  ! Local variables
    real     :: rho
  ! -------------------------
          
    rho = sqrt(1.+tan(x)**2.+tan(y)**2.)

    matrixG(1,1) = 1.+tan(x)**2.
    matrixG(1,2) = -tan(x)*tan(y)
    matrixG(2,1) = matrixG(1,2)
    matrixG(2,2) = 1.+tan(y)**2.
    
    matrixG      = radius**2 / (rho**4. * cos(x)**2. * cos(y)**2.) * matrixG

    return
  end subroutine calc_matrixG

  subroutine calc_matrixIG(matrixIG,x,y)

    implicit none

    real,intent(in)          :: x,y
    real,intent(out)         :: matrixIG(2,2)
  ! Local variables
    real                     :: rho
  ! -------------------------
    
    rho = sqrt(1.+tan(x)**2.+tan(y)**2.)

    matrixIG(1,1)=1.+tan(y)**2.
    matrixIG(1,2)=tan(x)*tan(y)
    matrixIG(2,1)=matrixIG(1,2)
    matrixIG(2,2)=1.+tan(x)**2.

    matrixIG = (rho**2.*cos(x)**2.*cos(y)**2.) / radius**2 * matrixIG

    return
  end subroutine calc_matrixIG

  subroutine calc_Jacobian(jab,x,y)

    implicit none

    real,intent(in)      :: x,y
    real,intent(out)     :: jab
  ! Local variables
    real                 :: rho
  ! -----------------------

    rho = sqrt(1 + tan(x)**2 + tan(y)**2)
    
    jab = 1. /( cos(x)**2 * cos(y)**2 * rho**3 )
    
    jab = jab * radius**2

    return
  end subroutine calc_Jacobian

  subroutine calc_matrixIA(matrixIA, lambda, theta, iPatch)

    implicit none

    integer,intent(in) :: iPatch     
    real,intent(in)    :: lambda,theta
    real,intent(out)   :: matrixIA(2,2)
  ! Local variables
    real   :: alambda,atheta,a,b,c,d,temp
  ! -------------------------
    
    if (iPatch <= 4) then
      alambda=lambda-dble(iPatch-1)*pi/2.
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=d*d*b*b+c*c
      matrixIA(1,1)=1./d
      matrixIA(1,2)=0.
      matrixIA(2,1)=a*c/temp
      matrixIA(2,2)=b/temp
    else if (iPatch==5) then
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=c+a*a*d*d/c
      matrixIA(1,1)=b/temp
      matrixIA(1,2)=-a/c/temp
      temp=c+b*b*d*d/c
      matrixIA(2,1)=a/temp
      matrixIA(2,2)=b/c/temp
    else
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=c+a*a*d*d/c
      matrixIA(1,1)=-b/temp
      matrixIA(1,2)=a/c/temp
      temp=c+b*b*d*d/c
      matrixIA(2,1)=a/temp
      matrixIA(2,2)=b/c/temp
    endif

    matrixIA = matrixIA/radius
          
    return
  end subroutine calc_matrixIA

  subroutine calc_matrixA(matrixA, lambda, theta, iPatch)

    implicit none

    integer,intent(in)    :: iPatch
    real,intent(in)       :: lambda,theta
    real,intent(out)      :: matrixA(2,2)
  ! Local variables
    real     :: alambda,atheta,a,b,c,d,temp,r
  ! -------------------------

    if (iPatch <= 4) then
      alambda=lambda-dble(iPatch-1)*pi/2.
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      matrixA(1,1)=d
      matrixA(1,2)=0.d0
      matrixA(2,1)=-c*d*a/b
      matrixA(2,2)=b*d*d+c*c/b
    else if (iPatch==5) then
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=1.+a*a*d*d/c/c
      matrixA(1,1)=b*c*temp
      matrixA(2,1)=-c*c*a*temp
      temp=1.+b*b*d*d/c/c
      matrixA(1,2)=a*c*temp
      matrixA(2,2)=b*c*c*temp
    else
      alambda=lambda
      atheta=theta
      a=sin(alambda)
      b=cos(alambda)
      c=sin(atheta)
      d=cos(atheta)
      temp=1.+a*a*d*d/c/c
      matrixA(1,1)=-b*c*temp
      matrixA(2,1)=c*c*a*temp
      temp=1.+b*b*d*d/c/c
      matrixA(1,2)=a*c*temp
      matrixA(2,2)=b*c*c*temp
    endif

    matrixA = matrixA*radius

    return
  end subroutine calc_matrixA
  
END MODULE projection_mod

