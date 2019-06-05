MODULE constants_mod
    
  implicit none
  INTEGER,PARAMETER :: DOF       = 4    ! Degree of Freedoms within a 1D element
  INTEGER,PARAMETER :: numVAR    = 3    ! numbers of variables
  
  REAL,PARAMETER    :: gravity   = 9.80616
  REAL,PARAMETER    :: pi        = 2.*asin(1.)
  
  REAL,PARAMETER    :: radius    = 6371220.
  REAL,PARAMETER    :: D2R       = PI/180.    ! convert degree into radian
  REAL,PARAMETER    :: R2D       = 180./PI    ! convert radian into degree
  REAL,PARAMETER    :: Omega     = 7.292E-5

  REAL,PARAMETER    :: FillValue = -9999999999999999.
END MODULE constants_mod
