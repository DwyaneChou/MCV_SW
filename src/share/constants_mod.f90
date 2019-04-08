MODULE constants_mod
    
  implicit none
  INTEGER,PARAMETER :: DOF       = 4    ! Degree of Freedoms within a 1D element
  INTEGER,PARAMETER :: numVAR    = 3    ! numbers of variables
  
  REAL,PARAMETER    :: gravity   = 9.80616
  REAL,PARAMETER    :: P0        = 1.E5
  REAL,PARAMETER    :: GAMMA     = 1.4
  REAL,PARAMETER    :: Rd        = 287.
  REAL,PARAMETER    :: Cv        = 717.5
  REAL,PARAMETER    :: Cp        = 1004.5
  REAL,PARAMETER    :: pi        = 2.*asin(1.)
  
  REAL,PARAMETER    :: radius    = 6371220.
  REAL,PARAMETER    :: D2R       = PI/180.    ! convert degree into radian
  REAL,PARAMETER    :: R2D       = 180./PI    ! convert radian into degree
  REAL,PARAMETER    :: OMEGA     = 7.292E-5
  REAL,PARAMETER    :: ALPHA_ROT = 0.         ! the flow orientation angle

END MODULE constants_mod
