MODULE constants
    
  implicit none
  INTEGER,PARAMETER :: DOF       = 3    ! Degree of Freedoms within a 1D element
  INTEGER,PARAMETER :: NumVAR    = 3    ! numbers of variables
  
  REAL,PARAMETER    :: g         = 9.80616d0
  REAL,PARAMETER    :: P0        = 1.E5
  REAL,PARAMETER    :: GAMMA     = 1.4d0
  REAL,PARAMETER    :: Rd        = 287.d0
  REAL,PARAMETER    :: Cv        = 717.5d0
  REAL,PARAMETER    :: Cp        = 1004.5d0
  REAL,PARAMETER    :: pi        = 2.d0*asin(1.d0)
  
  REAL,PARAMETER    :: radius    = 6371220.d0
  REAL,PARAMETER    :: D2R       = PI/180.d0    ! convert degree into radian
  REAL,PARAMETER    :: R2D       = 180.d0/PI    ! convert radian into degree
  REAL,PARAMETER    :: OMEGA     = 7.292E-5
  REAL,PARAMETER    :: ALPHA_ROT = 0.d0         ! the flow orientation angle

END MODULE constants
