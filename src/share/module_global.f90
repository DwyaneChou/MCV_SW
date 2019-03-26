MODULE global

  INTEGER,PARAMETER    :: NP=3        ! numbers of DOFS within a 1D element
  INTEGER,PARAMETER    :: NumVAR=5    ! numbers of variables
  INTEGER,PARAMETER    :: NumMOI=5    ! numbers of moisture variables

  REAL,PARAMETER       :: GRA=9.80616
  REAL,PARAMETER       :: P0=1.E5
  REAL,PARAMETER       :: GAMMA=1.4
  REAL,PARAMETER       :: RD=287.
  REAL,PARAMETER       :: CV=717.5
  REAL,PARAMETER       :: CP=1004.5
  REAL,PARAMETER       :: PI=2.*asin(1.)

#ifdef SPHERE
  REAL,PARAMETER       :: RA=6371220.
  REAL,PARAMETER       :: D2R=PI/180.    ! convert degree into radian
  REAL,PARAMETER       :: R2D=180./PI    ! convert radian into degree
  REAL,PARAMETER       :: OMEGA=7.292E-5
  !REAL,PARAMETER       :: OMEGA=0.       ! No rotation
  REAL,PARAMETER       :: ALPHA_ROT=0.   ! the flow orientation angle
  !REAL,PARAMETER       :: ALPHA_ROT=PI/4.   ! the flow orientation angle
  !REAL,PARAMETER       :: ALPHA_ROT=PI/2.   ! the flow orientation angle
#endif

#ifdef CARTESIAN
  REAL,PARAMETER       :: RA =1.
  REAL,PARAMETER       :: D2R=1.    ! convert degree into radian
  REAL,PARAMETER       :: R2D=1.    ! convert radian into degree
  REAL,PARAMETER       :: OMEGA=0.
#endif


END MODULE global
