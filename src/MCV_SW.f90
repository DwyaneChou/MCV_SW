    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multi-moment Constrained finite Volume shallow water dynamic core !
    ! Written by Lilong Zhou, March 26, 2019                            !
    ! Numerical Weather Prediction Center of CMA                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    program MCV_SW
      use parameters_mod, only: initParameters
      use stat_mod      , only: initStat
      use tend_mod      , only: initTend
      use mesh_mod      , only: initMesh
      use test_case_mod , only: initTestCase
      use output_mod    , only: write_netCDF
      implicit none
      
      integer :: timeStart,timeEnd
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      call initParameters
      call initMesh
      call initStat
      call initTend
      call initTestCase
      
      call write_netCDF
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',dble(timeEnd-timeStart)/10000.0,' seconds to run this program'
    end program MCV_SW
