!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_ReductionTest

!----------------------------------------------------------------------
!
!  this program tests POP global reduction operations, including
!  global sums, minval, maxval, etc.  This test code relies on the
!  redistribution mod working correctly for gather/scatter operations.
!
!----------------------------------------------------------------------

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_IOUnitsMod
   use POP_ConfigMod

   implicit none

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: &
      configUnit,      &! unit assigned to configuration file
      errorCode         ! error flag

   integer (POP_i4) :: & 
      i4Test, i4Expect

   real (POP_r4) :: &
      r4Test, r4Expect  ! test values

   real (POP_r8) :: &
      r8Test, r8Expect  ! test values

   logical (POP_logical) :: &
      logTest, logExpect

   character (POP_charLength) :: &
      charTest, charExpect

!----------------------------------------------------------------------
!
!  initialize communcation environment
!
!----------------------------------------------------------------------

   call POP_CommInitMessageEnvironment
   call POP_CommInit

   errorCode = POP_Success

!----------------------------------------------------------------------
!
!  set test and correct values
!
!----------------------------------------------------------------------

   i4Test     = -2345_POP_i4
   i4Expect   = 1234567890_POP_i4
   r4Test     = -2345._POP_r4
   r4Expect   = 123456._POP_r4
   r8Test     = -2345._POP_r8
   r8Expect   = 12345678.12345678_POP_r8
   logTest    = .false.
   logExpect  = .true.
   charTest   = 'unknownValue'
   charExpect = 'correctValue'

!----------------------------------------------------------------------
!
!  open default config file
!
!----------------------------------------------------------------------

   call POP_ConfigOpen(configUnit, errorCode)

   if (errorCode /= POP_Success) then
     call POP_ErrorSet(errorCode, &
        'ConfigTest: error opening config file')
   endif

!----------------------------------------------------------------------
!
!  read variables and test results
!
!----------------------------------------------------------------------

   write(POP_stdout,POP_delimFormat)
   write(POP_stdout,POP_blankFormat)
   write(POP_stdout,'(a)') 'Testing reads from correct input file'
   write(POP_stdout,'(a)') 'Should see no errors'
   write(POP_stdout,POP_blankFormat)
   write(POP_stdout,POP_delimFormat)

   call POP_ConfigRead(configUnit, 'testModule', 'iVariable', &
                       i4Test, i4Expect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading integer variable')
   endif

   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for integer variable')
   endif

   call POP_ConfigRead(configUnit, 'testModule', 'rVariable', &
                       r4Test, r4Expect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading real variable')
   endif

   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for real variable')
   endif

   call POP_ConfigRead(configUnit, 'testModule', 'dVariable', &
                       r8Test, r8Expect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading double variable')
   endif

   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for double variable')
   endif

   call POP_ConfigRead(configUnit, 'testModule', 'lVariable', &
                       logTest, logExpect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading logical variable')
   endif

   if (logTest /= logExpect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for logical variable')
   endif

   call POP_ConfigRead(configUnit, 'testModule', 'cVariable', &
                       charTest, charExpect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading character variable')
   endif

   if (charTest /= charExpect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for character variable')
   endif

!----------------------------------------------------------------------
!
!  close config file
!
!----------------------------------------------------------------------

   call POP_ConfigClose(configUnit, errorCode)

   if (errorCode /= POP_Success) then
     call POP_ErrorSet(errorCode, &
        'ConfigTest: error closing config file')
   endif

!----------------------------------------------------------------------
!
!  repeat by reading from alternate input file
!
!----------------------------------------------------------------------

   call POP_ConfigOpen(configUnit, errorCode, &
                       configFileName='altInFile')

   if (errorCode /= POP_Success) then
     call POP_ErrorSet(errorCode, &
        'ConfigTest: error opening alternate config file')
   endif

!----------------------------------------------------------------------
!
!  read variables and test results
!
!----------------------------------------------------------------------

   write(POP_stdout,POP_delimFormat)
   write(POP_stdout,POP_blankFormat)
   write(POP_stdout,'(a)') 'Testing reads from bad input file'
   write(POP_stdout,'(a)') 'Should see no errors for ints and reals'
   write(POP_stdout,'(a)') 'Logical should end up with default value'
   write(POP_stdout,'(a)') 'Character will fail with module not found'
   write(POP_stdout,POP_blankFormat)
   write(POP_stdout,POP_delimFormat)

   i4Test     = -2345_POP_i4
   r4Test     = -2345._POP_r4
   r8Test     = -2345._POP_r8
   logTest    = .false.
   charTest   = 'unknownValue'

   !*** use integer to test output strings
   call POP_ConfigRead(configUnit, 'testModule', 'iVariable',     &
                       i4Test, i4Expect, errorCode,               &
                       outStringBefore = 'The integer value is:', &
                       outStringAfter  = ' units')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading integer variable')
   endif

   !*** use real to test only before string
   if (i4Test /= i4Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for integer variable')
   endif

   call POP_ConfigRead(configUnit, 'testModule', 'rVariable', &
                       r4Test, r4Expect, errorCode,           &
                       outStringBefore = 'The real value is:')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading real variable')
   endif

   if (r4Test /= r4Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for real variable')
   endif

   !*** use double to test only after string
   call POP_ConfigRead(configUnit, 'testModule', 'dVariable', &
                       r8Test, r8Expect, errorCode,           &
                       outStringAfter  = ' units')

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading double variable')
   endif

   if (r8Test /= r8Expect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for double variable')
   endif

   !*** test missing variable name
   call POP_ConfigRead(configUnit, 'testModule', 'logicalVariable', &
                       logTest, logExpect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading logical variable')
   endif

   if (logTest /= logExpect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for logical variable')
   endif

   !*** test missing module name
   call POP_ConfigRead(configUnit, 'nonexistantModule', 'cVariable', &
                       charTest, charExpect, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: error reading character variable')
   endif

   if (charTest /= charExpect) then
      call POP_ErrorSet(errorCode, &
         'ConfigTest: bad value for character variable')
   endif

!----------------------------------------------------------------------
!
!  close config file
!
!----------------------------------------------------------------------

   call POP_ConfigClose(configUnit, errorCode)

   if (errorCode /= POP_Success) then
     call POP_ErrorSet(errorCode, &
        'ConfigTest: error closing alternate config file')
   endif

!----------------------------------------------------------------------
!
!  clean up
!
!----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode)
   call POP_CommExitMessageEnvironment

!----------------------------------------------------------------------

 end program POP_ReductionTest

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
