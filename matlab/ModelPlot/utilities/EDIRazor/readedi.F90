#include "fintrf.h"
!----------------------------------------------------------------------!
!--------------------------------------------------------! MEXFUNCTION
!----------------------------------------------------------------------!

	subroutine mexFunction(nlhs,plhs,nrhs,prhs)

!----------------------------------------------------------------------!
!                                                                      !
!  Gateway routine for reading *.edi file from matlab by calling PLTIO !
!  module.                                                             !
!                                                                      !
!  See also writeedi.                                                  !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com                                  !
!  Started       : 2014/02/14                                          !
!  Last modified : 2014/02/17 16:09:46.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2014-02-14 Created.                                                 !
!  2014-02-17 Added a step for checking the file existed or not.       !
!                                                                      !
!----------------------------------------------------------------------!

	use PLTIO, only: PreReadPLT,ReadPLT,FindPLTKW

	implicit none
!
!  mexFunction arguments.
!
	mwPointer plhs(*), prhs(*)
	integer(4) nlhs, nrhs
!
!  API function declarations.
!
	mwPointer mxGetPr
	mwPointer mxCreateDoubleMatrix
	mwSize    mxGetM, mxGetN
	mwPointer mxCreateString, mxGetString
	integer   status
!
!  Pointers and arguments for input/output mxArrays:
!
	mwSize    maxbuf,kwbuf
	parameter(maxbuf=80,kwbuf=10)
!
!  Array information:
!
	mwSize    mrows, ncols, size
!
!  Arugments for fortran routines.
!
	character(80) :: EDIFILE
	character(10) :: KW
	character(30) :: DATAID
	integer(4)    :: nr
	integer(4)    :: nf
	real(8)       :: x,y,z
	integer(4)    :: kr,kf

	integer(4)    :: err
	real(8),allocatable,dimension(:) :: val

	character(90) :: msg
	logical       :: existed

! -Variables- !--------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Check for proper number of arguments.
!
	if(nrhs .ne. 2) then
		msg = 'Only two inputs required, the EDI file name and the keyword!'
		call mexErrMsgIdAndTxt('MATLAB:dellatlon:nInput',trim(msg))
	!elseif(nlhs .gt. 1) then
		!msg = 'Too many output arguments.'
		!call mexErrMsgIdAndTxt ('MATLAB:dellatlon:nOutput',trim(msg))
	endif
!
!  Get the EDI file name.
!
	status = mxGetString(prhs(1),EDIFILE,maxbuf)
	if (status .NE. 0) then
		msg = 'Error reading EDI file name!'
		call mexErrMsgIdAndTxt('MATLAB:readedi:readError',msg)
	endif ! status.
!
!  Get the keyword.
!
	status = mxGetString(prhs(2),KW,kwbuf)
	if (status .NE. 0) then
		msg = 'Error reading the keyword!'
		call mexErrMsgIdAndTxt('MATLAB:readedi:readError',msg)
	endif ! status.
!
!  See the EDI file exited or not.
!
	inquire(file=EDIFILE,exist=existed)
	if (.NOT. existed) then
		msg = 'EDI file:'//trim(EDIFILE)//' does not existed!'
		call mexErrMsgIdAndTxt('MATLAB:readedi:readError',msg)
	endif ! edi file not existed.
!
!  Read the plt file.
!
	call PreReadPLT(EDIFILE,DATAID,nf,x,y,z)
!
!  Write out the frequency information and the 1st reciever's info.
!
	select case (KW)
	case ('head','HEAD','>head','>HEAD')
		plhs(1) = mxCreateString(DATAID)
		plhs(2) = mxCreateDoubleMatrix(1,1,0)
		plhs(3) = mxCreateDoubleMatrix(1,1,0)
		plhs(4) = mxCreateDoubleMatrix(1,1,0)
		plhs(5) = mxCreateDoubleMatrix(1,1,0)
		call mxCopyReal8ToPtr(dble(nf),mxGetPr(plhs(2)),1)
		call mxCopyReal8ToPtr(x,       mxGetPr(plhs(3)),1)
		call mxCopyReal8ToPtr(y,       mxGetPr(plhs(4)),1)
		call mxCopyReal8ToPtr(z,       mxGetPr(plhs(5)),1)
	case default
		existed = FindPLTKW(EDIFILE,trim(KW))
		if (.NOT. existed) then
			msg = 'Keyword: "'//trim(KW)//'" is not existed in '//   &
				   trim(EDIFILE)//' or is a wrong keyword.'
			call mexErrMsgIdAndTxt('MATLAB:readedi:readError',msg)
		endif ! if kw not existed.
		allocate(val(nf))
		call ReadPLT(EDIFILE,nf,adjustl(trim(KW)),val)
		plhs(1) = mxCreateDoubleMatrix(nf,1,0)
		call mxCopyReal8ToPtr(val,mxGetPr(plhs(1)),nf)
		deallocate(val)
	end select ! kw
!
!  End of the subroutine.
!
	return
	end subroutine mexFunction
