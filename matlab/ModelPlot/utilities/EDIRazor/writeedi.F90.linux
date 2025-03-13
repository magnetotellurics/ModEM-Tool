#include "fintrf.h"
!----------------------------------------------------------------------!
!--------------------------------------------------------! MEXFUNCTION
!----------------------------------------------------------------------!

	subroutine mexFunction(nlhs,plhs,nrhs,prhs)

!----------------------------------------------------------------------!
!                                                                      !
!  Gateway routine for writing *.edi file from matlab by calling PLTIO !
!  module                                                              !
!                                                                      !
!  See also XXX.                                                       !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com                                  !
!  Started       : 2014/02/15                                          !
!  Last modified : 2014/02/15 17:12:52.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2014-02-15 Created.                                                 !
!                                                                      !
!----------------------------------------------------------------------!

	use PLTIO, only: PreWritePLT,WritePLT

	implicit none
!
!  mexFunction arguments.
!
	mwPointer plhs(*), prhs(*)
	integer(4)   nlhs, nrhs
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
	mwSize    mrows, ncols
!
!  Arugments for fortran routines.
!
	character(80) :: EDIFILE
	character(10) :: KW
	character(30) :: DATAID
	character(20) :: temp
	integer(4)    :: nf
	real(8)       :: x,y,z
	integer(4)    :: kf
	integer(4)    :: err
	real(8),allocatable,dimension(:) :: val

	character(90) :: msg

! -Variables- !--------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Check for proper number of arguments.
!
	if(nrhs .ne. 3 .AND. nrhs .ne. 7) then
		msg = '3/7 inputs required, the EDI file name and the keyword!'
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
!  Write out the head information or the keyword information of EDI file.
!
	select case (KW)
	case ('head','HEAD','>head','>HEAD')
		status = mxGetString(prhs(3),DATAID,30)
		!call mxCopyPtrToInteger4(mxGetPr(prhs(4)),nf,1)
		call mxCopyPtrToReal8(mxGetPr(prhs(4)),x,1)
		nf = int(x)
		call mxCopyPtrToReal8(mxGetPr(prhs(5)),x,1)
		call mxCopyPtrToReal8(mxGetPr(prhs(6)),y,1)
		call mxCopyPtrToReal8(mxGetPr(prhs(7)),z,1)
		call PreWritePLT(EDIFILE,DATAID,nf,x,y,z)
	case default
		mrows = mxGetM(prhs(3))
		ncols = mxGetN(prhs(3))
		nf    = mrows*ncols
		allocate(val(nf))
		call mxCopyPtrToReal8(mxGetPr(prhs(3)),val,nf)
		call WritePLT(EDIFILE,nf,adjustl(trim(kw)),val)
		deallocate(val)
	end select ! kw
!
!  End of the subroutine.
!
	return
	end subroutine mexFunction

