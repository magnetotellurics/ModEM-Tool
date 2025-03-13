!----------------------------------------------------------------------!
!--------------------------------------------------------------! PLTIO
!----------------------------------------------------------------------!

	module PLTIO

!----------------------------------------------------------------------!
!                                                                      !
!  Module to IO the MT data in *.plt format of PHEONIX or EDI format.  !
!                                                                      !
!  See also EDIIO.                                                     !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com                                  !
!  Started       : 2012/11/16                                          !
!  Last modified : 2014/02/15 17:41:03.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2012-11-16 Created.                                                 !
!  2012-11-18 Released v1.0.                                           !
!  2012-11-20 I) Added a function for finding a keyword.               !
!             II) Updated the rules for finding a keyword in a line.   !
!  2012-12-06 Added subroutine PreWritePLT2 to output the XYZ          !
!             coordinate.                                              !
!  2013-12-21 Added keywords for supporting Phase Tensor.              !
!  2013-12-23 I) Fixed a bug in function d2dms - Output *** when lat<0.!
!             II) Moved d2dms & dms2d function to module LATLON.       !
!  2014-01-10 I) Added function GetPLTVAL for extracting some value by !
!                keywords.                                             !
!             II) Added subroutine CompZerrFlr for computing error     !
!                 floor of impedance.                                  !
!  2014-02-14 Deleted subroutine chkkw for better compatibility of     !
!             calling from MATLAB.                                     !
!  2014-02-15 I) Added DATAID I/O in PreReadPLT & PreWritePLT.         !
!             II) Fixed a bug of keyword matching in FindPLTKW.        !
!  2014-02-15 Released v1.1.                                           !
!                                                                      !
!----------------------------------------------------------------------!

	use LATLON, only: d2dms,dms2d

	implicit none

	public  PreReadPLT,ReadPLT,PreWritePLT,WritePLT,FindPLTKW
	public  PreWritePLT2,GetPLTVAL,CompZerrFlr
	private getsval,dms2d,d2dms,chkkw

	private sValue
	character(180) :: sValue
!
!  The contained subroutines.
!
	contains

!----------------------------------------------------------------------!
!---------------------------------------------------------! PREREADPLT !
!----------------------------------------------------------------------!

	subroutine PreReadPLT(PLTFILE,DATAID,nf,lat,lon,elv)
!
!  Routine to get the number of frequency and location of the site.
!
	implicit none
	character(*) :: PLTFILE
	character(*) :: DATAID
	integer(4)   :: nf
	real(8)      :: lat,lon,elv
	integer(4)   :: err
	integer(4)   :: k,j
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='old',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Get the data id.
!
	call getsval(101,'DATAID=',sValue)
	DATAID = sValue(1:len(DATAID))
!
!  Remove the double quotes if exised.
!
	k = index(DATAID,'"')
	j = index(DATAID,'"',BACK=.TRUE.)
	DATAID = trim(adjustl(DATAID(k+1:j-1)))
!
!  Get number of frequency.
!
	call getsval(101,'NFREQ=',sValue)
	read(sValue,*) nf
!
!  Get latitude or x-coordinate.
!
	call getsval(101,'LAT=',sValue)
	if (index(sValue,':') .NE. 0) then
		lat = dms2d(sValue)
	else
		read(sValue,*) lat
	endif ! X^o X' X"
!
!  Get longitude or y-coordinate.
!
	call getsval(101,'LONG=',sValue)
	if (index(sValue,':') .NE. 0) then
		lon = dms2d(sValue)
	else
		read(sValue,*) lon
	endif ! X^o X' X"
!
!  Get elevation or z-coordinate.
!
	call getsval(101,'ELEV=',sValue)
	read(sValue,*) elv
!
!  Close the PLT file.
!
	close(101)
!
!  end of the subroutine.
!
	end subroutine PreReadPLT

!----------------------------------------------------------------------!
!------------------------------------------------------------! READPLT !
!----------------------------------------------------------------------!

	subroutine ReadPLT(PLTFILE,nf,kw,val)
!
!  Routine to read the value of the keyword.
!
	implicit none
	character(*)   :: PLTFILE
	character(*)   :: kw
	integer(4)     :: kf,nf
	real(8)        :: val(nf)
	integer(4)     :: err
	character(180) :: sLine
!
!  Check the validity of the keyword.
!
	!call chkkw(kw)
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='old',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Find the requested block.
!
	read(101,'(a)') sLine
	do while (index(sLine,kw) .EQ. 0 .OR. index(sLine,kw//'.') .NE. 0)
		read(101,'(a)',end=90) sLine
	enddo ! kw
	read(101,*) (val(kf), kf=1,nf)
!
!  Close the PLT file.
!
	close(101)
	return
!
!  Error information.
!
90    write(0,*) 'ERROR: have not found the keyword: ',kw
	write(0,*) '       STOP'
	stop
!
!  end of the subroutine.
!
	end subroutine ReadPLT

!----------------------------------------------------------------------!
!--------------------------------------------------------! PREWRITEPLT !
!----------------------------------------------------------------------!

	subroutine PreWritePLT(PLTFILE,DATAID,nf,lat,lon,elv)
!
!  Routine to write the number of frequency and location of the site to
!  PLTFILE.
!
	implicit none
	character(*) :: PLTFILE
	character(*) :: DATAID
	integer(4)   :: nf
	real(8)      :: lat,lon,elv
	character(23):: str
	integer(4)   :: err
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='replace',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Write out the parameters.
!
	write(101,'(a5)') '>HEAD'
	write(101,*     ) ''

	write(101,'(a)') '    DATAID="'//adjustl(trim(DATAID))//'"'

	str(1 :4 ) = '    '
	str(5 :8 ) = 'LAT='
	str(9 :22) = d2dms(lat)
	write(101,'(a22)') str

	str(5 :9 ) = 'LONG='
	str(10:23) = d2dms(lon)
	write(101,'(a23)') str

	str(5 :9 ) = 'ELEV='
	write(101,'(a9,f8.2)') str(1:9),elv

	write(101,*     ) ''
	write(101,'(a8)') '>=MTSECT'
	write(101,*     ) ''

	str(1:6) = 'NFREQ='
	write(101,'(a6,i4)') str(1:6),nf
!
!  Close the PLT file.
!
	close(101)
!
!  end of the subroutine.
!
	end subroutine PreWritePLT

!----------------------------------------------------------------------!
!-------------------------------------------------------! PREWRITEPLT2 !
!----------------------------------------------------------------------!

	subroutine PreWritePLT2(PLTFILE,nf,lat,lon,elv)
!
!  Routine to write the number of frequency and location of the site to
!  PLTFILE.
!
	implicit none
	character(*) :: PLTFILE
	integer(4)   :: nf
	real(8)      :: lat,lon,elv
	character(22):: str
	integer(4)   :: err
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='replace',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Write out the parameters.
!
	write(101,'(a5)') '>HEAD'
	write(101,*     ) ''

	str(1 :4 ) = '    '
	str(5 :8 ) = 'LAT='
!	write(101,'(a8,f12.2)') str(1:8),lat
	write(str(9:22),'(f12.2)') lat
	write(101,'(a8,a14)') str(1:8),adjustl(str(9:22))

	str(5 :9 ) = 'LONG='
!	write(101,'(a9,f12.2)') str(1:9),lon
	write(str(10:22),'(f12.2)') lon
	write(101,'(a9,a13)') str(1:9),adjustl(str(10:22))

	str(5 :9 ) = 'ELEV='
!	write(101,'(a9,f12.2)') str(1:9),elv
	write(str(10:22),'(f12.2)') elv
	write(101,'(a9,a13)') str(1:9),adjustl(str(10:22))

	write(101,*     ) ''
	write(101,'(a8)') '>=MTSECT'
	write(101,*     ) ''

	str(1:6) = 'NFREQ='
	write(101,'(a6,i4)') str(1:6),nf
!
!  Close the PLT file.
!
	close(101)
!
!  end of the subroutine.
!
	end subroutine PreWritePLT2

!----------------------------------------------------------------------!
!-----------------------------------------------------------! WRITEPLT !
!----------------------------------------------------------------------!

	subroutine WritePLT(PLTFILE,nf,kw,val)
!
!  Routine to write out the value of the keyword.
!
	implicit none
	character(*)   :: PLTFILE
	character(*)   :: kw
	integer(4)     :: kf,nf
	real(8)        :: val(nf)
	integer(4)     :: err
!
!  Check the validity of the keyword.
!
	!call chkkw(kw)
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='old',position='append',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Write out the keyword related parameters.
!
	write(101,* ) ''
	write(101,'(a,a4,i4)') kw,' // ',nf
	!write(101,* ) val
	write(101,'(7es16.6)') val
!
!  Close the PLT file.
!
	close(101)
!
!  end of the subroutine.
!
	end subroutine WritePLT

!----------------------------------------------------------------------!
!---------------------------------------------------------! COMPERRFLR !
!----------------------------------------------------------------------!

	subroutine CompZerrFlr(perc,method,zxx,zxy,zyx,zyy,dxx,dxy,dyx,dyy)
!
!  Routine to compute the error floor for MT data.
!  Method = 1 : perc*sqrt(|Zxy*Zyx|) for all 4 components.
!  Method = 2 : perc*|Zxy| for Zxx & Zxy, perc*|Zyx| for Zyx & Zyy.
!  Method = 3 : method 2 for Zxy & Zyx, method 1 for Zxx & Zyy.
!  Method = 4 : sqrt(det(Z)) for all 4 components.
!
	implicit none
	real(8)    :: perc
	integer(4) :: method
	complex(8) :: zxx,zxy,zyx,zyy
	real(8)    :: dxx,dxy,dyx,dyy
	real(8)    :: sd

	select case (method)
	case (1)
		sd  = perc * dsqrt(cdabs(zxy*zyx))
		dxx = max(dxx,sd)
		dxy = max(dxy,sd)
		dyx = max(dyx,sd)
		dyy = max(dyy,sd)
	case (2)
		sd  = perc * cdabs(zxy)
		dxx = max(dxx,sd)
		dxy = max(dxy,sd)
		sd  = perc * cdabs(zyx)
		dyx = max(dyx,sd)
		dyy = max(dyy,sd)
	case (3)
		sd  = perc * cdabs(zxy)
		dxy = max(dxy,sd)
		sd  = perc * cdabs(zyx)
		dyx = max(dyx,sd)
		sd  = perc * dsqrt(cdabs(zxy*zyx))
		dxx = max(dxx,sd)
		dyy = max(dyy,sd)
	case (4)
		sd  = perc * dsqrt(cdabs(zxx*zyy-zxy*zyx))
		dxx = max(dxx,sd)
		dxy = max(dxy,sd)
		dyx = max(dyx,sd)
		dyy = max(dyy,sd)
	case default
		write(*,*) 'ERROR: unsupported error floor computing method!'
		write(*,*) '       STOP!'
		stop
	end select ! method
!
!  end of the subroutine.
!
	end subroutine CompZerrFlr

!----------------------------------------------------------------------!
!--------------------------------------------------------------! CHKKW !
!----------------------------------------------------------------------!

	subroutine chkkw(kw)
!
!  Routine to check the validity of keyword.
!
	implicit none
	character(*)   :: kw
!
!  Check if the keyword is valid.
!
	select case (trim(adjustl(kw)))
	case ('>FREQ', '>ZROT', '>ZXXR', '>ZXXI', '>ZXX.VAR',              &
		'>ZXYR', '>ZXYI', '>ZXY.VAR', '>ZYXR', '>ZYXI',                &
		'>ZYX.VAR', '>ZYYR', '>ZYYI', '>ZYY.VAR',                      &
		'>TZXR','>TZXI','>TZX.VAR','>TZYR','>TZYI','>TZY.VAR',         &
		'>TXR.EXP','>TXI.EXP','>TYR.EXP','>TYI.EXP',                   &
		'>RHOROT', '>RHOXX', '>RHOXX.ERR', '>RHOXY',                   &
		'>RHOXY.ERR', '>RHOYX', '>RHOYX.ERR', '>RHOYY',                &
		'>RHOYY.ERR', '>PHSXX', '>PHSXX.ERR', '>PHSXY',                &
		'>PHSXY.ERR', '>PHSYX', '>PHSYX.ERR', '>PHSYY',                &
		'>PHSYY.ERR', '>RES1DXY', '>DEP1DXY', '>RES1DYX', '>DEP1DYX',  &
		'>PTXX','>PTXY','>PTYX','>PTYY',&
		'>PTALPHA','>PTBETA','>PTMAX','>PTMIN')
		return

	case default
		write(0,*) 'ERROR: Unsupported keyword: ',kw
		write(0,*) '       Please check the PLT file!'
		write(0,*) '       STOP!'
		stop

	end select
!
!  end of the subroutine.
!
	end subroutine chkkw

!----------------------------------------------------------------------!
!------------------------------------------------------------! GETSVAL !
!----------------------------------------------------------------------!

	subroutine getsval(fid,kw,sVal)
!
!  Subroutine to get the Value string of one keyword from the PLT file.
!
	implicit none
	integer(4)     :: fid
	character(*)   :: kw
	character(*)   :: sVal
	character(180) :: sLine

	read(fid,'(a)') sLine
	do while (index(sLine,kw) .EQ. 0)
		read(fid,'(a)',end=90) sLine
	enddo ! kw
	sVal = sLine(index(sLine,'=')+1: )
	rewind(fid)

	return
!
!  Error information.
!
90  write(0,*) 'ERROR: have not found the keyword: ',kw
	write(0,*) '       STOP!'
	stop
!
!  End of the subroutine.
!
	end subroutine getsval

!----------------------------------------------------------------------!
!----------------------------------------------------------! GETPLTVAL !
!----------------------------------------------------------------------!

	function GetPLTVAL(PLTFILE,kw)
!
!  Function to get the character value of the KeyWord.
!
	implicit none
	character(*)   :: PLTFILE
	character(*)   :: kw
	integer(4)     :: err
	integer(4)     :: k,j
	character(60)  :: GetPLTVAL
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='old',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err

!
!  Get the keyword's string value.
!
	call getsval(101,trim(adjustl(kw)),GetPLTVAL)
!
!  Remove the double quotes if exised.
!
	k = index(GetPLTVAL,'"')
	j = index(GetPLTVAL,'"',BACK=.TRUE.)
	GetPLTVAL = trim(adjustl(GetPLTVAL(k+1:j-1)))
!
!  Close the PLT file.
!
	close(101)
!
!  End of the function.
!
	end function GetPLTVAL

!----------------------------------------------------------------------!
!----------------------------------------------------------! FINDPLTKW !
!----------------------------------------------------------------------!

	function FindPLTKW(PLTFILE,kw)
!
!  Function to find if the KeyWord is in the PLT file or not.
!
	implicit none
	character(*)   :: PLTFILE
	character(*)   :: kw
	character(180) :: sLine
	integer(4)     :: err
	logical        :: FindPLTKW
!
!  Check the validity of the keyword.
!
	!call chkkw(kw)
!
!  Open the PLT file.
!
	open(101,FILE=PLTFILE,status='old',iostat=err)
	if (err .NE. 0) then
		write(0,*) 'ERROR: Can not open PLTFILE: ',PLTFILE
		write(0,*) '       Please check!'
		write(0,*) '       STOP!'
		stop
	endif ! err
!
!  Loop over each line to find the keyword.
!
	read(101,'(a)') sLine
	do while (index(sLine,kw//' ') .EQ. 0 .OR. index(sLine,kw//'.') .NE. 0)
		read(101,'(a)',end=91) sLine
	enddo ! kw
!
!  If we have found the keyword, then close the PLT file and return.
!
	FindPLTKW = .TRUE.
	close(101)
	return
!
!  If we have not found the keyword.
!
91  FindPLTKW = .FALSE.
	close(101)
	return
!
!  End of the function.
!
	end function FindPLTKW
!
!  End of the module.
!
	end module PLTIO

