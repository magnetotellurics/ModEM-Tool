!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                     program EDI2LIST main frame                      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

	program edi2list

!----------------------------------------------------------------------!
!                                                                      !
!  Program to generate a list file for ModEM from a list of edi files. !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com.                                 !
!  Started       : 2014/01/03.                                         !
!  Version       : 1.0                                                 !
!  Released date : 2014/01/03.                                         !
!  Last modified : 2014/03/24 09:59:19.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2014-01-03 Created.                                                 !
!  2014-02-21 I) Added the extra-pos reading and frequency selection   !
!             support.                                                 !
!             II) Updated for compatibility of using PLTIO v1.1.       !
!                                                                      !
!----------------------------------------------------------------------!

	use PLTIO, only : PreReadPLT,ReadPLT,GetPLTVAL,CompZerrFlr,FindPLTKW
	use LATLON,only : ll2xy

	implicit none

	character(99) :: RUNFILE,EDILIST,OUTFILE,EDIFILE,FREQLIST,POSFILE
	character(40) :: temp,DTYPE
	character(60) :: DATAID
	integer(4)    :: err
	real(8)       :: lat0,lon0,lat,lon,elev
	real(8)       :: errflr,dxx,dxy,dyx,dyy
	integer(4)    :: errmethod
	complex(8)    :: ZXX,ZXY,ZYX,ZYY
	integer(4)    :: nsite,nfreq,nf,kf,kk
	real(8)       :: x,y,z
	real(8),allocatable,dimension(:) :: freq,efreq,fsel
	real(8),allocatable,dimension(:) :: zxxr,zxxi,zxyr,zxyi,zyxr,zyxi,zyyr,zyyi
	real(8),allocatable,dimension(:) :: dzxx,dzxy,dzyx,dzyy
	real(8),   parameter :: EPS= 1.d-6
	complex(8),parameter :: II = (0.d0,1.d0)

	character(56)   :: vers
	character(180)  :: sLine,sCode,sValue
	logical         :: bCommnet

	integer(4)    :: n,i
	real(8)       :: ZERO = 0

	logical       :: existed

	integer(4)     :: Linenum

! -Variables- !--------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Initialize all options.
!
	EDILIST   = ''
	FREQLIST  = ''
	OUTFILE   = ''
	DTYPE     = ''
	lat0      = ZERO/ZERO !NaN
	lon0      = ZERO/ZERO !NaN
	errflr    = ZERO/ZERO !NaN
	errmethod = 0
!
!  Deal with the command line options.
!
	n = iargc()
	select case (n)
	case (1)
		call getarg(1,RUNFILE)
!----------------------------------------------------------------------!
!  Read runfile                                                        !
!----------------------------------------------------------------------!
!
!  Open runfile.
!
		open(unit=100,file=RUNFILE,status='old',iostat=err)
		if (err .NE. 0) then
			write(*,*) 'Error opening RUNFILE: '//trim(RUNFILE)//' !'
			write(*,*) 'Please check the RUNFILE and try again.'
			write(*,*) 'STOP!'
			stop
		endif ! open runfile err.
!
!  Read in all parameters line by line.
!
		do while (.TRUE.)
!
!  Get the next code/value pair.
!
			read(100,'(A)',iostat=err) sLine
!
!  Check if is the end of the file.
!
			if (err .NE. 0) exit
!
!  Decode the current line.
!
			call ParseCode(len(sLine),sLine,sCode,sValue,bCommnet)
!
!  If is commnet, then skip to next line.
!
			if (bCommnet) cycle
!
!  Seeing what do we have.
!
			select case (trim(sCode))
			
			case ('version')
				vers = sValue(1:len(vers))
				write(*,*) 'The RUNFILE version is ',vers
				
			case ('edilist')
				EDILIST = trim(sValue(1:len(EDILIST)))
				
			case ('freqlist')
				FREQLIST = trim(sValue(1:len(FREQLIST)))
				
			case ('outfile')
				OUTFILE = trim(sValue(1:len(OUTFILE)))
				
			case ('datatype')
				DTYPE = trim(sValue(1:len(DTYPE)))
				
			case ('lat0lon0','origin','geoorigin')
				read(sValue,*) lat0,lon0
				
			case ('errorfloor')
				read(sValue,*) errflr
				
			case ('seterrorfloormethod')
				read(sValue,*) errmethod
				
			case default
				write(*,*) 'Error reading RUNFILE!'
				write(*,*) 'Unknown or unsupported code: ',sCode
				write(*,*) 'Please check the RUNFILE and try again.'
				stop
				
			end select ! sCode
!
!  move to next line.
!
		enddo ! while.
!
!  Close runfile.
!
		close(100)
!
!  Check all required parameters are read.
!
		if (EDILIST .EQ. '') then
			write(*,*) 'ERROR: can not find EDILIST file name! STOP!'
			stop
		endif ! EDIFILE.
		if (FREQLIST .EQ. '') then
			write(*,*) 'ERROR: can not find FREQLIST file name! STOP!'
			stop
		endif ! FREQFILE.
		if (OUTFILE .EQ. '') then
			write(*,*) 'WARNING: can not find OUTFILE file name!'
			write(*,*) '         use the default OUTFILE=edi2list.out'
			OUTFILE = 'edi2list.out'
		endif ! OUTFILE.
		if (DTYPE .EQ. '') then
			write(*,*) 'WARNING: can not find data type info!'
			write(*,*) '         use the default DTYPE=Full_Impedance'
			DTYPE = 'Full_Impedance'
		endif ! OUTFILE.
		if (isnan(lat0) .OR. isnan(lon0)) then
			write(*,*) 'ERROR: can not find origin info! STOP!'
			stop
		endif ! lat0/lon0
		if (isnan(errflr)) then
			write(*,*) 'WARNING: can not find ErrorFloor.'
			write(*,*) '         use the real error in each EDI file.'
			errflr = 0.d0
		endif ! errflr
		if (errmethod .EQ. 0) then
			write(*,*) 'WARNING: can not find SetErrorFloorMethod.'
			write(*,*) '         use the default ErrMethod=1'
			errmethod = 1
		endif ! errflr
!----------------------------------------------------------------------!
!  End of reading runfile.                                             !
!----------------------------------------------------------------------!
	case (8)
		i=0;
		i=i+1; call getarg(i,EDILIST)
		i=i+1; call getarg(i,FREQLIST)
		i=i+1; call getarg(i,OUTFILE)
		i=i+1; call getarg(i,DTYPE)
		i=i+1; call getarg(i,temp)
		read(temp,*) lat0
		i=i+1; call getarg(i,temp)
		read(temp,*) lon0
		i=i+1; call getarg(i,temp)
		read(temp,*) errflr
		i=i+1; call getarg(i,temp)
		read(temp,*) errmethod
	case default
		write(*,*) 'EDI2LIST v1.0'
		write(*,*) '- A tool to generate a list data file for ModEM from'
		write(*,*) '  a list of a PLT/EDI file.'
		write(*,*) '- Written by Yang Bo, IGG, CUG, 2014/01/02.'
		write(*,*) '- Usage: edi2list runfile'
		write(*,*) '         edi2list edilist.txt freq.list out.list datatype lat0 lon0 errflr errmethod'
		write(*,*) '         lon0 lat0: The lon/lat for the origin point'
		write(*,*) '                    of the survey area.'
		write(*,*) '         datatype: Off_Diagonal_Rho_Phase.'
		write(*,*) '                   Off_Diagonal_Impedance.'
		write(*,*) '                   Full_Impedance.'
		write(*,*) '                   Full_Vertical_Components.'
		stop
	end select ! n
!
!  get the frequency list.
!
	nfreq = Linenum(FREQLIST)
	allocate(freq(nfreq))
	open(300,file=FREQLIST,status='old',iostat=err)
	if (err .NE. 0) then
		write(*,*) 'ERROR: Can not open FREQ LIST file: ',FREQLIST
		write(*,*) '       Please check!'
		write(*,*) '       STOP!'
		stop
	endif ! err
	read(300,*) (freq(kf),kf=1,nfreq)
	close(300)
!
!  open the out file.
!
	open(200,file=OUTFILE,status='unknown',iostat=err)
	if (err .NE. 0) then
		write(*,*) 'ERROR: Can not open OUT file: ',OUTFILE
		write(*,*) '       Please check!'
		write(*,*) '       STOP!'
		stop
	endif ! err
!
!  Open the EDI list file.
!
	nsite = Linenum(EDILIST)
	open(100,file=EDILIST,status='old',iostat=err)
	if (err .NE. 0) then
		write(*,*) 'ERROR: Can not open EDI LIST file: ',EDILIST
		write(*,*) '       Please check!'
		write(*,*) '       STOP!'
		stop
	endif ! err
!
!  write out the header infomation of the output list file.
!
	write(200,'(A)') '#'
	write(200,'(A)') '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'
	write(200,'(A)') '> '//trim(adjustl(DTYPE))
	write(200,'(A)') '> exp(i\omega t)'
	write(200,'(A)') '> [mV/km]/[nT]'
	write(200,'(A)') '> 0.0'
	write(200,'(a2,f8.3,f8.3)') '> ',lat0,lon0
	write(200,'(a2,i4,i5)') '> ',nfreq,nsite
!
!  read the list file and loop over all sites.
!
	do while (.TRUE.)
!
!  read EDIFILE name one by one.
!
		!read(100,'(A)',iostat=err) EDIFILE
		read(100,*,iostat=err) EDIFILE
		if (err .NE. 0) exit
!
!  read the head information from a position file or the original EDI file.
!
		POSFILE = adjustl(trim(EDIFILE))//'.pos'
		inquire(file=POSFILE,exist=existed)
		if (existed) then
			call PreReadPLT(POSFILE,DATAID,nf,lat,lon,elev)
		else
			call PreReadPLT(EDIFILE,DATAID,nf,lat,lon,elev)
		endif ! position file is existed or not.
!
!  read frequency and it's selection info of current site.
!
		allocate(efreq(nf),fsel(nf))
		call ReadPLT(EDIFILE,nf,'>FREQ',efreq)
		if (FindPLTKW(EDIFILE,'>FSEL')) then
			call ReadPLT(EDIFILE,nf,'>FSEL',fsel)
		else
			fsel = 1;
		endif ! if.
!
!  read the impedance tensor.
!
		allocate(zxxr(nf),zxxi(nf),zyyr(nf),zyyi(nf))
		allocate(zxyr(nf),zxyi(nf),zyxr(nf),zyxi(nf))
		allocate(dzxx(nf),dzxy(nf),dzyx(nf),dzyy(nf))
		call ReadPLT(EDIFILE,nf,'>ZXXR',zxxr)
		call ReadPLT(EDIFILE,nf,'>ZXXI',zxxi)
		call ReadPLT(EDIFILE,nf,'>ZXYR',zxyr)
		call ReadPLT(EDIFILE,nf,'>ZXYI',zxyi)
		call ReadPLT(EDIFILE,nf,'>ZYXR',zyxr)
		call ReadPLT(EDIFILE,nf,'>ZYXI',zyxi)
		call ReadPLT(EDIFILE,nf,'>ZYYR',zyyr)
		call ReadPLT(EDIFILE,nf,'>ZYYI',zyyi)
		call ReadPLT(EDIFILE,nf,'>ZXX.VAR',dzxx)
		call ReadPLT(EDIFILE,nf,'>ZXY.VAR',dzxy)
		call ReadPLT(EDIFILE,nf,'>ZYX.VAR',dzyx)
		call ReadPLT(EDIFILE,nf,'>ZYY.VAR',dzyy)
!
!  write out type of data requested.
!
		select case (DTYPE)
		
		case ('Off_Diagonal_Impedance')
			! output
			do kf = 1,nfreq
				LOOP_FREQ1: do kk = 1,nf
					! see if current frequncy is removed or not
					if (dabs(fsel(kk)) .LT. EPS) then
						cycle LOOP_FREQ1
					endif ! if.
					if (dabs(efreq(kk)-freq(kf)) .LT. EPS) then
						! convert latlon to xy.
						call ll2xy(lat0,lon0,lat,lon,x,y)
						x = x*1000.d0
						y = y*1000.d0
						z = -elev
						! get site name.
						!DATAID = GetPLTVAL(EDIFILE,'DATAID=')
						! compute error floor.
						ZXX = zxxr(kk) + II*zxxi(kk)
						ZXY = zxyr(kk) + II*zxyi(kk)
						ZYX = zyxr(kk) + II*zyxi(kk)
						ZYY = zyyr(kk) + II*zyyi(kk)
						dxx=dzxx(kk);dxy=dzxy(kk);dyx=dzyx(kk);dyy=dzyy(kk)
						call CompZerrFlr(errflr,errmethod,ZXX,ZXY,ZYX, &
							             ZYY,dxx,dxy,dyx,dyy)
						! write out ZXY.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZXY',zxyr(kk),zxyi(kk),dxy
						! write out ZYX.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZYX',zyxr(kk),zyxi(kk),dyx
						exit LOOP_FREQ1
					endif ! find freq.
				enddo LOOP_FREQ1 ! kk
			enddo ! kf
			
		case ('Full_Impedance')
			! output
			do kf = 1,nfreq
				LOOP_FREQ2: do kk = 1,nf
					! see if current frequncy is removed or not
					if (dabs(fsel(kk)) .LT. EPS) then
						cycle LOOP_FREQ2
					endif ! if.
					if (dabs(efreq(kk)-freq(kf)) .LT. EPS) then
						! convert latlon to xy.
						call ll2xy(lat0,lon0,lat,lon,x,y)
						x = x*1000.d0
						y = y*1000.d0
						z = -elev
						! get site name.
						!DATAID = GetPLTVAL(EDIFILE,'DATAID=')
						! compute error floor.
						ZXX = zxxr(kk) + II*zxxi(kk)
						ZXY = zxyr(kk) + II*zxyi(kk)
						ZYX = zyxr(kk) + II*zyxi(kk)
						ZYY = zyyr(kk) + II*zyyi(kk)
						dxx=dzxx(kk);dxy=dzxy(kk);dyx=dzyx(kk);dyy=dzyy(kk)
						call CompZerrFlr(errflr,errmethod,ZXX,ZXY,ZYX, &
							             ZYY,dxx,dxy,dyx,dyy)
						! write out ZXY.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZXX',zxxr(kk),zxxi(kk),dxx
						! write out ZXY.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZXY',zxyr(kk),zxyi(kk),dxy
						! write out ZYX.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZYX',zyxr(kk),zyxi(kk),dyx
						! write out ZYX.
						write(200,90) 1.d0/freq(kf),trim(DATAID),lat,lon,&
							          x,y,z,'ZYY',zyyr(kk),zyyi(kk),dyy
						exit LOOP_FREQ2
					endif ! find freq.
				enddo LOOP_FREQ2 ! kk
			enddo ! kf	
			
		case ('Full_Vertical_Components')
			
		case ('Off_Diagonal_Rho_Phase')
			
		case default
		end select ! DTYPE
!
!  release momeries for current site.
!
		deallocate(efreq,fsel)
		deallocate(zxxr,zxxi,zxyr,zxyi,zyxr,zyxi,zyyr,zyyi)
		deallocate(dzxx,dzxy,dzyx,dzyy)
!
!  move to next site.
!
	enddo ! while.
!
!  close files.
!
	close(100)
	close(200)
!
!  release all memeroies.
!
	deallocate(freq)
!
!  data format.
!
90  format(es14.6,a10,f8.3,f8.3,f12.3,f12.3,f12.3,1x,a3,es14.6,es14.6,es14.6)
!
!  End of the program.
!
	end program edi2list

!----------------------------------------------------------------------!
!--------------------------------------------------------------! UPPER
!----------------------------------------------------------------------!

	function Linenum(FILENAME)

! Bo Yang, IGG/CUG, Wuhan, China.
! Dec 2012 - get the number of lines in text file.
	implicit none
	character(*)   :: FILENAME
	integer(4)     :: Linenum
	integer(4)     :: k
	integer(4)     :: err
	character(280) :: str
	open(100,file=FILENAME,status='old',iostat=err)
	if (err .NE. 0) then
		write(*,*) 'ERROR: Can not open: ',FILENAME
		write(*,*) '       Please check!'
		write(*,*) '       STOP!'
		stop
	endif ! err
	Linenum = 0
	do while(.TRUE.)
		read(100,'(A)',iostat=err) str
		if (err .NE. 0) exit
		Linenum = Linenum + 1
	enddo ! k
	close(100)
	end function Linenum

