#include "fintrf.h"
!----------------------------------------------------------------------!
!--------------------------------------------------------! MEXFUNCTION
!----------------------------------------------------------------------!

	subroutine mexFunction(nlhs, plhs, nrhs, prhs)

!----------------------------------------------------------------------!
!                                                                      !
!  Gateway routine for calling delaz.f from matlab.                    !
!                                                                      !
!  See also XXX.                                                       !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com                                  !
!  Started       : 2014/02/14                                          !
!  Last modified : 2014/02/18 10:16:19.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2014-02-14 Created.                                                 !
!                                                                      !
!----------------------------------------------------------------------!

	implicit none
!
!  mexFunction arguments.
!
	mwPointer plhs(*), prhs(*)
	integer   nlhs, nrhs
!
!  API function declarations.
!
	mwPointer mxGetPr
	mwPointer mxCreateDoubleMatrix
	integer   mxIsNumeric
	mwSize    mxGetM, mxGetN
!
!  Pointers to input/output mxArrays:
!
	mwPointer pelat,pelon,pslat,pslon,pdel
!
!  Array information:
!
	mwSize    mrows, ncols, size
!
!  Arguments for computational routine:
!
	real*8  elat,elon,slat,slon,del
	real*4  elat4,elon4,slat4,slon4,del4
	real*4  az,baz,delkm,saz,caz

	character(80) msg

! -Variables- !--------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!
!  Check for proper number of arguments.
!
	if(nrhs .ne. 4) then
		msg = 'Four input required!'
		call mexErrMsgIdAndTxt('MATLAB:dellatlon:nInput',trim(msg))
	elseif(nlhs .gt. 1) then
		msg = 'Too many output arguments.'
		call mexErrMsgIdAndTxt ('MATLAB:dellatlon:nOutput',trim(msg))
	endif
!
!  Create Fortran array from the input argument.
!
	pelat = mxGetPr(prhs(1))
	pelon = mxGetPr(prhs(2))
	pslat = mxGetPr(prhs(3))
	pslon = mxGetPr(prhs(4))
	call mxCopyPtrToReal8(pelat,elat,1)
	call mxCopyPtrToReal8(pelon,elon,1)
	call mxCopyPtrToReal8(pslat,slat,1)
	call mxCopyPtrToReal8(pslon,slon,1)
!
!  Create matrix for the return argument.
!
	plhs(1) = mxCreateDoubleMatrix(1,1,0)
	pdel    = mxGetPr(plhs(1))
!
!  Call the computational subroutine.
!
	elat4 = real(elat,4)
	elon4 = real(elon,4)
	slat4 = real(slat,4)
	slon4 = real(slon,4)
	call delaz(elat4,elon4,slat4,slon4,del4,az,baz,delkm,saz,caz)
!
!  debugging output.
!
	!write(msg,*) delkm
	!call mexErrMsgTxt(msg)
!
!  Load the delkm into pdel, which is the output to MATLAB.
!
	call mxCopyReal8ToPtr(dble(delkm),pdel,1)
!
!  End of the subroutine.
!
	return
	end subroutine mexFunction

!----------------------------------------------------------------------!
!  --------------------subrouitnes of delaz--------------------------  !
!----------------------------------------------------------------------!

        subroutine delaz(elat,elon,slat,slon,del,az,baz,delkm,saz,caz)
!-----
!     elat    R*4 Epicenter latitude (degrees)
!     elon    R*4 Epicenter longitude (degrees)
!     slat    R*4 Station latitude (degrees)
!     slon    R*4 Station longitude (degrees)
!     del R*4 Epicentral Distance (degrees)
!     az  R*4 Epicenter to Station Azimuth
!     baz R*4 Station to Epicenter Backazimuth
!     delkm   R*4 Epicentral Distance (km)
!     saz R*4 Sine of Azimuth
!     caz R*4 Cosine of Azimuth
!-----
        real*4 rad, e2, re
        data rad/0.017453292/,e2/0.9932315/, re/6371.003/
        call dircos(elat,elon,selat,celat,selon,celon,ea,eb,ec)
        call dircos(slat,slon,sslat,cslat,sslon,cslon,sa,sb,sc)

!-----
!     compute distance
!     Choose correct formula for short and large distances
!-----
        cdel = (ea*sa + eb*sb + ec*sc)
!-----
!     if DEL = [0,20)
!-----
        if(cdel .gt. 0.9396)then
            fac = (ea-sa)**2 + (eb-sb)**2 + (ec-sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*asin(fac)
!-----
!     if DEL = [20,160]
!-----
        else if(cdel .le. 0.9396 .and. cdel .ge. -0.9396)then
            del = acos(cdel)
!-----
!     if DEL = (160,180]
!-----
        else
            fac = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*acos(fac)
        endif
        
!-----
!     check for station or epicenter at pole
!-----
        if(elat.eq.90.0 .or. slat.eq.-90.0)then
            az = 180.0
            baz = 0.0
            saz =  0.0
            caz = -1.0
        else if(elat.eq.-90.0 .or. slat.eq.90.0)then
            az = 0.0
            baz = 180.0
            saz = 0.0
            caz = 1.0
        else
            saz = celat*(cslat * sin(rad*(slon - elon)))
            caz = (sslat - cdel*selat)
            fac = sqrt(saz**2 + caz**2)
            if(fac.gt.0.0)then
                saz = saz / fac
                caz = caz / fac
                az = atan2(saz,caz)
            
                sbz = - cslat*(celat * sin(rad*(slon - elon)))
                cbz = (selat - cdel*sslat)
                baz = atan2(sbz,cbz)
            else
                az = 0.0
                caz = 1.0
                saz = 0.0
                baz = 180.0
            endif
            az = az / rad
            baz = baz / rad
        endif
        delkm = del * re
        del = del / rad

!-----
!     put az and baz in the range [0,360)
!-----
        if(az .lt. 0.0)az = az + 360.0
        if(baz .lt. 0.0)baz = baz + 360.0
        return
        end
        
        subroutine dircos(lat,lon,slat,clat,slon,clon,aa,bb,cc)
        real*4 lat,lon,slat,slon,clat,clon, aa, bb, cc
        real*4 rad, e2
        data rad/0.017453292/,e2/0.993305615/
!-----
!     convert geographic latitude to geocentric
!     Use flattening of Chovitz (1981) f= 1/298.257 adopted by IUGG 1980
!     
!     The relation between geocentric and geographic latitude is
!     tan phi c = ( 1 - f)^2 tan phi g
!
!     To avoid problems at poles, define sin phi c and cos phi c
!     so that the relation holds ans also that s^2 + c^2 = 1
!
!     For geographic to geocentric use e2 = (1-f)^2 = 0.993395615
!     For geographic to equidistant use e2 =(1-f)^1.5=0.99776354
!     Brown, R. J. (1984). On the determination of source-receiver
!         distances using a new equidistant latitude,
!         Geophys. J. R. astr. Soc. 76, 445-459.
!     
!-----
        e2 = 0.99776354
        c = cos(rad*lat)
        s = sin(rad*lat)
        e4 = e2**2
        fac = sqrt(e4 + (1.0-e4)*c*c)
        slat = e2 * s /fac
        clat =      c /fac
        slon = sin(rad*lon)
        clon = cos(rad*lon)

        aa = clat * clon
        bb = clat * slon
        cc = slat
        
        return
        end
