!----------------------------------------------------------------------!
!-------------------------------------------------------------! LATLON
!----------------------------------------------------------------------!

	module LATLON

!----------------------------------------------------------------------!
!                                                                      !
!  Module to deal with latitude and longitude.                         !
!                                                                      !
!  See also XXX.                                                       !
!                                                                      !
!  Written by    : Yang Bo, Institute of Geophysics and Geomatics,     !
!                  China University of Geosciences, Wuhan, China.      !
!  E-mail        : yangbo.cug@163.com                                  !
!  Started       : 2013/12/24                                          !
!  Last modified : 2014/01/03 22:58:32.                                !
!  Copyright 2011-2016 Yang Bo, IGG, CUG.                              !
!  Comments, bug reports and questions, please send me an email.       !
!                                                                      !
!                             Modifications                            !
!                             =============                            !
!                                                                      !
!  2013-12-24 Created.                                                 !
!                                                                      !
!----------------------------------------------------------------------!

	implicit none

	real(8),parameter :: PI = 3.14159265358979d0

	public d2dms,dms2d

!
!  The contained subroutines.
!
	contains

!----------------------------------------------------------------------!
!--------------------------------------------------------------! DMS2D !
!----------------------------------------------------------------------!

	function dms2d(dms)
!
!  Function to convert longitude/latitude from Degrees,Minutes,Seconds
!  format to Decimal degrees format.
!
	implicit none
	character(*) :: dms
	real(8)      :: dms2d
	real(8)      :: dd,mm,ss
	integer(4)   :: idm,ims
	character(60):: sval

	idm = scan(dms,':')
	ims = scan(dms,':',back=.TRUE.)
	read(dms(1:idm-1),    *) dd
	read(dms(idm+1:ims-1),*) mm
	read(dms(ims+1: ),    *) ss
	dms2d = dabs(dd) + mm/60.d0 + ss/3600.d0
	if (dd .LT. 0) dms2d = - dms2d
!
!  End of the function.
!
	end function dms2d

!----------------------------------------------------------------------!
!--------------------------------------------------------------! D2DMS !
!----------------------------------------------------------------------!

	function d2dms(deg)
!
!  Function to convert longitude/latitude from Decimal degrees format
!  to Degrees,Minutes,Seconds format.
!
	implicit none
	character(15) :: d2dms
	real(8)       :: deg
	integer(4)    :: dd,mm
	real(8)       :: ss
	character(3)  :: sdd
	character(2)  :: smm
	character(8)  :: sss

	dd = int(dabs(deg))
	mm = int((dabs(deg)-dd)*60.d0)
	ss = (dabs(deg)-dd-dble(mm)/60.d0)*3600.d0
	sdd = ''
	smm = ''
	sss = ''
	write(sdd,'(i3.3)') dd
	write(smm,'(i2.2)') mm
	write(sss,'(f8.5)') ss
	d2dms = sdd//':'//smm//':'//adjustl(sss)
	if (deg .LT. 0.d0) then
		d2dms = '-'//d2dms
	endif ! deg < 0.
!
!  End of the function.
!
	end function d2dms

!----------------------------------------------------------------------!
!--------------------------------------------------------------! XY2LL !
!----------------------------------------------------------------------!

	subroutine xy2ll(lat0,lon0,x,y,lat,lon)
!
!  Subroutine for converting lon/lat to x/y by using Anna's algorithm.
!  (lat0,lon0) - origin.
!  x+ -> North, y+ -> East.
!  delx = KmPerDeg * (lat - lat0)
!  dely = KmPerDeg * (lon - lon0) * cos((lat+lat0)/2)
!
	implicit none
	real(8) :: lon0,lat0
	real(8) :: lon, lat
	real(8) :: x,   y
	real(8) :: KmPerDeg
	real(8) :: az,baz

	call delaz(lat0,lon0,lat0+1.d0,lon0,KmPerDeg)
	lat = x / KmPerDeg + lat0
	call delaz(lat0,lon0,lat0,lon0+1.d0,KmPerDeg)
	lon = y / KmPerDeg + lon0
!
!  End of the subroutine.
!
	end subroutine xy2ll

!----------------------------------------------------------------------!
!--------------------------------------------------------------! LL2XY !
!----------------------------------------------------------------------!

	subroutine ll2xy(lat0,lon0,lat,lon,x,y)
!
!  Subroutine for converting lon/lat to x/y by calling azdist.
!  (lat0,lon0) - origin.
!  x+ -> North, y+ -> East.
!
	implicit none
	real(8) :: lon0,lat0
	real(8) :: lon, lat
	real(8) :: KmPerDeg
	real(8) :: x,   y

	call delaz(lat0,lon0,lat,lon0,x)
	call delaz(lat0,lon0,lat0,lon,y)
	if (lat < lat0) x = -x
	if (lon < lon0) y = -y
!
!  End of the subroutine.
!
	end subroutine ll2xy

!----------------------------------------------------------------------!
        subroutine delaz(elat,elon,slat,slon,delkm)
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
        implicit none
        real(8) :: elat,elon,slat,slon,delkm
        real(8),parameter :: rad=0.017453292,e2=0.9932315,re=6371.003
        real(8) :: selat,celat,selon,celon,ea,eb,ec
        real(8) :: sslat,cslat,sslon,cslon,sa,sb,sc
        real(8) :: cdel,fac,del

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
        elseif (cdel .le. 0.9396 .and. cdel .ge. -0.9396) then
            del = acos(cdel)
!-----
!     if DEL = (160,180]
!-----
        else
            fac = (ea+sa)**2 + (eb+sb)**2 + (ec+sc)**2
            fac = sqrt(fac)/2.0
            del = 2.0*acos(fac)
        endif

        delkm = del * re
        del = del / rad

        return
        end subroutine delaz

!----------------------------------------------------------------------!
        subroutine dircos(lat,lon,slat,clat,slon,clon,aa,bb,cc)
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
        implicit none
        real(8) :: lat,lon,slat,slon,clat,clon, aa, bb, cc
        real(8) :: c,s,e4,fac
        real(8) :: rad,e2

        rad=0.017453292
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
        end subroutine dircos

!
!  End of the module.
!
	end module LATLON

