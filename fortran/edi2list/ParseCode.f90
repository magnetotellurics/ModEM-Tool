!==============================================================================!
!===================================================================! ParseCode
!==============================================================================!
	subroutine ParseCode( nLen, sLine, sCode, sValue, bComment )

! David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Nov 2006 - parse a line read from a file into a code & value.
! Force the code to be all lowercase with no ending colon.  Terminate
! the line at a '#' or '!' sign (these allow for user comments!)
!
!  Last modified : 2014/01/12 14:08:54.                                !
!  UPDATE LOG(by Bo Yang):
!  Jan 10, 2014: check if a whole line is empty, if so, return with
!                bComment=.True.
!
		implicit none

		! Args
		integer, intent(in)         :: nLen
		character(nLen)             :: sLine
		character(nLen), intent(out) :: sCode, sValue
		logical, intent(out)        :: bComment

		! Local vars
		integer :: iFrom, iTo

		! Init returns
		bComment = .false.
		sCode = ' '
		sValue = ' '

		! Convert all tab characters to spaces
!		forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '

		! Skip any beginning blanks
		do iFrom = 1,nLen
			if (sLine(iFrom:iFrom) .ne. ' ') exit
		enddo

		! If this whole line is empty, then return.(Added by Bo Yang)
		if( iFrom == nLen+1 ) then
			bComment = .true.
			return
		endif

		! If the first char is a comment char, then the whole line is a comment.
		if( sLine(iFrom:iFrom) == '#' .or. sLine(iFrom:iFrom) == '!' ) then
			bComment = .true.
			return
		endif

		! Pull off the code value.  Cvt to lowercase as we go.
		iTo = index(sLine,':') - 1
		if (iTo < iFrom) then
		    !write(*,*) 'Parsing Error: missing colon in line below:'
			!write(*,*) sLine
			bComment = .true. ! KWK Feb 2008 treat blank line like a comment 
			return
		endif
		sCode = sLine(iFrom:iTo)
		call Lower(sCode)

		! Skip spaces after the colon
		do iFrom = iTo+2,nLen
			if (sLine(iFrom:iFrom) .ne. ' ') exit
		enddo

		! Get the rest, up to any comment
		sValue = sLine(iFrom:)
		iTo = len_trim(sValue)

		!   iFrom = min(index(sValue,'%'), index(sValue,'!'))
		!    if (iFrom > 0 .and. iFrom < iTo) then
		!        sValue(iFrom:iTo) = ' '
		!    endif
		! KWK April 18,2008. Fixed the above three lines using this:
		iFrom = index(sValue,'#')
		if (iFrom > 0 .and. iFrom < iTo) then
			sValue(iFrom:iTo) = ' '
		endif
		iFrom = index(sValue,'!')
		if (iFrom > 0 .and. iFrom < iTo) then
			sValue(iFrom:iTo) = ' '
		endif
		!call Lower(sValue)     ! No: Some values are filenames which are case-sensitive on UNIX!

    end subroutine ParseCode


!==============================================================================!
!=======================================================================! LOWER
!==============================================================================!
    subroutine Lower( s )

! David Myer IGPP/SIO La Jolla CA 92093-0225
! DGM Nov 2006 - convert string to lower case
    implicit none
    character(*), intent(out)  :: s
    integer i

    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
      endif
    enddo

    end subroutine Lower
