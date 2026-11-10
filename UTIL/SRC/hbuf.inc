	
	integer hbuf_length, HBUF_MAX_LENGTH
	parameter(HBUF_MAX_LENGTH = 100000)  ! SAM-ASP writes thousands of profiles
	
	real hbuf(1000*1000)	
	character *20 namelist(HBUF_MAX_LENGTH)  ! as in SAM hbuffer.f90	
	character *80 deflist(HBUF_MAX_LENGTH)	
	character *10 unitlist(HBUF_MAX_LENGTH)	
	integer status(HBUF_MAX_LENGTH)
	integer average_type(HBUF_MAX_LENGTH)

	common /hbuffer/ hbuf,namelist,deflist,unitlist,
     &			 status,average_type,hbuf_length
