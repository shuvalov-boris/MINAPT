$NOTSTRICT
	INTERFACE TO INTEGER*2 FUNCTION system[C] (str[REFERENCE])
	  character*1 str
	end

	subroutine kmon(str)
	  character*1 str(1)
	  integer*2 system,i
	  external   system
	  
	  i=system(str)
	  
	  return
	end