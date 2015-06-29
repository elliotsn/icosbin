c******************************************************************************
c  PROGRAM: ptrilink
c
c  PURPOSE: This program is used after triangle numbers (bin numbers) have been added
c			using ptrinum. It can be used following ptrinum directly or following
c			ptrigather, which gathers points with the same triangle numbers.
c			The idea is that before ptrinum, the point stream originally came from 
c			pfootprint, where the output comprises the original Diviner point, followed by
c			n points representing that points effective field of view (EFOV). In ptrinum, the
c			original point is assigned a weight of zero and its component points that
c			represent its field of view should have weights that sum to 1.
c
c			In the group of points following the original point, many fields are
c			unnecessarily repeated because they are invariant between points in the same 
c			EFOV probability distribution. In order to reduce the 
c			amount of data that is output, this program splits the stream into two linked
c			lists. The original points, containing a large number of fields, are output
c			to one file. Whereas the gridded points (n per original point) have only the
c			variant values extracted and are written to a separate file. This vastly
c			reduces the space occupied by output data (by a factor typically between 9 and
c			 10) and makes much easier the building of the foundation database, because 
c			this rearranging can be distributed on the cluster, rather than on the build 
c			nodes.
c
c			This program does not output anything to stdout except reporting.
c
c			Fields that vary for each point in the input stream include:
c
c			talt,tslope,temis,tinc,tinsol,tphase,taz,trinum1,trinum2,trilat,trilon,weight
c
c			These fields will be written to the path in gfile. All fields that are not 
c			passed after the 'fields=' argument are assumed to be part of the original
c			point and are written to the path in ofile.
c
c  REQUIRED KEYWORDS: 
c
c		fields  =char -- Variant fields to store in gridded points file. Remaining fields
c						 in input pipe will be stored in original points file. This 
c						 keyword MUST always be the first argument to ptrilink.
c
c       ofile=char -- Filepath to write the original points to.
c
c		gfile=char -- Filepath to write the gridded points to.
c
c		weight  =char -- The name of the column that contains the weight of each point.
c						 Where weight==0 then this is the original point and all
c						 columns listed after the keyword 'fields' are ignored and 
c						 only the other columns are stored in the original points
c						 file (except weight).
c
c  OPTIONAL KEYWORDS:
c
c       des -- The descriptor file to describe the format and the
c                   fields of each record of the input data
c                   (the default is the previous pipe).
c		
c		nodes -- To prohibit the writing of the descriptor file. Under normal operation
c				 the des file is written as separate files created with names based
c				 on those provided in ofile and gfile.
c
c  Copyright (C) 2015  Elliot Sefton-Nash
c
c  LICENSE
c
c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.

c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License along
c    with this program; if not, write to the Free Software Foundation, Inc.,
c    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
c***********************************************************************************

	    program ptrilink
		implicit none

c  Variable declaration of variables common to all pipes
		include 'pipe1.inc'

	    integer MAXARG
        parameter (MAXARG=2048)
  
        integer*4 i,j,ifound,iunit,ndeslines,ios,ifirst
        character*80 cdesfile,junk
        character*8000 arg,key
		character*132 ofile, gfile
		character*10 odes, gdes

        character*48 gdesstitle(MAXCOL), odesstitle(MAXCOL)
		character*132 gdesltitle(MAXCOL), odesltitle(MAXCOL)
        
        integer*4 odeslines,gdeslines
        
        character*48 param(MAXCOL), weightp
        integer colnum(MAXCOL),gcolnum(MAXCOL),ocolnum(MAXCOL)
        integer narg,done

c		Have to be real so we can write to pipes format.
		integer*8 ocnt,gcnt

        character cjunk

c		We don't yet know how many columns will be in each of the output pipes.
c		gridded points are now stored as 4-byte reals to reduce disk consumption by
c		40-50% with no noticeable loss of precision.
		real*4, allocatable :: grdat(:)
		real*8, allocatable :: ordat(:)

        integer irecord, weighti

         
c 		Check for descriptor file in command line
        call getcmdchar('des',1,'optional',ifound,cdesfile)
        if( ifound .eq. 1 ) then
           iunit=7 ! If descriptor file found, read from unit 7
        else
           iunit=5 ! Else, get the descriptor from standard input
        endif
        call pdesread(iunit, cdesfile,  cdesheader,
     &               cdesstitle,  cdesltitle,
     &               ndeslines)

c 		Check the required keywords are present.
        call getcmdchar('ofile',1,'required',ifound, ofile)
        call getcmdchar('gfile',1,'required',ifound, gfile)
        
        call getcmdchar('weight',1,'required',ifound, weightp)
c		Get column number of this param in the input descriptor
        call findcol(weightp,cdesstitle,ndeslines, weighti)

c 		Find the column indices of the fields in the variant file.

c 		Retrieve all the column names that the user wants extracted
        do i=1,iargc()
           call getarg(i,arg)
           call getparam(arg,key)
           if(key.eq.'fields') then
c			  Get narg, the number of fields after the keyword 'fields'           
              call getnumval(arg,narg)
              call getcharval(arg,param)
              goto 5
           endif
           stop 'ptrilink error: fields= keyword is missing'
        enddo
        
c		Check that there aren't more fields requested than there are in the descriptor file
5       if(narg.gt.ndeslines) then
           write(0,*) 'ERROR in ptrilink: Trying to extract',
     &                narg,' values from a ',ndeslines,'-value dataset'
           call exit()
        endif

c 		Match the requested fields with their corresponding column numbers in descriptor file
c		Loop over all the fields in the input descriptor to check if each one matches one of the
c		ones after the keyword, 'fields'.

        gdeslines = 0	! number of columns in gridded point descriptor
		odeslines = 0	! number of columns in original point descriptor

c		For each column in the input descriptor        
        do i=1,ndeslines
        
c			Loop over each column listed after the 'fields' keyword.
            do j=1,narg
c			  If there's a match between the desired variant fields and this column
              if(param(j).eq.cdesstitle(i)) then

c				Increment the number of variant/gridded fields whose column number has
c				been found.
                gdeslines=gdeslines+1

c               Counter in list of fields is j, position in original cdesfile is i.
                gcolnum(gdeslines) = i
                gdesstitle(gdeslines) = cdesstitle(i)
                gdesltitle(gdeslines) = cdesltitle(i)
              	              	
c				Continue looping through until last field
              	goto 8           
              endif
            enddo
            
c			If a match wasn't found in the above loop, then we can assume that this
c			column is one to be added to the original point descriptor.
			odeslines=odeslines+1
			ocolnum(odeslines) = i ! i is the index of the column in the original descriptor 
			odesstitle(odeslines) = cdesstitle(i)
			odesltitle(odeslines) = cdesltitle(i)
            
c			If at any point there are less keywords left in the descriptor than are 
c			still required as listed after the 'fields' keywords, then report error and stop.
 8			if((narg-gdeslines).gt.(ndeslines-i)) then
				write (0,*)  
     & 'error in ptrilink: all columns after fields',param(j),
     & ' keyword not present in descriptor file'
                 stop
			endif     
        enddo

c		If not all the fields have been allocated, error.
		if((gdeslines + odeslines).ne.ndeslines) then
		write (0,*) 'error in ptrilink: Not all columns allocated.'
        stop
		endif

c		Report on descols
		write (0,*) '**********	FILE CONTENTS **********'
		write (0,*) ofile,' will contain: '
		do i=1,odeslines
			write (0,*) odesstitle(i)
		enddo
		write (0,*) '***********************************'
		write (0,*) gfile,' will contain: '
		do i=1,gdeslines
			write (0,*) gdesstitle(i)
		enddo
		write (0,*) 'opointnum'
    	write (0,*) '***********************************'

c		Number of columns for gridded points is one more than the number of extracted 
c		fields because we add one extra column to store the original point number		
        gdeslines = gdeslines + 1
		
c		Add original point number to gridded/variant output desc.
		gdesstitle(gdeslines)='opointnum'
		gdesltitle(gdeslines)='Original point number'

        cdesheader='Descripter file auto-generated by pipes'

c 		Check for 'nodes' in command line 
        call getcmdchar('nodes',0,'optional',ifound, cjunk)

		odes='fdsdbo.des'
	    gdes='fdsdbg.des'

c		Make sure the files don't exist, because pdeswrite returns an error if the file
c		already exists. Files are deleted on close.
		open (10, file=odes, form='unformatted', status='unknown')
		close (10, status='delete')
		open (11, file=gdes, form='unformatted', status='unknown')
		close (11, status='delete')

c		If nodes not found then write descriptors.
        if(ifound.eq.0) then

c		   For original points           
           call pdeswrite(10, odes,cdesheader,odesstitle,odesltitle,
     &			    odeslines)
c		   pdeswrite (in pipelibfort_g77.f doesn't add the 'end' marker to a 
c		   descriptor, so we append it manually.
   		   open (unit=10,file=odes,status='old',position='append')
		   write(10,*) "'",'end',"'"
           close(unit=10)
           
c		   For gridded/invariant points
           call pdeswrite(11, gdes,cdesheader,gdesstitle,gdesltitle,
     &			    gdeslines)
		   open (unit=11,file=gdes,status='old',position='append')
		   write(11,*) "'",'end',"'"
		   close(unit=11)
		   
        endif

        call chkcmdkey('fields')

c 		Check for bogus and line inputs
        call chkcmdkey ('bogus')

c		If we are here, the new descriptor files have been written (if requested)
c		and the output files are ready to be filled.

c		Counter for the original points.
        ocnt=0
		gcnt=0
		
c		Allocate the output arrays
		allocate(ordat(odeslines), stat=ios)		
		allocate(grdat(gdeslines), stat=ios)		
		
c		Open the output files for writing the binary pipes stream.
		open (unit=11,file=gfile,form='binary',status='unknown')
		open (unit=10,file=ofile,form='binary',status='unknown')
		
500     ios = read5(rdat, ndeslines) 

c		If end of file.        
        if(ios.eq.-1) goto 999

c		If this is an original point, write the invariant fields to the original point
c		file.
		if(rdat(weighti).EQ.0.) then

c		  Increment the original point counter.
		  ocnt=ocnt+1.
		  		  
c		  Get the original point fields from this record.		  
		  do i=1,odeslines
		  	ordat(i) = rdat(ocolnum(i))
		  enddo  
		 
		  write(10) ordat
		  
		else
		
		  gcnt=gcnt+1.
		
c		  Get the gridded point fields from this record.
		  do i=1,gdeslines-1
		  	grdat(i) = rdat(gcolnum(i))
		  enddo  
c		  Links this record to the original point.
		  grdat(gdeslines)=ocnt
		  		  
		  write(11) grdat
  
		endif
		
        goto 500

999		close(unit=10)
		close(unit=11)
		write (0,*) ocnt,' original points written to ',ofile
		write (0,*) gcnt,' gridded points written to ',gfile
		stop
		
		end

