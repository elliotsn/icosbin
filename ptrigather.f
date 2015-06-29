c******************************************************************************
c  PROGRAM: ptrigather
c
c  PURPOSE: This program is used after adding triangle numbers (bin numbers) to 
c			a stream of data points that has already been through pfootprint.
c			Each point is now represented by a much larger number of points describing 
c			a probability distribution of the effective field of view for that point.
c
c			It is very likely that many points from the same original measurement will 
c			have the same triangle number. This program averages all those points together
c			in order to reduce the data volume that is written out to disk.
c
c			The program adds a column to the pipe stream, 'weight', which contains
c	        the cumulative value of each measurement output. 'weight' is
c			equal to the number of data points averaged divided by the total number of
c			points used by pfootprint to model the EFOV.
c
c  REQUIRED KEYWORDS: 
c
c       points=char -- The number of points used to model each measurements EFOV.
c					   This is used to define the size of the bucket inside which
c					   points with the same triangle number are averaged. Note points
c					   must be 3 digits (100 <= points <= 999)
c
c  OPTIONAL KEYWORDS:
c
c       des -- The descriptor file to describe the format and the
c                   fields of each record of the input data
c                   (the default is the previous pipe).
c
c       nodes -- To prohibit the writing of the descriptor file to
c                the next pipe. This enables the unformatted output
c                to be redirected to an output file.
c
c  EXTRA OUTPUT COLUMNS:
c
c       weight -- The fraction of the original measurement that this point represents.
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

		program ptrigather
		implicit none

c  Variable declaration of variables common to all pipes
		include 'pipe1.inc'
	    character*48 p_trinum1, p_trinum2
c	    character*48 p_radiance, p_tb, p_talt, p_tslope, p_taz
		character*48 p_temis, p_tinc, p_tinsol, p_tphase
		character*48 p_tlat, p_tlon
		character*80 desfile, junk
        integer*4 ifirst, ifound, ndeslines, ios, iunit, ndeslinesO
		integer*4 bktcnt, points, outcnt
		
c		integer*4 ind_rad, ind_tb, ind_talt, ind_tslope, ind_taz
		integer*4 ind_temis, ind_tinc, ind_tinsol, ind_tphase
		integer*4 ind_tlat, ind_tlon
		integer*4 it1, it2, naind, nmind, stat
		
		integer*4 iin,iout,im,iav,irow,icol

		logical found, match
		
c	Add one column to output bucket because we add 'weight' to it.
		real*8, allocatable :: out(:,:)
		real*8, allocatable :: in(:,:)
		real*8, allocatable :: rdatO(:)
		
		integer*4, allocatable :: aind(:)
		integer*4, allocatable :: mind(:)
				
		logical done
		
        character*3 pointsc
         
c Check for descriptor file in command line
        call getcmdchar('des',1,'optional',ifound,desfile)
        if( ifound .eq. 1 ) then
           iunit=7 ! If descriptor file found, read from unit 7
        else
           iunit=5 ! Else, get the descriptor from standard input
        endif
        call pdesread(iunit, desfile,  cdesheader,
     &               cdesstitle,  cdesltitle,
     &               ndeslines)

c Check the required keywords are present.
        call getcmdchar('trinum1',1,'required',ifound, p_trinum1)
        call getcmdchar('trinum2',1,'required',ifound, p_trinum2)
        call getcmdchar('points',1,'required',ifound, pointsc)
        
c        call getcmdchar('radiance',1,'required',ifound, p_radiance)
c        call getcmdchar('tb',1,'required',ifound, p_tb)
c		call getcmdchar('talt',1,'required',ifound, p_talt)        
c        call getcmdchar('tslope',1,'required',ifound, p_tslope)
        call getcmdchar('temis',1,'required',ifound, p_temis)
        call getcmdchar('tinc',1,'required',ifound, p_tinc)
c        call getcmdchar('tinsol',1,'required',ifound, p_tinsol)
        call getcmdchar('tphase',1,'required',ifound, p_tphase)
c        call getcmdchar('taz',1,'required',ifound, p_taz)
        
        call getcmdchar('tlat',1,'required',ifound, p_tlat)
        call getcmdchar('tlon',1,'required',ifound, p_tlon)
                
c		Find the column numbers of trinum1 and trinum2 in the input stream
        call findcol(p_trinum1,cdesstitle,ndeslines, it1)
        call findcol(p_trinum2,cdesstitle,ndeslines, it2)

c		Find the column numbers of quantities that must be averaged in the input stream
c		call findcol(p_radiance,cdesstitle,ndeslines, ind_rad)
c		call findcol(p_tb,cdesstitle,ndeslines, ind_tb)
c		call findcol(p_talt,cdesstitle,ndeslines, ind_talt)
c		call findcol(p_tslope,cdesstitle,ndeslines, ind_tslope)
		call findcol(p_temis,cdesstitle,ndeslines, ind_temis)
		call findcol(p_tinc,cdesstitle,ndeslines, ind_tinc)
c		call findcol(p_tinsol,cdesstitle,ndeslines, ind_tinsol)
		call findcol(p_tphase,cdesstitle,ndeslines, ind_tphase)
c		call findcol(p_taz,cdesstitle,ndeslines, ind_taz)
		
		call findcol(p_tlat,cdesstitle,ndeslines, ind_tlat)
		call findcol(p_tlon,cdesstitle,ndeslines, ind_tlon)

c       Add weight to the descriptor file.
        ndeslinesO=ndeslines
        call adddesline(ndeslines+1,'weight',
     &            'Weight',cdesheader,
     &                  cdesstitle,cdesltitle,ndeslinesO)

c       Check for nodes keyword in command line
        call getcmdchar('nodes',0,'optional',ifound, junk)
        ! If 'nodes' is not present, send new descriptor file down pipe
        if (ifound.eq.0) then
          call pdeswrite(6, cdesname, cdesheader, 
     &               cdesstitle,  cdesltitle,
     &               ndeslinesO)
        endif

c       Check command line for bogus key words
        call chkcmdkey('bogus')
		
c	    Read integer version of points
		read(unit=pointsc,fmt='(I3)') points

c		Allocate the bucket arrays now we know how big they'll be
		allocate(out(points,ndeslinesO), stat=ios)
		allocate(in(points,ndeslines), stat=ios)

c		Allocate the output pipe array, one larger than the input because we add 'weight'
		allocate(rdatO(ndeslinesO), stat=ios)

c		Initialize the buckets. ndeslinesO is one more than the incoming pipe, because
c		we added the column, weight to the output pipe.
		do irow=1,points
			do icol=1,ndeslines
				in(irow,icol)=0.
			enddo
			do icol=1,ndeslinesO
				out(irow,icol)=0.
			enddo
		enddo

c		The number of points in the input bucket is zero.
		bktcnt=0

c		The number of collected points in the output bucket is zero
		outcnt=0

c		The indices of the columns to average for gathered points.
c		naind=9
		naind=5
		allocate(aind(naind), stat=ios)
		aind(1)=ind_temis
		aind(2)=ind_tinc
		aind(3)=ind_tphase
		aind(4)=ind_tlat
		aind(5)=ind_tlon
c		aind(6)=ind_talt
c		aind(7)=ind_tslope
c		aind(8)=ind_taz
c		aind(9)=ind_tinsol

		nmind=2
		allocate(mind(nmind), stat=ios)
		mind(1)=it1
		mind(2)=it2

c       Read data.
10      ios = read5(rdat, ndeslines)
c		   If there are no more points           
           if( ios .eq. -1 ) goto 999
			
c		   One more point in the bucket.	
		   bktcnt = bktcnt + 1

c	       The first point that's output by pfootprint is always the original Diviner 
c		   measurement. We can use this as a marker by setting the weight to 0, indicating
c		   that this point shouldn't be binned. All other gathered points will have 
c		   weight > 0.
		   if (bktcnt.EQ.1) then

c		     Set the weight for the first point to zero and send it down the pipe.
		     do icol=1,ndeslines
		       rdatO(icol) = rdat(icol)
		     enddo
		     rdatO(ndeslinesO)=0.
		     ios = write6(rdatO, ndeslinesO)
		   
		   else
c		     Add this point to the in bucket.
	  	     do icol=1,ndeslines
			   in(bktcnt-1,icol) = rdat(icol)
		     enddo
		   endif
		   
c		   Is the bucket full? If so, gather points in 'in' with equal trinums into single
c		   records then put them in out.
		   if(bktcnt .EQ. points+1) then
c		     Initialize the out bucket.		   
		     do irow=1,points
		       do icol=1,ndeslinesO
		         out(irow,icol)=0.
		       enddo
		     enddo
c
c	******* MERGE ALL 'IN' POINTS WITH SAME TRINUM AND PUT MEANS IN 'OUT' ********
c
c		   For each row in 'in'.
		   do iin=1,points
		    
c		     Check if the required columns match any rows in 'out'. They should never
c		     match more than one, so if a match is found it is assumed it's the only
c		     one.
		     iout=0
		     found=.FALSE.
	  	     do while ((.NOT. found) .AND. (iout.LE.outcnt))
		       iout=iout+1
c			   Test for this row in the out array if all the required columns match this
c			   point in the in array.
			   im=0
			   match=.TRUE.
			   do while ((im .LT. nmind).AND.(match .EQV. .TRUE.))
		         im=im+1
		         if(in(iin,mind(im)).NE.out(iout,mind(im))) then
			       match=.FALSE.
			     endif
			   enddo
c			   If we're here and match is still true, then we've found the one we're looking for.			
			   if(match) found=.TRUE.
		     enddo
		  
c		     For this row in 'in', we add it to the appropriate row in 'out'.
c		     We add to the running total if we've found a matching row in the out array.
		     if(found) then
c			   Increment the counter for this point in 'out', which is always in the final 
c			   column, at cols+1
		       out(iout,ndeslinesO)=out(iout,ndeslinesO)+1.
c			   Add running total to every column that we must.
			   do iav=1,naind
		out(iout,aind(iav))=out(iout,aind(iav))+in(iin,aind(iav))
			   enddo
		     else
c			   This is the first time this row has been encountered, 
c              copy it across to 'out'		  
		       outcnt=outcnt+1
		       do icol=1,ndeslines
		         out(outcnt,icol)=in(iin,icol)
		       enddo
c			   Increment the counter for this new row.
		  	   out(outcnt,ndeslinesO)=out(outcnt,ndeslinesO)+1
			 endif
		   enddo

c		   Now the gathered points should be copied across the out array, we still need to 
c		   divide the summed points by the number of original points, to get the average. 
           do iout=1,outcnt
			
c		     Divide running total by number of points contributing to get average.
		     do iav=1,naind
		       out(iout,aind(iav))=out(iout,aind(iav))/out(iout,ndeslinesO)
             enddo

c		     Turn counts into weight based on the number of points in the inbucket, which is
c		     held in 'points'.
		     out(iout,ndeslinesO)=out(iout,ndeslinesO)/points
           enddo
				
c		   Finally, write each gathered point in out to the pipe stream.
		   do irow=1,outcnt 
c		     Send this point in out down the pipe, put it in a tmp out variable, 
c		     rdatO, first.
		     do icol=1,ndeslinesO
		       rdatO(icol) = out(irow,icol)
		     enddo
		     ios = write6(rdatO, ndeslinesO)
	       enddo

c		     Reset the counters and start again.
		     bktcnt=0
       	     outcnt=0
           endif
       	goto 10
		  
999     continue

        endfile (unit=6)
        ios=close5()
        stop
        end
        
c				  All times COULD lie on boundaries, BUT if all from one point should have
c				  same housekeeping data. This relies on pfootprint outputting data in
c				  point sequential order.
c				  So shouldn't need to average.
c				   date			- Should be same
c				   month		- Should be same
c				   year 		- Should be same
c				   hour 		- Should be same
c				   minute		- Should be same
c				   second		- Should be same
c				      
c				   jdate		- Should be same
c				   orbit		- Should be same
c				   
c				   sundst		- Should be same
c				   sclk			- Should be same
c				   scrad		- Should be same
c				   scalt		- Should be same
c				   el_cmd		- Should be same
c				   az_cmd		- Should be same
c				   c			- Should be same
c				   det			- Should be same
c				   radiance		- Average
c				   tb			- Average
c				   cemis		- Should be same
c				   csunzen		- Should be same
c				   csunazi		- Should be same
c				   cloctime		- Should be same
c				   qca			- Should be same
c				   qge			- Should be same
c				   qmi			- Should be same
c				   af			- Should be same
c				   talt			- Average OR don't include because it should be in DEM
c				   tslope		- Average OR don't include because it should be in DEM
c				   temis		- Average
c				   tinc			- Average
c				   tinsol		- Average
c				   tphase		- Average
c				   taz			- Average OR don't include because it should be in DEM
c				   tlat			- Average - because within one triangle it's more correct to
c											save the average geographic position of them
c											rather than the center of the triangle.
c				   tlon			- Average -  "
c				   trinum1
c				   trinum2		
c				   trilat       - Don't need - Geographic center of the triangle not 
c											   representative of geographic center of 
c											   point spread in that triangle.
c				   trilon		- Don't need - 	"