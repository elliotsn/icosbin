c******************************************************************************
c  PROGRAM: ptrinum
c
c  PURPOSE: For each data point, using lat and lon, this program calculates 
c           the triangle number on a hierarchical triangular mesh and adds it to
c           the pipe data stream. To be represented in double precision format 
c
c  REQUIRED KEYWORDS: 
c
c       lev=char -- The level of the triangular grid (assuming division by 2
c                   of sides at each new level). Level 0 is an icosahedron.
c
c       lon=char  -- The column containing east longitude
c
c       lat=char  -- The column containing latitude
c
c  OPTIONAL KEYWORDS:
c
c       des -- The descriptor file to describe the format and the
c                   fields of each record of the input data
c                   (the default is the previous pipe).
c
c.......nodes -- To prohibit the writing of the descriptor file to
c                the next pipe.  This enables the unformatted output
c                to be redirected to an output file.
c
c.......gridlonlat -- Once the final triangle number is found, set the
c.....................lat and lon parameters to the mid-point of the triangle.
c.....................This is used to spatially grid points rather than average their 
c.....................location within triangles.
c  OUTPUT COLUMNS:
c
c    trinum1, - the triangle number of intersection between a ray drawn from the
c    trinum2    sphere's origin to the centre of the triangle. For a level 14 
c               grid the number will be 15 digits long. The number has no zeros in
c               it, so that leading zeros are not cropped if the first triangle is
c               numbered < 10. trinum1 contains the digits of the triangle number up to
c               level 5, so 8 digits. trinum2 contains the remainder.
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

		program ptrinum
		implicit none

c  Variable declaration of variables common to all pipes
		include 'pipe1.inc'
	    character*48 latparam, lonparam, levparam
		character*80 desfile, junk
        integer*4 ifirst, ifound, ndeslines, ios, iunit

c       Triangle gridding vars        
        character*2 levc
c        integer newcol
        integer lonind,latind,levind,lev,trinumlen
        integer ticos(20,3), ilev, itri, numtri, tri(4,3)
        integer i, ix, iv
        real*8 lat,lon,triaddr10,base,intri(3,3)
        real*8 vert(6,3), vicos(12,3), ray(2,3), x(3)
        real*8 firsttri, tridec, tridecdig
        real*8 trinum, trinum1, trinum2
        real*8 center(3), trilat, trilon
        logical found, intersects
         
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

        call getcmdchar('lon',1,'required',ifound, lonparam)
        call getcmdchar('lat',1,'required',ifound, latparam)
        call getcmdchar('lev',1,'required',ifound, levc)
        
c		Find the column numbers of lon and lat in the input stream
        call findcol(lonparam,cdesstitle,ndeslines, lonind)
        call findcol(latparam,cdesstitle,ndeslines, latind)

c       Add trinum1 and trinum2 to the descriptor file.
        call adddesline(ndeslines+1,'trinum1',
     &            'Triangle number 1',cdesheader,
     &                  cdesstitle,cdesltitle,ndeslines)     
        
        call adddesline(ndeslines+1,'trinum2',
     &            'Triangle number 2',cdesheader,
     &                  cdesstitle,cdesltitle,ndeslines)
        
c		Add two new fields to represent the center lat-lon of the 
c		final triangle that is determined for this point.
		call adddesline(ndeslines+1,'trilat',
     &            'Center latitude of triangle',cdesheader,
     &                  cdesstitle,cdesltitle,ndeslines)
        
        call adddesline(ndeslines+1,'trilon',
     &            'Center longitude of triangle',cdesheader,
     &                  cdesstitle,cdesltitle,ndeslines)
        
c       Check for nodes keyword in command line
        call getcmdchar('nodes',0,'optional',ifound, junk)
        ! If 'nodes' is not present, send new descriptor file down pipe
        if (ifound.eq.0) then
          call pdeswrite(6, cdesname, cdesheader, 
     &               cdesstitle,  cdesltitle,
     &               ndeslines)
        endif

c       Check command line for bogus key words
        call chkcmdkey('bogus')
		
c	    Read integer version of lev
		read(unit=levc,fmt='(I2)') lev
		
c		Check level boundaries
		if ((lev .GT. 16) .OR. (lev .LT. 0)) then	
			write(0,*) 'ptrinum: ERROR - Please ensure 0<=lev<=16'
			goto 999
	    endif	
		
c		We are dividing triangle sides by 2, so triangle addresses on the grid are 
c		determined with 4 choices at each new level.

c		Set the starting triangle grid. The starting
c		vertices (list of xyz coordinates) and triangles 
c       (list of abc vertices) of an icosahedron.
        call geticos(vicos, ticos)
        	    
c       Read data. Very important that -4, because we are adding 4 new columns
c		to the pipe.
10      ios = read5(rdat, ndeslines-4)
           if( ios .eq. -1 ) goto 999
            
c        	Get lat and lon from the input pipe
           	lat = rdat(latind)
            lon = rdat(lonind)
			
c			Set up the triangle number = 9 as the marker to avoid leading zeros when the 
c			number of the first triangle in the icosahedron is < 10.
			trinum=900.
			
c			Level 0 icosahedral grid is an icosahedron.			           
        	ilev=0

c			Convert the point in lat lon space to a 3d vector
		    call sph2cart(lat,lon, x)

c			Set up a ray travelling from the origin to this point on the surface
			do i=1,3
				ray(1,i) = 0.
				ray(2,i) = x(i)
			enddo

c			Icosahedra have 20 faces
			numtri=20
			
c			For each level up until the maximum level number.
        	do while (ilev.LE.lev)

        		if (ilev.GT.0) then
c				  Set numtri = 4, because now that the 
c				  first triangle in the starting icosahedron has beed found, all
c				  subsequent tests will be on a subgrid of 4 triangles.
				  numtri=4
				endif
        		
        		found=.FALSE.
        	    itri=1
        	            	    
c			    For each triangle in this level     	    
        	    do while ( (.NOT. found) .AND. (itri.LE.numtri) )
       	              
c                  Get this triangle into a 3x3 array to pass to the intersect check
				   if (ilev.EQ.0) then
				   
				     do ix=1,3
				   	    do iv=1,3
				           intri(iv,ix) = vicos( ticos(itri,iv) ,ix)
				        enddo
				     enddo
				   else
				     do ix=1,3
				   	    do iv=1,3
				           intri(iv,ix) = vert( tri(itri,iv) ,ix)
				        enddo
				     enddo
				   endif

c        	       Test if ray intersects with triangle
        	       call intersect3dRayTriangle(ray,intri, intersects)

        	       if (intersects .EQV. .TRUE.) then
c				      If there's an intersection, add this triangle number to the output.			

        	          if (ilev.EQ.0) then
c					    If we are in level zero then trinum=900 and we just add the first
c					    triangle number in the icosahedron to this.
        	            trinum=trinum+DBLE(itri)
        	          else
c						For all levels > 0 we just multiply the triangle number by 10 to
c						make a new trailing zero, then add the current triangle number
c						(which will always be < 10) to it.
        	            trinum=(trinum*10.)+DBLE(itri)
					  endif
					   
c 					  Store the first 5 levels (8 digits) in trinum1, and the remainder 
c					  in trinum2.
					  if (ilev.EQ.5) then
					  	trinum1=trinum
					    trinum=0.	  
					  endif
        	       
c				      Get the new triangles
        	          call divTriBy2(intri,  tri,vert)
        	       
        	          found=.TRUE.
        	       endif

       			   itri=itri+1
        	    enddo

        	    ilev=ilev+1
        	enddo
			
c 			If the trinum is <= 5 then trinum2 is empty
			if (lev .LT. 5) then
				trinum1=trinum
				trinum2=0.
			else
c				The remainder of the triangle number is stored in trinum2.
				trinum2=trinum
			endif

c			Find the center point of this triangle
       		do ix=1,3 
       			center(ix)=(vert(1,ix)+vert(2,ix)+vert(3,ix))/3
       		enddo

c			Convert it to a lat-lon on the sphere.
			call cart2sph(center, trilat,trilon)			

c			Now convert trinum to base 10 from base 4
c			call base2dec(trinum,base, tridec)
			
c			How many digits is tridec?
c			tridecdig =	aint(log10(tridec))+1
			
c			Get the final triangle number. Digits are:
c			Position	Value			Meaning 
c			1			9				Leading 9 so that no leading zeros are lost.
c			2-3			0-9				Number in icosahedron	
c			4-end		0-9				Decimal representation of base 4 triangle number.
c										Convert the number formed by these digits to base 
c										4 in order to retrieve the actual triangle number.
c			trinum = (firsttri*(10**tridecdig)) + tridec
c            trinum = tridecdig
c            trinum=tridec
c          Multiply first triangle number by the appropriate power of 10 to add the 
c		   required number of trailing zeros to accomodate triaddr10.
c           trinum=dble(firsttri)*10**(aint(log10(triaddr10))+1)
c     &	+triaddr10
           
c		   Send the trinums and triangle center lat-lon down the pipe.
           rdat(ndeslines-3) = trinum1
           rdat(ndeslines-2) = trinum2

           rdat(ndeslines-1) = trilat
           rdat(ndeslines)   = trilon
           
           ios = write6(rdat, ndeslines)
	
        goto 10
999     continue
        endfile (unit=6)
        ios=close5()
        stop
        end

c	Converts a real*8 integer in base b to base 10.
c	Works when 1 < b < 10
c	Note that this is not true base b, because there are not allowed to be any zeros
c   in the input number. We subtract 1 from thisn in order to make it true base b before
c	converting to base 10.
	  subroutine base2dec(in,base, out)

	  real*8 in, base, out, thisn
	  integer n, i
	  
      n=aint(log10(in))+1.
      out=0.
      do i=1,n
          thisn=aint(in/10.**(n-i)) - aint(in/10.**(n-i+1))*10. - 1.
          out=out+thisn*base**(n-i)
      enddo
	  return
	  end

c Program to convert longitude and latitude into cartesian coords
c on a sphere.
c
c Input: 
c    lon,lat    latitude and longitude
c Output:
c    x(3)      cartesian coords, x, y, z
	  subroutine sph2cart(lat,lon,x)
	  real*8 x(3),lat,lon,pi,deg2rad,latr,lonr
	  pi=acos(0.0)
      deg2rad=pi/180.
c calc cartesian coords  
	  latr = lat*deg2rad
  	  lonr = lon*deg2rad
  	  x(1) = cos(latr)*cos(lonr)
  	  x(2) = cos(latr)*sin(lonr)
  	  x(3) = sin(latr)
	  return
	  end

c convert spherical to cartesian coordinates
	  subroutine cart2sph(x, theta,phi)
c  This routine is part of the International Astronomical Union's
c  SOFA (Standards of Fundamental Astronomy) software collection.
c
c  Status:  vector/matrix support routine.
c
c  Given:
c     x(3)     3D vector
c
c  Returned:
c     theta    latitude angle  (degrees)
c     phi      longitude angle (degrees)
c
c  Notes:
c
c  1) x can have any magnitude; only its direction is used.
c  2) If x is null, zero theta and phi are returned.
c  3) At either pole, zero theta is returned.
c-----------------------------------------------------------------------
      real*8 x(3), theta, phi, d, pi, rad2deg

	  pi=acos(0.0)
      rad2deg=180./pi

      d = x(1)*x(1) + x(2)*x(2)

      if ( d .EQ. 0D0 ) then
         phi = 0D0
      else
         phi = atan2(x(2),x(1))
      end if

      if ( x(3) .EQ. 0D0 ) then
         theta = 0D0
      else
         theta = atan2(x(3),sqrt(d))
      end if
      
      phi=phi*rad2deg
      theta=theta*rad2deg
            
	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C cross product of two vectors
      subroutine cross(x, y, a)
      real*8 x(3), y(3), a(3)
      a(1) = x(2)*y(3) - x(3)*y(2)
      a(2) = x(3)*y(1) - x(1)*y(3)
      a(3) = x(1)*y(2) - x(2)*y(1)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dot product of two vectors
      subroutine dot(x, y, a)
      real*8 x(3), y(3), a
      a = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  normalize vectors
      subroutine normalize(x, a)
      real*8 x(3), a(3), b
      b = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      a(1) = x(1)/b
      a(2) = x(2)/b
      a(3) = x(3)/b
      return
      end

c
c Function to divide a 3D triangle into 4 smaller triangles and normalise the
c resulting triangles onto the unit sphere. 
c
	  subroutine divTriBy2(intri, outtri,outvert)
      
      real*8 intri(3,3), outvert(6,3), tmp(3), tmpnew(3)
      integer outtri(4,3)
      
c     Put 3 old vertices into new array, these are keepers.    
      do iv=1,3
         do ix=1,3 
         	outvert(iv,ix)=intri(iv,ix)
         enddo
      enddo
      
c     Find 3 midpoints of input triangle based on old vertices.
      do ix=1,3
         outvert(4,ix)=(intri(1,ix)+intri(2,ix))/2.
         outvert(5,ix)=(intri(2,ix)+intri(3,ix))/2.
         outvert(6,ix)=(intri(3,ix)+intri(1,ix))/2.
      enddo
      
c     New vertices need to to be normalized onto the unit sphere. Old ones already are.
      do iv=4,6
      	 do ix=1,3
	        tmp(ix) = outvert(iv,ix)
		 enddo
		 call normalize(tmp, tmpnew)   
		 do ix=1,3
	        outvert(iv,ix) = tmpnew(ix)
		 enddo
      enddo
      
c     4 new triangles defined by 6 vertices:
      outtri(1,1) = 1
      outtri(1,2) = 4
      outtri(1,3) = 6
      
      outtri(2,1) = 4
      outtri(2,2) = 2
      outtri(2,3) = 5
      
      outtri(3,1) = 6
      outtri(3,2) = 5
      outtri(3,3) = 3
      
      outtri(4,1) = 5
      outtri(4,2) = 6
      outtri(4,3) = 4
      
      return
      end

c	  Subroutine to return the starting vertices and triangle number
c	  matrices of an icosahedron	
      subroutine geticos(v, t)
	  	
	    real*8 v(12,3), pi, phi
	    integer t(20,3)
	  
        pi=acos(0.0)
c       Golden ratio
        phi=2.*cos(pi/5.)
	  	  
c		These methods don't work in F77.
c	    vertex(1,1:3)=(0.,phi,1.) 
c        vertex(2,1:3)=(0.,-phi,1.)
c        vertex(3,1:3)=(0.,phi,-1.)
c        vertex(4,1:3)=(0.,-phi,-1.)
c        vertex(5,1:3)=(1.,0.,phi)
c        vertex(6,1:3)=(-1.,0.,phi)
c        vertex(7,1:3)=(1.,0.,-phi)
c        vertex(8,1:3)=(-1.,0.,-phi)
c        vertex(9,1:3)=(phi,1.,0.)
c        vertex(10,1:3)=(-phi,1.,0.)
c        vertex(11,1:3)=(phi,-1.,0.)
c        vertex(12,1:3)=(-phi,-1.,0.)
c	    data vertex/0.,0.,0.,0.,1.,-1.,1.,-1.,phi,-phi,phi,-phi,
c	 &  phi,-phi,phi,-phi,0.,0.,0.,0.,1.,1.,-1.,-1.,
c	 &  1.,1.,-1.,-1.,phi,phi,-phi,-phi,0.,0.,0.,0./
	  
c       12 vertices of the starting icosahedron, x, y, z.
        v(1,1)=0.
        v(1,2)=phi
        v(1,3)=1.
        
        v(2,1)=0.
        v(2,2)=-phi
        v(2,3)=1.
        
        v(3,1)=0.
        v(3,2)=phi
        v(3,3)=-1.
        
        v(4,1)=0.
        v(4,2)=-phi
        v(4,3)=-1.
        
        v(5,1)=1.
        v(5,2)=0.
        v(5,3)=phi
        
        v(6,1)=-1.
        v(6,2)=0.
        v(6,3)=phi
        
        v(7,1)=1.
        v(7,2)=0.
        v(7,3)=-phi
        
        v(8,1)=-1.
        v(8,2)=0.
        v(8,3)=-phi
        
        v(9,1)=phi
        v(9,2)=1.
        v(9,3)=0.
        
        v(10,1)=-phi
        v(10,2)=1.
        v(10,3)=0.
        
        v(11,1)=phi
        v(11,2)=-1.
        v(11,3)=0.
        
        v(12,1)=-phi
        v(12,2)=-1.
        v(12,3)=0.
        

c 		The 20 starting triangles, defined by thier vertices
    	t(1,1) = 2
		t(1,2) = 4
		t(1,3) = 11
		
		t(2,1) = 5
		t(2,2) = 2
		t(2,3) = 11
		
		t(3,1) = 9
		t(3,2) = 5
		t(3,3) = 11
			
		t(4,1) = 7
		t(4,2) = 9
		t(4,3) = 11
		
		t(5,1) = 11
		t(5,2) = 7
		t(5,3) = 4
		
		t(6,1) = 4
		t(6,2) = 2
		t(6,3) = 12
		
		t(7,1) = 6
		t(7,2) = 12
		t(7,3) = 2
		
		t(8,1) = 2
		t(8,2) = 5
		t(8,3) = 6
		
		t(9,1) = 1
		t(9,2) = 6
		t(9,3) = 5
		
		t(10,1) = 5
		t(10,2) = 9
		t(10,3) = 1

		t(11,1) = 3
		t(11,2) = 1
		t(11,3) = 9

		t(12,1) = 9
		t(12,2) = 7
		t(12,3) = 3

		t(13,1) = 8
		t(13,2) = 3
		t(13,3) = 7

		t(14,1) = 7
		t(14,2) = 4
		t(14,3) = 8

		t(15,1) = 12
		t(15,2) = 8
		t(15,3) = 4

		t(16,1) = 12
		t(16,2) = 6
		t(16,3) = 10

		t(17,1) = 6
		t(17,2) = 1
		t(17,3) = 10
		
		t(18,1) = 1
		t(18,2) = 3
		t(18,3) = 10

		t(19,1) = 3
		t(19,2) = 8
		t(19,3) = 10

		t(20,1) = 10
		t(20,2) = 12
		t(20,3) = 8
		
      return
      end
      
c Function to determine if a ray intersects with a triangle in 3D space, 
c using parametric coordinate system.
c
c Original in C++ by Dan Sunday
c Adapted for F77 by Elliot Sefton-Nash 2013/04/08
c
c Input:
c
c   ray - real*8 array of two points defining a vector a to b in the form
c			x y z
c	      a * * *
c 		  b * * *
c
c   tri - a real*8 3x3 array with row number (first index) denoting
c         vertex and column number (second index) referring to unit 
c         vector:
c                        x y z
c                      a * * *
c                      b * * *
c                      c * * *
c 
c Output:
c
c   intersects - logical returned to indicate intersect or not
c
	  subroutine intersect3dRayTriangle(ray, tri, intersects)

	  real*8 a, b, u(3), v(3), tri(3,3), w0(3), direc(3), ray(2,3)
	  real*8 uu, vv, uv, wv, wu, D, s, t, n(3), w(3), Ip(3)
	  integer i, count
	  logical intersects

c     Vector    dir, w0, w;            ray vectors
c     float     r, a, b;               params to calc ray-plane intersect
c
c	  u and v are vectors in the triangle along vertices v2 and v3 from v1.
c
c     Get triangle edge vectors.
      do i=1,3
    	u(i) = tri(2,i) - tri(1,i)
    	v(i) = tri(3,i) - tri(1,i) 
	  enddo
c	  Get plane normal	  
      call cross(u,v, n)
        
c     If the triangle is degenerate
	  count=0
      do i=1,3
         if (n(i).EQ.0.) then
        	count = count + 1    
         endif
      enddo
      if (count.EQ.3) then
        intersects = .FALSE.
	    return
	  endif
					      
      do i=1,3
c       Ray direction vector
        direc(i) = ray(2,i) - ray(1,i)       
        w0(i) = ray(1,i) - tri(1,i)
      enddo
      
      call dot(n,w0, a)
      a = -a
      call dot(n,direc, b)
      
      if (abs(b).LT.0.000000000001) then
c		    % Ray is disjoint from or parallel to triangle plane
		    intersects = .FALSE.
            return         
      endif
    
c 	  Get intersect point of ray with triangle plane
      r = a/b
      
      if (r.LT.0) then
c		 % Ray goes away from triangle         
         intersects = .FALSE.
         return
      endif
      
c     For a segment, also test if (r > 1.0) => no intersect
c     Intersect point of ray and plane
      do i=1,3
         Ip(i) = ray(1,i) + r*direc(i)
      enddo
    
c     Is the intersect point inside the triangle?
      call dot(u,u, uu)
      call dot(u,v, uv)
	  call dot(v,v, vv)
    
      do i=1,3
         w(i) = Ip(i) - tri(1,i)
      enddo
      
      call dot(w,u, wu)
      call dot(w,v, wv)
      
      D = uv * uv - uu * vv
    
c 	  Get and test parametric coords
      s = (uv * wv - vv * wu) / D
    
      if ((s .LE. 0.) .OR. (s .GT. 1.)) then
c        I is outside T, i.e. off the end of the triangle side.
         intersects = .FALSE.
         return
      endif
    
      t = (uv * wu - uu * wv) / D
      if ((t.LT.0.) .OR. ((s + t).GT.1.)) then
c        I is outside T
         intersects = .FALSE.
         return
      endif
      
c     If we are here, intersection point I is inside the triangle.
      intersects = .TRUE.
      return
      end