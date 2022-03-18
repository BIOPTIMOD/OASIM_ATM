      program dat2NetCDF

!  Computes monthly mean non-diurnal spectral irradiance W/m2.
!  Accounts for surface reflectance so this is just below the
c  surface (0-).  

      parameter(im=360,jm=180)
      parameter(alon0=-180.0,alat0=-90.0)
      parameter(dlon=1.0,dlat=1.0)
      parameter(nlt=33)
      CHARACTER(len=100) :: filein
      CHARACTER(len=100) :: fileout
      CHARACTER(len=17) :: date_start
      CHARACTER(len=17) :: date_end
      real taua(im,jm,nlt),asymp(im,jm,nlt),ssalb(im,jm,nlt)   !
      real myarray(im,jm,nlt,3)
      byte msk
      common /bmsk/ msk(im,jm)
      data flagl /999.0E9/

      call getarg(1, filein)
      write (*,*) 'filein ', filein
      call getarg(2, fileout)
      write (*,*) 'fileout ', fileout
      call getarg(3, date_start)
      write (*,*) 'date_start', date_start
      call getarg(4, date_end)
      write (*,*) 'date_end', date_end

      open(3,file=TRIM(filein),status='old',form='unformatted')
      read(3) myarray(:,:,:,:)
      close(3)
      do nl = 1,nlt
        do j = 1,jm
          do i = 1,im
               taua(i,j,nl)  = myarray(i,j,nl,1)
               asymp(i,j,nl) = myarray(i,j,nl,2)
               ssalb(i,j,nl) = myarray(i,j,nl,3)
          enddo
        enddo
      enddo
      call WRITE_AVE_2D(TRIM(fileout),'taua','asymp','ssalb'
     &                 ,date_start
     &                 ,date_end, taua, asymp, ssalb)

      end
      
