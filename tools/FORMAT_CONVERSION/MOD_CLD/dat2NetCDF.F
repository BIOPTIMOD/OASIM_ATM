      program dat2NetCDF

!  Computes monthly mean non-diurnal spectral irradiance W/m2.
!  Accounts for surface reflectance so this is just below the
c  surface (0-).  

      parameter(im=360,jm=180)
      parameter(alon0=-180.0,alat0=-90.0)
      parameter(dlon=1.0,dlat=1.0)
      parameter(time=12)
      CHARACTER(len=14) :: filein
      CHARACTER(len=13) :: fileout
      CHARACTER(len=17) :: date_start
      CHARACTER(len=17) :: date_end
      integer month
      real*4 ccovm(im,jm,time),cldtcm(im,jm,time)  
      real*4 rlwpm(im,jm,time),cdrem(im,jm,time) 
      real*4 myarray(im,jm,4,time)
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
      read(3)imc,jmc
      do month =1,12
        read(3) myarray(:,:,:,month)
      enddo
      close(3)

      do month = 1,12
        do j = 1,jm
            do i = 1,im
               ccovm(i,j,month)  = myarray(i,j,1,month)
               cldtcm(i,j,month) = myarray(i,j,2,month)
               rlwpm(i,j,month)  = myarray(i,j,3,month)
               cdrem(i,j,month)  = myarray(i,j,4,month)
            enddo
        enddo
       enddo

      call WRITE_AVE_2D(fileout,'ccovm','cldtcm', 'rlwpm','cdrem'
     &                 ,date_start
     &                 ,date_end, ccovm, cldtcm, rlwpm, cdrem)

      end
      
