       SUBROUTINE EXISTVAR(fileNetCDF, varname,B)
       USE netcdf
       IMPLICIT NONE
       character, INTENT(IN) :: fileNetCDF*(*) ,varname*(*)
       LOGICAL,   INTENT(OUT)::B
       ! local
       integer stat,ncid,VARid,counter

       counter=0
       B=.false.
         stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)  
       call handle_err1(stat, counter,FileNetCDF)
         stat = nf90_inq_varid (ncid, varname, VARid)

         if(stat .ne. nf90_NoErr)  then
           !write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
           B = .false.
         else
            B=.true.
        endif
        stat = nf90_close(ncid)                           
       call handle_err1(stat, counter,FileNetCDF)
       END SUBROUTINE EXISTVAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE WRITE_AVE_2D(fileNetCDF, VAR1, VAR2, VAR3, VAR4, VAR5, datefrom, dateTo, M1, M2, M3, M4, M5)
       USE netcdf
       IMPLICIT NONE
       integer, parameter :: jpiglo = 360 
       integer, parameter :: jpjglo = 180
       integer, parameter :: ntime  = 12

       character*(*),intent(in) :: fileNetCDF
       character(LEN=17),intent(in) :: datefrom, dateTo
       real,intent(in),dimension(jpiglo, jpjglo, ntime) :: M1, M2, M3, M4, M5

       character(LEN=20) :: VAR1, VAR2, VAR3, VAR4, VAR5
       integer :: istart,iend

       integer :: ji,jj
       integer :: s, nc, counter
       integer :: timid, yid, xid, nid
       integer :: idvartime, idphit, idlamt, idVAR1, idVAR2, idVAR3, idVAR4, idVAR5

       real :: lat_actual_range(2), lon_actual_range(2)
       REAL :: totglamt(jpiglo,jpjglo),  totgphit(jpiglo,jpjglo)
         lon_actual_range=(/-180.0 , 180.0   /)
         lat_actual_range=(/-90.0  , 90.0    /)


        do jj=1,jpjglo
            do ji=1,jpiglo
                totglamt(ji,jj) = -180.0 + REAL(ji) -1.0
                totgphit(ji,jj) = -90.0  + REAL(jj) -1.0
            enddo
        enddo

        s = nf90_create(fileNetCDF, or(nf90_clobber,NF90_HDF5), nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, nf90_global, 'Convenctions'   ,'COARDS')
        s = nf90_put_att(nc, nf90_global, 'DateStart'     , datefrom)
        s = nf90_put_att(nc, nf90_global, 'Date__End'     ,   dateTo)

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'lon'           , jpiglo,  xid)
        s= nf90_def_dim(nc,'lat'           , jpjglo,  yid)
        s= nf90_def_dim(nc,'time'          ,  ntime,timid)

        ! ********** VARIABLES *****************
        !s = nf90_def_var(nc,'time',         nf90_double,(/timid/),       idvartime)
        s = nf90_def_var(nc,'lat'   ,       nf90_float, (/yid/),            idphit)
        s = nf90_def_var(nc,'lon'   ,       nf90_float, (/xid/),            idlamt)

        s = nf90_def_var(nc,trim(VAR1) ,        nf90_float, (/xid,yid,timid/),  idVAR1)

        s = nf90_def_var(nc,trim(VAR2) ,        nf90_float, (/xid,yid,timid/),  idVAR2)

        s = nf90_def_var(nc,trim(VAR3) ,        nf90_float, (/xid,yid,timid/),  idVAR3)

        s = nf90_def_var(nc,trim(VAR4) ,        nf90_float, (/xid,yid,timid/),  idVAR4)

        s = nf90_def_var(nc,trim(VAR5) ,        nf90_float, (/xid,yid,timid/),  idVAR5)

        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idphit,'actual_range' ,lat_actual_range)

        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')
        s = nf90_put_att(nc,idlamt,'actual_range' ,lon_actual_range)

        s = nf90_put_att(nc,idVAR1, 'long_name'    ,trim(VAR1))
        s = nf90_put_att(nc,idVAR1, 'missing_value' ,-1)

        s = nf90_put_att(nc,idVAR2, 'long_name'    ,trim(VAR2))
        s = nf90_put_att(nc,idVAR2, 'missing_value' ,-1)

        s = nf90_put_att(nc,idVAR3, 'long_name'    ,trim(VAR3))
        s = nf90_put_att(nc,idVAR3, 'missing_value' ,-1)

        s = nf90_put_att(nc,idVAR4, 'long_name'    ,trim(VAR4))
        s = nf90_put_att(nc,idVAR4, 'missing_value' ,-1)

        s = nf90_put_att(nc,idVAR5, 'long_name'    ,trim(VAR5))
        s = nf90_put_att(nc,idVAR5, 'missing_value' ,-1)

        s =nf90_enddef(nc)

        counter=0
        s = nf90_put_var(nc, idlamt,  REAL(totglamt(:,jpjglo),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idphit,  REAL(totgphit(jpiglo,:),4) )
        call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR1  ,  M1) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR2 ,  M2) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR3  ,  M3) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR4 ,  M4) 
       call handle_err1(s,counter,fileNetCDF)
        s = nf90_put_var(nc, idVAR5 ,  M5) 
       call handle_err1(s,counter,fileNetCDF)


        s =nf90_close(nc)

       END SUBROUTINE WRITE_AVE_2D


      !****************************************************************************
        subroutine handle_err1(status,mycount, fileNetCDF)
        USE netcdf
        integer status,mycount
        character fileNetCDF*(*)
        mycount =mycount+1
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'netcdf call',mycount,'with status = ',status
           write(*,*)  'file :', fileNetCDF
           write(*,*) nf90_strerror(status)
           write(*,*) 'Stopped'
           STOP 1
        endif

        end subroutine handle_err1


      ! **************************************************************************
        subroutine handle_err2(status,fileNetCDF,varname)
        USE netcdf
        integer status
        character fileNetCDF*(*) ,varname*(*)

        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2
