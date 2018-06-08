program ConvertRIO


! how to run:  ./ConvertRIO 2016 05 12
! command line arguments: year month day (f.e. 2016 05 12)

implicit none
real :: id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20,tmp21,tmp22,tmp23,tmp24
real :: xlamb(2098),ylamb(2098),values_pm10(2098,24),values_no2(2098,24),values_o3(2098,24),values_pm25(2098,24),values_bc(2098,24),values_no(2098,24)
character*4 :: jaar
character*2 :: maand,dag
character(len=8) :: fmt
character*2 :: str
character*240 :: heading
integer :: i,j,uur,lamx_au,lamy_au,k,i2,j2,k2,xlam_tmp,ylam_tmp
real :: values_aur(81,58,6),lat_au(81,58),lon_au(81,58),lat,lon,lamx_aud,lamy_aud
CHARACTER *100 buffer  

! convert command line arguments
CALL getarg(1,buffer)
READ(BUFFER,*) jaar
CALL getarg(2,buffer)
READ(BUFFER,*) maand
CALL getarg(3,buffer)
READ(BUFFER,*) dag


write(*,*) jaar
write(*,*) maand
write(*,*) dag



OPEN(UNIT=1,FILE='lambert4x4.dat')
READ(1,*)
do i=1,2098
READ(1,*) id_tmp,xlam_tmp,ylam_tmp
 xlamb(i)=xlam_tmp
 ylamb(i)=ylam_tmp
enddo
CLOSE(1)

OPEN(UNIT=1,FILE='pm10_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_pm10(i,1)=tmp1
 values_pm10(i,2)=tmp2
 values_pm10(i,3)=tmp3
 values_pm10(i,4)=tmp4
 values_pm10(i,5)=tmp5
 values_pm10(i,6)=tmp6
 values_pm10(i,7)=tmp7
 values_pm10(i,8)=tmp8
 values_pm10(i,9)=tmp9
 values_pm10(i,10)=tmp10
 values_pm10(i,11)=tmp11
 values_pm10(i,12)=tmp12
 values_pm10(i,13)=tmp13
 values_pm10(i,14)=tmp14
 values_pm10(i,15)=tmp15
 values_pm10(i,16)=tmp16
 values_pm10(i,17)=tmp17
 values_pm10(i,18)=tmp18
 values_pm10(i,19)=tmp19
 values_pm10(i,20)=tmp20
 values_pm10(i,21)=tmp21
 values_pm10(i,22)=tmp22
 values_pm10(i,23)=tmp23
 values_pm10(i,24)=tmp24
enddo

CLOSE(1)


OPEN(UNIT=1,FILE='bc_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_bc(i,1)=tmp1
 values_bc(i,2)=tmp2
 values_bc(i,3)=tmp3
 values_bc(i,4)=tmp4
 values_bc(i,5)=tmp5
 values_bc(i,6)=tmp6
 values_bc(i,7)=tmp7
 values_bc(i,8)=tmp8
 values_bc(i,9)=tmp9
 values_bc(i,10)=tmp10
 values_bc(i,11)=tmp11
 values_bc(i,12)=tmp12
 values_bc(i,13)=tmp13
 values_bc(i,14)=tmp14
 values_bc(i,15)=tmp15
 values_bc(i,16)=tmp16
 values_bc(i,17)=tmp17
 values_bc(i,18)=tmp18
 values_bc(i,19)=tmp19
 values_bc(i,20)=tmp20
 values_bc(i,21)=tmp21
 values_bc(i,22)=tmp22
 values_bc(i,23)=tmp23
 values_bc(i,24)=tmp24
enddo

CLOSE(1)

OPEN(UNIT=1,FILE='pm25_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_pm25(i,1)=tmp1
 values_pm25(i,2)=tmp2
 values_pm25(i,3)=tmp3
 values_pm25(i,4)=tmp4
 values_pm25(i,5)=tmp5
 values_pm25(i,6)=tmp6
 values_pm25(i,7)=tmp7
 values_pm25(i,8)=tmp8
 values_pm25(i,9)=tmp9
 values_pm25(i,10)=tmp10
 values_pm25(i,11)=tmp11
 values_pm25(i,12)=tmp12
 values_pm25(i,13)=tmp13
 values_pm25(i,14)=tmp14
 values_pm25(i,15)=tmp15
 values_pm25(i,16)=tmp16
 values_pm25(i,17)=tmp17
 values_pm25(i,18)=tmp18
 values_pm25(i,19)=tmp19
 values_pm25(i,20)=tmp20
 values_pm25(i,21)=tmp21
 values_pm25(i,22)=tmp22
 values_pm25(i,23)=tmp23
 values_pm25(i,24)=tmp24
enddo


CLOSE(1)

OPEN(UNIT=1,FILE='o3_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_o3(i,1)=tmp1
 values_o3(i,2)=tmp2
 values_o3(i,3)=tmp3
 values_o3(i,4)=tmp4
 values_o3(i,5)=tmp5
 values_o3(i,6)=tmp6
 values_o3(i,7)=tmp7
 values_o3(i,8)=tmp8
 values_o3(i,9)=tmp9
 values_o3(i,10)=tmp10
 values_o3(i,11)=tmp11
 values_o3(i,12)=tmp12
 values_o3(i,13)=tmp13
 values_o3(i,14)=tmp14
 values_o3(i,15)=tmp15
 values_o3(i,16)=tmp16
 values_o3(i,17)=tmp17
 values_o3(i,18)=tmp18
 values_o3(i,19)=tmp19
 values_o3(i,20)=tmp20
 values_o3(i,21)=tmp21
 values_o3(i,22)=tmp22
 values_o3(i,23)=tmp23
 values_o3(i,24)=tmp24
enddo


CLOSE(1)


OPEN(UNIT=1,FILE='no2_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_no2(i,1)=tmp1
 values_no2(i,2)=tmp2
 values_no2(i,3)=tmp3
 values_no2(i,4)=tmp4
 values_no2(i,5)=tmp5
 values_no2(i,6)=tmp6
 values_no2(i,7)=tmp7
 values_no2(i,8)=tmp8
 values_no2(i,9)=tmp9
 values_no2(i,10)=tmp10
 values_no2(i,11)=tmp11
 values_no2(i,12)=tmp12
 values_no2(i,13)=tmp13
 values_no2(i,14)=tmp14
 values_no2(i,15)=tmp15
 values_no2(i,16)=tmp16
 values_no2(i,17)=tmp17
 values_no2(i,18)=tmp18
 values_no2(i,19)=tmp19
 values_no2(i,20)=tmp20
 values_no2(i,21)=tmp21
 values_no2(i,22)=tmp22
 values_no2(i,23)=tmp23
 values_no2(i,24)=tmp24
enddo


CLOSE(1)

OPEN(UNIT=1,FILE='no_1h_'//jaar//maand//dag//'.dat')
do i=1,2098
    READ(1,'(a240)') heading
    CALL splitsmaarop(heading,id_tmp,date_tmp1,extra,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10 &
 ,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20 &
 ,tmp21,tmp22,tmp23,tmp24)

 values_no(i,1)=tmp1
 values_no(i,2)=tmp2
 values_no(i,3)=tmp3
 values_no(i,4)=tmp4
 values_no(i,5)=tmp5
 values_no(i,6)=tmp6
 values_no(i,7)=tmp7
 values_no(i,8)=tmp8
 values_no(i,9)=tmp9
 values_no(i,10)=tmp10
 values_no(i,11)=tmp11
 values_no(i,12)=tmp12
 values_no(i,13)=tmp13
 values_no(i,14)=tmp14
 values_no(i,15)=tmp15
 values_no(i,16)=tmp16
 values_no(i,17)=tmp17
 values_no(i,18)=tmp18
 values_no(i,19)=tmp19
 values_no(i,20)=tmp20
 values_no(i,21)=tmp21
 values_no(i,22)=tmp22
 values_no(i,23)=tmp23
 values_no(i,24)=tmp24
enddo


CLOSE(1)




do uur=1,24
values_aur=0.

 fmt = '(I2.2)' ! an integer of width 5 with zeros at the left


write (str,fmt) uur-1 ! converting integer to string using a 'internal file'

 
do i2=1,81
 do j2=1,58

  lamx_au=(i2-1)*4000+0000
  lamy_au=(j2-1)*4000+18000
  do k2=1,2098
   if (lamx_au .eq. xlamb(k2) .and. lamy_au .eq. ylamb(k2)) then
    values_aur(i2,j2,1)=1.e-9*values_o3(k2,uur)
    values_aur(i2,j2,2)=1.e-9*values_no2(k2,uur)
    values_aur(i2,j2,3)=1.e-9*values_no(k2,uur)
    values_aur(i2,j2,4)=1.e-9*values_pm10(k2,uur)
    values_aur(i2,j2,5)=1.e-9*values_pm25(k2,uur)
    values_aur(i2,j2,6)=1.e-9*values_bc(k2,uur)
    
   endif
  enddo
 enddo
enddo
OPEN(unit=11, FILE=''//jaar//maand//dag//str//'Aurora_conc',form='unformatted')
 WRITE(11) values_aur
CLOSE(11)
enddo

CONTAINS

SUBROUTINE splitsmaarop(heading,id,date_tmp1,extra,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24)
implicit none
CHARACTER*240 :: heading
integer :: vorige,j,p
real :: g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,zevenentwintig(27)
real :: date_tmp1,extra,id
vorige=0
j=0
do p=1,240
	if (heading(p:p) .eq. ";") then
		j=j+1
		zevenentwintig(j)=strtoreal(heading(vorige+1:p-1))
		vorige=p
        if (j .eq. 26) then           
            zevenentwintig(27)=strtoreal(heading(p+1:240))
            exit
        endif
	endif
enddo



20 continue
id=zevenentwintig(1)
date_tmp1=zevenentwintig(2)
extra=zevenentwintig(3)
g1=zevenentwintig(4)
g2=zevenentwintig(5)
g3=zevenentwintig(6)
g4=zevenentwintig(7)
g5=zevenentwintig(8)
g6=zevenentwintig(9)
g7=zevenentwintig(10)
g8=zevenentwintig(11)
g9=zevenentwintig(12)
g10=zevenentwintig(13)
g11=zevenentwintig(14)
g12=zevenentwintig(15)
g13=zevenentwintig(16)
g14=zevenentwintig(17)
g15=zevenentwintig(18)
g16=zevenentwintig(19)
g17=zevenentwintig(20)
g18=zevenentwintig(21)
g19=zevenentwintig(22)
g20=zevenentwintig(23)
g21=zevenentwintig(24)
g22=zevenentwintig(25)
g23=zevenentwintig(26)
g24=zevenentwintig(27)



END SUBROUTINE


REAL function strtoreal(str)

implicit none

character*(*), intent(in):: str
REAL :: result

read(str,*) result
strtoreal= result

end function strtoreal



END PROGRAM

 


