program THALASSA_DRIVER_TOL
! Description:
!    Driver program for Thalassa. It launches a batch of propagations using
!    distinct values of the integration tolerance, and measures the CPU time
!    for each propagation.
!
! ==============================================================================

! MODULES
use KINDS, only: dk

! VARIABLES
implicit none ! "you got it baby"
! Tolerance intervals
real(dk)  ::  tolmin,tolmax,delta
real(dk),allocatable  ::  tol(:)
integer   ::  itol,ntol
! Input file stuff
integer              ::  tic,toc,rate
integer,parameter    ::  id_input = 10, tol_line = 28, nlines = 40
integer              ::  iline,ioflag
character(len=256)   ::  input_file(1:nlines)
character(len=8)     ::  date
character(len=6)     ::  time

! ==============================================================================

! Assign tolerance interval
write(*,'(a)') 'THALASSA - DRIVER (TOL)'
write(*,'(a)') 'Assign tol_max:'
read(*,*) tolmax
write(*,'(a)') 'Assign tol_min:'
read(*,*) tolmin
write(*,'(a)') 'Assign ntol:'
read(*,*) ntol

! Read the whole input file
open(unit=id_input,file='in/input.txt',status='old',action='readwrite')
ioflag = 0
iline  = 0
do while (ioflag == 0)
    iline = iline + 1
    read(unit=id_input,iostat=ioflag,fmt='(a256)') input_file(iline)

end do
rewind(unit=id_input)  ! for later

! Space the tolerance vector logarithmically
allocate(tol(1:ntol))
delta = (log10(tolmin) - log10(tolmax))/(ntol - 1)
tol(1) = tolmax

call DATE_AND_TIME(date,time)

! Assign tolerance value to input file and run
do itol = 1,ntol
  tol(itol) = 10._dk**(log10(tolmax) + (itol - 1)*delta)

  ! Input file: modify tolerance and output dir
  write(input_file(tol_line),'(a10,es22.15)') 'tol:      ',tol(itol)
  write(input_file(nlines-1),'(4(a),i2.2,a)')&
  &'out: /home/davide/Documents/PhD/02.Progetti/THESS/data/Thalassa/',&
  &date,time,'/t',itol,'/'
  write(input_file(nlines),'(4(a))')&
  &'stat: /home/davide/Documents/PhD/02.Progetti/THESS/data/Thalassa/',&
  &date,time,'/'

  write(id_input,'(a)') (input_file(iline), iline = 1,nlines )
  rewind(id_input)

  ! Call Thalassa main program
  call SYSTEM('./thalassa.x')



end do

deallocate(tol)
close(id_input)

end program THALASSA_DRIVER_TOL
