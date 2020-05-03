
program option

 use iso_c_binding
 implicit none

 logical :: option1, option2, option3
 integer :: option1_read, option2_read, option3_read
 integer :: opt

 option1 = .false.
 option2 = .false.
 option3 = .true.

 ! set options
 opt= 0
 if (option1) opt = ibset(opt,0)
 if (option2) opt = ibset(opt,1)
 if (option3) opt = ibset(opt,2)

 write(*,*) opt

 ! read options
 option1_read = ibits(opt,0,1)
 option2_read = ibits(opt,1,1)
 option3_read = ibits(opt,2,1)

 write(*,*) option1_read, option2_read, option3_read

end program option
