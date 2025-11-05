module test_module
    implicit none
    contains

    subroutine test_module_hello_world()
        write(*,*) 'hello world'
    end subroutine


end module test_module

program main
    implicit none

end program main