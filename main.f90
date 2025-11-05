module test_module
    implicit none
    contains

    subroutine test_module_hello_world()
        write(*,*) 'hello world'
    end subroutine


end module test_module

module parameters_for_simulation
    implicit none

    !dimention: シミュレーションの次元を選択する。2->二次元、３->三次元
    !pargicle_distance: 初期の粒子間距離[m]
    !time_interval: 計算ステップの時間幅[s]
    !output_interval: 何ステップごとに粒子の位置を出力するか

    !------------change allowed--------------!
    integer,parameter :: dimention=3
    integer,parameter :: paricle_distance=0.025
    integer,parameter :: time_interval=0.001
    integer, parameter :: output_interval=20
    !------------change allowed--------------!

end module parameters_for_simulation



program main 
    use test_module
    implicit none
    
end program