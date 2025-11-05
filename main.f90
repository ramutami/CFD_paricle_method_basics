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
    !particle_limit: 粒子数の上限を与える
    !finish_time: 終了時刻[s]

    !------------change allowed--------------!
    integer,parameter :: dimention=3
    real(8),parameter :: paricle_distance=0.025
    real(8),parameter :: time_interval=0.001
    integer, parameter :: output_interval=20
    integer,parameter :: particle_limit=5000
    real(8),parameter :: finish_time
    
    !------------change allowed--------------!

end module parameters_for_simulation





program main 
    use test_module
    implicit none
    
end program