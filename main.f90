!------------change allowed--------------!
    !module parameters_for_simulationを適宜変更すること。
    !module inital_condition下に初期条件について定めるsubroutineを作る必要がある。
    
!------------change allowed--------------!





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
    integer,parameter :: particle_limit=1
    real(8),parameter :: finish_time = 2.0
    !------------change allowed--------------!

end module parameters_for_simulation

module global_variables
    use parameters_for_simulation
    implicit none

    !global変数：すべてのmoduleのsubroutineで参照可能な変数
    !allocatはsubroutineで行う。

    !particle_position: 粒子iのx,y,z座標をparticle_position(i,x=1/y=2/z=3)に収納する。
    !particle_velocity: 粒子iのx,y,z方向の速度ををparticle_velocity(i,x=1/y=2/z=3)に収納する。
    !accleration: 粒子iに関する加速度をacceleration(i,x=1/y=2/z=3)に収納する。
    !number_density: 重みつき粒子数密度Σw(ri-rj)

    real(8),allocatable :: particle_position(:,:)
    real(8),allocatable :: particle_velocity(:,:)
    real(8),allocatable :: acceleration(:,:)
    real(8),allocatable :: number_density(:)

    contains
    subroutine allocate_global_variables
        implicit none
        allocate(particle_position(particle_limit,3))
        allocate(particle_velocity(particle_limit,3))
        allocate(acceleration(particle_limit,3))
        allocate(number_density(particle_limit))
    end subroutine

end module global_variables

module initial_particle_position_velocity_particle_type
    use global_variables
    implicit none

    !以下のsubroutineで初期の粒子の配置、壁の配置、ダミー壁の配置、粒子速度を定める。
    !壁、ダミー壁はすべて粒子として扱う。
    !粒子についての位置particle_positions(dim,n)を定めた上で、
    !その粒子のtype:particle_type(n)が液体粒子なのか壁なのかダミー壁なのかを(後から)定めることで、
    !液体粒子、壁、ダミー壁の位置を定めることができる。
    !粒子速度については、particle_velocity(dim,n)を定める。

    contains
    subroutine water_tank_and_water_column()
        implicit none

        particle_position(1,1) = 1.0
    
    end subroutine

end module initial_particle_position_velocity_particle_type


module output_module
    use global_variables
    implicit none

    contains
    subroutine writedatainvtuformat(file_number)
        !VTKファイルに粒子についての情報を出力するサブルーチン。
        !file_numberは（連番）VTKファイルの番号。
        implicit none
        character(128) :: file_name
        character(256) :: temp_char
        integer,intent(in) :: file_number
        integer :: i

        !./output_vtuファイル下の、output_(file_number).vtuに出力する。
        !現在、file_numberは６桁まで対応している。
        write(file_name,'(A,I6.6,A)') "./output_vtu/output_",file_number,".vtu"
        
        !出力部分
        open(10,file=file_name,action='write',status='replace',position='append')

        !地の文
        write(10,'(A)') "<?xml version='1.0',encoding='UTF-8'?>"
        write(10,'(A)')"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
        write(10,'(A)') "   <UnstructuredGrid>"
        write(temp_char,'(A,I0,A,I0,A)')"      <Piece NUmberOfcells='",particle_limit,"' NumberOfPoints='",particle_limit,"'>"
        write(10,'(A)') temp_char
        write(10,*) 
        write(10,'(A)') "         <Points>"
        write(10,'(A)') "            <DataArray NumberOfComponents='3' type='Float32' Name='Particle_Position' format='ascii'>"
        
        !粒子の座標出力
        do i = 1,particle_limit
            write(10,'(A)',advance='no')"            " 
            write(10,*)particle_position(i,:)
        end do

        !地の文
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Points>"
        write(10,*)
        write(10,'(A)') "      </Piece>"
        write(10,'(A)') "   </UnstructuredGrid>"
        write(10,'(A)') "</VTKFile>"



        close(10)


    end subroutine writedatainvtuformat
end module output_module





program main 
    !---------modules----------!
    use test_module
    use global_variables
    use output_module
    !---------modules----------!

    use initial_particle_position_velocity_particle_type
    implicit none


    !-------calling subroutines-------!
    call allocate_global_variables()
    call water_tank_and_water_column()
    call writedatainvtuformat(2)
    !-------calling subroutines-------!




    write(*,*) particle_position(1,1)

end program