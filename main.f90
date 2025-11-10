!------------change allowed--------------!
    !module parameters_for_simulationを適宜変更すること。
    !module inital_condition下に初期条件について定めるsubroutineを作る必要がある。
    !*all global parameters and variables must be defind in module:parameters_and_variables_for_simulation

!------------change allowed--------------!




module parameters_and_variables_for_simulation
    implicit none


    !gloval_parameters
        !dimention: シミュレーションの次元を選択する。2->二次元、３->三次元
        !pargicle_distance: 初期の粒子間距離[m]
        !time_interval: 計算ステップの時間幅[s]
        !output_interval: 何ステップごとに粒子の位置を出力するか
        !finish_time: 終了時刻[s]

        !------------change allowed--------------!
        integer,parameter :: dimention=3
        real(8),parameter :: paricle_distance=0.025
        real(8),parameter :: time_interval=0.001
        integer, parameter :: output_interval=20
        real(8),parameter :: finish_time = 2.0
        !------------change allowed--------------!
    !

    !global_variables
        !global変数：すべてのmoduleのsubroutineで参照可能な変数
        !allocatはsubroutineで行う。

        !particle_position: 粒子iのx,y,z座標をparticle_position(i,x=1/y=2/z=3)に収納する。
        !particle_velocity: 粒子iのx,y,z方向の速度ををparticle_velocity(i,x=1/y=2/z=3)に収納する。
        !accleration: 粒子iに関する加速度をacceleration(i,x=1/y=2/z=3)に収納する。
        !number_density: 重みつき粒子数密度Σw(ri-rj)
        !number_of_particles：粒子総数。粒子の初期配置を決定した段階で、particle_position等はこの数にallocateする。
        !particle_type：粒子の属性。液体粒子か、壁粒子か、ダミー壁粒子か
        !particle_prssure：各粒子位置での圧力
        !original_layer：各粒子が最初、水柱の中でどの高さに属しているか。

        real(8),allocatable :: particle_position(:,:)
        real(8),allocatable :: particle_velocity(:,:)
        real(8),allocatable :: acceleration(:,:)
        real(8),allocatable :: number_density(:)
        real(8),allocatable :: particle_type(:)
        real(8),allocatable :: particle_pressure(:)
        real(8),allocatable :: Original_layer(:)
        integer :: number_of_particles

    !

end module parameters_and_variables_for_simulation

module initial_particle_position_velocity_particle_type
    use parameters_and_variables_for_simulation
    implicit none

    !以下のsubroutineで初期の粒子の配置、壁の配置、ダミー壁の配置、粒子速度を定める。
    !壁、ダミー壁はすべて粒子として扱う。
    !粒子についての位置particle_positions(dim,n)を定めた上で、
    !その粒子のtype:particle_type(n)が液体粒子なのか壁なのかダミー壁なのかを(後から)定めることで、
    !液体粒子、壁、ダミー壁の位置を定めることができる。

    contains
    subroutine water_tank_and_water_column_2d(x_watertank,y_watertank,x_column,y_column)

        !二次元水柱崩壊を計算する場合の初期条件
        !水槽の大きさはx_watertank,y_watertank、水柱の大きさはx_column,y_column
        implicit none
        integer,intent(in)::x_watertank,y_watertank,x_column,y_column

        !count_num_of_particles
        !暫定的に
        number_of_particles = 5

        !メモリの確保
        allocate(particle_position(number_of_particles,3))
        allocate(particle_velocity(number_of_particles,3))
        allocate(acceleration(number_of_particles,3))
        allocate(number_density(number_of_particles))
        allocate(particle_type(number_of_particles))
        allocate(particle_pressure(number_of_particles))
        allocate(Original_layer(number_of_particles))

        !
    
    end subroutine

end module initial_particle_position_velocity_particle_type


module output_module
    use open
    use parameters_and_variables_for_simulation
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
        open(10,file=file_name,action='write',status='replace')
        !地の文
        write(10,'(A)') "<?xml version='1.0' encoding='UTF-8'?>"
        write(10,'(A)')"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>"
        write(10,'(A)') "   <UnstructuredGrid>"
        write(temp_char,'(A,I0,A,I0,A)')"      <Piece NumberOfcells='",number_of_particles,"' NumberOfPoints='",number_of_particles,"'>"
        write(10,'(A)') temp_char
        !粒子の座標出力
        write(10,*) 
        write(10,'(A)') "         <Points>"
        write(10,'(A)') "            <DataArray NumberOfComponents='3' type='Float32' Name='Particle_Position' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)particle_position(i,:)
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Points>"
        write(10,*)
        write(10,'(A)') "         <PointData>"
        !particle_typeの出力
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='INT32' Name='Particle_type' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) particle_type(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !絶対速度の出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='abs_velocity' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) sqrt(particle_velocity(i,1)**2 + particle_velocity(i,2)**2 + particle_velocity(i,3)**2)
        end do
        write(10,'(A)') "            </DataArray>"
        !圧力の出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) particle_pressure(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !初期の液体の高さの出力
        write(10,*)
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Float32' Name='Original_layer' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*) Original_layer(i)
        end do
        write(10,'(A)') "            </DataArray>"
        !その他記述(cellなど)
        write(10,'(A)') "         </PointData>"
        write(10,*)
        write(10,'(A)') "         <Cells>"
        write(10,'(A)') "            <DataArray type='INT32' Name='connectivity' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i-1
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='INT32' Name='offsets' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='INT32' Name='types' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)1
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Cells>"
        write(10,*)
        write(10,'(A)') "      </Piece>"
        write(10,'(A)') "   </UnstructuredGrid>"
        write(10,'(A)') "</VTKFile>"
        close(10)
    end subroutine writedatainvtuformat
end module output_module

module maincalculation
    use parameters_and_variables_for_simulation
    use output_module
    implicit none

    contains

    subroutine mainloopofsimulation()
        implicit none
        integer outputstep


        !outputstep = 0
        !do while not end sim, call writedatainvtu(outputstep),
        !then, outputstep+=1

        do  outputstep= 1,5
            call writedatainvtuformat(outputstep)
        end do

    end subroutine
    

end module maincalculation




program main 
    !---------modules----------!
    use omp_lib
    use initial_particle_position_velocity_particle_type
    use maincalculation
    !---------modules----------!
    implicit none


    !-------calling subroutines-------!
    call water_tank_and_water_column_2d(15,6,5,5)
    call mainloopofsimulation
    !-------calling subroutines-------!




    write(*,*) particle_position(1,1)

end program