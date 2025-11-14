!------------change allowed--------------!
    !module parameters_for_simulationを適宜変更すること。
    !module inital_condition下に初期条件について定めるsubroutineを作る必要がある。
    !*all global parameters and variables must be defind in module:parameters_and_variables_for_simulation
!------------change allowed--------------!

!------------Caution!--------------!
    !i,j,j…などはsubroutine内部でのみ使うようにすること。
!------------Caution!--------------!






module parameters_and_variables_for_simulation
    implicit none


    !gloval_parameters
        !dimention: シミュレーションの次元を選択する。2->二次元、３->三次元
        !pargicle_distance: 初期の粒子間距離[m]
        !time_interval: 計算ステップの時間幅[s]
        !output_interval: 何ステップごとに粒子の位置を出力するか
        !finish_time: 終了時刻[s]

        !------------change allowed--------------!
        real(8),parameter :: particle_distance=0.025
        real(8),parameter :: time_interval=0.001
        integer, parameter :: output_interval=20
        real(8),parameter :: finish_time = 2.0
        !------------change allowed--------------!


        integer,parameter :: fluid = 1
        integer,parameter :: wall = 2
        integer,parameter :: dummywall = 3
        real(8),parameter :: g_x = 0.0
        real(8),parameter :: g_y = -9.80665
        real(8),parameter :: g_z = 0.0

    !

    !global_variables
        !global変数：すべてのmoduleのsubroutineで参照可能な変数
        !allocatはsubroutineで行う。

        !particle_position: 粒子iのx,y,z座標をparticle_position(i,x=1/y=2/z=3)に収納する。
        !particle_velocity: 粒子iのx,y,z方向の速度ををparticle_velocity(i,x=1/y=2/z=3)に収納する。
        !particle_acceleration: 粒子iに関する加速度をacceleration(i,x=1/y=2/z=3)に収納する。
        !number_density: 重みつき粒子数密度Σw(ri-rj)
        !number_of_particles：粒子総数。粒子の初期配置を決定した段階で、particle_position等はこの数にallocateする。
        !particle_type：粒子の属性。液体粒子か、壁粒子か、ダミー壁粒子か
        !particle_prssure：各粒子位置での圧力
        !original_layer：各粒子が最初、水柱の中でどの高さに属しているか。

        real(8),allocatable :: particle_position(:,:)
        real(8),allocatable :: particle_velocity(:,:)
        real(8),allocatable :: particle_acceleration(:,:)
        real(8),allocatable :: number_density(:)
        integer,allocatable :: particle_type(:)
        real(8),allocatable :: particle_pressure(:)
        real(8),allocatable :: Original_layer(:)
        integer :: number_of_particles
        real(8) :: n0_for_laplacian

    !

end module parameters_and_variables_for_simulation

module initial_particle_position_velocity_particle_type
    use omp_lib
    use parameters_and_variables_for_simulation
    implicit none

    !以下のsubroutineで初期の粒子の配置、壁の配置、ダミー壁の配置、粒子速度を定める。
    !壁、ダミー壁はすべて粒子として扱う。
    !粒子についての位置particle_positions(dim,n)を定めた上で、
    !その粒子のtype:particle_type(n)が液体粒子なのか壁なのかダミー壁なのかを(後から)定めることで、
    !液体粒子、壁、ダミー壁の位置を定めることができる。

    contains
    subroutine water_tank_and_water_column_2d(x_watertank,y_watertank,x_watercolumn,y_watercolumn,wallthickness,dummywallthickness)

        !二次元水柱崩壊を計算する場合の初期条件
        !水槽の大きさは,壁の内側では買った大きさがx_watertank*y_watertank、水柱の大きさはx_watercolumn*y_watercolumn [m]

        use parameters_and_variables_for_simulation
        implicit none

        real(8),intent(in)::x_watertank,y_watertank,x_watercolumn,y_watercolumn
        integer,intent(in)::wallthickness,dummywallthickness
        !水槽や水柱のx,y方向の粒子数、例えばnx_watertank = x_watertank/particle_distance
        integer :: nx_watertank,ny_watertank,nx_watercolumn,ny_watercolumn
        real(8) :: epsilon = 0.01
        !計算ようの変数
        integer :: i,ix,iy
        integer :: temp_int_1,temp_int_2,temp_int_3

        !nx,nyの計算。epsilon*particl_distanceをつけるのは、real(8)の不安定さ故。
        nx_watertank = floor((x_watertank+particle_distance*epsilon)/particle_distance)
        ny_watertank = floor((y_watertank+particle_distance*epsilon)/particle_distance)
        nx_watercolumn = floor((x_watercolumn+particle_distance*epsilon)/particle_distance)
        ny_watercolumn = floor((y_watercolumn+particle_distance*epsilon)/particle_distance)
        write(*,*) "water tank size :", nx_watertank,ny_watertank
        write(*,*) "water column size : ",nx_watercolumn,ny_watercolumn

        !num_of_particle
        !必要な粒子数のメモリの数：
        !壁については　(nxtank+2*(wallthickness+dummywallthickness))*(nytank+wallthicnkness+dummywallthickness)-nx_tank*ny_tank
        !水柱については　nx_watercolumn*ny_watercolumn
        number_of_particles = (nx_watertank+2*(wallthickness+dummywallthickness))*(ny_watertank+wallthickness+dummywallthickness)-&
        & nx_watertank*ny_watertank
        number_of_particles = number_of_particles + nx_watercolumn*ny_watercolumn
        write(*,*) 'number of particles : ', number_of_particles

        !メモリの確保
        !二次元三次元によらず、三次元の座標を確保する。（計算の内容は変わらないため。）
        allocate(particle_position(number_of_particles,3))
        allocate(particle_velocity(number_of_particles,3))
        allocate(particle_acceleration(number_of_particles,3))
        allocate(number_density(number_of_particles))
        allocate(particle_type(number_of_particles))
        allocate(particle_pressure(number_of_particles))
        allocate(Original_layer(number_of_particles))

        !初期粒子位置の設定        
        !y座標を固定し、xについてのloopを回しながら初期条件を設定してゆく。

        !以下はdummywall＋wallの厚さが4の場合

        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |                                                     |           |
        !|            |○____________________________________________________|           | iY=+1
        !|                                                                              | iY=-0
        !|                                                                              | iY=-1
        !|                                                                              | iY=-2 
        !|______________________________________________________________________________| iY=-3
        ! -3,-2,-1, 0,+1 =iX

        !○の位置が(ix,iy) = 1,1の座標となる。
        !粒子番号iについては、上の図で左下を1として、
        !xを走査しながらyを上げていってiを順次定めていく。（water_tankに対し）
        !また、water_tankの最後の粒子をi=i_last_tankとした時、
        !水柱の左下に来る粒子をi=i_last_tank+1として、同じようにxを走査しながらyを上げていってiを順次定めていく。

        !計算量削減のため、毎回計算する変数はループの外側でtemp_intで計算してしまう。

        !床dummy
        temp_int_1 = 1-wallthickness-dummywallthickness
        temp_int_2 = nx_watertank+2*(wallthickness+dummywallthickness)
        temp_int_3 = 1-wallthickness-dummywallthickness
        !$omp parallel do private(i,iX)
        do iY = 1-wallthickness-dummywallthickness,1-wallthickness-1
            do iX = 1-wallthickness-dummywallthickness,nx_watertank+wallthickness+dummywallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do 
        end do
        !end $omp parallel do

        !床wall/dummy
        !$omp parallel do private(i,iX)
        do iY = 1-wallthickness,1-1
            do iX = 1-wallthickness-dummywallthickness,1-wallthickness-1    
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = 1-wallthickness,nx_watertank+wallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+wallthickness+1,nx_watertank+wallthickness+dummywallthickness
                i = (iY-temp_int_1)*temp_int_2+(iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do

        !壁
        temp_int_1 = (wallthickness+dummywallthickness)*(nx_watertank+2*(wallthickness+dummywallthickness))
        temp_int_2 = 2*(wallthickness+dummywallthickness)
        temp_int_3 = 1-wallthickness-dummywallthickness
        !$omp parallel do private(i,iX)
        do iY = 1,ny_watertank-wallthickness
            do iX = 1-wallthickness-dummywallthickness,1-wallthickness-1
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = 1-wallthickness,1-1
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+1,nx_watertank+wallthickness
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = wall
                Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+wallthickness+1,nx_watertank+wallthickness+dummywallthickness
                i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
                particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
                particle_velocity(i,:) = real(0.0,kind=8)
                particle_acceleration(i,:) = real(0.0,kind=8)
                particle_type(i) = dummywall
                Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do
        !$omp parallel do private(i,iX)
        do iY = ny_watertank-wallthickness+1,ny_watertank
            do iX = 1-wallthickness-dummywallthickness,1-1
            i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = wall
            Original_layer(i) = iY*particle_distance
            end do
            do iX = nx_watertank+1,nx_watertank+wallthickness+dummywallthickness
            i = temp_int_1 + (iY-1)*temp_int_2 + (iX-temp_int_3+1)-nx_watertank
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = wall
            Original_layer(i) = iY*particle_distance
            end do
        end do
        !end $omp parallel do
        
        temp_int_1 = (nx_watertank+2*(wallthickness+dummywallthickness))*(ny_watertank+wallthickness+dummywallthickness)-&
        & nx_watertank*ny_watertank
        do iY = 1,ny_watercolumn   
            do iX = 1,nx_watercolumn
            i = temp_int_1 + (iY-1)*nx_watercolumn+iX
            particle_position(i,:) = (/iX*particle_distance,iY*particle_distance,real(0.0,kind=8)/)
            particle_velocity(i,:) = real(0.0,kind=8)
            particle_acceleration(i,:) = real(0.0,kind=8)
            particle_type(i) = fluid
            Original_layer(i) = iY*particle_distance
            end do
        end do

    end subroutine

end module initial_particle_position_velocity_particle_type

module function_module  
    implicit none
    contains 

    real(8) function weight_function(distance,Re)
        implicit none
        real(8),intent(in)::distance,Re
        if (distance<Re) then
            weight_function =(Re/distance)-1.0
        else if (distance>=Re) then
            weight_function = 0.0
        end if
        return 
    end function weight_function

end module function_module

module calculation_of_parameters
end module calculation_of_parameters

module calculation_module
    use parameters_and_variables_for_simulation
    implicit none
    contains

    subroutine calgravity()
        !重力項による粒子にかかる加速度
        implicit none
        !内部変数
        integer :: i

        !$omp parallel do
        do i = 1,number_of_particles
            if (particle_type(i)== fluid) then
                particle_acceleration(i,1) = g_x
                particle_acceleration(i,2) = g_y
                particle_acceleration(i,3) = g_z
            end if
        end do
        !end $omp parallel do

    end subroutine calgravity

    subroutine calviscosity()
        implicit none


    end subroutine calviscosity

    subroutine moveparticle()
        !加速度を受けて粒子の位置を更新する
        implicit none
        !内部変数
        integer :: i

        !$omp parallel do
        do i= 1,number_of_particles
            if (particle_type(i) == fluid) then
                !速度の更新
                particle_velocity(i,1) = particle_velocity(i,1)+particle_acceleration(i,1)*time_interval
                particle_velocity(i,2) = particle_velocity(i,2)+particle_acceleration(i,2)*time_interval
                particle_velocity(i,3) = particle_velocity(i,3)+particle_acceleration(i,3)*time_interval
                !位置の更新
                particle_position(i,1) = particle_position(i,1)+particle_velocity(i,1)*time_interval
                particle_position(i,2) = particle_position(i,2)+particle_velocity(i,2)*time_interval
                particle_position(i,3) = particle_position(i,3)+particle_velocity(i,3)*time_interval
            end if
        end do
        !end $omp parallel do

    end subroutine moveparticle

end module calculation_module

module output_module
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
        write(temp_char,'(A,I0,A,I0,A)')"      <Piece NumberOfCells='",number_of_particles,"' NumberOfPoints='",number_of_particles,"'>"
        write(10,'(A)') temp_char
        !粒子の座標出力
        write(10,*) 
        write(10,'(A)') "         <Points>"
        write(10,'(A,I0,A)') "            <DataArray NumberOfComponents='3' type='Float32' Name='Particle_Position' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)particle_position(i,:)
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,'(A)') "         </Points>"
        write(10,*)
        write(10,'(A)') "         <PointData>"
        !particle_typeの出力
        write(10,'(A)') "            <DataArray NumberOfComponents='1' type='Int32' Name='Particle_type' format='ascii'>"
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
        write(10,'(A)') "            <DataArray type='Int32' Name='connectivity' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i-1
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='Int32' Name='offsets' format='ascii'>"
        do i = 1,number_of_particles
            write(10,'(A)',advance='no')"            " 
            write(10,*)i
        end do
        write(10,'(A)') "            </DataArray>"
        write(10,*)
        write(10,'(A)') "            <DataArray type='Int32' Name='types' format='ascii'>"
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

module mainloop
    use omp_lib
    use parameters_and_variables_for_simulation
    use output_module
    implicit none

    contains

    subroutine mainloopofsimulation()
        use calculation_module
        implicit none
        integer timestep,outputstep


        !outputstep = 0
        !do while not end sim, call writedatainvtu(outputstep),
        !then, outputstep+=1

        outputstep=0
        mainloop : do  timestep= 0,20000
            !初回の処理
            if(timestep == 0)then
                call writedatainvtuformat(0)
                cycle mainloop
            end if
            call calgravity()
            call moveparticle()
            if (mod(timestep,20)==0) then
                outputstep=outputstep+1
                call writedatainvtuformat(outputstep)
                write(*,*) time_interval*outputstep
            end if
        end do mainloop 

    end subroutine
    

end module mainloop


program main 
    !---------modules----------!
    use initial_particle_position_velocity_particle_type
    use mainloop
    !---------modules----------!
    implicit none


    !-------calling subroutines-------!
    call water_tank_and_water_column_2d(real(2.0,8),real(1.2,8),real(1.0,8),real(1.2 ,8),3,2)
    call mainloopofsimulation
    !-------calling subroutines-------!


end program