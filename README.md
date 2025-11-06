# CFD_paricle_method_basics
粒子法を用いた簡単な水柱崩壊のプログラム

## このプログラムについて
このプログラムでは、粒子法を用いた水柱崩壊のシミュレーションを行う。以下の参考書の付録のc++で書かれたプログラムを、FORTRANで使えるようにすることを主な目的としている。。

### 参考書
粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜,丸善出版

## 粒子法の概要
微分方程式を計算機で解く場合と基本的なアイデアは同じである。すなわち、微分方程式を時間ステップについて離散化して、各ステップごとに粒子の位置を更新＋出力していく。

### Navier-stokes equation
ラグランジュ形式で書いたNavier-stokes equationは以下のよう

$$\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}$$

この式を時間について離散化すれば、

$$\Delta \boldsymbol{u} = f\cdot\Delta{t}$$
$$\Delta \boldsymbol{r} = \Delta\boldsymbol{u}\cdot\Delta{t}$$

のようにして、一ステップ後（すなわち
$$\Delta t$$
後）の粒子の位置
$$\boldsymbol{r}$$
と粒子の速度
$$\boldsymbol{u}$$
を得ることができる。

### プログラムの概要



```fortran:
program main
|-initial_particle_position_velocity_particle_type
|-calConstantParameter
|-mainLoopOfSimuation   
  |-calGravity          
  |-calViscosity        
  |-moveParticle        
  |-collision
  |-calPressure
    |-calNumberDensity
    |-setBoundaryCondition
    |-setSourceTerm
    |-setMatrix
    |-solveSimultaniousEquationByGaussEliminationMethod
    |-removeNegativePressure
    |-setMinimumPressure
  |-calPressureGradient
  |-moveParticleUsingPressureGradient
  |-if (timestep = outputstep) {-writeData_inVtuFormat}
```

### それぞれのsubroutineの説明

<ins>mainLoopOfSimuation</ins><br>
k~k+1ステップに粒子の情報を更新するルーチン

<ins>calGravity</ins><br>
重力による粒子の加速を計算するルーチン

<ins>calViscosity</ins><br>
粘性項による粒子の加速を計算するルーチン

<ins>moveParticle</ins><br>
CalGravity,calViscosityで計算した加速度をつかって粒子の位置・速度を更新するルーチン





    



## 参考文献
[Lee et al] (https://www.sciencedirect.com/science/article/pii/S0045782519305067)