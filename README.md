# CFD_paricle_method_basics
粒子法を用いた簡単な水柱崩壊のプログラム

# このプログラムについて
このプログラムでは、粒子法を用いた水柱崩壊のシミュレーションを行う。以下の参考書の付録のc++で書かれたプログラムを、FORTRANで使えるようにすることを主な目的としている。

## 参考書
粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜,丸善出版

# 粒子法の概要
　以下は、「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」の内容の要点を抜き出したものである。
　
　基本的には、微分方程式を計算機で解く場合となアイデアは同じである。すなわち、微分方程式を時間ステップについて離散化して、各ステップごとに粒子の位置を更新＋出力していく。

## Navier-stokes equation
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

## ラグランジュ微分の離散化

## 粘性項の離散化

## 圧力項の計算

# プログラムの概要
実際のプログラムで行われる処理（サブルーチン）を下に示す。

```fortran:
program main

|-initial_particle_position_velocity_particle_type
|-calConstantParameter
|-mainLoopOfSimulation   
  |-calGravity          
  |-calViscosity        
  |-moveParticle        
  |-collision
  |-calPressure
    |-calNumberDensity
    |-setBoundaryCondition
    |-setSourceTerm
    |-setMatrix
      |-exceptionalProcessingForBoundaryCondition
        |-checkBoundaryCondition
        |-increaseDiagonalTerm
    |-solveSimultaniousEquationByGaussEliminationMethod
    |-removeNegativePressure
    |-setMinimumPressure
  |-calPressureGradient
  |-moveParticleUsingPressureGradient
  |-if (timestep = outputstep) {writeData_inVtuFormat}
  |-if (time>finish time){exit mainloopOfSimulation}

end program
```

## それぞれのsubroutineの説明

### mainLoopOfSimuation

　k~k+1ステップに粒子の情報を更新するルーチン。以下のcalGravity~movePariceUsingPressureGradientで構成される。

　今、改めてナビエストークス方程式を見てみると、

$$\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}$$

のようになっており、流体の加速度
$$D\boldsymbol{u}/Dt$$は
$$-1/\rho\cdot\nabla P$$
、
$$\nu\nabla^2\mathbf{u}$$、
$$\boldsymbol{g}$$
の三つの項目で構成されていることがわかり、それぞれ圧力項、粘性項、重力項と呼ばれる。
　mainLoopOfSimuationにおいては、まず重力項と粘性項をcalGravity、calViscosityによって計算し、それら加速度を用いて一旦粒子の情報をアップデートする。次に、圧力項をcalPressureGradientによって計算し、それを用いて再度粒子の情報をアップデートする。また、calPressureGradientを計算するためにはまず、calPressureを用いて粒子の圧力を計算する。また、これらの計算は壁粒子に対しては行わないように気を付ける。



<ins>calGravity</ins>：重力項による粒子の加速を計算するルーチン

<ins>calViscosity</ins>：粘性項による粒子の加速を計算するルーチン

<ins>moveParticle</ins>：calGravity,calViscosityで計算した加速度をつかって粒子の位置・速度を更新するルーチン

<ins>collision</ins>：例外処理。異常接近した粒子があった場合に、粒子間距離を広げる。

<ins>calPressure</ins>：各粒子位置での圧力を計算するルーチン。以下のcalNumberDensity~setMinimumPressureで構成される。

　今、改めて圧力は

$$-\dfrac{1}{\rho_0}\dfrac{2d}{\lambda^0 n^0}\displaystyle\sum_{j\neq i}\left(P_j^{k+1}-P_i^{k+1}\right)w(|\boldsymbol{r}_j^{\*}-\boldsymbol{r}_i^{\*}|)=\dfrac{1}{\Delta t^2}\dfrac{n_i^{\*}-n^0}{n^0}$$

によって計算されるのであった。（詳しくは「<ins>粒子法入門〜流体シミュレーションの基礎から並列計算と可視化まで〜</ins>」参照。）よって、圧力を計算するために必要なのは、
$$
n_i^{*}
$$
、

- <ins>calnumberDensity</ins>：粒子数密度を計算するルーチン

- <ins>setBoundaryCondition</ins>：


<ins>calPressureGradient</ins>：圧力項による粒子の加速を計算するルーチン

<ins>moveParticleusingPressureGradient</ins>：calPressureGradientで計算した加速度を使って粒子の位置・速度を更新するルーチン。ついでに、加速度をゼロにリセットする。


# 関数

## weight関数

重み関数。

# 各ルーチンの詳細について

## writeData_inVtuFormat
VTKフォーマットは、paraviewを用いた可視化のためのフォーマットで、以下の様な構造になっている.

```xml
<?xml version='1.0',encoding='UTF-8'?>
<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0,1' type='UnstructuredGrid'>
    <unstructuredGrid>
    
    </unstructuredGrid>
</VTKFile>
```
構造だけ抜き出すと以下のよう。[参考](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html?utm_source=chatgpt.com)
```xml
<VTKFile type="UnstructuredGrid" ...>
  <UnstructuredGrid>
    <Piece NumberOfPoints="#" NumberOfCells="#">
        <PointData>...<!-あ->...</PointData>
        <CellData>...</CellData>
        <Points>...<!-座標->...</Points>
        <Cells>...</Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>

```

### VTKFile 条件子

<ins>byte_order</ins>：エンディアン。（データを前から並べるか、後ろから並べるか、くらいの認識。）macはLittleEndian。

<ins>version</ins>：VTKフォーマットの仕様のバージョン。1.0でいいのかな。

<ins>type</ins>：stucturedは格子データのような規則正しく並んだデータ、unstructuredは粒子データのような不規則に並んだデータ。色々あるっぽいが、粒子法では```unstructuredGrid```を用いる。unstructuredGridは拡張子```.vtu```に対応しているので、このコードでは```.vtu```ファイルに出力している。

### points
各粒子の座標情報が入る。例えば三次元なら、以下のよう。
```datastudio
0.0 0.0 0.0
0.0 0.0 0.1
0.0 0.0 0.2
.
.
.
```

### 




