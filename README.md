# CFD_paricle_method_basics
粒子法を用いた簡単な水柱崩壊のプログラム

# このプログラムについて
このプログラムでは、粒子法を用いた水柱崩壊のシミュレーションを行う。以下の参考書の付録のc++で書かれたプログラムを、FORTRANで使えるようにすることを主な目的としている。また、この文章は自分の頭の中を整理するために書いていると言う節もある。

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

のようにして、一ステップ後（すなわち $\Delta t$ 後）の粒子の位置 $\boldsymbol{r}$ と粒子の速度 $\boldsymbol{u}$ を得ることができる。

## ラグランジュ微分の離散化

## 粘性項の離散化

## 圧力項の計算

# プログラムの概要
実際のプログラムで行われる処理（サブルーチン）を下に示す。

```fortran:
program main

|-water_tank_and_water_column_2d
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
### water_tank_and_water_column_2d
　初期の粒子配置を設定するルーチン。
```txt
|            |                                                     |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |oooooooooooo                                         |           |
|            |_____________________________________________________|           | 
|                                                                              | 
|                                                                              | 
|                                                                              |  
|______________________________________________________________________________| 
```
こんな感じで水槽と水柱を配置する。具体的には、必要になる粒子の数を```numberofparticles```として、i=1~numberofparticlesの各粒子に対して、粒子の位置を表す```particleposition(i)```や粒子の種類(壁など)を表す```particletype(i)```などを、メモリを確保しながら（allocateしながら）入力していく。

### mainLoopOfSimuation

　k〜k+1ステップに粒子の情報を更新するルーチン。以下のcalGravity~movePariceUsingPressureGradientで構成される。

　今、改めてナビエストークス方程式を見てみると、

$$\dfrac{D\boldsymbol{u}}{Dt}  = -\dfrac{1}{\rho}\nabla P +\nu\nabla^2\mathbf{u}+\boldsymbol{g}$$

のようになっており、流体の加速度$D\boldsymbol{u}/Dt$は$-1/\rho\cdot\nabla P$、$\nu\nabla^2\mathbf{u}$、$\boldsymbol{g}$の三つの項目で構成されていることがわかり、それぞれ圧力項、粘性項、重力項と呼ばれる。
　mainLoopOfSimuationにおいては、まず重力項と粘性項をcalGravity、calViscosityによって計算し、それら加速度を用いて一旦粒子の情報をアップデートする。次に、圧力項をcalPressureGradientによって計算し、それを用いて再度粒子の情報をアップデートする。また、calPressureGradientを計算するためにはまず、calPressureを用いて粒子の圧力を計算する。また、これらの計算は壁粒子に対しては行わないように気を付ける。



<ins>calGravity</ins>：重力項による粒子の加速を計算するルーチン

<ins>calViscosity</ins>：粘性項による粒子の加速を計算するルーチン

<ins>moveParticle</ins>：calGravity,calViscosityで計算した加速度をつかって粒子の位置・速度を更新するルーチン

<ins>collision</ins>：例外処理。異常接近した粒子があった場合に、粒子間距離を広げる。

<ins>calPressure</ins>：各粒子位置での圧力を計算するルーチン。以下のcalNumberDensity~setMinimumPressureで構成される。

　今、改めて圧力$P$は

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

## water_tank_and_column_2d
２次元の水柱に関しての初期状態を定めるルーチン。ルーチンの概要は以下のよう。

```fortran
subroutine water_tank_and_water_column_2d(x_watertank,y_watertank,x_watercolumn,y_watercolumn,wallthickness,dummywallthickness)

    nx_watertank = ...
    number_of_particles = ...
    allocate(particle_position(number_of_particles,3))
        .
        .
        .
    
    do iY = ...
        do iX =...
        !particle_position(i,:)=...
        end do 
    end do

end subroutine
```

まず、引数として
```fortran
x_watertank,y_watertank,x_watercolumn,y_watercolumn,wallthickness,dummywallthickness
```
の6つを受け入れる。```x_watertank,y_watertank```で水槽の（内壁の）大きさを定め、```x_watercolumn,y_watercolumn```で水柱の大きさを定め、```wallthickness,dummywallthickness```で各種壁の厚みを定める。いずれも受け入れ引数の単位は[m]。
　で、あらかじめ設定しておいた初期粒子間距離[m]を用いて、トータルで必要な粒子数```numberofparticles```を求め、```particle_position```等に必要粒子数分のメモリを確保していく。
 　その後、```iX,iY```によるループを用いて
$$
(\text{iX}\times\text{初期粒子間距離},\text{iY}\times\text{初期粒子間距離},)
$$
に位置する粒子についての```particletype(i)```などの情報を入れていく。このループは並列化して行う。


## writeData_inVtuFormat
VTKフォーマットは可視化のためのフォーマット。空行と空白を同様に扱う。VTKファイルは基本は以下の様な構造になっている.

```xml
<?xml version='1.0',encoding='UTF-8'?>
<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>
    <unstructuredGrid>
        <Piece NUmberOfcells='50000' NumberOfPoints='50000'>

            <Points>
                <Dataarray NumberOfComponents='2or3' type='Float32' Name='Position' format='ascii'>
                    <!座標データ、n行2/3列、nは粒子数、2/3は空間の次元>
                </Dataarray>
            </Points>

            <PointData>
                <Dataarray NumberOfComponents='i' type='Int32/Float32/etc' Name='hoge' format='ascii'>
                    <!粒子のhoge属性に関するデータ、n行i列、nは粒子数iはデータの次元>
                </Dataarray>
                <Dataarray NumberOfComponents='j' type='Int32/Float32/etc' Name='Foo' format='ascii'>
                    <!粒子のFoo属性に関するデータ、n行j列、nは粒子数iはデータの次元>
                </Dataarray>
                    .
                    .
                    .
            </PointData>

            <Cells>
                <Dataarray　type='Int32' Name='connectivity' format='ascii'>
                    <!後ほど解説>
                </Dataarray>
                <Dataarray　type='Int32' Name='offsets' format='ascii'>
                    <!後ほど解説>
                </Dataarray>
                <Dataarray type='Int32' Name='types 'format='ascii'>
                    <!後ほど解説>
                </Dataarray>
            </Cells>
            
            <Celldata>
                <Dataarray NumberOfComponents='i' type='Int32/Float32/etc' Name='hage' format='ascii'>
                    <!cellのhage属性に関するデータ、m行i列、mはcellの数、iはデータの次元>
                </Dataarray>
                <Dataarray NumberOfComponents='j' type='Int32/Float32/etc' Name='Fooo' format='ascii'>
                    <!cellのFooo属性に関するデータ、m行j列、mはcellの数、jはデータの次元>
                </Dataarray>
                    .
                    .
                    .
            </Celldata>
        </Piece>
    </unstructuredGrid>
</VTKFile>
```
主要な構造だけ抜き出すと以下のよう。[参考](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html?utm_source=chatgpt.com)
```xml
<VTKFile type="UnstructuredGrid" ... >
  <UnstructuredGrid>
    <Piece NumberOfPoints="#" NumberOfCells="#">
        <Points>...<!-座標->...</Points>
        <PointData>...<!-あ->...</PointData>
        <Cells>...</Cells>
        <CellData>...</CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>

```
<br>

### VTKFile 条件子
```xml
<VTKFile xmlns='VTK' byte_order='LittleEndian' version='1.0' type='UnstructuredGrid'>
```
で指定される、VTKファイルに関する諸々のパラメタ。

<ins>byte_order</ins>：エンディアン。（データを前から並べるか、後ろから並べるか、くらいの認識。）macはLittleEndian。

<ins>version</ins>：VTKフォーマットの仕様のバージョン。1.0でいいのかな。

<ins>type</ins>：stucturedは格子データのような規則正しく並んだデータ、unstructuredは粒子データのような不規則に並んだデータ。色々あるっぽいが、粒子法では```unstructuredGrid```を用いる。unstructuredGridは拡張子```.vtu```に対応しているので、このコードでは```.vtu```ファイルに出力している。
<br>

### piece条件子
```NumberOfPoints```は粒子数、```NumberOfCell```はセルの数。先に下のpointとcellに関する解説を読んだ方がわかりやすいと思う。
<br>

### points
```xml
<Points>
    <Dataarray NumberOfComponents='2or3' type='Float32' Name='Position' format='ascii'>
        0.0 0.0 0.0
        0.0 0.0 0.1
        0.0 0.0 0.2
            .
            .
            .
    </Dataarray>
</Points>
```
各粒子の座標情報が入る。```NumberOfComponents```は空間の次元。xmlは改行と空白を区別しないので、
```xml
0.0 0.0 0.0
0.0 0.0 0.1
0.0 0.0 0.2
    .
    .
    .
```
の部分は
```xml
0.0 0.0 0.0 0.0 0.0 0.1 ...
```
のように認識される。よって、```NumberOfComponents='3'```で、「この座標は三次元やで。やから３つごとにデータを読み取りなさいな」と言う指定を入れてやる必要がある。```type```はデータの型（INTなど）、```name```はまあ、自由に設定していいそのDataarrayの名前。```Format```はasciiかbinary。大規模計算とかだとbinaryの方がいいらしい。
<br>

### pointdata
各粒子にデータを載せる。圧力だったり温度だったり粒子の種類だったり。例えば温度という一次元データを各粒子に持たせることを考えると、
```xml
<Dataarray NumberOfComponents='1' type='Int32/Float32/etc' Name='temperature' format='ascii'>
    273
    280
    250
     .
     .
     .
</Dataarray>
```
ここでも、```NumberOfComponetnsは次元。paraviewなら三次元データ（速度場とか）も可視化できるそうな。すごいッピ！
### cells
ざっくばらんにいって仕舞えば、メッシュ。一つのメッシュのことをセルと言う。例えば下の i、j、k、l 番目の粒子は四角形のセルを作っている。
```
   i .__. k
     |  |
   j .__. l
```
そんな感じに、いろんな形のセルをいくらでも作ることができて、それを指定するのが```<cells>```の項目。（つまり、粒子集合の部分集合がセル）必ず、以下の三つをもつ。（それぞれ```Name```固定。）
```xml
<Cells>
    <Dataarray type='Int32' Name='connectivity' format='ascii'>
    <Dataarray type='Int32' Name='offsets' format='ascii'>
    <Dataarray type='Int32' Name='types 'format='ascii'>
</Cells>
```
<br>

<ins>**connectivity**</ins>：まずconnectivityだが、これはcellの頂点に関するデータ。
```xml
<Dataarray type='Int32' Name='connectivity' format='ascii'>
    0 3 5
    3 5 7 8
      .
      .
      .
</Datarray>
```
上のコードなら、1,4,6番目の粒子で三角形のセルを作り、4,6,8,9番目で４角形のセルを作り...という風。（vtkファイルの配列は0始まりなので、```<points>```で定めた粒子において、１番目の粒子のindexが0となる。~~0から始まるの直感的じゃなさすぎてほんと嫌い。Fortranを見習ってほしい。~~）ただ、何度も述べるように、xmlは改行を認識しないので実際には
```xml
<Dataarray type='Int32' Name='connectivity' format='ascii'>
    0 3 5 3 5 7 8...
</Datarray>
```
のように認識されている。よって、どこでcellを区切るのかを明示してやる必要がある。それをするのが以下の```offset```
<br>

<ins>**offset**</ins>
```xml
<Dataarray type='Int32' Name='offsets' format='ascii'>
    3
    7
    .
    .
    .
</Dataarray>
```
上のコードは、```connectivity```において「三番目までの粒子番号が最初のセル、（四番目から）七番目までの粒子番号が二つ目のセル、...」という指定を行う。
<br>

<ins>**types**</ins>：それぞれのcellが三角形なのか四角形なのか、みたいなことを明示してやる必要がある。それぞれの形状とtype番号との対応は[リンク](https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html)にある。今、
```cpp
  VTK_TRIANGLE = 5,
  VTK_QUAD = 9,
```
なので、```types```は以下のように書けば良いとわかる。
```xml
<Dataarray type='UInt8' Name='types 'format='ascii'>
    5
    9
    .
    .
    .
</Dataarray>
```
ちなみに、```types```の種類的に、```UINT8``` (符号なし8ビット)で十分だったりする。

### celldata
pointdataと同様、cellに温度だったり圧力だったりの情報を載せる。
```xml
<Dataarray NumberOfComponents='1' type='Int32/Float32/etc' Name='pressure' format='ascii'>
    10
    15
    .
    .
    .
</Dataarray>
```
cellの数だけデータ数（行数）が存在することになる。
<br>

### 粒子法におけるcellの指定
粒子法においては、cellは（多分）使う必要がない。なので、各点が一つのcellを作るとしてcellデータを記入すれば良い。この時、typesは```VTK_VERTEX=1```なので、1と記入すれば良い。
```hml
<cells>
    <Dataarray type='Int32' Name='connectivity' format='ascii'>
        0
        1
        .
        .
        .
        n-1
    </Dataarray>
    <Dataarray type='Int32' Name='offsets' format='asciss'>
        1
        2
        .
        .
        .
        n
    </Dataarray>
    <Dataarray type='UInt8' Name='types' format='ascii'>
        1
        1
        1
        .
        .
        .
    </Dataarray>        
</cells>
```
