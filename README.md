# CFD_paricle_method_basics
粒子法を用いた簡単な水柱崩壊のプログラム

## このプログラムについて
このプログラムでは、粒子法を用いた水柱崩壊のシミュレーションを行う。以下の参考文献の付録のc++で書かれたプログラムを、FORTRANに書き直したものである。

### 参考文献
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
$$ \boldsymbol{u}$$
を得ることができる。

### プログラムの概要