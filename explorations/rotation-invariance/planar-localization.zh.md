# 平面点插入：用交换涂抹直接得到算子范数极限

2026-10-18（Asia/Shanghai）· ChatGPT 生成的局部推导，尚未经独立数学审阅

## 范围与结果

本笔记接续 [上一份探索](local-note.zh.md)，检查原稿 Appendix A.3–A.5 的点插入步骤。原稿为 Sussman–Alpoge 项目中的 [Rotation invariance from lattice symmetry and reflection positivity](https://github.com/ethanwsussman/Yang--Mills-LLM-experiment/blob/785a68fdeef78a1da854c4807fe1f17b597d05e7/QED3_slop/YM_Rotation_Invariance_simplified.pdf)，固定版本 `785a68fdeef78a1da854c4807fe1f17b597d05e7`，尤其公式 (50)、(53)、(55)–(60)。

得到的局部改进是：**在原稿的有界插入与协变恒等式成立时，平面正则化插入直接在算子范数中收敛，并有二次误差界。** 不必先取弱算子子列，再用 Montel 与唯一性识别。

这是一种可替换 A.5 中间论证的证明，没有发现足以推翻原稿的缺陷，也不声称具有文献新颖性。这里只建立局部解析模块；Section 3 的谱比较、全局旋转拼接、具体场论模型的假设与 YM4 / 质量间隙均未因此完成。

## 1. 输入与记号

沿用原稿的 Hilbert 空间 \(\mathcal H\)，以及

\[
q(\xi)=\operatorname{Re}\xi_t-|\operatorname{Im}\xi_x|,
\qquad W=\{q>0\},\qquad e=(1,0).
\]

需要以下输入，均局限在当前平面，横向变量仍保留 \(L^2\) 涂抹：

1. \(V(\xi)\) 在 \(W\) 内算子范数全纯，满足半群律，在 \(q\ge0\) 上为压缩算子。
2. 对支集位于 \(|t|\le\delta\) 的复 \(L^2\) 测试，\(B_\delta\) 线性且
   \[
   \|B_\delta(\chi)\|\le C_0\|\chi\|_2.
   \tag{1}
   \]
3. 对实 \(z=(t,x)\)、\(r\ge |t|\)，有
   \[
   V(re+z)B_\delta(\chi)V(re-z)
   =B_{\delta+r}(\tau_z\chi).
   \tag{2}
   \]

\(C_0\) 与 slab 宽度无关。原稿 A.3 的指数源估计与 A.4 的 Cauchy 估计提供 (1)；(2) 是其 (53)。这两点是后文的明确前提，而非由正则化极限倒推。

令 \(f\in L^2(\mathbb R^{d-2};\Lambda^2\mathbb C^d)\)，

\[
j_\varepsilon(z)=\frac{1_{|z|<\varepsilon}}{\pi\varepsilon^2},
\qquad
A_\varepsilon(\xi,\eta)
=V(\xi-\varepsilon e)B_\varepsilon(j_\varepsilon f)V(\eta-\varepsilon e).
\tag{3}
\]

省略下标 \(f\)。所有范数均为 \(\mathcal B(\mathcal H)\) 算子范数；设
\(g=\min(q(\xi),q(\eta))>0\)。

## 2. 先核实一致界：常数可明确保留

原稿 (55) 的证明只用到 (1)、(2) 与全纯性，未使用点插入极限。因此可以先取为输入，也可以如下复核。

令 \(r=(g-\varepsilon)/2\)，固定单位向量 \(u,v\)，定义

\[
h_\varepsilon(\zeta)
=\langle u,A_\varepsilon(\xi+\zeta,\eta-\zeta)v\rangle.
\]

它在 \(|\operatorname{Re}\zeta_t|<r\)、\(|\operatorname{Im}\zeta_x|<r\) 内全纯。写 \(\zeta=z+ib\)，当 \(|b_x|<r\) 时，可将外侧因子取为

\[
V(\xi-(r+\varepsilon)e+ib),\qquad
V(\eta-(r+\varepsilon)e-ib).
\]

它们均为压缩算子且不依赖实变量 \(z\)。对支于 \((-r,r)^2\) 的 \(\alpha\)，用 (2) 积分，再用 (1) 与 Young 不等式，得

\[
\left|\int\alpha(z)h_\varepsilon(z+ib)\,dz\right|
\le C_0\|f\|_2\|\alpha*j_\varepsilon\|_2
\le C_0\|f\|_2\|\alpha\|_2.
\]

故每个这样的虚平移上，实正方形内的 \(L^2\) 范数至多 \(C_0\|f\|_2\)。对 \(\mathbb C^2\) 中半径 \(r'<r\) 的复多圆盘使用均值不等式，并将虚部积分扩大到 \((-r,r)^2\)，得到

\[
|h_\varepsilon(0)|^2
\le \frac{4r^2}{\pi^2r'^4}C_0^2\|f\|_2^2.
\]

令 \(r'\uparrow r\)，再对 \(u,v\) 取上确界，即得

\[
\boxed{\ \|A_\varepsilon(\xi,\eta)\|
\le \frac{C_1\|f\|_2}{g-\varepsilon},\quad C_1=4C_0/\pi.\ }
\tag{4}
\]

常数是否最优并不重要；关键是损失一幂间距，且不依赖正则化半径。原稿以同一个 \(K\) 吸收这些数值常数，是可行的。

## 3. 两个圆盘的涂抹可以交换

定义实平移平均

\[
(\mathcal S_\sigma A_\varepsilon)(\xi,\eta)
=\int j_\sigma(z)A_\varepsilon(\xi+z,\eta-z)\,dz.
\]

当 \(\varepsilon+\sigma<\delta<g\) 时，由 (2)，

\[
\mathcal S_\sigma A_\varepsilon
=V(\xi-\delta e)
B_\delta((j_\sigma*j_\varepsilon)f)
V(\eta-\delta e).
\tag{5}
\]

逐点使用 (2) 合法，因为 \(|z_t|\le\sigma<\delta-\varepsilon\)。测试函数的时间支集也位于 \(|t|<\delta\)。积分是 Bochner 积分：平移在 \(L^2\) 连续，\(B_\delta\) 有界，外侧因子固定。

卷积交换立即给出

\[
\boxed{\mathcal S_\sigma A_\varepsilon
=\mathcal S_\varepsilon A_\sigma.}
\tag{6}
\]

这一步完全发生在正则化插入之间，尚未使用任何极限场。

## 4. 范数 Cauchy 估计与直接构造

设 \(0<\varepsilon,\sigma\le g/4\)。固定 \(|z|\le\sigma\)，记 \(L=|z_t|+|z_x|\)。当 \(L>0\) 时，考虑一复变量算子值函数

\[
G(w)=A_\varepsilon(\xi+wz,\eta-wz),
\qquad |w|\le R=\frac{g}{2L}.
\]

在该圆盘上，两侧的 \(q\) 均至少为 \(g/2\)；减去 \(\varepsilon\) 后仍至少为 \(g/4\)。由 (4)，

\[
\|G(w)\|\le M=4C_1\|f\|_2/g.
\]

令 \(s=R^{-1}=2L/g\)。由 \(L\le\sqrt2\sigma\) 可知 \(s^2\le1/2\)。Cauchy 系数估计和 \(w\mapsto-w\) 对称性给出

\[
\left\|\frac{G(1)+G(-1)}2-G(0)\right\|
\le M\frac{s^2}{1-s^2}
\le \frac{32C_1\|f\|_2}{g^3}L^2.
\]

\(z=0\) 的情形平凡。又因均匀圆盘满足

\[
\int j_\sigma(z)L^2\,dz
\le2\int j_\sigma(z)|z|^2\,dz=\sigma^2,
\]

且 \(j_\sigma(z)=j_\sigma(-z)\)，积分得到

\[
\|\mathcal S_\sigma A_\varepsilon-A_\varepsilon\|
\le32C_1\|f\|_2\sigma^2/g^3.
\tag{7}
\]

将 (6) 夹在两次这样的估计之间，得

\[
\boxed{\ \|A_\varepsilon-A_\sigma\|
\le\frac{32C_1\|f\|_2}{g^3}(\varepsilon^2+\sigma^2).\ }
\tag{8}
\]

因此 \(A_f(\xi,\eta)=\lim_{\varepsilon\downarrow0}A_\varepsilon(\xi,\eta)\) 直接存在于完备的算子范数空间，并且

\[
\|A_\varepsilon-A_f\|\le32C_1\|f\|_2\varepsilon^2/g^3,
\qquad
\|A_f\|\le C_1\|f\|_2/g.
\tag{9}
\]

收敛在 \(W\times W\) 的每个紧集上一致。每个 \(A_\varepsilon\) 本已联合范数全纯，所以极限联合范数全纯；不需要先经弱算子极限，也不需要另作半群夹乘来证明范数全纯性。线性、平移与伴随恒等式可直接传到范数极限。

## 5. 恢复涂抹与更细的展开

对紧支 \(\alpha\in L^2(\mathbb R^2)\)，其时间支集严格位于 \(|t|<\delta\)，由正则化版 (5) 及 \(\alpha*j_\varepsilon\to\alpha\) 的 \(L^2\) 收敛，得到原稿 (60)：

\[
B_\delta(\alpha f)
=\int\alpha(z)A_f(\delta e+z,\delta e-z)\,dz.
\tag{10}
\]

严格位于 slab 内提供一致的正间距；\(\alpha\) 紧支且 \(L^2\)，因而也 \(L^1\)，足以控制这个算子范数积分。这里 (10) 是范数极限构造的结果，不是 (8) 的前提。

取 \(\varepsilon<\delta<g\)，将 (10) 用于 \(j_\varepsilon\)，再夹乘两侧半群，便有精确恒等式

\[
A_\varepsilon(\xi,\eta)
=\int j_\varepsilon(z)A_f(\xi+z,\eta-z)\,dz.
\tag{11}
\]

令 \(h(\zeta)=A_f(\xi+\zeta,\eta-\zeta)\)。圆盘的一阶矩为零，二阶矩为

\[
\mathbb E[z_t^2]=\mathbb E[z_x^2]=\varepsilon^2/4,
\qquad \mathbb E[z_tz_x]=0.
\]

于是

\[
\boxed{\ A_\varepsilon
=A_f+\frac{\varepsilon^2}{8}
(\partial_{\zeta_t}^2+\partial_{\zeta_x}^2)h(0)+\mathcal R_\varepsilon.\ }
\tag{12}
\]

当 \(\varepsilon\le g/4\) 时，同样的一复变量偶次幂估计，现在使用 \(\|h(wz)\|\le2C_1\|f\|_2/g\)，给出

\[
\|A_\varepsilon-A_f\|\le16C_1\|f\|_2\varepsilon^2/g^3,
\qquad
\|\mathcal R_\varepsilon\|
\le\frac{256}{3}C_1\|f\|_2\varepsilon^4/g^5.
\tag{13}
\]

后一个常数来自 \(\mathbb E L^4\le4\mathbb E|z|^4=4\varepsilon^4/3\)。这些导数是有界算子值全纯函数的导数；无需把它们改写成可能涉及无界生成元定义域的形式对易子。

更一般的圆盘平均级数为

\[
\int j_\varepsilon(z)h(z)\,dz
=\sum_{m\ge0}\frac{\varepsilon^{2m}}{4^m m!(m+1)!}
(\partial_{\zeta_t}^2+\partial_{\zeta_x}^2)^m h(0),
\tag{14}
\]

在足够小半径内范数收敛。前三个系数是 \(1,1/8,1/192\)。这也给出一个可选的加速近似：

\[
\widehat A_\varepsilon=(4A_{\varepsilon/2}-A_\varepsilon)/3,
\qquad
\|\widehat A_\varepsilon-A_f\|
\le\frac{320}{9}C_1\|f\|_2\varepsilon^4/g^5.
\tag{15}
\]

它只是在固定正间距下消掉二阶偏差的线性组合。它带有负权重，不能据此声称仍有原正则化的正性性质，也不能据此降低全局频率增长指数。

## 6. 对后续构造的作用与边界

对有限个插入，若每个两侧间距至少为 \(r>0\)，且 \(\varepsilon\le r/4\)，用望远镜展开即可控制乘积：

\[
\left\|\prod_{j=1}^N A_{f_j,\varepsilon}
-\prod_{j=1}^N A_{f_j}\right\|
\le32N C_1^N(4/3)^{N-1}
\frac{\varepsilon^2}{r^{N+2}}\prod_{j=1}^N\|f_j\|_2.
\tag{16}
\]

这避免了“弱算子收敛能否穿过乘积”的疑问。它在 Appendix A.6/A.7 的紧支、严格分离配置上提供可用的加强；原稿已有其他识别积分的办法，故不把这项加强描述为修复已确认的错误。

但 \(r\downarrow0\) 时误差界发散。虚旋转下可用间距也可能随角参数衰减。因此本结果：

- 适用于固定间距或保持正下界的紧参数集；
- 不能直接用于同时取碰撞极限、无穷虚角极限；
- 不能代替原稿的闭合性无接触项论证；
- 不改变原稿每次插入的一幂间距损失，也未进一步改善上份笔记的频率阈值。

## 7. 检查与下一步

[check_planar_averaging.py](check_planar_averaging.py) 用精确有理数、极坐标偶次积分递推与多项式运算核对圆盘矩、(14) 的低阶系数、两次涂抹交换和 (15) 的二阶消去。有限维生成元例子只说明展开可有非零二阶项，不构造满足全部场论假设的模型。这些是公式检查，不能替代上述算子范数证明。

这次让我更在意的是：两个近似之间可能有一个精确的交换关系，比先构造极限再辨认它更有用。下一步值得检查的地方是正文 §4.2–§4.4：全局角度拼接、不同排序区域之间的一致性，以及继续后的核如何重新识别为原来的矩。这里仍没有做完整审计。
