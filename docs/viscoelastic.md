# 黏弹性开发第一阶段：材料与独立干摩擦

这是独立干摩擦任务，不是完整黏弹混合润滑。稳定分支保持rc10。

## 两个Persson输入文件

主曲线：首行可为材料名，后面三列 f[Hz], storage[MPa], loss[MPa]。
温移表：首行可为材料名，后面三列 T[C], aT, bT（线性值，均为正）。
空行及#或!开头注释可用；频率/温度严格递增，不静默排序、平滑或补尾。
默认标准主曲线表示Young模量；shear或reduced必须显式选择，分别换算
E=2(1+nu)G及E=(1-nu²)E_reduced，然后内部只执行一次接触约化。

E(omega,T)=bT E(aT omega,Tref)。参考温度须与温移表的单位移位点一致。
可以不填温移文件且在参考温度计算；也可手填shift_factors={'aT':...,'bT':...}，
与shift_factors_file二选一。参考温度为必填元数据，程序不从材料名称猜测。
温移表的log(aT)、log(bT)按温度线性插值，不外推。
模量在log频率上用保形插值；正模量/损耗用log值，含零损耗时直接插值损耗。
默认frequency_policy='error'：移位后频率覆盖不足则停止。
显式constant_endpoint仅用于诊断，超出范围时记录out_of_band，不能当作测得的外推。
有限带宽数据的独立插值不保证Kramers--Kronig因果一致性。

## 运行已知解析例子

    python examples/viscoelastic/run_example.py

默认results/run须为空。可用 --output 指定另一空目录，--temperature 10 改温度。
示例材料是明确标记的合成标准线性固体，不是用户橡胶：
E0=1MPa, Einf=10MPa, tau=0.01s，参考温度20C。
E'=E0+(Einf-E0)(omega*tau)^2/(1+(omega*tau)^2)，
E''=(Einf-E0)omega*tau/(1+(omega*tau)^2)。配套温移表为演示假设，非实测。
master_curve.dat、shift_factors.dat及roughness.dat均在示例目录中。

## 日常参数与开关

main.py打开RUN_VISCOELASTIC_FRICTION，关闭所有其他任务及独立流动开关。
material.model='viscoelastic'，填写master_curve_file、poisson_ratio；
删除young_modulus_pa/effective_modulus_pa或设为None，不允许同时指定固定模量。
宏观构型中显式设置deformable_body='lower'或'upper'，粗糙度必须位于相反的刚体侧。
此任务中上下位置只标明材料与粗糙度归属；速度是两体相对速度大小。
不求宏观圆柱压力/变形，所以输入的是viscoelastic.pressures_pa，不是总力或线载荷。
参考示例脚本中的完整参数字典，可直接换成自己的文件与温度。

## 理论与摩擦分解

采用Persson等温、稳态、各向同性单侧刚性粗糙面对光滑黏弹半空间的公式：
https://doi.org/10.1063/5.0293616 的式17--20。

G(q)=1/8 integral[q0,q] dq q³ C(q) integral[0,2pi] dphi
 |E(q v cos(phi))/((1-nu²)p)|²；P(q)=erf(1/(2sqrt(G)))。
S(q)=kappa+(1-kappa)P(q)²，默认kappa=0.5（与旧间隙能量gamma独立）。
mu_visco=1/2 integral dq q³ C(q) P(q)S(q) integral dphi
 cos(phi) Im[E(qv cos(phi))/((1-nu²)p)]。
四象限对称性在代码中显式处理。损耗输入为非负幅值；复数API采用exp(+i omega t)，
负频率为正频率的共轭，与原版另一傅里叶号约定的物理耗散相同。

mu_adhesive=tau*P(q1)/p，代表真实接触面积上的界面剪切，不是法向黏附力。
本任务为干态，因此mu_fluid=0；mu_total=mu_viscoelastic+mu_adhesive。
这些系数均按同一名义压力归一化；不用接触面积内平均应力替代名义压力。

谱积分拆分每个PSD折点，用各小区间精确的q³C谱矩及中点材料响应，累加G。
角积分使用Gauss节点，分别加倍谱分辨率和角分辨率核查mu和面积；
核查不改变主结果。只有压力、速度变化时同一速度的材料积分在内存复用。

## 验证及限制

可独立运行 `python examples/viscoelastic/verify.py` 重复7项检查。

测试覆盖SI转换、温移、模量类型、频率越界、无耗散弹性极限、正反向对称、
全接触角积分系数、以及与独立DOP853连续ODE积分对照。
原NBR格式文件可读取，但未将其与同工况Fortran的摩擦输出完成对照；不能声称原版复现。
本阶段未实现宏观黏弹变形/滚动损耗、液膜耦合、热、应变软化、法向黏附或双侧粗糙。
软体自身粗糙不能直接套用刚性粗糙扫过软体的损耗公式，因此该组合明确拒绝，而不把所有黏弹损耗设为0。
原混合润滑和tau拟合仍拒绝viscoelastic；禁止把独立干摩擦结果不加说明地加到旧液膜解。

后续混合润滑的目标分解为mu_viscoelastic+mu_adhesive+mu_fluid，但还需要
自洽地更新接触面积、平均间隙/输运闭合及宏观材料响应；不能只把E换成复数或加一个摩擦常数。
