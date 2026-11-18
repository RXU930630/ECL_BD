# ECL_BD 0.2.0rc4 — 本次硅油实验

本分支包含完整运行程序、已填写的parameters.py，以及用户要求上传的27点Excel实验数据（40、100、160N各9个速度）。没有历史测试结果或已计算的拟合结果。独立剪切强度拟合开关已开启。

rc4已将积分自动加密上限改为4次，修正40N临界节点的收敛中断，保持0.002相对容差。输出目录改为新的silicone_fit_excel_02，避免覆盖原失败诊断。

## 下载与运行

在本分支页面点击绿色Code按钮，再点击Download ZIP；解压完整目录。在解压目录打开终端，执行：

```bash
python -m pip install -r requirements.txt
python main.py
```

已填好的数据位于inputs/experiment.xlsx，可以直接用Excel打开编辑，保持第1行表头不变。无需再转换CSV或复制其他文件。

拟合结果保存在results/silicone_fit_excel_02/，用Excel打开fit_results.xlsx。拟合概况页显示剪切强度Pa和MPa，各实验点页显示固体/流体分项、面积和残差。计算可能需要数分钟；输出目录非空时请在parameters.py中换一个新名称。

## 已填写的工况

- 成形外半径0.17m，轴向长度0.10m，流向投影长度0.03m，厚度0.02m仅记录。
- E=1.8MPa，泊松比0.5，动力黏度1Pa·s，无滑移，两端储液表压0。
- 法向载荷逐行读取Excel，自动按轴向长度换算线载荷。速度取Excel中的实验速度。
- 接触剪切强度保持未知，由拟合任务求出；尚未进行实验精度或迁移性验证。

程序仍使用平面应变半空间弹性近似。详细使用与物理范围见docs/local_run.md、docs/finite_geometry.md及docs/shear_strength_fit.md。这些通用说明中的“纯净包”指local-trial-v0.2.0rc4分支；本分支额外包含本次实验输入。
