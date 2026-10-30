# 用户提供的历史输入

从本次会话上传的compare.zip中按原字节复制，仅缩短文件名，未修改内容。

| 当前文件 | 原始文件名 |
| --- | --- |
| roughness_psd.dat | IN_ROADTYPE.q0.qr.q1.h0.H.1logWavevector.2logPowerSpectrum |
| pressure_profile.dat | IN.1SyringeName.1FluidPressure_MPa.1x_meter.2ContactPressure_MPa |

PSD名称为BD_Libertas_numerical_2，共4个谱点；q0=10、qr=1、q1=9.8e7 m⁻¹。数值谱按文件插值，不用qr替代q0。压力文件坐标单位m、接触压力单位MPa，头部外侧流体压力0.414 MPa不是额外施加的界面参考压力。

此工况的E=2.07 MPa、ν=0.5、液体黏度0.001 Pa·s。参数快照保存在同目录的parameters.json，仅由开发验证读取；日常入口只使用根目录的parameters.py。

SHA256：

```
roughness_psd.dat: abcc26e1d55e576148c30e123b232763fcc4d881683638121e06c0af9f3580d4
pressure_profile.dat: 59c5f8a6a5d6f24e694fab536888afe515804d378c26269fcad9ac979dcb1e94
```
