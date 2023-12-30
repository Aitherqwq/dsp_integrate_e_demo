# dsp_integrate_e_demo
数字信号处理大作业，语音处理

## 作业要求
- 采集一段语音信号，听一听播放效果；对其进行频谱分析，画出信号的时域波形图和频谱图。
- 生成2个频率分别为50Hz和1.8KHz的正弦信号，作为噪声加入语音信号中，听一听播放效果；画出信号的时域波形图和频谱图；比较加噪前后的语音信号，分析发生的变化。
- 设计合适的滤波器，对带有噪声的语音信号进行滤波，听一听播放效果；画出滤波前后的时域波形图和频谱图，比较滤波前后的语音信号，分析发生的变化。
- 把正弦型噪声换成高斯随机噪声，重新做一下。分析两次实验结果的异同和原因。
- 混响，利用电脑或手机的录音功能，采集一段自己的音频，并进行频谱分析；给这段语音叠加上延迟，产生混响（回声）效果，再进行频谱分析。
- 语速改变，利用电脑或手机的录音功能，采集一段自己的音频，通过处理，改变语速。
- 变音，采集一段自己的声音，实现男声变女声，女声变男声的效果。

## 功能实现
加噪、FIR滤波器滤波、维纳滤波器滤波、混响、变音。
dsp_integrate_e_demo.mlapp是有GUI的主程序。包含的音频文件可以作为示例音频使用。

## 开源代码
wiener_as函数由开源代码修改   Copyright (c) 2006 by Philipos C. Loizou, MIT License
