代码包括
convection-diffusion例子, 见红宝书11.2节

格式
1. 对流项：Upwind, CD, Power-law, SOU (to be implemented), Quick (to be implemented)
2. 粘性项：CD (继承自Python_heat_conduction的例子)

单独运行代码
python3 convection_diffusion.py
python3 convection_diffusion2.py

脚本运行代码
./job_Pe_L_center (该脚本内部实际运行文件convection_diffusion.py，用以复现红宝书上的Fig. 11.7)
./job_Pe_L (该脚本内部实际运行文件convection_diffusion2.py，用以复现红宝书上的Fig. 11.2)
