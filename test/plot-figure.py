import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# ===========================================================
# 创建数据点
x = np.linspace(0, 1, 100)  # 在0到1之间生成100个x坐标点
y = np.linspace(0, 1, 100)  # 在0到1之间生成100个y坐标点
x, y = np.meshgrid(x, y)    # 创建坐标网格

z = np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y)  # 计算z坐标

# 绘制3D图形
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot_surface(x, y, z, cmap='viridis')  # 绘制3D曲面

# 设置坐标轴标签
ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$')
ax1.set_zlabel('$z$')


# ===========================================================
# 读取CSV文件
data = pd.read_csv("data/numerical_result_4.csv")

# 提取x、y、z列的数据
x = data['x']
y = data['y']
z = data['z']
# 创建三维图形对象
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

# 绘制三维表面图
ax2.plot_trisurf(x, y, z, cmap='viridis')

# 设置坐标轴标签
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

# 显示图形
plt.show()