import matplotlib.pylab as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams["figure.autolayout"] = True

fig = plt.figure(figsize=(10,10))

data = pd.read_csv("/home/toberoi/hpc/project-3-winter-2023-tinaoberoi/nbody.csv", header = None, names=['iter', 'x', 'y', 'z'])

all_iterations = data.iter.unique()

data_dict = {elem: pd.DataFrame() for elem in all_iterations}

for key in data_dict.keys():
	data_dict[key] = data[:][data.iter==key]

# for iter_key, data_val in data_dict.items():
#     fig = plt.figure(figsize=(10, 10))
#     ax = plt.axes(projection='3d')
#     ax.scatter3D(data_val['x'], data_val['y'], data_val['z'])
#     plt.title(f"Iteration {iter_key}")
#     plt.show()
    
def animate(i):
	data_val = data_dict[i+1]
	plt.clf()

	ax = plt.axes(projection='3d')

	ax.set_xlim3d([-10.0, 10.0])
	ax.set_xlabel('X')

	ax.set_ylim3d([-10.0, 10.0])
	ax.set_ylabel('Y')

	ax.set_zlim3d([-10.0, 10.0])
	ax.set_zlabel('Z')

	ax.scatter3D(data_val['x'], data_val['y'], data_val['z'])


ani = FuncAnimation(fig, animate, repeat=True, frames=20)

ani.save("thing.gif", writer=PillowWriter(fps=25))