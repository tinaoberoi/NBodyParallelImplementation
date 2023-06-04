import matplotlib.pyplot as plt

avg_time = [17.129681, 8.567937, 4.296844, 2.215315, 1.726375, 1.166894]
tot_time = [325.463935, 162.790811, 81.640037, 42.090989, 32.801116, 22.170987]
cores = [4, 8, 16, 64, 128, 256]

plt.xlim(0, 260)
plt.ylim(1.166, 1300)

plt.plot(cores, avg_time, 'o:r', label="Average Time")
plt.plot(cores, tot_time, 'o:b', label="Total Time"),

plt.title('Scaling Plot for 100000 body configurations')
plt.xlabel("Cores")
plt.ylabel("Time Taken (sec)")
plt.legend(loc="upper right")
plt.savefig('openmp.png')