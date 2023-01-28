import matplotlib.pyplot as plt


Core = [1, 4, 8, 12, 16]

Ideal_2 = [2, 8, 16, 24, 32]
Ideal_4 = [4, 16, 32, 48, 64]
Ideal_6 = [6, 24, 48, 72, 96]
Ideal_8 = [8, 32, 64, 96, 128]

fig = plt.figure()

plt.subplot_tool()

N2 = fig.add_subplot(2, 2, 1)
N4 = fig.add_subplot(2, 2, 2)
N6 = fig.add_subplot(2, 2, 3)
N8 = fig.add_subplot(2, 2, 4)


N2_32 = [1.55, 2.58, 2.06, 1.32, 1.11]
N2_1000 = [1.90, 6.84, 11.22, 14.59, 16.26]
N2_10000 = [1.90, 7.32, 13.11, 18.36, 23.38]
#N4 = [2.18, 2.04, 1.50, 1.25, 0.89]
#N6 = [2.47, 1.87, 1.37, 0.99, 0.79]
#N8 = [2.52, 1.52, 1.07, 0.89, 0.62]

N2.scatter(Core,Ideal_2)
N2.scatter(Core,N2_32)
N2.scatter(Core,N2_1000)
N2.scatter(Core,N2_10000)
#plt.scatter(Core,N8)


N2.plot(Core,Ideal_2,label='Linear ')
N2.plot(Core,N2_32,label='32 Body')
N2.plot(Core,N2_1000,label='1000 Body')
N2.plot(Core,N2_10000,label='10000 Body')
#plt.plot(Core,N8,label='8 Node')

N2.set_title('2 node')
N2.set_xlabel("Cores")
N2.set_ylabel("Speedup")
N2.legend()

#N4
N4_32 = [2.18, 2.04, 1.50, 1.25, 0.89]
N4_1000 = [3.73, 11.14, 14.97, 15.98, 14.41]
N4_10000 = [3.81, 14.04, 24.69, 32.96, 37.70]


N4.scatter(Core,Ideal_4)
N4.scatter(Core,N4_32)
N4.scatter(Core,N4_1000)
N4.scatter(Core,N4_10000)
#plt.scatter(Core,N8)


N4.plot(Core,Ideal_4,label='Linear ')
N4.plot(Core,N4_32,label='32 Body')
N4.plot(Core,N4_1000,label='1000 Body')
N4.plot(Core,N4_10000,label='10000 Body')
#plt.plot(Core,N8,label='8 Node')

N4.set_title('4 node')
N4.set_xlabel("Cores")
N4.set_ylabel("Speedup")
N4.legend()

#N6
N6_32 = [2.47, 1.87, 1.37, 0.99, 0.79]
N6_1000 = [5.39, 12.92, 15.21, 13.26, 10.76]
N6_10000 = [5.65, 20.26, 34.93, 43.38, 48.60]


N6.scatter(Core,Ideal_6)
N6.scatter(Core,N6_32)
N6.scatter(Core,N6_1000)
N6.scatter(Core,N6_10000)
#plt.scatter(Core,N8)

N6.plot(Core,Ideal_6,label='Linear ')
N6.plot(Core,N6_32,label='32 Body')
N6.plot(Core,N6_1000,label='1000 Body')
N6.plot(Core,N6_10000,label='10000 Body')
#plt.plot(Core,N8,label='8 Node')

N6.set_title('6 node')
N6.set_xlabel("Cores")
N6.set_ylabel("Speedup")
N6.legend()

#N8
N8_32 = [2.52, 1.52, 1.07, 0.89, 0.62]
N8_1000 = [6.99, 14.57, 13.48, 11.21, 8.28]
N8_10000 = [7.53, 25.73, 41.43, 48.68, 46.39]


N8.scatter(Core,Ideal_8)
N8.scatter(Core,N8_32)
N8.scatter(Core,N8_1000)
N8.scatter(Core,N8_10000)
#plt.scatter(Core,N8)

N8.plot(Core,Ideal_8,label='Linear ')
N8.plot(Core,N8_32,label='32 Body')
N8.plot(Core,N8_1000,label='1000 Body')
N8.plot(Core,N8_10000,label='10000 Body')
#plt.plot(Core,N8,label='8 Node')

N8.set_title('8 node')
N8.set_xlabel("Cores")
N8.set_ylabel("Speedup")
N8.legend()

plt.show()