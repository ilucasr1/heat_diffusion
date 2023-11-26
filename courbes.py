import matplotlib.pyplot as plt
import numpy as np

# Générer des données pour trois courbes
# x = np.array([i*5 for i in range(21)])
# fast   = np.array([12, 5, 6, 8, 8, 8, 9, 10, 11, 12, 11, 11, 11, 12, 13, 12, 14, 14, 14, 15, 16])
# medium = np.array([136, 38, 29, 25, 29, 19, 20, 19, 20, 23, 19, 20, 22, 20, 20, 23, 23, 22, 25, 22, 24])
# normal = np.array([0, 0, 270, 220, 170, 140, 133, 128, 162, 185, 148, 150, 166, 146, 140, 151, 150, 121, 132, 118, 125])

# # Tracer les trois courbes sur le même graphique
# plt.plot(x, fast, label='FAST')
# plt.plot(x, medium, label='MEDIUM')
# plt.plot(x, normal, label='NORMAL')

# # Ajouter des titres et légendes
# plt.title('3D simulation on FAST, MEDIUM and NORMAL mode')
# plt.xlabel('time (s)')
# plt.ylabel('processus')
# plt.legend()

# # Afficher le graphique
# plt.show()





# # Générer des données pour trois courbes
# x = np.array([i*5 for i in range(21)])
# y_1D   = np.array([])
# y_2D = np.array([])
# y_3D = np.array([0, 38, 29, 25, 29, 19, 20, 19, 20, 23, 19, 20, 22, 20, 20, 23, 23, 22, 25, 22, 24])

# # Tracer les trois courbes sur le même graphique
# plt.plot(x, y_1D, label='1D')
# plt.plot(x, y_2D, label='2D')
# plt.plot(x, y_3D, label='3D')

# # Ajouter des titres et légendes
# plt.title('1D, 2D and 3D simulation on MEDIUM mode')
# plt.xlabel('time (s)')
# plt.ylabel('processus')
# plt.legend()

# # Afficher le graphique
# plt.show()

# Générer des données pour trois courbes
x = np.array([i*5 for i in range(21)])
y_1D = np.array([137, 38, 36, 28, 29, 31, 34, 28, 28, 28, 30, 28, 29, 29, 29, 27, 28, 28, 29, 30, 31])
y_2D = np.array([130, 37, 28, 28, 31, 23, 24, 22, 23, 22, 23, 24, 23, 24, 24, 24, 23, 25, 27, 25, 25])
y_3D = np.array([132, 38, 29, 25, 29, 19, 20, 19, 20, 23, 19, 20, 22, 20, 20, 23, 23, 22, 25, 22, 24])

# Tracer les trois courbes sur le même graphique
plt.plot(x, y_1D, label='1D')
plt.plot(x, y_2D, label='2D')
plt.plot(x, y_3D, label='3D')

# Ajouter des titres et légendes
plt.title('1D, 2D and 3D simulation on MEDIUM mode')
plt.xlabel('time (s)')
plt.ylabel('processus')
plt.legend()

# Afficher le graphique
plt.show()

#medium de 100 à 150 [124, 135, 121, 125, 150, 125, 133, 130, 140, 147]


#medium 1D de 1 à 150

#medium 2D de 1 à 150[136, 37, 28, 28, 31, 23, 24, 22, 23, 22, 23, 24, 23, 24, 23, 25, 31, 25, 29, 29, 30, 31, 32, 34]