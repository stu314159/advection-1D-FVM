function phi = van_Leer(theta)

phi =(theta + abs(theta))/(1 + abs(theta));