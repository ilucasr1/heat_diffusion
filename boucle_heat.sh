#!/bin/bash

# Affichage de la commande
echo "Lancement de la commande avec mpiexec --n 50 --hostfile \$OAR_NODEFILE heatsink_v1"

# Exécution de la commande
mpiexec --n 50 --hostfile $OAR_NODEFILE heatsink_v1

# Ajout d'une ligne vide pour la clarté
echo ""

# Boucle de 10 itérations
for ((i=8; i<=20; i++)); do
    # Calcul de la valeur de X
    X=$(((i==0)?1:i * 5))
    
    # Affichage de la commande
    echo "Lancement de la commande avec mpiexec --n $X --hostfile \$OAR_NODEFILE heatsink_v1"
    
    # Exécution de la commande
    mpiexec --n $X --hostfile $OAR_NODEFILE heatsink_v1

    # Ajout d'une ligne vide pour la clarté
    echo ""
done

