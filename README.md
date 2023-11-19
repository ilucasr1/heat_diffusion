# heat_diffusion

AMELIORATIONS :

1) dans do_halo
on attend de tout recevoir pour traiter.
essayer de traiter des qu'on recoit

2) tout reecrire dans un unique tableau T au lieu de contenir dans T et reorganiser dans R (pour le gather)


REPRENDRE : 

1) Il faut modifier le remplissage de chunk pour les gather (cas dim >= 2) en reorganisant

puis faire de meme pour dim >= 1