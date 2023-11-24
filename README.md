# heat_diffusion

AMELIORATIONS :

1) dans do_halo
on attend de tout recevoir pour traiter.
essayer de traiter des qu'on recoit

2) tout reecrire dans un unique tableau T au lieu de contenir dans T et reorganiser dans R (pour le gather)

3) implementer checkpoint

4) changer association divider <-> dim_chunk (pb si div[1] > 2)

5) tenir compte des pb aux bords si divider ne divise pas n/m/o



REPRENDRE : 

1) comprendre d'ou vient le decalage de 0.1

2) NETTOYER LE CODE

3) copier code de grid

puis faire de meme pour dim >= 1




FORMAT save.txt

0 : no save / 1 : save
time && n_steps
my_rank
z0
.
.
.
zn
my_rank
.
. 
.