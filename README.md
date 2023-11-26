# heat_diffusion
## Compilation
Compile with 
```
make
``` 
then execute with 
```
make run NB_PROC=p
``` 
where p is the number of processus you want to run the simulation on

## Display
To print the image, follow the instructions printed in the terminal at the end of the execution

## Save Format

first line : 0 = no save OR p = save with p the number of processus we used
seconde line : time and n_steps
then 1 blocks per processus (in the order of their rank) :
    1st line : my_rank
    then 1 block per step on z axe :
        1st line : position on z axe
        then 1 line per step on y axe :
            1 value per step on x axe


example :
```
4
144.000 36000
0
z = 0
293.150000 293.150000 ... 293.150000 293.150000
.
.
.
293.150000 293.150000 ... 293.150000 293.150000
z = 0.002
.
.
.
z = 0.056
307.774269 307.815516 ... 303.161329 303.048255
.
.
.
307.774268 307.815516 ... 303.161329 303.048255
1
.
. 
.
3
.
.
.
```