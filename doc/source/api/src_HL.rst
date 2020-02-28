******************
source codes by HL
******************

These code files are stored in *src1* directory


Codes by Heliana Luparello:

location:

/mnt/is2/mcb/correlations

/mnt/is2/mcb/profiles

source codes are located in *src1* directory.


Correlations
============

EN ESTE DIRECTORIO ESTÁN LOS CÓDIGOS PARA CORRER LAS
CORRELACIONES DE LOS MAPAS DE TEMPERATURA EN EN CMB,
adaptados para clemente.

* 180deg contiene:

cmb_corr_clem_All.py ---> hace las correlaciones en todo el rango de escalas del mapa
cmb_corr_all.sh   ---> es el script para mandar a correr el .py en clemente.

* 7deg contiene:

cmb_corr_clem_All_7deg.py ---> hace las correlaciones en el rango de escalas hasta 7grados
cmb_corr_highres.sh ---> es el script para mandar a correr el .py en clemente.

En todos los casos, para correrlo hay que modificar los paths de entrada y salida de datos,
porque va a escribir en un directorio en el que solo yo tengo permisos, y se va a romper.

Estos códigos también se pueden correr en mirta3 directamente, sin sistema de cola,
para hacer pruebas. para que funcionen rápido hay que cambiar el número de randoms
y el número de cores que se usan.




Profiles
========


ESTOS SON LOS CÓDIGOS PARA HACER LOS PERFILES RADIALES DE TEMPERATURA ALREDEDOR DE GALAXIAS.

en particular, este es para las galaxias late que están edge--on. eso se cambia al principio
del codigo, eligiendo convenientemente las galaxias.

esto no está paralelizado, y corre en mirta3.

cmb_functions.py ---> son las funciones para hacer los perfiles

cmb_profL05.py  ---> se encarga de levantar las rutinas, setea parámetros y corre los perfiles.




