******************
source codes by EB
******************

These code files are stored in *src2* directory



File CMB_MeandT.py:

Computes the Mean of the CMB Temperature fluctuations within a sky
region encompassing the union of areas of disks of sizes
r_disk = A_k * 10**(r_ext) centred at the position of the galaxies.
The samples employed are galaxies of 2MASS catalog separated by
morphological type Sa, Sb, and Sc using the four CMB maps with
fourground substraction, namely SMICA, SEVEM, NILC and Commander.


File CMB_MeandT_random.py:

Computes the Mean of the CMB Temperature fluctuations within a sky
region encompassing the union of areas of disks of sizes
r_disk = A_k * 10**(r_ext) centred at the position of random galaxies.
Random samples have the same number of objects than the samples
Sa, Sb, and Sc with the same distribution of angular sizes
r_disk = A_k * 10**(r_ext).
The four CMB maps with fourground substraction, namely SMICA,
SEVEM, NILC and Commander are also employed.

File CMB_MeandT_rings.py:

Computes the Mean of the CMB Temperature fluctuations within rings
of external radius of sizes r_disk = A_k * 10**(r_ext) and width
"bins_sz" centred at the position of the galaxies.
The samples employed are galaxies of 2MASS catalog separated by
morphological type Sa, Sb, and Sc using the four CMB maps with
fourground substraction, namely SMICA, SEVEM, NILC and Commander.

File CMB_MeandT_rings_random.py:
   
Computes the Mean of the CMB Temperature fluctuations within rings
of external radius of sizes r_disk = A_k * 10**(r_ext) and width
"bins_sz" centred at the position of random galaxies.
The samples employed are galaxies of 2MASS catalog separated by
morphological type Sa, Sb, and Sc using the four CMB maps with
fourground substraction, namely SMICA, SEVEM, NILC and Commander.

File MeanT_of_disks_ring.py:
   
Computes the Mean of the Mean CMB Temperature fluctuations within rings
of external radius of sizes r_disk = A_k * 10**(r_ext) and width
"bins_sz" centred at the position of random galaxies.
The samples employed are galaxies of 2MASS catalog separated by
morphological type Sa, Sb, and Sc using the four CMB maps with
fourground substraction, namely SMICA, SEVEM, NILC and Commander.

