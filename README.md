# Protein-glycan docking

## Software Requirements
If you are a University of Mississippi member, the following software are required to complete this tutorial:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php)\
\
Other users are required to have:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/), [AutoDock Vina](https://vina.scripps.edu/), [Vina-carb](http://legacy.glycam.org/docs/othertoolsservice/downloads/downloads-software/index.html) and [GlycoTorch Vina](https://github.com/EricBoittier/GlycoTorch-Vina)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php), [AutoDock Vina](https://vina.scripps.edu/)

Vina-carb and GlycoTorch Vina software are not available for windows OS.

## Protein-glycan docking
Docking of Lewis Y Tetrasaccharide to Humanized Fab using AutoDock Vina and Vina-carb. Structure of the complaex is availabe in the protein data bank under PDB ID 1S3K(https://www.rcsb.org/structure/1S3K)

```
1S3K.pdb
receptor.pdb
LeY-xray.pdb
LeY-glycam.pdb
config_vc.txt
config_vina.txt
```
```
module load mgltools/v2.1.5.7 
prepare_ligand4.py -l LeY-xray.pdb -A hydrogens
prepare_ligand4.py -l LeY-glycam.pdb -A hydrogens



prepare_ligand4.py -l ligand.pdb -o ligand.pdbqt -A hydrogens
prepare_ligand4.py -l ligand.pdb -o ligand.pdbqt -A hydrogens


prepare_receptor4.py -r receptor.pdb -o receptor_rigig.pdb -A "hydrogens"
prepare_flexreceptor4.py -r receptor.pdbqt -s receptor:H:TYR32_TYR33_TYR50_TRP105 

```

```
module load vina-carb/v1.2 
module load autodock-vina
```



```
[sushil@fucose flex_lig]$ vina --config config_vina.txt
#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#                                                               #
# O. Trott, A. J. Olson,                                        #
# AutoDock Vina: improving the speed and accuracy of docking    #
# with a new scoring function, efficient optimization and       #
# multithreading, Journal of Computational Chemistry 31 (2010)  #
# 455-461                                                       #
#                                                               #
# DOI 10.1002/jcc.21334                                         #
#                                                               #
# Please see http://vina.scripps.edu for more information.      #
#################################################################

Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 0
Performing search ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -7.6      0.000      0.000
   2         -7.4      2.842      8.087
   3         -7.2      2.292      6.803
   4         -7.1      2.377      4.747
   5         -7.1      2.242      6.845
   6         -6.6      2.034      6.669
   7         -6.6      2.217      7.389
   8         -6.4      2.350      4.405
   9         -6.4      1.241      4.338
  10         -6.2      2.007      6.492
Writing output ... done.


```

```
[sushil@fucose flex_lig]$ vina-carb --config config_vc.txt 
#################################################################
#	 		Vina-Carb Results	 	     	#
#################################################################

WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs
Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 0
Performing search ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
done.
Refining results ... 	BEFORE corresponding chi_energy: 0
	BEFORE corresponding chi_energy: 0
	BEFORE corresponding chi_energy: 0
	.
  .
  .
	BEFORE corresponding chi_energy: 8.41328
	BEFORE corresponding chi_energy: 10.815
	BEFORE corresponding chi_energy: 14.0039
	BEFORE corresponding chi_energy: 16.9191
	BEFORE corresponding chi_energy: 11.793
done.

mode |   affinity |   chi   |  affinity | dist from best mode
     | (kcal/mol) |  energy |  - chi    | rmsd l.b.| rmsd u.b.
-----+------------+---------+--------------------+------------
   1         -7.6      0.0      -7.6      0.000      0.000
   2         -6.9      0.0      -6.9      2.397      7.545
   3         -6.8      0.0      -6.8      1.919      5.885
   4         -6.7      0.0      -6.7      2.204      6.004
   5         -6.5      0.0      -6.5      1.829      5.924
   6         -6.4      0.0      -6.4      2.322      6.022
   7         -6.4      0.0      -6.4      2.862      8.269
   8         -6.2      0.0      -6.2      2.376      8.255
   9         -6.1      0.0      -6.1      2.840      6.246
  10         -6.1      0.0      -6.1      2.176      8.020
Writing output ... done.
```



```
[sushil@fucose flex_rec]$ vina --config config_vina.txt 
#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#                                                               #
# O. Trott, A. J. Olson,                                        #
# AutoDock Vina: improving the speed and accuracy of docking    #
# with a new scoring function, efficient optimization and       #
# multithreading, Journal of Computational Chemistry 31 (2010)  #
# 455-461                                                       #
#                                                               #
# DOI 10.1002/jcc.21334                                         #
#                                                               #
# Please see http://vina.scripps.edu for more information.      #
#################################################################

Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 0
Performing search ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -7.4      0.000      0.000
   2         -7.4      2.218      6.436
   3         -7.1      2.679      6.334
   4         -7.1      1.627      6.702
   5         -6.9      2.672      5.984
   6         -6.5      3.311      6.055
   7         -6.5      2.336      5.568
   8         -6.3      3.296      5.312
   9         -6.2      2.155      5.361
  10         -6.2      2.393      5.291
Writing output ... done.

```

```
[sushil@fucose flex_rec]$ vina-carb --config config_vc.txt 
#################################################################
#	 		Vina-Carb Results	 	     	#
#################################################################

WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs
Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 0
Performing search ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
done.
Refining results ... 	BEFORE corresponding chi_energy: 0
	BEFORE corresponding chi_energy: 0
	BEFORE corresponding chi_energy: 0
	.
  .
  .
  BEFORE corresponding chi_energy: 10.1698
	BEFORE corresponding chi_energy: 19.245
done.

mode |   affinity |   chi   |  affinity | dist from best mode
     | (kcal/mol) |  energy |  - chi    | rmsd l.b.| rmsd u.b.
-----+------------+---------+--------------------+------------
   1         -7.6      0.0      -7.6      0.000      0.000
   2         -7.0      0.0      -7.0      1.409      4.556
   3         -7.0      0.0      -7.0      1.848      5.843
   4         -6.6      0.0      -6.6      1.633      4.717
   5         -6.5      0.0      -6.5      1.828      3.665
   6         -6.4      0.0      -6.4      1.873      4.753
   7         -6.3      0.0      -6.3      1.809      5.279
   8         -6.2      0.0      -6.2      1.817      5.066
   9         -6.2      0.0      -6.2      1.394      4.693
  10         -6.1      0.0      -6.1      2.455      6.065
Writing output ... done.

```

