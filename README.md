## Software Requirements
If you are a University of Mississippi member, the following software are required to complete this tutorial:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php)\
\
Other users are required to have:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/), [AutoDock Vina](https://vina.scripps.edu/), [Vina-carb](http://legacy.glycam.org/docs/othertoolsservice/downloads/downloads-software/index.html) and [GlycoTorch Vina](https://github.com/EricBoittier/GlycoTorch-Vina)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php), [AutoDock Vina](https://vina.scripps.edu/)

Vina-carb and GlycoTorch Vina software are not available for windows OS.

# Protein-glycan docking
https://autodock-vina.readthedocs.io/en/latest/docking_flexible.html

https://www.click2drug.org/index.php#Binding%20free%20energy%20estimation 

https://casfaculty.fiu.edu/David.Chatfield/workshop/materials/drug-design/fiu-docking-tutorial.pdf
https://personal.utdallas.edu/~son051000/comp/EdelmiroMoman.pdf

https://www.bu.edu/chemistry/files/2010/11/Docking.20101122.pdf

Docking of Lewis Y Tetrasaccharide to Humanized Fab using AutoDock Vina and Vina-carb. Structure of the complaex is availabe in the protein data bank under PDB ID 1S3K(https://www.rcsb.org/structure/1S3K)


Open PuTTY and connect to Fucose workstation using your credencials:
```
ssh -X guest01@fucose.pharmacy.olemiss.edu
```
##Connect to remote workstation:

All input files have been prepared and placed under ./practice directory for your convenience . 

Login to fucose using the instructions provided during the lecture. Linux/Mac users can use terminal to connect to ccbrc workstation, whereas windows users should use PyTTY to connect to the ccbrc workstation.
```
leo:~ sushil$ ssh -X guestXX@machine.host.name
guestXX@machine.host.name's password: 
Last login: Tue Jul 19 13:33:49 2022 from idose.pharmacy.olemiss.edu

##########################################################################
##									##
##    Computational Chemistry and Bioinformatics Research Core (CCBRC)	##
##									##
##			Support: sushil@olemiss.edu			##
##									##
##----------------------------------------------------------------------##
## 	Access to this machine is strictly for research and to  	##
##	authorized users only. 						## 
##									##
##########################################################################
-bash-4.2$ pwd

```
Now if you will type "ls -l" and hit enter, there should be two directories "practice and tutorial" available to everyone.
```
-bash-4.2$ ls -l
total 8
drwxr-xr-x. 2 sushil cgw 4096 Dec 15  2021 practice
drwxrwxr-x. 2 sushil cgw 4096 Dec 16  2021 tutorial
```

Directory "tutorial has all the precalculated data if you want to look into the correct output files. Another directory named "practice" has only input files for docking, and you can run calculations under this directory. To do so, type cd "./practice" and hit enter. Then type "ls -l," and it should list two direcotiries:
```
flex_lig  : input files for rigid receptor + flexible ligand docking  
flex_rec  " input files for flexible receptor + flexible ligand docking  
```
To use any software in this machines you need to load them as module. All teh available software can be listed using follwoing command:
```
module avail

------------------------------------ /usr/share/Modules/modulefiles ------------------------------------
dot         module-git  module-info modules     null        use.own

------------------------------------------- /etc/modulefiles -------------------------------------------
amber/20                     fesetup                      netpbm
**autodock-vina**                glycotorch-vina              pymol/v2
basecalling/filtlong         guppy/gpu-6.0.1              rosetta/2020.50.61505
basecalling/flye-2.9         maxcluster/0.6.6             schrodinger/2020.4
basecalling/miniasm-0.3-r179 **mgltools/v2.1.5.7**            sire/2020.1
basecalling/polypolish       modeller/10.0                smina/1.2.2
basecalling/porechop         modeller/9.25                spicker/v1.0
basecalling/pycoqc           mpi/mpich-3.0-x86_64         **vina-carb/v1.2**
boost                        mpi/mpich-x86_64             vmd/1.9.3
cresset/Flare                msub/v1.0                    xscore/v1.2.1
cresset/Forge                naccess1/v2.1.1
cresset/Spark                namd/3.0alpha

```
## Prepare flexible ligand and rigid receptor input files for docking
Change direcotry to 'flex_lig' and run follwoing commands to perform docking:

```
module load autodock-vina 				# Load AutoDock Vina to your working environment
module load vina-carb/v1.2 				# Load Vina-carb to your working environment
module load mgltools/v2.1.5.7				# Load MGL Tools to your working environment

prepare_ligand4.py -l LeY-xray.pdb -A hydrogens 	# Convert LeY-xray.pdb into LeY-xray.pdbqt (Vina Input file)
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A "hydrogens"	# Convert receptor.pdb in receptor.pdbqt (Vina Input file)
```
If you see an error message "**init.c(556):ERROR:161: Cannot initialize TCL**" run follwoing command to fix it
```
unset LD_LIBRARY_PATH

```
## Perform docking using AutoDock Vina andvina-carb
Now run docking using Vina and vina-carb as below:

```
$ vina --config config_vina.txt   	#this will run vina using input parameters in 'config_vina.txt' file

Output:
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
This run should finish within a minute and file "LeY-xray_vina_out.pdbqt" containing all the docked poses will be generated. Now repeate docking using Vina-carb and input parameters present in "config_vc.txt" file (as below). It should generate a file "LeY-xray_vc_out.pdbqt" taht contains all the docking poses. 

```
$ vina-carb --config config_vc.txt   	#this will run vina-carb using input parameters in 'config_vc.txt' 

Output:
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
## Visulizing docking poses:

Now, copy _receptor.pdb_, _LeY-xray_vc_out.pdbqt_, and _LeY-xray_vina_out.pdbqt_ into back into your computer and load them in PyMOL.


# Flexible receptor docking 

Change direcotory to 
```
$ cd ~/practice/flex_rec
$ ls -l 

total 860
-rw-r--r--. 1 sushil cgw 338094 Jul 19 13:45 1S3K.pdb
-rw-r--r--. 1 sushil cgw    619 Jul 19 13:45 config_vc.txt
-rw-r--r--. 1 sushil cgw    580 Jul 19 13:45 config_vina.txt
-rw-r--r--. 1 sushil cgw   7371 Jul 19 13:45 LeY-xray.pdb
-rw-r--r--. 1 sushil cgw 524003 Jul 19 13:45 receptor.pdb
```


## Prepare flexible ligand and rigid receptor input files for docking

```
$ prepare_ligand4.py -l LeY-xray.pdb -A hydrogens
$ prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A hydrogens
$ prepare_flexreceptor4.py -r receptor.pdbqt -s receptor:H:TYR32_TYR33_TYR50_TRP105 #prepare flexible and rigid receptor input files
$ ls -l
total 1508
-rw-r--r--. 1 sushil cgw 338094 Jul 19 13:45 1S3K.pdb
-rw-r--r--. 1 sushil cgw    619 Jul 19 13:45 config_vc.txt
-rw-r--r--. 1 sushil cgw    580 Jul 19 13:45 config_vina.txt
-rw-r--r--. 1 sushil cgw   7371 Jul 19 13:45 LeY-xray.pdb
-rw-r--r--. 1 sushil cgw   6483 Jul 19 15:04 LeY-xray.pdbqt
-rw-r--r--. 1 sushil cgw   4767 Jul 19 15:05 receptor_flex.pdbqt
-rw-r--r--. 1 sushil cgw 524003 Jul 19 13:45 receptor.pdb
-rw-r--r--. 1 sushil cgw 323496 Jul 19 15:04 receptor.pdbqt
-rw-r--r--. 1 sushil cgw 320080 Jul 19 15:05 receptor_rigid.pdbqt

```
Please not that '-s' option in prepare_flexreceptor4.py script allows user to procvide residues which you want to keep flexible during the docking. In this case we are keeping four amino-acid residues _Tyr3 Tyr33 Tyr50 and Trp105_ flexible during the docking. In general one should keep only those residue to be flexible whose conformational changes can affect ligand binding. Output file _receptor_flex.pdbqt_ and _receptor_rigid.pdbqt_ contains flexible and rigid part of the receptor for docking.

## Perform docking using AutoDock Vina andvina-carb
Now run docking using Vina and vina-carb as below:

```
$ vina --config config_vina.txt   	#this will run vina using input parameters in 'config_vina.txt' file
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
   2         -7.4      2.228      6.434
   3         -7.1      2.702      6.362
   4         -7.1      1.622      6.723
   5         -6.9      2.668      5.994
   6         -6.6      3.388      5.927
   7         -6.4      2.305      6.497
   8         -6.3      2.856      6.635
   9         -6.3      3.289      5.325
  10         -6.1      2.488      5.993
Writing output ... done.
```
Flexible ligand docking should take longer time compared to rigid receptor docking. Repeat the same claulcation using vina-carb and you hsould get finally two outfiles _LeY-xray_vina_flex_out.pdbqt_ and _LeY-xray_vc_flex_out.pdbqt_. Copy both these files in your local computer and analyze the binding modes. 















## Protein-glycan docking

```
1S3K.pdb
receptor.pdb
LeY-xray.pdb
LeY-glycam.pdb
config_vc.txt
config_vina.txt
```
```

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
[sushil@fucose flex_lig]$ 

```

```
[sushil@fucose flex_lig]$ 
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

