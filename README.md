# Protein-glycan docking tutorial:

## Software Requirements

Following software are required to complete this tutorial:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php)\
\
Other users are required to have:\
Linux/Mac Users: [PyMOL](https://pymol.org/2/), [AutoDock Vina](https://vina.scripps.edu/), [Vina-carb](http://legacy.glycam.org/docs/othertoolsservice/downloads/downloads-software/index.html) and [GlycoTorch Vina](https://github.com/EricBoittier/GlycoTorch-Vina)\
Windows Users: [PyMOL](https://pymol.org/2/), [PuTTY](https://www.putty.org/), [WinSCP](https://winscp.net/eng/download.php), [AutoDock Vina](https://vina.scripps.edu/)

Vina-carb and GlycoTorch Vina software are not available for windows OS. All these programs have been preinstalled on your desktop or we will access them remotely in one of the CCBRC Workstations.

This tutorial aims to dock a tetrasaccharide to FAb using the [AutoDock Vina](https://vina.scripps.edu), and [Vina-carb](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00834). You can download all the input files by clicking on Code --> Download Zip. Unzip this file and go inside the docking directory. All these files have been placed on the computer you will be working in. 

We will break down this tutorial into five steps.

## 1. Obtain protein and glycan complex for docking: 
In this tutorial, we will be docking a Lewis Y Tetrasaccharide with Humanized Fab. First of all, download the X-ray structure of the complex from PDB [PDB ID: 1S3K](https://www.rcsb.org/structure/1s3k). The following figure shows the LeY binding to Fab as seen in its crystal structure. 
![alt text](https://github.com/glycodynamics/gly-docking/blob/main/images/Input_Fig.png)

 Now, open the PDB Structure in PyMOL and prepare the structure of the ligand-protein for docking. Perform the following structure manipulations:

#### Remove crystal waters: 
Open 1S3K.pdb file in PyMOL and click non on Action --> remove waters

#### Save ligand and protein separately
Now split the complex in protein and glycan and save them separately in two separate PDB files. \
Select LeY tetrasaccharide by clicking left mouse button on each monosaccharide. Then click on File --> export molecule --> Selection (sele) --> Save File name:  LeY-xray --> Files of type: pdb --> save. A [Ley-xray.pdb](https://github.com/glycodynamics/gly-docking/blob/main/LeY-xray.pdb) file will be saved in your computer 

Now select all the co-crystalized ligands and remove them (Action --> remove atoms). \
Then save protein: File --> export molecule --> Selection (1s3k) --> Save File name: receptor --> Files of type: pdb --> save. \
A [receptor.pdb](https://github.com/glycodynamics/gly-docking/blob/main/receptor.pdb) file will be saved to the current working directory of your computer.

These files have been prepared and placed under ./practice directory under your account in fucose. 

## 2. Connect to remote workstation:

Login to fucose using the instructions provided during the lecture. Linux/Mac users can use the terminal to connect to ccbrc workstation, whereas windows users should use PyTTY to connect to the CCBRC workstation.

```
$ ssh -X guestXX@fucose.pharmacy.olemiss.edu ## replace XX with your serial number

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
```
Now if you will type "ls -l" and hit enter, there should be two directories "practice and tutorial" available to everyone.
```
-bash-4.2$ ls -l
total 8
drwxr-xr-x. 2 sushil cgw 4096 Dec 15  2021 practice
drwxrwxr-x. 2 sushil cgw 4096 Dec 16  2021 tutorial
```

Directory "tutorial has all the precalculated data if you want to look into the correct output files. Another directory named "practice" has only input files for docking, and you can run calculations under this directory. To do so, type cd "./practice" and hit enter. Then type "ls -l," and it should list two directories:
```
flex_lig  : input files for rigid receptor + flexible ligand docking  
flex_rec  " input files for flexible receptor + flexible ligand docking  
```
To use any software in this machine you need to load them as a module. All the available software can be listed using the follwoing command:
```
$ module avail

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
## 3. Prepare flexible ligand and rigid receptor input files for docking
Change directory to 'flex_lig' and run the follwoing commands to perform docking:

```
module load autodock-vina 				# Load AutoDock Vina to your working environment
module load vina-carb/v1.2 				# Load Vina-carb to your working environment
module load mgltools/v2.1.5.7				# Load MGL Tools to your working environment

prepare_ligand4.py -l LeY-xray.pdb -A hydrogens 	# Convert LeY-xray.pdb into LeY-xray.pdbqt (Vina Input file)
prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A "hydrogens"	# Convert receptor.pdb in receptor.pdbqt (Vina Input file)
```
Nore: If you see an error message "**init.c(556):ERROR:161: Cannot initialize TCL**" run the following command to fix it
```
unset LD_LIBRARY_PATH

```

## 4. Perform docking using AutoDock Vina and vina-carb
In addition to receptor and ligand input files in PDBQT format, a configuration containing all the docking parameters is needed for docking. There are two configuration files _config_vc.txt_ (for vina_carb) and _config_vina.txt_ (for Vina) present in the practice directory. These configuration files have mainly following input parameters:

```
receptor=receptor.pdbqt		# name of receptor input file
ligand=LeY-xray.pdbqt		# name of the ligand input file
out=LeY-xray_vina_out.pdbqt	# name of the docking output file

center_x=-32.79          	# x, y and z coordinate of the search box
center_y=40.56          
center_z=14.89          
                         
size_x=30                	# Size of the search space box in Angstrom in x, y and z direction
size_y=30                
size_z=25                
                         
exhaustiveness=8   		# exhaustiveness of the global search (roughly proportional to time)
seed=0      			# explicit random seed - reproducibelity          
num_modes=10            	# number of docking solutions
energy_range=5          	# within enery raange 
cpu=8				# number of CPU cores used for docking
chi_coeff=1			# Chi coefficient energy (used in vina-carb only)              
chi_cutoff=2			# Chi cutoff energy (used in vina-carb only)
```


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
This run should finish within a minute and file "LeY-xray_vina_out.pdbqt" containing all the docked poses will be generated. Now repeat docking using Vina-carb and input parameters present in "config_vc.txt" file (as below). It should generate a file "LeY-xray_vc_out.pdbqt" that contains all the docking poses. 

**Perform docking using Vina-carb:**

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
## 5. Analyzing Docking Results
Since have docked LeY to Fab using two different software, AutoDock Vina, Vina-carb. These three programs differ in their approach to sample glycosidic linkage between monosaccharide units. Now visualize docking poses and see which program predicted LeY binding similar to what has been seen in the crystal structure of the complex obtained from PDB [PDB ID: 1S3K](https://www.rcsb.org/structure/1s3k). Copy _receptor.pdb_, _LeY-xray_vc_out.pdbqt_, and _LeY-xray_vina_out.pdbqt_ from fucose to your computer using WinSCP and load them in PyMOL.Now start PyMOL and load these files and analyze all docking poses by pressing the right arrow key of the keyboard. 


![alt text](https://github.com/glycodynamics/gly-docking/blob/main/images/Vina_docking_Fig.png)

You will see that the top scoring pose from AutoDock is very different from the crystal structure binding pose of LeY. However, the Vina-carb docking pose superimposes very much over the crystal structure binding pose. This shows that Chi intrinsic energy function introduced in vina-carb produces accurate results compared to Vina. However the 3rd docking pose of very accurate but the scoring function failed to rank this pose as the top ranking pose. Suppose if we did not know how LeY binds to Fab, we would never that docking results by vina are not correct. Therefore, for the cases where the docking protocol cannot be validated by reproducing the crystal structures of some ligands, it is advised to analyze many top scoring poses carefully instead of relying on the just top scoring pose. 


## 6. Flexible receptor docking 

Change directory to lex_rec and repeat the calculations for rigid receptor docking.
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

## 7. Prepare flexible ligand and rigid receptor input files for docking
In flexible receptor docking, you need to prepare an input file for the flexible and rigid parts of the receptor. We will use _prepare_flexreceptor4.py_ script of MGL Tools to prepare flexible receptor where _Tyr3 Tyr33 Tyr50 and Trp105_ are allowed to move during the docking. 

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
Please note that '-s' option in prepare_flexreceptor4.py script allows the user to provide residues that you want to keep flexible during the docking. In this case, we are keeping four amino-acid residues _Tyr3 Tyr33 Tyr50 and Trp105_ flexible during the docking. In general, one should keep only those residues to be flexible whose conformational changes can affect ligand binding. Output file _receptor_flex.pdbqt_ and _receptor_rigid.pdbqt_ contains flexible and rigid part of the receptor for docking.

## 8. Perform docking using AutoDock Vina and vina-carb

The configuration file for flexible receptor docking is nearly identical to rigid receptor docking except for the input file names. Since the receptor has been split into flexible and rigid parts, the configuration file should have the follwoing keywords

```
receptor=receptor_rigid.pdbqt	#name of file containing rigid part of the receptor
flex=receptor_flex.pdbqt	#name of file containing flexible part of the receptor

_REST ALL WILL REMAIN SAME_

```

Now run docking using Vina and vina-carb as earlier for rigid receptor docking:

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
Flexible ligand docking should take a longer time compared to rigid receptor docking. Repeat the same calculation using vina-carb and you should get finally two output files _LeY-xray_vina_flex_out.pdbqt_ and _LeY-xray_vc_flex_out.pdbqt_. Copy both these files to your local computer and analyze the binding modes. 


## 9. Analyzing flexible docking results

Copy output files _LeY-xray_vc_out.pdbqt_ and _LeY-xray_vina_out.pdbqt_ to your local computer and load them in PyMOL for visual inspection. You will see that Vina and vina-carb show similar behavior in flexible receptor docking and no significant change in the docking conformations has been observed. This shows that flexible receptor docking may not always be needed as rigid receptor docking can dock LeY well into Fab. 



## Useful resources:


https://www.click2drug.org/index.php#Binding%20free%20energy%20estimation 
