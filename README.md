# TTSCP_L102C_DSB-MBnThz-vs-TTSCP_L102C_DDSB-MBnThz

"This documentation represents the individual steps taken for the molecular dynamics (MD) simulation analysis for the following paper: XXX

In the MD analysis, our goal was to infer whether the presence/absence of a disulphide-bond (DSB) between 13th and 60th residues of TTSCP would yield different structural behaviour for the protein. Experimental studies showed that in the absence of DSB, the enzyme's catalytic performance dropped. It is often an observed phenomenon that removal of DSBs can cause conformational changes in the proteins, and we wanted to investigate whether such changes occur in our system and whether they can explain the lower catalytic performance of the DSB mutant. We performed this analysis under the presence of the MbnThz cofactor ligated to C102, which is the catalytic unit of the enzyme.

Firstly, we predicted the protein structures of TTSCP with both the DSB (C13, C60, and C102) and no DSB (A13, A60, C102) using AlphaFold2 (AF2). We will call these TTSCP_L102C and TTSCP_ΔDSB_L102C respectively from now. Protein structures can be predicted using ColabFold website (especially for small number of structures, one may find it easier to use a locally installed AF2 if the number of structures to process is large).  https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=KK7X9T44pWb7. ColabFold combines Google Colab and AF2 together and enables users to use a free-cloud based GPU providing system for protein structure prediction. In the experimental part of the study, the following protein sequences were used.

TTSCP_L102C

MELFTEAWAQAYCRKLNESEAYRKAASTWEGSLALAVRPDPKAGFPKGVAVVLDLWHGACRGAKAVEGEAEADFVIEADLATWQEVLEGRLEPLSALMRGLCELKKGTIAALAPYAQAAQELVKVAREVAENLYFQ

TTSCP_ΔDSB_L102C

MELFTEAWAQAYARKLNESEAYRKAASTWEGSLALAVRPDPKAGFPKGVAVVLDLWHGAARGAKAVEGEAEADFVIEADLATWQEVLEGRLEPLSALMRGLCELKKGTIAALAPYAQAAQELVKVAREVAENLYFQ

We can therefore put these sequences to the "query_sequence" section of ColabFold. We are defining our jobname as TTSCP_L102C and TTSCP_L102C_DDSB for their respective sequences. We would also like to perform AMBER relaxation on the predictions to achieve a more accurate geometry. For this, we are changing the “num_relax” variable to 5 to perform the relaxation on all of the predictions. We are not changing any other variable.

After AF2 job is finished, the files are ready to be downloaded. Results are usually saved in the contents folder, so you can click on the folder symbol on the left hand side of the page, and find the content folder with your job name. Your files would be saved in a .zip extension, so once you see your job_name.zip, you can click on the three dots on the right side of it and click “Download” to Download your files.

 
AF2 names the files based on the ranking, which is based on AF2's confidence in its prediction per-residue basis. We are choosing the structure file (with pdb extension) that has a label rank_001 in it for further processing. This file contains hydrogens added after the AF2 prediction, but we would like to add the hydrogens and determine the protonation states of residues via a dedicated software for a better accuracy. For this reason, we open our structure in a molecular visualization tool (such as PyMOL), and run the following command

remove hydro

This command removes the hydrogen atoms from the structure. We then save this molecule via File-Export Molecule, and name it as TTSCP_L102C_DDSB_nohydrogen.pdb. Then, we go to PDB2PQR website for the protonation of our molecule at a desired pH.  
https://server.poissonboltzmann.org/pdb2pqr. Here, we can choose “Upload a PDB file” option to load our file for protonation. We are changing the pH from 7.0 to 8.0 under the pKa options (This is the pH at which the catalysis reaction was carried out in the paper). We leave the “Use PROPKA to assign protonation states at provided pH” option as ticked. “For the Forcefield Options”, we are choosing AMBER and for the “Please choose an output naming scheme to use” we are choosing AMBER also. Then, we click “Start Job”. Within a couple of minutes, we see that we our job has been completed, and a “Use results with APBS” option appears on the right hand side of the screen. We click it, and click “Start Job” again from the page that is opened. Again, within a couple of minutes, an option called “View in 3Dmol” appears on the right side of the page. We click this, and a new page with our molecule is visible. We would like to download this molecule to our system, therefore we choose “Export as: PyMol” on the menu on the right side, and click Export. This downloads a .zip file with a folder, and several files within that folder. We open the file with the .pqr extension in PyMOL, and save it using  File-Export Molecule, where we name the file as TTSCP_L102C_DDSB_hydrogenadded.pdb. Importantly, we should manually delete the “HG” atom of the cysteine on the 102th position. This is because this hydrogen should not be there when the MbnThz cofactor is covalently attached to the cysteine. 

We now need to prepare our MBnThz cofactor. Such chemicals can be drawn with the help of online servers such as https://cactus.nci.nih.gov/translate/. In this, we can click on the “Start Structure Editor” button, which should open a new window with a JSME Structure Editor. We then right-click and find the “Paste MOL or SDF” on “SMILES” option. This opens a new small windows to enable us to put our chemical in its SMILES representation. We know that the SMILES representation of MBnThz is Cc2c(CCN1C(=O)CCC1=O)sc[n+]2Cc3ccccc3. However, we will need the cysteine residue that we are covalently attaching our MBnThz to , since we are going to align the cysteine (CYS) part of this molecule we create onto the CYS that is already present in the AF2 prediction. To obtain the CYS+MBnThz, we are using Cc2c(CCN1C(=O)CC(SC[C@@H](N)C(=O)O)C1=O)sc[n+]2Cc3ccccc3 SMILES code, and click “Accept”. It is important to note that one hydrogen in the thiazole ring is missing here, as the molecule is expected to have a +1 charge. This should show our cysteine+MBnThz together. We can then click “Submit Molecule”, and go back to https://cactus.nci.nih.gov/translate/. We are choosing "PDB" "Kekule representation" and "3D coordinates" for SMILES output format and click “Translate”. Once a new page opens that says “User-defined exchange format file: Click Here”, we click on “Click Here”. This downloads CYS+MBnThz together. We name this file MBnThz_CYS.pdb.
Next, we are going to put the cofactor into the protein. Since experimentally we demonstrated that cofactor is conjugated to cysteine on the 102th residue, we can ligate the cofactor into the protonated protein structure we obtained out of PDB2PQR server. To do this, we are going to take advantage of the “cysteine-like” part of the cofactor, and we are going to align this bit of the molecule onto the cysteine of the protein itself. 

For this purpose, we are first creating a .mol2 file that is generated from the .pdb after RESP charge determination. Our molecule after RESP charge determination has capping groups, so we should first remove these capping groups. We can manually choose the atoms for these capping groups and delete them in PyMOL, only to be left with the cysteine+cofactor. After that, we save this file as MBNThz_CYS.pdb by doing “File – Export Molecule”. Importantly, our molecule has a formal charge of +1, so we are taking that into consideration.

The cofactor molecule is first modified with an acetyl group at the N-terminus and a methylamide group at the C-terminus, in order to mimic the protein environment (presence of upstream and downstream amino acid  residues). The resulting SMILES string is: “O=C(NC)[C@@H](NC(C)=O)CSC1C(N(CCC2=C(C)[N+](CC3=CC=CC=C3)=CS2)C(C1)=O)=O”. This SMILES can then be used as input in the resp_charges_generation.ipynb Jupyter notebook.

In the Jupyter notebook, the cofactor molecule is built from the modified SMILES string using rdkit, which also adds hydrogen atoms. Based on the visualisation of the molecule, the atom indices are split into capping indices and backbone and side chain indices. This is essential as the capping acetyl and methylamide groups will need to be removed at the end in order to insert the residual cofactor into the protein.

A desired number of conformers (which can be modified with the n_confs variable) is then embedded into the cofactor molecule, and the conformers are aligned and energy minimised. The cofactor molecule is then passed to psiresp, along with the following options

    ​ Constraints: The acetyl and methylamide groups are constrained to have a charge sum of 0. This will ensure that the cofactor molecule will retain an integer charge once the capping groups are removed. The overall charge sum of the molecule is set to +1.
    ​ 
    ​ Geometry Optimisation: The geometry optimisation is set to True, and the level of theory is set to HF/6-31G* (in order to make the parameterised cofactor compatible with the Amber protein force field parameters).
    ​ 
    ​ Electrostatic Potential: The level of theory is set to HF/6-31G* (for the same reason as above).
The QM calculations (first geometry optimisation, followed by single point calculations) are run in the background with psi4. Once finished, the partial charges are output as a numpy array, with the array indices corresponding to the atom indices in the molecule visualisation. 

antechamber -fi pdb -fo mol2 -i MBNThz_CYS.pdb -o MBNThz_CYS.mol2 -pf y -c bcc -nc 1

This operation created a file with atom types and partial charges determined by AM1-BCC method. We are manually changing the charges in this mol2 file (by opening the file as a text editor) based on the output of the RESP calculation operation. 

Later, we are opening this mol2 file in PyMOL and then doing “File – Export Molecule”, and saving it as MBNThz_CYS_frommol2.pdb This step ensures avoiding problems while performing the alignment. It is also important that while we are aligning the cofactor onto the cysteine residue of the protein, the atoms that will be used in the alignment are named the same. For this reason, we are making sure that the atoms to be aligned to the cysteine in the protein are renamed as C, CA, CB, SG and O, N. This can be done manually by opening the MBNThz_CYS_frommol2.pdb file in a text editor. We should also change the atom types for these atoms in the  MBNThz_CYS.mol2 file and we should give them C, CX, CT, S, O and N respectively to have an entry like this in the file.

      1 N            8.033     -4.576      4.686  N        102 MBN       0.0508721369056092
      2 CA           7.606     -3.197      4.503  CX       102 MBN       0.05185889196762
      3 C            8.139     -2.346      5.656  C        102 MBN      -0.1199508359830749
      4 O            8.143     -2.769      6.811  O        102 MBN      -0.2697978250720955
      5 CB           6.078     -3.102      4.4140 CT       102 MBN       0.0832908019003133
      6 SG           5.455     -3.973      2.9480 S        102 MBN      -0.1527596148450697

We can then perform the alignment. To do this, we are opening both  MBNThz_CYS_frommol2.pdb and TTSCP_L102C_DDSB_hydrogenadded.pdb simultaneously in PyMOL, using “File – Open”. After opening these two files, we are individually selecting the atoms to be aligned (C, CA, CB, SG, O, N) in both the ligand and the protein itself by first switching “Selecting – Residues “ to “Selecting – Atoms” and then selecting individual atoms. As we click the atoms individually, PyMOL will automatically open a new selection as can be seen on the right side of the screen as (sele). We are renaming the ligand selections as lig, and protein selections as prot. Renaming can be achieved by clicking the “A” button near (sele), and then clicking rename selection and entering the desired names on the screen appearing. Then, in the command line of PyMOL we are running:

align lig, prot

After this command, our ligand is sticking out of the protein. Since the initial docking results showed that the cofactor should be in the protein’s core, we are manually changing the cofactor’s position. To do this, we are going to the “3-button editing mode” in PyMOL, and bringing the cofactor back into the protein core. Ideally, we are placing the cofactor such that there are not many major clashes with the cartoon representation of the protein, and also with not many atoms. Since we do not want to change anything from the protein structure itself, we are removing the C, CA, CB, SG, O, N from the ligand structure. It might be tricky to do it in the sticks representation, so we are first changing the lig into a line representation by clicking “S-lines”. We can then remove the atoms (also hydrogen atoms associated with them) by selecting all the atoms and “right click – remove”. We are saving this structure as “TTSCP_L102C_DDSB_hydrogenadded_conjugated.pdb”. Since we know  now that the cysteine at position 102 and our cofactor are ligated to each other, we can think of them as one big unnatural residue. Therefore, the CYS 102 entries when we open this pdb file in a text editor can be replaced with MBN 102! Similarly, the atoms related to MBN could be renamed MBN 102 if they are not so. We are also copy-pasting all the atoms from the actual MBN (without cysteine), beneath the cysteine atom entries, as this will be useful during renumbering, We are saving this file as “TTSCP_L102C_DDSB_hydrogenadded_conjugated_prerenumber.pdb”. Also, if the MBN atoms are logged in as HETATM in this file, we are changing them to ATOM, and we make sure that we keep all the coordinate columns aligned with each other and rest of the file. Then, using pdb4amber, we run the following command in our terminal:

pdb4amber -i TTSCP_L102C_DDSB_hydrogenadded_conjugated_prenumber.pdb -o TTSCP_L102C_DDSB_hydrogenadded_conjugated_renum.pdb

Now, we have all of our atoms in a numeric order (check one more time to make sure that none of the MBN atoms are returned to HETATM). 

Now, it is time for us to define this newly created non-proteinogenic MBN residue. To do this, we will need to create a .lib file for it. Also, we would like to create force field parameters for this residue. We are running the following bash scripts in our commands:

./generate_force.sh

This creates both the forcefield parameters for our new MBN residue, as well as it creates a lib file for it. In this lib file, we are checking if our partial charges are assigned correctly to each individual atom. We also need to make sure we have the identical coordinates in it as our pdb file. We can extract the coordinates directly from the .pdb file TTSCP_L102C_DDSB_hydrogenadded_conjugated_renum.pdb.

awk '/^ATOM/ && substr($0, 18, 3) == "MBN" { print substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8) }' TTSCP_L102C_DDSB_hydrogenadded_conjugated_renum.pdb > coordinates.pdb 

The coordinates.pdb file should then contain the coordinates of the MBN residue. We are copy-pasting this to the coordinates section of the .lib file. 

In the frcmod file we created, there is one parameter that is likely missing. We make sure that we are manually entering the following line into the ANGLE section in the frcmod file we created

C -N-hn    48.300     117.55

Our simulation preparation files are finally done! We can now generate both the topology and coordinate files that will be used by our simulations. For this, we are running the following command:

tleap -s -f tleap_preminimization.in > tleap_preminimization.out

We are checking the tleap_preminimization.out file manually using a pdb file. If you go to the end of the file, you should see 4 warnings in total and no errors. The warnings that appear in this way are likely fine. 

Note: The tleap_preminimization.in file for the TTSCP_L102C_DSB will differ from TTSCP_L102C_DDSB due to the presence of a disulphide bond. We should explicitly make this bond in tleap, if you inspect the provided .in file, you will see an extra line for it.

Now, we have the simulation files ready. We are first running the run_energy_minimization.py script to get a reasonable structure to start our simulation with. This should help us avoid having visually weird bonds after the simulations (they sometimes happen if the atoms are too close, PyMOL automatically makes them look like they are bonds and then they are visually bad).

python3 run_energy_minimization.pdb

We do not have to wait for this simulation to finish, and it does not have to be a long simulation. What matters is that the structure is relaxed enough such that the atoms that may make weird bonds with the cofactor are now away from it. This script gives us a .dcd output, which can be used to visually inspect the simulation over time. To do this, we first open TTSCP_L102C_DDSB_preenergyminimization.pdb in  PyMOL. Then, we also open the TTSCP_L102C_DDSB_energyminimization.dcd file in it. We are removing the water and sodium atoms and saving this file as TTSCP_L102C_DDSB_energyminimized.pdb. 

We also need to make sure that the MBN.lib file is updated with the new coordinates. To do this, we extract the coordinates from the pdb file we created 

awk '/^ATOM/ && substr($0, 18, 3) == "MBN" { print substr($0, 31, 8), substr($0, 39, 8), substr($0, 47, 8) }' TTSCP_L102C_DDSB_energyminimized.pdb > coordinates.pdb 

We then run the following command:

tleap -s -f tleap_postminimization.in > tleap_postminimization.out

We are checking the tleap_preminimization.out file manually using a pdb file. If you go to the end of the file, you should see 4 warnings in total and no errors. The warnings that appear in this way are likely fine. 

Note: The tleap_preminimization.in file for the TTSCP_L102C_DSB will differ from TTSCP_L102C_DDSB due to the presence of a disulphide bond. We should explicitly make this bond in tleap, if you inspect the provided .in file, you will see an extra line for it.

Now, we have the simulation files ready. We are running TTSCP_L102C_DSB_1.py up to 5, and TTSCP_L102C_DDSB_1.py up to 5 to run the simulations in parallel.
