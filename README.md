# CEPLA: Comparative Epitope Landscape Analysis
Pairwise comparisons of proteins based on their epitope landscape profiles

# How to install

The best way is to create an anaconda \citep{anaconda} environment that contains all specific collection of modules that the user have to install for running CEPLA.

Download and install Conda

    Introductions will be found on the website:
    https://conda.io/projects/conda/en/latest/user-guide/install/index.html

Create a Conda environment with Python 3.9 or higher

    conda create -n cepla python=3.9

Activate the Conda environment. The user will need to activate the Conda environment in each terminal in which the user want to use Cepla.

    conda activate cepla

Install require modules

    conda install -c anaconda pyqt
    conda install -c anaconda pathlib
    conda install -c anaconda biopython
    conda install -c conda-forge biotite
    conda install -c conda-forge matplotlib
    conda install -c anaconda numpy

Download CEPLA script from GitHub


# Usage

Go into the folder where the *Gui.py* script is located and run the main program.
The Graphical User Interface will open in a seperate window in which the user can select the relevant data by using the *browse-button*.

Please note that for the reference and for the query sequence the user have to choose a fasta file format as an input.
For the T-cell epitope data the user have to choose the path to the epitope folder the user generated with the IEDB tool before. The exact description of how the epitope data should be saved will be described in the next chapter.
The B-cell epitopes are also previously generated with either the Bepipred2.0 tool from IEDB or Epidope. 
For the epitope analyse a csv file has to be uploaded. When using Epidope B-cell epitopes, the file *epidope\scores.csv* have to be used. Bepipred2.0 generates just one csv-file.
Unfortunately, the program is error-prone here. Please note that the epitopes are stored in the individual folders as described in the next section. Only then the analysis can work. 
It must also be ensured that the individual epitope data have the same file name as the sequence file. 

There exist mandatory fields and extended fields.


| mandatory field | what to do | datatyp |
| -------|------|------|
| reference sequence | choose path to sequence | fasta |
| query sequence | choose path to sequence | fasta |
| reference T-cell epitope data | choose path to folder | folder with subfolder |
| query T-cell epitope data | choose path to folder | folder with subfolder |
| Job name | give your analysis a name | - |
| mhc class | choose mhcI, mhcII or both | - |
       
   

 | extended field | what to do | datatyp |
 |-------|-----|-----|
 | domain architecture | choose path to domainfile | txt |
 | population coverage | choose one of the options | - |
 | show plots | typ in the plotnumber | - |
 | outputpath | typ in outputpath | - |
 | additional T-cell epitope data  | choose path to folder | folder with subfolder |
 | reference B-cell epitope data | choose path to file | csv |
 | query B-cell epitope data | choose path to file | csv |
 


The extended fields have to be confirmed with the *use-button*.
The analysis starts by clicking on *Go-button* and a new window will be open.
Screen 2 shows in an output window an overview of the selected data.

When the analysis is finished, the alignment is displayed in the lower window, which was generated with muscle3 version 3.7.
With a click on *show alignment-button* the alignment will be pop up in an extra window. Certainly the process takes a moment. 
Clicking on *show colored alignment–button* the alignment with colored amino acid regions will be shown. The intensity of the color indicates the similarity of the sequences.
The *output-button* shows the result of the analysis in a new window and the profile can be simply saved by clicking on the save-symbol left in the top bar.
The result is saved as a svg in the output folder automatically after generating the epitope landscape profile.

If all mandatory fields and extended fields are selected, this is the presentation of the results in different subplots in the appropriate order. The numbering starts with 0.

* Subplot 0 and 1: predicted T-cell epitopes (created with IEDB tool)
* Subplot 2: epitopescore (1 - percentile value)
* Subplot 3: difference plot (shows changes between reference and query)
* Subplot 4: binding numbers
* Subplot 5: identity plot (shows mutations (black) and gaps (green))
* Subplot 6: additional T-cell epitope data (created with IEDB tool)
* Subplot 7: population coverage
* Subplot 8 and 9: predicted B-cell epitopes (created with IEDB)
* Subplot 10: domain architecture

If only certain subplots are to be displayed, these are also the numbers of the subplots that must be specified in the extended field *show plots*.
MCHI and MHCII are displayed in separate windows if both results were desired.
There are two sliders on the right side. With the upper slider the user can set the epitope score of the epitopes to filter them. They will be then just visualized in subplot 0 as colored bars.
The default value is 0.9 and shows the epitopes with the highest binding affinity (0.9-1.0).
With the slider below the slider can also set the epitope score for epitopes. However, the filtered epitopes then only refer to subplots 1, 4, 7.
The default value is 0, so every epitope with an epitope score until 1.0 will be display.
These two sliders in combination allow viewing of all epitopes up to an epitope score of 1.0 and filtering and visualizing epitopes with high binding affinity.





# Prediction of T-cell epitopes

To predict linear T-cell epitopes, the IEDB's T-cell epitope – [MHC binding prediction tool v2.24](https://downloads.iedb.org/tools/) was used. It is a consensus approach and predicts IC50 values for peptides that bind to specific MHC molecules. There is a website and a download version. In the master's thesis, the download version was used and the epitopes were predicted using the script. Linux was used as the operating system and the scripts *pred\binding.py* for MHC I and *mhc\II\binding.py* for MHC II were called in an Anaconda environment with Python 2.7.15 installed.
In general, the tool consists of two subprograms: Peptide Binding to MHC Class I Molecules and Peptide Binding to MHC Class II Molecules. Both scripts must be downloaded and used separately.

Shell command for MHC-I binding predictions: 

    python predict_binding.py IEDB_recommended HLA-A*01:01,HLA-A*02:01,HLA-A*02:03,HLA-A*02:06,
    HLA-A*03:01,HLA-A*11:01,HLA-A*23:01,HLA-A*24:02,HLA-A*26:01,
    HLA-A*30:01,HLA-A*30:02,HLA-A*31:01,HLA-A*32:01,HLA-A*33:01,HLA-A*68:01,
    HLA-A*68:02,HLA-B*07:02,HLA-B*08:01,HLA-B*15:01,HLA-B*35:01,HLA-B*40:01,
    HLA-B*44:02,HLA-B*44:03,HLA-B*51:01,HLA-B*53:01,HLA-B*57:01,HLA-B*58:01
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8 sequence.fa > output.out;

The exact code call was repeated accordingly for the other peptide lengths: 9, 10 and 11.

Shell command for MHC-II binding predictions: 

    python mhc_II_binding.py IEDB_recommended HLA-DRB1*01:01,HLA-DRB1*03:01,
    HLA-DRB1*04:01,HLA-DRB1*04:05,HLA-DRB1*07:01,HLA-DRB1*08:02,HLA-DRB1*09:01,
    HLA-DRB1*11:01,HLA-DRB1*12:01,HLA-DRB1*13:02,HLA-DRB1*15:01,HLA-DRB3*01:01,
    HLA-DRB3*02:02,HLA-DRB4*01:01,HLA-DRB5*01:01,HLA-DQA1*05:01/DQB1*02:01,
    HLA-DQA1*05:01/DQB1*03:01,HLA-DQA1*03:01/DQB1*03:02,HLA-DQA1*04:01/DQB1*04:02,
    HLA-DQA1*01:01/DQB1*05:01,HLA-DQA1*01:02/DQB1*06:02,HLA-DPA1*02:01/DPB1*01:01,
    HLA-DPA1*01:03/DPB1*02:01,HLA-DPA1*01:03/DPB1*04:01,HLA-DPA1*03:01/DPB1*04:02,
    HLA-DPA1*02:01/DPB1*05:01,HLA-DPA1*02:01/DPB1*14:01 sequence.fa > output.out;
    
    
    
  
The MHCI/MHCII predicting tool, was used to predict linear T-cell epitopes. The analysis has now been specially programmed for this tool.
Only epitopes of a certain length can bind to MHCI. Therefore, different epitope lengths (8, 9, 10, 11 amino acids) are predicted with recommended alleles.
To use the Comparative Epitope Landscape Analysis, the predicted epitope data must be saved in appropriate folders. For example, epitopes with length 8 go in folder 8 and epitopes with length 9 in folder 9.
Epitopes can bind to MHCII regardless of their length. The results are saved in a folder named 2.



# Prediction of B-cell epitopes 

The online *Antibody Epitope prediction tool*  from the IEDB website was used to generate B-cell epitopes. It is a method for predicting continuous antibody epitopes from protein sequences. This tool offers various operations to predict the localization of continuous epitopes on parameters such as hydrophilicity, flexibility, accessibility, whorls, exposed surface, polarity and antigenic propensity of polypeptide chains.
In this work, *BepiPred-2.0* was selected from the various offered prediction methods. It is a sequential B-cell epitope predictor. The tool is using a random forest algorithm to predict B-cell epitopes from a protein structure. This algorithm is trained on epitopes and non-epitopes amino acids obtained from crystal structures.
The tool creates a csv-file as an output file, which can then be used in CEPLA.
This file contains 4 columns: Position, Residue, Score and Assignmnet.
For each amino acid in the protein sequence, a antigenicity scores is calculated that indicates whether or not it might be part of an epitope. 
The antigenicity score describe the ability to bind with a specific antibody.
A threshold value is formulated for this purpose (default value 0.5). If the value is above the threshold, it is predicted to be part of an epitope and gets an 'E' in the assignment-column of the output file.

The *EpiDope* software offers another prediction possibility. "EpiDope has been shown to be the best-performing among currently available B-cell epitope prediction tools." [collatz2021epidope] For this reason, we wanted to give the user the opportunity to also display epitope profiles of B-cells predicted with this tool.
The Software uses a deep neural network to predict linear B-cell epitopes and is supported on Linux and Mac. A detailed installation specification and condition guide can be found on GitHub.

The only input file EpiDope needs is a Fasta file of the protein sequence to be scanned for epitopes.
EpiDope produces different output csf-files. The file *epidope\scores.csv* lists the predicted score per amino acid in a similar way like *BepiPred-2.0* and *predicted\epitopes.csv* lists all regions with a antigenicity scores higher than the defined threshold (>0.8) and the last output *predicted\epitopes\sliced.faa* is a multi-fasta file of potential epitopes.

Shell command for B-cell epitope prediction with EpiDope: 

    epidope −i /path_to/fasta.fa −o /path_to/results/
    
* -i: Path to the protein sequence in fasta format
* -o: Output path for the result datas
