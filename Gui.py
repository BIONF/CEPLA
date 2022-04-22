from PyQt5.QtWidgets import QApplication, QDialog, QFileDialog, QMessageBox,QLineEdit
from PyQt5 import QtWidgets
from PyQt5.uic import loadUi
from pathlib import Path
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import biotite.sequence as bseq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics
from Bio import SeqIO
import sys, os, shutil


# Module:
import epitope_editing
import align_process
import subplot_creater





class MainWindow(QDialog):

    """
    Class for Screen 1 that provides all the functionality behind the button s and calls the main method to start the analys.

    """
    def __init__(self):

        """
        Loads the file gui2.ui to display the graphical user interface and calls the button functions.

        :param paths: list included: refSeq, specSeq, epitopeRef, epitopeSp, domain, addfile, epitopeRefB, epitopeSpB
        :param lines: list included: popCov, showplots, outputpath
        :param comboBox1: Scroll down menu for population coverage
        :param Browse_*number*: Browse-button for different selections.
        """
        super(MainWindow, self).__init__()
        loadUi("gui2.ui",self)
        self.paths = ["","","","","","","",""]
        self.lines = ["","",""]


        self.comboBox1.addItem("No") 
        self.comboBox1.addItem("World")
        self.comboBox1.addItem("East Asia")
        self.comboBox1.addItem("Northeast Asia")
        self.comboBox1.addItem("South Asia")
        self.comboBox1.addItem("Southeast Asia")
        self.comboBox1.addItem("Southwest Asia")
        self.comboBox1.addItem("Europe")
        self.comboBox1.addItem("East Africa")
        self.comboBox1.addItem("West Africa")
        self.comboBox1.addItem("Central Africa")
        self.comboBox1.addItem("North Africa")
        self.comboBox1.addItem("South Africa")
        self.comboBox1.addItem("West Indies")
        self.comboBox1.addItem("North America")
        self.comboBox1.addItem("Central Amercia")
        self.comboBox1.addItem("South America")
        self.comboBox1.addItem("Oceania")



        self.Browse1.clicked.connect(lambda: self.browsefiles(0)) # reference sequence
        self.Browse_2.clicked.connect(lambda: self.browsefiles(1)) # query sequence
        self.Browse_3.clicked.connect(lambda: self.browsefolder(2)) # reference epitope data (T-cell)
        self.Browse_4.clicked.connect(lambda: self.browsefolder(3)) # query epitope data (T-cell)
        self.Browse_5.clicked.connect(lambda: self.browsefiles(4)) # domain architecture
        self.Browse_6.clicked.connect(lambda: self.browsefiles(5)) # additional epitope data
        self.Browse_7.clicked.connect(lambda: self.browsefiles(6)) # reference epitope data (B-cell)
        self.Browse_8.clicked.connect(lambda: self.browsefiles(7)) # query epitope data (B-cell)
        self.useBut.clicked.connect(self.clickMethod)
        self.goBut.clicked.connect(self.goButtonFunc)

    def goButtonFunc(self):

        """
        Function behind the Go button. Starts the analysis.

        """
        print(type(self.filenameSp),self.filenameSp.text())
        if self.filenameSp.text() != '' or self.filenameRef.text() != '' or self.epitopedataSp.text() != '' or self.epitopedataRef.text() != '': # check if the User selected the data
            if not os.path.exists(self.outputP.text() + "/" + self.jobname.text()): # check if the jobname already exist
                self.startProgram(self.lines[2],self.paths[2],self.paths[3],self.paths[1],self.paths[0],self.lines[0],self.paths[5],self.lines[1],self.paths[4], self.jobname.text(), self.paths[6], self.paths[7], self.refname.text(), self.queryname.text())
                self.gotoScreen2()
            else:
                self.show_popup()
        else:
            self.show_popup()



    def startProgram(self,outputpath,epitopeRef,epitopeSp,specSeq,refSeq,popCov,addfile,subplot,domain,jobname, epitopeRefB, epitopeSpB, refName, quName):

        """
        It is checked which mhc class it is and the variables are provided for the analysis.
        In addition, it is checked whether the job name already exists. If so, a pop up warning window will open.

        :param outputpath: Path to output folder
        :param epitopeRef: Path to T-cell epitopes of the reference
        :param epitopeSp: Path to T-cell epitopes of the query
        :param specSeq: Path to query sequence
        :param refSeq: Path to reference sequence
        :param popCov: Defined population or area
        :param addfile: Path of extra file
        :param subplot: Plotnumber for the subplots
        :param domain: Path to domain file
        :param jobname: Name of the Job for the output folder
        :param epitopeRefB: Path to B-cell epitopes of the reference
        :param epitopeSpB: Path to B-cell epitopes of the query
        :param refName: Name of Reference for plottitel or labels
        :param quName: Name of query for plottitel or labels
        :param mhc: MHC class

        :return: listResult with all parameters for starting the analysis
        """


        mhc="both"
        self.specSeq = specSeq
        self.refSeq = refSeq
        self.epitopeRef = epitopeRef
        self.epitopeSp = epitopeSp
        self.popCov = popCov
        self.addfile = addfile
        self.subplot = subplot
        self.domain = domain
        self.jobname = jobname
        self.refName = refName
        self.quName = quName
        self.outputpath = outputpath
        self.epitopeRefB = epitopeRefB
        self.epitopeSpB = epitopeSpB
        subplot = "all"


        if self.I.isChecked():
            mhc = 'I'
        if self.II.isChecked():
            mhc = 'II'
        if self.both.isChecked():
            mhc= self.both.text()

        self.listResult = epitopeanalyse(outputpath,epitopeSp,epitopeRef,specSeq,refSeq,popCov,mhc,addfile,subplot,domain,jobname,epitopeRefB,epitopeSpB, refName, quName)



    def show_popup(self):

        """
        Pop up warning window for still existing jobnames.

        """
        msg = QMessageBox()
        msg.setWindowTitle("Error job name")
        msg.setText('Falsch!!!!!!!')

        x = msg.exec_()

    def gotoScreen2(self):

        """
        Switches to Screen 2

        """
        widget.setCurrentIndex(widget.currentIndex()+1)

    def clickMethod(self):

        """
        If the Use-button is pressed, the advanced settings are applied.

        """

        self.lines[0] = self.comboBox1.currentText()
        self.lines[1] = self.showplots.text()
        self.lines[2] = self.outputP.text()
        

    def browsefiles(self,n):

        """
        Function behind the browse button. Opens the home directory from which you can then search for the relevant files.

        :param fname: Path to file

        """
        fname =QFileDialog.getOpenFileName(self, 'Open file', r'C:\Users\stefa\Desktop')[0]
        self.paths[n] =fname
        self.updateTexts()


    def browsefolder(self,n):

        """
        Function behind the browse button. Opens the home directory from which you can then search for the relevant folder.

        :param fname:  Path to folder

        """
        fname =str(QFileDialog.getExistingDirectory(self, 'Open file', r'C:\Users\stefa\Desktop'))
        self.paths[n] =fname
        self.updateTexts()

    def updateTexts(self):

        """
        Updates the text in the input boxes

        """
        self.filenameRef.setText(self.paths[0])
        self.filenameSp.setText(self.paths[1])
        self.epitopedataRef.setText(self.paths[2])
        self.epitopedataSp.setText(self.paths[3])
        self.filenameDom.setText(self.paths[4])
        self.addData.setText(self.paths[5])
        self.refBcell.setText(self.paths[6])
        self.quBcell.setText(self.paths[7])



#______________________________________________________________________________________________________________________#




class Screen2(QDialog,QtWidgets.QGraphicsView):
    def __init__(self):
        super(Screen2, self).__init__()
        loadUi("gui3.1.ui",self)
        self.le = QLineEdit()

        self.alignBut.clicked.connect(self.openFile)
        self.colorBut.clicked.connect(self.showAlignment2)
        self.outBut.clicked.connect(self.showPlot)
        self.percenSlider.valueChanged.connect(self.v_change)
        self.percenSlider_2.valueChanged.connect(self.epiV_change)



    def openFile(self):

        """
        Opens the alignment file of the reference and query in a texteditor.

        """
        path = mainwindow.outputpath + "/" + mainwindow.jobname+"/"+mainwindow.jobname+"_aligningData/aligned/TestSpecie/" + os.path.splitext(os.path.basename(mainwindow.specSeq))[0] + "_" + os.path.splitext(os.path.basename(mainwindow.refSeq))[0] + "_alignedseqs.txt"
        os.system('chmod u+x '+ path)
        os.system('gedit ' + path)


    def showAlignment(self):

        """
        Load the alignment file of the reference and query and shows it formatted on Screen 2.

        """
        alignment = []
        path = mainwindow.outputpath + "/"+ mainwindow.jobname+"/"+mainwindow.jobname+"_aligningData/aligned/TestSpecie/"
        files1 = os.listdir(path)
        for i in files1:
            fullpath = os.path.join(path, i)
        record = list(SeqIO.parse(fullpath,'fasta'))
        seq1 = record[0].seq
        seq2 = record[1].seq
        alignment.append([seq1,seq2,0,len(seq1)])
        if screen2.alignField.toPlainText() == '':
            screen2.alignField.append(format_alignment(*alignment[0]))


    def showAlignment2(self):

        """
        Creates BLOSUM62 matrix and perform pairwise sequence alignment with affine gap penalty. Terminal gaps are not penalized.
        Draw first and only alignment and the color intensity indicates the similiarity.
        The colored alignment is saved as png in the output folder.

        """
        path = mainwindow.outputpath + "/" +mainwindow.jobname+"/"+mainwindow.jobname+"_aligningData/Aligndata/sequences_aligned/" + os.path.splitext(os.path.basename(mainwindow.specSeq))[0] + "_" + os.path.splitext(os.path.basename(mainwindow.refSeq))[0] + ".txt"

        splited_list = []
        
        splited_list.append(mainwindow.refName)
        splited_list.append(mainwindow.quName)

        record = list(SeqIO.parse(path, 'fasta'))
        seq1 = record[0].seq
        seq2 = record[1].seq

        seq3 = bseq.ProteinSequence(seq1)
        seq4 = bseq.ProteinSequence(seq2)

        matrix = align.SubstitutionMatrix.std_protein_matrix()

        alignments = align.align_optimal(seq3, seq4, matrix, gap_penalty=(-10, -1), terminal_penalty=False)

        self.fig = plt.figure(figsize=(42, 16))
        self.ax = self.fig.add_subplot(111)
        graphics.plot_alignment_similarity_based(
              self.ax, alignments[0], matrix=matrix, labels=[splited_list[0], splited_list[1]],
              show_numbers=True, show_line_position=True, spacing=1)

        self.fig.tight_layout()
        plt.savefig(mainwindow.outputpath + "/" +mainwindow.jobname+"/"+mainwindow.jobname+"_aligningData/"+ mainwindow.refName + '_'+ mainwindow.quName + "_coloredAlignment.png")  # figure will be saved here


    def showPlot(self):

        """
        Plots will be generated by using the current epitope score.

        """
        currentPercentile(self.percenSlider.value(), self.percenSlider_2.value())


    def v_change(self):

        """
        Epitope score slider for subplot 3, 5, 7

        """
        my_value = self.percenSlider.value()
        scaledValue = float(my_value)/100
        self.valueBox.setText(str(scaledValue))
        

    def epiV_change(self):

        """
        Epitope score slider for subplot 1, 2

        """
        my_value = self.percenSlider_2.value()
        scaledValue = (float(my_value)/100)
        self.lineEdit.setText(str(scaledValue))



#______________________________________________________________________________________________________________________#



def epitopeanalyse(outputpath,epitopeSp,epitopeRef,specSeq,refSeq,popCov,mhc, addfile,subplot,domain,jobname, epitopeRefB, epitopeSpB, refName, quName):

    """
    This is the main program. The folders and subfolders are generated in the output path and then the analysis is carried out with several subprograms.
    The input list contains the MHC class. The analysis is carried out separately for MHC class I and II.

    :param outputpath: Path to output folder
    :param epitopeRef: Path to T-cell epitopes of the reference
    :param epitopeSp: Path to T-cell epitopes of the query
    :param specSeq: Path to query sequence
    :param refSeq: Path to reference sequence
    :param popCov: Defined population or area
    :param addfile: Path of extra file
    :param subplot: Plotnumber for the subplots
    :param domain: Path to domain file
    :param jobname: Name of the Job for the output folder
    :param epitopeRefB: Path to B-cell epitopes of the reference
    :param epitopeSpB: Path to B-cell epitopes of the query
    :param refName: Name of Reference for plottitel or labels
    :param quName: Name of query for plottitel or labels
    :param mhc: MHC class
    :param epivalue: selected epitope score from the second slider


    """

    if domain == '':
        answ = 'No'
    else:
        answ = 'Yes'



    # Output overview will be printed in the textbrowser at Screen 2

    screen2.textBrowser.append('Overview:')
    screen2.textBrowser.append('Output file is '+ outputpath)
    screen2.textBrowser.append('Query epitope file is '+ epitopeSp)
    screen2.textBrowser.append('Reference epitope file is '+ epitopeRef)
    screen2.textBrowser.append('Query sequence is '+ specSeq)
    screen2.textBrowser.append('Reference sequence is '+ refSeq)
    screen2.textBrowser.append('addfile is '+ addfile)
    screen2.textBrowser.append("population selected: "+ popCov)
    screen2.textBrowser.append('Class is '+ mhc)
    screen2.textBrowser.append('Subplots will be print '+ subplot)
    screen2.textBrowser.append('Domain architecture will be print: '+ answ)
    screen2.textBrowser.append('Reference B-cell epitope data: ' + epitopeRefB)
    screen2.textBrowser.append('Query B-cell epitope data: ' + epitopeSpB)
    screen2.textBrowser.append('_________________________________________________________')
    screen2.textBrowser.append('')
    screen2.textBrowser.append('The calculation may take a while ...')


    # ------------------------------------------------------------------------------------------------------------------#


    # Folder and subfolder will be generated in the output path

    folderlist=[2,8,9,10,11]
    for i in folderlist:
        os.makedirs(outputpath + "/"+jobname+"/seq1/filtered/"+str(i))
        os.makedirs(outputpath + "/"+jobname+"/seq2/filtered/"+str(i))

    os.makedirs(outputpath + "/"+jobname+"/seq1/results/mhc_i")
    os.makedirs(outputpath + "/"+jobname+"/seq1/results/mhc_ii")
    os.makedirs(outputpath + "/"+jobname+"/seq1/results/modified")

    os.makedirs(outputpath + "/"+jobname+"/seq2/results/mhc_i")
    os.makedirs(outputpath + "/"+jobname+"/seq2/results/mhc_ii")
    os.makedirs(outputpath + "/"+jobname+"/seq2/results/modified")

    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/Aligndata/sequences_aligned")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/aligned/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_i/results_seqs_ALL/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_i/filtered_results_seqs_ALL/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_i/result_for_plot/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_ii/results_seqs_ALL/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_ii/filtered_results_seqs_ALL/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_ii/result_for_plot/TestSpecie")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/seq")
    os.makedirs(outputpath + "/"+jobname+"/"+jobname+"_aligningData/svg-Data")

    # ------------------------------------------------------------------------------------------------------------------#


    inputList = []
    epivalue = screen2.percenSlider_2.value()

    if mhc == 'both':
        inputList.append('I')
        inputList.append('II')
    else:
        inputList.append(mhc)


    for hla in inputList:

        # scripts for analysis 1:
        epitope_editing.crop_annotation(epitopeSp, outputpath + "/"+jobname+"/seq1/filtered/", specSeq, "mhc_"+hla.lower())
        epitope_editing.crop_annotation(epitopeRef, outputpath + "/"+jobname+"/seq2/filtered/", refSeq, "mhc_"+hla.lower())
        epitope_editing.combine_epitopes(outputpath + "/"+jobname+"/seq1/filtered", outputpath + "/"+jobname+"/seq1/megafile.txt")
        epitope_editing.combine_epitopes(outputpath + "/" + jobname + "/seq2/filtered", outputpath + "/" + jobname + "/seq2/megafile.txt")
        epitope_editing.get_epitopes(outputpath + "/"+jobname+"/seq1/megafile.txt", outputpath + "/"+jobname+"/seq1/results")
        epitope_editing.get_epitopes(outputpath + "/"+jobname+"/seq2/megafile.txt", outputpath + "/"+jobname+"/seq2/results")

        screen2.textBrowser.append("Analysis 1 was successful")


        # scripts for analysis 2:
        if popCov != "No":
            os.system("python popCov/calculate_population_coverage.py -p " + popCov +" -c "+ hla +" -f " + outputpath + "/"+jobname+"/seq1/results")
            os.system("python popCov/calculate_population_coverage.py -p " + popCov + " -c " + hla + " -f " + outputpath + "/"+jobname+"/seq2/results")

            screen2.textBrowser.append("Population coverage was calculated successfully")
        else:
            dict1 =[]
            dict2 =[]

            with open(outputpath + "/"+jobname+"/seq1/results/"+Path(specSeq).stem+".txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.replace('\n','')
                    line += '\t0.0\n'
                    dict1.append(line)
            with open(outputpath + "/"+jobname+"/seq1/results/modified/"+Path(specSeq).stem+".txt", "w") as o:
                for i in dict1:
                    o.write(i)

            with open(outputpath + "/"+jobname+"/seq2/results/"+Path(refSeq).stem+".txt", "r") as g:
                lines2 = g.readlines()
                for line2 in lines2:
                    line2 = line2.replace('\n', '')
                    line2 += '\t0.0\n'
                    dict2.append(line2)
            with open(outputpath + "/"+jobname+"/seq2/results/modified/"+Path(refSeq).stem+".txt", "w") as a:
                for j in dict2:
                    a.write(j)

        shutil.copy(specSeq, outputpath + "/"+jobname+"/"+jobname+"_aligningData/seq")
        shutil.copy(refSeq, outputpath + "/"+jobname+"/"+jobname+"_aligningData/seq")

        destFile = outputpath + "/"+jobname+"/"+jobname+"_aligningData/seq/" + os.path.splitext(os.path.basename(specSeq))[
            0] + "_" + os.path.splitext(os.path.basename(refSeq))[0] + ".txt"
        shutil.copyfile(specSeq, destFile)

        with open(destFile, 'a') as resultfile:
            resultfile.write("\n\n")
            with open(refSeq, 'r') as readfile:
                for line in readfile:
                    resultfile.write(line)

        shutil.copy(destFile, outputpath + "/"+jobname+"/"+jobname+"_aligningData/Aligndata/sequences_aligned")
        align_process.seqfilter(outputpath + "/"+jobname+"/"+jobname+"_aligningData/Aligndata/sequences_aligned", outputpath,jobname)
        screen2.textBrowser.append("The two sequences were successfully aligned")

        align_process.startaligner(outputpath + "/"+jobname+"/"+jobname+"_aligningData/aligned/TestSpecie", outputpath + "/"+jobname+"/seq1/results/modified", outputpath + "/"+jobname+"/seq2/results/modified", outputpath, hla, jobname, screen2.percenSlider.value(), screen2.percenSlider_2.value())
        align_process.filterGaps(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_"+hla.lower()+"/results_seqs_ALL", outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_"+hla.lower()+"/filtered_results_seqs_ALL", specSeq, refSeq,screen2.percenSlider.value(), screen2.percenSlider_2.value())
        align_process.median(outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_"+hla.lower()+"/filtered_results_seqs_ALL", outputpath + "/"+jobname+"/"+jobname+"_aligningData/mhc_"+hla.lower()+"/result_for_plot", specSeq, refSeq,screen2.percenSlider.value(),screen2.percenSlider_2.value())
        screen2.textBrowser.append("The files were processed based on references")

        screen2.showAlignment()
    return inputList



def format_alignment(align1, align2, begin, end, full_sequences=False):

    """Format the alignment prettily into a string.

    IMPORTANT: Gap symbol must be "-" (or ['-'] for lists)!

    Since Biopython 1.71 identical matches are shown with a pipe
    character, mismatches as a dot, and gaps as a space.

    Prior releases just used the pipe character to indicate the
    aligned region (matches, mismatches and gaps).

    Also, in local alignments, if the alignment does not include
    the whole sequences, now only the aligned part is shown,
    together with the start positions of the aligned subsequences.
    The start positions are 1-based; so start position n is the
    n-th base/amino acid in the *un-aligned* sequence.

    NOTE: This is different to the alignment's begin/end values,
    which give the Python indices (0-based) of the bases/amino acids
    in the *aligned* sequences.

    If you want to restore the 'historic' behaviour, that means
    displaying the whole sequences (including the non-aligned parts),
    use ``full_sequences=True``. In this case, the non-aligned leading
    and trailing parts are also indicated by spaces in the match-line.


    Note: This function was taken from Biopython.pairwise2.format_alignment and modified and adapted for this program.

    """

    align_begin = begin
    align_end = end
    start1 = start2 = ""
    start_m = begin  # Begin of match line (how many spaces to include)
    # For local alignments:
    if not full_sequences and (begin != 0 or end != len(align1)):
        # Calculate the actual start positions in the un-aligned sequences
        # This will only work if the gap symbol is '-' or ['-']!
        start1 = str(len(align1[:begin]) - align1[:begin].count("-") + 1) + " "
        start2 = str(len(align2[:begin]) - align2[:begin].count("-") + 1) + " "
        start_m = max(len(start1), len(start2))
    elif full_sequences:
        start_m = 0
        begin = 0
        end = len(align1)

    if isinstance(align1, list):
        # List elements will be separated by spaces, since they can be
        # of different lengths
        align1 = [a + " " for a in align1]
        align2 = [a + " " for a in align2]

    s1_line = ["{:>{width}}".format(start1, width=start_m)]  # seq1 line
    m_line = [" " * start_m]  # match line
    s2_line = ["{:>{width}}".format(start2, width=start_m)]  # seq2 line

    for n, (a, b) in enumerate(zip(align1[begin:end], align2[begin:end])):
        # Since list elements can be of different length, we center them,
        # using the maximum length of the two compared elements as width
        m_len = max(len(a), len(b))
        s1_line.append("{:^{width}}".format(a, width=m_len))
        s2_line.append("{:^{width}}".format(b, width=m_len))
        if full_sequences and (n < align_begin or n >= align_end):
            m_line.append("{:^{width}}".format(" ", width=m_len))  # space
            continue
        if a == b:
            m_line.append("{:^{width}}".format("|", width=m_len))  # match
        elif a.strip() == "-" or b.strip() == "-":
            m_line.append("{:^{width}}".format(" ", width=m_len))  # gap
        else:
            m_line.append("{:^{width}}".format(".", width=m_len))  # mismatch

    return "\n".join(["".join(s1_line), "".join(m_line), "".join(s2_line)])


def currentPercentile(panelValue, panelValue_2):

    """
    If one of the two epitope score pannel is changed, Analysis2 is called up again for the calculation and filtered according to the corresponding epitopes.
    At the end the result will be shown.

    :param panelValue: epitope score for subplot 3,5,7
    :param panelValue_2: epitope score for subplot 1,2

    """

    for hla in mainwindow.listResult:
        align_process.startaligner(mainwindow.outputpath + "/" + mainwindow.jobname + "/" + mainwindow.jobname + "_aligningData/aligned/TestSpecie",
                                    mainwindow.outputpath + "/" + mainwindow.jobname + "/seq1/results/modified",
                                    mainwindow.outputpath + "/" + mainwindow.jobname + "/seq2/results/modified", mainwindow.outputpath, hla, mainwindow.jobname,panelValue, panelValue_2)

        align_process.filterGaps(
            mainwindow.outputpath + "/" + mainwindow.jobname + "/" + mainwindow.jobname + "_aligningData/mhc_" + hla.lower() + "/results_seqs_ALL",
            mainwindow.outputpath + "/" + mainwindow.jobname + "/" + mainwindow.jobname + "_aligningData/mhc_" + hla.lower() + "/filtered_results_seqs_ALL",
            mainwindow.specSeq, mainwindow.refSeq, panelValue, panelValue_2)

        align_process.median(
            mainwindow.outputpath + "/" + mainwindow.jobname + "/" + mainwindow.jobname + "_aligningData/mhc_" + hla.lower() + "/filtered_results_seqs_ALL",
            mainwindow.outputpath + "/" + mainwindow.jobname + "/" + mainwindow.jobname + "_aligningData/mhc_" + hla.lower() + "/result_for_plot", mainwindow.specSeq,
            mainwindow.refSeq, panelValue, panelValue_2)

        subplot_creater.subplotter(
            mainwindow.outputpath + "/" + mainwindow.jobname + "/seq2/results/modified/" + name2 + ".txt",
            mainwindow.outputpath + "/" + mainwindow.jobname + "/seq1/results/modified/" + name1 + ".txt", name1, name2, mainwindow.addfile,
            mainwindow.subplot, mainwindow.outputpath, mainwindow.domain, hla, mainwindow.jobname, panelValue, mainwindow.epitopeRefB, mainwindow.epitopeSpB, mainwindow.popCov, mainwindow.refName, mainwindow.quName, panelValue_2)

    subplot_creater.plt.show()





# main
app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
screen2 = Screen2()
widget.addWidget(mainwindow)
widget.addWidget(screen2)
widget.setGeometry(50,50,1000,600)
widget.setWindowTitle("Comparative Epitope Landscapeprofile")
widget.show()
sys.exit(app.exec_())



