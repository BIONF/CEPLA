import shutil
from Bio.Align.Applications import MuscleCommandline
import numpy as np
import os
from Bio import SeqIO

def seqfilter(seq_aligndict, outputpath,jobname):
    """
    Gets the two sequences in fasta format and aligns them with muscle. After the process the sequences will copy in the outputfolder.

    :param seq_aligndict: Both sequences in one fasta file
    :param outputpath: Path to output folder
    :param jobname: Name of the Job for the output folder

    """
    files1 = os.listdir(seq_aligndict)

    for i in files1:
        filename= i.split('.txt')[0]
        muscle_cline = MuscleCommandline(input= outputpath+ "/"+jobname+"/"+jobname+"_aligningData/Aligndata/sequences_aligned/"+i, out = outputpath+ "/"+jobname+"/"+jobname+"_aligningData/aligned/TestSpecie/" + filename + "_" + "alignedseqs.txt")
        muscle_cline()
    shutil.copy(outputpath+ "/"+jobname+"/"+jobname+"_aligningData/aligned/TestSpecie/" + filename + "_" + "alignedseqs.txt", outputpath+ "/"+jobname+"/"+jobname+"_aligningData/seq")


def startaligner(alignedDict,epitopeDict_spezie, epitopeDict_spezie2,outpath,mhc, jobname, pcv, epivalue):

    """
    Reads all the data and stores it in variables. Then the main align function is called.

    :param alignedDict: Aligned sequences in fasta format
    :param epitopeDict_spezie: Filtered and combined epitope data from the query with eventual population coverage data.
    :param epitopeDict_spezie2: filtered and combined epitope data from the reference with eventual population coverage data.
    :param outpath: Path to output folder
    :param mhc: MHC class
    :param jobname: Name of the job for the output folder
    :param pcv: epitope score value from the second slider for subplot 3, 5, 7

    """

    alignfiles = os.listdir(alignedDict)
    epitopefiles = os.listdir(epitopeDict_spezie)
    epitopefiles2 = os.listdir(epitopeDict_spezie2) #ref

    pathToEpitopefile = os.path.join(epitopeDict_spezie,epitopefiles[0])
    pathToEpitopefile2 = os.path.join(epitopeDict_spezie2,epitopefiles2[0]) #ref

    speziename1 = epitopefiles[0].split(".txt")[0]
    speziename2 = epitopefiles2[0].split(".txt")[0] #ref


    for i in alignfiles:
        path = os.path.join(alignedDict, i)
        file = open(path, "r")
        sequencesResult = seqReader(file)
        aligner(pathToEpitopefile2, pathToEpitopefile, sequencesResult , speziename1 , speziename2, outpath+"/"+jobname+"/"+jobname+"_aligningData/mhc_"+mhc.lower()+"/results_seqs_ALL/", mhc, pcv, epivalue)


def aligner(file_h,file_b,sequencesList, spezie1, spezie2, outputPath, mhc, pcv, epivalue):

    """
    In this method, for each amino acid in the sequence, the epitope values, population coverage values, and
    binding numbers are determind by assigning the maximum.

    :param file_h: Path to reference epitope data
    :param file_b: Path to query epitope data
    :param sequencesList: List with both sequences [query, reference]
    :param spezie1: Name of query
    :param spezie2: Name of reference
    :param outputPath: Path to output folder
    :param mhc: MHC class
    :param pcv: epitope score value of second slider for subplot 3, 5, 7

    """

    speziename = 'TestSpecie'

    file = open(file_h, "r") # ref
    file2 = open(file_b, "r")  # query

    pcv = float(pcv/100)
    epivalue = float(epivalue/100)


    # ref
    seq1_string = sequencesList[1]
    gap_pos1 = gap_searcher(seq1_string)

    # ref
    epitopPosList= fileReader(file,mhc, pcv)
    epitopelist = epitopeScoreCalc(epitopPosList,len(seq1_string))

    # query
    seq2_string = sequencesList[0]
    gap_pos2 = gap_searcher(seq2_string)

    # query
    epitopPosList2 = fileReader(file2,mhc, pcv)
    epitopelist2 = epitopeScoreCalc(epitopPosList2,len(seq2_string))


    # EpitopeScore
    epitopelist_new = insert_gap(gap_pos1,epitopelist,0)

    zIndex = 0

    for z in epitopelist_new:
        if z == 100:
            epitopelist_new[zIndex] = float(1)
        zIndex += 1

    poslist = []
    for j in epitopelist_new:
        value = float(1 - j)
        poslist.append(value)


    epitopelist2_new = insert_gap(gap_pos2, epitopelist2,0)

    zIndex =0
    for z in epitopelist2_new:
        if z == 100:
            epitopelist2_new[zIndex] = float(1)
        zIndex +=1

    poslist2 = []
    for j in epitopelist2_new:
        value = float(1 - j)
        poslist2.append(value)


    file.close()



# here to create the seq-output: (AA     epitopescore    populationCoV   bindingnumber)
    fileOutputCreater(file_h, sequencesList[1], spezie2, speziename, spezie1, outputPath,mhc,pcv, epivalue)
    fileOutputCreater(file_b, sequencesList[0], spezie1, speziename, spezie2, outputPath,mhc,pcv, epivalue)


def fileReader(file,mhc,pcv):

    """
    Reads in the epitope data and filters the epitopes according to the mhc class and the set percentile value.

    :param file: epitope data
    :param mhc: MHC class
    :param pcv: set percentile value
    :return: Returns a list of filtered rows from the epitope data. [start, end, population coverage, mhc class ]

    """
    x1_axesList = []
    for line in file:
        line = line.rstrip('\n').split('\t')
        if line[4] == "mhc_"+mhc.lower():
            if line[9] != "":
                scoreS = 1 - float(line[5]) # line[5] = percentile value
                print(scoreS)
                print(pcv)
                if scoreS >= float(pcv) : # pcv = epitope score
                    x1_axesList.append([int(line[2]),int(line[3]),float(line[9]),int(line[8]), float(line[5])])
            else:

                scoreS = 1 - float(line[5])
                if scoreS >= float(pcv):
                    x1_axesList.append([int(line[2]),int(line[3]),float(0),int(line[8]), float(line[5])])

    epitopPosList= sorted(x1_axesList)

    return epitopPosList

def epitopeScoreCalc(epitopPosList,length):

    """
    Calculates the Epitope Score. (1pcv)

    :param epitopPosList: The filtered epitope data.
    :param length: sequence length
    :return: Returns the list with the calculated Epitope scores.

    """

    #Epitopescore
    filteredList1 = [(i[0], i[1], i[4]) for i in epitopPosList]
    epitopelist = [100 for i in range(length)]
    for i in filteredList1:
        for j in range(i[0], i[1] + 1):
            if epitopelist[j - 1] > i[2]:
                epitopelist[j - 1] = i[2]

    return epitopelist


def popCovCalc(epitopPosList,length):

    """
    Calculates the Population Coverage.

    :param epitopPosList: The filtered epitope data
    :param length: sequence length
    :return: Returns the list with population coverage values

    """
    #populationCov
    filteredList2 = [(i[0], i[1], i[2]) for i in epitopPosList]
    popCovlist = [0 for i in range(length)]
    for i in filteredList2:
        for j in range(i[0], i[1] + 1):
            if popCovlist[j - 1] < i[2]:
                popCovlist[j - 1] = i[2]

    return popCovlist


def bindingNumCalc(epitopPosList,length):

    """
    Calculates the Binding Numbers.

    :param epitopPosList: The filtered epitope data
    :param length: sequence length
    :return: Returns the list with binding numbers

    """
    # bindingNum
    filteredList5 = [(i[0], i[1], i[3]) for i in epitopPosList]
    bindinglist = [0 for i in range(length)]
    for i in filteredList5:
        for j in range(i[0], i[1] + 1):
            if bindinglist[j - 1] < i[2]:
                bindinglist[j - 1] = i[2]
    return bindinglist

def seqReader(file):

    """
    Reads the fasta file and filters out the sequence.

    :param file: Path to sequence
    :return: Returns two lists with the parsed sequences from the fasta file

    """

    seq1_string = []
    seq2_string = []

    sequences = list(SeqIO.parse(file,'fasta'))
    spezierecord= sequences[0]
    humanrecord = sequences[1]
    seq1_spezie = spezierecord.seq
    seq2_human = humanrecord.seq

    for i in seq1_spezie:
        seq1_string.append(i)

    for j in seq2_human:
        seq2_string.append(j)


    return seq1_string,seq2_string


def gap_searcher(seq):

    """
    Finds the gaps in the sequence and remembers the position(s) in a list.

    :param seq: Sequence string
    :return: list with the gap position(s)

    """

    lst = np.array(seq)
    pos = np.where(lst =='-')

    return pos


def insert_gap(gap_pos,list,flag):

    """
    Completes the columns of the gaps with 0.0 or 1.0.

    :param gap_pos: list with the gap position(s)
    :param list: epitope list
    :param flag: Sign for value 0.0 or 1.0
    :return: revised list

    """

    if flag == 1:
        for i in gap_pos:
            for j in i:
                list.insert(j,0.0)
        return list
    else:
        for i in gap_pos:
            for j in i:
                list.insert(j,1.0)
        return list


def fileOutputCreater(file2epitope, seq1_string, name, spezie, alignSeqName, outputPath, mhc, pcv, epivalue):

    """
    Generates the aligned file.

    :param file2epitope: Path to epitope data
    :param seq1_string: sequnce string
    :param name: name of reference or query
    :param spezie: Name of the folder 'Testspecie'
    :param alignSeqName: Path to aligned and filtered fasta file
    :param outputPath: Path to output folder
    :param mhc: MHC class
    :param pcv: epitope score from the second slider
    :param epivalue: epitope score from the first slider

    """

    file1 = open(file2epitope, "r")
    epitopPosList = fileReader(file1,mhc, pcv)
    file1.close()


    gap_pos1 = gap_searcher(seq1_string)
    epitopelist = epitopeScoreCalc(epitopPosList,len(seq1_string))
    popCovlist = popCovCalc(epitopPosList,len(seq1_string))
    bindinglist = bindingNumCalc(epitopPosList,len(seq1_string))

    # EpitopeScore
    epitopelist_new = insert_gap(gap_pos1, epitopelist,0)

    zIndex = 0
    for z in epitopelist_new:
        if z == 100:
            epitopelist_new[zIndex] = float(1)
        zIndex += 1

    poslist = []
    for j in epitopelist_new:
        value = float(1 - j)
        poslist.append(round(value,2))

    # PopulationCoV
    popCovlist_new = insert_gap(gap_pos1, popCovlist,1)

    # BindingNum
    bindinglist_new = insert_gap(gap_pos1, bindinglist,1)

    if not os.path.exists(outputPath + "/" + spezie + "/" + str(pcv) + "/" + str(epivalue)):
        os.makedirs(outputPath + "/"+ spezie + "/" + str(pcv) + "/" + str(epivalue))
    full_path = os.path.join(outputPath, spezie, str(pcv), str(epivalue), name + '_' + alignSeqName + ".txt")
    if os.path.isfile(full_path) == False:
        out = open(full_path, 'w')
        for i in range(len(seq1_string)):
            out.write(str(seq1_string[i]) + '\t' + str(poslist[i]) +'\t' + str(popCovlist_new[i]) + '\t' + str(
                bindinglist_new[i]) + '\n')
        out.close()


def filterGaps(inpath, outpath, name1, name2, pcv, epivalue):

    """
    In this function, the gaps created by the align process are removed. The position of the gaps in the reference is
    noted and deleted at the same point in the query sequence. This is how the reference-based alignment
    process is created.

    :param inpath:
    :param outpath:
    :param name1:
    :param name2:
    :param pcv:
    :return:
    """

    pcv = float(pcv/100)
    epivalue = float(epivalue / 100)

    dir_in = os.listdir(inpath)
    for dir in dir_in:
        path_in = os.path.join(inpath, dir)
        path_out = os.path.join(outpath, dir)
        files = os.listdir(os.path.join(path_in, str(float(pcv)), str(epivalue)))
        data = {}
        if not os.path.exists(path_out + "/" + str(float(pcv)) + "/" + str(float(epivalue))):
            os.makedirs(path_out + "/" + str(float(pcv))+ "/" + str(float(epivalue)))
        for file in files:
            speciename1 = os.path.splitext(os.path.basename(name1))[0]
            speciename2 = os.path.splitext(os.path.basename(name2))[0] # ref

            if file.startswith(name2):
                if speciename2 in data:
                    data[speciename2].append(file)
                else:
                    data[speciename2] = [file]
            else:
                if speciename1 in data:
                    data[speciename1] = [file] + data[speciename1]
                else:
                    data[speciename1] = [file]

        for key in data:
            path_h = os.path.join(path_in, str(float(pcv)), str(float(epivalue)), data[key][1])
            path_a = os.path.join(path_in, str(float(pcv)), str(float(epivalue)), data[key][0])
            file_h = open(path_h, "r").readlines()
            file_a = open(path_a, "r").readlines()

            outpath_h = os.path.join(path_out, str(float(pcv)), str(float(epivalue)), data[key][1])
            outpath_a = os.path.join(path_out, str(float(pcv)), str(float(epivalue)), data[key][0])
            outfile_h = open(outpath_h, "a")
            outfile_a = open(outpath_a, "a")

            for i in range(len(file_h)):
                if not file_h[i].startswith("-"):
                    outfile_h.write(file_h[i])
                    outfile_a.write(file_a[i])


def median(inpath, outpath, name1, name2,pcv, epivalue):

    """
    I had more than one isolate-sequence. In this case, I calculated the median of all sequences values.

    :param inpath: Path to filtered Amino Acid file
    :param outpath: Path to output folder
    :param name1: Name of query
    :param name2: Name of reference
    :param pcv: percentile value of second slider

    """

    pcv = float(pcv/100)
    epivalue = float(epivalue / 100)
    specie = "TestSpecie"

    speciename1 = os.path.splitext(os.path.basename(name1))[0]
    speciename2 = os.path.splitext(os.path.basename(name2))[0]  # ref


    data = {speciename2: [], speciename1: [], "Results": {speciename2:[], speciename1:[]}}
    # load files
    file_inpath = os.path.join(inpath, specie, str(float(pcv)), str(float(epivalue)))
    for file in os.listdir(file_inpath):
        file_path = os.path.join(file_inpath, file)
        content = open(file_path, "r").readlines()
        if file.startswith(speciename2):
            data[speciename2].append(content)
        else:
            data[speciename1].append(content)
    # calculate median
    for key in [speciename2,speciename1]:
        for i in range(0,len(data[key][0])):
            median = []
            values = [[],[],[]]
            for file in data[key]:
                content = file[i].replace("\n","").split("\t")
                values[0].append(float(content[1]))
                values[1].append(float(content[2]))
                values[2].append(float(content[3]))
                if len(median) == 0:
                    median.append(content[0])
            for j in range(len(values)):
                values[j].sort()
                half = len(values[j]) // 2
                if len(values[j]) % 2 == 0:
                    median.append((values[j][half-1]+values[j][half])/2.0)
                else:
                    median.append(values[j][half])
            data["Results"][key].append(median)
    #get files
    for key in data["Results"]:
        if not os.path.exists(outpath + "/" + specie + "/" + str(float(pcv)) + "/" + str(epivalue)):
            os.makedirs(outpath + "/" + specie + "/" + str(float(pcv)) + "/" + str(epivalue))
        out_path = os.path.join(outpath, specie, str(pcv), str(epivalue), key + ".txt")
        if os.path.isfile(out_path) == False:
            outfile = open(out_path, "a")
            for line in data["Results"][key]:
                outfile.write(str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\n")



def startlistMaker(list):

    """
    Generates a list from 1 to length of a template list. [0,1,2,3,4,5,6,7,8....]

    :param list: list with values
    :return: list with ascending numbers

    """
    startlist = []
    for i in range(1, len(list) + 1):
        startlist.append(int(i))
    return startlist


def plotter(dictpath, name1, name2, pcv):

    """
    List all files and returns a long list of all value with corresponding start list

    :param dictpath: Amino acid file 'results_for_plot'
    :param name1: Name of query
    :param name2: Name of reference
    :param pcv: percentile value of second slider
    :return: list with all values and startlist

    """

    speciename1 = os.path.splitext(os.path.basename(name1))[0]  #ref
    speciename2 = os.path.splitext(os.path.basename(name2))[0]
    scorelist_human = []
    popValuelist_human = []
    bindNumlist_human = []

    scorelist_spezie = []
    popValuelist_spezie = []
    bindNumlist_spezie = []

    filedict = os.listdir(dictpath)
    for j in filedict:
        if j.startswith(speciename1): #ref

            fullpath = os.path.join(dictpath, j)
            file1 = open(fullpath, "r")
            for i in file1:
                i = i.rstrip('\n').split('\t')
                if float(i[1]) >= float(pcv):
                    scorelist_human.append(float(i[1]))
                    popValuelist_human.append(float(i[2]))
                    bindNumlist_human.append(float(i[3]))
                else:
                    scorelist_human.append(float(0.0))
                    popValuelist_human.append(float(0.0))
                    bindNumlist_human.append(float(0.0))



        if j.startswith(speciename2):

            fullpath = os.path.join(dictpath, j)
            file1 = open(fullpath, "r")
            for i in file1:
                i = i.rstrip('\n').split('\t')
                if float(i[1]) >= float(pcv):
                    scorelist_spezie.append(float(i[1]))
                    popValuelist_spezie.append(float(i[2]))
                    bindNumlist_spezie.append(float(i[3]))
                else:
                    scorelist_spezie.append(float(0.0))
                    popValuelist_spezie.append(float(0.0))
                    bindNumlist_spezie.append(float(0.0))


    startlist_EScore_h = startlistMaker(scorelist_human)
    startlist_PScore_h = startlistMaker(popValuelist_human)
    startlist_BScore_h = startlistMaker(bindNumlist_human)


    startlist_EScore_b = startlistMaker(scorelist_spezie)
    startlist_PScore_b = startlistMaker(popValuelist_spezie)
    startlist_BScore_b = startlistMaker(bindNumlist_spezie)


    return startlist_EScore_h, scorelist_human, startlist_PScore_h, popValuelist_human, startlist_BScore_h, bindNumlist_human, startlist_EScore_b, scorelist_spezie, startlist_PScore_b, popValuelist_spezie, startlist_BScore_b, bindNumlist_spezie



