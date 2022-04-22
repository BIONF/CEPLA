import matplotlib.pyplot as plt
from Bio import SeqIO
import os
from matplotlib.collections import LineCollection
import numpy as np
import matplotlib.patches as mpatches
import csv
import align_process


def subplotter(file1, file2, sname1, sname2, extra, plotNum, outputpath, domainA, mhc, jobname, pcv, epitopeRefB,
               epitopeSpB, popCov, refName, quName, epivalue):

    """
    Generates the subplots in the output.

    :param file1: filtered epitope data of reference
    :param file2: filtered epitope data of query
    :param sname1:  Name of query
    :param sname2: Name of reference
    :param extra: Path to extra file
    :param plotNum: list with all subplotnumbers
    :param outputpath: Path to output folder
    :param domainA: Path to domain architecture
    :param mhc: MHC class
    :param jobname: Name of the Job for the output folder
    :param pcv: percentile value of the second slider
    :param epitopeRefB: Path to the reference epitope data (B-cell)
    :param epitopeSpB: Path to the query epitope data (B-cell)
    :param popCov: Defined population or area
    :param refName: entered reference name for label or titel
    :param quName: entered query name for label or titel
    :param epivalue: epitope score from first slider

    """

    extraName = os.path.splitext(os.path.basename(extra))[0]
    x1_axesList = []
    x2_axesList = []
    x3_axesList = []
    y3_axeslist = []
    y4_axeslist = []
    y5_axeslist = []
    speciesfile1 = open(file1, "r")  # ref
    speciesfile2 = open(file2, "r")

    pcv = float(pcv/100)
    epivalue = float(epivalue/100)
    print("Jetzt umgewandelt episcore: ", epivalue)
    print("Jetzt umgewandelt pcv: ", pcv)

    # for B-cell Dataset
    if not epitopeRefB == '' and not epitopeSpB == '':
        liste = csvReader(epitopeRefB)
        bScoresRef = liste[0]

        liste = csvReader(epitopeSpB)
        bScoresSp = liste[0]

        path = outputpath + "/" + jobname + "/" + jobname + "_aligningData/aligned/TestSpecie/"
        gapPositions = alignReader(path)
        for i in gapPositions[0]:
            for j in i:
                start = float(bScoresRef[j])
                end = float(bScoresRef[j + 1])
                avg = float((end + start) / 2)
                bScoresRef.insert(j, avg)
        for j in gapPositions[1]:
            for k in j:
                start = float(bScoresSp[k])
                end = float(bScoresSp[k + 1])
                avg = float((end + start) / 2)
                bScoresSp.insert(k, avg)


    # for domain architecture
    if not domainA == '':
        with open(domainA, 'r') as domainfile:
            next(domainfile)
            domaindict = []
            namelist = []
            for line in domainfile:
                line = line.rstrip('\n').split(' ')
                start = int(line[0])
                end = int(line[1])
                name = line[2]
                domaindict.append([start, end])
                namelist.append(name)

    for line in speciesfile1:  # ref
        line = line.rstrip('\n').split('\t')
        if line[9] != "":
            if line[4] == "mhc_" + mhc.lower():
                if 1 - float(line[5]) <= float(pcv):
                    x1_axesList.append([int(line[2]), int(line[3]), float(line[9]), int(line[8]), float(line[5])])

        else:
            if line[4] == "mhc_" + mhc.lower():
                if 1 - float(line[5]) <= float(pcv):
                    x1_axesList.append([int(line[2]), int(line[3]), float(0), int(line[8]), float(line[5])])

    epitopPosList = sorted(x1_axesList)
    speciesfile1.close()

    for line2 in speciesfile2:
        line2 = line2.rstrip('\n').split('\t')
        if line2[9] != "":
            if line2[4] == "mhc_" + mhc.lower():
                if 1 - float(line2[5]) <= float(pcv):
                    x2_axesList.append([int(line2[2]), int(line2[3]), float(line2[9]), int(line2[8]), float(line2[5])])
        else:
            if line2[4] == "mhc_" + mhc.lower():
                if 1 - float(line2[5]) <= float(pcv):
                    x2_axesList.append([int(line2[2]), int(line2[3]), float(0), int(line2[8]), float(line2[5])])

    epitopPosList2 = sorted(x2_axesList)
    speciesfile2.close()



    filtered_epitopesRef = epitopeFilter(file1, mhc, epivalue)
    filtered_epitopesSp = epitopeFilter(file2, mhc, epivalue)


    # for extra data
    if not extra == '':
        extrafile = open(extra, "r")
        for line3 in extrafile:
            line3 = line3.rstrip('\n').split('\t')
            x3_axesList.append([int(line3[1]), int(line3[2])])
        epitopPosList3 = sorted(x3_axesList)
        extrafile.close()

    for z in epitopPosList:  # percentil
        z.pop()

    for j in epitopPosList:  # binding
        j.pop()

    for i in epitopPosList:  # population
        i.pop()

    for c in epitopPosList2:  # percentil
        c.pop()

    for a in epitopPosList2:  # binding
        a.pop()

    for b in epitopPosList2:  # population
        b.pop()


#----------------------------------------------------------------------------------------------------------------------#

    # Generation of the subplots

    if plotNum == 'all':
        my_plotting_list = list(range(0, 6))
        if extra != '':
            my_plotting_list.append(6)
        if popCov != 'No':
            my_plotting_list.append(7)
        if epitopeRefB != '' and epitopeSpB != '':
            my_plotting_list.append(8)
        if epitopeRefB != '' and epitopeSpB != '':
            my_plotting_list.append(9)
        if domainA != '':
            my_plotting_list.append(10)
    else:
        my_plotting_list = list(map(int, plotNum.split(",")))

    fig, ax = plt.subplots(len(my_plotting_list), 1, sharex='all')

    fig.suptitle(quName + " and " + refName + " MHC class " + mhc, fontsize=18, fontweight="bold")
    if not isinstance(ax, list):
        ax = [ax]
    if len(my_plotting_list) > 1:
        ax = ax[0]

    startlist = []
    counter = 0
    plt.xlabel('sequence [per AA position]', fontsize=11, fontweight="bold")



    if 0 in my_plotting_list:  # ref
        for i in filtered_epitopesRef:
            ax[counter].plot(i, [0, 0], linewidth=8)
            startlist.append(i[0])
        for i in range(len(startlist)):
            y4_axeslist.append(int("1"))
        ax[counter].plot(startlist, y4_axeslist, color='white')
        ax[counter].set_ylabel('predicted epitopes \n' + refName, fontsize=11, fontweight="bold", rotation=0,
                               labelpad=80)
        ax[counter].axes.yaxis.set_ticks([])

        startlist = []
        counter += 1



    if 1 in my_plotting_list:
        for i in filtered_epitopesSp:
            ax[counter].plot(i, [0, 0], linewidth=8)
            startlist.append(i[0])
        for i in range(len(startlist)):
            y5_axeslist.append(int("1"))
        ax[counter].plot(startlist, y5_axeslist, color='white')
        ax[counter].set_ylabel('predicted epitopes \n' + quName, fontsize=11, fontweight="bold", rotation=0,
                               labelpad=80)
        ax[counter].axes.yaxis.set_ticks([])

        startlist = []
        counter += 1

    result = align_process.plotter(
        outputpath + "/" + jobname + "/" + jobname + "_aligningData/mhc_" + mhc.lower() + "/result_for_plot/TestSpecie/" + str(
            float(pcv)) + "/"+ str(float(epivalue)), sname2, sname1, pcv)
    identitiy_string = seq_identity(
        outputpath + "/" + jobname + "/" + jobname + "_aligningData/aligned/TestSpecie/" + sname1 + "_" + sname2 + "_alignedseqs.txt")



    if 2 in my_plotting_list:
        # The x and y data to plot
        y1 = np.array(result[1])  # ref
        x1 = np.arange(len(y1))
        y2 = np.array(result[7])
        x2 = np.arange(len(y2))

        # Create line segments: 1--2, 2--17, 17--20, 20--16, 16--3, etc.
        segments_x1 = np.r_[x1[0], x1[1:-1].repeat(2), x1[-1]].reshape(-1, 2)
        segments_y1 = np.r_[y1[0], y1[1:-1].repeat(2), y1[-1]].reshape(-1, 2)
        segments_x2 = np.r_[x2[0], x2[1:-1].repeat(2), x2[-1]].reshape(-1, 2)
        segments_y2 = np.r_[y2[0], y2[1:-1].repeat(2), y2[-1]].reshape(-1, 2)

        # Assign colors to the line segments
        lc1 = ["red" if k >= len(segments_y1) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                        segments_y2[k][1] else "grey" for k in range(len(segments_y2))]  # sname1
        lc2 = ["blue" if k >= len(segments_y2) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                         segments_y2[k][1] else "grey" for k in range(len(segments_y1))]  # sname2

        segments2 = [list(zip(x, y)) for x, y in zip(segments_x1, segments_y1)]
        segments1 = [list(zip(x, y)) for x, y in zip(segments_x2, segments_y2)]

        ax[counter].add_collection(LineCollection(segments1, colors=lc1))
        ax[counter].add_collection(LineCollection(segments2, colors=lc2))

        ax[counter].set_ylabel('epitope score', fontsize=11, fontweight="bold", rotation=0, labelpad=60)
        ax[counter].set_ylim(-0.1, max(result[1]) * 1.1)

        counter += 1

    bluelistx = []  # sname2
    bluelisty = []
    redlistx = []  # sname1
    redlisty = []
    blacklisty = [] # for plotting
    if 3 in my_plotting_list:
        for i in range(len(result[1])):
            redlistx.append(result[6][i])
            redlisty.append(result[7][i])
            bluelisty.append(result[1][i])
            bluelistx.append(result[0][i])

        for j in range(len(redlisty)):
            value = bluelisty[j] - redlisty[j]
            blacklisty.append(value)

        ax[counter].plot(bluelistx, blacklisty, color='black')
        ax[counter].set_ylabel('∆ score', fontsize=11, fontweight="bold", rotation=0, labelpad=60)

        counter += 1



    if 4 in my_plotting_list:
        # The x and y data to plot
        y1 = np.array(result[5])  # sname2
        x1 = np.arange(len(y1))
        y2 = np.array(result[11])  # sname1
        x2 = np.arange(len(y2))

        # Create line segments: 1--2, 2--17, 17--20, 20--16, 16--3, etc.
        segments_x1 = np.r_[x1[0], x1[1:-1].repeat(2), x1[-1]].reshape(-1, 2)
        segments_y1 = np.r_[y1[0], y1[1:-1].repeat(2), y1[-1]].reshape(-1, 2)
        segments_x2 = np.r_[x2[0], x2[1:-1].repeat(2), x2[-1]].reshape(-1, 2)
        segments_y2 = np.r_[y2[0], y2[1:-1].repeat(2), y2[-1]].reshape(-1, 2)

        # Assign colors to the line segments
        lc1 = ["red" if k >= len(segments_y1) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                        segments_y2[k][1] else "grey" for k in range(len(segments_y2))]  # sname1
        lc2 = ["blue" if k >= len(segments_y2) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                         segments_y2[k][1] else "grey" for k in range(len(segments_y1))]  # sname2

        segments2 = [list(zip(x, y)) for x, y in zip(segments_x1, segments_y1)]
        segments1 = [list(zip(x, y)) for x, y in zip(segments_x2, segments_y2)]

        ax[counter].add_collection(LineCollection(segments1, colors=lc1))
        ax[counter].add_collection(LineCollection(segments2, colors=lc2))

        ax[counter].set_ylabel('binding number', fontsize=11, fontweight="bold", rotation=0, labelpad=60)
        ax[counter].set_ylim(-0.1, max(result[5]) * 1.1)

        counter += 1



    if 5 in my_plotting_list:
        for i in range(0, len(identitiy_string)):  # identity diagram
            ax[counter].plot([i, i + 1], [-0.16, 0.17], linewidth=1, color='black')
        for i in range(1, len(identitiy_string) - 1):
            ax[counter].plot([i, i + 1], [-0.15, 0.16], linewidth=1, color='white')
            ax[counter].plot([i, i + 1], [1, 1], linewidth=1, color='white')
            ax[counter].plot([i, i + 1], [-1, -1], linewidth=1, color='white')
            ax[counter].set_ylabel('mutational differences', fontsize=11, fontweight="bold", rotation=0,
                                   labelpad=80)
            ax[counter].axes.yaxis.set_ticks([])

        for i in range(0, len(identitiy_string)):
            if identitiy_string[i] == 1:
                pass
            elif identitiy_string[i] == 2:
                ax[counter].plot([i, i], [-0.16, 0.17], linewidth=1, color='black')
            elif identitiy_string[i] == 3:
                ax[counter].plot([i, i], [-0.16, 0.17], linewidth=1, color='green')

        counter += 1



    if 6 in my_plotting_list:  # extrafile
        for i in epitopPosList3:
            ax[counter].plot(i, [0, 0], linewidth=8)
            startlist.append(i[0])
        for i in range(len(startlist)):
            y3_axeslist.append(int("1"))
        ax[counter].plot(startlist, y3_axeslist, color='white')
        ax[counter].set_ylabel('epitopes of \n' + extraName, fontsize=11, fontweight="bold", rotation=0,
                               labelpad=60)

        counter += 1



    if 7 in my_plotting_list:  # population coverage

        # The x and y data to plot
        y1 = np.array(result[3])  # sname2
        x1 = np.arange(len(y1))
        y2 = np.array(result[9])  # sname1
        x2 = np.arange(len(y2))

        # Create line segments: 1--2, 2--17, 17--20, 20--16, 16--3, etc.
        segments_x1 = np.r_[x1[0], x1[1:-1].repeat(2), x1[-1]].reshape(-1, 2)
        segments_y1 = np.r_[y1[0], y1[1:-1].repeat(2), y1[-1]].reshape(-1, 2)
        segments_x2 = np.r_[x2[0], x2[1:-1].repeat(2), x2[-1]].reshape(-1, 2)
        segments_y2 = np.r_[y2[0], y2[1:-1].repeat(2), y2[-1]].reshape(-1, 2)

        # Assign colors to the line segments
        lc1 = ["red" if k >= len(segments_y1) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                        segments_y2[k][1] else "grey" for k in range(len(segments_y2))]  # sname1
        lc2 = ["blue" if k >= len(segments_y2) or segments_y1[k][0] != segments_y2[k][0] or segments_y1[k][1] !=
                         segments_y2[k][1] else "grey" for k in range(len(segments_y1))]  # sname2

        segments2 = [list(zip(x, y)) for x, y in zip(segments_x1, segments_y1)]
        segments1 = [list(zip(x, y)) for x, y in zip(segments_x2, segments_y2)]

        ax[counter].add_collection(LineCollection(segments1, colors=lc1))
        ax[counter].add_collection(LineCollection(segments2, colors=lc2))

        ax[counter].set_ylabel('population coverage', fontsize=11, fontweight="bold", rotation=0, labelpad=60)
        ax[counter].set_ylim(-0.1, max(result[3]) * 1.1)

        counter += 1



    if 8 in my_plotting_list:
        threshold = 0.8
        startlist = []
        for i in range(0, len(bScoresRef)):
            startlist.append(i)
        liste = goodPlot(startlist, bScoresRef, 0.8)
        ax[counter].plot(liste[0], liste[1], linewidth=1, color='blue')
        for j in range(len(startlist)):
            y3_axeslist.append(float("0.8"))
        ax[counter].plot(startlist, y3_axeslist, linewidth=1, color='green')
        ax[counter].fill_between(liste[0], liste[1], [threshold for i in range(len(liste[0]))],
                                 where=[liste[1][i] >= threshold for i in range(len(liste[0]))], color='yellow')
        ax[counter].set_ylabel('predicted b-cell \n epitopes', fontsize=11, fontweight="bold", rotation=0, labelpad=60)

        threshold = 0.8
        startlist2 = []
        for k in range(0, len(bScoresSp)):
            startlist2.append(k)
        liste2 = goodPlot(startlist2, bScoresSp, 0.8)
        ax[counter].plot(liste2[0], liste2[1], linewidth=1, color='red')
        ax[counter].fill_between(liste2[0], liste2[1], [threshold for i in range(len(liste2[0]))],
                                 where=[liste2[1][i] >= threshold for i in range(len(liste2[0]))], color='yellow')

        counter += 1



    bluelistx = []  # sname2
    blacklisty = []  # zum plotten
    if 9 in my_plotting_list:  # heatbeatplot
        for j in range(len(bScoresRef)):
            bluelistx.append(j)

        for k in range(len(bScoresSp)):
            value = bScoresRef[k] - bScoresSp[k]
            blacklisty.append(value)

        ax[counter].plot(bluelistx, blacklisty, color='black')  # epitopescoreBcell
        ax[counter].set_ylabel('∆ score', fontsize=11, fontweight="bold", rotation=0, labelpad=60)

        counter += 1



    if 10 in my_plotting_list:
        cnt = 0
        for i in domaindict:
            n = namelist[cnt]
            pos = i[0]
            ax[counter].plot(i, [0, 0], linewidth=16)
            ax[counter].annotate(n, xy=(pos, 0), xycoords='data')
            ax[counter].set_ylabel('domain architecture', fontsize=11, fontweight="bold", rotation=0,
                                   labelpad=80)
            ax[counter].axes.yaxis.set_ticks([])
            cnt += 1

        counter += 1

    red_patch = mpatches.Patch(color='red', label=quName, linewidth=0.1)
    blue_patch = mpatches.Patch(color='blue', label=refName, linewidth=0.1)
    green_patch = mpatches.Patch(color='green', label='threshold', linewidth=0.1)
    plt.legend(handles=[red_patch, blue_patch, green_patch], bbox_to_anchor=(0.9, 0), loc="lower right",
               bbox_transform=fig.transFigure, ncol=3, fontsize='large')

    figure = plt.gcf()  # get current figure
    figure.set_size_inches(16, 6)  # set figure's size manually to your full screen (16x6)
    fig.savefig(
        outputpath + "/" + jobname + "/" + jobname + "_aligningData/svg-Data/" + sname1 + "_" + sname2 + "_MHC" + mhc + "_" + str(
            pcv) + "_" + str(epivalue) + ".svg", format='svg', bbox_inches='tight', dpi=200)  # here you save it as svg


def alignReader(path):

    files1 = os.listdir(path)
    for i in files1:
        fullpath = os.path.join(path, i)
    record = list(SeqIO.parse(fullpath, 'fasta'))
    seq1 = record[0].seq  # sp
    seq2 = record[1].seq  # ref

    gapPosRef = gap_searcher(seq2)
    gapPosSp = gap_searcher(seq1)

    return gapPosRef, gapPosSp

def gap_searcher(seq):

    lst = np.array(seq)
    pos = np.where(lst =='-')

    return pos


def seq_identity(file):

    """
    Reads the alignment and compares AA for AA where there are mutations.
    If the AA are equal = 1
    If the AA != 2
    If there is a gap = 3

    :param file: Path to Alignment file
    :return: list with numbers (1,2,3) to express the mutation informations

    """
    seq1_string = []
    seq2_string = []
    identity_string = []

    sequences = list(SeqIO.parse(file, 'fasta'))
    spezierecord = sequences[0]
    humanrecord = sequences[1]
    seq1_spezie = spezierecord.seq
    seq2_human = humanrecord.seq

    for i in seq1_spezie:
        seq1_string.append(i)

    for j in seq2_human:
        seq2_string.append(j)

    counter = 0
    for i in seq1_string:
        if i == seq2_string[counter]:
            identity_string.append(int(1))
            counter += 1
        elif i != seq2_string[counter] and i != '-':
            identity_string.append(int(2))
            counter += 1
        elif i == '-':
            identity_string.append(int(3))
            counter += 1

    return identity_string


def csvReader(path):

    """
    Reads the csv from Epidope or BepiPred-2.0 and filters for the scores and epitope position.

    :param path: Path to B-cell Epitope csv file
    :return: list with epitopescore, epitope positions and full sequence

    """
    scoreList = []
    fillListPos = []
    fillListVal = []
    seq = []
    with open(path, newline='') as csvfile:
        for line in csvfile:
            if line.startswith("Position,"):
                spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in spamreader:
                    seq.append(row[1])
                    scoreList.append(float(row[2]))
                    if float(row[2]) >= 0.5:
                        fillListPos.append(row[0])
                        fillListVal.append(float(row[2]))
            if line.startswith("position/header"):
                spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
                next(spamreader)
                for row in spamreader:
                    scoreList.append(float(row[2]))
                    if float(row[2]) >= 0.8:
                        fillListPos.append(row[0])
                        fillListVal.append(float(row[2]))

    fullSeq = "".join(seq)  # führt die einzelnen Buchstaben zusammen zu einer Sequenz

    return scoreList, fillListPos, fillListVal, fullSeq


def goodPlot(x, y, t):

    """
    For existing List of x and y values, create and insert a point (a,b) into the lists,
    where the resulting line plot frm x and y intersects with y = t.

    :param x: x achses value
    :param y: y achses value
    :param t: threshold value


    """
    for i in range(1, len(x)):
        if (y[i - 1] < t and y[i] > t) or (y[i - 1] > t and y[i] < t):
            a = ((t - y[i - 1]) * (x[i] - x[i - 1])) / (y[i] - y[i - 1]) + x[i - 1]
            x = [x[j] for j in range(i)] + [a] + [x[j - 1] for j in range(i + 1, len(x) + 1)]
            y = [y[j] for j in range(i)] + [t] + [y[j - 1] for j in range(i + 1, len(y) + 1)]

    return x, y



def epitopeFilter(epifile, mhc, epivalue):

    """
    It only returns the epitopes with a epitope score >= 0.9 in subplot 1 and 2

    :param epifile: filtered epitope data of reference or query
    :param mhc: MHC class
    :param epivalue: epitope score for first slider. Just for subplot 1,2
    :return: list with epitopes with epitopscore > 0.9 (default)

    """
    filtered_epitopes = []

    epitopefile = open(epifile, "r")
    for line in epitopefile:
        line = line.rstrip('\n').split('\t')
        if line[9] != "":
            if line[4] == "mhc_" + mhc.lower():
                if 1-float(line[5]) >= float(epivalue):  # default 0.9
                    filtered_epitopes.append([int(line[2]), int(line[3])])

    return filtered_epitopes