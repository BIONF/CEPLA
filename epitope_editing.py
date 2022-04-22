import os

def crop_annotation(inpath, outpath, seqpath, tool):

    """
    This Skript filtered out lines of interests from the Epitopedata-set: HLA seqname start stop epitope class percentilVal seqlength
    The default value of the percentile value is 1.0.

    :param inpath: Path of Epitopdata
    :param outpath: Path to output folder
    :param seqpath: Path to sequence of reference or query
    :param tool: mhc class
    :return file with filtered epitopes

    """
    seqs = read_seqs(seqpath)
    for i in seqs:
        name = os.path.splitext(os.path.basename(seqpath))[0]
        if tool == "mhc_i":
            for j in ["8", "9", "10", "11"]:
                full_path = os.path.join(outpath, j, name + ".txt")
                input_path = os.path.join(inpath,j, name+".out")
                if os.path.exists(os.path.join(outpath, j)) and os.path.exists(input_path):
                    out = open(full_path, 'a')
                    if os.path.isfile(input_path):
                        infile = open(input_path, "r")
                        line = infile.readline()
                        line = infile.readline()
                        cells = line.split('\t')
                        while len(cells) > 7:
                            if float(cells[7]) <= 1.0:
                                out.write(cells[0] + '\t' + name + '\t' + cells[2] + '\t' + cells[3] + '\t' + cells[5] + '\t' + tool + '\t' + cells[7] + '\t' + str(i[1]) + '\n')
                                out.flush()
                            line = infile.readline()
                            cells = line.split('\t')
        else:
            for j in ["2"]:
                full_path = os.path.join(outpath, "2", name + ".txt")
                input_path = os.path.join(inpath, j,name+".out")
                if os.path.exists(os.path.join(outpath, "2")) and os.path.exists(input_path):
                    out = open(full_path, 'a')
                    if os.path.isfile(input_path):
                        infile = open(input_path, "r")
                        line = infile.readline()
                        line = infile.readline()
                        cells = line.split('\t')
                        while len(cells) > 7:
                            if float(cells[7]) <= 1.0:
                                out.write(cells[0] + '\t' + name + '\t' + cells[2] + '\t' + cells[3] + '\t' + cells[6] + '\t' + tool + '\t' + cells[7] + '\t' + str(i[1]) + '\n')
                                out.flush()
                            line = infile.readline()
                            cells = line.split('\t')

def read_seqs(seqs):

    """
    Reads the sequence

    :param seqs: Path to the sequence
    :return: returns the sequence in a list

    """
    all_seqs = []

    if os.path.isdir(seqs):
        for file in os.listdir(seqs):
            full_path = os.path.join(seqs,file)
            if os.path.isfile(full_path):
                infile = open(full_path, "r")
                lines = infile.readlines()
                all_seqs.append((lines[0].lstrip(">").rstrip("\n"), len(lines[1])))
        return all_seqs

    else:
        file = open(seqs, "r")
        lines = file.readlines()
        all_seqs.append((lines[0].lstrip(">").rstrip("\n"), len(lines[1])))
        return all_seqs


def combine_epitopes(inpath, outpath):

    """
    This function deletes duplicate entries in the filtered epitopes, combines overlapping epitopes and assigns them the new start and end values.

    :param inpath: Path to the filtered epitopes
    :param outpath: Path to outputfile named 'megafile.txt'
    :return: megafile.txt : all filtered and combined epitopes

    """
    data = {}
    for sub_dir in [os.path.join(inpath, d) for d in os.listdir(inpath)]:
        for file in os.listdir(sub_dir):
            file_path = os.path.join(sub_dir, file)
            if os.path.isfile(file_path):
                infile = open(file_path, "r")
                li = infile.readline()
                li = li.strip("\n")
                cells = li.split("\t")
                k = 0
                while len(cells) > 6:
                    if cells[1] not in data:
                        data[cells[1]] = []
                    if cells not in data[cells[1]]:
                        overlapping = [i for i in range(len(data[cells[1]])) if data[cells[1]][i][0] == cells[0] and data[cells[1]][i][5:] == cells[5:]]
                        if len(overlapping) > 0:
                            combined = []
                            for i in overlapping:
                                start = int(data[cells[1]][i][2])
                                end = int(data[cells[1]][i][3])
                                if not (int(cells[3]) < start or int(cells[2]) > end):
                                    combined.append(i)
                                    b = [[start, end], data[cells[1]][i][4]]
                                    a = [[int(cells[2]), int(cells[3])], cells[4]]
                                    cells[4] = combineTech(a, b)
                                    cells[2] = min(start, int(cells[2]))
                                    cells[3] = max(end, int(cells[3]))
                            combined.sort(reverse=True)
                            for i in combined:
                                data[cells[1]].pop(i)
                            data[cells[1]].append(cells)
                        else:
                            data[cells[1]].append(cells)
                    li = infile.readline()
                    li = li.strip("\n")
                    cells = li.split("\t")
                    k += 1
    outfile = open(outpath, "a")
    for seq in data:
        for line in data[seq]:
            outfile.write(
                line[0] + "\t" + line[1] + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + line[4] + "\t" + line[
                    5] + "\t" + line[6] + "\t" + line[7] + "\n" )
            outfile.flush()


def combineTech(a,b):

    """
    Merges overlapping epitopes into one large epitope.

    :param a: Lukas nochmal fragen
    :param b: Lukas nochmal fragen
    :return: new start end end value

    """
    if a[0][0] < b[0][0]:
        if a[0][1] > b[0][1]:
            return a[1]
        start = a[1][:b[0][0]-a[0][0]]
        return start + b[1]
    else:
        if b[0][1] > a[0][1]:
            return b[1]
        start = b[1][:a[0][0]-b[0][0]]
        return start + a[1]


def get_epitopes(inpath, outpath):

    """
    Filters out the epitope sequence and the associated HLA allele.

    :param inpath: Path to megafile.txt : all filtered and combined epitopes
    :param outpath: Path to output folder. Results will be saved in subfolder 'results'
    :return: filtered file

    """
    infile = open(inpath, "r")
    cnt=0
    megalist = []
    epitoplist = []
    epitopdict = {}
    merkerseq = ""
    k=0

    for line in infile:
        k +=1
        line = line.strip('\n')
        cells = line.split('\t')
        if merkerseq =='' or (cells[1] == megalist[-1][1] and cells[5] == megalist[-1][4]):
            allel = cells[0]
            seq = cells[1]
            start = cells[2]
            end = cells[3]
            peptid = cells[4]
            mclass = cells[5]
            value = cells[6]
            length = cells[7]
            merkerseq = cells[1]
            epitoplist.append(peptid)
            megalist.append([allel,seq,start,end,mclass,value,length,peptid,cnt])

        else:
            full_path1 = os.path.join(outpath, str(megalist[0][1]) + ".txt")
            full_path2 = os.path.join(outpath, megalist[0][4] + "/" + str(megalist[0][1]) + "_" + megalist[0][4] + ".txt")
            out1 = open(full_path1,'a')
            out2 = open(full_path2, 'a')
            l = list(set(epitoplist))
            for i in l:
                counter = epitoplist.count(i)
                for j in megalist:
                    if i in j:
                        index = megalist.index(j)
                        megalist[index][8] = megalist[index][8] + counter
                        out1.write(
                            str(j[0]) + "\t" + str(j[1]) + "\t" + str(j[2]) + "\t" + str(j[3]) + "\t" + str(j[4]) + "\t" + str(
                                j[5]) + "\t" + str(j[6]) + "\t" + str(j[7]) + "\t" + str(j[8]) + "\n")
                        out1.flush()
                        if i not in epitopdict:
                            epitopdict[i] = []
                        epitopdict[i].append(str(j[0]))
            for i in epitopdict:
                allNames = epitopdict[i][0].replace("/", ",HLA-")
                for j in range(1,len(epitopdict[i])):
                    allNames += "," + epitopdict[i][j].replace("/", ",HLA-")
                out2.write(str(i) + "\t" + str(allNames) + "\n")

            epitopdict = {}
            megalist = []
            epitoplist = []
            allel = cells[0]
            seq = cells[1]
            start = cells[2]
            end = cells[3]
            peptid = cells[4]
            mclass = cells[5]
            value = cells[6]
            length = cells[7]
            merkerseq = cells[1]
            epitoplist.append(peptid)
            megalist.append([allel, seq, start, end, mclass, value, length, peptid, cnt])

    full_path1 = os.path.join(outpath, str(megalist[0][1]) + ".txt")
    full_path2 = os.path.join(outpath, megalist[0][4] + "/" + str(megalist[0][1]) + "_" + megalist[0][4] + ".txt")
    out1 = open(full_path1, 'a')
    out2 = open(full_path2, 'a')
    l = list(set(epitoplist))
    for i in l:
        counter = epitoplist.count(i)
        for j in megalist:
            if i in j:
                index = megalist.index(j)
                megalist[index][8] = megalist[index][8] + counter
                out1.write(
                    str(j[0]) + "\t" + str(j[1]) + "\t" + str(j[2]) + "\t" + str(j[3]) + "\t" + str(j[4]) + "\t" + str(
                        j[5]) + "\t" + str(j[6]) + "\t" + str(j[7]) + "\t" + str(j[8]) + "\n")
                out1.flush()
                if i not in epitopdict:
                    epitopdict[i] = []
                epitopdict[i].append(str(j[0]))
    for i in epitopdict:
        allNames = epitopdict[i][0].replace("/", ",HLA-")
        for j in range(1, len(epitopdict[i])):
            allNames += "," + epitopdict[i][j].replace("/", ",HLA-")
        out2.write(str(i) + "\t" + str(allNames) + "\n")

