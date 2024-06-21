#*****************************位点筛选
from itertools import islice

def site_screen_len(file_bed,file_bed_out,STRLen=None):

    if STRLen == None:
        STRLen = [10,100]

    bed_list = []
    with open(file_bed) as f:
        for line in islice(f, 0,1):
            line_0 = line.split()
            bed_list.append(line_0)
        for line in islice(f, 1, None):
            line = line.split()
            repeat_len,repeat_num = int(line[2]),int(line[4])
            temp_bed = []
            if repeat_len*repeat_num >= STRLen[0] and repeat_len*repeat_num <= STRLen[1]:
                for i in line:
                    temp_bed.append(i)
                bed_list.append(temp_bed)
    f.close()
    print('筛选后的STR位点数：',len(bed_list))

    with open(file_bed_out, 'w') as f:
        for i in range(len(bed_list)):
            for j in range(len(bed_list[i])):
                f.write(str(bed_list[i][j])+'\t')
            f.write('\n')
    f.close()

def site_screen_motifLen(file_bed,file_bed_out,motifLen):

    bed_list = []

    if motifLen == 'all':

        with open(file_bed) as f:
            for line in islice(f, 0, None):
                line = line.split()
                temp_bed = []
                for i in line:
                    temp_bed.append(i)
                bed_list.append(temp_bed)
        f.close()
    elif isinstance(motifLen,list) == True:

        motifLen_min,motifLen_max = motifLen[0],motifLen[1]
        if motifLen_min == '': motifLen_min = 1
        if motifLen_max == '': motifLen_max = 100
        motifLen_list = [i for i in range(motifLen_min,motifLen_max+1,1)]

        with open(file_bed) as f:
            for line in islice(f, 0, 1):
                line_0 = line.split()
                bed_list.append(line_0)
            for line in islice(f, 1, None):
                line = line.split()
                repeat_len = int(line[2])
                temp_bed = []
                if repeat_len in motifLen_list:
                    for i in line:
                        temp_bed.append(i)
                    bed_list.append(temp_bed)
        f.close()

    else:
        with open(file_bed) as f:
            for line in islice(f, 0, 1):
                line_0 = line.split()
                bed_list.append(line_0)
            for line in islice(f, 1, None):
                line = line.split()
                repeat_len = int(line[2])
                temp_bed = []
                if repeat_len == motifLen:
                    for i in line:
                        temp_bed.append(i)
                    bed_list.append(temp_bed)
        f.close()

    print('筛选后的STR位点数：',len(bed_list))

    with open(file_bed_out, 'w') as f:
        for i in range(len(bed_list)):
            for j in range(len(bed_list[i])):
                f.write(str(bed_list[i][j])+'\t')
            f.write('\n')
    f.close()