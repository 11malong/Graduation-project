from itertools import islice
import copy
import numpy as np

def del_site(del_index, file_MSIbed, file_out):

    bed_list = []
    with open(file_MSIbed) as f:
        for line in islice(f, 1, None):
            bed_list.append(line.split())
    f.close()
    print(len(bed_list))
    for i in reversed(list(del_index)):
        del bed_list[i]
    print(len(bed_list))

    with open(file_out, 'w') as f:
        f.write(
            'chromosome' + '\t' + 'location' + '\t' + 'repeat_unit_length' + '\t' + 'repeat_unit_binary' + '\t' + 'repeat_times' + '\t' + 'left_flank_binary' + '\t' +
            'right_flank_binary' + '\t' + 'repeat_unit_bases' + '\t' + 'left_flank_bases' + '\t' + 'right_flank_bases' + '\n')

        for i in range(len(bed_list)):
            for j in range(len(bed_list[i])):
                if j == len(bed_list[i]) - 1:
                    f.write(bed_list[i][j] + '\n')
                else:
                    f.write(bed_list[i][j] + '\t')
    f.close()

def get_MSIbed_startPos(file_MSIbed):
    startPos_MSIbed = []
    endPos_MSIbed = []
    with open(file_MSIbed) as f:
        for line in islice(f, 1, None):
            startPos_MSIbed.append(line.split()[1])
            endPos_MSIbed.append(int(line.split()[1]) + int(line.split()[2]) * int(
                line.split()[4]) - 1)
    f.close()

    return startPos_MSIbed, endPos_MSIbed

def get_complexBed_pos(file_complexBed):
    startPos_complexBed = []
    endPos_complexBed = []
    with open(file_complexBed) as f:
        for line in islice(f, 0, None):
            startPos_complexBed.append(line.split()[1])
            endPos_complexBed.append(line.split()[2])
    f.close()

    return startPos_complexBed, endPos_complexBed

def get_delIndex(file_MSIbed, file_complexBed):
    startPos_MSIbed, endPos_MSIbed = get_MSIbed_startPos(file_MSIbed)
    startPos_complexBed, endPos_complexBed = get_complexBed_pos(file_complexBed)
    del_index = set()

    for i in range(len(startPos_MSIbed)):
        for j in range(len(startPos_complexBed)):
            if int(endPos_MSIbed[i]) < int(startPos_complexBed[j]):
                break
            if int(startPos_MSIbed[i]) > int(startPos_complexBed[j]) and int(startPos_MSIbed[i]) < int(
                    endPos_complexBed[j]):
                del_index.add(i)
            elif int(endPos_MSIbed[i]) > int(startPos_complexBed[j]) and int(endPos_MSIbed[i]) < int(
                    endPos_complexBed[j]):
                del_index.add(i)

    return del_index

def get_delFile(file_MSIbed, file_complexBed_list, file_out):
    startPos_MSIbed, endPos_MSIbed = get_MSIbed_startPos(file_MSIbed)

    file_temp = copy.deepcopy(file_MSIbed)
    file_temp = open(file_temp)
    lines = file_temp.readlines()
    print(len(startPos_MSIbed))

    for file_complexBed in file_complexBed_list:
        startPos_complexBed, endPos_complexBed = get_complexBed_pos(file_complexBed)
        print(file_complexBed)

        for i in range(len(startPos_MSIbed)):
            for j in range(len(startPos_complexBed)):
                if int(startPos_MSIbed[i]) > int(startPos_complexBed[j]) and int(startPos_MSIbed[i]) < int(
                        endPos_complexBed[j]):
                    lines.pop(i + 1)
                elif int(endPos_MSIbed[i]) > int(startPos_complexBed[j]) and int(endPos_MSIbed[i]) < int(
                        endPos_complexBed[j]):
                    lines.pop(i + 1)

    file_temp.close()

    file_out = open(file_out, 'w')
    file_out.writelines(lines)
    file_out.close()

    return None

def read_MSIbed(file_MSIbed):
    chr_name = []
    startPos_MSIbed = []
    endPos_MSIbed = []
    with open(file_MSIbed) as f:
        for line in islice(f, 1, None):
            chr_name.append(line.split()[0])
            startPos_MSIbed.append(line.split()[1])
            endPos_MSIbed.append(int(line.split()[1]) + int(line.split()[2]) * int(
                line.split()[4]) - 1)
    f.close()

    return chr_name, startPos_MSIbed, endPos_MSIbed

def read_complexBed(file_complexBed):
    chr_name = []
    startPos_complexBed = []
    endPos_complexBed = []
    with open(file_complexBed) as f:
        for line in islice(f, 0, None):
            chr_name.append(line.split()[0])
            startPos_complexBed.append(line.split()[1])
            endPos_complexBed.append(line.split()[2])
    f.close()

    return chr_name, startPos_complexBed, endPos_complexBed

def clearSite(file_MSIbed, file_complexBed_list, file_out):

    bed_dict = {}

    chrName_MSIbed, startPos_MSIbed, endPos_MSIbed = read_MSIbed(file_MSIbed)

    this_start, last_end = 0, 0
    for i in range(len(chrName_MSIbed)):
        chr_name = chrName_MSIbed[i]
        if chr_name not in bed_dict.keys():
            bed_dict[chr_name] = []
            if i == 0:
                this_start = int(startPos_MSIbed[i])
                continue
            else:
                last_end = int(endPos_MSIbed[i - 1])
            bed_dict[chr_name].append(np.zeros((1, last_end - this_start + 1)))
            bed_dict[chr_name].append(np.array(range(this_start, last_end + 1, 1)))

            this_start = int(startPos_MSIbed[i])

    if len(bed_dict.keys()) == 1:
        this_start, last_end = int(startPos_MSIbed[0]), int(endPos_MSIbed[-1])
        key_name = list(bed_dict.keys())[0]
        bed_dict[key_name] = []
        bed_dict[key_name].append(np.zeros((1, last_end - this_start + 1)))
        bed_dict[key_name].append(np.array(range(this_start, last_end + 1, 1)))

    for complexBed in file_complexBed_list:
        if complexBed == 'None':
            continue
        chrName_complexBed, startPos_complexBed, endPos_complexBed = read_complexBed(complexBed)

        for i in range(len(chrName_complexBed)):
            chr_name = chrName_complexBed[i]
            if chr_name not in bed_dict.keys():
                continue
            else:
                dict_startPos, dict_endPos = bed_dict[chr_name][1][0], bed_dict[chr_name][1][-1]
                startPos_complex, endPos_complex = int(startPos_complexBed[i]), int(
                    endPos_complexBed[i])

                if startPos_complex <= dict_startPos and endPos_complex >= dict_endPos:
                    bed_dict[chr_name][0] = 1
                elif (startPos_complex < dict_startPos and endPos_complex < dict_startPos) or \
                        (startPos_complex > dict_endPos and endPos_complex > dict_endPos):
                    continue
                elif startPos_complex <= dict_startPos and endPos_complex < dict_endPos:
                    bed_dict[chr_name][0][0: endPos_complex - dict_startPos + 1] = 1
                elif startPos_complex > dict_startPos and endPos_complex >= dict_endPos:
                    bed_dict[chr_name][0][startPos_complex - dict_startPos: dict_endPos + 1] = 1
                elif startPos_complex > dict_startPos and endPos_complex < dict_endPos:
                    index_startPos_complexBed = startPos_complex - dict_startPos
                    len_complex = endPos_complex - startPos_complex
                    bed_dict[chr_name][0][index_startPos_complexBed: index_startPos_complexBed + 1 + len_complex] = 1
                else:
                    print('error:get_clearSite.py clearSite()!')

    bed_list = []
    f = open(file_MSIbed)
    for line in f.readlines():
        bed_list.append(line)
    f.close()

    with open(file_out, 'w') as f:
        f.write(
            'chromosome' + '\t' + 'location' + '\t' + 'repeat_unit_length' + '\t' + 'repeat_unit_binary' + '\t' + 'repeat_times' + '\t' + 'left_flank_binary' + '\t' +
            'right_flank_binary' + '\t' + 'repeat_unit_bases' + '\t' + 'left_flank_bases' + '\t' + 'right_flank_bases' + '\n')

        for i in range(len(chrName_MSIbed)):
            print(round(i / len(chrName_MSIbed), 1))
            chr_name = chrName_MSIbed[i]
            num_array, label_array = bed_dict[chr_name][0], bed_dict[chr_name][1]
            index_start, index_end = np.where(label_array == int(startPos_MSIbed[i]))[0][0], \
                                     np.where(label_array == int(endPos_MSIbed[i]))[0][0]
            if np.sum(num_array[index_start:index_end + 1]) > 0:
                continue
            else:
                f.writelines(bed_list[i] + '\n')
    f.close()

def getClearSite_FP(file_bed,file_out):

    out_file=open(f"{file_out}","w")
    out_file.write(
        'chromosome' + '\t' + 'location' + '\t' + 'repeat_unit_length' + '\t' + 'repeat_unit_binary' + '\t' + 'repeat_times' + '\t' + 'left_flank_binary' + '\t' +
        'right_flank_binary' + '\t' + 'repeat_unit_bases' + '\t' + 'left_flank_bases' + '\t' + 'right_flank_bases' + '\n')
    count_filter,count_nofilter=0,0
    with open(file_bed, 'r') as f:
        for line in f:
            if line.startswith('chromosome'):
                continue
            else:
                motif,leftFlank,rightFlank = line.split()[7],line.split()[8],line.split()[9]
                if (leftFlank[-1] in list(motif)) or (rightFlank[0] in list(motif)):
                    count_filter += 1
                    continue
                if len(set(list(leftFlank[-3:])))==1 or len(set(list(rightFlank[:3])))==1:
                    count_filter += 1
                    continue
                out_file.write(line)
                count_nofilter+=1
    print('被过滤的位点数：',count_filter,'没被过滤的位点数：',count_nofilter)

def getClearSite(file_MSIbed, filter_beds, file_out):
    filter_beds_dict = {}
    ms_beds_dict = {}
    chrom_len = {}
    out_file=open(f"{file_out}","w")
    out_file.write(
        'chromosome' + '\t' + 'location' + '\t' + 'repeat_unit_length' + '\t' + 'repeat_unit_binary' + '\t' + 'repeat_times' + '\t' + 'left_flank_binary' + '\t' +
        'right_flank_binary' + '\t' + 'repeat_unit_bases' + '\t' + 'left_flank_bases' + '\t' + 'right_flank_bases' + '\n')
    for item in filter_beds:
        for line in open(item):
            chrom, start, end = line[:-1].split("\t")[:3]
            start, end = int(start), int(end)
            if chrom not in filter_beds_dict:
                filter_beds_dict[chrom] = []
                chrom_len[chrom] = 0
            filter_beds_dict[chrom].append([start, end])
            chrom_len[chrom] = end if chrom_len[chrom] > end else chrom_len[chrom]
    for line in open(file_MSIbed):
        if line[:5]=="chrom":continue
        lineinfo=line[:-1].split("\t")
        chrom,start=lineinfo[:2]
        start,repeat_unit_length,repeat_times=int(start),int(lineinfo[2]),int(lineinfo[4])
        end=start+repeat_unit_length*repeat_times
        if chrom not in ms_beds_dict:
            ms_beds_dict[chrom] = []
        if chrom not in chrom_len:
            chrom_len[chrom] = 0
        ms_beds_dict[chrom].append([start, end,line])
        chrom_len[chrom] = end if chrom_len[chrom] > end else chrom_len[chrom]
    for chrom, chrom_l in chrom_len.items():
        print(f"[Info] Processing {chrom}")
        filter_dots = np.zeros(chrom_l + 1)
        if chrom in filter_beds_dict:
            for region in filter_beds_dict[chrom]:
                filter_dots[region[0]:region[1]] = 1
        if chrom  in ms_beds_dict:
            for region in ms_beds_dict[chrom]:
                if filter_dots[region[0]:region[1]].sum() > 0: continue
                out_file.write(f"{region[2]}")
    out_file.close()

if __name__ == '__main__':

    file_bed = '/mnt/d/science/simulateData/remakeData/realDataToMshunter/2023year/chr1_filter_areaRid.list'
    file_out = '/mnt/d/science/simulateData/remakeData/realDataToMshunter/2023year/chr1_filter_FP.list'
    getClearSite_FP(file_bed, file_out)
