def readFasta(file_fa):

    fa_str = ''

    with open(file_fa) as f:
        for line in islice(f, 0, None):
            line = line.strip('\n')
            if ">" in line:
                fa_dict_header = line
            else:
                fa_str = fa_str + line
    f.close()

    return fa_str

def readFasta_mulChr(file_fa):

    fa_str_list = []

    fasta = pysam.FastaFile(file_fa)
    chrName_list = fasta.references
    for chrName in chrName_list:
        fa_str = fasta.fetch(chrName)
        fa_str_list.append(fa_str)

    return fa_str_list,chrName_list

def filter_condition_1(sequence,base,num=2):

    if str(base*num) in sequence:
        return True
    else:
        return False

def remakeData_len_mulChr(file_fa,file_bed,variationLen,str_remake_div=0.1,str_remake_gap=2,method_num = 5):

    variationLen = abs(variationLen)

    if type(variationLen) == set:
        variationLen = int(list(variationLen)[0])

    fa_str_list,chrName_list = readFasta_mulChr(file_fa)

    fa_make_list = copy.deepcopy(fa_str_list)
    change = [0 for i in range(len(chrName_list))]
    remakeAfter = []

    gap = 0

    with open(file_bed) as f:
        for line in f:

            if line.startswith('chrom'):
                continue

            if 0 < gap <= str_remake_gap:
                gap += 1
                continue

            if gap == 0 or gap > str_remake_gap:

                gap = 0

                chrName,repeat_start, repeat_len, repeat_num, repeat_base = line.split()[0],int(line.split()[1]),int(line.split()[2]),float(line.split()[4]),line.split()[7]

                if filter_condition_1(fa_str_list[repeat_start - 10:repeat_start], line.split()[7], num=2):
                    continue

                if method_num == 1:

                    if variationLen > repeat_len*repeat_num:
                        variationLen = repeat_len*repeat_num
                    variationLen=int(variationLen)

                    chrName_index = chrName_list.index(chrName)
                    fa_make_list[chrName_index] = fa_make_list[chrName_index][:repeat_start - change[
                        chrName_index] + 1] + fa_make_list[chrName_index][variationLen + repeat_start - change[chrName_index] + 1:]
                    change[chrName_index] += variationLen
                    repNum_after = repeat_num - (variationLen / repeat_len)

                    remakeAfter.append([chrName, repeat_start, repeat_base, method_num, repeat_num, repNum_after])

                elif method_num == 5:
                    variationLen=int(variationLen)

                    if variationLen % repeat_len == 0:
                        insert_repeatUnit_num = int(variationLen / repeat_len)
                        chrName_index = chrName_list.index(chrName)

                        fa_make_list[chrName_index] = fa_make_list[chrName_index][:repeat_start + change[chrName_index] +1 ] + repeat_base * insert_repeatUnit_num + \
                                                      fa_make_list[chrName_index][repeat_start + change[chrName_index]+1 :]
                        change[chrName_index] += variationLen
                        repNum_after = repeat_num + insert_repeatUnit_num

                    else:
                        insert_repeatUnit_num = int(variationLen // repeat_len)
                        insert_leftBase_len = int(variationLen % repeat_len)
                        insert_leftBase = repeat_base[0:insert_leftBase_len]
                        chrName_index = chrName_list.index(chrName)
                        if repeat_len==4:
                            fa_make_list[chrName_index] = fa_make_list[chrName_index][:repeat_start + change[chrName_index] + 2] + repeat_base * insert_repeatUnit_num + insert_leftBase + \
                                                          fa_make_list[chrName_index][repeat_start + change[chrName_index] + 2:]
                        elif repeat_len == 6:
                            fa_make_list[chrName_index] = fa_make_list[chrName_index][:repeat_start + change[chrName_index] +2] + repeat_base * insert_repeatUnit_num + insert_leftBase + \
                                                          fa_make_list[chrName_index][repeat_start + change[chrName_index] +2:]
                        else:
                            fa_make_list[chrName_index] = fa_make_list[chrName_index][:repeat_start + change[chrName_index]+1 ] + repeat_base * insert_repeatUnit_num + insert_leftBase + \
                                                      fa_make_list[chrName_index][repeat_start + change[chrName_index]+1 :]
                        change[chrName_index] += variationLen
                        repNum_after = repeat_num + (variationLen / repeat_len)

                    remakeAfter.append([chrName,repeat_start,repeat_base,method_num,repeat_num,repNum_after])

                gap += 1
    f.close()

    return fa_make_list,chrName_list,remakeAfter

def get_fileLineCount(file):

    return len(open(file).readlines())

def writeFasta_remake(fa_make,file_out,header = '>chr21'):

    temp = 1
    with open(file_out, 'w') as f:
        f.write(header + '\n')
        for i in range(int(len(fa_make) / 50)):
            f.write(fa_make[0 + i * 50:50 + i * 50])
            f.write('\n')
            temp = i + 1
        f.write(fa_make[temp * 50:])
    f.close()

def writeFile_remake(remake_output,file_out,chrNum = 'chr21'):

    [repStar_list, repEnd_list, repNum_before_list, repBase_before_list, repNum_after_list,
                             repBase_after_list, method_num_list] = remake_output
    with open(file_out,'w') as f:
        f.write('chrNum'+'\t'+'repeat_startPos'+'\t'+'repeat_endPos'+'\t'+'repNum_before[repBase_before]'+'\t'+'repNum_after[repBase_after]'+'\t'+'method_num'+'\n')
        for i in range(len(repStar_list)):
            f.write(chrNum+'\t'+str(repStar_list[i])+'\t'+str(repEnd_list[i])+'\t'+str(repNum_before_list[i])+'['+str(repBase_before_list[i])+']'+'\t'+
                    str(repNum_after_list[i])+'['+str(repBase_after_list[i])+']'+'\t'+str(method_num_list[i]))
            f.write('\n')

    f.close()

def write_out_mulChr(fa_make_list, chrName_list, remakeAfter,file_out_fa,file_out_remake):

    with open(file_out_remake,'w') as f:
        f.write('chrName'+'\t'+'repeat_startPos'+'\t'+'repeat_endPos'+'\t'+'repNum_before[repBase_before]'+'\t'+'repNum_after[repBase_after]'+'\t'+'method_num'+'\n')
        for i in range(len(remakeAfter)):
            f.write(remakeAfter[i][0]+'\t'+str(remakeAfter[i][1])+'\t'+str(remakeAfter[i][1]+remakeAfter[i][4]*len(remakeAfter[i][2]))+'\t'+
                    str(remakeAfter[i][4])+'['+str(remakeAfter[i][2])+']'+'\t'+str(remakeAfter[i][-1])+'['+str(remakeAfter[i][2])+']'+'\t'+str(remakeAfter[i][3]))
            f.write('\n')
    f.close()

    with open(file_out_fa, 'w') as f:
        for i in range(len(chrName_list)):
            f.write('>'+str(chrName_list[i])+'\n')
            print('Processing No ',i,'. chr')
            for j in fa_make_list[i]:
                f.writelines(j)
                f.write('\n')
    f.close()

def remakeData_len_homo_mulChr(file_fa,file_bed,variationLen,file_out_fa,file_out_remake,str_remake_div=0.1,str_remake_gap=2,method_num = 5):

    fa_make_list, chrName_list, remakeAfter = remakeData_len_mulChr(file_fa,file_bed,variationLen,str_remake_div,str_remake_gap,method_num,)

    print('改造了的位点有',len(remakeAfter),'个')

    write_out_mulChr(fa_make_list, chrName_list, remakeAfter,file_out_fa,file_out_remake)

if __name__ == '__main__':

    file_fa = '/mnt/d/science/simulateData/remakeData/chr21.fa'
    file_bed = '/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/hg38Chr21_clear.list'
    file_out_fa = '/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/hahahahaha.fa'
    file_out_remake = '/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/hahahahaha.txt'
    remakeData_len_homo(file_fa,file_bed,40,file_out_fa,file_out_remake)

