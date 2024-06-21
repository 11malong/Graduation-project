def readFile_accuracy_hete(file_accuracy_hete):

    variaLen_dad,variaLen_mum,RMSE = [],[],[]
    with open(file_accuracy_hete) as f:
        for line in islice(f, 2, None):
            line = ' '.join(line.split())
            variaLen_dad.append(float(line.split(' ')[0]))
            variaLen_mum.append(float(line.split(' ')[1]))
            RMSE.append(float(line.split(' ')[2]))
    f.close()

    return variaLen_dad,variaLen_mum,RMSE

def draw_accuracy_hete(file_accuracy_hete):

    variaLen_dad, variaLen_mum, accu = readFile_accuracy_hete(file_accuracy_hete)
    variaLen_gap = [variaLen_mum[i] - variaLen_dad[i] for i in range(len(variaLen_dad))]

    plt.figure()
    plt.scatter(variaLen_gap, accu, c='red', s=10, label='legend')
    plt.xlabel("variaLen_gap", fontdict={'size': 10})
    plt.ylabel("accuracy(0/1)", fontdict={'size': 10})
    plt.legend(loc='best')
    for a, b in zip(variaLen_gap, accu):
        plt.text(round(a,3), round(b,3), round(b,3), ha='center', va='center', fontsize=10)

    path = '/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/accuracy_hete.jpg'
    plt.savefig(path)

    plt.show()

def readFile_trueValue_homo(file_true):

    repeatNum_true_list = []
    Pos_list = []

    with open(file_true) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            repeatNum_true_list.append(num)
            Pos_list.append(line[1])
    f.close()

    return repeatNum_true_list,Pos_list

def readFile_trueValueLength_homo(file_true):

    repeatNum_true_list = []
    Pos_list = []

    chrName_list = []

    with open(file_true) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            motifLen = len(line[4].split('[')[1].split(']')[0])
            repeatNum_true_list.append(num*motifLen)
            Pos_list.append(line[1])

            chrName_list.append(line[0])
    f.close()

    return repeatNum_true_list,Pos_list,chrName_list

def readFile_trueValue_hete(file_true_dad,file_true_mum):

    repeatNum_true_list = []
    repeatNum_true_dad,repeatNum_true_mum = [],[]
    Pos_list = []

    with open(file_true_dad) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            repeatNum_true_dad.append(num)
            Pos_list.append(line[1])
    f.close()

    with open(file_true_mum) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            repeatNum_true_mum.append(num)
    f.close()

    for i in range(len(repeatNum_true_dad)):
        temp = []
        temp.append(repeatNum_true_dad[i])
        temp.append(repeatNum_true_mum[i])
        repeatNum_true_list.append(temp)

    return repeatNum_true_list,Pos_list

def readFile_trueLength_hete(file_true_dad,file_true_mum):

    repeatNum_true_list = []
    repeatNum_true_dad,repeatNum_true_mum = [],[]
    Pos_list = []

    with open(file_true_dad) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            motifLen = len(line[4].split('[')[1].split(']')[0])
            repeatNum_true_dad.append(num*motifLen)
            Pos_list.append(line[1])
    f.close()

    with open(file_true_mum) as f:
        for line in islice(f, 1, None):
            line = line.split()
            num = int(float(line[4].split('[')[0]))
            motifLen = len(line[4].split('[')[1].split(']')[0])
            repeatNum_true_mum.append(num*motifLen)
    f.close()

    for i in range(len(repeatNum_true_dad)):
        temp = []
        temp.append(repeatNum_true_dad[i])
        temp.append(repeatNum_true_mum[i])
        repeatNum_true_list.append(temp)

    return repeatNum_true_list,Pos_list

def readFile_predictValue(file_predict):

    AL_list,POS_list = [],[]
    with open(file_predict) as f:
        for line in islice(f, 0, None):
            if line.startswith('#'):
                continue
            temp = line.split('\t')[-1]
            temp = temp.split(':')[-1].strip('\n')
            if temp.split('/')[0].isdecimal() and temp.split('/')[1].isdecimal():
                temp_list = [int(temp.split('/')[0]),int(temp.split('/')[1])]
                AL_list.append(temp_list)
                POS_list.append(line.split('\t')[1])
    f.close()

    return AL_list,POS_list

def readFile_predictValue_2(file_predict):

    AL_list,POS_list = [],[]
    with open(file_predict) as f:
        for line in f:
            if line.startswith('#'):
                continue
            motifLen=len(line.split()[3].split('[')[1].split(']')[0])
            if line.split()[4]=='.':
                continue
            elif ',' in line.split()[4]:
                if motifLen==3 or motifLen==4 or motifLen==6:
                    AL_list.append([int(float(line.split()[4].split(',')[0].split('[')[0])) * motifLen , int(float(line.split()[4].split(',')[1].split('[')[0])) * motifLen])
                else:
                    AL_list.append([round(float(line.split()[4].split(',')[0].split('[')[0]),0) * motifLen , round(float(line.split()[4].split(',')[1].split('[')[0]),0) * motifLen])
            else:
                if motifLen==3 or motifLen==4 or motifLen==6:
                    AL_list.append([int(float(line.split()[4].split('[')[0])) * motifLen , int(float(line.split()[4].split('[')[0])) * motifLen])
                else:
                    AL_list.append([round(float(line.split()[4].split('[')[0]),0) * motifLen , round(float(line.split()[4].split('[')[0]),0) * motifLen])
            POS_list.append(line.split()[1]) #['100002']
    f.close()

    return AL_list,POS_list

def get_Dataframe(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,mode='homo'):

    maxValue,minValue = 0,0

    if isinstance(repeatNum_true_list[0],list):
        true_temp, predict_temp = sum(repeatNum_true_list,[]), sum(repeatNum_predict_list, [])
    else:
        true_temp,predict_temp = repeatNum_true_list,sum(repeatNum_predict_list,[])
    maxValue,minValue = max([max(true_temp),max(predict_temp)]),min([min(true_temp),min(predict_temp)])

    maxValue,minValue=int(maxValue),int(minValue)
    df = pd.DataFrame(np.zeros([maxValue-minValue+1, maxValue-minValue+1]), index=range(minValue,maxValue+1,1),columns=range(minValue,maxValue+1,1))
    df_sum = 0

    for i in range(len(pos_true_list)):
        if pos_true_list[i] in pos_predict_list:

            predict_index = pos_predict_list.index(pos_true_list[i])

            if mode == 'homo':
                R1 = repeatNum_true_list[i]
                P1, P2 = int(repeatNum_predict_list[predict_index][0]), int(repeatNum_predict_list[predict_index][1])

                df.iloc[R1 - minValue,P1 - minValue] = df.iloc[R1 - minValue,P1 - minValue] + 1
                df.iloc[R1 - minValue,P2 - minValue] = df.iloc[R1 - minValue,P2 - minValue] + 1
                df_sum += 2
            elif mode == 'hete':
                R1,R2 = repeatNum_true_list[i][0],repeatNum_true_list[i][1]
                P1,P2 = int(repeatNum_predict_list[predict_index][0]),int(repeatNum_predict_list[predict_index][1])

                if abs(R1-P1) + abs(R2-P2) < abs(R1-P2) + abs(R2-P1):
                    df.iloc[R1 - minValue,P1 - minValue] = df.iloc[R1 - minValue,P1 - minValue] + 1
                    df.iloc[R2 - minValue,P2 - minValue] = df.iloc[R2 - minValue,P2 - minValue] + 1
                else:
                    df.iloc[R1 - minValue, P2 - minValue] = df.iloc[R1 - minValue, P2 - minValue] + 1
                    df.iloc[R2 - minValue, P1 - minValue] = df.iloc[R2 - minValue, P1 - minValue] + 1

                df_sum += 2
            else:
                print('ERROR:no such mode!')
    return df

def get_heatDataframe(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,motifLen,mode='homo'):

    maxValue,minValue = 0,0

    if isinstance(repeatNum_true_list[0],list):
        true_temp, predict_temp = sum(repeatNum_true_list,[]), sum(repeatNum_predict_list, [])
    else:
        true_temp,predict_temp = repeatNum_true_list,sum(repeatNum_predict_list,[])
    maxValue,minValue = max([max(true_temp),max(predict_temp)]),min([min(true_temp),min(predict_temp)])

    maxValue,minValue=int(maxValue),int(minValue)
    df = pd.DataFrame(np.zeros([maxValue-minValue+1, maxValue-minValue+1]), index=range(minValue,maxValue+1,1),columns=range(minValue,maxValue+1,1))
    df_sum = 0

    for i in range(len(pos_true_list)):
        if pos_true_list[i] in pos_predict_list:

            predict_index = pos_predict_list.index(pos_true_list[i])

            if mode == 'homo':
                R1 = repeatNum_true_list[i]
                P1, P2 = int(repeatNum_predict_list[predict_index][0]), int(repeatNum_predict_list[predict_index][1])

                if abs(R1-P1)<=abs(R1-P2):
                    df.iloc[R1 - minValue,P1 - minValue] = df.iloc[R1 - minValue,P1 - minValue] + 1
                    df.iloc[R1 - minValue,P1 - minValue] = df.iloc[R1 - minValue,P1 - minValue] + 1
                else:
                    df.iloc[R1 - minValue,P2 - minValue] = df.iloc[R1 - minValue,P2 - minValue] + 1
                    df.iloc[R1 - minValue,P2 - minValue] = df.iloc[R1 - minValue,P2 - minValue] + 1

                df_sum += 2
            elif mode == 'hete':
                R1,R2 = repeatNum_true_list[i][0],repeatNum_true_list[i][1]
                P1,P2 = round(repeatNum_predict_list[predict_index][0],0),round(repeatNum_predict_list[predict_index][1],0)

                if abs(R1-P1) + abs(R2-P2) < abs(R1-P2) + abs(R2-P1):
                    df.iloc[R1 - minValue,P1 - minValue] = df.iloc[R1 - minValue,P1 - minValue] + 1
                    df.iloc[R2 - minValue,P2 - minValue] = df.iloc[R2 - minValue,P2 - minValue] + 1
                else:
                    df.iloc[R1 - minValue, P2 - minValue] = df.iloc[R1 - minValue, P2 - minValue] + 1
                    df.iloc[R2 - minValue, P1 - minValue] = df.iloc[R2 - minValue, P1 - minValue] + 1

                df_sum += 2
            else:
                print('ERROR:no such mode!')

    index_list = df.columns.values
    list_ToBeDeleted= []
    row_ToBeDeleted,col_ToBeDeleted = [],[]
    flua_para =1

    for i in range(len(index_list)):
        if int(df.iloc[i,:].sum()) == 0:
            row_ToBeDeleted.append(i)
            continue
        for i in range(len(index_list)):
            if int(df.iloc[:,i].sum()) == 0:
                col_ToBeDeleted.append(i)
                continue
            break
        break
    if len(row_ToBeDeleted) <= len(col_ToBeDeleted):
        list_ToBeDeleted.append(row_ToBeDeleted)
    else:
        list_ToBeDeleted.append(col_ToBeDeleted)

    row_ToBeDeleted,col_ToBeDeleted = [],[]
    flua_para = 3

    for i in range(len(index_list)-1,-1,-1):
        if int(df.iloc[i,:].sum()) < flua_para:
            row_ToBeDeleted.append(i)
            continue
        for i in range(len(index_list)-1,-1,-1):
            if int(df.iloc[:,i].sum()) < flua_para:
                col_ToBeDeleted.append(i)
                continue
            break
        break
    if len(row_ToBeDeleted) <= len(col_ToBeDeleted):
        list_ToBeDeleted.append(row_ToBeDeleted)
    else:
        list_ToBeDeleted.append(col_ToBeDeleted)

    list_ToBeDeleted = sum(list_ToBeDeleted,[])
    df = df.drop(df.index[list_ToBeDeleted])
    df.drop(df.columns[list_ToBeDeleted], axis=1,inplace=True)

    return df

def test_haha(pos_predict_list,pos_true_list,repeatNum_true_list,repeatNum_predict_list):

    pos_list = []
    repeatNum_true,repeatNum_predict = [],[]

    for i in range(len(pos_true_list)):
        if pos_true_list[i] in pos_predict_list:

            pos_list.append(pos_true_list[i])
            repeatNum_true.append(repeatNum_true_list[i])
            predict_index = pos_predict_list.index(pos_true_list[i])
            repeatNum_predict.append(repeatNum_predict_list[predict_index])

    dataTest = {'pos':pos_list,
                'repeatNum_true':repeatNum_true,
                'repeatNum_predict':repeatNum_predict}
    df2 = pd.DataFrame(dataTest)

    return df2

def test_output_errSTR(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,variaLen,chrName_list=None,mode='homo'):

    pos_output,repeatNum_true,repeatNum_predict = [],[],[]
    fluctuat_para = 5

    chr_output = []

    if mode == 'homo':
        for i in range(len(pos_true_list)):
            if pos_true_list[i] not in pos_predict_list:
                continue
            predict_index = pos_predict_list.index(pos_true_list[i])
            if (abs(repeatNum_true_list[i] - repeatNum_predict_list[predict_index][0]) > variaLen - fluctuat_para) or \
                    (abs(repeatNum_true_list[i] - repeatNum_predict_list[predict_index][1]) > variaLen - fluctuat_para):
                pos_output.append(pos_true_list[i])
                repeatNum_true.append(repeatNum_true_list[i])
                repeatNum_predict.append(repeatNum_predict_list[predict_index])

                chr_output.append(chrName_list[i])
        dataTest = {'chrName':chr_output,
                    'pos(error STR)': pos_output,
                    'repeatNum_true': repeatNum_true,
                    'repeatNum_predict': repeatNum_predict,
                    }
    elif mode == 'hete':
        countTemp = 0
        for i in range(len(pos_true_list)):
            if pos_true_list[i] not in pos_predict_list:
                countTemp += 1
                continue
            predict_index = pos_predict_list.index(pos_true_list[i])
            true_temp = sorted(repeatNum_true_list[i])
            predict_temp = sorted(repeatNum_predict_list[predict_index])
            if (abs(true_temp[0] - predict_temp[0]) > variaLen - fluctuat_para) or \
                    (abs(true_temp[1] - predict_temp[1]) > variaLen - fluctuat_para):
                pos_output.append(pos_true_list[i])
                repeatNum_true.append(repeatNum_true_list[i])
                repeatNum_predict.append(repeatNum_predict_list[predict_index])

        print('draw_all.py:the number of true STR not in predict STR is: ',countTemp)

        dataTest = {'pos(error STR)':pos_output,
                    'repeatNum_true': repeatNum_true,
                    'repeatNum_predict': repeatNum_predict
                    }
    else:
        dataTest = {}
    df = pd.DataFrame(dataTest)
    df.to_csv('/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/df_errorSTR.csv')

def test_output_accSTR(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,variaLen,mode='homo'):

    pos_output,repeatNum_true,repeatNum_predict = [],[],[]
    fluctuat_para = 3
    countTemp = 0

    if mode == 'hete':

        for i in range(len(pos_true_list)):

            predict_index = pos_predict_list.index(pos_true_list[i])
            true_temp = sorted(repeatNum_true_list[i])
            predict_temp = sorted(repeatNum_predict_list[predict_index])
            if (abs(true_temp[0] - predict_temp[0]) < variaLen - fluctuat_para) and \
                    (abs(true_temp[1] - predict_temp[1]) < variaLen - fluctuat_para):
                countTemp += 1
                pos_output.append(pos_true_list[i])
                repeatNum_true.append(repeatNum_true_list[i])
                repeatNum_predict.append(repeatNum_predict_list[predict_index])

    print('the number of accurate STR site is: ', countTemp)

    dataTest = {'pos(accurate STR)':pos_output,
                'repeatNum_true': repeatNum_true,
                'repeatNum_predict': repeatNum_predict
                }
    df = pd.DataFrame(dataTest)
    df.to_csv('/mnt/d/science/simulateData/remakeData/remakeDataToMshunter/bigScale_test/df_accSTR.csv')

def draw_heatmap(file_true,file_predict,image_path,variaLen,motifLen='all',mode='homo'):

    chrName_list = None

    if mode == 'homo':
        repeatNum_true_list,pos_true_list,chrName_list = readFile_trueValueLength_homo(file_true)
    elif mode == 'hete':
        file_true_dad, file_true_mum = file_true[0],file_true[1]
        repeatNum_true_list, pos_true_list = readFile_trueLength_hete(file_true_dad, file_true_mum)
    else:
        repeatNum_true_list,pos_true_list = [],[]
        print('ERROR:no such mode!')
    repeatNum_predict_list,pos_predict_list = readFile_predictValue_2(file_predict)

    print('predictLen:',len(pos_predict_list))
    print('trueLen:',len(pos_true_list))
    if len(repeatNum_predict_list)==0 or len(pos_predict_list)==0:
        print('vcf file got no result.')
        return
    df = get_heatDataframe(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,motifLen,mode)

    sns.set()
    sns.heatmap(df, center=0, cmap="RdBu_r")
    plt.xlabel('predict value',fontsize=15)
    plt.ylabel('true value',fontsize=15)
    string_text = 'variaLen = '+str(variaLen)+' '+', motif length = all'

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    plt.text(0.1, -0.5, string_text,color='r')
    plt.savefig(image_path)
    plt.close()

def draw_histogram(file_true,file_predict,image_path_real,image_path_predict,variaLen,motifLen,mode='homo'):

    if mode == 'homo':
        repeatNum_true_list,pos_true_list,chrName_list = readFile_trueValueLength_homo(file_true)
        print(len(repeatNum_true_list),len(pos_true_list),len(chrName_list))
    elif mode == 'hete':
        file_true_dad, file_true_mum = file_true[0],file_true[1]
        repeatNum_true_list, pos_true_list = readFile_trueLength_hete(file_true_dad, file_true_mum)
    else:
        repeatNum_true_list,pos_true_list = [],[]
        print('ERROR:no such mode!')
    repeatNum_predict_list,pos_predict_list = readFile_predictValue_2(file_predict)
    print('predictLen:',len(pos_predict_list))
    if len(repeatNum_predict_list)==0 or len(pos_predict_list)==0:
        print('vcf file got no result.')
        return
    df = get_heatDataframe(repeatNum_true_list,repeatNum_predict_list,pos_true_list,pos_predict_list,motifLen,mode)
    index_list = df.columns.values

    df_real = pd.DataFrame(columns=['repeat num','STR site num'])
    repeat_num_real,site_num_real = [],[]
    df_predict = pd.DataFrame(columns=['repeat num','STR site num'])
    repeat_num_predict,site_num_predict = [],[]
    STR_num_all_real,STR_num_all_predict = 0,0


    for i in range(len(index_list)):
        temp = int(df.iloc[i,:].sum()) / 2
        site_num_real.append(temp) #行
        repeat_num_real.append(int(index_list[i]))
    repeat_num_processed,site_num_processed = [],[]
    for i in range(len(repeat_num_real)):
        if int(site_num_real[i]) != 0:
            repeat_num_processed.append(repeat_num_real[i])
            site_num_processed.append(site_num_real[i])
            STR_num_all_real += site_num_real[i]
    df_real['repeat num'] = repeat_num_processed
    df_real['STR site num'] = site_num_processed


    df_real.plot.bar(x="repeat num", y="STR site num",label='variaLen is %s' % variaLen)
    plt.xlabel('True value')
    plt.ylabel('Number of STR')
    x = MultipleLocator(10)

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    plt.legend()
    plt.savefig(image_path_real)
    plt.close()

    for i in range(len(index_list)):
        site_num_predict.append(int(df.iloc[:,i].sum()))
        repeat_num_predict.append(int(index_list[i]))
    repeat_num_processed, site_num_processed = [], []
    for i in range(len(repeat_num_predict)):
        if int(site_num_predict[i]) != 0:
            repeat_num_processed.append(repeat_num_predict[i])
            site_num_processed.append(site_num_predict[i])
            STR_num_all_predict += site_num_predict[i]
    df_predict['repeat num'] = repeat_num_processed
    df_predict['STR site num'] = site_num_processed

    df_predict.plot.bar(x="repeat num", y="STR site num", label='variaLen is %s' % variaLen)
    plt.xlabel('Predict value')
    plt.ylabel('Number of STR')
    x = MultipleLocator(10)

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    plt.legend()
    plt.savefig(image_path_predict)
    plt.close()

    return STR_num_all_real,STR_num_all_predict

def draw_variaLen_accu(data_variaLen,data_accu):

    plt.bar(data_variaLen,data_accu)
    plt.scatter(data_variaLen,data_accu,color='black')
    for a, b in zip(data_variaLen, data_accu):
        plt.text(a, b, b, ha='center', va='bottom')
    plt.xlabel('Variation Length (bp)')
    plt.ylabel('RMSE')
    plt.savefig('/mnt/d/science/simulateData/remakeData/draw_varia_RMSE.jpg')




