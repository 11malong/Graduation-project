def get_df(file_true,file_predict,mode='homo'):

    if mode == 'homo':
        repeatNum_true_list,pos_true_list,chrName_list = readFile_trueValueLength_homo(file_true)
    elif mode == 'hete':
        file_true_dad, file_true_mum = file_true[0],file_true[1]
        repeatNum_true_list, pos_true_list = readFile_trueLength_hete(file_true_dad, file_true_mum)
    else:
        repeatNum_true_list,pos_true_list = [],[]
        print('ERROR:no such mode!')

    repeatNum_predict_list,pos_predict_list = readFile_predictValue_2(file_predict)

    if len(repeatNum_predict_list)==0 or len(pos_predict_list)==0:
        print('vcf file got no result.')
        return

    df = get_heatDataframe(repeatNum_true_list, repeatNum_predict_list, pos_true_list, pos_predict_list,
                            mode)

    return df

def get_error(file_true,file_predict,image_error_path,variaLen,mode='homo',error_way='normal'):

    df = get_df(file_true,file_predict,mode)
    error_list = []

    P_mean = sum(df.columns.values.tolist())/len(df.columns.values.tolist())
    P_num = len(df.columns.values.tolist())

    row_index_start = int(list(df)[0])

    R_list = []

    for row_index in range(len(df.values)):
        if sum(df.values[row_index]) == 0:
            continue
        R = row_index_start + row_index
        if R == 0:
            continue
        R_list.append(R)

    R_mean = sum(R_list)/len(R_list)

    for row_index in range(len(df.values)):
        if sum(df.values[row_index]) == 0:
            continue
        R = row_index_start + row_index
        if R == 0:
            continue

        error_temp = 0
        error_P,error_R = 0,0

        if error_way == 'normal':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = min([abs(P-R),R]) / R
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / (sum(df.values[row_index]) + 0.001))
        elif error_way == 'RMSE':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = (P-R)**2
                error_temp += error * df.values[row_index][i]
            error_list.append(np.sqrt(error_temp/(sum(df.values[row_index] * P_num))))
        elif error_way == '1/2RMSE':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = 0.5*(P-R)**2
                if abs(P-R) > abs(variaLen)*0.8:continue
                error_temp += error * df.values[row_index][i]
            error_list.append(np.sqrt(error_temp/(sum(df.values[row_index] * P_num))))
        elif error_way == 'log-cosh':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = np.log(np.cosh(P-R))
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / sum(df.values[row_index]))
        elif error_way == 'huber':
            para = 5
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                if abs(P-R) <= para:
                    error = 0.5*(P-R)**2
                else:
                    error = para*abs(P-R)-1*(para)**2
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / sum(df.values[row_index]))
        elif error_way == 'insensitive':
            para = 3
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                if abs(P-R) <= para:
                    error = 0
                else:
                    error = abs(P-R)-para
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / sum(df.values[row_index]))
        elif error_way == 'KL':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = R*(np.log(R)-np.log(P))
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / sum(df.values[row_index]))
        elif error_way == 'sqrt-normal':
            for i in range(len(df.values[row_index])):
                P = row_index_start + i
                error = np.sqrt(min([abs(P-R),R]))
                error_temp += error * df.values[row_index][i]
            error_list.append(error_temp / sum(df.values[row_index]))
        else:
            print('ERROR:no such error_way input!')

    print('均值误差:',sum(error_list)/len(error_list))

    plt.figure()
    plt.plot([R_list[i] for i in range(len(R_list))],error_list,color='g',label='variaLen is %s, motif length = all' % variaLen)
    plt.scatter([R_list[i] for i in range(len(R_list))],error_list,color='r')
    plt.xlabel('Number of repeat units')
    if error_way == '1/2RMSE':
        plt.ylabel('RMSE')
    else:
        plt.ylabel('error rate')
    plt.legend(loc='best',fontsize=8)
    plt.savefig(image_error_path)
    plt.close()

    return error_list,R_list








