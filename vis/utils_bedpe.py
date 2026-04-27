import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from util_chrom_start import get_chrom_starts

######### Load Pred Mat ###################
def get_preds(chroms, path):

    y_z_hat_list = []
    tmp = []
    for chrom in chroms:
        tmp.append(np.load(path + '{}.npz'.format(chrom))['arr_0'])
    tmp = np.concatenate([y for y in tmp], axis = 0)
    y_z_hat_list.append(tmp)
    y_z_hat_list = np.concatenate([np.expand_dims(y,1) for y in y_z_hat_list], axis = 1)

    return y_z_hat_list

def get_combined_yhat(y_hat_list, start_ind, end_ind, offset = 200, avg_stripe = False): 
    pred_len = 200
    y_hat_list_reshaped = np.concatenate([x.reshape(1,-1) for x in y_hat_list[start_ind:end_ind]], axis = 0)
    chrom_length = y_hat_list_reshaped.shape[0]
    mat = []

    for i in range(chrom_length):
        mat.append(np.insert(np.zeros(chrom_length+offset+1), i, np.insert(y_hat_list_reshaped[i],pred_len,0)))
    summed = pd.DataFrame(
    np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
        ).fillna(0).T.iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values
    
    if avg_stripe:
        summed = (pd.DataFrame(np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
                ).fillna(0).iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values+pd.DataFrame(
            np.array(mat)).reindex(np.arange(-1*pred_len,chrom_length,1)
                ).fillna(0).T.iloc[0:chrom_length+pred_len,0:chrom_length+pred_len].values)/2
    # if in inference use offset -2MB    
    summed = summed[200:-200,200:-200] # remove padded region
    # if in inference use offset 0
    #summed = summed[:-200,:-200] # remove padded region

    return summed
    
def load_pred_mat(chrom, path):
    y_hat = get_preds([chrom], path)
    pred_mat = get_combined_yhat(y_hat[:,0,:], start_ind = 0, end_ind = y_hat.shape[0], avg_stripe=True)

    return pred_mat

def get_starts(chrom_list, genome, step = 10e3):
    startl = []
    chroml = []
    for chrom in chrom_list:
        chrom = f'chr{chrom}'
        cur_starts = list(np.arange(get_chrom_starts(genome)[0][chrom],get_chrom_starts(genome)[1][chrom]-5000000, step).astype(int))
        startl = startl + cur_starts
        chroml = chroml + list(np.repeat(chrom, len(cur_starts)))
    return startl, chroml
    
def matrix_pairs(chrom, pred_mat, genome):
    st, cl = get_starts([chrom], genome)
    coords = coo_matrix(pred_mat)

    startI = ((coords.row)*10000)
    startJ = ((coords.col)*10000)
    score = coords.data
    str = np.repeat(0, len(coords.data))
    chr = np.repeat(int(chrom), len(coords.data))
    fragI = np.repeat(0, len(coords.data))
    fragJ = np.repeat(1, len(coords.data))
    
    # Create a three-column coordinate value matrix
    coordinate_matrix = np.column_stack((str, chr, startI, fragI, str, chr, startJ, fragJ, score))
    
    pred_df = pd.DataFrame(coordinate_matrix, columns=['strI','chrI','startI', 'fragI', 'strJ', 'chrJ', 'startJ', 'fragJ', 'score'])
    
    int_cols = ['strI','chrI','startI', 'fragI', 'strJ', 'chrJ', 'startJ', 'fragJ']
    for i in int_cols:
        pred_df[i] = pred_df[i].astype(int)

    return pred_df

######### Making Bedpe Files ###################
def makeBedpe(chrom, cell_types, pred_path, bedpe_path, cutoff, genome, loop_size_cutoff, percent_co = True):
    for ct in cell_types:
        pred_mat = load_pred_mat(chrom, f"{pred_path}/{ct}/prediction_{ct}_chr")

        pred_df = matrix_pairs(chrom, pred_mat, genome)
        pred_df = pred_df[np.abs(pred_df['startI']-pred_df['startJ']) > loop_size_cutoff] # Loop length must be larger than 

        if percent_co == False:
            pred_df1 = pred_df[pred_df["score"] > cutoff] # If want to do z-value cutoff instead 
        else:
            pred_df1 = pred_df[pred_df["score"]>np.percentile(pred_df.loc[:,"score"].to_list(), cutoff)]
        #print(np.percentile(pred_df1.loc[:,"score"].to_list(),cutoff))

        
        bedpe_df = pd.DataFrame(np.array([
            [chrom] * pred_df1.shape[0],
            [int(i) for i in pred_df1.loc[:,"startI"].to_list()],
            [int(i+10000) for i in pred_df1.loc[:,"startI"]],
            [chrom] * pred_df1.shape[0],
            [int(i) for i in pred_df1.loc[:,"startJ"].to_list()],
            [int(i+10000) for i in pred_df1.loc[:,"startJ"]],
            ["significant"] * pred_df1.shape[0],
            pred_df1.loc[:,"score"],
            ["."] * pred_df1.shape[0],
            ["."] * pred_df1.shape[0],
        ]).T)

       
        bedpe_df.to_csv(bedpe_path+'/bedpe_files/'+ct+'_chrom'+str(chrom)+'_co'+str(cutoff)+'_10kb.bedpe', sep='\t', header=None, index=None)