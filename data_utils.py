from numpy import isnan
from numpy import nan
from numpy import log2
from statistics import median
from scipy.stats import ttest_ind
from statistics import mean
import pandas as pd

def normalize(data_raw):
    data_log2 = log2(data_raw.replace(0,nan))
    data = data_log2.apply(lambda series: series-median(series.dropna()))
    return data

def get_of_type(cell_type, names, option='contains'):
    if option=='contains':
        cells_of_type = list(s for s in names if cell_type in s)
    elif option=='start':
        cells_of_type = list(s for s in names if s.startswith(cell_type))
    else:
        print ('INVALID OPTION')
        return ('INVALID OPTION')
    return cells_of_type
        
def check_present(row, null_value='nan'):
    bool_row = []
    for i in row:
        if null_value == 'nan':
            bool_row.append(not isnan(i))
        elif null_value== 0:
            bool_row.append(i !=0)
        else:
            bool_row.append(i != null_value)
    return sum(bool_row)

def check_n_of_each_type(row, cell_types=["1_B_", "1_T_"], option='contains', min_reps=3):
    present_in_types = []
    for i in cell_types:
        cells_of_type=get_of_type(i, row.index, option=option)
        data_by_type = row.loc[cells_of_type]
        in_type = check_present(data_by_type)
        three_in_type = bool(in_type >= min_reps)
        present_in_types.append(three_in_type)
    if sum(present_in_types) == len(cell_types):
        return True
    else: return False
    
def presence_by_type(row, cell_types=["1_B_", "1_T_"], option='contains', min_reps=3, null_value='nan'):
    present_in_types = {}
    for i in cell_types:
        cells_of_type=get_of_type(i, row.index, option=option)
        data_by_type = row.loc[cells_of_type]
        in_type = check_present(data_by_type, null_value=null_value)
        three_in_type = bool(in_type >= min_reps)
        present_in_types[i] = three_in_type
    return pd.Series(present_in_types)
    
def check_presence_absence(row, cell_types=["1_B_", "1_T_"], min_reps=3, option='contains'):
    present_in_types = {}
    for i in cell_types:
        cells_of_type=get_of_type(i, row.index, option=option)
        data_by_type = row.loc[cells_of_type]
        in_type = check_present(data_by_type)
        present_in_types[i] = in_type
    if 0 in list(present_in_types.values()):#absent in one type
        if present_in_types[cell_types[0]] >= min_reps:
            return cell_types[0]
        elif present_in_types[cell_types[1]]>= min_reps:
            return cell_types[1]
        
def ttest_wrapper(row, cell_types = ["1_B_", "1_T_"], option='contains'):
    split_row = []
    for i in cell_types:
        cells_of_type=get_of_type(i, row.index, option=option)
        split_row.append(row.loc[cells_of_type])
        
    tstat = ttest_ind(split_row[0],split_row[1])
    tstat = pd.Series(dict(statistic=tstat[0], pvalue=tstat[1]))
    return tstat
        
def get_fold_changes(row, cell_types=["1_B_", "1_T_"], option='contains'):
    means = {}
    for i in cell_types:
        cells_of_type=get_of_type(i, row.index, option=option)
        
        data_by_type = row.loc[cells_of_type]
        means[i] = mean(data_by_type.dropna())
    return (means[cell_types[0]]-means[cell_types[1]])


#get higher-in-B and higher-in-T proteins
def is_altered(tscore, pvalue=.01,change_factor=2, cell_types=["B cells", "T cells"], change_col ='log2(B)-log2(T)'):
    if change_factor > 0:
        log2_fold_change=log2(change_factor)
    else: log2_fold_change=0
    if tscore['pvalue'] < pvalue:
        if tscore[change_col] > log2_fold_change:
            #first type is statistically bigger
            return cell_types[0]
        elif tscore[change_col] < -log2_fold_change:
            #second type is statistically bigger
            return cell_types[1]
        
        
