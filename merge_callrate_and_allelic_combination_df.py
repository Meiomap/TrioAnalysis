'''Calculate statistics (z-score and modified z-score) for callrates and allelic combinations
'''

import pandas as pd
from collections import defaultdict

import pickle





'''Allelic Combinations
'''
def procedure_process_allelic_combinations():
    d = defaultdict(list)
    df1=pd.read_csv('Families_Dump/Allelic_combinations_thr075.csv')
    df1.fillna(0,inplace=True)

    #calculate zscore:
    for col in df1.columns[4:]:
        col_zscore = col + '_zscore'
        col_madscore= col + '_mzscore'
        df1[col_zscore] = (df1[col] - df1[col].mean()) / df1[col].std(ddof=1)
        df1[col_madscore] = (0.6745*(df1[col] - df1[col].median(skipna=True)))/ df1[col].mad(skipna=True)
       #(0.6745 * (df1[col] - df1[col].median()) / MAD)


    for i,row in df1.iterrows():
        for indi in [row['father'],row['mother'],row['child']]:
            d[indi].append(row['familyid'])

    df1.sort_values(by='child', inplace=True)
    return(df1,d)


#######################################################


'''Call rates
'''

def procedure_process_callrates():
    df2=pd.read_csv('Callrates_thr075_ddof1.csv')

    for col in df2.columns[1:5]:
        col_madscore = col + '_mzscore'
        df2[col_madscore] = (0.6745 * (df2[col] - df2[col].median(skipna=True))) / df2[col].mad(skipna=True)

    df2_series_adj_ar=[]

    for i2,row2 in df2.iterrows():
        familyids=d[row2['individual']]
        if len(familyids)==1:
            df2_series_adj_ar.append(row2)
        else:
            for f in familyids:
                _s=pd.Series(row2.copy(),name=row2.name)
                _s['familyid']=f
                df2_series_adj_ar.append(_s)


    df2_adj=pd.concat(df2_series_adj_ar,axis=1,keys=[s.name for s in df2_series_adj_ar]).transpose()

    #df2.set_index('individual',inplace=True)

    df2_adj.sort_values(by='individual',inplace=True)

    return(df2_adj)

##########################################################
'''Helper function for creating multiindex
'''

def create_multi_col_rec(c):
    if 'mzscore' in c:
        return (''.join(c.split('_')[0:-1]), 'mzscore')
    elif 'zscore' in c:
        return (''.join(c.split('_')[0:-1]),'zscore')
    else:
        return (c,'data')

######################################################

zscore_left=-1.0
zscore_right=2
print('#Zscore thresholds:')
print(zscore_left,zscore_right)


df1,d=procedure_process_allelic_combinations()
df2_adj=procedure_process_callrates()


#merged=pd.merge_asof(left=df2,right=df1,right_on='child',left_on='individual',right_by=['father','mother'],left_by=['individual','individual'])
maindf=pd.merge(left=df2_adj.sort_values(by='familyid'),right=df1.sort_values(by='familyid'),on='familyid',how='outer')
maindf.to_excel('Callrates_and_Allelic_combinations_batch1-5_v3_075.xlsx')
maindf.columns=pd.MultiIndex.from_tuples([create_multi_col_rec(i) for i in maindf.columns])

allfamilies_list=maindf[('familyid','data')].unique()

#for zscore_left,zscore_right in [(-1,1),(-1.5,1.5),(-2,2)]:
#for zscore_left,zscore_right in [(-1.5,1.5)]:



#hqfamilies_df=maindf[((maindf.loc[:,(slice(None),'mzscore')]<zscore_right) & (maindf.loc[:,(slice(None),'mzscore')]>zscore_left)).all(axis=1)]
#hqfamilies_list=set(hqfamilies_df.familyid.iloc[:,0])

hqfamilies_list_adv=[]

for famid,df_fam in maindf.groupby(by=('familyid','data')):
    mother_id=df_fam.loc[:,('mother','data')].iloc[0]
    father_id=df_fam.loc[:,('father','data')].iloc[0]
    child_id=df_fam.loc[:,('child','data')].iloc[0]
    df_parents=df_fam.loc[df_fam.individual.isin([mother_id, father_id]).iloc[:, 0], :]
    columns=[c for c in df_fam.columns if c[1]=='mzscore' and (len(c[0])==4 or 'adj' in c[0] or 'Total' in c[0])]

    zscore_sum=df_fam.loc[(df_fam.individual==child_id).iloc[:,0],(slice(None),'mzscore')].sum(axis=1)

    if not df_parents[((df_parents.loc[:, columns] < zscore_right) & (
                df_parents.loc[:, columns] > zscore_left)).all(axis=1)].empty:
        print(famid,df_parents.loc[:, columns])
        hqfamilies_list_adv.append((famid,abs(zscore_sum.iloc[0])))

ar=[]
for familyid,zscore_sum in hqfamilies_list_adv:
    ar.append((familyid,zscore_sum,familyid.split('.')[0]))

ar2=[]

hq_families_df=pd.DataFrame(ar,columns=['familyid','zscore_sum','familyidBIG'])

for famid, fam_df in hq_families_df.groupby(by='familyidBIG'):
    ar2.append(fam_df[fam_df.zscore_sum==fam_df.zscore_sum.min()].familyid.iloc[0])




print(len(ar2))

print(ar2)


#print(len(hqfamilies_list))
#print(len(hqfamilies_list_adv))


#count=Counter([i.split('.')[0] for i in hqfamilies_list_adv])
#print(len(count.keys()))

#for k in count.items():
#    print(k)




features=['b_allele_freq',
         'combination_trio',
         'log_r_ratio',
         'log_r_ratio_1000',
         'log_r_ratio_500',
         'rank',
         'rank_size',
          'rank_size_nnc',
         'score']

first=True

#nameto_repeat='gn3056.2'

#ar2=['fd14680.1']
#first=False

with open('HQ_families_40families_mzscoreMin075_28_1_2021.csv',newline='',mode='a') as fout:
  for fam in ar2:
      new_fam=pickle.load(open('PickledFamilies/' + fam + '075.pkl','rb'))
      print(fam,new_fam._children)
      new_fam.smooth_logr(1000)
      new_fam.smooth_logr(500)
      new_fam.calculate_supporting_snps()
      dfout=new_fam._dataobject.df[new_fam._children].stack(level='individual').reset_index(level=[1, 2, 3]).reset_index(level=0, drop=True)
      dfout['familyid']=fam
      if first:
          dfout.to_csv(fout,mode='a',index_label=False,index=False,header=True)
          first=False
      else:
          dfout.to_csv(fout, mode='a', index_label=False, index=False, header=False)




