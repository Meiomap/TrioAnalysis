import pandas as pd
import SureTypeSC as sc
import os
import pickle
import Family as mainclass
from collections import defaultdict
import sys


STATE=2

batch_names=['batch2',
             'batch3',
             'batch4',
             'batch5',
             'batch6']

samplesheet_path_ar=[r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_batch2\IDAT\204634470174\CYTOSNP_CHIP_1.csv',
                  r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_batch3\IDAT\CYTOSNP_CHIP_2.csv',
                  r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_batch4\IDAT\CYTOSNP_CHIP_MERGED.csv',
                  r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_batch5\IDAT\CYTOSNP_CHIP_8_13.csv',
                  r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_batch6\IDAT\204634470129\Cytosnp_chip_14.csv']


manifest_path_hg38=r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_development\manifest_cluster\HumanCytoSNP-12v2-1_NS550_A1.bpm'
manifest_path_hg19=r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_development\manifest_cluster\hg19\HumanCytoSNP-12v2-1_NS550.bpm'
cluster_path_ar=r'C:\Users\gqc954\Documents\__WORK__\IGENOMIX_development\manifest_cluster\HumanCytoSNP-12v2-1_NS550.egt'


#for  in zip(batch_names,samplesheet_path,[manifest_path],[cluster_path]):


def call_test():
    return(0.1)

#check if there are any duplicates in the samplesheet, if
def check_samplesheet_consistency():
    ar=[]
    new_samplesheet_ar=[]
    for samplesheet,batch_nr in zip(samplesheet_path_ar,batch_names):
        _df=pd.read_csv(samplesheet, skiprows=8)
        _df['batch']=batch_nr
        ar.append(_df)
    _a=pd.concat(ar)
    _a['key']=_a['Sample_ID'].str.lower()
    _a['version']=_a.groupby('key').cumcount()+1
    #_a['Sample_ID']=_a['Sample_ID'].str.cat(_a['version'].astype(str),sep='.')
    _a['Sample_ID'] = _a.apply(lambda x: x['Sample_ID'] if x['version']==1 else x['Sample_ID']+'.' + str(x['version']),axis=1)
    #_a['SentrixBarcode_A']=_a['SentrixBarcode_A'].apply(int)

    for samplesheet_orig,batch_nr in zip(samplesheet_path_ar,batch_names):
        with open(samplesheet_orig) as myfile:
            head = [next(myfile) for x in range(7)]
        print(head)
        new_samplesheet=samplesheet_orig.split('.')[0]+'_VERSIONING.csv'
        new_samplesheet_ar.append(new_samplesheet)
        with open(new_samplesheet,'w', newline='') as sout:
            for line in head: sout.write(line)
            _a[_a['batch']==batch_nr][[c for c in _a.columns if c not in ['key','version','batch']]].to_csv(sout,mode='a',index_label=False,index=False)
    return new_samplesheet_ar



DEFAULT_DIC=os.getcwd()

ASSEMBLY='hg19'
PREFIX=''
if ASSEMBLY=='hg38':
    SUFFIX=''
    manifest_path_ar=manifest_path_hg38
else:
    SUFFIX='hg19'
    manifest_path_ar=manifest_path_hg19

new_samplesheet_ar=check_samplesheet_consistency()

while(STATE<3):
    if STATE==1:
      _loaded_pd_samplesheets=[]
      # load data and serialize

      for batch, samplesheet_path, manifest_path, cluster_path in zip(batch_names,
                                                       new_samplesheet_ar,
                                                       len(samplesheet_path_ar) * [manifest_path_ar],
                                                       len(samplesheet_path_ar) * [cluster_path_ar]):
        #if batch=='batch6' or batch=='batch5':
        #if batch == 'batch5':
          _cur_wd= '\\'.join(samplesheet_path.split('\\')[:-1])
          os.chdir(_cur_wd)
          data_df = sc.basic(manifest_path, cluster_path, samplesheet_path)
          os.chdir(DEFAULT_DIC)
          data_df.to_pickle('Serialized/%s_%s_df%s.pkl' % (str(STATE),str(batch),str(SUFFIX)))

          #_df = pd.read_csv(samplesheet, skiprows=8, sep=',')
          #assert (len(_df.columns) == 5 and _df.columns[0] == 'Sample_ID')
          #_loaded_pd_samplesheets.append(_df)


      #samplesheets_df=pd.concat(_loaded_pd_samplesheets)
      #print(samplesheets_df)

    elif STATE==2:
        _loaded_pd_samplesheets=[]
        _loaded_dataobjects=[]
        for batch, samplesheet_path, manifest_path, cluster_path in zip(batch_names,
                                                                        samplesheet_path_ar,
                                                                        len(samplesheet_path_ar) * [manifest_path_ar],
                                                                        len(samplesheet_path_ar) * [cluster_path_ar]):
            _df = pd.read_csv(samplesheet_path, skiprows=8, sep=',')
            assert (len(_df.columns) == 5 and _df.columns[0] == 'Sample_ID')
            #_loaded_pd_samplesheets.append(_df)

            _d=pd.read_pickle('Serialized/%s_%s_df%s.pkl' % (str(STATE-1), str(batch),str(SUFFIX)))
            _loaded_dataobjects.append(sc.Data.create_from_frame(_d))
              #dfs=sc.Data.create_from_frame(_d)




        dfs=mainclass.merge(_loaded_dataobjects,axis=1)
        #dfs=_loaded_dataobjects[0]
        #print(dfs.df)

        #sample_pd_merged=pd.concat(_loaded_pd_samplesheets)
        #print(sample_pd_merged)

        identity_ar=[]
        ar=[]
        #dfs.apply_NC_threshold_3(0.15,inplace=True)

        dfs.df.columns = dfs.df.columns.set_names(['individual', 'feature'])
        #print(dfs.df.loc[:,(slice(None),'gtype')].stack(level=1).corr(call_test))
        #dfs.df.reset_index(drop=True).stack(level='individual').reset_index()[
        #    ['gtype', 'individual', 'score']].to_csv('Families_Dump/Families_Gtypes_threshold015.stacked.csv')
        #sys.exit()

        #generate identity matrices first

        dfs.apply_NC_threshold_3(0.15, inplace=True)

        for n, dataobject in mainclass.FamilyIterator(dfs).loop_families():
                new_fam = mainclass.Family(n, dataobject)
                new_fam.parse_members()
                new_fam.print_family_members()
                identity_ar.append(new_fam.measure_identity())

        idendities_df=pd.concat(identity_ar)
        #idendities_df.set_index('familyid',append=True).stack(level=0).to_csv('Families_Dump/Identities_long.csv')
        idendities_df.to_csv('Families_Dump/Identities_long_thr015.csv')

        def test(x):
            s=sum(x[['AA','BB','AB']])
            if s==0:
                s+=0.0001
            return x[['NC'!=i for i in x.index]]/s




        d = defaultdict()

        dfs.apply_NC_threshold_3(0.75, inplace=True)

        for n,dataobject in mainclass.FamilyIterator(dfs).loop_families_separate_versions():
            #if n=='17481':#swap father and child
            #print(set(dataobject.df.columns.get_level_values(0)))
            new_fam = mainclass.Family(n, dataobject)
            new_fam.parse_members()
            for i in new_fam.get_array_of_members(): d[i]=n
            #new_fam.print_family_members()
            if new_fam.isIncomplete():
              print('Incomplete families:')
              new_fam.print_family_members()
            #  print(':')

            if not new_fam.isIncomplete():# and new_fam.oneChild():
              #identity_ar.append(new_fam.measure_identity())
              new_fam.print_family_members()
              new_fam._dataobject.df.columns=new_fam._dataobject.df.columns.set_names(['individual', 'feature'])
              #ar.append(new_fam._dataobject.df.reset_index(drop=True).stack(level='individual').reset_index()[['gtype','individual','score']])
              new_fam.calculate_supporting_snps()
              ar.append(new_fam.get_allelic_combinations_df())
              pickle.dump(new_fam, open('PickledFamilies/' + n+ '075.pkl', "wb"))
              #new_fam.dump_result('Families_Dump/%s%s.csv' % (n, SUFFIX))

              #new_fam._dataobject.df.stack(level='individual')


              #new_fam.dump_result('Families_Dump/%s%s.csv' % (n,SUFFIX) )
              #new_fam.export_to_tripod_format('Families_Dump/%s%s.tripod.csv' % (n,SUFFIX))
              print('-----')

        outdf = pd.concat(ar)
        outdf.to_csv('Families_Dump/Allelic_combinations_thr075.csv')


        callrates=dfs.get_call_rates_per_sample().fillna(0)#.transform(lambda x: x / 293552.0).fillna(0)
        callrates_adj=callrates.apply(test,axis=1)
        callrates_adj.columns = pd.Index([i + '_adj' for i in callrates_adj.columns])
        callrates_adj['Total_CR']=dfs.get_call_rates_per_sample().transform(lambda x: x / 293552.0).fillna(0).apply(lambda x: x[['AA','BB','AB']].sum(),axis=1)
        for col in callrates_adj.columns:
            col_zscore = col + '_zscore'
            callrates_adj[col_zscore] = (callrates_adj[col] - callrates_adj[col].mean()) / callrates_adj[col].std(ddof=1)

        #for i,row in callrates_adj.iterrows():
        callrates_adj['familyid']=[d[i] for i in callrates_adj.index.values]
        callrates_adj.to_excel('Callrates_thr075_ddof1.xlsx')
        callrates_adj.to_csv('Callrates_thr075_ddof1.csv')

        #outdf.stack(level=0).to_csv('Families_Dump/Families_Gtypes_threshold015.stacked.csv')

    STATE+=1



























