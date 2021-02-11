import SureTypeSC as sc
import pandas as pd
import itertools
import numpy as np
from collections import Counter
from functools import reduce
from scipy import signal

def group_family(s):
    def member(m):
        m=m.lower().replace('_','')
        if 'm' in m[-1] or 'f' in m[-1]:
            return m[:-1]
        else:
            return m

    if 'gr9746_' in s[0]:
        retval = 'gr9746'
    elif '.' in s[0]:
        _n=s[0].split('.')[0]
        retval=member(_n)
    else:
          retval=member(s[0])
    return retval
    '''
    if 'm' == s[0][-1] or 'f' == s[0][-1]:
      return s[0][:-1]

    elif '.' in s[0]:
        return s[0].split('.')[0]
    else:
      return s[0]
    '''


def merge(dataobjects,axis=1):
    assert(type(dataobjects) == type([]))
    list_of_individuals=set()
    versions=dict()

    for d in dataobjects:
        cur_col_index_vals=d._data.columns.get_level_values(0).values
        cur_individuals=set(d._data.columns.get_level_values(0).values)
        intersect=cur_individuals.intersection(list_of_individuals)
        list_of_individuals=list_of_individuals.union(cur_individuals)
        if len(intersect)>0:
            for k in intersect:
                if k not in versions.keys():
                    versions[k]=2
                else:
                    versions[k]+=1

        new_level_vals=[str(i) + '.'+ str(versions[i]) if i in versions.keys() else i for i in d._data.columns.get_level_values(0)]
        d._data.columns=pd.MultiIndex.from_tuples([t for t in zip(new_level_vals,d._data.columns.get_level_values(1))])

    dataobject_dataframes_ar=[d._data for d in dataobjects]
    #nr_of_samples_with_key=Counter(reduce(lambda x,y: x+y,[list(set(c.columns.get_level_values(0).values)) for c in dataobject_dataframes_ar]))

    if axis==1:
      return sc.Data(df=pd.concat(dataobject_dataframes_ar, axis=1), type=dataobjects[0]._type,
                   transl=dataobjects[0].transl, container=dataobjects[0].container)
    else:
        return sc.Data(df=pd.concat(dataobject_dataframes_ar, axis=0), type=dataobjects[0]._type,
                transl=dataobjects[0].transl, container=dataobjects[0].container)

class FamilyIterator():
    def __init__(self, data):
        """
        Args:
            data:
        """
        self._dataobject = data

    def loop_families_separate_versions(self):
        # return Data object with families assuming pgd_familyid
        for g, d in self._dataobject.df.groupby(by=group_family, axis=1):
            members=set(d.columns.get_level_values(0))
            if any(['.' in i for i in members]):#versioning on, generate a subfamily
                if g!='gr9746':
                  mother_versions=set(sorted([i  for i in members if 'm'  in i.split('.')[0][-1]]))
                  father_versions=set(sorted([i  for i in members if 'f' in i.split('.')[0][-1]]))
                  child_versions=sorted(list(members - (mother_versions | father_versions)))
                else:
                    mother_versions = set(sorted([i for i in members if '2' in i.split('.')[0][-1]]))
                    father_versions = set(sorted([i for i in members if '3' in i.split('.')[0][-1]]))
                    child_versions = sorted(list(members - (mother_versions | father_versions)))


                for i,(mother,father,child) in enumerate(itertools.product(mother_versions, father_versions, child_versions)):
                    newfamid=g + '.' + str(i+1)
                    yield (newfamid, sc.DataLoader.Data(d.loc[:,([mother,father,child])], type=self._dataobject._type, transl=self._dataobject.transl,
                                                 container=self._dataobject.container))
            else:
              yield (g, sc.DataLoader.Data(d, type=self._dataobject._type, transl=self._dataobject.transl,
                                         container=self._dataobject.container))



    def loop_families(self):
        # return Data object with families assuming pgd_familyid
        for g, d in self._dataobject.df.groupby(by=group_family, axis=1):
              yield (g, sc.DataLoader.Data(d, type=self._dataobject._type, transl=self._dataobject.transl,
                                         container=self._dataobject.container))


class Family():
    def __init__(self, familyid, data):
        """
        Args:
            familyid:
            data:
        """
        self._familyid = familyid
        self._fatherid = []
        self._motherid = []
        self._children=[]
        self._dataobject = data
        self.mat=None
        self.pat=None


    def isIncomplete(self):
        ret=False
        if len(self._fatherid)==0 or len(self._motherid)==0 or len(self._children)==0:
            ret=True
        return ret

    def oneChild(self):
        return len(self._children)==1


    def print_family_members(self):
        for key,val in zip(['familyid','father','mother','child'],
                           [self._familyid,','.join(self._fatherid),','.join(self._motherid),','.join(self._children)]):
            #assert(val==''),'Incomplete member %s : %s' % (key,val)
            print('%s:%s'%(key,val))

    def parse_members(self):
        """

        """
        samples = self._dataobject.df.columns.get_level_values(level=0)
        for _s in set(samples):
            if '.' in _s:
                n=_s.split('.')[0]
            else: n=_s

            if "m" == n[-1] or n=='gr9746_1':
                self._motherid.append(_s)
            elif "f" == n[-1] or n=='gr9746_2':
                self._fatherid.append(_s)
            else:
                self._children.append(_s)

        #if self._familyid=='17387':
        #    tmp=self._fatherid
        #    self._fatherid=self._children[0]
        #    self._children[0]=tmp



    def calculate_heterozygosity(self):
        self._dataobject.df.columns.set_names(['individual', 'feature'], inplace=True)

        self._dataobject.df.reset_index(level=2, inplace=True)
        self._dataobject.df.set_index('Position', drop=False, append=True, inplace=True)
        self._dataobject.df.sort_index(level=[1, 2], inplace=True)





    def calculate_supporting_snps(self):
        '''
        condition_ab=(self._dataobject.df.loc[:,(self._motherid,'gtype')]=='AA') & \
                     (self._dataobject.df.loc[:,(self._fatherid,'gtype')]=='BB')
        condition_ba=(self._dataobject.df.loc[:,(self._motherid,'gtype')]=='BB') & \
                     (self._dataobject.df.loc[:,(self._fatherid,'gtype')]=='AA')
        '''

        self._dataobject.df.columns.set_names(['individual','feature'],inplace=True)
        if 'Position' not in self._dataobject.df.columns.get_level_values(0):
          self._dataobject.df.reset_index(level=2,inplace=True)
          self._dataobject.df.set_index('Position',drop=False,append=True,inplace=True)
          self._dataobject.df.sort_index(level=[1, 2], inplace=True)
        #for i,(_mother,_father,_child) in enumerate(itertools.product(self._motherid,self._fatherid,self._children)):
        for _child in self._children:
           '''
           self._dataobject.df[(_child, 'inf')]=np.nan
           self._dataobject.df.loc[:,(_child, 'inf')][condition_ab] = 'ab'
           self._dataobject.df.loc[:,(_child, 'inf')][condition_ba] = 'ba'
           '''
           #self._dataobject.df.loc[:,(_child, 'combination_trio')]=self._dataobject.df.loc[:,[(self._motherid, 'gtype'),
           #                        (self._fatherid, 'gtype'),
           #                        (_child, 'gtype')]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
           self._dataobject.df.loc[:,(_child, 'combination_trio')]=\
               self._dataobject.df.loc[:, (self._motherid + self._fatherid + [_child], 'gtype')].apply(lambda row: ''.join(row.values.astype(str)), axis=1)


           #for g,d in self._dataobject.df.reset_index(level=1).groupby(by=[(_child, 'combination'),('Chr','')]):
           #   print(1)

           #self._dataobject.df.loc[:, (_child, 'combination_parents')]=self._dataobject.df.loc[:,[(self._motherid, 'gtype'),
           #                        (self._fatherid, 'gtype')]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

           self._dataobject.df.loc[:, (_child, 'combination_parents')]=\
               self._dataobject.df.loc[:, (self._motherid + self._fatherid, 'gtype')].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

           #self._dataobject.df.reset_index(level=2,inplace=True)
           #self._dataobject.df.set_index('Position',drop=False,append=True,inplace=True)
           self._dataobject.df.set_index((_child, 'combination_trio'), append=True, inplace=True,drop=False)
           self._dataobject.df.sort_index(level=[1, 2], inplace=True)
           self._dataobject.df[(_child,'rank')]=self._dataobject.df.groupby(level=[1, 3])['Position'].rank()
           #rank het genotypes
           self._dataobject.df.loc[self._dataobject.df.loc[:, (_child, 'gtype')] == 'AB', :].groupby(level=1)['Position'].rank()

           self._dataobject.df.reset_index(level=3,drop=True)

           self._dataobject.df.set_index((_child, 'combination_parents'), append=True, inplace=True,drop=False)
           self._dataobject.df.sort_index(level=[1, 2, 3], inplace=True)
           #self._dataobject.df[(_child, 'rankmax_parents')] = self._dataobject.df.groupby(level=[1, 4])['Position'].rank().cummax()
           mask=['NC' not in i for i in self._dataobject.df.index.get_level_values(3)]
           self._dataobject.df[(_child, 'rank_size')]=self._dataobject.df.groupby(level=[1, 4])['Position'].transform('size')
           self._dataobject.df[(_child, 'rank_size_nnc')] = self._dataobject.df[mask].groupby(level=[1, 4])[
               'Position'].transform('size')
           print('updated_method')


           self._dataobject.df.reset_index(level=[3, 4], drop=True, inplace=True)
           #print(1)


           #print(1)

    def dump_result(self,filename):
        '''
        def _indtype(name):
            if 'f' == name[-1]:
                return 'Father'
            elif 'm' == name[-1]:
                return 'Mother'
            elif '.' in name:
                return 'Child' + name.split('.')[-1]
            else:
                return 'Child'
        '''

        def _indtype(name):
            if name==self._fatherid:
                return 'Father'
            elif name==self._motherid:
                return 'Mother'
            elif '.' in name:
                return 'Child' + name.split('.')[-1]
            elif name in self._children:
                return 'Child'
            else:
                return name



        outputfeatures=['gtype',
                        'combination_trio',
                        'combination_parents',
                        'rank',
                        'rank_size',
                        'individual',
                        'score',
                        'log_r_ratio',
                        'b_allele_freq',
                        'log_r_ratio_1000']
        out_df=self._dataobject.df
        chroms=out_df.index.get_level_values('Chr')
        out_df=out_df.loc[chroms!='0',:]
        out_df.columns = pd.MultiIndex.from_tuples([(_indtype(ind), feat) for ind, feat in out_df.columns])
        out_df.columns.set_names(['individual', 'feature'], inplace=True)
        out_df.stack(level=0).reset_index(level=-1)[outputfeatures].to_csv(filename)





        #self._dataobject.df.loc[:, (self._children,slice(None))].stack(level=0).reset_index(level=-1)[['gtype','combination_trio','combination_parents','rank','rank_size','individual','score','log_r_ratio','b_allele_freq','log_r_ratio_1000']].to_csv(filename)

    def serialize(self,filename):
        pass


    def get_allelic_combinations_df(self):
        seldf=self._dataobject.df.loc[:,(self._children,slice(None))].reset_index(drop=True).stack(level='individual').reset_index()[['combination_trio','combination_parents', 'score']]
        seldf['combination_parents_filt']=seldf.loc[~seldf['combination_parents'].str.contains('NC'), 'combination_parents']
        seldf['combination_trio_filt']=seldf.loc[~seldf['combination_trio'].str.contains('NC'), 'combination_trio']

        #for label, outdf=seldf.loc[:, [i for i in seldf.columns if 'combination' in i]].apply(lambda x: x.value_counts(normalize=True))

        outdf=pd.concat(
            [seldf['combination_trio_filt'].value_counts(normalize=True).to_frame().transpose().reset_index(drop=True),
             seldf['combination_parents_filt'].value_counts(normalize=True).to_frame().transpose().reset_index(drop=True)], axis=1)
        outdf.index=pd.MultiIndex.from_tuples([(self._motherid[0],self._fatherid[0],self._children[0],self._familyid)],names=['mother','father','child','familyid'])
        #outdf=seldf.value_counts(normalize=True).reset_index(level=2)
        #outdf=seldf.loc[:, [i for i in seldf.columns if 'combination' in i]].apply(lambda x: x.value_counts(normalize=True))
        #outdf.columns=outdf.columns.str.replace('combination',self._familyid)
        #outdf.columns=pd.Index(['score',self._familyid])

        return outdf

    def smooth_logr(self,winsize=1000):
        print('running updated smooth_logr')
        for _child in self._children:
            mask=self._dataobject.df.loc[:,[(_child, 'gtype')]]!='NC'
            subject=self._dataobject.df.loc[mask.values,[(_child, 'log_r_ratio')]]
            subject_smoothed=subject.\
                groupby(level=1, as_index=False).\
                rolling(winsize,center=True,min_periods=int(winsize/10)).\
                mean().\
                reset_index(0,drop=True)

            self._dataobject.df[[(_child, 'log_r_ratio_'+str(winsize))]]=subject_smoothed




    def get_array_of_members(self):
        return self._fatherid + self._motherid + self._children


    def measure_identity(self):
        ar=[]
        for (ind1,ind2) in itertools.combinations(self._fatherid + self._motherid + self._children,2):
            df1=self._dataobject.df.loc[:, (ind1, 'gtype')]
            df2=self._dataobject.df.loc[:, (ind2, 'gtype')]
            assert(len(df1)==len(df2))
            condition=(df1!='NC') & (df2!='NC') & (df1==df2)
            totalnnc=sum((df1!='NC') & (df2!='NC'))
            total=len(df1)
            if totalnnc==0:
               val=sum(condition)/(totalnnc+0.0001)
            else:
               val = sum(condition) / (totalnnc + 0.0001)
            ar.append((ind1,ind2,val,totalnnc,totalnnc/total,total))
        retdf=pd.DataFrame.from_records(ar,columns=['ind1','ind2','identity','totalnnc','nncratio','totalsnps'])
        #retdf=retdf.set_index(['ind1', 'ind2']).unstack('ind2')
        retdf['familyid']=self._familyid
        return retdf








    def export_to_tripod_format(self,filename):
        def _indtype(name):
            if 'f' == name[-1]:
                return 'Father'
            elif 'm' == name[-1]:
                return 'Mother'
            else:
                return 'Child'
        rename_dict={'gtype':'GType',
                     'b_allele_freq':'BAF',
                     'log_r_ratio':'LRR'}

        #d = dict(zip(self._dataobject.df.columns.levels[1], ["b1", "c1", "f1"]))
        out_df=self._dataobject.df.loc[:, (slice(None), ['gtype', 'b_allele_freq', 'log_r_ratio'])]
        chroms=out_df.index.get_level_values('Chr')
        out_df = out_df.loc[chroms!='0',:]

        out_df.columns=pd.MultiIndex.from_tuples([(_indtype(ind),rename_dict[feat]) for ind,feat in out_df.columns])
        out_cols = out_df.loc[:, (['Father', 'Mother', 'Child'], ['GType', 'BAF', 'LRR'])]
        out_df.loc[:,(['Father','Mother','Child'],['GType','BAF','LRR'])].to_csv(filename,sep='\t',header=['.'.join(list(c)) for c in out_cols],index=True)
        #SNP_Name
        #Chromosome
        #Position
        #Father.GType
        #Father.BAF
        #Father.LRR
        #Mother.GType
        #Mother.BAF
        #Mother.LRR
        #Child.GType
        #Child.BAF
        #Child.LRR










    # def postpro(self,embryoid,chr='1'):
    #     '''
    #     This function is the calculate the posterior probability of all the embryos of one certain chromosome
    #
    #     '''
    #     if self._motherid!="" and self._fatherid!="":
    #
    #         df = self._dataobject.df.copy();df = df.reset_index(); df = df.sort_values(by = ['Chr','Position'])
    #
    #         datagroup = df.groupby('Chr'); data = datagroup.get_group(chr)
    #
    #         mother_gt=(data[self._motherid, "gtype"]); father_gt = (data[self._fatherid, "gtype"]); parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #         loc = np.where(  (parent_gt == "AABB") | (parent_gt == 'BBAA') |(parent_gt == 'BBBB') | (parent_gt == 'AAAA')); parent_gt = parent_gt[loc] #print parent_gt
    #         #all_p = {}; all_ratios = {};com_seq = {}
    #         #for i in self._embryoids:
    #
    #
    #         i_gt = (data[embryoid,"gtype"]); i_gt = np.asarray(i_gt)[loc]
    #
    #
    #         i_all =  list(parent_gt + i_gt)
    #         #com_seq[i] = i_all
    #         i_prob,i_ratio = prob_cal(i_all,i_gt)
    #         #all_p[i] = i_prob
    #         #all_ratios[i] = i_ratio
    #     return i_prob,i_ratio,i_all#all_p,all_ratios,com_seq
    #
    # def male_or_female(self,embryoid):
    #     df = self._dataobject.df.copy().reset_index(); df = df.sort_values(by = ['Chr','Position'])
    #     datagroup = df.groupby('Chr').get_group('Y')[embryoid,'gtype']
    #     ratio_list = [(Counter(datagroup)['AA'] +  Counter(datagroup)['BB'])/(1.0*datagroup.shape[0]),(Counter(datagroup)['AB'])/(1.0*datagroup.shape[0]),(Counter(datagroup)['NC'])/(1.0*datagroup.shape[0])]
    #     print('Callrate of y chromosome: (based the filter threshold on the Data object): ', ratio_list[0] + ratio_list[1])#dict(zip(['Hom','Het','NC'],ratio_list))
    #     if (ratio_list[0] + ratio_list[1])  > 0.3 :
    #         return 'male'
    #     else: return 'female'
    #
    #
    #
    #
    # def postpro_b_l(self,embryoid,chr='1'):
    #     '''
    #     This function is the calculate the posterior probability of all the embryos of one certain chromosome
    #
    #     '''
    #     #if chr == 'X':
    #     #self._dataobject.df.sort_index(level=1,inplace=True)
    #     #self._dataobject.df.sort_index(level=2,inplace=True)
    #     df = self._dataobject.df.copy();df = df.reset_index(); df = df.sort_values(by = ['Chr','Position'])
    #     datagroup = df.groupby('Chr'); data = datagroup.get_group(chr)
    #     mother_gt=(data[self._motherid, "gtype"]); father_gt = (data[self._fatherid, "gtype"]); parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #     i_gt = np.asarray(data[embryoid,"gtype"]); i_all = parent_gt + i_gt; pos = np.asarray(data['Position'])
    #
    #
    #     loc = np.where((i_all == "AABBAA") | (i_all == 'BBAABB')|(i_all == 'BBAAAB')|(i_all == "AABBAB")|(i_all == "AABBBB")|(i_all == "BBAAAA"));
    #
    #     #parent_gt = parent_gt[loc] #print parent_gt
    #
    #     #i_gt = np.asarray(data[embryoid,"gtype"])[loc]
    #     logr = np.asarray(data[embryoid,"log_r_ratio"])[loc]
    #     baf = np.asarray(data[embryoid,'b_allele_freq'])[loc]
    #     gscore = np.asarray(data[embryoid,'score'])[loc]
    #     print(embryoid,chr)
    #     sc_score = np.asarray(data[embryoid,'rf-gda_ratio:1.0_prob'])[loc]
    #
    #     #print np.mean(baf),np.std(baf),np.mean(gscore),np.std(gscore),np.mean(sc_score),np.std(sc_score)
    #
    #     pos = pos[loc]
    #
    #     i_all =  list(i_all[loc])
    #
    #     to_save = (len(i_gt))*['-1']
    #     loc = loc[0]
    #
    #
    #     #i_prob,i_ratio = prob_cal(i_all,i_gt)
    #     #all_p[i] = i_prob
    #     #all_ratios[i] = i_ratio
    #     return sc_score,gscore,baf,logr,i_all,pos,loc,to_save
    #
    # def postpro_persample(self):
    #     '''
    #     This function is the calculate the posterior probability of all the embryos of all chromosomes
    #     '''
    #     if self._motherid != "" and self._fatherid != "":
    #         df = self._dataobject.df.copy(); df = df.reset_index(); df = df.sort_values(by=['Chr', 'Position'])
    #         mother_gt = (df[self._motherid, "gtype"]); father_gt = (df[self._fatherid, "gtype"])
    #         parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         loc = np.where((parent_gt == 'AAAA') | (parent_gt == "AABB") | (parent_gt == "BBBB") | (parent_gt == 'BBAA'))
    #         parent_gt = parent_gt[loc]
    #         all_p = {}
    #         for i in self._embryoids:
    #             i_gt = (df[i, "gtype"]);i_gt = np.asarray(i_gt)[loc]
    #             i_all = list(parent_gt + i_gt) ; i_prob = prob_cal(i_all,i_gt)
    #             all_p[i] = i_prob
    #
    #         return all_p
    #     else: return 'no parents info'
    #
    #
    #
    #
    # def filter(self,embryoid,seed,chr='2',type='Mat'):
    #     '''
    #     This function is to calculate the emission probability of the method that takes the embryo as seed
    #
    #     Args:
    #         seed: the name of the seed embryo
    #         chr: the name of the chromosome
    #         type: Mat or Pat
    #
    #     '''
    #     if self._motherid != "" and self._fatherid != "":
    #
    #         datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']); datagroup = datagroup.groupby('Chr');data = datagroup.get_group(chr)
    #         mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #
    #
    #
    #
    #         if type =='Mat':
    #             parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         else:
    #             parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #         seed_gt = np.asarray(data[seed,'gtype'])
    #         #all_genotype = {}; all_ratios = {} ; loc_dict = {}
    #         #for i in self._embryoids:
    #         if embryoid != seed:
    #             i_gt = data[embryoid,'gtype'];i_gt = np.asarray(i_gt)
    #             parent_i = np.asarray(parent_gt)
    #             i_all = parent_i + seed_gt + i_gt
    #             final_loc = np.where((i_all == 'ABBBABAB')|(i_all == 'ABBBBBBB')|(i_all == 'ABAAABAB')|(i_all == 'ABAAAAAA')|(i_all == 'ABABBBBB')|(i_all == 'ABABAAAA')
    #                                        |(i_all == 'ABBBBBAB')|(i_all == 'ABBBABBB')|(i_all == 'ABAAAAAB')|(i_all == 'ABAAABAA')|(i_all == 'ABABAABB')|(i_all == 'ABABBBAA'))
    #
    #             i_all = i_all[final_loc]
    #             #loc_dict[i] = final_loc
    #             length = i_all.shape[0]
    #             #all_genotype[i] = i_all
    #
    #             dic = Counter(i_all); ratio =  dict(list(zip(list(dic.keys()),[int(j)/(1.0*length) for j in list(dic.values())])));
    #
    #
    #
    #             #all_ratios[i] = ratio
    #     return i_all,ratio,final_loc
    #
    #
    #
    #
    #
    #
    # def logr_baf(self,chr='2',type='Mat'):
    #     if self._motherid != "" and self._fatherid != "":
    #         datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']); datagroup = datagroup.groupby('Chr')
    #         data = datagroup.get_group(chr); position = np.asarray(data['Position'])
    #         #mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #         if type =='Mat':
    #             mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #             parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         else:
    #             mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #             parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #
    #
    #         for i in self._embryoids:
    #             bafi = data[i,'b_allele_freq']; bafi = np.asarray(bafi)
    #             logri = data[i,'log_r_ratio']; logri = np.asarray(logri)
    #             i_gt = data[i,'gtype'];i_gt = np.asarray(i_gt)
    #             parent_i = np.asarray(parent_gt)
    #             i_all = parent_i  + i_gt
    #             #print i_all.shape
    #             final_loc = np.where((i_all == 'ABAAAB')|(i_all == 'BAAAAB'))#(i_all == 'ABAAAB')|(i_all == 'ABAAAA')|(i_all == 'ABABBB')|(i_all == 'ABABAA')
    #                                 #|(i_all == 'BABBAB')|(i_all == 'BABBBB')|(i_all == 'BAAAAA')|(i_all == 'BAAAAB')|(i_all == 'BAABAA')|(i_all == 'BAABBB'))
    #
    #
    #             np.save(type+i+'_chr'+chr+'_geno_ab.npy',i_all[final_loc])
    #             np.save(type+i+'_chr'+chr+'_baf_ab.npy',bafi[final_loc])
    #             np.save(type+i+'_chr'+chr+'_logr_ab.npy',logri[final_loc])
    #             np.save(type+i+'_chr'+chr+'_pos_ab.npy',position[final_loc])
    #
    #
    #         return 'Finished'
    #
    #
    # def logr_baf_aabb(self,chr='2',type='Mat'):
    #     if self._motherid != "" and self._fatherid != "":
    #         datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']); datagroup = datagroup.groupby('Chr')
    #         data = datagroup.get_group(chr); position = np.asarray(data['Position'])
    #         #mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #         if type =='Mat':
    #             mother_gt = (data[self._motherid, "m_phased"]);father_gt = (data[self._fatherid, "gtype"])
    #             parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         else:
    #             mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "f_phased"])
    #             parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #
    #
    #         for i in self._embryoids:
    #             bafi = data[i,'b_allele_freq']; bafi = np.asarray(bafi)
    #             logri = data[i,'log_r_ratio']; logri = np.asarray(logri)
    #             i_gt = data[i,'gtype'];i_gt = np.asarray(i_gt)
    #             parent_i = np.asarray(parent_gt)
    #             i_all = parent_i  + i_gt
    #             #print i_all.shape
    #             final_loc = np.where((i_all == 'ABAAAB'))#|(i_all == 'ABAABB'))#(i_all == 'ABAAAA'))
    #
    #
    #             np.save(type+i+'_chr'+chr+'_geno_1.npy',i_all[final_loc])
    #             np.save(type+i+'_chr'+chr+'_baf_1.npy',bafi[final_loc])
    #             np.save(type+i+'_chr'+chr+'_logr_1.npy',logri[final_loc])
    #             np.save(type+i+'_chr'+chr+'_pos_1.npy',position[final_loc])
    #
    #
    #         return 'Finished'
    #
    # def mat_male_x(self,seed,embryoid):
    #     datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']);datagroup = datagroup.groupby('Chr');
    #     data = datagroup.get_group('X')
    #     mother_gt = np.asarray(data[self._motherid, "gtype"]); i_gt = np.asarray(data[embryoid,'gtype']); seed_gt = np.asarray(data[seed,'gtype'])
    #     all_i = mother_gt + seed_gt + i_gt ;all_i = list(map(malex, list(all_i)))
    #     return all_i
    #
    #
    # def phased_mat_male_x(self,embryoid):
    #     datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']);datagroup = datagroup.groupby('Chr');
    #     data = datagroup.get_group('X')
    #     mother_gt = np.asarray(data[self._motherid, "m_phased"]); i_gt = np.asarray(data[embryoid,'gtype'])
    #     all_i = mother_gt + i_gt ;all_i = list(map(g_malex, list(all_i)))
    #     return all_i
    #
    #
    #
    #
    #
    #
    #
    #
    # def filter_grand(self,embryoid,chr='2',type='Mat'):
    #     '''
    #     This function is to calculate the emission probablility of the method with grandparent for one certain chromosome
    #     '''
    #
    #
    #     if self._motherid != "" and self._fatherid != "":
    #         datagroup = self._dataobject.df.copy().reset_index() ; datagroup = datagroup.sort_values(by=['Chr','Position']); datagroup = datagroup.groupby('Chr')
    #         data = datagroup.get_group(chr)
    #         #mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "gtype"])
    #         if type =='Mat':
    #             mother_gt = (data[self._motherid, "m_phased"]);father_gt = (data[self._fatherid, "gtype"])
    #             parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         else:
    #             mother_gt = (data[self._motherid, "gtype"]);father_gt = (data[self._fatherid, "f_phased"])
    #             parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #         #all_genotype = {}; all_ratios = {} ; loc_dict = {}
    #
    #         #for i in self._embryoids:
    #
    #         i_gt = data[embryoid,'gtype'];i_gt = np.asarray(i_gt)
    #         parent_i = np.asarray(parent_gt)
    #         i_all = parent_i  + i_gt
    #
    #         final_loc = np.where((i_all == 'ABBBAB')|(i_all == 'ABBBBB')|(i_all == 'ABAAAB')|(i_all == 'ABAAAA')|(i_all == 'ABABBB')|(i_all == 'ABABAA')
    #                                 |(i_all == 'BABBAB')|(i_all == 'BABBBB')|(i_all == 'BAAAAA')|(i_all == 'BAAAAB')|(i_all == 'BAABAA')|(i_all == 'BAABBB'))
    #
    #         i_all = i_all[final_loc]; #print i_all[-1]
    #         #loc_dict[i] = final_loc
    #         length = i_all.shape[0]
    #         #all_genotype[i] = i_all
    #         dic = Counter(i_all); ratio =  dict(list(zip(list(dic.keys()),[int(j)/(1.0*length) for j in list(dic.values())])));
    #         #all_ratios[i] = ratio; #bafi = np.asarray(data[i,'b_allele_freq']), logri = np.asarray(data[i,'log_r_ratio'])
    #             #print bafi.shape, logri.shape
    #             #baf[i] = bafi[final_loc]
    #             #logr[i] = logri[final_loc]
    #
    #     return i_all,ratio,final_loc#all_genotype,all_ratios,loc_dict
    #
    # def filter_persample(self,seed,chr='2',type='Mat'):
    #     '''
    #     This function is to calculate the emission probability based on all chromosomes
    #     '''
    #
    #     if self._motherid != "" and self._fatherid != "":
    #         df = self._dataobject.df.copy().reset_index(); df.sort_values(by=['Chr', 'Position']); datagroup = df.groupby('Chr')
    #         data = datagroup.get_group(chr)
    #         mother_gt = (df[self._motherid, "gtype"]); father_gt = (df[self._fatherid, "gtype"])
    #         if type =='Mat':
    #             parent_gt = np.asarray(mother_gt) + np.asarray(father_gt)
    #         else:
    #             parent_gt = np.asarray(father_gt) + np.asarray(mother_gt)
    #         seed_gt = np.asarray(df[seed, 'gtype'])
    #         all_genotype = {}; all_ratios = {}; loc_dict = {}
    #
    #         for i in self._embryoids:
    #             if i != seed:
    #                 i_gt = df[i,'gtype'];i_gt = np.asarray(i_gt)
    #                 parent_i = np.asarray(parent_gt)
    #                 i_all = parent_i + seed_gt + i_gt
    #                 final_loc = np.where((i_all == 'ABBBABAB')|(i_all == 'ABBBBBBB')|(i_all == 'ABAAABAB')|(i_all == 'ABAAAAAA')|(i_all == 'ABABBBBB')|(i_all == 'ABABAAAA')
    #                                        |(i_all == 'ABBBBBAB')|(i_all == 'ABBBABBB')|(i_all == 'ABAAAAAB')|(i_all == 'ABAAABAA')|(i_all == 'ABABAABB')|(i_all == 'ABABBBAA'))
    #
    #                 i_all = i_all[final_loc] #print i_all[-1]
    #                 #loc_dict[i] = final_loc
    #                 length = i_all.shape[0]
    #                 #all_genotype[i] = i_all
    #
    #                 dic = Counter(i_all); ratio =  dict(list(zip(list(dic.keys()),[int(j)/(1.0*length) for j in list(dic.values())])));
    #                 all_ratios[i] = ratio
    #
    #                 if type == 'Mat':
    #                     i_all_chr = np.asarray(data[self._motherid, "gtype"]) + np.asarray(data[self._fatherid, "gtype"]) + np.asarray(data[seed,'gtype']) + np.asarray(data[i,'gtype'])
    #                 else:
    #                     i_all_chr = np.asarray(data[self._fatherid, "gtype"]) + np.asarray(data[self._motherid, "gtype"]) + np.asarray(data[seed, 'gtype']) + np.asarray(
    #                         data[i, 'gtype'])
    #                 final_loc = np.where((i_all_chr == 'ABBBABAB') | (i_all_chr == 'ABBBBBBB') | (i_all_chr == 'ABAAABAB') | (
    #                                 i_all_chr == 'ABAAAAAA') | (i_all_chr == 'ABABBBBB') | (i_all_chr == 'ABABAAAA')
    #                                          | (i_all_chr == 'ABBBBBAB') | (i_all_chr == 'ABBBABBB') | (i_all_chr == 'ABAAAAAB') | (
    #                                                      i_all_chr == 'ABAAABAA') | (i_all_chr == 'ABABAABB') | (
    #                                                      i_all_chr == 'ABABBBAA'))
    #
    #                 loc_dict[i] = final_loc
    #                 all_genotype[i] = i_all_chr[final_loc]
    #
    #
    #
    #         return all_genotype,all_ratios,loc_dict
    #
    # def calculate_ppc_error(self):
    #     '''calculate parent-parent-child error'''
    #     if self._motherid != "" or self._fatherid != "":
    #         mother_gt = self._dataobject.df[self._motherid, "gtype"]
    #         father_gt = self._dataobject.df[self._fatherid, "gtype"]
    #
    #         self._dataobject.df[self._motherid, "gtype_bin"] = mother_gt.replace(conversion, inplace=False)
    #         self._dataobject.df[self._fatherid, "gtype_bin"] = father_gt.replace(conversion, inplace=False)
    #         # father_gt.replace(conversion,inplace=True)
    #         for embryo in self._embryoids:
    #             self._dataobject.df[embryo, "gtype_bin"] = self._dataobject.df[embryo, "gtype"].replace(conversion,
    #                                                                                                     inplace=False)
    #             self._dataobject.df[embryo, "ppce"] = self._dataobject.df[
    #                 [(cell, "gtype_bin") for cell in [self._motherid, self._fatherid, embryo]]].apply(func, axis=1)
    #             self._dataobject.container.append("ppce")
    #             self._dataobject.container.append("gtype_bin")
    #             # convert embryos back to gtypes
    #             # self._dataobject.df[embryo, "gtype"].replace(deconversion, inplace=True)
    #         # convert parent back to gtypes
    #         # mother_gt.replace(deconversion, inplace=True)
    #         # father_gt.replace(deconversion, inplace=True)
    #
    # def calculate_mat_pat_strand(self):
    #     if self._motherid != "" or self._fatherid != "":
    #         mother_gt = self._dataobject.df[self._motherid, "gtype"]
    #         father_gt = self._dataobject.df[self._fatherid, "gtype"]
    #
    #         for embryo in self._embryoids:
    #             x = self._dataobject.df[[(cell, "gtype") for cell in [self._motherid, self._fatherid, embryo]]].apply(
    #                 calcstrands, axis=1)
    #             self._dataobject.df[embryo, "mat_strand"] = x.values[:, 0]
    #             self._dataobject.df[embryo, "pat_strand"] = x.values[:, 1]
    #             self._dataobject.container.append("mat_strand")
    #             self._dataobject.container.append("pat_strand")
    #             # self._dataobject.df[embryo, ["mat_strand", "pat_strand"]]
    #
    # def compare_embryo_strands(self):
    #     for embryoseed in self._embryoids:
    #         seed = embryoseed
    #         for embryo in self._embryoids:
    #             if embryo != seed:
    #                 # self._dataobject.df[embryo, "mat_"+seed]=self._dataobject.df[embryo, "mat_strand"].dropna().eq(self._dataobject.df[seed, "mat_strand"].dropna()).apply(int)
    #                 self._dataobject.df[embryo, "mat_" + seed] = self._dataobject.df[
    #                     [(em, "mat_strand") for em in [embryo, seed]]].apply(cmpfunc, axis=1)
    #                 self._dataobject.df[embryo, "pat_" + seed] = self._dataobject.df[
    #                     [(em, "pat_strand") for em in [embryo, seed]]].apply(cmpfunc, axis=1)
    #                 # self._dataobject.df[embryo, "pat_" + seed] = self._dataobject.df[embryo, "pat_strand"].dropna().eq(self._dataobject.df[seed, "pat_strand"].dropna()).apply(int)
    #                 self._dataobject.container.append("mat_" + seed)
    #                 self._dataobject.container.append("pat_" + seed)
    #
    # def phase_parents(self):
    #     '''
    #     This function is to phase parents
    #     '''
    #     # phased = {}
    #     mother_gt = self._dataobject.df[self._motherid, "gtype"].copy()
    #     father_gt = self._dataobject.df[self._fatherid, "gtype"].copy()
    #     if self._paternalgrandfatherid != '' and self._paternalgrandmotherid != '':
    #         pgf_gt = self._dataobject.df[self._paternalgrandfatherid, "gtype"]
    #         pgm_gt = self._dataobject.df[self._paternalgrandmotherid, "gtype"]
    #         heter = np.where(father_gt == 'AB');
    #         f_phased = father_gt
    #         f_phased[heter[0]] = list(map(het_c, np.sum((pgf_gt[heter[0]], pgm_gt[heter[0]]), axis=0)))
    #         self._dataobject.df[self._fatherid, 'f_phased'] = f_phased
    #         self._dataobject.container.append("f_phased")
    #         print(self._familyid + 'Using paternalgrandparents info')
    #         self.pat=True
    #     # else:
    #     # print self._familyid + ": No paternalgranparents info"
    #
    #     if self._maternalgrandfatherid != '' and self._maternalgrandmotherid != '':
    #         pmf_gt = self._dataobject.df[self._maternalgrandfatherid, "gtype"]
    #         pmm_gt = self._dataobject.df[self._maternalgrandmotherid, "gtype"]
    #         heter = np.where(mother_gt == 'AB');
    #         m_phased = mother_gt
    #         m_phased[heter[0]] = list(map(het_c, np.sum((pmf_gt[heter[0]], pmm_gt[heter[0]]), axis=0)))
    #         self._dataobject.df[self._motherid, 'm_phased'] = m_phased
    #         self._dataobject.container.append("m_phased")
    #         print(self._familyid + 'Using Maternalgrandparents info')
    #         self.mat=True
    #     # else:
    #     # print self._familyid + ": No maternalgranparents info"
    #
    #
    #     if self._siblingid != '' and self._maternalgrandfatherid == '' and self._maternalgrandmotherid == '':
    #         #we have sibling - these flags indicate we can phase both maternal and paternal strand
    #         self.mat=True
    #         self.pat=True
    #         sbl_gt = self._dataobject.df[self._siblingid, 'gtype']
    #         m_phased = mother_gt.copy();
    #         f_phased = father_gt.copy()
    #         heter_f = np.where(father_gt == 'AB');
    #         heter_m = np.where(mother_gt == 'AB')
    #
    #         m_phased[heter_m[0]] = list(
    #             map(sibling, np.sum((m_phased[heter_m[0]], f_phased[heter_m[0]], sbl_gt[heter_m[0]]), axis=0)))
    #         self._dataobject.df[self._motherid, 'm_phased'] = m_phased
    #         self._dataobject.container.append("m_phased")
    #
    #
    #         print(self._familyid + ': Using siblings info')
    #
    #
    #     if self._siblingid != '' and self._paternalgrandfatherid == '' and self._paternalgrandmotherid == '':
    #         sbl_gt = self._dataobject.df[self._siblingid, 'gtype']
    #         m_phased = mother_gt.copy();
    #         f_phased = father_gt.copy()
    #         heter_f = np.where(father_gt == 'AB');
    #         heter_m = np.where(mother_gt == 'AB')
    #
    #         f_phased[heter_f[0]] = list(
    #             map(sibling, np.sum((f_phased[heter_f[0]], m_phased[heter_f[0]], sbl_gt[heter_f[0]]), axis=0)))
    #         self._dataobject.df[self._fatherid, 'f_phased'] = f_phased
    #         self._dataobject.container.append("f_phased")
    #
    #         print(self._familyid + ': Using siblings info')
    #         self.sib=True
    #
    #     if self._maternalgrandfatherid != '' and self._maternalgrandmotherid == '':
    #         m_phased = mother_gt;
    #         heter_m = np.where(mother_gt == 'AB')
    #         mgf_gt = self._dataobject.df[self._maternalgrandfatherid, "gtype"]
    #         m_phased[heter_m[0]] = list(map(f_single, mgf_gt[heter_m[0]]))
    #         self._dataobject.df[self._motherid, 'm_phased'] = m_phased
    #         print(self._familyid + ': Using single maternalgrandfather info')
    #         self.mat=True
    #
    #     if self._paternalgrandfatherid != '' and self._paternalgrandmotherid == '':
    #         f_phased = father_gt;
    #         heter_f = np.where(father_gt == 'AB')
    #         mgf_gt = self._dataobject.df[self._paternalgrandfatherid, "gtype"]
    #         f_phased[heter_f[0]] = list(map(f_single, mgf_gt[heter_f[0]]))
    #         self._dataobject.df[self._fatherid, 'f_phased'] = f_phased
    #         print(self._familyid + ': Using single paternalgrandfather info')
    #         self.pat=True
    #
    #     if self._maternalgrandmotherid != '' and self._maternalgrandfatherid == '':
    #         m_phased = mother_gt;
    #         heter_m = np.where(mother_gt == 'AB')
    #         mgm_gt = self._dataobject.df[self._maternalgrandmotherid, "gtype"]
    #         m_phased[heter_m[0]] = list(map(m_single, mgm_gt[heter_m[0]]))
    #         self._dataobject.df[self._motherid, 'm_phased'] = m_phased
    #         self._dataobject.container.append("m_phased")
    #
    #         print(self._familyid + ': Using single maternalgrandmother info')
    #         self.mat=True
    #
    #     if self._paternalgrandmotherid != '' and self._paternalgrandfatherid == '':
    #         f_phased = father_gt;
    #         heter_f = np.where(father_gt == 'AB')
    #         mgm_gt = self._dataobject.df[self._paternalgrandmotherid, "gtype"]
    #         f_phased[heter_f[0]] = list(map(m_single, mgm_gt[heter_f[0]]))
    #         self._dataobject.df[self._fatherid, 'f_phased'] = f_phased
    #         self._dataobject.container.append("f_phased")
    #         print(self._familyid + ': Using single paternalgrandmother info')
    #         self.pat=True
    #
    # def haplotype_embryos(self, seed):
    #
    #     '''
    #     This function is to phase embryos by using one embryo as a seed
    #     '''
    #     new_val_col_mat='Hap_strand:Mat_alg:{}_seed:{}'.format('IBD-SEED', seed)
    #     new_val_col_pat='Hap_strand:Pat_alg:{}_seed:{}'.format('IBD-SEED', seed)
    #
    #     mother_gt = self._dataobject.df[self._motherid, "gtype"].copy();
    #     father_gt = self._dataobject.df[self._fatherid, "gtype"].copy()
    #     self._siblingid = seed;
    #     print('Seed', self._siblingid);
    #     sbl_gt = self._dataobject.df[self._siblingid, 'gtype'].copy()
    #     m_phased = mother_gt.copy();
    #     f_phased = father_gt.copy()
    #     heter_f = np.where(father_gt == 'AB');
    #     heter_m = np.where(mother_gt == 'AB')
    #     f_phased[heter_f[0]] = list(
    #         map(sibling, np.sum((f_phased[heter_f[0]], mother_gt[heter_f[0]], sbl_gt[heter_f[0]]), axis=0)))
    #     m_phased[heter_m[0]] = list(
    #         map(sibling, np.sum((m_phased[heter_m[0]], father_gt[heter_m[0]], sbl_gt[heter_m[0]]), axis=0)))
    #     self._dataobject.df[self._motherid, 'm_phased'] = m_phased
    #     self._dataobject.df[self._fatherid, 'f_phased'] = f_phased
    #     self._dataobject.container.append("m_phased")
    #     self._dataobject.container.append("f_phased")
    #
    #     df_empty = pd.DataFrame({'Name': self._dataobject.df.index.get_level_values(0),
    #                              'Chr': self._dataobject.df.index.get_level_values(1),
    #                              'Position': self._dataobject.df.index.get_level_values(2)}
    #                             )
    #
    #     df_empty= pd.DataFrame().reindex_like(self._dataobject.df)
    #     for j in self._embryoids:
    #         em_p = haplotype(self._dataobject.df[self._fatherid, 'f_phased'].copy(), mother_gt.copy(),
    #                          self._dataobject.df[j, 'gtype'].copy())
    #         em_m = haplotype(self._dataobject.df[self._motherid, 'm_phased'].copy(), father_gt.copy(),
    #                          self._dataobject.df[j, 'gtype'].copy())
    #
    #         df_empty[j,new_val_col_mat] = list(map(int, em_m))
    #         df_empty[j,new_val_col_pat] = list(map(int, em_p))
    #
    #         self._dataobject.df[j, new_val_col_mat] = list(map(int, em_m))
    #         self._dataobject.df[j, new_val_col_pat] = list(map(int, em_p))
    #
    #
    #     #df_empty.to_csv(
    #     #    'Seed{}_{}_hom_{}_het_{}_phasing.hap'.format(str(seed), str(alg), str(hom_threshold), str(het_threshold)),
    #     #    header=True, sep='\t', index=None)
    #     return df_empty
    #
    # def get_haplotype(self):
    #     '''
    #     This function is to phase all the embryos by using phased parents
    #
    #     Return: dataframe of all phased embryos
    #     '''
    #     #df_empty = pd.DataFrame({'Name': self._dataobject.df.index.get_level_values(0),
    #     #                         'Chr': self._dataobject.df.index.get_level_values(1),
    #     #                         'Position': self._dataobject.df.index.get_level_values(2)}
    #     df_empty=pd.DataFrame().reindex_like(self._dataobject.df)
    #
    #     all_items = list(self._dataobject.df.stack(level=0))
    #
    #     new_val_col_mat='Hap_strand:{}_alg:{}_seed:{}'.format('Mat', 'IBD-GRAND', 'parent')
    #     new_val_col_pat='Hap_strand:{}_alg:{}_seed:{}'.format('Pat', 'IBD-GRAND', 'parent')
    #     if 'm_phased' in all_items:
    #         phased_par = self._dataobject.df[self._motherid, 'm_phased']
    #         ot_par = self._dataobject.df[self._fatherid, 'gtype']
    #         for embryo in self._embryoids:
    #             em_gt = self._dataobject.df[embryo, 'gtype']
    #             hap_em = haplotype(phased_par, ot_par, em_gt)
    #             # print len(hap_em),(df_empty.values).shape
    #             df_empty[embryo,new_val_col_mat] = list(map(int, hap_em))
    #             self._dataobject.df[embryo,new_val_col_mat]=df_empty[embryo,new_val_col_mat]
    #             self._dataobject.container.append(new_val_col_mat)
    #     if 'f_phased' in all_items:
    #         phased_par = self._dataobject.df[self._fatherid, 'f_phased']
    #         ot_par = self._dataobject.df[self._motherid, 'gtype']
    #         for embryo in self._embryoids:
    #             em_gt = self._dataobject.df[embryo, 'gtype']
    #             hap_em = haplotype(phased_par, ot_par, em_gt)
    #             df_empty[embryo,new_val_col_pat] = list(map(int, hap_em))
    #             self._dataobject.df[embryo,new_val_col_pat]=df_empty[embryo,new_val_col_pat]
    #             self._dataobject.container.append(new_val_col_pat)
    #     df_empty = df_empty.sort_values(by=['Chr', 'Position'])
    #     return df_empty

