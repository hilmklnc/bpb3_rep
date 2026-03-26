#import stored_tfs as stf
from bpb3.script.script_high.other import stored_tfs as stf
#import re
import pandas as pd
import glob
import os
import numpy as np
import argparse
#exec(open("check_accuracy.py" ).read())

def my_parser(parser):
   required= parser.add_argument_group("Required ")
   required.add_argument("--in_bpb3_ranking_path", help="File path for bpb3 ranking results, such as ../data/out/sep8_run1/test_67snps/foreground/test_seq/rankings/",
                          required=True, type=str, metavar="FILE PATH")   
   
   optional= parser.add_argument_group("Option parameters , have default values")
   optional.add_argument("--in_topRank_cutoff", help="Top rank cutof value for selected TFs for comparison, default=10 ",
                       type=int, metavar="NUMBER", default=10)
   return parser


def load_results(in_files, rankCutoff=0):
 #load all snp to dataframe based on bpb3 ranking file
 record_snp2rankings={}
 for fi in in_files:
    tmp_df0=pd.read_csv(fi,sep='\t',skiprows=[0])
    if rankCutoff> 0:
      #print('Ranked results selected by top , ', rankCutoff)
      #do filter data based on cutoff
      #assume the second column is 'rank' order from bpb3 ranking results
      #tmp_df0.iloc[:,1]=pd.to_numeric(tmp_df0.iloc[:,1])
      #tmp_df=tmp_df0[tmp_df0.iloc[:,1]<=rankCutoff].copy()
      #or use column name
      tmp_df0['rank']=pd.to_numeric(tmp_df0['rank'])
      tmp_df=tmp_df0[tmp_df0['rank']<=rankCutoff].copy()
    else:
      #print('Input all ranked results')
      tmp_df=tmp_df0.copy()

    infile=open(fi,'r')
    firstline=infile.readline()
    tmp_snp=firstline.split('\t')[-1].strip()

    record_snp2rankings[tmp_snp]=tmp_df.copy()
 return record_snp2rankings

def Intersection(lst1, lst2):
    return set(lst1).intersection(lst2)

def isTfCorrect(tfs, x, tfNameMap):
  out_TF=False
  for tf in tfs:
    if tf not in tfNameMap.keys():
       print('Unknown TF: ', tf)
       return
    for tfName in tfNameMap[tf]:
      #if re.search(tfName, x, re.IGNORECASE):
      if  isinstance(x,list):
        for xx in x:
          #print(xx,x)
          if tfName.lower() in xx.lower():
            out_TF= True
      else:
          if tfName.lower() in x.lower():
            out_TF=True
  return out_TF

def find_best_rank(record_snp2rankings_fi,fi,file_name,tfs, change,tfNameMap,out_df_columns):
  '''first look for the correc change, if not find then look for the remaining results with different changes'''
  tmp_out=record_snp2rankings_fi[record_snp2rankings_fi[file_name].apply(lambda x: isTfCorrect(tfs,x,tfNameMap))].copy()
  if 'gain' in change.lower():
        out=tmp_out[tmp_out.change=='Positive'].copy()
        out_rev=tmp_out[tmp_out.change=='Negative'].copy()
  elif 'loss' in change.lower():
        out=tmp_out[tmp_out.change=='Negative'].copy()
        out_rev=tmp_out[tmp_out.change=='Positive'].copy()
  else:
        #print(change)
        out=tmp_out.copy()
        out_rev=pd.DataFrame()

  if out.shape[0]>0:
       idx2min=out['rank'].idxmin()
       #print(out.loc[idx2min,:])
       if len(out_df_columns)==3:
         new_row={'mutation': fi ,'bestRank': out.loc[idx2min,['rank']].to_list()[0] ,
                 'bestScore': out.loc[idx2min,['score_pca']].to_list()[0]}
       else:
         #print(out.loc[idx2min,['mlps']].to_list())
         new_row={'mutation': fi ,'bestRank': out.loc[idx2min,['rank']].to_list()[0] ,
                 'bestScore': out.loc[idx2min,['score_pca']].to_list()[0],
                  'file_name': out.loc[idx2min,['file_name']].to_list()[0],
                   'mlps': out.loc[idx2min,['mlps']].to_list()[0]}
 
  elif out.shape[0]==0 and out_rev.shape[0]>0:
       idx2min=out_rev['rank'].idxmin()
       if len(out_df_columns)==3:
         new_row={'mutation': fi ,'bestRank': out_rev.loc[idx2min,['rank']].to_list()[0] ,
                 'bestScore': out_rev.loc[idx2min,['score_pca']].to_list()[0]}
       else:
         new_row={'mutation': fi ,'bestRank': out_rev.loc[idx2min,['rank']].to_list()[0] ,
                 'bestScore': out_rev.loc[idx2min,['score_pca']].to_list()[0],
                  'file_name': out_rev.loc[idx2min,['file_name']].to_list()[0],
                   'mlps': out_rev.loc[idx2min,['mlps']].to_list()[0]}
  else:
       if len(out_df_columns)==3:
         new_row={'mutation':fi, 'bestRank': np.NaN, 'bestScore': np.NaN }
       else:
         new_row={'mutation':fi, 'bestRank': np.NaN, 'bestScore': np.NaN , 'file_name': '','mlps':''}
  return new_row

def get_tfs_and_change(fi):
    correct_tfs=fi.split('|')
    correct_tf=correct_tfs[0:-1]
    change=correct_tfs[-1]
    tfs=[]
    for ii in correct_tf:
       tfs.append(ii.split('_')[1])
    return tfs, change 

def best_rank_in_result(stf_snps,tfNameMap,record_snp2rankings,out_df):
  #out_df=pd.DataFrame(columns=['mutation','bestRank','bestScore'])
  out_df_list=[]
  for fi in stf_snps:
    #print(fi)
    tfs, change= get_tfs_and_change(fi)

    new_row=find_best_rank(record_snp2rankings[fi],fi,'file_name',tfs, change,tfNameMap,out_df.columns)
    #jbw 2026
    #out_df=out_df.append(new_row,ignore_index=True)
    out_df_list.append(new_row.copy())
  #jbw 2026
  out_df=pd.DataFrame(out_df_list)
  return out_df.copy()

def accuracy(out_df, top_rank_cutoff):
  total_df=out_df.shape[0]
  print(total_df)
  out2df=out_df[out_df.bestRank<=top_rank_cutoff].copy()
  fd_num= out2df.shape[0]
  accuracy=fd_num/total_df *100
  median_rank=out2df.bestRank.median()
  mean_rank=out2df.bestRank.mean()
  print(fd_num, ' of ', total_df, ' are correct (', accuracy, '% )')
  print('Mean rank: ', mean_rank, ' Median rank:', median_rank)
  return total_df, fd_num, out2df

def run(args):
  in_result_folder=args.in_bpb3_ranking_path

  in_files=glob.glob(os.path.join(in_result_folder,'*ranking.tsv'))

  total_snps=len(in_files)
  print('Total :', total_snps)

  #load all snp to dataframe
  record_snp2rankings=load_results(in_files)

  #check accuracyi in fredriksson snps
  out_df=pd.DataFrame(columns=['mutation','bestRank','bestScore'])
  fd_out_df= best_rank_in_result(stf.fredrikssonMuts,stf.tfNameMap,record_snp2rankings,out_df)
  hg_out_df= best_rank_in_result(stf.hgmdMuts, stf.tfNameMap, record_snp2rankings,out_df)
  ad_out_df= best_rank_in_result(stf.andersonEpsteinMuts,stf.tfNameMap, record_snp2rankings,out_df)

  fd_out_df.bestRank=pd.to_numeric(fd_out_df.bestRank)
  hg_out_df.bestRank=pd.to_numeric(hg_out_df.bestRank)
  ad_out_df.bestRank=pd.to_numeric(ad_out_df.bestRank)

  #calculate accuracy
  top_rank_cutoff=args.in_topRank_cutoff
  print('Considering top ', top_rank_cutoff, ' TFs , change direction: False')

  print('fredrikssonMuts')
  fd_total_df, fd_num, fd_out2df=accuracy(fd_out_df,top_rank_cutoff)

  print('hgmdMuts')
  hg_total_df, hg_num, hg_out2df=accuracy(hg_out_df,top_rank_cutoff)

  print('andersonEpsteinMuts')
  ad_total_df, ad_num, ad_out2df=accuracy(ad_out_df,top_rank_cutoff)

  print('\n')
  total_mut=fd_total_df+hg_total_df+ad_total_df
  good_num=fd_num+ hg_num+ ad_num
  percent=good_num/total_mut*100

  all_df=pd.concat([fd_out_df,hg_out_df,ad_out_df])
  #replace low rank as NA
  all_df.loc[all_df.bestRank>top_rank_cutoff,['bestRank','bestScore']]=np.nan
  print('Total: ',good_num, ' of ',total_mut, ' are correct (',percent, '% )')
  print('All: mean rank=',all_df.bestRank.mean(), ' median rank=',all_df.bestRank.median())

  all_df.bestRank=all_df.bestRank.fillna('NA')
  all_df.bestScore=all_df.bestScore.fillna('NA')
  out_file=os.path.join(in_result_folder,'../ranks.tsv')
  all_df.to_csv(out_file,sep='\t',index=False)
  print('Export at: ', out_file)
  return all_df, fd_total_df, hg_total_df, ad_total_df


####################################
if __name__=='__main__':
 args= my_parser(argparse.ArgumentParser('python check_accuracy.py')).parse_args()
 run(args)


