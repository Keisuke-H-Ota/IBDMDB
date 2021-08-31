# mmvec 用の入力データを作成する

import os
import pandas as pd
import numpy as np
import random
import streamlit as st
import altair as alt
from altair_saver import save
from streamlit_vega_lite import vega_lite_component, altair_component
from altair import limit_rows, to_values
import toolz
t = lambda data: toolz.curried.pipe(data, limit_rows(max_rows=100000000), to_values)
alt.data_transformers.register('custom', t)
alt.data_transformers.enable('custom')

# メタデータ（md）ファイルの読み取り
cwd = os.getcwd()
gwd = cwd.replace('/InputData','')
md_path = ['RawData','MetaData','hmp2_metadata.csv']
md_file = os.path.join(gwd, *md_path)
md = pd.read_table(md_file, delimiter=',') # データ型を指定できずにエラーが出るが、使わない行なので無視してよい

# メタボロミクス（mbx）ファイルの読み取り
mbx_path = ['RawData','Stool','MBX','HMP2_metabolomics.csv.gz']
mbx_file = os.path.join(gwd, *mbx_path)
mbx = pd.read_table(mbx_file, delimiter=',', index_col=6) # データ型を指定できずにエラーが出るが使わない行なので無視

# mbx データの集計
mbx = mbx.dropna(subset=['Metabolite']).groupby('Metabolite').sum()
mbx = mbx.iloc[:,3:].T
mbx['External ID'] = mbx.index

# site sub coll を抜き出す
mbx_md = md[md['data_type']=='metabolomics'][['External ID','site_sub_coll']]
mbx = pd.merge(mbx_md,mbx)
mbx = mbx.set_index('site_sub_coll')
mbx = mbx.drop('External ID',axis=1)
# print(mbx)

# MGX データのメタデータを取り出す
mgx_md = md[md['data_type']=='metagenomics'][['External ID','site_sub_coll']]

# MGX データ（細菌）の読取
Stool_Microbes_path = ['RawData','Stool','MGX','taxonomic_profiles_3.tsv.gz']
Stool_Microbes_file = os.path.join(gwd, *Stool_Microbes_path)
Stool_Microbes_df = pd.read_table(Stool_Microbes_file, sep='\t' ,header=0, index_col=0)
Stool_Microbes_df.columns = list(map(lambda x:x[:-8],Stool_Microbes_df.columns.tolist()))

# 分類階級ごとに集計
L_lst = ['p','c','o','f','g','s']
Level_lst = ['phylum','class','order','family','genus','species']
tax_lst = []
for i in Stool_Microbes_df.index.tolist():
	tax_lst.append(i[(i.rfind('__'))-1])
Stool_Microbes_df['Level'] = tax_lst

input_data = pd.DataFrame()
for L,Level in zip(L_lst,Level_lst):
	# print('################ {0} {1} ################'.format(L,Level))
	mgx = Stool_Microbes_df[Stool_Microbes_df['Level'].isin(['W',L])].drop('Level',axis=1)
	microbes_lst = sorted(list(set(mgx.index.tolist())))
	microbes_lst = [microbes_lst[0]]+list(map(lambda x:x[x[:x.rfind('__')].rfind('__')-1:],microbes_lst[1:]))
	mgx.index = microbes_lst
	mgx = mgx.T
	mgx['External ID'] = mgx.index

	# site sub coll を抜き出す
	mgx = pd.merge(mgx_md,mgx)
	mgx = mgx.set_index('site_sub_coll')
	mgx = mgx.drop('External ID',axis=1)
	mgx = mgx.loc[~mgx.index.duplicated(),:]

	# 両方に共通するサンプル名を取得
	sn = list(set(mgx.index)&set(mbx.index))
	# 生データにカンマが含まれるので tsv で保存する。ここで保存される tsv は mmvec の入力データとなる
	tmp_df = mgx.T[sn]
	tmp_df.to_csv('mmvec/MGX-MBX/input/{}.tsv'.format(Level), header=True, index=True, sep='\t')

	# 相関係数を計算し、可視化する生データを用意する
	tmp_df['genre'] = Level
	# print(tmp_df)
	input_data = pd.concat([input_data,tmp_df],axis=0)
	del tmp_df

# 生データにカンマが含まれるので tsv で保存する
tmp_df = mbx.T[sn]
tmp_df.to_csv('mmvec/MGX-MBX/input/mbx.tsv', header=True, index=True, sep='\t')
# 相関係数を計算し、可視化する生データを用意する
tmp_df['genre'] = 'metabolites'
# print(tmp_df)
input_data = pd.concat([input_data,tmp_df],axis=0)
input_data.index = input_data.index.set_names('feature')

# 入力データを保存
input_data.to_csv('mmvec/MGX-MBX/input/MGX-MBX.tsv', header=True, index=True, sep='\t')
