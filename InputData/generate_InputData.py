import os
import numpy as np
import pandas as pd
import numpy as np

L_lst = ['p','c','o','f','g','s']
Level_lst = ['phylum','class','order','family','genus','species']

# メタデータ読取・保存

cwd = os.getcwd()
gwd = cwd.replace('/InputData','')
md_path = ['RawData','MetaData','hmp2_metadata.csv']
md_file = os.path.join(gwd, *md_path)
md_df = pd.read_table(md_file, delimiter=',' ,sep=',' ,header=0, index_col=0) # データ型を指定できずにエラーが出るが、使わない行なので無視してよい
md_df.to_csv(os.path.join(cwd, *['MetaData','hmp2_metadata.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

def create_metadata(md_df, data_type):
	md_df = md_df[md_df['data_type']==data_type]
	md_df = md_df.set_index('External ID')
	return md_df

mgx_md_df = create_metadata(md_df, 'metagenomics')
# mgx_md_df.to_csv(os.path.join(cwd, *['MetaData','mgx_md_df.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')
mbx_md_df = create_metadata(md_df, 'metabolomics')
# mbx_md_df.to_csv(os.path.join(cwd, *['MetaData','mbx_md_df.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')
mpx_md_df = create_metadata(md_df, 'proteomics')
# mpx_md_df.to_csv(os.path.join(cwd, *['MetaData','mpx_md_df.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')
vx_md_df = create_metadata(md_df, 'viromics')
# vx_md_df.to_csv(os.path.join(cwd, *['MetaData','vx_md_df.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')
participant_lst = mgx_md_df['Participant ID'].values.tolist() + mbx_md_df['Participant ID'].values.tolist() + mpx_md_df['Participant ID'].values.tolist() + vx_md_df['Participant ID'].values.tolist()
participant_lst = sorted(list(set(participant_lst)))
pd.DataFrame(participant_lst).to_csv(os.path.join(cwd, *['MetaData','participant.tsv']), header = True, index = True, sep = '\t')

# MGX データ（細菌）の読取
Stool_Microbes_path = ['RawData','Stool','MGX','taxonomic_profiles_3.tsv.gz']
Stool_Microbes_file = os.path.join(gwd, *Stool_Microbes_path)
Stool_Microbes_df = pd.read_table(Stool_Microbes_file, sep='\t' ,header=0, index_col=0)
Stool_Microbes_df.columns = list(map(lambda x:x[:-8],Stool_Microbes_df.columns.tolist()))

# 分類階級ごとに集計
tax_lst = []
for i in Stool_Microbes_df.index.tolist():
	tax_lst.append(i[(i.rfind('__'))-1])
Stool_Microbes_df['Level'] = tax_lst
Stool_Microbes_all_df = pd.DataFrame()
for L,Level in zip(L_lst,Level_lst):
	tmp_df = Stool_Microbes_df[Stool_Microbes_df['Level'].isin(['W',L])].drop('Level',axis=1)
	microbes_lst = sorted(list(set(tmp_df.index.tolist())))
	microbes_lst = [microbes_lst[0]]+list(map(lambda x:x[x[:x.rfind('__')].rfind('__')-1:],microbes_lst[1:]))
	tmp_df.index = microbes_lst

	# メタデータと MGX データ（細菌）の結合
	tmp_df = tmp_df.T.join(mgx_md_df[['Participant ID','week_num']])
	tmp_df = tmp_df.groupby(['Participant ID','week_num']).mean().reset_index()
	tmp_df = pd.melt(tmp_df,id_vars=['Participant ID','week_num'],var_name="Microbes",value_name="Composition")
	tmp_df['Level'] = Level
	Stool_Microbes_all_df = pd.concat([Stool_Microbes_all_df,tmp_df])
Stool_Microbes_all_df = Stool_Microbes_all_df.set_index('Participant ID')
Stool_Microbes_all_df.to_csv(os.path.join(cwd, *['Stool','MGX','Stool_MGX_Microbes.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

# 代謝経路データ読取
Stool_PWY_path = ['RawData','Stool','MGX','pathabundances_3.tsv.gz']
Stool_PWY_file = os.path.join(gwd, *Stool_PWY_path)
Stool_PWY_df = pd.read_table(Stool_PWY_file, sep='\t' ,header=0, index_col=0)
Stool_PWY_df.columns = list(map(lambda x:x[:-18],Stool_PWY_df.columns.tolist()))

## 代謝経路の集計データを抽出
PWY_agg_lst = [p for p in Stool_PWY_df.index.tolist() if not '|' in p][2:]
Stool_PWY_agg_df = Stool_PWY_df.T[PWY_agg_lst]
PWY_agg_lst = sorted(list(map(lambda x:x[x.rfind(':')+1:],PWY_agg_lst)))
Stool_PWY_agg_df.columns = PWY_agg_lst
pd.DataFrame(PWY_agg_lst).to_csv(os.path.join(cwd, *['Stool','MGX','Stool_MGX_PWY_lst.tsv']), header = True, index = True, sep = '\t')

## メタデータと代謝経路（集計）データの結合
Stool_PWY_agg_df = Stool_PWY_agg_df.join(mgx_md_df[['Participant ID','week_num']])
Stool_PWY_agg_df = Stool_PWY_agg_df.groupby(['Participant ID','week_num']).mean().reset_index()
Stool_PWY_agg_df = pd.melt(Stool_PWY_agg_df,id_vars=['Participant ID','week_num'],var_name="PWY",value_name="RPKs(CPM)")
Stool_PWY_agg_df = Stool_PWY_agg_df.set_index('Participant ID')
Stool_PWY_agg_df.to_csv(os.path.join(cwd, *['Stool','MGX','Stool_MGX_PWY.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

# メタボロームデータの読取
Stool_mbx_path = ['RawData','Stool','MBX','HMP2_metabolomics.csv.gz']
Stool_mbx_file = os.path.join(gwd, *Stool_mbx_path)
Stool_mbx_df = pd.read_table(Stool_mbx_file, delimiter=',', index_col=6) # データ型を指定できずにエラーが出るが使わない行なので無視
Stool_mbx_df = Stool_mbx_df.dropna(subset=['Metabolite']).groupby('Metabolite').sum()
Stool_mbx_df = Stool_mbx_df.iloc[:,3:].T
## メタデータとの結合
Stool_mbx_df = Stool_mbx_df.join(mbx_md_df[['Participant ID','week_num']])
Stool_mbx_df = Stool_mbx_df.groupby(['Participant ID','week_num']).mean().reset_index()
Stool_mbx_df = pd.melt(Stool_mbx_df,id_vars=['Participant ID','week_num'],var_name="Metabolites",value_name="Abundance")
Stool_mbx_df = Stool_mbx_df.set_index('Participant ID')
Stool_mbx_df.to_csv(os.path.join(cwd, *['Stool','MBX','Stool_MBX.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

# プロテオームデータの読取
Stool_mpx_path = ['RawData','Stool','MPX','HMP2_proteomics_ecs.tsv']
Stool_mpx_file = os.path.join(gwd, *Stool_mpx_path)
Stool_mpx_df = pd.read_table(Stool_mpx_file, delimiter='\t', index_col=0).T # データ型を指定できずにエラーが出るが使わない行なので無視
Stool_mpx_df.columns = list(map(lambda x:x[x.rfind(':')+1:],Stool_mpx_df.columns.tolist()))
Stool_mpx_df = Stool_mpx_df.drop('UNGROUPED',axis=1)
## メタデータとの結合
Stool_mpx_df = Stool_mpx_df.join(mpx_md_df[['Participant ID','week_num']])
Stool_mpx_df = Stool_mpx_df.groupby(['Participant ID','week_num']).mean().reset_index()
Stool_mpx_df = pd.melt(Stool_mpx_df,id_vars=['Participant ID','week_num'],var_name="Gene",value_name="Abundance")
Stool_mpx_df = Stool_mpx_df.set_index('Participant ID')
Stool_mpx_df.to_csv(os.path.join(cwd, *['Stool','MPX','Stool_MPX.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

# ウイルスデータの読取
Stool_vx_path = ['RawData','Stool','VX','taxonomic_profiles.tsv']
Stool_vx_file = os.path.join(gwd, *Stool_vx_path)
Stool_Viruses_df = pd.read_table(Stool_vx_file, sep='\t' ,header=0, index_col=0)

# 分類階級の選択（p,c,o,g,s から選択）
Level_lst = ['species','genus','family','order','class','phylum']
taxon_lst = Stool_Viruses_df.index.tolist()
for Level in Level_lst:
	tax_lst = []
	tmp_lst = []
	for i in taxon_lst:
	    tax_lst.append(i[i.rfind('__')+2:])
	    tmp_lst.append(i[:i.rfind('__')-2])
	taxon_lst = tmp_lst
	Stool_Viruses_df[Level] = tax_lst

Stool_Viruses_df = pd.melt(Stool_Viruses_df,id_vars=Level_lst,var_name="Participant ID",value_name="Composition")
Stool_Viruses_df = Stool_Viruses_df.set_index('Participant ID')
Stool_Viruses_df = Stool_Viruses_df.join(vx_md_df[['Participant ID','week_num']])
Stool_Viruses_df = pd.melt(Stool_Viruses_df,id_vars=['Participant ID','week_num','Composition'],var_name="Level",value_name="Virus")
Stool_Viruses_df = Stool_Viruses_df.set_index('Participant ID')
Stool_Viruses_df.to_csv(os.path.join(cwd, *['Stool','VX','Stool_VX.tsv.gz']), header = True, index = True, sep = '\t', compression='gzip')

quit()

# キャッシュ受取
mgx_md_df,mbx_md_df,mpx_md_df,vx_md_df,participant_lst = create_metadata_and_participant_lst()
Stool_Microbes_all_df = create_Stool_Microbes_all_df(mgx_md_df)
PWY_agg_lst, Stool_PWY_agg_df = create_Stool_PWY_agg_df(mgx_md_df)
Stool_mbx_df = create_Stool_mbx_df(mbx_md_df)
Stool_mpx_df = create_Stool_mpx_df(mpx_md_df)
Stool_Viruses_df = create_Stool_Viruses_df(vx_md_df)
participant = st.selectbox('参加者を選んでください',participant_lst)

# 細菌分布の可視化
Level_radio = alt.binding_radio(options=Level_lst)
Level_select = alt.selection_single(fields=['Level'], bind=Level_radio, name="Level", init={'Level': Level_lst[0]})

microbes_select = alt.selection_multi(empty='all', fields=['Microbes'])

bars = alt.Chart().mark_bar().encode(
    x='week_num:O',
    y='sum(Composition)',
    color=alt.condition(microbes_select,
                        'Microbes',
                        alt.value('lightgray'),
                        legend=None
                        ),
    tooltip=['Microbes','Composition']
).properties(
    width=500, height=500
).add_selection(
    microbes_select
)

# Create a selection that chooses the nearest point & selects based on x-value
nearest = alt.selection(type='single', nearest=True, on='mouseover',
                        fields=['week_num'], empty='none')

# The basic line
line = alt.Chart().mark_line().encode(
    x='week_num:Q',
    y='Composition:Q',
    color = 'Microbes'
)

# Transparent selectors across the chart. This is what tells us
# the x-value of the cursor
selectors = alt.Chart().mark_point().encode(
    x='week_num:Q',
    opacity=alt.value(0),
).add_selection(
    nearest
)

# Draw points on the line, and highlight based on selection
points = line.mark_point().encode(
    opacity=alt.condition(nearest, alt.value(1), alt.value(0))
)

# Draw text labels near the points, and highlight based on selection
text = line.mark_text(align='left', dx=5, dy=-5).encode(
    text=alt.condition(nearest, 'Composition:Q', alt.value(' '))
)

# Draw a rule at the location of the selection
rules = alt.Chart().mark_rule(color='gray').encode(
    x='week_num:Q',
).transform_filter(
    nearest
)

# Put the five layers into a chart and bind the data
lines = alt.layer(
    line, selectors, points, rules, text
).properties(
    width=500, height=500
).transform_filter(
    microbes_select
)

Stool_Microbes_chart = alt.hconcat(bars, lines, data = Stool_Microbes_all_df[Stool_Microbes_all_df['Participant ID']==participant]).add_selection(
    Level_select
).transform_filter(
    Level_select
)

# 代謝酵素分布の可視化

PWY_dropdown = alt.binding_select(options=PWY_agg_lst)
PWY_select = alt.selection_single(fields=['PWY'], bind=PWY_dropdown, name="microbes", init={'PWY': PWY_agg_lst[0]})

# Create a selection that chooses the nearest point & selects based on x-value
nearest = alt.selection(type='single', nearest=True, on='mouseover',
                        fields=['week_num'], empty='none')

# The basic line
line = alt.Chart().mark_line().encode(
    x='week_num:Q',
    y='RPKs(CPM):Q',
    color='PWY'
)

# Transparent selectors across the chart. This is what tells us
# the x-value of the cursor
selectors = alt.Chart().mark_point().encode(
    x='week_num:Q',
    opacity=alt.value(0),
).add_selection(
    nearest
)

# Draw points on the line, and highlight based on selection
points = line.mark_point().encode(
    opacity=alt.condition(nearest, alt.value(1), alt.value(0))
)

# Draw text labels near the points, and highlight based on selection
text = line.mark_text(align='left', dx=5, dy=-5).encode(
    text=alt.condition(nearest, 'RPKs(CPM):Q', alt.value(' '))
)

# Draw a rule at the location of the selection
rules = alt.Chart().mark_rule(color='gray').encode(
    x='week_num:Q',
).transform_filter(
    nearest
)

# Put the five layers into a chart and bind the data
lines = alt.layer(
    line, selectors, points, rules, text,
).properties(
    width=500, height=500
).add_selection(
    PWY_select
).transform_filter(
    PWY_select
)

PWY_selector = alt.selection_multi(empty='all', fields=['PWY'])

bars = alt.Chart().mark_bar().encode(
    x='week_num:O',
    y='sum(RPKs(CPM))',
    color=alt.condition(PWY_selector,
                        'PWY',
                        alt.value('lightgray'),
                        legend=None
                        ),
    tooltip=['PWY','RPKs(CPM)','week_num']
).properties(
    width=500, height=500
).add_selection(
    PWY_selector
)

Stool_PWY_agg_chart = alt.hconcat(bars, lines, data = Stool_PWY_agg_df[Stool_PWY_agg_df['Participant ID'] == participant])

# メタボローム分布可視化

metabolites_select = alt.selection_multi(empty='all', fields=['Metabolites'])

bars = alt.Chart().mark_bar().encode(
    x='week_num:O',
    y='sum(Abundance)',
    color=alt.condition(metabolites_select,
                        'Metabolites',
                        alt.value('lightgray'),
                        legend=None
                        ),
    tooltip=['Metabolites','Abundance']
).properties(
    width=500, height=500
).add_selection(
    metabolites_select
)

# Create a selection that chooses the nearest point & selects based on x-value
nearest = alt.selection(type='single', nearest=True, on='mouseover',
                        fields=['week_num'], empty='none')

# The basic line
line = alt.Chart().mark_line().encode(
    x='week_num:Q',
    y='Abundance:Q',
    color = 'Metabolites'
)

# Transparent selectors across the chart. This is what tells us
# the x-value of the cursor
selectors = alt.Chart().mark_point().encode(
    x='week_num:Q',
    opacity=alt.value(0),
).add_selection(
    nearest
)

# Draw points on the line, and highlight based on selection
points = line.mark_point().encode(
    opacity=alt.condition(nearest, alt.value(1), alt.value(0))
)

# Draw text labels near the points, and highlight based on selection
text = line.mark_text(align='left', dx=5, dy=-5).encode(
    text=alt.condition(nearest, 'Abundance:Q', alt.value(' '))
)

# Draw a rule at the location of the selection
rules = alt.Chart().mark_rule(color='gray').encode(
    x='week_num:Q',
).transform_filter(
    nearest
)

# Put the five layers into a chart and bind the data
lines = alt.layer(
    line, selectors, points, rules, text
).properties(
    width=500, height=500
).transform_filter(
    metabolites_select
)

Stool_mbx_chart = alt.hconcat(bars, lines, data = Stool_mbx_df[Stool_mbx_df['Participant ID'] == participant])

# プロテオーム分布可視化
gene_select = alt.selection_multi(empty='all', fields=['Gene'])

bars = alt.Chart().mark_bar().encode(
    x='week_num:O',
    y='sum(Abundance)',
    color=alt.condition(gene_select,
                        'Gene',
                        alt.value('lightgray'),
                        legend=None
                        ),
    tooltip=['Gene','Abundance']
).properties(
    width=500, height=500
).add_selection(
    gene_select
)

# Create a selection that chooses the nearest point & selects based on x-value
nearest = alt.selection(type='single', nearest=True, on='mouseover',
                        fields=['week_num'], empty='none')

# The basic line
line = alt.Chart().mark_line().encode(
    x='week_num:Q',
    y='Abundance:Q',
    color = 'Gene'
)

# Transparent selectors across the chart. This is what tells us
# the x-value of the cursor
selectors = alt.Chart().mark_point().encode(
    x='week_num:Q',
    opacity=alt.value(0),
).add_selection(
    nearest
)

# Draw points on the line, and highlight based on selection
points = line.mark_point().encode(
    opacity=alt.condition(nearest, alt.value(1), alt.value(0))
)

# Draw text labels near the points, and highlight based on selection
text = line.mark_text(align='left', dx=5, dy=-5).encode(
    text=alt.condition(nearest, 'Abundance:Q', alt.value(' '))
)

# Draw a rule at the location of the selection
rules = alt.Chart().mark_rule(color='gray').encode(
    x='week_num:Q',
).transform_filter(
    nearest
)

# Put the five layers into a chart and bind the data
lines = alt.layer(
    line, selectors, points, rules, text
).properties(
    width=500, height=500
).transform_filter(
    gene_select
)

Stool_mpx_chart = alt.hconcat(bars, lines, data = Stool_mpx_df[Stool_mpx_df['Participant ID'] == participant])

# ウイルス分布可視化
Level_radio = alt.binding_radio(options=Level_lst)
Level_select = alt.selection_single(fields=['Level'], bind=Level_radio, name="Level", init={'Level': Level_lst[0]})

viruses_select = alt.selection_multi(empty='all', fields=['Virus'])

bars = alt.Chart().mark_bar().encode(
    x='week_num:O',
    y='sum(Composition)',
    color=alt.condition(viruses_select,
                        'Virus',
                        alt.value('lightgray'),
                        legend=None
                        ),
    tooltip=['Virus','Composition']
).properties(
    width=500, height=500
).add_selection(
    viruses_select
)

# Create a selection that chooses the nearest point & selects based on x-value
nearest = alt.selection(type='single', nearest=True, on='mouseover',
                        fields=['week_num'], empty='none')

# The basic line
line = alt.Chart().mark_line().encode(
    x='week_num:Q',
    y='Composition:Q',
    color = 'Virus'
)

# Transparent selectors across the chart. This is what tells us
# the x-value of the cursor
selectors = alt.Chart().mark_point().encode(
    x='week_num:Q',
    opacity=alt.value(0),
).add_selection(
    nearest
)

# Draw points on the line, and highlight based on selection
points = line.mark_point().encode(
    opacity=alt.condition(nearest, alt.value(1), alt.value(0))
)

# Draw text labels near the points, and highlight based on selection
text = line.mark_text(align='left', dx=5, dy=-5).encode(
    text=alt.condition(nearest, 'Composition:Q', alt.value(' '))
)

# Draw a rule at the location of the selection
rules = alt.Chart().mark_rule(color='gray').encode(
    x='week_num:Q',
).transform_filter(
    nearest
)

# Put the five layers into a chart and bind the data
lines = alt.layer(
    line, selectors, points, rules, text
).properties(
    width=500, height=500
).transform_filter(
    viruses_select
)

Stool_Viruses_chart = alt.hconcat(bars, lines, data = Stool_Viruses_df[Stool_Viruses_df['Participant ID'] == participant]).add_selection(
    Level_select
).transform_filter(
    Level_select
)


# 可視化
st.write(Stool_Microbes_chart)
st.write(Stool_PWY_agg_chart)
st.write(Stool_mbx_chart)
st.write(Stool_mpx_chart)
st.write(Stool_Viruses_chart)
# altair のデータサイズが膨大なため 
# /Users/ootakeisuke/.pyenv/versions/3.8.4/lib/python3.8/site-packages/streamlit/server/server_util.py で
# MESSAGE_SIZE_LIMIT を 200 MB に変更


