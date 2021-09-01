import os
import streamlit as st
import altair as alt
from altair_saver import save
import numpy as np
import pandas as pd
from skbio.stats.composition import ancom
from skbio.stats.composition import multiplicative_replacement # ANCOM の psudo-count
from collections import defaultdict
from streamlit_vega_lite import vega_lite_component, altair_component # Altair の異なる図の interactive 性をつなげる
# https://github.com/domoritz/streamlit-vega-lite
from altair import limit_rows, to_values
import toolz
t = lambda data: toolz.curried.pipe(data, limit_rows(max_rows=100000000), to_values)
alt.data_transformers.register('custom', t)
alt.data_transformers.enable('custom')

# ページ全体のレイアウトを決める。
# 参考：https://discuss.streamlit.io/t/where-to-set-page-width-when-set-into-non-widescreeen-mode/959

max_width = 1111
padding_top = 0
padding_right = 0
padding_left = 0
padding_bottom = 0

st.markdown(
        f"""
<style>
    .reportview-container .main .block-container{{
        max-width: {max_width}px;
        padding-top: {padding_top}rem;
        padding-right: {padding_right}rem;
        padding-left: {padding_left}rem;
        padding-bottom: {padding_bottom}rem;
    }}
    .reportview-container .main {{
        primaryColor="#F63366";
        backgroundColor="#FFFFFF";
        secondaryBackgroundColor="#F0F2F6";
        textColor="#262730";
        font="sans serif";
    }}
</style>
""",
        unsafe_allow_html=True,
    )


Level_lst = ['species','genus','family','order','class','phylum']

# データ読取
class Read:
    
    cwd = os.getcwd()
    gwd = cwd + '/InputData'

    @st.cache(allow_output_mutation=True, max_entries=10, ttl=3600)
    def Metadata(self):
        Metadata_path = os.path.join(self.gwd, *['MetaData','hmp2_metadata.tsv.gz'])
        Metadata_df = pd.read_table(Metadata_path, sep='\t' ,header = 0, index_col = 1)
        return Metadata_df

    @st.cache(allow_output_mutation=True, max_entries=10, ttl=3600)
    def Stool_MGX_Microbes(self):
        Stool_MGX_Microbes_path = os.path.join(self.gwd, *['Stool','MGX','Stool_MGX_Microbes.tsv.gz'])
        Stool_MGX_Microbes_df = pd.read_table(Stool_MGX_Microbes_path, sep='\t' ,header=0)
        return Stool_MGX_Microbes_df

    @st.cache(allow_output_mutation=True, max_entries=10, ttl=3600)
    def Stool_MGX_PWY(self):
        Stool_MGX_PWY_path = os.path.join(self.gwd, *['Stool','MGX','Stool_MGX_PWY.tsv.gz'])
        Stool_MGX_PWY_df = pd.read_table(Stool_MGX_PWY_path, sep='\t' ,header=0)
        return Stool_MGX_PWY_df

    @st.cache(allow_output_mutation=True, max_entries=10, ttl=3600)
    def Stool_MBX(self):
        Stool_MBX_path = os.path.join(self.gwd, *['Stool','MBX','Stool_MBX.tsv.gz'])
        Stool_MBX_df = pd.read_table(Stool_MBX_path, sep='\t' ,header=0)
        return Stool_MBX_df

    @st.cache(allow_output_mutation=True, max_entries=10, ttl=3600)
    def Stool_MGX_MBX(self):
        Stool_MGX_MBX_path = os.path.join(self.gwd, *['mmvec','MGX-MBX','input','MGX-MBX.tsv'])
        Stool_MGX_MBX_df = pd.read_table(Stool_MGX_MBX_path, delimiter='\t',index_col=0,header=0)
        return Stool_MGX_MBX_df

read = Read()
Metadata_df = read.Metadata()
Stool_MGX_Microbes_df = read.Stool_MGX_Microbes()
Stool_MGX_PWY_df = read.Stool_MGX_PWY()
Stool_MBX_df = read.Stool_MBX()
Stool_MGX_MBX_df = read.Stool_MGX_MBX()

# タイトル

st.title('炎症性腸疾患を対象にしたマルチオミクス解析')
st.subheader('Multi-omics analysis of inflammatory bowel disease')

text = '''

このアプリケーションでは、**炎症性腸疾患（Inflammatory Bowel Disease:IBD）**を対象にした公共オミクスデータベース **[IBD Multi-omics DataBase](https://ibdmdb.org/)** に登録されている糞便の**ショットガンメタゲノムデータ**および**メタボロームデータ**を可視化・解析することができます。\n

'''

st.markdown(text)


with st.expander('IBD とは'):
    text = '''
    * 慢性腸炎を引き起こす原因不明の疾患
        - 若年層の発症率が高い
        - 症状に個人差が多い

    * クローン病（CD）と潰瘍性大腸炎（UC）に分類
        - CD: 結腸や小腸末端部（回腸）を中心に断続的な病変
        - UC: 結腸や直腸を中心に連続的な病変

    * 腸内の細菌叢や代謝物の変化に伴う腸管免疫の破綻が発症・増悪に関与
        - 腸内細菌叢の乱れ：Dysbiosis
        - 慢性腸炎に伴う吸収障害によるアミノ酸の増加   
    '''
    st.markdown(text)


# 被験者の選択
st.sidebar.markdown('# 被験者の選択')

text = '''

``Suggestion``にチェックを入れると実施される nonIBD と IBD の比較解析で使用される被験者の条件を選択してください。

'''

st.sidebar.markdown(text)

Education_Level = ['Some college', 'Some high school','no degree','7th grade or less','High school graduate or GED',"Bachelor's degree","Master's degree",'Unknown/Not Reported','Professional/Doctoral degree','Associate degree']
Occupation_lst = ['Paid','Retired','Student','Unknown/Not Reported','Unpaid/volunteer']
diagnosis_lst = ['CD','UC','nonIBD']

select_sex_exec = st.sidebar.checkbox('性別でフィルタリングしますか？')
if select_sex_exec == True:
    select_sex = st.sidebar.multiselect('性別を選択してください',['Male','Female'],['Male','Female'])
    Metadata_df = Metadata_df[Metadata_df['sex'].isin(select_sex)]

select_age_exec = st.sidebar.checkbox('年齢でフィルタリングしますか？')
if select_age_exec == True:
    select_age = st.sidebar.slider('年齢を選択してください',min_value=0 , max_value=100, value=(0, 100))
    Metadata_df = Metadata_df[Metadata_df['consent_age'] >= select_age[0]]
    Metadata_df = Metadata_df[Metadata_df['consent_age'] <= select_age[1]]

select_EL_exec = st.sidebar.checkbox('学歴でフィルタリングしますか？')
if select_EL_exec == True:
    select_EL = st.sidebar.multiselect('最終学歴を選択してください',Education_Level,Education_Level)
    Metadata_df = Metadata_df[Metadata_df['Education Level'].isin(select_EL)]

select_OL_exec = st.sidebar.checkbox('職歴でフィルタリングしますか？')
if select_OL_exec == True:
    select_OL = st.sidebar.multiselect('職歴を選択してください',Occupation_lst,Occupation_lst)
    Metadata_df = Metadata_df[Metadata_df['Occupation'].isin(select_OL)]

select_diag_exec = st.sidebar.checkbox('病名でフィルタリングしますか？')
if select_diag_exec == True:
    select_diag = st.sidebar.multiselect('病名を選択してください',diagnosis_lst,diagnosis_lst)
    if len(select_diag) <= 1 :
        st.info("病名は 2 つ以上選択してください")
        st.stop()
    if not 'nonIBD' in select_diag:
        st.info("nonIBD は必ず選択してください")
        st.stop()
    Metadata_df = Metadata_df[Metadata_df['diagnosis'].isin(select_diag)]
else:
    select_diag = diagnosis_lst


text = ''' 

このプロジェクトでは**クローン病（CD）**患者 67 名、**潰瘍性大腸炎（UC）**患者 38 名、**非炎症性腸疾患（nonIBD）**患者 27 名から約 1 年間にわたって糞便試料が提供されています。
アプリケーションでは被験者ごとに約 1 年間分のメタゲノムデータおよびメタボロームデータを追跡・可視化することができます。**次のプルダウンから着目したい被験者を選択しましょう。**\n

'''
st.markdown(text)

Metadata_dist_df = Metadata_df.copy()
Metadata_dist_df = Metadata_dist_df.set_index('Participant ID')
Metadata_dist_df = Metadata_dist_df.loc[~Metadata_dist_df.index.duplicated(),:]
Participant_lst = sorted((Metadata_dist_df.index.tolist()))
# 個人の選択
Participant = st.selectbox('被験者を選択してください',Participant_lst)

# サンプリングデータ
Sampling_df = Metadata_df[['data_type','week_num','Participant ID']]
Sampling_chart = alt.Chart(Sampling_df[(Sampling_df['Participant ID']==Participant)&(Sampling_df['data_type'].isin(['metagenomics','metabolomics']))]).mark_point(
        filled=True, 
        size=200,
        opacity=0.7
    ).encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('data_type',
            sort=['metagenomics','metabolomics'],
            axis=alt.Axis(labelAngle=-45,labelFontSize=15, ticks=True, titleFontSize=18, title='サンプルの種類')
            ),
        color=alt.Color('data_type',legend=None, scale=alt.Scale(domain=['metagenomics','metabolomics'],range=['steelblue', 'darkorange']),),
        tooltip='week_num'
    ).properties(
        width=500, height=200
    )

Participant_info = Metadata_df[Metadata_df['Participant ID']==Participant][['diagnosis','sex','Age at diagnosis','consent_age','Occupation','Education Level']]
# Participant_info = Participant_info.astype({'Age at diagnosis': int, 'consent_age': int})
Participant_info = Participant_info.astype({'Age at diagnosis': str, 'consent_age': str})
Participant_info = Participant_info.iloc[0]
Participant_info.index = ['病名','性別','発症時の年齢','年齢','職業','学歴']
Participant_info.name = Participant

md_col1, md_col2 = st.columns(2)
with md_col1:
    st.markdown('### 患者の情報')
    st.table(Participant_info)
with md_col2:
    st.markdown('### サンプリング時期')

    st.write(Sampling_chart)
    
    text = '''

    初診日を 0 週として、メタゲノムデータおよびメタボロームデータがどのタイミングで採取されたのかを表示しています。

    '''

    st.markdown(text)

text = '''

性質の似ている被験者を調べたい場合は [PLS-ROG](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/cem.2883) の解析結果を参考にしましょう。
[こちら](https://keisukeota.shinyapps.io/IBDMDB/)を参照してください。

'''

st.markdown(text)

text = '''

また性別や年齢層、学歴、職業といった患者の属性ごとに IBD の特徴となる菌叢や代謝物を探索することができます。
サイドバーで被験者の選択条件を設定しましょう。選択条件を満たした被験者の度数分布が以下に表示されます。

'''

st.markdown(text)

@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def Metadata_dist_chart(select_column,axis_title):
    if not select_column == 'consent_age':
        Metadata_dist_chart = alt.Chart(Metadata_dist_df).mark_bar(opacity=0.5).encode(
            x=alt.X(select_column, 
                axis=alt.Axis(labelAngle=-45, labelFontSize=15, ticks=True, titleFontSize=18, title=axis_title)
                ),
            y=alt.Y('count()',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='人数'),
                stack=None
                ),
            color=alt.Color('diagnosis',legend=None),
            column=alt.Column('diagnosis:N',
                header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                )
            ).properties(
            width=1000/len(select_diag),
            height=300
            ).interactive()
    else:
        Metadata_dist_chart = alt.Chart(Metadata_dist_df).mark_bar(opacity=0.5).encode(
            x=alt.X(select_column, 
                bin=alt.Bin(maxbins=6),
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title=axis_title)
                ),
            y=alt.Y('count()',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='人数'),
                stack=None
                ),
            color=alt.Color('diagnosis',legend=None),
            column=alt.Column('diagnosis:N',
                header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                )
            ).properties(
            width=1000/len(select_diag),
            height=300
            ).interactive()
    return Metadata_dist_chart

with st.expander('被験者の度数分布'):

    Participant_info_lst = ['性別','年齢','職業','学歴']
    select_Participant_info = st.selectbox('項目を選んでください',Participant_info_lst, index=0,key='select_Participant_info_lst')

    if select_Participant_info == '性別':
        Metadata_dist_chart = Metadata_dist_chart('sex','性別')
    elif select_Participant_info == '年齢':
        Metadata_dist_chart = Metadata_dist_chart('consent_age','年齢')
    elif select_Participant_info == '学歴':
        Metadata_dist_chart = Metadata_dist_chart('Education Level','学歴')
    elif select_Participant_info == '職業':
         Metadata_dist_chart = Metadata_dist_chart('Occupation','職業')

    st.write(Metadata_dist_chart)

# メタゲノムデータの可視化

st.markdown('# 菌叢の細菌組成')

text = '''

被験者の時系列細菌組成データがプルダウンで選択された分類階級にしたがって表示されます。着目したい種をクリックすると他の種が灰色で表示されます。2 種目以降をクリックする場合は、shift キーを押しながらクリックしましょう。
``Suggestion`` にチェックを入れると nonIBD との比較解析を準備できます。

'''

st.markdown(text)

# 分類階級を選択
Level_lst = ['species','genus','family','order','class','phylum']
Level = st.selectbox('分類階級を選んでください',Level_lst, index=5,key='Stool_MGX_Microbes_select_Level')
Stool_MGX_Microbes_Analytics = st.checkbox('Suggestion', key='Stool_MGX_Microbes_Suggestion')
Stool_MGX_Microbes_Analytics_exec = False
if Stool_MGX_Microbes_Analytics == True:

    st.markdown('``Exec`` にチェックを入れると解析が実行されます。')
    Stool_MGX_Microbes_Analytics_exec = st.checkbox('Exec', key='Stool_MGX_Microbes_Analytics_exec')

    with st.expander('各患者のサンプル抽出数'):
        text = '''

        1 年間で糞便を提供した回数は被験者によって様々です。nonIBD との比較解析において、各被験者から何回ずつサンプルを抽出するのかをスライダーで選択してください。サンプル提供回数がスライダーで選択された回数に満たない被験者は解析から除外されます。\n
        ``各患者のサンプリング回数`` をクリックすると、左側に各被験者のサンプル提供回数が表示され、右側にスライダーで選択された回数以上サンプルを提供した被験者の人数が表示されます。
        
        '''
        
        st.markdown(text)

    min_ns = st.slider('各患者のサンプル抽出数を選んでください',min_value=1, max_value=25,step=1, key='Stool_MGX_Microbes_min_ns') # minimum number of samples
    
    # if Stool_MGX_Microbes_Analytics_exec == True:
    #     Stool_MGX_Microbes_Dysbiotic_CD = st.sidebar.checkbox('Dysbiotic CD サンプルを抽出する')
    #     if Stool_MGX_Microbes_Dysbiotic_CD == True:
    #         select_CD_sample = st.sidebar.slider('CD 患者の dysbiosis score を選択してください',0.00, 1.00, (0.00, 1.00))
    #     Stool_MGX_Microbes_Dysbiotic_UC = st.sidebar.checkbox('Dysbiotic UC サンプルを抽出する')
    #     if Stool_MGX_Microbes_Dysbiotic_UC == True:
    #         select_UC_sample = st.sidebar.slider('UC 患者の dysbiosis score を選択してください',0.00, 1.00, (0.00, 1.00))

@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def Stool_MGX_Microbes_chart(): 

    # 細菌分布の可視化
    Level_radio = alt.binding_radio(options=Level_lst)
    Level_select = alt.selection_single(fields=['Level'], bind=Level_radio, name="select", init={'Level': Level_lst[0]})

    Microbes_select = alt.selection_multi(empty='all', fields=['Microbes'])

    bars = alt.Chart().mark_bar().encode(
        x=alt.X('week_num:O',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('sum(Composition)',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='組成（％）')
            ),
        color=alt.condition(Microbes_select,'Microbes',alt.value('lightgray'),legend=None),
        tooltip=['Microbes','Composition']
    ).properties(width=500, height=500).add_selection(Microbes_select)

    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection(type='single', nearest=True, on='mouseover',fields=['week_num'], empty='none')

    # The basic line
    line = alt.Chart().mark_line().encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('Composition:Q',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='組成（％）')
            ),
        color = 'Microbes'
    )

    # Transparent selectors across the chart. This is what tells us
    # the x-value of the cursor
    selectors = alt.Chart().mark_point().encode(
        x=alt.X(
            'week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        opacity=alt.value(0),
    ).add_selection(nearest)

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
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
    ).transform_filter(nearest)

    # Put the five layers into a chart and bind the data
    lines = alt.layer(
        line, selectors, points, rules, text
    ).properties(
        width=500, height=500
    ).transform_filter(
        Microbes_select
    )

    Stool_MGX_Microbes_chart = alt.hconcat(bars, lines, data = Stool_MGX_Microbes_df[(Stool_MGX_Microbes_df['Participant ID'] == Participant)&(Stool_MGX_Microbes_df['Level']==Level)])

    # Stool_MGX_Microbes_Analytics = st.checkbox('Microbes Suggestion by Analytics')
    if Stool_MGX_Microbes_Analytics == True :

        # メタゲノムデータ提供者の処理
        member_df = Metadata_df[['Participant ID','data_type','week_num','diagnosis']]
        Stool_MGX_member_df = member_df[member_df['data_type']=='metagenomics'][['Participant ID','week_num','diagnosis']]
        Stool_MGX_member_df['External ID'] = Stool_MGX_member_df.index.tolist()

        # Defaultdict を用いて、各患者が 1 年間で何回サンプルを提供してくれたのかをカウント
        Stool_MGX_member_counter = defaultdict(int)
        for key in Stool_MGX_member_df['Participant ID'].values.tolist():
            Stool_MGX_member_counter[key] += 1

        # 最低回数以上サンプルを提供してくれた患者をリストアップ
        # min_ns = st.slider('患者の最低サンプル提供回数を選んでください（菌叢解析）',min_value=1, max_value=25,step=1) # minimum number of samples
        Stool_MGX_input_members = [k for k, i in Stool_MGX_member_counter.items() if i >= min_ns]
        tmp = []
        for P in Participant_lst:
            if P in Stool_MGX_input_members:
                tmp.append(P)
        Stool_MGX_input_members = tmp

        # 各患者から何回目に採取されたサンプルなのかをカウント
        Participant_memo = Stool_MGX_member_df['Participant ID'].values.tolist()[0]
        cnt = 0
        cnt_lst = []
        for row in Stool_MGX_member_df.values.tolist():
            if Participant_memo == row[0]:
                cnt += 1
                cnt_lst.append(cnt)
            else:
                Participant_memo = row[0]
                cnt = 1
                cnt_lst.append(cnt)
        Stool_MGX_member_df['Sampling_counter'] = cnt_lst

        # 各被験者の採取回数を表示

        histgram = alt.Chart(Stool_MGX_member_df).mark_bar(opacity=0.5).encode(
            y=alt.Y("Participant ID:O", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='被験者'),
                sort='-x'
                ),
            x=alt.X('count(Participant ID)',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='サンプリング回数'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(width=500,height=1600).interactive()

        bar = alt.Chart(Stool_MGX_member_df[Stool_MGX_member_df['Sampling_counter']==min_ns]).mark_bar(opacity=0.5).encode(
            y=alt.Y("count(Participant ID)", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='人数'),
                sort='-x'
                ),
            x=alt.X('diagnosis',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='診断'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(
            width=150,height=1600
        ).interactive()

        Stool_MGX_Microbes_sample_dist_chart = histgram|bar

        # Stool_MGX_Microbes_Analytics_exec = st.checkbox('Exec（菌叢解析）')
        if Stool_MGX_Microbes_Analytics_exec == True:

            # リストアップされた患者の菌叢データを集める
            Stool_MGX_Microbes_Level_df = Stool_MGX_Microbes_df[Stool_MGX_Microbes_df['Level']==Level]
            Stool_MGX_Microbes_Level_messy_df = pd.pivot_table(Stool_MGX_Microbes_Level_df, index=['Participant ID', 'week_num'], columns='Microbes', values='Composition')
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df.T[Stool_MGX_input_members].T
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df.reset_index()
            Stool_MGX_Microbes_Level_messy_df = pd.merge(Stool_MGX_Microbes_Level_messy_df,Stool_MGX_member_df,on=['Participant ID','week_num'],how='left')
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df[Stool_MGX_Microbes_Level_messy_df['Sampling_counter']<=min_ns]
            
            # Dysbiosis score の計算
            # まずは Bray-Curtis 距離の計算
            Stool_MGX_Microbes_Level_messy_config_df = Stool_MGX_Microbes_Level_messy_df[['Participant ID','week_num','diagnosis','External ID']].copy()
            Stool_MGX_Microbes_Level_BC_input_df = Stool_MGX_Microbes_Level_messy_df.drop(['Participant ID','week_num','diagnosis','External ID','Sampling_counter'],axis=1)
            from skbio.diversity import beta_diversity
            Stool_MGX_Microbes_Level_BC_output_df = beta_diversity("braycurtis", Stool_MGX_Microbes_Level_BC_input_df.values, Stool_MGX_Microbes_Level_messy_config_df['External ID'].values.tolist())
            Stool_MGX_Microbes_Level_BC_output_df = pd.DataFrame(Stool_MGX_Microbes_Level_BC_output_df[0:len(Stool_MGX_Microbes_Level_messy_config_df['External ID'].values)])
            
            # ここから Dysbiosis score の計算
            DS_lst = []
            i = 0
            for row1 in Stool_MGX_Microbes_Level_messy_config_df.values.tolist():
                cnt = 0
                # 自分以外の nonIBD サンプルをメモするリスト row_lst
                row_lst = []
                for row2 in Stool_MGX_Microbes_Level_messy_config_df.values.tolist():
                    if not row1[0] == row2[1]:
                        if row2[2] == 'nonIBD':
                            row_lst.append(cnt)
                    cnt += 1
                # row_lst にメモされている相手との BC 距離を Dysbiosis score とする
                DS_lst.append(np.median(Stool_MGX_Microbes_Level_BC_output_df.values[i,row_lst]))
                i += 1

            Stool_MGX_Microbes_Level_messy_config_df['Dysbiosis_Score'] = DS_lst

            # Dysbiosis score 描画

            # Create a selection that chooses the nearest point & selects based on x-value
            nearest = alt.selection(type='single', nearest=True, on='mouseover',fields=['week_num'], empty='none')

            # The basic line
            line = alt.Chart().mark_line().encode(
                x=alt.X('week_num',
                    scale=alt.Scale(domain=[-5,60]),
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週')
                    ),
                y=alt.Y('Dysbiosis_Score:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='Dysbiosis score'),
                    scale=alt.Scale(domain=[0,1])
                    ),
            )

            # Transparent selectors across the chart. This is what tells us
            # the x-value of the cursor
            selectors = alt.Chart().mark_point().encode(
                x=alt.X('week_num',
                    scale=alt.Scale(domain=[-5,60]),
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週')
                    ),
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
                text=alt.condition(nearest, 'Dysbiosis_Score:Q', alt.value(' '))
            )

            # Draw a rule at the location of the selection
            rules = alt.Chart().mark_rule(color='gray').encode(
                x=alt.X('week_num',
                    scale=alt.Scale(domain=[-5,60]),
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週')
                    ),
            ).transform_filter(
                nearest
            )

            # Put the five layers into a chart and bind the data
            Stool_MGX_Microbes_DS_Participant_chart = alt.layer(
                line, selectors, points, rules, text, data = Stool_MGX_Microbes_Level_messy_config_df[Stool_MGX_Microbes_Level_messy_config_df['Participant ID']==Participant]
            ).properties(
                width=500, height=500
            )

            jitter = alt.Chart(Stool_MGX_Microbes_Level_messy_config_df).mark_circle(size=20).encode(
                x=alt.X('jitter:Q',
                    title=None,
                    axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                    scale=alt.Scale(),
                    ),
                y=alt.Y('Dysbiosis_Score:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='Dysbiosis score')
                    ),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                    ),
                tooltip=['Participant ID','week_num']
            ).transform_calculate(
                # Generate Gaussian jitter with a Box-Muller transform
                jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
            ).properties(
                width=40,
                height=500
            )

            boxplot = alt.Chart(Stool_MGX_Microbes_Level_messy_config_df).mark_boxplot(
                opacity=0.5,
                size=30,
                ticks=alt.MarkConfig(width=30), 
                median=alt.MarkConfig(color='black',size=30),
                outliers=alt.MarkConfig(opacity=0)
            ).encode(
                y=alt.Y('Dysbiosis_Score:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title=None)
                    ),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                    )
            ).properties(width=40,height=500)

            Stool_MGX_Microbes_DS_chart = Stool_MGX_Microbes_DS_Participant_chart | jitter | boxplot

            # ここから ANCOM 前処理
            # # Dysbiosis score で解析サンプルをトリミングすることができる
            # Stool_MGX_Microbes_Level_messy_df = pd.merge(Stool_MGX_Microbes_Level_messy_df,Stool_MGX_Microbes_Level_messy_config_df[['Participant ID','week_num','Dysbiosis_Score']],on=['Participant ID','week_num'],how='left')
            # if Stool_MGX_Microbes_Dysbiotic_CD == True:
            #     Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df[((Stool_MGX_Microbes_Level_messy_df['Dysbiosis_Score'] > select_CD_sample[0])&(Stool_MGX_Microbes_Level_messy_df['Dysbiosis_Score'] < select_CD_sample[1])&(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'CD'))|(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'UC')|(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'nonIBD')]
            # if Stool_MGX_Microbes_Dysbiotic_UC == True:
            #     Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df[((Stool_MGX_Microbes_Level_messy_df['Dysbiosis_Score'] > select_UC_sample[0])&(Stool_MGX_Microbes_Level_messy_df['Dysbiosis_Score'] < select_UC_sample[1])&(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'UC'))|(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'CD')|(Stool_MGX_Microbes_Level_messy_df['diagnosis'] == 'nonIBD')]
            # Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df.drop('Dysbiosis_Score',axis=1)
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df.set_index(['External ID','diagnosis'])

            ## 可視化する細菌の種類をリストにまとめてメモしておく
            Stool_MGX_Microbes_Level_sum_series = Stool_MGX_Microbes_Level_messy_df[Stool_MGX_Microbes_Level_messy_df['Participant ID']==Participant].drop(['Participant ID','week_num','Sampling_counter'],axis=1).sum(axis=0)
            # st.write(Stool_MGX_Microbes_Level_sum_series)
            Stool_MGX_Microbes_lst = Stool_MGX_Microbes_Level_sum_series[Stool_MGX_Microbes_Level_sum_series > 0].index.tolist()
            
            ## Psudo-count
            Stool_MGX_Microbes_ANCOM_input_df = Stool_MGX_Microbes_Level_messy_df.drop(['Participant ID','week_num','Sampling_counter'],axis=1)
            Stool_MGX_Microbes_ANCOM_input_df = pd.DataFrame(
                data = multiplicative_replacement(Stool_MGX_Microbes_ANCOM_input_df.values),
                index = Stool_MGX_Microbes_ANCOM_input_df.index, 
                columns = Stool_MGX_Microbes_ANCOM_input_df.columns
                )

            # ANCOM 実施
            Stool_MGX_Microbes_ANCOM_input_df = Stool_MGX_Microbes_ANCOM_input_df.reset_index()
            Stool_MGX_Microbes_ANCOM_input_df = Stool_MGX_Microbes_ANCOM_input_df.set_index('External ID')
            Stool_MGX_Microbes_ANCOM_input_series = Stool_MGX_Microbes_ANCOM_input_df['diagnosis']
            Stool_MGX_Microbes_ANCOM_results = ancom(Stool_MGX_Microbes_ANCOM_input_df.drop(['diagnosis'],axis=1), Stool_MGX_Microbes_ANCOM_input_series)[0]
            Stool_MGX_Microbes_ANCOM_results = Stool_MGX_Microbes_ANCOM_results[Stool_MGX_Microbes_ANCOM_results['Reject null hypothesis']]['W']
            
            # ANCOM 結果の描画
            Stool_MGX_Microbes_ANCOM_results = pd.DataFrame(Stool_MGX_Microbes_ANCOM_results)
            Stool_MGX_Microbes_ANCOM_results = Stool_MGX_Microbes_ANCOM_results.sort_values('W',ascending=False).head(20)
            Stool_MGX_Microbes_ANCOM_results['Microbes'] = Stool_MGX_Microbes_ANCOM_results.index.tolist()

            # ANCOM で CD、nonIBD、UC 間で統計的有意差のある細菌の生データを可視化
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df[list(set(Stool_MGX_Microbes_ANCOM_results['Microbes'].values.tolist()+Stool_MGX_Microbes_lst))]
            Stool_MGX_Microbes_Level_messy_df = Stool_MGX_Microbes_Level_messy_df.reset_index()
            Stool_MGX_Microbes_Level_messy_df = pd.merge(Stool_MGX_Microbes_Level_messy_df,Stool_MGX_member_df.drop(['diagnosis','Sampling_counter','External ID'],axis=1),on=['External ID'],how='left')
            Stool_MGX_Microbes_Level_tidy_df = pd.melt(Stool_MGX_Microbes_Level_messy_df,id_vars=['External ID','diagnosis','Participant ID','week_num'],var_name="Microbes",value_name="Composition" )

            selector = alt.selection_single(empty='all', fields=['Microbes'])

            Stool_MGX_Microbes_ANCOM_results_chart = alt.Chart(Stool_MGX_Microbes_ANCOM_results).mark_bar().encode(
                x=alt.X('Microbes:N',
                    sort='-y',
                    axis=alt.Axis(labelAngle=-45,labelFontSize=6, ticks=True, titleFontSize=18, title='細菌')
                    ),
                y=alt.Y('W',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='W（変数重要度）')
                    ),
                color=alt.condition(Microbes_select,'Microbes:N',alt.value('lightgray'),legend=None),
                tooltip='Microbes'
            ).properties(width=500, height=500).add_selection(Microbes_select)

            jitter = alt.Chart(Stool_MGX_Microbes_Level_tidy_df).mark_circle(size=20).encode(
                x=alt.X('jitter:Q',
                    title=None,
                    axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                    scale=alt.Scale(),
                ),
                y=alt.Y('Composition:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='組成（％）')
                    ),
                color=alt.Color('Microbes',legend=None),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                    ),
                tooltip=['Participant ID','week_num']
            ).transform_calculate(
                # Generate Gaussian jitter with a Box-Muller transform
                jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
            ).add_selection(
                Microbes_select
            ).transform_filter(
                Microbes_select
            ).properties(
                width=59,height=500
            )

            boxplot = alt.Chart(Stool_MGX_Microbes_Level_tidy_df).mark_boxplot(
                opacity=0.5,
                size=30,
                ticks=alt.MarkConfig(width=30), 
                median=alt.MarkConfig(color='black',size=30),
                outliers=alt.MarkConfig(opacity=0)
            ).encode(
                y=alt.Y('Composition:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title=None)
                    ),
                color=alt.Color('Microbes',legend=None),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, labelColor='gray',titleFontSize=18, title='病名', titleColor='gray')
                    )
            ).add_selection(
                Microbes_select
            ).transform_filter(
                Microbes_select
            ).properties(
                width=59,height=500
            )

            Stool_MGX_Microbes_ANCOM_results_chart = Stool_MGX_Microbes_ANCOM_results_chart | jitter | boxplot
            Stool_MGX_Microbes_chart = Stool_MGX_Microbes_chart & Stool_MGX_Microbes_ANCOM_results_chart

            return Stool_MGX_Microbes_chart, Stool_MGX_Microbes_DS_chart, Stool_MGX_Microbes_sample_dist_chart
        else: 
            return Stool_MGX_Microbes_chart, Stool_MGX_Microbes_sample_dist_chart
    else:
        return Stool_MGX_Microbes_chart


if Stool_MGX_Microbes_Analytics == True :
    if Stool_MGX_Microbes_Analytics_exec == True:
        Stool_MGX_Microbes_chart, Stool_MGX_Microbes_DS_chart, Stool_MGX_Microbes_sample_dist_chart = Stool_MGX_Microbes_chart()
        with st.expander('各患者のサンプリング回数'):
            st.write(Stool_MGX_Microbes_sample_dist_chart)

        st.write(Stool_MGX_Microbes_chart)

        text = '''

        解析結果が下に追加されます。下左図の各プロットは細菌の種類に対応しています。カーソルを合わせると細菌名を確認できます。
        縦軸は **[ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/)** で計算される変数重要度 **W** を表しています。
        すなわち **W** の大きな細菌が CD と UC と nonIBD を区別する際に重要であるということです。重要な代謝経路が見つからなかった場合は何も表示されません。\n
        下右図では、細菌ごとに組成の分布が箱ひげ図で表示されます。ANCOM の W では群間の大小関係はわからないので、CD、UC、nonIBD の大小比較に用いるとよいでしょう。\n
        着目したい細菌をクリックすると他の種が灰色で表示されます。上の時系列データとも連動しています。色々クリックしてみましょう。

        '''
        st.markdown(text)

        st.write(Stool_MGX_Microbes_DS_chart)

        text = '''

        また菌叢組成の乱れ度を表す **Dysbiosis score** も計算されます。Dysbiosis score は nonIBD の各被験者に対する Bray-Curtis 距離の中央値として定義されます。
        すなわち、nonIBD の集団と菌叢の組成が異なるほど値は 1 近づきます。反対に nonIBD の集団と菌叢の組成が似ているほど値は 0 に近づきます。

        '''
        st.markdown(text)
        
    else:
        Stool_MGX_Microbes_chart, Stool_MGX_Microbes_sample_dist_chart = Stool_MGX_Microbes_chart()
        with st.expander('各患者のサンプリング回数'):
            st.write(Stool_MGX_Microbes_sample_dist_chart)
        st.write(Stool_MGX_Microbes_chart)
else:
    Stool_MGX_Microbes_chart = Stool_MGX_Microbes_chart()
    st.write(Stool_MGX_Microbes_chart)

# 代謝酵素分布の可視化

st.markdown('# 菌叢の機能組成')

text = '''

菌叢の機能組成データが表示されます。菌叢の機能とは、菌叢が発現する酵素によって行われる**代謝**のことです。ショットガンメタゲノムデータには細菌がもつ酵素遺伝子の read が含まれています。
ここでは **[HUMAnN3.0](https://qiita.com/keisuke-ota/items/4e34ceffb33f2e983f09)** によって推定される酵素を発現するポテンシャルをみていることに注意しましょう。\n
なお代謝経路の相対存在量は、** Copies Per Million :CPM ** という 100 万 bp あたりの read 数に正規化されています。``Suggestion`` にチェックを入れると nonIBD との比較解析を準備できます。

'''

st.markdown(text)

Stool_MGX_PWY_Analytics = st.checkbox('Suggestion', key='Stool_MGX_PWY_Suggestion')
Stool_MGX_PWY_Analytics_exec = False
if Stool_MGX_PWY_Analytics == True:

    st.markdown('``Exec`` にチェックを入れると解析が実行されます。')
    Stool_MGX_PWY_Analytics_exec = st.checkbox('Exec', key='Stool_MGX_PWY_Analytics_exec')

    with st.expander('各患者のサンプル抽出数'):
        text = '''

        1 年間で糞便を提供した回数は被験者によって様々です。nonIBD との比較解析において、各被験者から何回ずつサンプルを抽出するのかをスライダーで選択してください。サンプル提供回数がスライダーで選択された回数に満たない被験者は解析から除外されます。\n
        ``各患者のサンプリング回数`` をクリックすると、左側に各被験者のサンプル提供回数が表示され、右側にスライダーで選択された回数以上サンプルを提供した被験者の人数が表示されます。 
        
        '''
        
        st.markdown(text)

    min_ns = st.slider('各患者のサンプル抽出数を選んでください',min_value=1, max_value=25,step=1, key='Stool_MGX_PWY_min_ns') # minimum number of samples
    
@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def Stool_MGX_PWY_chart(): 

    PWY_select = alt.selection_multi(empty='all', fields=['PWY'])

    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection(type='single', nearest=True, on='mouseover',
                            fields=['week_num'], empty='none')

    # The basic line
    line = alt.Chart().mark_line().encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('RPKs(CPM):Q',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='CPM')
            ),
        color='PWY'
    )

    # Transparent selectors across the chart. This is what tells us
    # the x-value of the cursor
    selectors = alt.Chart().mark_point().encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        opacity=alt.value(0),
    ).add_selection(nearest)

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
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
    ).transform_filter(nearest)

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


    bars = alt.Chart().mark_bar().encode(
        x=alt.X('week_num:O',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('sum(RPKs(CPM))',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='CPM')
            ),
        color=alt.condition(PWY_select,'PWY',alt.value('lightgray'),legend=None),
        tooltip=['PWY','RPKs(CPM)','week_num']
    ).properties(
        width=500, height=500
    ).add_selection(
        PWY_select
    )

    Stool_MGX_PWY_chart = alt.hconcat(bars, lines, data = Stool_MGX_PWY_df[Stool_MGX_PWY_df['Participant ID']==Participant])

    # Stool_MGX_PWY_Analytics = st.checkbox('Pathway Suggestion by Analytics')
    if Stool_MGX_PWY_Analytics == True :

        # メタゲノムデータ提供者の処理
        member_df = Metadata_df[['Participant ID','data_type','week_num','diagnosis']]
        Stool_MGX_member_df = member_df[member_df['data_type']=='metagenomics'][['Participant ID','week_num','diagnosis']]
        Stool_MGX_member_df['External ID'] = Stool_MGX_member_df.index.tolist()
        # st.table(Stool_MGX_member_df.head(20))

        # Defaultdict を用いて、各患者が 1 年間で何回サンプルを提供してくれたのかをカウント
        Stool_MGX_member_counter = defaultdict(int)
        for key in Stool_MGX_member_df['Participant ID'].values.tolist():
            Stool_MGX_member_counter[key] += 1

        # 最低回数以上サンプルを提供してくれた患者をリストアップ
        # min_ns = st.slider('患者の最低サンプル提供回数を選んでください（代謝経路解析）',min_value=1, max_value=25,step=1) # minimum number of samples
        Stool_MGX_input_members = [k for k, i in Stool_MGX_member_counter.items() if i >= min_ns]
        tmp = []
        for P in Participant_lst:
            if P in Stool_MGX_input_members:
                tmp.append(P)
        Stool_MGX_input_members = tmp

        # 各患者から何回目に採取されたサンプルなのかをカウント
        Participant_memo = Stool_MGX_member_df['Participant ID'].values.tolist()[0]
        cnt = 0
        cnt_lst = []
        for row in Stool_MGX_member_df.values.tolist():
            if Participant_memo == row[0]:
                cnt += 1
                cnt_lst.append(cnt)
            else:
                Participant_memo = row[0]
                cnt = 1
                cnt_lst.append(cnt)
        Stool_MGX_member_df['Sampling_counter'] = cnt_lst
        # st.table(Stool_MGX_member_df.head(20))

        # 各被験者の採取回数を表示

        histgram = alt.Chart(Stool_MGX_member_df).mark_bar(opacity=0.5).encode(
            y=alt.Y("Participant ID:O", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='被験者'),
                sort='-x'
                ),
            x=alt.X('count(Participant ID)',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='サンプリング回数'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(
            width=500, height=1600
        ).interactive()

        bar = alt.Chart(Stool_MGX_member_df[Stool_MGX_member_df['Sampling_counter']==min_ns]).mark_bar(opacity=0.5).encode(
            y=alt.Y("count(Participant ID)", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='人数'),
                sort='-x'
                ),
            x=alt.X('diagnosis',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='診断'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(
            width=150,height=1600
        ).interactive()

        Stool_MGX_PWY_sample_dist_chart = histgram|bar

        # Stool_MGX_PWY_Analytics_exec = st.checkbox('Exec（代謝経路解析）')
        if Stool_MGX_PWY_Analytics_exec == True:

            # リストアップされた患者の代謝経路データを集める
            Stool_MGX_PWY_messy_df = pd.pivot_table(Stool_MGX_PWY_df, index=['Participant ID', 'week_num'], columns='PWY', values='RPKs(CPM)')
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.T[Stool_MGX_input_members].T
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.reset_index()
            Stool_MGX_PWY_messy_df = pd.merge(Stool_MGX_PWY_messy_df,Stool_MGX_member_df,on=['Participant ID','week_num'],how='left')
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df[Stool_MGX_PWY_messy_df['Sampling_counter']<=min_ns]
            
            # st.table(Stool_MGX_PWY_messy_df.head(20))
            ## 可視化する代謝経路の種類をリストにまとめてメモしておく
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.set_index(['External ID','diagnosis'])
            Stool_MGX_PWY_sum_series = Stool_MGX_PWY_messy_df[Stool_MGX_PWY_messy_df['Participant ID']==Participant].drop(['Participant ID','week_num','Sampling_counter'],axis=1).sum(axis=0)
            Stool_MGX_PWY_lst = Stool_MGX_PWY_sum_series[Stool_MGX_PWY_sum_series > 0].index.tolist()
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.reset_index()
            ##

            # # Dysbiosis score で解析サンプルをトリミングすることができる
            # if Stool_MGX_Microbes_Analytics_exec == True:
            #     Stool_MGX_PWY_messy_df = pd.merge(Stool_MGX_PWY_messy_df,Stool_MGX_Microbes_Level_messy_config_df[['Participant ID','week_num','Dysbiosis_Score']],on=['Participant ID','week_num'],how='left')
            #     if Stool_MGX_Microbes_Dysbiotic_CD == True:
            #         Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df[((Stool_MGX_PWY_messy_df['Dysbiosis_Score'] > select_CD_sample[0])&(Stool_MGX_PWY_messy_df['Dysbiosis_Score'] < select_CD_sample[1])&(Stool_MGX_PWY_messy_df['diagnosis'] == 'CD'))|(Stool_MGX_PWY_messy_df['diagnosis'] == 'UC')|(Stool_MGX_PWY_messy_df['diagnosis'] == 'nonIBD')]
            #     if Stool_MGX_Microbes_Dysbiotic_UC == True:
            #         Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df[((Stool_MGX_PWY_messy_df['Dysbiosis_Score'] > select_UC_sample[0])&(Stool_MGX_PWY_messy_df['Dysbiosis_Score'] < select_UC_sample[1])&(Stool_MGX_PWY_messy_df['diagnosis'] == 'UC'))|(Stool_MGX_PWY_messy_df['diagnosis'] == 'CD')|(Stool_MGX_PWY_messy_df['diagnosis'] == 'nonIBD')]
            #     Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df[Stool_MGX_PWY_messy_df['Dysbiosis_Score'] >= 0] # 同じタイミングで MGX が計測されていない場合は Dysbiosis score を計算できないので削除
            #     # st.table(Stool_MGX_PWY_messy_df.head())
            #     Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.drop('Dysbiosis_Score',axis=1)


            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.set_index(['External ID','diagnosis'])
            Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.drop(['Participant ID','week_num','Sampling_counter'],axis=1)            # ## 全ての代謝経路の存在量が 0 のサンプルを除く
            # st.table(Stool_MGX_PWY_messy_df.T[[('HSM6XRR9', 'UC'),('HSMA33IO', 'UC'),('MSM6J2HB', 'UC'),('PSM6XBTD', 'UC')]].head())
            Stool_MGX_PWY_sum_series = Stool_MGX_PWY_messy_df.sum(axis=1)
            Stool_MGX_PWY_zero_lst = Stool_MGX_PWY_sum_series[Stool_MGX_PWY_sum_series == 0].index.tolist()
            if len(Stool_MGX_PWY_zero_lst) > 0:
                Stool_MGX_PWY_messy_df = Stool_MGX_PWY_messy_df.drop(Stool_MGX_PWY_zero_lst,axis=0)
            # st.table(Stool_MGX_PWY_messy_df.head(20))
            ## Psudo-count
            Stool_MGX_PWY_ANCOM_input_df = pd.DataFrame(
                data = multiplicative_replacement(Stool_MGX_PWY_messy_df.values),
                index = Stool_MGX_PWY_messy_df.index, 
                columns = Stool_MGX_PWY_messy_df.columns
                )

            # ANCOM 実施
            Stool_MGX_PWY_ANCOM_input_df = Stool_MGX_PWY_ANCOM_input_df.reset_index()
            Stool_MGX_PWY_ANCOM_input_df = Stool_MGX_PWY_ANCOM_input_df.set_index('External ID')
            Stool_MGX_PWY_ANCOM_input_series = Stool_MGX_PWY_ANCOM_input_df['diagnosis']
            Stool_MGX_PWY_ANCOM_results = ancom(Stool_MGX_PWY_ANCOM_input_df.drop(['diagnosis'],axis=1), Stool_MGX_PWY_ANCOM_input_series)[0]
            Stool_MGX_PWY_ANCOM_results = Stool_MGX_PWY_ANCOM_results[Stool_MGX_PWY_ANCOM_results['Reject null hypothesis']]['W']

            # ANCOM 結果の描画
            Stool_MGX_PWY_ANCOM_results = pd.DataFrame(Stool_MGX_PWY_ANCOM_results)
            Stool_MGX_PWY_ANCOM_results = Stool_MGX_PWY_ANCOM_results.sort_values('W',ascending=False).head(20)
            Stool_MGX_PWY_ANCOM_results['PWY'] = Stool_MGX_PWY_ANCOM_results.index.tolist()

            # ANCOM で CD、nonIBD、UC 間で統計的有意差のある細菌の生データを可視化
            Stool_MGX_PWY_Significant_messy_df = Stool_MGX_PWY_messy_df[list(set(Stool_MGX_PWY_ANCOM_results['PWY'].values.tolist()+Stool_MGX_PWY_lst))]
            Stool_MGX_PWY_Significant_messy_df = Stool_MGX_PWY_Significant_messy_df.reset_index()
            Stool_MGX_PWY_Significant_messy_df = pd.merge(Stool_MGX_PWY_Significant_messy_df,Stool_MGX_member_df.drop(['diagnosis','Sampling_counter','External ID'],axis=1),on=['External ID'],how='left')
            Stool_MGX_PWY_Significant_tidy_df = pd.melt(Stool_MGX_PWY_Significant_messy_df,id_vars=['External ID','diagnosis','Participant ID','week_num'],var_name="PWY",value_name="RPKs(CPM)" )

            Stool_MGX_PWY_ANCOM_results_chart = alt.Chart(Stool_MGX_PWY_ANCOM_results).mark_bar().encode(
                x=alt.X('PWY:N',
                    sort='-y',
                    axis=alt.Axis(labelAngle=-45,labelFontSize=6, ticks=True, titleFontSize=18, title='代謝経路')
                    ),
                y=alt.Y('W',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='W')
                    ),
                color=alt.condition(PWY_select,'PWY:N',alt.value('lightgray'),legend=None),
                tooltip='PWY'
            ).properties(
                width=500, height=500
            ).add_selection(
                PWY_select
            )

            jitter = alt.Chart(Stool_MGX_PWY_Significant_tidy_df).mark_circle(size=20).encode(
                x=alt.X(
                    'jitter:Q',
                    title=None,
                    axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                    scale=alt.Scale(),
                ),
                y=alt.Y('RPKs(CPM):Q'),
                color=alt.Color('PWY',legend=None),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, titleFontSize=18,title='病名')
                    ),
                tooltip=['Participant ID','week_num']
                ).transform_calculate(
                    # Generate Gaussian jitter with a Box-Muller transform
                    jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
                ).add_selection(
                    PWY_select
                ).transform_filter(
                    PWY_select
                ).properties(
                    width=59,height=500
                )

            boxplot = alt.Chart(Stool_MGX_PWY_Significant_tidy_df).mark_boxplot(
                opacity=0.5,
                size=30,
                ticks=alt.MarkConfig(width=30), 
                median=alt.MarkConfig(color='black',size=30),
                outliers=alt.MarkConfig(opacity=0)
            ).encode(
                y='RPKs(CPM):Q',
                color=alt.Color('PWY',legend=None),
                column=alt.Column('diagnosis:N', 
                header=alt.Header(labelFontSize=15, titleFontSize=18,title='病名'))
            ).add_selection(
                PWY_select
            ).transform_filter(
                PWY_select
            ).properties(
                width=59,height=500
            )

            Stool_MGX_PWY_ANCOM_results_chart = Stool_MGX_PWY_ANCOM_results_chart | jitter | boxplot
            Stool_MGX_PWY_chart = Stool_MGX_PWY_chart & Stool_MGX_PWY_ANCOM_results_chart

        return Stool_MGX_PWY_chart, Stool_MGX_PWY_sample_dist_chart
    else:
        return Stool_MGX_PWY_chart

if Stool_MGX_PWY_Analytics == True :
    Stool_MGX_PWY_chart, Stool_MGX_PWY_sample_dist_chart = Stool_MGX_PWY_chart()
    with st.expander('各患者のサンプリング回数'):
        st.write(Stool_MGX_PWY_sample_dist_chart)
    
    if Stool_MGX_PWY_Analytics_exec == True: 

        text = '''

        解析結果が下に追加されます。下左図の各プロットは代謝経路に対応しています。カーソルを合わせると代謝経路名を確認できます。
        縦軸は **[ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/)** で計算される変数重要度 **W** を表しています。
        すなわち **W** の大きな代謝経路が CD と UC と nonIBD を区別する際に重要であるということです。重要な代謝経路が見つからなかった場合は何も表示されません。\n
        下右図では、代謝経路ごとに組成の分布が箱ひげ図で表示されます。ANCOM の W では群間の大小関係はわからないので、CD、UC、nonIBD の大小比較に用いるとよいでしょう。\n
        着目したい代謝経路をクリックすると他の代謝経路が灰色で表示されます。上の時系列データとも連動しています。色々クリックしてみましょう。

        '''

    else:

        text = ''

    st.write(Stool_MGX_PWY_chart)
    st.markdown(text)

else:
    Stool_MGX_PWY_chart = Stool_MGX_PWY_chart()
    st.write(Stool_MGX_PWY_chart)

# メタボローム分布可視化

st.markdown('# メタボロームデータ')

text = '''

被験者の時系列メタボロームデータが表示されます。各代謝物は LC-MS によって定量されています。``Suggestion`` にチェックを入れると nonIBD との比較解析を準備できます。

'''

st.markdown(text)
Stool_MBX_Analytics = st.checkbox('Suggestion', key='Stool_MBX_Suggestion')
Stool_MBX_Analytics_exec = False
if Stool_MBX_Analytics == True:

    st.markdown('``Exec`` にチェックを入れると解析が実行されます。')
    Stool_MBX_Analytics_exec = st.checkbox('Exec', key='Stool_MBX_Analytics_exec')

    with st.expander('各患者のサンプル抽出数'):
        text = '''

        1 年間で糞便を提供した回数は被験者によって様々です。nonIBD との比較解析において、各被験者から何回ずつサンプルを抽出するのかをスライダーで選択してください。サンプル提供回数がスライダーで選択された回数に満たない被験者は解析から除外されます。\n
        ``各患者のサンプリング回数`` をクリックすると、左側に各被験者のサンプル提供回数が表示され、右側にスライダーで選択された回数以上サンプルを提供した被験者の人数が表示されます。\n
        
        '''
        
        st.markdown(text)

    min_ns = st.slider('各患者のサンプル抽出数を選んでください',min_value=1, max_value=7,step=1, key='Stool_MBX_min_ns') # minimum number of samples
    
@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def Stool_MBX_chart(): 

    metabolites_select = alt.selection_multi(empty='all', fields=['Metabolites'])

    bars = alt.Chart().mark_bar().encode(
        x=alt.X('week_num:O',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('sum(Abundance)',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='存在量')
            ),
        color=alt.condition(metabolites_select,'Metabolites',alt.value('lightgray'),legend=None),
        tooltip=['Metabolites','Abundance']
    ).properties(
        width=500, height=500
    ).add_selection(
        metabolites_select
    )

    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection(type='single', nearest=True, on='mouseover',fields=['week_num'], empty='none')

    # The basic line
    line = alt.Chart().mark_line().encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
        y=alt.Y('Abundance:Q',
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='存在量')
            ),
        color = 'Metabolites'
    )

    # Transparent selectors across the chart. This is what tells us
    # the x-value of the cursor
    selectors = alt.Chart().mark_point().encode(
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
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
        x=alt.X('week_num',
            scale=alt.Scale(domain=[-5,60]),
            axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='週数')
            ),
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

    Stool_MBX_chart = alt.hconcat(bars, lines, data = Stool_MBX_df[Stool_MBX_df['Participant ID']==Participant])

    # Stool_MBX_Analytics = st.checkbox('Metabolites Suggestion by Analytics')
    if Stool_MBX_Analytics == True :

        # メタゲノムデータ提供者の処理
        member_df = Metadata_df[['Participant ID','data_type','week_num','diagnosis']]
        Stool_MBX_member_df = member_df[member_df['data_type']=='metabolomics'][['Participant ID','week_num','diagnosis']]
        Stool_MBX_member_df['External ID'] = Stool_MBX_member_df.index.tolist()
        # st.table(Stool_MBX_member_df.head(20))

        # Defaultdict を用いて、各患者が 1 年間で何回サンプルを提供してくれたのかをカウント
        Stool_MBX_member_counter = defaultdict(int)
        for key in Stool_MBX_member_df['Participant ID'].values.tolist():
            Stool_MBX_member_counter[key] += 1

        # 最低回数以上サンプルを提供してくれた患者をリストアップ
        # min_ns = st.slider('患者の最低サンプル提供回数を選んでください（メタボローム解析）',min_value=1, max_value=7,step=1) # minimum number of samples
        Stool_MBX_input_members = [k for k, i in Stool_MBX_member_counter.items() if i >= min_ns]
        tmp = []
        for P in Participant_lst:
            if P in Stool_MBX_input_members:
                tmp.append(P)
        Stool_MBX_input_members = tmp

        # 各患者から何回目に採取されたサンプルなのかをカウント
        Participant_memo = Stool_MBX_member_df['Participant ID'].values.tolist()[0]
        cnt = 0
        cnt_lst = []
        for row in Stool_MBX_member_df.values.tolist():
            if Participant_memo == row[0]:
                cnt += 1
                cnt_lst.append(cnt)
            else:
                Participant_memo = row[0]
                cnt = 1
                cnt_lst.append(cnt)
        Stool_MBX_member_df['Sampling_counter'] = cnt_lst
        # st.table(Stool_MBX_member_df.head(20))

        # 各被験者の採取回数を表示

        histgram = alt.Chart(Stool_MBX_member_df).mark_bar(opacity=0.5).encode(
            y=alt.Y("Participant ID:O", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='被験者'),
                sort='-x'
                ),
            x=alt.X('count(Participant ID)',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='サンプリング回数'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(
            width=500,height=1600
        ).interactive()

        bar = alt.Chart(Stool_MBX_member_df[Stool_MBX_member_df['Sampling_counter']==min_ns]).mark_bar(opacity=0.5).encode(
            y=alt.Y("count(Participant ID)", 
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='人数'),
                sort='-x'
                ),
            x=alt.X('diagnosis',
                axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='診断'),
                stack=None
                ),
            color=alt.Color('diagnosis'),
        ).properties(
            width=150,height=1600
        ).interactive()

        Stool_MBX_sample_dist_chart = histgram|bar

        # Stool_MBX_Analytics_exec = st.checkbox('Exec（メタボローム解析）')
        if Stool_MBX_Analytics_exec == True:

            # リストアップされた患者の代謝経路データを集める
            Stool_MBX_messy_df = pd.pivot_table(Stool_MBX_df, index=['Participant ID', 'week_num'], columns='Metabolites', values='Abundance')
            Stool_MBX_messy_df = Stool_MBX_messy_df.T[Stool_MBX_input_members].T
            Stool_MBX_messy_df = Stool_MBX_messy_df.reset_index()
            Stool_MBX_messy_df = pd.merge(Stool_MBX_messy_df,Stool_MBX_member_df,on=['Participant ID','week_num'],how='left')
            Stool_MBX_messy_df = Stool_MBX_messy_df[Stool_MBX_messy_df['Sampling_counter']<=min_ns]

            # ## 可視化する代謝経路の種類をリストにまとめてメモしておく
            Stool_MBX_messy_df = Stool_MBX_messy_df.set_index(['External ID','diagnosis'])
            Stool_MBX_sum_series = Stool_MBX_messy_df[Stool_MBX_messy_df['Participant ID']==Participant].drop(['Participant ID','week_num','Sampling_counter'],axis=1).sum(axis=0)
            Stool_MBX_lst = Stool_MBX_sum_series[Stool_MBX_sum_series > 0].index.tolist()
            Stool_MBX_messy_df = Stool_MBX_messy_df.reset_index()
            # st.write(Stool_MBX_lst)

            # # Dysbiosis score で解析サンプルをトリミングすることができる
            # if Stool_MGX_Microbes_Analytics_exec == True:
            #     Stool_MBX_messy_df = pd.merge(Stool_MBX_messy_df,Stool_MGX_Microbes_Level_messy_config_df[['Participant ID','week_num','Dysbiosis_Score']],on=['Participant ID','week_num'],how='left')
            #     if Stool_MGX_Microbes_Dysbiotic_CD == True:
            #         Stool_MBX_messy_df = Stool_MBX_messy_df[((Stool_MBX_messy_df['Dysbiosis_Score'] > select_CD_sample[0])&(Stool_MBX_messy_df['Dysbiosis_Score'] < select_CD_sample[1])&(Stool_MBX_messy_df['diagnosis'] == 'CD'))|(Stool_MBX_messy_df['diagnosis'] == 'UC')|(Stool_MBX_messy_df['diagnosis'] == 'nonIBD')]
            #     if Stool_MGX_Microbes_Dysbiotic_UC == True:
            #         Stool_MBX_messy_df = Stool_MBX_messy_df[((Stool_MBX_messy_df['Dysbiosis_Score'] > select_UC_sample[0])&(Stool_MBX_messy_df['Dysbiosis_Score'] < select_UC_sample[1])&(Stool_MBX_messy_df['diagnosis'] == 'UC'))|(Stool_MBX_messy_df['diagnosis'] == 'CD')|(Stool_MBX_messy_df['diagnosis'] == 'nonIBD')]
            #     Stool_MBX_messy_df = Stool_MBX_messy_df[Stool_MBX_messy_df['Dysbiosis_Score'] >= 0] # 同じタイミングで MGX が計測されていない場合は Dysbiosis score を計算できないので削除
            #     # st.table(Stool_MBX_messy_df.head())
            #     Stool_MBX_messy_df = Stool_MBX_messy_df.drop('Dysbiosis_Score',axis=1)

            #
            Stool_MBX_messy_df = Stool_MBX_messy_df.set_index(['External ID','diagnosis','Participant ID','week_num','Sampling_counter'])
            # Stool_MBX_messy_df = Stool_MBX_messy_df.drop(['Participant ID','week_num','Sampling_counter'],axis=1)
            # st.table(Stool_MBX_messy_df.head())

            # PLS 解析で入力する応答変数を作成
            Stool_MBX_messy_index_lst = Stool_MBX_messy_df.index.tolist()
            Stool_MBX_messy_diagnosis_array = np.array(Stool_MBX_messy_index_lst)[:,1]
            Stool_MBX_messy_diagnosis_dummy_array = pd.get_dummies(Stool_MBX_messy_diagnosis_array).values

            # PLS 実施
            from sklearn.cross_decomposition import PLSRegression
            pls = PLSRegression(n_components=2)
            pls_fit = pls.fit(X=Stool_MBX_messy_df.values,Y=Stool_MBX_messy_diagnosis_dummy_array)

            # VIP 計算 https://github.com/scikit-learn/scikit-learn/issues/7050#issuecomment-345208503
            def calc_VIP(model):
                t = model.x_scores_
                w = model.x_weights_
                q = model.y_loadings_
                p, h = w.shape
                vips = np.zeros((p,))
                s = np.diag(t.T @ t @ q.T @ q).reshape(h,-1)
                total_s = np.sum(s)
                for i in range(p):
                    weight = np.array([ (w[i,j]/np.linalg.norm(w[:,j]))**2 for j in range(h) ])
                    vips[i] = np.sqrt( p*(s.T @ weight)/total_s )
                return vips

            pls_vip = calc_VIP(pls_fit)
            Stool_MBX_ana_res_df = pd.DataFrame({'Metabolites': Stool_MBX_messy_df.columns.tolist(),'VIP': pls_vip.tolist()}) 

            # ランダムフォレスト実施
            from sklearn.ensemble import RandomForestClassifier
            clf_rf = RandomForestClassifier()
            clf_rf.fit(Stool_MBX_messy_df.values, Stool_MBX_messy_diagnosis_dummy_array)
            Stool_MBX_ana_res_df['MDI'] = clf_rf.feature_importances_  # 変数重要度 mean decrease impurity （https://qh73xebitbucketorg.readthedocs.io/ja/latest/1.Programmings/python/LIB/scikit-learn/EnsembleMethods/randomizedForest/）
            
            # 可視化するのは選択された個人がもっている代謝物のみ
            Stool_MBX_ana_res_df = Stool_MBX_ana_res_df[Stool_MBX_ana_res_df['Metabolites'].isin(Stool_MBX_lst)]

            # 可視化
            Stool_MBX_ana_res_chart = alt.Chart(Stool_MBX_ana_res_df).mark_point(
                filled=True, 
                size=200,
                opacity=0.7
            ).encode(
                x=alt.X('VIP:Q',
                    axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='VIP（変数重要度）'),
                    ),
                y=alt.Y('MDI:Q',
                    axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='MDI（変数重要度）'),
                    ),
                color=alt.condition(metabolites_select,'Metabolites:N',alt.value('lightgray'),legend=None),
                tooltip='Metabolites'
            ).properties(
                width=500, height=500
            ).add_selection(
                metabolites_select
            )

            # 可視化するのは選択された個人がもっている代謝物のみ
            Stool_MBX_messy_df = Stool_MBX_messy_df[Stool_MBX_lst]
            Stool_MBX_tidy_df = pd.melt(Stool_MBX_messy_df.reset_index(),id_vars=['External ID','diagnosis','Participant ID','week_num','Sampling_counter'],var_name="Metabolites",value_name="Abundance")

            jitter = alt.Chart(Stool_MBX_tidy_df).mark_circle(size=20).encode(
                x=alt.X('jitter:Q',
                    title=None,
                    axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                    scale=alt.Scale(),
                ),
                y=alt.Y('Abundance:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title='存在量')
                    ),
                color=alt.Color('Metabolites',legend=None),
                column=alt.Column('diagnosis:N', 
                    header=alt.Header(labelFontSize=15, titleFontSize=18,title='病名')
                    ),
                tooltip=['Participant ID','week_num']
            ).transform_calculate(
                # Generate Gaussian jitter with a Box-Muller transform
                jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
            ).add_selection(
                metabolites_select
            ).transform_filter(
                metabolites_select
            ).properties(
                width=59,height=500
            )

            boxplot = alt.Chart(Stool_MBX_tidy_df).mark_boxplot(
                opacity=0.5,
                size=30,
                ticks=alt.MarkConfig(width=30), 
                median=alt.MarkConfig(color='black',size=30),
                outliers=alt.MarkConfig(opacity=0)
            ).encode(
                y=alt.Y('Abundance:Q',
                    axis=alt.Axis(labelAngle=0,labelFontSize=15, ticks=True, titleFontSize=18, title=None)
                    ),
                color=alt.Color('Metabolites',legend=None),
                column=alt.Column('diagnosis:N', 
                header=alt.Header(labelFontSize=15, titleFontSize=18,title='病名')
                )
            ).add_selection(
                metabolites_select
            ).transform_filter(
                metabolites_select
            ).properties(
                width=59,height=500
            )

            Stool_MBX_ana_res_chart = Stool_MBX_ana_res_chart | jitter | boxplot
            Stool_MBX_chart = Stool_MBX_chart & Stool_MBX_ana_res_chart

        return Stool_MBX_chart, Stool_MBX_sample_dist_chart
    else:
        return Stool_MBX_chart

if Stool_MBX_Analytics == True :
    Stool_MBX_chart, Stool_MBX_sample_dist_chart = Stool_MBX_chart()
    with st.expander('各患者のサンプリング回数'):
        st.write(Stool_MBX_sample_dist_chart)

    if Stool_MBX_Analytics_exec == True: 

        text = '''

            解析結果が下に追加されます。下左図の各プロットは代謝物に対応しています。カーソルを合わせると代謝物名を確認できます。
            横軸は **PLS-DA** で計算される変数重要度 **VIP（Variable Importance in Projection）** を表し、縦軸は Random Forest で計算される変数重要度 **MDI（Mean Decrease in Impurity）** を表しています。
            いずれの指標も大きいほど CD と UC と nonIBD を区別する際に重要であるということです。\n
            下右図では、代謝経路ごとに存在量の分布が箱ひげ図で表示されます。変数重要度では群間の大小関係はわからないので、CD、UC、nonIBD の大小比較に用いるとよいでしょう。\n
            着目したい代謝物をクリックすると他の種が灰色で表示されます。上の時系列データとも連動しています。色々クリックしてみましょう。

            '''
    else:
        text = ''

    st.write(Stool_MBX_chart)
    st.markdown(text)

else:
    Stool_MBX_chart = Stool_MBX_chart()
    st.write(Stool_MBX_chart)

st.markdown('# 細菌と代謝物の関連について')

text = '''

糞便に含まれる細菌と代謝物の関連を調べます。ここでは**スピアマンの相関係数**と [**mmvec**](https://www.nature.com/articles/s41592-019-0616-3?proof=t) によって推定される**共起確率**が表示されます。
共起確率は全被験者のデータに基づきあらかじめ算出された結果が表示されますが、相関係数はサイドバーで選択された被験者から算出された結果が表示されます。

'''
st.markdown(text)

with st.expander('共起確率と相関係数について'):
    text = '''
    共起確率とは、細菌と代謝物のどちらか一方が観察されたとき、もう一方も観察される確率のことです。
    したがって細菌と代謝物のどちらか一方、あるいは両方の存在量が少ない場合、共起確率は低くなります。
    すなわち共起確率は、存在量の多いメジャーな細菌や代謝物であるほど高くなる傾向があります。
    そのため、共起確率が高いメジャーな細菌と代謝物の組み合わせでも相関係数が低かったり、
    共起確率が低いマイナーな細菌と代謝物の組み合わせでも相関係数が高かったりすることもあります。

    '''
    st.markdown(text)

Level = st.selectbox('細菌の分類階級を選んでください',Level_lst, index=5,key='Stool_MGX_MBX_select_Level')

with st.expander('各患者のサンプル抽出数'):
    text = '''

    1 年間で糞便を提供した回数は被験者によって様々です。各被験者から何回ずつ抽出されたサンプルで解析するのかをスライダーで選択してください。サンプル提供回数がスライダーで選択された回数に満たない被験者は解析から除外されます。\n
    ``各患者のサンプリング回数`` をクリックすると、左側に各被験者のサンプル提供回数が表示され、右側にスライダーで選択された回数以上サンプルを提供した被験者の人数が表示されます。\n
    ここではメタゲノムデータとメタボロームデータが揃っているサンプルのみを数えています。
    
    '''

min_ns = st.slider('各患者のサンプル抽出数を選んでください',min_value=1, max_value=7,step=1, key='Stool_MGX_MBX_min_ns') # minimum number of samples

@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def Stool_MGX_MBX():
    Stool_MGX_MBX_Level_df = Stool_MGX_MBX_df[(Stool_MGX_MBX_df['genre']==Level)|(Stool_MGX_MBX_df['genre']=='metabolites')]

    # メタデータを付与（診断、性別、年齢など sidebar で選択できるもの）
    sn = Stool_MGX_MBX_Level_df.columns.tolist()[:-1] # Sample Name
    select_md = Metadata_df[Metadata_df['site_sub_coll'].isin(sn)][['Participant ID','site_sub_coll','diagnosis','sex','consent_age','Education Level','Occupation']]
    select_md = select_md.set_index('site_sub_coll')
    select_md = select_md.loc[~select_md.index.duplicated(),:]

    # 各患者から何回目に採取されたサンプルなのかをカウント
    Participant_memo = select_md['Participant ID'].values.tolist()[0]
    cnt = 0
    cnt_lst = []
    for row in select_md['Participant ID'].values.tolist():
        if Participant_memo == row:
            cnt += 1
            cnt_lst.append(cnt)
        else:
            Participant_memo = row
            cnt = 1
            cnt_lst.append(cnt)
    select_md['Sampling_counter'] = cnt_lst

    # サイドバーに基づいてトリミング
    if select_sex_exec == True:
        select_md = select_md[select_md['sex'].isin(select_sex)]
    if select_age_exec == True:
        select_md = select_md[select_md['consent_age'] >= select_age[0]]
        select_md = select_md[select_md['consent_age'] <= select_age[1]]
    if select_EL_exec == True:
        select_md = select_md[select_md['Education Level'].isin(select_EL)]
    if select_OL_exec == True:
        select_md = select_md[select_md['Occupation'].isin(select_OL)]

    # サンプリング回数の描画   
    histgram = alt.Chart(select_md).mark_bar(opacity=0.5).encode(
        y=alt.Y("Participant ID:O", 
            axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='被験者'),
            sort='-x'
            ),
        x=alt.X('count(Participant ID)',
            axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='サンプリング回数'),
            stack=None
            ),
        color=alt.Color('diagnosis'),
    ).properties(
        width=500,height=1600
    ).interactive()

    bar = alt.Chart(select_md[select_md['Sampling_counter']==min_ns]).mark_bar(opacity=0.5).encode(
        y=alt.Y("count(Participant ID)", 
            axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18, title='人数'),
            sort='-x'
            ),
        x=alt.X('diagnosis',
            axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18,title='診断'),
            stack=None
            ),
        color=alt.Color('diagnosis'),
    ).properties(
        width=150,height=1600
    ).interactive()

    Stool_MGX_MBX_sample_dist_chart = histgram|bar

    # スライドバーで選択したサンプリング回数だけ各被験者からサンプルを抽出して相関係数を計算
    select_md = select_md[select_md['Sampling_counter']<=min_ns]
    Stool_MGX_MBX_Level_df = Stool_MGX_MBX_Level_df[select_md.index.tolist()+['genre']]

    # 細菌と代謝物の相関関係を調べる
    corr = Stool_MGX_MBX_Level_df.drop('genre',axis=1).T.corr(method='spearman')
    mgx_lst = Stool_MGX_MBX_Level_df[Stool_MGX_MBX_Level_df['genre']==Level].index.tolist()
    mbx_lst = Stool_MGX_MBX_Level_df[Stool_MGX_MBX_Level_df['genre']=='metabolites'].index.tolist()
    corr = corr[mbx_lst].T[mgx_lst]
    # print(corr)

    # mmvec データ読取
    mmvec_path = ['InputData','mmvec','MGX-MBX','output',Level,'latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ranks.txt']
    cwd = os.getcwd()
    mmvec_file = os.path.join(cwd, *mmvec_path)
    mmvec = pd.read_table(mmvec_file, delimiter='\t',index_col=0)

    # corr データから mmvec データのない細菌と代謝物を除く
    mgx_list = mmvec.columns.tolist()
    mbx_list = mmvec.index.tolist()
    corr = corr[mgx_list].T[mbx_list].T
    # print(corr)

    # tidy データへの変換
    corr['MBX'] = corr.index
    corr_tidy = pd.melt(corr,id_vars="MBX",var_name="MGX",value_name="corr" )

    # corr データ と mmvec データの結合
    corr = pd.DataFrame()
    corr_list = []
    mm_list = []
    metabo_list = []
    micro_list = []
    for row in corr_tidy.values.tolist():
        metabo_list.append(row[0])
        micro_list.append(row[1])
        mm_list.append('{}__AND__{}'.format(row[0],row[1]))
        corr_list.append(row[2])
    corr['mm'] = mm_list
    corr['corr'] = corr_list
    corr['MBX'] = metabo_list
    corr['MGX'] = micro_list

    mmvec['MBX'] = mmvec.index
    mmvec_tidy = pd.melt(mmvec,id_vars="MBX",var_name="MGX",value_name="corr" )
    mmvec = pd.DataFrame()
    mmvec_list = []
    mm_list = []
    for row in mmvec_tidy.values.tolist():
        mm_list.append('{}__AND__{}'.format(row[0],row[1]))
        mmvec_list.append(row[2])
    mmvec['mm'] = mm_list
    mmvec['mmvec'] = mmvec_list
    corr_mmvec = pd.merge(mmvec,corr,how='left',on='mm')
    return corr_mmvec, select_md, Stool_MGX_MBX_Level_df, mgx_list, mbx_list, Stool_MGX_MBX_sample_dist_chart

# mmvec score と相関係数の散布図を作成
@st.cache(allow_output_mutation=True, suppress_st_warning=True, max_entries=10, ttl=3600)
def corr_mmvec_scatter(data):
    selected = alt.selection_single(on="click", empty="none", fields=['mm'], init={'mm': corr_mmvec['mm'].values.tolist()[0]})

    return alt.Chart(data).mark_circle(size=30).encode(
            x=alt.X('mmvec:Q',axis=alt.Axis(labelFontSize=15, titleFontSize=18, labels=True,ticks=True,grid=False,title='mmvec で推定される細菌と代謝物の共起確率の CLR 変換値')),
            y=alt.Y('corr:Q',axis=alt.Axis(labelFontSize=15, titleFontSize=18, labels=True,ticks=True,grid=False,title='細菌と代謝物のスピアマン相関係数'),scale=alt.Scale(domain=[-1, 1])),
        color=alt.condition(selected, alt.value("steelblue"), alt.value("lightgray")),
        tooltip=['MBX','MGX']
    ).add_selection(selected).properties(
        width=400,
        height=400,
    )

st.markdown('``Visualize`` にチェックを入れると解析結果の可視化を準備できます。')
Stool_MGX_MBX_mmvec_exec = st.checkbox('Visualize')
if Stool_MGX_MBX_mmvec_exec == True: 

    corr_mmvec, select_md, Stool_MGX_MBX_Level_df, mgx_list, mbx_list, Stool_MGX_MBX_sample_dist_chart = Stool_MGX_MBX()

    with st.expander('各患者のサンプリング回数'):
        st.write(Stool_MGX_MBX_sample_dist_chart)

    # 可視化する細菌と代謝物を選択する

    st.markdown('着目したい細菌または代謝物があれば、ここで選択することができます。')

    select_MGX = st.multiselect('細菌を選んでください',sorted(mgx_list))
    if len(select_MGX) == 0:
        select_MGX = mgx_list
    select_MBX = st.multiselect('代謝物を選んでください',sorted(mbx_list))
    if len(select_MBX) == 0:
        select_MBX = mbx_list

    st.markdown('``Exec`` にチェックを入れると解析結果が可視化されます。')

    Stool_MGX_MBX_Analytics_exec = st.checkbox('Exec', key='Stool_MGX_MBX_Analytics_exec')
    if Stool_MGX_MBX_Analytics_exec == True:
        mmvec_col1, mmvec_col2 = st.columns(2)

        with mmvec_col1:
            # Altair が描画されると共に Altair 内で選択した行の情報が dict 型で返される。
            # species を選択した際に図の容量が大きくなりすぎるのを防ぐために mmvec score が 0 以上のものを描画

            Stool_MGX_MBX_mmvec = altair_component(altair_chart=corr_mmvec_scatter(corr_mmvec[((corr_mmvec['MGX'].isin(select_MGX))&(corr_mmvec['MBX'].isin(select_MBX)))]))

            text = '''

            上図が解析結果です。各プロットは細菌と代謝物のペアに対応しています。カーソルを合わせると細菌名と代謝物名が表示されます。
            縦軸はスピアマンの相関係数を表し、横軸は mmvec で推定される共起確率の CLR（Center logratio）変換値を表しています。
            CLR 変換値は平均で 0 になる対数値です。ここでは CLR 変換値が高いペアの共起確率は他のペアよりも高いといえます。\n
            着目したい細菌と代謝物のペアをクリックすると、右に分布が表示されます。上記のプルダウンから細菌または代謝物を選択すると、探索しやすくなります。

            '''
            
            st.markdown(text)

        if 'mm' in list(Stool_MGX_MBX_mmvec.keys()): # 何か選択されると 'name' key が加わる
            select_mm = list(Stool_MGX_MBX_mmvec['mm'])[0].split('__AND__')
            Stool_MGX_MBX_Level_df = select_md[['diagnosis','Participant ID']].join(Stool_MGX_MBX_Level_df.T)


            # 選択された細菌と代謝物の組み合わせで生データを可視化
            Stool_MGX_MBX_chart = alt.Chart(Stool_MGX_MBX_Level_df).mark_circle(
                size=30
            ).encode(
                x=alt.X(select_mm[0],
                    axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18)
                    ),
                y=alt.Y(select_mm[1],
                    axis=alt.Axis(labelFontSize=15, ticks=True, titleFontSize=18)
                    ),
                color=alt.Color('diagnosis'),
                tooltip=['Participant ID',select_mm[0],select_mm[1]]
            ).properties(
                width=500,height=500
            ).interactive()

            with mmvec_col2:

                st.write(Stool_MGX_MBX_chart)

                text = '''

                左図で選択した細菌と代謝物のペアの相関を見ることができます。カーソルを合わせると被験者を判別できます。

                '''
                st.markdown(text)