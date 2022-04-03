import pandas as pd
import numpy as np
import os
import re
import plotly.express as px
import tqdm
import random

from structuremap.plotting import scale_pvals

def test_identical_ids(directory_cif, directory_pae):
    p_list_cif = list()
    for file in os.listdir(directory_cif):
        if file.endswith("cif"):
            protein_id = re.sub(r'.cif', '', file)
            p_list_cif.append(protein_id)

    p_list_pae = list()
    for file in os.listdir(directory_pae):
        if file.endswith("hdf"):
            protein_id = re.sub(r'.hdf', '', file)
            protein_id = re.sub(r'pae_', '', protein_id)
            p_list_pae.append(protein_id)

    p_list_cif = np.unique(p_list_cif)
    p_list_pae = np.unique(p_list_pae)

    np.testing.assert_equal(len(p_list_cif),len(p_list_pae))
    np.testing.assert_equal(p_list_cif,p_list_pae)

    print('Number of unique proteins with cif and pae file: ',len(p_list_cif))

def plot_enrichment_david(data, fdr_threshold=0.01, fold_enrichment_threshold=2, count_threshold=10):
    df = data.copy(deep=True)
    df = df[df['Fold Enrichment'] >= fold_enrichment_threshold]
    df = df[df['Count'] >= count_threshold]
    df = df[df['FDR'] <= fdr_threshold]

    df['Term'] = [re.sub('GO:.*~','',x) for x in df['Term']]

    df['neg_log_FDR'] = -np.log10(df.FDR)

    category_dict = {}

    df['neg_log_FDR_round'] = scale_pvals(df.neg_log_FDR)

    category_dict['neg_log_FDR_round'] = list(reversed(['> 100','> 50','> 10','> 5','> 2','> 0']))
    color_dict = {'> 100':'rgb(177, 63, 100)',
                  '> 50':'rgb(221, 104, 108)',
                  '> 10':'rgb(241, 156, 124)',
                  '> 5':'rgb(245, 183, 142)',
                  '> 2':'rgb(246, 210, 169)',
                  '> 0':'grey'}

    fig = px.bar(df,
                 x='Fold Enrichment',
                 y='Term',
                 text='Count',
                 orientation='h',
                 color='neg_log_FDR_round',
                 hover_data=['Fold Enrichment','FDR','Count'],
                 #category_orders=category_dict,
                 category_orders=category_dict,
                 color_discrete_map = color_dict,
                 template="simple_white",
                 #width=900, height=600)
                )

    fig.update_yaxes(title_text='')
    fig.update_layout(margin = dict(l=400, r=150, b=10), legend=dict(title='-log10(FDR)'))

    config={'toImageButtonOptions': {'format': 'svg', 'filename':'structure ptm enrichment'}}

    return(fig.show(config=config))

def get_loop_positions(x : str) -> list:
    x_split = re.sub(' ', '', x)
    x_split = re.sub('\\[N/A\\]', 'nan', x)
    x_split = x_split.split(',')
    x_split = [y.split('...') for y in x_split]
    return x_split

def test_get_loop_positions():
    assert get_loop_positions("195...205, 366...378") == [['195', '205'], [' 366', '378']]
    assert get_loop_positions("195...205") == [['195', '205']]
    assert get_loop_positions("[N/A]") == [['nan']]

test_get_loop_positions()

def extract_loop_annotation(
    loop_df : pd.DataFrame,
) -> list:

    df = loop_df.copy(deep=True)

    loop_list = list()

    for col in ['gloop','aloop','achelix']:
        df[col] = df[col].apply(get_loop_positions)
        df[col] = df[col].apply(lambda x: [list(np.arange(int(p[0]),int(p[1])+1)) if p[0] != 'nan' else [] for p in x])
        df[col] = df[col].apply(lambda x: [item for sublist in x for item in sublist])

        position_list = list()
        protein_list = list()
        df.apply(lambda x: position_list.append(x[col]), axis=1)
        df.apply(lambda x: protein_list.append(list(np.repeat(x['uniprot_id'],len(x[col])))), axis=1)

        position_list = [item for sublist in position_list for item in sublist]
        protein_list = [item for sublist in protein_list for item in sublist]

        loop_df = pd.DataFrame({'protein_id': protein_list, 'position': position_list, col : 1})
        loop_list.append(loop_df)

    return loop_list


def annotate_ptm_data(alphafold_data: pd.DataFrame) -> pd.DataFrame:

    alphafold_data_annotated = alphafold_data.copy(deep=True)
    ptm_dir = 'data/ptm_data/'

    for file in tqdm.tqdm(os.listdir(ptm_dir)):
        print(file)
        if file == "phosphositeplus_annotation.csv":
            ptm_file = pd.read_csv(ptm_dir+file, skiprows=1)
        else:
            ptm_file = pd.read_csv(ptm_dir+file)
        if 'AA' in ptm_file.columns:
            alphafold_data_annotated = alphafold_data_annotated.merge(ptm_file, how='left', on=['protein_id','AA','position'])
        else:
            alphafold_data_annotated = alphafold_data_annotated.merge(ptm_file, how='left', on=['protein_id','position'])
        alphafold_data_annotated = alphafold_data_annotated.fillna(0)

    return(alphafold_data_annotated)

def generate_ptm_site_dict(alphafold_df):
    all_ptm_datasets = ['ub_treated_only', 'ac', 'ac_reg', 'ga', 'gl', 'gl_reg', 'm', 'm_reg',
       'p', 'p_reg', 'sm', 'sm_reg', 'ub', 'ub_reg',
       'p_functional_0', 'p_functional_5', 'p_functional_6', 'p_functional_7',
       'p_functional_8', 'p_functional_9', 'p_stukalov',
       'ub_shared', 'ub_untreated_only', 'p_sugiyama', 'p_sugiyama_psp',
       'p_sugiyama_ochoa','p_sugiyama_stukalov']
    ptm_dict = {}
    for d in all_ptm_datasets:
        df_d = alphafold_df[alphafold_df[d] == 1]
        unique_aa = list(np.unique(df_d.AA.values))
        ptm_dict.update({d: unique_aa})

    return(ptm_dict)

def plot_shortIDR_activationLoop(df, kinase):

    ripk2 = df[df.protein_id==kinase].reset_index(drop=True)
    ripk2 = ripk2.assign(
        pep_sym='square',
        p_sym='circle',
        idr_y = 0)

    ripk2['IDR_col'] = np.where((ripk2['IDR']==1), 'idr', 'structured')
    ripk2['p_reg_col'] = np.where((ripk2['p_reg']==1), 'regulatory', 'unknown')
    ripk2['p_reg_y'] = np.where((ripk2['p_reg']==1), 0.2, 0.12)

    ripk2['idr_count'] = 0
    idr_count=0
    for i in range(ripk2.shape[0]-1):
        if (ripk2['flexible_pattern_extended_5'].values[i]==0) & (ripk2['flexible_pattern_extended_5'].values[i+1]==1):
            idr_count=idr_count+1
        if (ripk2['flexible_pattern_extended_5'].values[i]==1):
            ripk2['idr_count'].values[i] = idr_count

    symbol_dict = {'square': 1, 'circle':0}
    color_map ={'regulatory': "#DC143C", 'unknown':"#FF9A9F", 'idr':"#B5B5B5", 'structured': "#367BC3"}

    fig_idr0 = px.scatter(ripk2[ripk2.IDR==0], x='position', y='idr_y', color='IDR_col', symbol='pep_sym',
                          symbol_map=symbol_dict, opacity=0.8, color_discrete_map=color_map)
    fig_idr1 = px.scatter(ripk2[ripk2.IDR==1], x='position', y='idr_y', color='IDR_col', symbol='pep_sym',
                          symbol_map=symbol_dict, opacity=0.8, color_discrete_map=color_map)

    fig = fig_idr0.add_traces(list(fig_idr1.select_traces()))

    df_plot_ptm = ripk2[ripk2.p==1].reset_index(drop=True)
    fig_p = px.scatter(df_plot_ptm,
                       x='position', y='p_reg_y', color='p_reg_col',
                       symbol='p_sym', symbol_map=symbol_dict,
                       color_discrete_map=color_map)
    fig = fig.add_traces(list(fig_p.select_traces()))

    for i in range(0, df_plot_ptm.shape[0]):
        fig = fig.add_shape(
                dict(
                    type="line",
                    x0=df_plot_ptm.position.values[i],
                    y0=0.02,
                    x1=df_plot_ptm.position.values[i],
                    y1=df_plot_ptm.p_reg_y.values[i] - 0.02, #0.18,
                    line=dict(
                        color='black',
                        width=1
                    )
                )
        )

    if ripk2[ripk2.flexible_pattern_extended_5==1].shape[0] > 0:
        for c in np.arange(1,idr_count+1):
            fig = fig.add_shape(
                type="line",
                x0=np.min(ripk2[(ripk2.idr_count==c)].position),
                x1=np.max(ripk2[(ripk2.idr_count==c)].position),
                y0=-0.1,
                y1=-0.1,
                line=dict(color="#66985E",width=4),
                opacity=1
            )

            fig = fig.add_annotation(dict(font=dict(color='#66985E', size=9),
                                                    x=np.mean(ripk2[ripk2.idr_count==c].position),
                                                    y=-0.2,
                                                    showarrow=False,
                                                    text="Short IDR",
                                                    textangle=0,
                                                    xanchor='center'))


    if ripk2[ripk2.aloop==1].shape[0] > 0:
        fig = fig.add_shape(
            type="line",
            x0=np.min(ripk2[ripk2.aloop==1].position),
            x1=np.max(ripk2[ripk2.aloop==1].position),
            y0=0.3,
            y1=0.3,
            line=dict(color="#003500",width=4),
            opacity=1
        )

        fig = fig.add_annotation(dict(font=dict(color='#003500', size=9),
                                                x=np.mean(ripk2[ripk2.aloop==1].position),
                                                y=0.4,
                                                showarrow=False,
                                                text="A-loop",
                                                textangle=0,
                                                xanchor='center'))

    fig = fig.update_layout(title=ripk2.protein_id[0],
                            template="simple_white",
                            xaxis={'visible': False, 'showticklabels': False},
                            yaxis={'visible': False, 'showticklabels': False},
                            legend=dict(
                                yanchor="bottom",
                                y=-1,
                                xanchor="left",
                                x=0,
                                orientation='h'),
                            yaxis_range=[-0.3,0.5],
                            width=900, height=250,
                           )

    config={'toImageButtonOptions': {'format': 'svg', 'filename':'short_IDR_plot_'+ripk2.protein_id[0]}}

    return(fig.show(config=config))

def format_for_3Dviz(df: pd.DataFrame,
                     ptm_dataset: str) -> pd.DataFrame:
    df_mod = df[["protein_id","AA","position",ptm_dataset]]
    df_mod = df_mod.rename(columns={"protein_id": "unique_protein_id",
                                    "AA": "modified_sequence",
                                    "position": "start"})
    df_mod["modified_sequence"] = [mod+"_"+str(i) for i,mod in enumerate(df_mod["modified_sequence"])]
    df_mod["all_protein_ids"] = df_mod["unique_protein_id"]
    df_mod["PTMsites"] = 0
    df_mod["start"] = df_mod["start"]-1
    df_mod["end"] = df_mod["start"]
    df_mod["PTMsites"] = [[i] for i in df_mod["PTMsites"]]
    df_mod = df_mod[df_mod[ptm_dataset] == 1]
    df_mod["marker_symbol"] = 1
    df_mod["PTMtypes"] = [[ptm_dataset] for i in df_mod["PTMsites"]]
    df_mod = df_mod.dropna(subset=['PTMtypes']).reset_index(drop=True)
    return df_mod

def string_to_motif(start_motif : str) -> str:
    motif = re.sub('X','A-Z',start_motif)
    motif = ''.join(['['+i+']' for i in motif.split(',')])
    motif = re.sub('\\[nan\\]','',motif)
    return motif

def test_string_to_motif():
    assert string_to_motif('R,X,R,X,X') == '[R][A-Z][R][A-Z][A-Z]'
    assert string_to_motif('nan') == ''
    assert string_to_motif('S,T,nan') == '[S][T]'

test_string_to_motif()


def get_kinase_substrates(alphafold_df,
                          kinase_df,
                          kinase,
                          max_pPSE = None,
                          write_output=True,
                          name="kinase_substrates",
                          window_size=6,
                          random_seed=44):
    random.seed(random_seed)
    kinase_list = list()
    inside_list = list()
    cutoff_list = list()
    detected_count=0
    missed_count=0
    short_count=0
    high_ppse_count=0
    kinase_sub = kinase_df[kinase_df.Kinase==kinase].reset_index(drop=True)
    for i in tqdm.tqdm(np.arange(0, kinase_sub.shape[0])):
        df_sub = alphafold_df[
            (alphafold_df.protein_id==kinase_sub.protein_id[i]) &
            (alphafold_df.position.isin(range(kinase_sub.position[i]-window_size,
                                              kinase_sub.position[i]+window_size+1)))
        ]
        sub = df_sub.AA.values
        if len(sub) > 0:
            if sub[6] == kinase_sub.AA[i]:
                if len(sub) == (2*window_size)+1:
                    sub = ''.join(sub)
                    kinase_list.append(sub)
                    if max_pPSE:
                        ppse = df_sub.nAA_12_70_pae.values[window_size]
                        if ppse <= max_pPSE:
                            inside_list.append(sub)
                            detected_count+=1
                        else:
                            cutoff_list.append(sub)
                            high_ppse_count+=1
                    else:
                        detected_count+=1
                else:
                    short_count+=1
            else:
                missed_count+=1
        else:
            missed_count+=1

    print("Detected: ",detected_count,
          " Missed: ", missed_count,
          " Too short: ", short_count,
          " pPSE cutoff: ",high_ppse_count)

    if write_output:
        textfile = open(name+"_"+kinase+"_all.txt", "w")
        for element in kinase_list:
            textfile.write(element + "\n")

        textfile = open(name+"_"+kinase+"_maxPPSE_"+str(max_pPSE)+"_in.txt", "w")
        for element in inside_list:
            textfile.write(element + "\n")

        if len(cutoff_list)>0:
            textfile = open(name+"_"+kinase+"_maxPPSE_"+str(max_pPSE)+"_out.txt", "w")
            for element in cutoff_list:
                textfile.write(element + "\n")

            random_inside_subset = random.sample(inside_list, len(cutoff_list))
            textfile = open(name+"_"+kinase+"_maxPPSE_"+str(max_pPSE)+"_in_random_sub.txt", "w")
            for element in random_inside_subset:
                textfile.write(element + "\n")

    return [kinase_list, inside_list, cutoff_list, random_inside_subset]

def extract_short_IDR_info(df):
    shortIDR_sub = df[(df.flexible_pattern_extended_5==1)][["protein_id","AA","position"]].reset_index(drop=True)
    proteins = list()
    short_idr_extended_start = list()
    short_idr_extended_end = list()
    short_idr_sequence = list()
    seq_idx = 0
    for i in range(0,shortIDR_sub.shape[0]):
        if i==0:
            proteins.append(shortIDR_sub.protein_id[i])
            short_idr_extended_start.append(shortIDR_sub.position[i])
            short_idr_sequence.append(shortIDR_sub.AA[i])
        elif i==shortIDR_sub.shape[0]-1:
            short_idr_extended_end.append(shortIDR_sub.position[i])
            short_idr_sequence[seq_idx] = short_idr_sequence[seq_idx]+shortIDR_sub.AA[i]
        else:
            if shortIDR_sub.protein_id[i] != shortIDR_sub.protein_id[i-1]:
                proteins.append(shortIDR_sub.protein_id[i])
                short_idr_extended_end.append(shortIDR_sub.position[i-1])
                short_idr_extended_start.append(shortIDR_sub.position[i])
                seq_idx +=1
                short_idr_sequence.append(shortIDR_sub.AA[i])
            elif shortIDR_sub.protein_id[i] == shortIDR_sub.protein_id[i-1]:
                if shortIDR_sub.position[i] != shortIDR_sub.position[i-1]+1:
                    proteins.append(shortIDR_sub.protein_id[i])
                    short_idr_extended_end.append(shortIDR_sub.position[i-1])
                    short_idr_extended_start.append(shortIDR_sub.position[i])
                    seq_idx += 1
                    short_idr_sequence.append(shortIDR_sub.AA[i])
                else:
                    short_idr_sequence[seq_idx] = short_idr_sequence[seq_idx]+shortIDR_sub.AA[i]
    res_df = pd.DataFrame({'protein_id':proteins,
                           'IDR_start':short_idr_extended_start,
                           'IDR_end':short_idr_extended_end,
                           'sequence':short_idr_sequence})
    return res_df
