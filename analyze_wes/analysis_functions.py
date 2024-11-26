import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd
import seaborn as sns
from tqdm import tqdm

colors_40 = ['#a9a9a9', '#2f4f4f', '#556b2f', '#a0522d', '#228b22', '#191970', '#8b0000', '#808000', '#4682b4', '#9acd32',
             '#20b2aa', '#00008b', '#32cd32', '#daa520', '#8fbc8f', '#800080', '#b03060', '#ff4500', '#ff8c00', '#6a5acd',
             '#ffff00', '#0000cd', '#00ff00', '#9402d3', '#ba55d3', '#dc143c', '#00bfff', '#f4a460', '#adff2f', '#ff00ff',
             '#1e90ff', '#f0e68c', '#fa8072', '#dda0dd', '#afeeee', '#98fb98', '#7fffd4', '#ffdab9', '#ff69b4', '#ffb6c1']

colors_80 = ['#696969', '#a9a9a9', '#d3d3d3', '#2f4f4f', '#556b2f', '#6b8e23', '#a0522d', '#2e8b57', '#228b22', '#7f0000', 
             '#191970', '#006400', '#808000', '#483d8b', '#b22222', '#5f9ea0', '#778899', '#3cb371', '#bc8f8f', '#663399',
             '#008080', '#bdb76b', '#4682b4', '#d2691e', '#9acd32', '#20b2aa', '#cd5c5c', '#00008b', '#4b0082', '#32cd32',
             '#daa520', '#8fbc8f', '#800080', '#b03060', '#66cdaa', '#9932cc', '#ff4500', '#ff8c00', '#ffa500', '#ffd700',
             '#ffff00', '#c71585', '#0000cd', '#deb887', '#40e0d0', '#7fff00', '#00ff00', '#9400d3', '#ba55d3', '#00fa9a',
             '#00ff7f', '#4169e1', '#dc143c', '#00ffff', '#00bfff', '#f4a460', '#9370db', '#0000ff', '#adff2f', '#d8bfd8',
             '#b0c4de', '#ff7f50', '#ff00ff', '#db7093', '#fa8072', '#eee8aa', '#ffff54', '#6495ed', '#dda0dd', '#b0e0e6',
             '#ff1493', '#7b68ee', '#ffa07a', '#ee82ee', '#98fb98', '#87cefa', '#7fffd4', '#ffdab9', '#ff69b4', '#ffc0cb']

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

colorblind = sns.color_palette("colorblind", 9).as_hex()
colorblind_hex_codes = ['#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133']

linestyles = ['-', ':', '--']


# display_colors(colors_40) 
# display_colors(colors_80)
# display_colors(default_colors)
# display_colors(colorblind)


def display_colors(colors_list):

    num_colors = len(colors_list)
    
    cmap = ListedColormap(colors_list)
    vals = range(len(colors_list))
    bounds = np.append(vals, vals[-1] + 1)
    norm = BoundaryNorm(bounds, ncolors=num_colors)
    
    fig, ax = plt.subplots(figsize=(0.35 * num_colors, 1))
    fig.subplots_adjust(bottom=0.5)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                 cax=ax, orientation='horizontal')

    for i in range(num_colors):
        ax.text(i + 0.35 / 2, 0.5, i)

    ax.axes.get_xaxis().set_visible(False)
    
    plt.show()



def pickle_that(obj, fp):
    
    with open(fp, 'wb') as out_handle:
        pickle.dump(obj, out_handle)


def unpickle_that(fp):
    
    with open(fp, 'rb') as in_handle:
        return pickle.load(in_handle)
        

def build_paths(paths):
    
    for path in paths:
        
        if not os.path.exists(path):
            build_path(path)
            print(f"Created directory: {path}")


def build_path(path):
    os.makedirs(path)


def build_oncoprint_df(mutations_dict, gene_set, pt_sample_ordering):

    track_types = ['frameshift_indel', 'missense_mutation', 'nonsense_mutation', 'silent_intron_igr', 'splice_utr', 'inframe_indel']

    # midx = pd.MultiIndex.from_product([g])
    idx = list(range(len(gene_set) * len(track_types)))

    op_df = pd.DataFrame('', index=idx, columns=['track_name', 'track_type'] + pt_sample_ordering)

    i = 0
    for track_type in track_types:
        
        for gene in gene_set:
            
            op_df.loc[i, 'track_name'] = gene
            op_df.loc[i, 'track_type'] = track_type
            
            i += 1
            
    op_df.index = op_df['track_name'] + '_' + op_df['track_type']

    for gene in tqdm(gene_set):

        for pt in mutations_dict:
        
            for sample in mutations_dict[pt]:
                
                sample_gene_mutations_df = mutations_dict[pt][sample].loc[mutations_dict[pt][sample]['Gene_Name'] == gene]
                
                for _, row in sample_gene_mutations_df.iterrows():
                    
                    found = False
                                
                    if any([var_type in row['Mutation_Type'] for var_type in ['intron_variant', 'synonymous_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'sequence_feature', 'protein_protein_contact', 'intergenic_region']]):
                        op_df.loc[f"{gene}_silent_intron_igr", sample] = 'Silent/Intron/IGR'
                        found = True
                        
                    if any([var_type in row['Mutation_Type'] for var_type in ['splice_region_variant', '5_prime_UTR_premature_start_codon_gain_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_region_variant&intron_variant']]):
                        op_df.loc[f"{gene}_splice_utr", sample] = 'Splice/UTR'
                        found = True

                    if any([var_type in row['Mutation_Type'] for var_type in ['frameshift_variant']]):
                        op_df.loc[f"{gene}_frameshift_indel", sample] = 'Frameshift Indel'
                        found = True
                        
                    if any([var_type in row['Mutation_Type'] for var_type in ['conservative_inframe_deletion', 'conservative_inframe_insertion', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion']]):
                        op_df.loc[f"{gene}_inframe_indel", sample] = 'Inframe Indel'
                        found = True
                        
                    if any([var_type in row['Mutation_Type'] for var_type in ['missense_variant']]):
                        op_df.loc[f"{gene}_missense_mutation", sample] = 'Missense Mutation'
                        found = True
                        
                    if any([var_type in row['Mutation_Type'] for var_type in ['stop_gained', 'stop_lost']]):
                        op_df.loc[f"{gene}_nonsense_mutation", sample] = 'Nonsense Mutation'
                        found = True
                        
                    if any([var_type in row['Mutation_Type'] for var_type in ['structural_interaction_variant']]):
                        
                        if len(row['REF']) == len(row['ALT']):
                            op_df.loc[f"{gene}_silent_intron_igr", sample] = 'Silent/Intron/IGR'
                            
                        else:
                            op_df.loc[f"{gene}_inframe_indel", sample] = 'Inframe Indel'
        
                        found = True
                        
                    if not found:
                        print(sample_gene_mutations_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Gene_Name', 'Mutation_Type']].sort_values(by=['CHROM', 'POS']).values)

    op_df = op_df.drop(columns=['track_type']).set_index('track_name').fillna('')

    return op_df


def create_mutation_markers(op_df):

    mutation_markers = {}

    for i, category in enumerate(np.sort(np.unique(op_df.values.flatten()))):
        
        if category == '':
            continue
            
        mutation_markers[category] = dict(color=colorblind[i - 1])
        
        if category in ['Missense Mutation', 'Frameshift Indel']:
            mutation_markers[category]['marker'] = 'rect'
            mutation_markers[category]['zindex'] = 0
            
        elif category in ['Silent/Intron/IGR', 'Nonsense Mutation']:
            mutation_markers[category]['marker'] = 'fill'
            mutation_markers[category]['height'] = 0.5
            mutation_markers[category]['width'] = 0.5
            mutation_markers[category]['zindex'] = 1
            
        elif category in ['Splice/UTR', 'Inframe Indel']:
            mutation_markers[category]['marker'] = 'x'
            mutation_markers[category]['zindex'] = 2

    return mutation_markers


def build_oncoprint_df_2(mutations_dict, gene_set, pt_sample_ordering):

    mutation_types_dict = {'intron_variant': 'Silent/Intron/IGR',
                           'synonymous_variant': 'Silent/Intron/IGR',
                           'upstream_gene_variant': 'Silent/Intron/IGR',
                           'downstream_gene_variant': 'Silent/Intron/IGR',
                           'sequence_feature': 'Silent/Intron/IGR',
                           'protein_protein_contact': 'Silent/Intron/IGR',
                           'intergenic_region': 'Silent/Intron/IGR',
                            'splice_region_variant': 'Splice/UTR',
                            '5_prime_UTR_premature_start_codon_gain_variant': 'Splice/UTR',
                            '3_prime_UTR_variant': 'Splice/UTR',
                            '5_prime_UTR_variant': 'Splice/UTR',
                            'splice_region_variant&intron_variant': 'Splice/UTR',
                            'splice_acceptor_variant&intron_variant': 'Splice/UTR',
                            'splice_region_variant&synonymous_variant': 'Splice/UTR',
                            'splice_donor_variant&intron_variant': 'Splice/UTR',
                            'frameshift_variant': 'Frameshift Indel',
                            'frameshift_variant&splice_region_variant': 'Frameshift Indel',
                            'conservative_inframe_deletion': 'Inframe Indel',
                            'conservative_inframe_insertion': 'Inframe Indel',
                            'disruptive_inframe_deletion': 'Inframe Indel',
                            'disruptive_inframe_insertion': 'Inframe Indel',
                            'missense_variant': 'Missense Mutation',
                            'missense_variant&splice_region_variant': 'Missense Mutation',
                            'stop_gained': 'Nonsense Mutation',
                            'stop_lost': 'Nonsense Mutation',
                            'stop_gained&splice_region_variant': 'Nonsense Mutation'}
    
    mutation_codes_dict = {'Silent/Intron/IGR': 0,
                           'Splice/UTR': 1,
                           'Frameshift Indel': 2,
                           'Inframe Indel': 3,
                           'Missense Mutation': 4,
                           'Nonsense Mutation': 5,
                           'Multiple Mutations': 6}

    op_df = pd.DataFrame(-1, index=gene_set, columns=pt_sample_ordering)

    for gene in tqdm(gene_set):

        for pt in mutations_dict:

            for sample in mutations_dict[pt]:

                sample_gene_mutations_df = mutations_dict[pt][sample].loc[mutations_dict[pt][sample]['Gene_Name'] == gene]

                for _, row in sample_gene_mutations_df.iterrows():
                    
                    if row['Mutation_Type'] in mutation_types_dict:

                        mutation_type = mutation_types_dict[row['Mutation_Type']]

                        # FIX ME: Multiple mutations where at least one is a productive mutation 
                        if op_df.loc[gene, sample] == -1:
                            op_df.loc[gene, sample] = mutation_codes_dict[mutation_type]
                            
                        else:
                            op_df.loc[gene, sample] = mutation_codes_dict['Multiple Mutations']

                    else:
                        print(sample_gene_mutations_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Gene_Name', 'Mutation_Type']].sort_values(by=['CHROM', 'POS']).values)
                    
    return op_df


def calculate_tmb(op_df, exome_size=32000000):

    tmb_df = pd.DataFrame(0, index=op_df.columns, columns=['mut_count', 'tmb'])

    for sample in op_df.columns:

        tmb_df.loc[sample, 'mut_count'] += op_df[sample].value_counts().get(2, 0)
        tmb_df.loc[sample, 'mut_count'] += op_df[sample].value_counts().get(3, 0)
        tmb_df.loc[sample, 'mut_count'] += op_df[sample].value_counts().get(4, 0)
        tmb_df.loc[sample, 'mut_count'] += op_df[sample].value_counts().get(5, 0)
        tmb_df.loc[sample, 'mut_count'] += op_df[sample].value_counts().get(6, 0)

        tmb_df.loc[sample, 'tmb'] = tmb_df.loc[sample, 'mut_count'] / exome_size
     
    return tmb_df

# def create_oncoprint(op_df):


