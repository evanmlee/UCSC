import SSfasta
import pandas as pd
import numpy as np
from IPython.display import display
import directoryUtility
import os
import re
from collections import OrderedDict
import fastaUtility
from collections import Counter
from Bio.SubsMat.MatrixInfo import blosum62

def gen_blos_df():
    global aas, blosum62_bg
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-']
    bg_probs = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]#, 0.000001]
    blosum62_bgdict = dict(zip(aas,bg_probs))
    blosum62_bg = bg_probs
    blos_df = pd.DataFrame(index=aas[:-1],columns=aas[:-1])
    for pair in blosum62:
        val = blosum62[pair]
        first, second = pair[0],pair[1]
        if first in aas and second in aas:
            blos_df.loc[first,second] = val
            blos_df.loc[second,first] = val
    sim_matrix = blos_df.values
    return aas, blosum62_bg, blos_df, sim_matrix

def find_uniques(align_df, sub_freq_threshold, test_species_idx,display_uniques):
    #align
    uniques = pd.DataFrame(index=align_df.index)
    for pos in align_df:
        col = align_df[pos]
        sub_counts = col.value_counts()
        c_counts = sub_counts.value_counts()
        for index,value in sub_counts.iteritems():
    #         if(value == 1):# and c_counts.loc[value] == 1:
            if value <= sub_freq_threshold:
                if col[test_species_idx] == index and index != 'X':
                    uniques.loc[:,pos] = col
                    break
    if display_uniques:
        print("Threshold Number of Sequences: "+str(int(sub_freq_threshold)))#/len(ordered)))
        display(uniques)
    return uniques


def gap_ratio(col):
    return sum([c == '-' for c in col]) / len(col)


def sub_ratio(col, sub):
    return sum([c == sub for c in col]) / len(col)


def species_gaps_dict(align_df, gap_char='-'):
    gaps_dict = OrderedDict()
    unique_gaps = OrderedDict()
    unique_subs = OrderedDict()
    align_nd = align_df.values
    for j, spec in enumerate(align_df.index):
        gaps_dict[spec] = []
        unique_gaps[spec] = []
        unique_subs[spec] = []
        spec_row = align_df.loc[spec, :].values
        #         row = align_nd[j]
        for i, entry in enumerate(spec_row):
            #             outgroup_col_ = align_df.iloc[:,i].copy()
            outgroup_col_ = align_nd.T[i].copy()
            #             outgroup_col_.drop(spec,inplace=True)
            outgroup_col_ = np.delete(outgroup_col_, obj=j)
            ratio = gap_ratio(outgroup_col_)
            if entry == gap_char:
                gaps_dict[spec].append(i + 1)
                # Check if remainder of positions have gap; if none do, add to unique gaps.
                if ratio <= 0.1:
                    unique_gaps[spec].append(i + 1)
            else:
                if ratio >= 0.9:
                    unique_subs[spec].append(i + 1)

    return gaps_dict, unique_gaps, unique_subs


def seq_to_align_pos(seq_pos, species, gap_dict):
    spec_gaps = gap_dict[species]
    for gap in spec_gaps:
        if gap <= seq_pos:
            seq_pos += 1
    return seq_pos


def outgroup_seq_pos(gaps, test_pos, test_species, out_species):
    test_gaps = gaps[test_species]
    out_gaps = gaps[out_species]
    out_pos = test_pos
    for test_gap in test_gaps:
        if test_gap <= out_pos and test_gap not in out_gaps:
            out_pos += 1
    for out_gap in out_gaps:
        if out_gap <= out_pos and out_gap not in test_gaps:
            out_pos -= 1
            #             print("hurray")
    return out_pos


def force_noshift(align_df, test_species, subs):
    test_row = align_df.loc[test_species, :]
    for pos in subs:
        sub = subs[pos]
        if test_row.loc[pos] != sub:
            test_row.loc[pos] = sub


def convert_sub_str(sub_str):
    # Converts a substitution string (i.e. P39T) into a tuple of (position,outgroup variant,test variant)
    # for example 'P39T' -> (39,'P','T')
    try:
        groups = re.search("([A-Z])(\d+)([A-Z])", sub_str).groups()
        out, pos, test = groups[0], int(groups[1]), groups[2]
        return (pos, out, test)
    except AttributeError as e:
        #         raise UserWarning("oops")
        print(sub_str + " did not match expected pattern i.e. P39T, A79V")
        print("Will not attempt to realign for this input string")
        return None


def force_with_shifts(align_df, test_species, sub_str, gap_dict, unique_gaps):
    # subs is a dictionary from positions to tuples of (outgroup substitution, test variant)
    sub = convert_sub_str(sub_str)
    if sub:
        seq_pos, out, test = sub[0], sub[1], sub[2]
        # Check test is at seq_pos
        align_pos = seq_to_align_pos(seq_pos, test_species, gap_dict)
        test_variant = align_df.loc[test_species, align_pos]
        if test_variant != test:
            print("Input string: " + sub_str)
            print(
                "Could not find specified substitution at specified position. Alignment will not be changed Continuing...")
            pass
        else:
            # Idea: generate outgroup seq pos for each outgroup species individually. Check consensus of position, use mode to select for predicted position. Check substitution at that position. Threshold.
            outgroup = align_df.index.copy().drop(test_species)
            for out_species in outgroup:
                out_seq_pos = outgroup_seq_pos(gap_dict, seq_pos, test_species, out_species)
                out_align_pos = seq_to_align_pos(out_seq_pos, out_species, gap_dict)
    else:
        pass


def align_pos_to_native_pos(align_df, test_species, align_pos):
    pre_pos = align_df.loc[test_species, :align_pos - 1]
    pre_pos_vc = pre_pos.value_counts()
    native_pos = align_pos
    if '-' in pre_pos_vc:
        native_pos -= pre_pos_vc['-']
    return native_pos


def calculate_sequence_weights(msa, aas):
    """ Calculate the sequence weights using the Henikoff '94 method
    for the given msa. """
    n = len(msa)
    alignlen = len(msa[0])
    seq_weights = np.zeros((n))
    for col in msa.T:
        freqs = Counter(col)
        del freqs['-']
        num_observed = len(freqs)

        #         d_ = [freqs[aa]*num_observed if aa != '-' else 0 for aa in col]
        d_ = [freqs[aa] * num_observed for aa in col]
        inverse = [1 / d if d > 0 else 0 for d in d_]
        seq_weights += inverse
    seq_weights /= alignlen
    return seq_weights

def weighted_fcpseudo(col, seq_weights, pc_amount=0.000001):
    if len(seq_weights) != len(col):
        return [1]*len(col)
    else:
        freqs = np.array([pc_amount]*len(aas))
        for i,aa in enumerate(aas):
#             print(col)
#             print([seq_weights[j] if aa == entry else 0 for j,entry in enumerate(col)])
            freqs[i] += sum([seq_weights[j] if aa == entry else 0 for j,entry in enumerate(col)])
        for i,entry in enumerate(freqs):
#             if entry > pc_amount:
#                  freqs[i] -= pc_amount
            freqs[i] /= (sum(seq_weights) + len(aas)*pc_amount)
        return freqs

def weighted_gap_penalty(col,weights):
    #Return simple gap penalty multiplier for column
    gap_weights = np.array([w if c == "-" else 0 for c,w in zip(col,weights)])
    gap_sum = sum(gap_weights)
    return 1 - (gap_sum/sum(weights))
# print(weighted_gap_penalty(['A','A','-'],[0.3333,0.3333,0.3333]))
def relative_entropy(pC,bg,gap_penalty=0):
    """Calculate the relative entropy of the column distribution with a
    background distribution specified in bg_distr. This is similar to the
    approach proposed in Wang and Samudrala 06."""
    if len(bg) == 20 and len(pC) == 21:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    ratio = pC/bg
    log_ratio = np.log(ratio)/np.log(2)
    RE = sum([a*b for a,b in zip(pC,log_ratio)])
    if gap_penalty:
        RE*= gap_penalty
    return RE

def col_RE(col,bg,weights,gap_penalty=1):
    pC = weighted_fcpseudo(col,weights)
    if len(bg) == 20 and len(pC) == 21:
        pC = pC[:-1]
        pC /= sum(pC)
    ratio = pC/bg
    log_ratio = np.log(ratio)/np.log(2)
    RE = sum([a*b for a,b in zip(pC,log_ratio)])
    if gap_penalty:
        RE*= gap_penalty
    return RE
# relative_entropy(np.array([0.25,0.25,0.25,0.25]),np.array([0.1,0.1,0.1,0.7]))

def JSD(col,bg,weights,gap_penalty=1,jsd_lam=0.5):
    pC = weighted_fcpseudo(col,weights)
    if len(bg) == 20:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    gap_penalty = weighted_gap_penalty(col,weights)
    r_dist = np.array([jsd_lam*aa for aa in pC]) + np.array([(1-jsd_lam)*aa for aa in bg])

    RE_pCr = relative_entropy(pC,r_dist,gap_penalty)
    RE_qr = relative_entropy(bg,r_dist, gap_penalty)
    jsd = jsd_lam*RE_pCr + (1-jsd_lam)*RE_qr
    return jsd

#Generate a Jensen Shannon Difference Series for position indices within the MSA
def generate_jsd_series(test_species,align_df, alignlen, aas,blosum62_bg):
    jsd_srs = pd.Series(index=range(1,alignlen))
    jsd_lam = 0.5
    gap_threshold = 0.3
    gap_positions = []
    jsd_df = align_df.drop(test_species,axis=0)
    # jsd_df = align_df.copy()
    jsd_nd = jsd_df.values
    weights = calculate_sequence_weights(jsd_nd,aas)
    # print(weights)
    for i,col in enumerate(jsd_nd.T):
        gap_penalty = weighted_gap_penalty(col,weights)
        jsd = JSD(col,blosum62_bg,weights,gap_penalty)
        jsd_srs[i+1] = jsd
        if gap_ratio(col) > gap_threshold:
            gap_positions.append(i+1)
    # jsd_zscores = calc_z_scores(jsd_srs)
    with pd.option_context("display.max_rows",None,"display.max_columns", None):
#         display(jsd_srs)
        pass
    # return jsd_srs, jsd_zscores, gap_positions
    return jsd_srs, gap_positions

def summary_table2(align_df, pos_jsds,test_species,blos_df, display_summary):
    jsd_nd = align_df.copy().drop(index=test_species).values
    unique_positions = pos_jsds.index
    unique_positions.name = "MSA Position"
    n_seq = len(align_df.index)
    cols_list = ["Test Species Position", "Test Variant","Test Variant Instances","Outgroup Variant",\
                 "Outgroup Variant Instances","Aligned Sequences", "Gap Fraction", "JSD",\
                 "Test vs Outgroup Blosum62", "Outgroup Pairwise Blosum62"]
    summary_df = pd.DataFrame(index=unique_positions,columns=cols_list)
    for pos in unique_positions:
        native_pos = align_pos_to_native_pos(align_df, test_species,pos)
        jsd = pos_jsds[pos]
        align_col_srs = align_df.loc[:,pos]
        align_col = align_col_srs.values
        outgroup_variant = align_col_srs.mode().iloc[0]
        gap_fraction = gap_ratio(align_col)
        test_variant = align_df.loc[test_species,pos]
        vc = align_df[pos].value_counts()
        og_var_freq = vc[outgroup_variant]
        test_var_freq = vc[test_variant]
        col = jsd_nd.T[pos-1]
        if not test_variant == '-':
            blos_mean = np.mean(blos_df[test_variant][col])
        else:
            blos_mean = np.nan
        #pairwise calculations
        others_col = jsd_nd.T[pos-1]
        pairwise_scores = []
        skip_chars = ['-','X','Z']
        for i,first in enumerate(others_col):
            for j,second in enumerate(others_col):
                if i<j:
                    if not first in skip_chars and not second in skip_chars:
                        pairwise_scores.append(blos_df[first][second])
                    else:
                        pass
        if pairwise_scores:
#             blosum62_pairwise[pos] = np.mean(pairwise_scores)
            blos_pw = np.mean(pairwise_scores)
        else:
            blos_pw = np.nan
        row_values = [native_pos, test_variant, test_var_freq,outgroup_variant,og_var_freq,n_seq,gap_fraction,jsd,blos_mean,blos_pw]
        summary_row = pd.Series(name=pos,data=dict(zip(cols_list,row_values)))
#         summary_df = summary_df.append(summary_row)
        summary_df.loc[pos] = summary_row
    #Output conditional on config option
    if display_summary:
        print(test_species+" Semi-Unique Substitutions Summary")
        print("Number of Species: "+str(len(align_df.index)))
        with pd.option_context("display.max_rows",None,"display.max_columns", None,"display.max_colwidth",200):
            display(summary_df)
#             pass
    return summary_df

def write_output4(records_fpath,summary_fpath,records_df=pd.DataFrame(),summary_df=pd.DataFrame(),write_seq=False):
    """.fasta file output writing (both pre and post alignment) are now handled during
    process_input (specifically in seq_srs_to_align_df by calling construct_id_dm)"""
    if not records_df.empty:
        if not write_seq and "seq" in records_df:
            drop_seq = records_df.drop(columns="seq",inplace=False)
            drop_seq.to_csv(records_fpath,float_format='%.5f')
        else:
            records_df.to_csv(records_fpath,float_format='%.5f')
    if not summary_df.empty:
#         summary_path = "{0}/{1}_summary.csv".format(gene_dir_path,gene_name)
        summary_df.to_csv(summary_fpath,float_format='%.5f')


def write_SS_analysis(dir_vars,xref_table,taxid_dict):

    aas, blosum62_bg, blos_df, sim_matrix = gen_blos_df()
    overall_summary_cols = ["Gene", "MSA Position", "Test Species Position", "Test Variant", "Test Variant Instances", \
                            "Outgroup Variant", "Outgroup Variant Instances", \
                            "Aligned Sequences", "Gap Fraction", "JSD", \
                            "Test vs Outgroup Blosum62", "Outgroup Pairwise Blosum62"]
    updated_summary_ordered_cols = ["Test Species Position", "Test Variant", "Test Variant Instances", \
                                    "Outgroup Variant", "Outgroup Variant Instances", \
                                    "Aligned Sequences", "Gap Fraction", "JSD", \
                                    "Test vs Outgroup Blosum62", "Outgroup Pairwise Blosum62"]
    numeric_cols = ["Gap Fraction", "JSD", \
                    "Test vs Outgroup Blosum62", "Outgroup Pairwise Blosum62"]
    UNIQUE_THRESH = 0.05


    # dtype_specification = zip(numeric_cols,[np.float64]*len(numeric_cols))
    overall_summary = pd.DataFrame(columns=overall_summary_cols)

    config,taxid_dict = directoryUtility.config_initialization()
    dataset_config, taxonomy_config, directory_config = config["DATASET"], config["TAXONOMY"], config["DIRECTORY"]

    output_dir = "summary/{0}".format(dataset_identifier)
    directoryUtility.create_directory(output_dir)

    length_metrics_fpath = "{0}/length_metrics.tsv".format(combined_run_dir)
    lm_columns = ["NCBI_gid","euth_mean","euth_median","allncbi_mean","allncbi_median","bestncbi_mean","bestncbi_median"]
    if os.path.exists(length_metrics_fpath):
        field_types = [str]
        field_types.extend([float]*6)
        field_conv_dict = dict(zip(lm_columns,field_types))
        lm_df = pd.read_csv(length_metrics_fpath,sep='\t',index_col="UCSC_transcript_id",dtype=field_conv_dict)

    ts_id_list = lm_df.index

    for col in numeric_cols:
        overall_summary[col] = overall_summary[col].astype(np.float64)
    no_uniques = {"9999":[],"29073":[],"10181":[],"9994":[]}
    overall_summary = {}
    for taxid in taxid_dict:
        overall_summary[taxid] = pd.DataFrame(columns=overall_summary_cols)

    for ts_id in ts_id_list[:]:

        xref_row = xref_table.loc[ts_id,:]
        NCBI_gid = xref_row["NCBI_gid"]
        HGNC_symbol = xref_row["HGNC_symbol"]

        gene_output_dir = "{0}/{1}".format(output_dir,NCBI_gid)
        directoryUtility.create_directory(gene_output_dir)

        best_NCBI_fasta = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI_parent_dir,NCBI_gid)
        NCBI_fasta_df = fastaUtility.NCBI_fasta_df(best_NCBI_fasta,taxid_dict)

        MSA_fpath = "{0}/combined/{1}.fasta".format(bestNCBI_parent_dir,ts_id)
        # print(os.path.exists(MSA_fpath))
        MSA_srs = SSfasta.fasta_to_srs(MSA_fpath)
        MSA_df = SSfasta.align_srs_to_df(MSA_srs)
        # all_records_df = pd.read_csv(all_records_fpath, index_col="Unnamed: 0")

        NCBI_id_df = MSA_df.loc[MSA_df.index.str.contains('XP')]
        # display(NCBI_id_df)
        res_id_strs = ['hetGla2','speTri2','XP']
        res_pat = '|'.join(res_id_strs)
        res_index = MSA_df.index[MSA_df.index.str.contains(res_pat)]
        # print(res_index)
        outgroup_MSA_df = MSA_df.drop(index=res_index)
        n, align_len = len(outgroup_MSA_df.index), len(outgroup_MSA_df.columns)
        #         #Unique substitutions must be present in < UNIQUE_THRESH*100 % of sequences (or min of 1)
        sub_freq_threshold = max(int(n * UNIQUE_THRESH), 1)
        for spec_idx, record_row in NCBI_fasta_df.iterrows():

            align_spec_row = MSA_df.loc[spec_idx,:]
            NCBI_taxid = record_row["NCBI_taxid"]
            gene_tax_summary_fpath = "{0}/{1}_summary.tsv".format(gene_output_dir,NCBI_taxid)
            uniques = find_uniques(MSA_df, sub_freq_threshold, spec_idx, display_uniques=False)
            spec_no_uniques = no_uniques[NCBI_taxid]
            spec_overall_summary = overall_summary[NCBI_taxid]
            spec_overall_summary_fpath = "{0}/{1}_summary.tsv".format(output_dir,NCBI_taxid)
            if len(uniques.columns) == 0:
                #             write_output3(gene_name,run_name,final_records_df)
                spec_no_uniques.append(ts_id)
                no_uniques[NCBI_taxid] = spec_no_uniques
                continue
            else:
                JSD_df = outgroup_MSA_df.append(align_spec_row)
                align_nd = outgroup_MSA_df.values
                seq_weights = calculate_sequence_weights(align_nd, aas)
                # Unique substitution identification, metrics calculations and individual gene output writing
                gap_dict, unique_gaps, unique_subs = species_gaps_dict(JSD_df)
                jsd_srs, gap_positions = generate_jsd_series(spec_idx, JSD_df, align_len, aas, blosum62_bg)
                unique_JSDs = jsd_srs.loc[uniques.columns]
                # sw_metrics = calc_sw_metrics(jsd_srs, uniques, gap_positions, align_len, config)
                #             if window_calculations:
                #                 summary_df = summary_table(final_align_df,sw_metrics,test_species_id,blos_df,display_summary)
                #             else:
                summary_df = summary_table2(JSD_df, unique_JSDs, spec_idx, blos_df, display_summary=False)
                # write_output4(all_records_fpath, prev_summary_fpath, summary_df=summary_df)
                if not summary_df.empty:
                    #         summary_path = "{0}/{1}_summary.csv".format(gene_dir_path,gene_name)
                    summary_df.to_csv(gene_tax_summary_fpath, float_format='%.5f',sep='\t')
                # MSA_pos_idx = summary_df.index
                # AGS_var_srs = AGS_row.loc[MSA_pos_idx]
                # AGS_13LGS_consensus_srs = AGS_var_srs == summary_df["Test Variant"]
                # updated_summary_df = summary_df.copy()
                # updated_summary_df.loc[:, "AGS Variant"] = AGS_var_srs
                # updated_summary_df.loc[:, "AGS 13LGS Consensus"] = AGS_13LGS_consensus_srs

                # updated_summary_df = updated_summary_df[updated_summary_ordered_cols]
                # write_output4(all_records_fpath, updated_summary_fpath, summary_df=updated_summary_df)

                summary_df.reset_index(drop=False, inplace=True)
                if not "Gene" in summary_df.columns:
                    summary_df.insert(0, "Gene", [HGNC_symbol] * len(summary_df))
                spec_overall_summary = spec_overall_summary.append(summary_df, sort=False)
                overall_summary[NCBI_taxid] = spec_overall_summary
                spec_overall_summary.to_csv(spec_overall_summary_fpath, float_format='%.5f')

