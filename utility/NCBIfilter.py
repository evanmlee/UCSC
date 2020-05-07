import pandas as pd
import numpy as np
import importlib
from utility import directoryUtility
importlib.reload(directoryUtility)
from IPython.display import display
from utility.fastaUtility import UCSC_fasta_df, NCBI_fasta_df, filter_fasta_infile, min_dist_spec_record, construct_id_dm
import os
from query import NCBI_homology

def load_transcript_table(transcripts_fpath):
    transcript_table = pd.read_csv(transcripts_fpath,sep='\t',index_col="UCSC_transcript_id")
    return transcript_table

def select_NCBI_records(dir_vars,taxid_dict,UCSC_tid,NCBI_gid,selection="identity"):
    from IPython.display import display
    UCSC_parent,allNCBI,bestNCBI = [dir_vars[k] for k in ["UCSC_raw_parent","allNCBI_parent","bestNCBI_parent"]]

    raw_UCSC_fpath = "{0}/{1}.fasta".format(UCSC_parent,UCSC_tid)
    NCBI_all_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(allNCBI,NCBI_gid)
    combined_all_aln_fpath = "{0}/combined/{1}.fasta".format(allNCBI,UCSC_tid)

    UCSC_df = UCSC_fasta_df(raw_UCSC_fpath)
    allNCBI_df = NCBI_fasta_df(NCBI_all_aln_fpath,taxid_dict=taxid_dict)

    CLOSEST_EVO = {'9999':['mm10','rn6','speTri2'],'29073':['ailMel1','canFam3'],
                   '10181':['hetGla2','mm10','rn6'], '9994':['mm10','rn6','speTri2']}

    allseq_df = UCSC_df.append(allNCBI_df,sort=False)
    allseq_df.loc[allNCBI_df.index,"NCBI_taxid"] = allNCBI_df["NCBI_taxid"]
    id_dm, align_srs = construct_id_dm(allseq_df,combined_all_aln_fpath,aligned=True)
    bestNCBI_ids = []
    bestNCBI_dict = {}
    for taxid in taxid_dict:
        NCBI_tid_df = allNCBI_df.loc[allNCBI_df["NCBI_taxid"]==taxid,:]
        #Skip if no taxid records in allNCBI_df (ie no records for that species)
        if not NCBI_tid_df.empty:
            spec_record_ids = NCBI_tid_df.index
            closest_evos = CLOSEST_EVO[taxid]
            ce_pat = '|'.join(closest_evos)
            id_calc_records = UCSC_df.index[UCSC_df.index.str.contains(ce_pat)]
            #Check if records in id_calc_records are all empty strings
            if sum([len(seq) for seq in UCSC_df.loc[id_calc_records,"sequence"]]) == 0:
                #TODO: Figure out filtering if no UCSC records for comparison (likely only for polar bear NCBI)
                continue
            md_row, min_dist = min_dist_spec_record(id_dm,align_srs.index,spec_record_ids,
                                                                 id_calc_records,allseq_df)
            bestNCBI_ids.append(md_row.name)
            bestNCBI_dict[taxid] = md_row

    bestNCBI_df = allNCBI_df.loc[bestNCBI_ids]
    bestNCBI_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI,NCBI_gid)
    bestNCBI_all_fpath = "{0}/combined/{1}.fasta".format(bestNCBI,UCSC_tid)
    if True:#not os.path.exists(bestNCBI_aln_fpath): (os.path check now in outer function in UCSC_filter)
        filter_fasta_infile(bestNCBI_ids,NCBI_all_aln_fpath,bestNCBI_aln_fpath)
        final_align_ids = list(UCSC_df.index)
        final_align_ids.extend(bestNCBI_ids)
        # print(final_align_ids)
        filter_fasta_infile(final_align_ids,combined_all_aln_fpath,bestNCBI_all_fpath)
    #UCSC_euth: Records which taxonomically belong to euarchontoglires or laurasiatheria
    UCSC_euth_df = UCSC_df.iloc[:51, :]
    UCSC_euth_mean_len, UCSC_euth_median_len = UCSC_euth_df["length"].mean(), UCSC_euth_df["length"].median()
    allNCBI_mean_len, allNCBI_median_len = allNCBI_df["length"].mean(),allNCBI_df["length"].median()
    bestNCBI_mean_len, bestNCBI_median_len = bestNCBI_df["length"].mean(), bestNCBI_df["length"].median()
    length_labels = ["euth_mean","euth_median","allncbi_mean","allncbi_median","bestncbi_mean","bestncbi_median"]
    length_vals = [UCSC_euth_mean_len,UCSC_euth_median_len,
                   allNCBI_mean_len,allNCBI_median_len,
                   bestNCBI_mean_len,bestNCBI_median_len]
    length_metrics = dict(zip(length_labels,length_vals))
    return length_metrics

def filter_allNCBI_data(force_overwrite_tid_subset=[],alt_lm_fpath=""):
    # Config initialization, taxonomy dictionaries for NCBI species
    config, taxid_dict, dir_vars = directoryUtility.config_initialization()
    spec_name_dict = {taxid_dict[tid]: tid for tid in taxid_dict}
    # Directory name initialization
    dataset_config, taxonomy_config, directory_config = config["DATASET"], config["TAXONOMY"], config["DIRECTORY"]
    dataset_name, dataset_identifier = dir_vars["dataset_name"], dir_vars["dataset_identifier"]
    reorg_parent_dir, UCSC_raw, xref_xml, xref_summary = [dir_vars[k] for k in ["reorg_parent", "UCSC_raw_parent",
                                                                                "xref_xml", "xref_summary"]]
    combined_run, allNCBI, bestNCBI,orthologs_dir = [dir_vars[k] for k in ["combined_run",
                                                                "allNCBI_parent", "bestNCBI_parent","orthologs_parent"]]
    # Read stable ID table, stable ID xref table.
    transcripts_fpath = "{0}/ensembl_stable_ids_{1}.tsv".format(reorg_parent_dir, dataset_name)
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(xref_summary)
    transcript_table = load_transcript_table(transcripts_fpath)
    xref_table = NCBI_homology.load_NCBI_xref_table(xref_fpath,gid_dtype=int)
    orthologs_fpath = "{0}/orthologs_final.tsv".format(orthologs_dir)
    ortholog_id_table = NCBI_homology.load_orthologs_table(orthologs_fpath)
    NCBI_gid_table = xref_table.loc[:, ["stable_gene_id", "NCBI_gid", "HGNC_symbol"]]
    if alt_lm_fpath:
        length_metrics_fpath = alt_lm_fpath
    else:
        length_metrics_fpath = "{0}/length_metrics.tsv".format(combined_run)
    lm_columns = ["NCBI_gid", "euth_mean", "euth_median", "allncbi_mean", "allncbi_median", "bestncbi_mean",
                  "bestncbi_median"]
    if not os.path.exists(length_metrics_fpath):
        lm_df = pd.DataFrame(columns=lm_columns,dtype=float)
        lm_df = lm_df.astype(dtype={'NCBI_gid':str})
        lm_df.index.name = "UCSC_transcript_id"
    else:
        field_types = ['str']
        field_types.extend(['float'] * 6)
        field_conv_dict = dict(zip(lm_columns, field_types))
        lm_df = pd.read_csv(length_metrics_fpath, sep='\t', index_col="UCSC_transcript_id", dtype=field_conv_dict)

    non_empty_orthoIDs = ortholog_id_table.dropna(how='all')
    some_NCBI_IDs = NCBI_gid_table.loc[(NCBI_gid_table["NCBI_gid"].isin(non_empty_orthoIDs.index)), :]
    for UCSC_tid, row in some_NCBI_IDs.iterrows():
        NCBI_gid = row["NCBI_gid"]
        best_NCBI_fpath = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI, NCBI_gid)
        combined_best_tid_fpath = "{0}/combined/{1}.fasta".format(bestNCBI,UCSC_tid)
        # if (not os.path.exists(best_NCBI_fpath)) or UCSC_tid in force_overwrite_tid_subset:
        if (not os.path.exists(combined_best_tid_fpath)) or UCSC_tid in force_overwrite_tid_subset:
            lm_df.loc[UCSC_tid, "NCBI_gid"] = NCBI_gid
            try:
                length_metrics = select_NCBI_records(dir_vars, taxid_dict, UCSC_tid, NCBI_gid)
            except Exception as e:
                print("UCSC_tid: {0}".format(UCSC_tid))
                display("NCBI_gid: {0}".format(NCBI_gid))
                raise e

            lm_df.loc[UCSC_tid, lm_columns[1:]] = pd.Series(length_metrics)
            lm_df.to_csv(length_metrics_fpath, sep='\t', float_format='%.4f')

    #Reorder to index order from xref table
    lm_df = lm_df.loc[some_NCBI_IDs.index,:]
    lm_df.to_csv(length_metrics_fpath, sep='\t', float_format='%.4f')

def main():
    pd.options.display.max_columns = None
    from utility import UCSCerrors
    from utility.directoryUtility import config, taxid_dict, dir_vars
    ortholog_errors_fpath = "{0}/ortholog_errors.tsv".format(dir_vars['summary_run'])
    check_oe, ortholog_errors_df = UCSCerrors.load_errors(ortholog_errors_fpath)
    ambig_errors= ortholog_errors_df.loc[ortholog_errors_df['error_code']==3,:]
    alt_lm_fpath = "{0}/length_metrics_allncbi_rerun.tsv".format(dir_vars['combined_run'])
    filter_allNCBI_data()#force_overwrite_tid_subset=ambig_errors['tid'].unique())

if __name__ == '__main__':
    main()




