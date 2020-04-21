import pandas as pd
import numpy as np
import importlib
import directoryUtility
importlib.reload(directoryUtility)
from IPython.display import display
pd.options.display.max_columns = None
from fastaUtility import UCSC_fasta_df, NCBI_fasta_df, filter_fasta_infile, min_dist_spec_record, construct_id_dm
import os
import NCBI_homology

def load_transcript_table(transcripts_fpath):
    transcript_table = pd.read_csv(transcripts_fpath,sep='\t',index_col="UCSC_transcript_id")
    return transcript_table

def select_NCBI_records(dir_vars,tax_id_dict,UCSC_tid,NCBI_gid,selection="identity"):
    from IPython.display import display
    UCSC_parent,allNCBI_parent,bestNCBI_parent = [dir_vars[k] for k in ["UCSC_raw_parent","allNCBI_parent",
                                                                        "bestNCBI_parent"]]

    raw_UCSC_fpath = "{0}/{1}.fasta".format(UCSC_parent,UCSC_tid)
    NCBI_all_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(allNCBI_parent,NCBI_gid)
    combined_all_aln_fpath = "{0}/combined/{1}.fasta".format(allNCBI_parent,UCSC_tid)

    UCSC_df = UCSC_fasta_df(raw_UCSC_fpath)
    allNCBI_df = NCBI_fasta_df(NCBI_all_aln_fpath,tax_id_dict)

    CLOSEST_EVO = {'9999':['mm10','rn6','speTri2'],'29073':['ailMel1','canFam3'],
                   '10181':['hetGla2','mm10','rn6'], '9994':['mm10','rn6','speTri2']}

    allseq_df = UCSC_df.append(allNCBI_df,sort=False)
    allseq_df.loc[allNCBI_df.index,"NCBI_taxid"] = allNCBI_df["NCBI_taxid"]
    # display(allseq_df["sequence"])
    id_dm, align_srs = construct_id_dm(allseq_df,combined_all_aln_fpath,aligned=True)
    bestNCBI_ids = []
    bestNCBI_dict = {}
    for taxid in tax_id_dict:
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
    # display(bestNCBI_df)
    bestNCBI_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI_parent,NCBI_gid)
    bestNCBI_all_fpath = "{0}/combined/{1}.fasta".format(bestNCBI_parent,UCSC_tid)
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

def filter_allNCBI_data():
    # Config initialization, taxonomy dictionaries for NCBI species
    config, taxid_dict, dir_vars = directoryUtility.config_initialization()
    spec_name_dict = {taxid_dict[tid]: tid for tid in taxid_dict}
    # Directory name initialization
    dataset_config, taxonomy_config, directory_config = config["DATASET"], config["TAXONOMY"], config["DIRECTORY"]
    dataset_name, dataset_identifier = dir_vars["dataset_name"], dir_vars["dataset_identifier"]
    reorg_parent_dir, UCSC_raw, xref_xml, xref_summary = [dir_vars[k] for k in ["reorg_parent", "UCSC_raw_parent",
                                                                                "xref_xml", "xref_summary"]]
    orthologs_dir = dir_vars["orthologs_parent"]
    combined_run_dir, allNCBI_parent, bestNCBI_parent = [dir_vars[k] for k in
                                                         ["combined_run", "allNCBI_parent", "bestNCBI_parent"]]
    # Read stable ID table, stable ID xref table.
    transcripts_fpath = "{0}/ensembl_stable_ids_{1}.tsv".format(reorg_parent_dir, dataset_name)
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(xref_summary)
    transcript_table = load_transcript_table(transcripts_fpath)
    # xref_table = load_xrefs_table(xref_fpath)
    xref_table = NCBI_homology.load_NCBI_xref_table(xref_fpath,gid_dtype=int)
    orthologs_fpath = "NCBI_orthologs/orthologs_final.tsv"
    ortholog_id_table = NCBI_homology.load_orthologs_table(orthologs_fpath)
    NCBI_gid_table = xref_table.loc[:, ["stable_gene_id", "NCBI_gid", "HGNC_symbol"]]

    length_metrics_fpath = "{0}/length_metrics.tsv".format(combined_run_dir)
    lm_columns = ["NCBI_gid", "euth_mean", "euth_median", "allncbi_mean", "allncbi_median", "bestncbi_mean",
                  "bestncbi_median"]
    if not os.path.exists(length_metrics_fpath):
        field_types = [str]
        field_types.extend([float] * 6)
        field_conv_dict = dict(zip(lm_columns, field_types))
        lm_df = pd.DataFrame(columns=lm_columns, dtype=field_conv_dict)
        lm_df.index.name = "UCSC_transcript_id"
    else:
        field_types = [str]
        field_types.extend([float] * 6)
        field_conv_dict = dict(zip(lm_columns, field_types))
        lm_df = pd.read_csv(length_metrics_fpath, sep='\t', index_col="UCSC_transcript_id", dtype=field_conv_dict)

    non_empty_orthoIDs = ortholog_id_table.dropna(how='all')
    # print(len(ortholog_id_table))
    # print(len(non_empty_orthoIDs))
    # no_NCBI_IDs = NCBI_gid_table.loc[~(NCBI_gid_table["NCBI_gid"].isin(non_empty_orthoIDs.index)), :]
    some_NCBI_IDs = NCBI_gid_table.loc[(NCBI_gid_table["NCBI_gid"].isin(non_empty_orthoIDs.index)), :]

    for UCSC_tid, row in some_NCBI_IDs.iterrows():  # NCBI_gid_table.iterrows():
        NCBI_gid = row["NCBI_gid"]
        best_NCBI_fpath = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI_parent, NCBI_gid)
        if not os.path.exists(best_NCBI_fpath):
            # print(best_NCBI_fpath)
            lm_df.loc[UCSC_tid, "NCBI_gid"] = NCBI_gid
            try:
                length_metrics = select_NCBI_records(dir_vars, taxid_dict, UCSC_tid, NCBI_gid)
            except KeyError as ke:
                print("UCSC_tid: {0}".format(UCSC_tid))
                display("NCBI_gid: {0}".format(NCBI_gid))
                raise ke

            lm_df.loc[UCSC_tid, lm_columns[1:]] = pd.Series(length_metrics)
            lm_df.to_csv(length_metrics_fpath, sep='\t', float_format='%.4f')
    pass

def main():
    filter_allNCBI_data()

if __name__ == '__main__':
    main()