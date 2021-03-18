import pandas as pd
import os
from utility import fastaUtility as fasta
from utility.directoryUtility import config, dir_vars, taxid_dict, ucsc_taxid_dict
from IPython.display import display
from utility import UCSCerrors as errors

from conservation import conservation_summary as cs

def read_transcript_length_info(tid):

    ucsc_raw_dir = "{0}/{1}/{2}".format(dir_vars["reorg_parent"],dir_vars["dataset_identifier"],dir_vars["clade_level"])
    best_ncbi_raw = "{0}/best_raw".format(dir_vars["bestNCBI_parent"])
    ucsc_fasta_fpath = "{0}/{1}.fasta".format(ucsc_raw_dir,tid)
    ucsc_df = fasta.load_UCSC_fasta_df(ucsc_fasta_fpath)

    lengths_row = {}
    for taxid in taxid_dict:
        label = "{0}_ncbi_length".format(taxid)
        ncbi_fasta_fpath = "{0}/{1}/{2}.fasta".format(best_ncbi_raw,taxid,tid)
        if os.path.exists(ncbi_fasta_fpath):
            ncbi_df = fasta.load_NCBI_fasta_df(ncbi_fasta_fpath,taxid_dict)
            ncbi_length = ncbi_df["length"].iloc[0]
            lengths_row[label] = ncbi_length
        else:
            #Force value to 0 to prevent float format in output .tsv file
            lengths_row[label] = 0
    for taxid in ucsc_taxid_dict:
        label = "{0}_ucsc_length".format(taxid)
        lengths_row[label] = ucsc_df.loc[ucsc_df["NCBI_taxid"]==taxid,"length"].iloc[0]
    for taxid in [9606,10090]:
        label = "{0}_ucsc_length".format(taxid)
        lengths_row[label] = ucsc_df.loc[ucsc_df["NCBI_taxid"] == taxid, "length"].iloc[0]
    return lengths_row

def transcript_length_table(xref_fpath):
    from query import orthologUtility as orutil
    xref_df = orutil.load_NCBI_xref_table(xref_fpath)
    ncbi_test_columns = ["{0}_ncbi_length".format(taxid) for taxid in taxid_dict]
    ucsc_test_columns = ["{0}_ucsc_length".format(taxid) for taxid in ucsc_taxid_dict]
    ucsc_og_columns = ["{0}_ucsc_length".format(taxid) for taxid in [9606,10090]]
    table_cols = ["Gene"]+ncbi_test_columns+ucsc_test_columns+ucsc_og_columns
    tl_df = pd.DataFrame(columns=table_cols)
    COUNT = 0
    for tid,row in xref_df.iterrows():
        lengths_row = read_transcript_length_info(tid)
        lengths_row = pd.Series(lengths_row,name=tid)
        tl_df = tl_df.append(lengths_row)
        tl_df.loc[tid, "Gene"] = row["HGNC_symbol"]
        COUNT+=1
    with pd.option_context("display.max_columns",None):
        display(tl_df)
    tl_outpath = "{0}/normalization_record_lengths.tsv".format(dir_vars["combined_run"])
    tl_df.to_csv(tl_outpath,sep='\t',float_format='%.12g')

def __append__gene_symbols_tl_table():
    tl_table_fpath = "{0}/normalization_record_lengths.tsv".format(dir_vars["combined_run"])
    tl_df = pd.read_csv(tl_table_fpath, sep='\t', index_col=0)
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars["xref_summary"])
    from query import orthologUtility as orutil
    xref_df = orutil.load_NCBI_xref_table(xref_fpath)
    symbol_srs = xref_df.loc[tl_df.index, "HGNC_symbol"]
    tl_df.insert(loc=0, column="gene", value=symbol_srs)
    display(tl_df)
    tl_df.to_csv(tl_table_fpath, sep='\t')

def cumulative_metric_table(taxid,xref_fpath,table_name="nongaps_suppl",source_db="NCBI",lc_fpath="",
                            metrics=["JSD","Norm. BLOSUM62","JSD BLOSUM Product"],rewrite_length_table=False,
                            rewrite_clm_table=False):
    from query import orthologUtility as orutil
    xref_df = orutil.load_NCBI_xref_table(xref_fpath)
    cm_outpath = "{0}/conservation/{1}/{2}_ncm_table.tsv".format(dir_vars["summary_run"], source_db, taxid)
    if os.path.exists(cm_outpath) and not rewrite_clm_table:
        print("Cumulative Metrics Table already exists for taxid: "+str(taxid)+
              ", specify 'rewrite_clm_table=True to recalculate")
        return
    if rewrite_length_table:
        transcript_length_table(xref_fpath)
    tl_outpath = "{0}/normalization_record_lengths.tsv".format(dir_vars["combined_run"])
    tl_df = pd.read_csv(tl_outpath,sep='\t',index_col=0)
    tl_col_label,hs_tl_col_label = "{0}_{1}_length".format(taxid,source_db.lower()),"9606_ucsc_length"
    tl_col, hs_tl_col = tl_df[tl_col_label], tl_df[hs_tl_col_label]

    cons_parent_dir = "{0}/conservation/{1}".format(dir_vars["summary_run"],source_db)
    summary_table_fpath = "{0}/{1}_{2}.tsv".format(cons_parent_dir,taxid,table_name)
    summary_df = pd.read_csv(summary_table_fpath,sep='\t',index_col=0)

    if source_db == "NCBI":
        analysis_errors_fpath = "{0}/{1}_analysis_errors.tsv".format(dir_vars["summary_run"],taxid)
    elif source_db == "UCSC":
        analysis_errors_fpath = "{0}/UCSC/{1}/{1}_analysis_errors.tsv".format(cons_parent_dir,taxid)
    check_aes, analysis_error_df = errors.load_errors(analysis_errors_fpath, error_type="SequenceAnalysisError")
    misc_cols = ["Gene","Test Record Length","Human Record Length","n Uniques"]
    cm_cols = ["Cumulative {0}".format(m) for m in metrics]
    lncm_cols= ["Length Normalized Cumulative {0}".format(m) for m in metrics]
    uncm_cols = ["Uniques Normalized Cumulative {0}".format(m) for m in metrics]
    cm_df = pd.DataFrame(index=xref_df.index,columns=misc_cols+cm_cols+lncm_cols+uncm_cols)
    # COUNT = 0
    sum_func = dict(zip(metrics,["sum"]*len(metrics)))
    for tid,row in xref_df.iterrows():
        # COUNT += 1
        # if COUNT > 20:
        #     break
        cm_df.loc[tid,"Gene"] = row["HGNC_symbol"]
        cm_df.loc[tid, "Human Record Length"] = hs_tl_col[tid]
        if tid not in summary_df["Transcript ID"].unique():
            if check_aes and tid in analysis_error_df["tid"].unique():
                cm_df.loc[tid,cm_cols] = [0]*len(cm_cols)
                cm_df.loc[tid, lncm_cols] = [0] * len(lncm_cols)
                cm_df.loc[tid, uncm_cols] = [0] * len(uncm_cols)
                cm_df.loc[tid, "Test Record Length"] = tl_col[tid]
                cm_df.loc[tid, "n Uniques"] = 0
            else:
                continue
        else:
            tid_sub_df = summary_df.loc[summary_df["Transcript ID"]==tid,:]
            sum_agg = tid_sub_df.agg(func=sum_func)
            cm_df.loc[tid,cm_cols] = sum_agg.values
            tlen = tl_col[tid]
            cm_df.loc[tid,lncm_cols] = sum_agg.values*100/tlen
            cm_df.loc[tid, uncm_cols] = sum_agg.values / len(tid_sub_df)
            cm_df.loc[tid, "Test Record Length"] = tl_col[tid]
            cm_df.loc[tid, "n Uniques"] = len(tid_sub_df)
    cm_df.to_csv(cm_outpath, sep='\t', float_format="%.6f")

def metric_ranked_gene_sets(taxid,source_db="NCBI",lc_filter=True):
    metrics = ["Length Normalized Cumulative JSD BLOSUM Product", "Uniques Normalized Cumulative JSD BLOSUM Product"]
    gs_prefixes = ["length_norm","uniques_norm"]
    cm_inpath = "{0}/conservation/{1}/{2}_ncm_table.tsv".format(dir_vars["summary_run"], source_db, taxid)
    gene_sets_dir = "{0}/conservation/gene_sets/{1}/".format(dir_vars["summary_run"],taxid)
    cm_df = pd.read_csv(cm_inpath,sep='\t',index_col=0)

    if lc_filter:
        if source_db == "NCBI":
            lc_fpath = "{0}/length_checks.tsv".format(dir_vars["combined_run"])
        else:
            lc_fpath = "{0}/ucsc_length_checks.tsv".format(dir_vars["combined_run"])
        lc_df = pd.read_csv(lc_fpath,sep='\t',index_col=0)
        lc_col_label = "{0}_check".format(taxid)
        passed_df = lc_df.loc[lc_df[lc_col_label] == True, :]
        anal_set = cm_df.loc[passed_df.index, :]
    else:
        anal_set = cm_df.dropna(axis=0, how='any')
    # print(len(anal_set))
    # print(len(anal_set["Gene"].unique()))
    display(anal_set)
    metrics_dict = dict(zip(metrics,gs_prefixes))
    for metric in metrics_dict:
        sorted = anal_set.sort_values(by=metric,ascending=False)
        gs = sorted["Gene"].unique()
        gs_fpath = "{0}/{1}_{2}_genes.txt".format(gene_sets_dir,taxid,metrics_dict[metric])
        # display(gs)
        # display(sorted)
        cs.write_gene_set_txt(gs,gs_fpath)



def main():
    # read_transcript_length_info("ENST00000367772.8")
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars["xref_summary"])
    write_transcript_lengths = False
    if write_transcript_lengths:
        transcript_length_table(xref_fpath)
    #CLM Table Calculation
    rewrite_clm_tables = False
    for taxid in taxid_dict:
        cumulative_metric_table(taxid,xref_fpath,rewrite_clm_table=rewrite_clm_tables)
        metric_ranked_gene_sets(taxid, source_db="NCBI")
    for taxid in ucsc_taxid_dict:
        cumulative_metric_table(taxid, xref_fpath,source_db="UCSC",rewrite_clm_table=rewrite_clm_tables)
        metric_ranked_gene_sets(taxid, source_db="UCSC")



if __name__ == "__main__":
    main()

