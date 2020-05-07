import pandas as pd
import os
from utility.UCSCerrors import DataFileError

def SQL_transcript_table(dir_vars,tid_index,write_table=False):
    """Returns and optionally writes to file an int-indexed table of UCSC transcript IDs, and corresponding Ensembl
    transcript stable IDs and versions.

    :param (dict) dir_vars: contains path variables for run directories
    :param (pandas Index) tid_index: Index of UCSC_transcript_ids (ENST############.##)
    :param (boolean) write_table: If true, write to sql_transcript_fpath
    :return: DataFrame containing int-indexed UCSC_transcript IDs as well as stable transcript ID and version numbers
    which will be used to query Ensembl SQL tables.
    """
    reorg_parent,dataset_name = dir_vars["reorg_parent"],dir_vars["dataset_name"]
    sql_transcript_fpath = "{0}/{1}_transcript_table.tsv".format(reorg_parent, dataset_name)
    sql_transcript_table_index = tid_index.to_series()
    sql_transcript_table = pd.DataFrame(index=sql_transcript_table_index,
                                        columns=['stable_transcript_id', 'version'], dtype='str')
    sql_transcript_table.loc[:, 'stable_transcript_id'] = sql_transcript_table_index.str.extract("(ENST\d+)\.\d+")
    sql_transcript_table.loc[:, 'version'] = sql_transcript_table_index.str.extract("ENST\d+\.(\d+)")
    int_indexed = sql_transcript_table.reset_index(drop=False)
    int_indexed.index.name = "UCSC_transcript_index"
    if write_table:
        int_indexed.to_csv(sql_transcript_fpath,sep='\t',header=False)
    return int_indexed

def format_orthologs_sql(dir_vars,taxid_dict,formatted_outpath):
    """Reformats NCBI ortholog gid information from the SQL database tables to store appropriate orthologs in separate
    columns for each species.

    :param dir_vars: dictionary of directory paths based on run parameters, see directoryUtility
    :param taxid_dict: dictionary mapping NCBI taxids of species for which NCBI data will be fetched to species name
    :return: formatted_df: DataFrame indexed on NCBI human gid, columns named to [taxid]_gid and populated with
            NCBI gids of orthologs from appropriate species, based on data present in orthologs_sql_fpath.
    :return: absent_taxids: list of taxids that are in taxid_dict but not present in the ortholog_sql table
    """
    from IPython.display import display
    orthologs_parent,dataset_name = dir_vars["orthologs_parent"],dir_vars["dataset_name"]
    orthologs_sql_fpath = "{0}/{1}_orthologs_sql.tsv".format(orthologs_parent,dataset_name)
    if not os.path.exists(orthologs_sql_fpath):
        raise DataFileError(0,"Missing run orthologs file at path: {0}".format(orthologs_sql_fpath))
    col_names = ['ortholog_id','tax_id','GeneID','relationship','Other_tax_id','Other_GeneID']
    orthologs_sql_df = pd.read_csv(orthologs_sql_fpath,sep='\t',dtype=str,index_col=0,names=col_names)
    NCBI_gid_index = pd.Index(orthologs_sql_df["GeneID"].drop_duplicates())
    formatted_columns = ["{0}_gid".format(taxid) for taxid in taxid_dict]
    formatted_df = pd.DataFrame(index=NCBI_gid_index,columns=formatted_columns)
    for _,row in orthologs_sql_df.iterrows():
        hgid = row["GeneID"]
        spec_tax_id = row["Other_tax_id"]
        spec_gid = row["Other_GeneID"]
        spec_col = "{0}_gid".format(spec_tax_id)
        formatted_df.loc[hgid,spec_col] = spec_gid
    formatted_df.index.name = "NCBI_hgid"
    formatted_df.to_csv(formatted_outpath, sep='\t')

    #Check for species from taxid_dict not represented in sql table
    other_tid_unique = orthologs_sql_df["Other_tax_id"].unique()
    absent_taxids = [taxid for taxid in taxid_dict if taxid not in other_tid_unique]
    return formatted_df, absent_taxids