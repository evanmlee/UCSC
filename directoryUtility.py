import os
import shutil
import glob
import pandas as pd

def create_directory(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def remove_thing(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


def empty_directory(path):
    for i in glob.glob(os.path.join(path, '*')):
        remove_thing(i)


def write_errors(gene_name, message, errors_fpath="tmp/errors.tsv"):
    if not os.path.exists(errors_fpath):
        errors_f = open(errors_fpath, 'wt')
        errors_f.write("gene\tmessage\n")
        errors_f.close()
    f_line = "{0}\t{1}\n".format(gene_name, message)
    # Check if f_line already in errors_f
    errors_f = open(errors_fpath, 'rt')
    errors_flines = errors_f.readlines()
    if f_line not in errors_flines:
        errors_f = open(errors_fpath, 'at')
        errors_f.write(f_line)
        errors_f.close()

# Funcitons for reading config files, including run parameters, gene symbols list, and species lists.
def parse_config(config_file="config/config.txt"):
    """Parse config text file (INI format) to establish paramters for the run

    config_file: path to the config file ("config/config.txt" by default)
    """
    import configparser

    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def create_run_directories(dir_vars,taxid_dict):
    #Creates all directories from directory variables as well as intermediate directories
    for dir_label in dir_vars:
        dir_var = dir_vars[dir_label]
        create_directory(dir_var)
    #tmp file storage (distance matrix calculations, etc)
    create_directory("tmp")
    #creates defined directory structure for orthologs, allNCBI, bestNCBI parent directories
    orthologs_children = [k.strip() for k in taxid_dict.keys()]
    for child in orthologs_children:
        childpath = "{0}/{1}".format(dir_vars["orthologs_parent"],child)
        create_directory(childpath)
    allNCBI_children = ["combined","NCBI_alignments","NCBI_raw"]
    for child in allNCBI_children:
        childpath = "{0}/{1}".format(dir_vars["allNCBI_parent"], child)
        create_directory(childpath)
    bestNCBI_children = ["combined","NCBI_alignments"]
    for child in bestNCBI_children:
        childpath = "{0}/{1}".format(dir_vars["bestNCBI_parent"], child)
        create_directory(childpath)


def directory_variables(config):
    """Generates appropriately named directory paths for run based on parameters in config.

    Briefly:
    "reorg_parent": reorganized data parent directory. Contains split UCSC raw data, Ensembl xref xml data
    and summary tables
    "dataset": knownCanonical, knownGene, or refSeq. Specified in config
    "dataset_identifier": full dataset name (corresponds to UCSC raw Conservation track file names)
    "orthologs_parent": Parent directory for NCBI ortholog sequence data. Contains child dirs named by species taxonom ID
    "xref_xml": Raw ensembl xref xml data from REST API (subdir of reorganized_data)
    "xref_summary": Contains summary information for transcripts/ external refernces within greater UCSC dataset
    (subdir of reorganized_data)
    "combined_run": combined_alignments run-specific subdirectory to store allNCBI/ bestNCBI alignments
    "UCSC_raw_parent": parent directory for split UCSC raw (unsupplemented) transcript alignments
    "allNCBI_parent": parent directory for alignment outputs for 1) combined UCSC, NCBI and 2) just NCBI
    for all NCBI records per species
    "bestNCBI_parent": parent directory for alignment outputs for 1) combined UCSC, NCBI and 2) just NCBI
    for best NCBI record per species
    "summary_run": output parent directory for analysis outputs (subdirectory of UCSC/summary)

    :param config: configparser object from config/config.txt
    :return: dir_vars: dictionary of directory variable labels to directory variable paths
    """
    dataset_config, taxonomy_config, directory_config = config["DATASET"], config["TAXONOMY"], config["DIRECTORY"]
    hg_version,seq_format,dataset_name = [dataset_config[k] for k in dataset_config.keys()]
    reorg_parent = directory_config["ReorganizedDataDir"]
    dataset_identifier = "{0}{1}_{2}".format(hg_version, seq_format, dataset_name)

    combined_parent = directory_config["CombinedAlignmentDir"]
    combined_run_name = directory_config["CombinedRunSubDir"]

    orthologs_parent = "NCBI_orthologs"
    combined_run = "{0}/{1}/{2}".format(combined_parent, dataset_identifier, combined_run_name)
    UCSC_raw_parent = "{0}/{1}".format(reorg_parent, dataset_identifier)
    xref_summary = "{0}/xref_summary_{1}".format(reorg_parent, dataset_name)
    xref_xml = "{0}/xref_xml".format(reorg_parent)
    allNCBI_parent = "{0}/allNCBI".format(combined_run)
    bestNCBI_parent = "{0}/bestNCBI".format(combined_run)
    summary_run = "summary/{0}".format(dataset_identifier)

    dir_variable_labels = ["reorg_parent","dataset_name","dataset_identifier","orthologs_parent","xref_xml","xref_summary",
                          "combined_run","UCSC_raw_parent","allNCBI_parent","bestNCBI_parent","summary_run"]
    dir_values = [reorg_parent,dataset_name,dataset_identifier,orthologs_parent,xref_xml,xref_summary,
                  combined_run,UCSC_raw_parent,allNCBI_parent,bestNCBI_parent,summary_run]
    dir_vars = dict(zip(dir_variable_labels,dir_values))
    return dir_vars

def config_initialization():
    """Reads config file, generates directories based on run parameters in config

    :return: config: configparser object created from text file at config/config.txt
    :return: taxid_dict: dictionary from NCBI tax id (as string) to species name (with capitalization as standard)
    :return: dir_vars: dictionary from directory path labels to paths
    """
    config = parse_config()
    dataset_config, taxonomy_config, directory_config = config["DATASET"], config["TAXONOMY"], config["DIRECTORY"]
    taxids = list(taxonomy_config.keys())
    spec_names = [taxonomy_config[k].strip() for k in taxonomy_config.keys()]
    taxid_dict = dict(zip(taxids,spec_names))
    dir_vars = directory_variables(config)
    create_run_directories(dir_vars,taxid_dict)
    return config, taxid_dict, dir_vars


def load_ortholog_table(ortholog_fpath):
    import pandas as pd
    ortholog_table = pd.read_csv(ortholog_fpath,sep='\t',index_col="NCBI_gid",dtype='str')
    return ortholog_table

def ortholog_file_cleanup(ortholog_fpath):
    from IPython.display import display
    ortholog_table = ortholog_table = pd.read_csv(ortholog_fpath,sep='\t',index_col=0,dtype=str)
    # display(ortholog_table)
    print("Ortholog entries: {0}".format(len(ortholog_table)))
    for col in ortholog_table:
        col_dropna = ortholog_table[col].dropna()
        valid_idx = col_dropna.str.contains("^10|^11")
        col_valid = col_dropna.loc[valid_idx]
        display(col_valid)
        col_sus = col_dropna.loc[~(valid_idx)]
        display(col_sus)
        # suspect_vals = ortholog_table.loc[valid_idx,col]
        # print(len(suspect_vals))
        # suspect_vals = suspect_vals.dropna()
        # display(suspect_vals)

# import importlib
# import UCSCfilter
# importlib.reload(UCSCfilter)
# ortholog_fpath = "NCBI_orthologs/NCBI_orthologs.tsv"
# print(os.path.exists(ortholog_fpath))
# ortholog_file_cleanup(ortholog_fpath)


