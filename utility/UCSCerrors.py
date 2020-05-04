class Error(Exception):
    pass
class DataFileError(Error):

    error_type = "DataFileError"
    def __init__(self, code, message):
        self.code = code
        self.message = message

class SequenceDataError(Error):
    error_type = "SequenceDataError"
    def __init__(self, code, message):
        self.code = code
        self.message = message

class NCBIQueryError(Error):

    def __init__(self, code, message):
        self.code = code
        self.message = message

def write_errors(errors_fpath,gene_symbol,error):
    """Maintains a tsv/ DataFrame of gene symbols and errors generated during a run.
    """
    etype,ecode,emsg = error.error_type,error.code,error.message
    if os.path.exists(errors_fpath):
        check_ef,errors_df = load_errors(errors_fpath)
        if gene_symbol in errors_df['gene_symbol'].unique():
            gene_error_df = errors_df.loc[errors_df['gene_symbol']==gene_symbol,:]
            if gene_error_df['error_message'].str.contains(emsg).any() or \
                (etype == 'SequenceDataError' and gene_error_df['error_type'].str.contains(etype).any()):
                print_errors(errors_df,gene_symbol,message=emsg)
                return
    else:
        error_columns = ['gene_symbol','error_type','error_code','error_message']
        errors_df = pd.DataFrame(columns=error_columns)
    error_row = pd.Series({'gene_symbol':gene_symbol,'error_type':etype,'error_code':ecode,'error_message':emsg})
    errors_df = errors_df.append(error_row,ignore_index=True)
    print_errors(errors_df,gene_symbol,emsg)
    errors_df.to_csv(errors_fpath,sep='\t')

def print_errors(errors_df,gene_symbol,message=None,error_type=None):
    """

    :param errors_df: DataFrame of errors information
    :param gene_symbol: symbol string used to filter errors_df
    :param message: Optional, if provided only prints errors with text matching message
    :param error_type: Optional, if provided, only prints errors of specified error_type
    :return: N/A. Writes error information to output
    """
    symbol_df = errors_df.loc[errors_df['gene_symbol']==gene_symbol,:]
    error_rows = symbol_df
    if message and error_rows['error_message'].str.contains(message).any():
        error_rows = error_rows.loc[error_rows['error_message']==message,:]
    elif error_type and error_rows['error_type'].str.contains(error_type).any():
        error_rows = error_rows.loc[error_rows['error_type'] == error_type, :]
    for idx,error_row in error_rows.iterrows():
        genename, error_type, error_code, error_msg = error_row.values
        print("{0}\t{1}\t{2}\t{3}".format(genename, error_type, error_code, error_msg))

def load_errors(errors_fpath,error_type=""):
    """Loads errors_df from specified file path. Used to prevent repeat operations that generate errors (ie failed
    OrthoDB or NCBI queries, lack of available GeneCards alias data

    :param errors_fpath: file path to errors.tsv which logs errors raised in acquisition/ analysis process. If path
    doesn't exist, returns check_error_file as False and an empty DataFrame (which won't be accessed).
    :param error_type: if provided, will only check against logged errors of given type
    :return: check_error_file: boolean on whether errors_df exists
    :return: errors_df: if check_error, returns file loaded from
    """
    if os.path.exists(errors_fpath):
        errors_df = pd.read_csv(errors_fpath, delimiter='\t',index_col=0)
        if error_type:
            errors_df = errors_df.loc[errors_df['error_type'] == error_type, :]
        check_error_file = True
        return check_error_file, errors_df
    else:
        check_error_file = False
        return check_error_file, pd.DataFrame(columns=['gene_symbol','error_type','error_code','error_message'])