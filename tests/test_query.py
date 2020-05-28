import pandas as pd
import os
import unittest
from Bio import SeqIO
from IPython.display import display

from query import orthologUtility as orutil
from utility.directoryUtility import dir_vars, taxid_dict
from query import api_orthologs as ao

class NCBI_API_Test(unittest.TestCase):

    def test_mrna_record_query(self):

        pid_tsv_fpath = "{0}/NCBI_protein_IDs.tsv".format(dir_vars['orthologs_parent'])
        pid_df = pd.read_csv(pid_tsv_fpath,sep='\t',index_col=0)
        xrefs_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
        xref_df = orutil.load_NCBI_xref_table(xrefs_fpath,gid_dtype='int')
        test_row = xref_df.iloc[0,:]

        test_tid, test_gid = test_row.name, test_row['NCBI_gid']
        # display(pid_df)
        # str_test_gid = str(test_gid)
        # pid_row = pid_df.loc[str_test_gid,:]
        pid_row = pid_df.loc[test_gid,:]
        test_pid_str = pid_row["9999_pids"]
        print(test_pid_str)

        taxid = 9999
        best_record_fpath = "{0}/best_raw/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'],taxid,test_tid)
        self.assertTrue(os.path.exists(best_record_fpath))
        best_fastas = SeqIO.parse(best_record_fpath,'fasta')
        best_record_id = [fasta.id for fasta in best_fastas][0]
        # print(best_record_id)
        test_mrna_fpath = "{0}/{1}/{2}.fasta".format(dir_vars['orthologs_nuc'], taxid, test_tid)
        if not os.path.exists(test_mrna_fpath):
            ao.mrna_ortholog_request(dir_vars,test_tid,test_gid,taxid,best_record_id,test_pid_str)
        self.assertTrue(os.path.exists(test_mrna_fpath))

        mrna_fastas = SeqIO.parse(test_mrna_fpath,'fasta')
        for fasta in mrna_fastas:
            mrna_fasta_seq = fasta.seq

        best_fastas = SeqIO.parse(best_record_fpath, 'fasta')
        for fasta in best_fastas:
            ortholog_aa_seq = fasta.seq

        from Bio import Seq
        #Translate mrna sequence, ignore termination character
        mrna_trans = Seq.translate(mrna_fasta_seq)[:-1]
        self.assertTrue(mrna_trans == ortholog_aa_seq)
        # print(mrna_trans)
        # print(ortholog_aa_seq)
