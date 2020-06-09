import unittest
from utility import NCBIfilter, fastaUtility, directoryUtility, UCSCerrors
import pandas as pd
import os
from IPython.display import display

class FastaDataFrameTest(unittest.TestCase):

    def test_load_combined_df(self):
        config,taxid_dict,dir_vars = directoryUtility.config_initialization()
        test_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        df = fastaUtility.load_UCSC_NCBI_df(test_fpath,taxid_dict)
        # UCSC_euth: Records which taxonomically belong to euarchontoglires or laurasiatheria
        euth_tax_df = df.iloc[:51, :]
        euth_tax_df = euth_tax_df.loc[:,['species_name','NCBI_taxid']]
        noneuth_tax_df = df.iloc[51:100,:]
        noneuth_tax_df = noneuth_tax_df.loc[:,['species_name','NCBI_taxid']]
        all_ncbi_tids,all_ncbi_names = list(taxid_dict.keys()),list(taxid_dict.values())
        tid_sets = [all_ncbi_tids,euth_tax_df['NCBI_taxid'],noneuth_tax_df['NCBI_taxid']]
        for tid_set in tid_sets:
            for tid in tid_set:
                self.assertTrue(tid in df['NCBI_taxid'].unique())
        #Test UCSC subset
        df = fastaUtility.load_UCSC_NCBI_df(test_fpath,taxid_dict,UCSC_subset=euth_tax_df.index)
        self.assertTrue(len(df) == len(euth_tax_df)+len(all_ncbi_tids))
        for tid in noneuth_tax_df:
            self.assertFalse(tid in df['NCBI_taxid'].unique())
        #Test NCBI subset
        test_ncbi_subset = df.loc[df['NCBI_taxid']==9999,:].index
        df = fastaUtility.load_UCSC_NCBI_df(test_fpath,NCBI_subset=test_ncbi_subset)
        self.assertTrue(len(df)==101)
        self.assertFalse('XP_015337450.1' in df.index)

        display(df)

        #Test taxid_dict doesn't correspond to table file in config, should cause filtering even if ncbi_subset not provided
        test_tid_dict = {9999:'Urocitellus parryii'}
        df = fastaUtility.load_UCSC_NCBI_df(test_fpath,ncbi_taxid_dict=test_tid_dict)

    def test_UCSC_df(self):
        test_path = "reorganized_data/hg38AA_knownCanonical/all/ENST00000002596.6.fasta"
        df = fastaUtility.load_UCSC_fasta_df(test_path)
        self.assertTrue(df.index.name == 'record_id')
        self.assertTrue(len(df)==100)
        self.assertTrue(len(df['NCBI_taxid'].dropna()) == 100)

    def test_tax_partitions(self):
        ucsc_test_path = "reorganized_data/hg38AA_knownCanonical/all/ENST00000002596.6.fasta"
        df = fastaUtility.load_UCSC_fasta_df(ucsc_test_path)
        bor_df,rest_df = fastaUtility.partition_UCSC_by_clade(df,'boreoeutheria')
        self.assertTrue("Condylura cristata" in bor_df["species_name"].unique())
        self.assertFalse("Condylura cristata" in rest_df["species_name"].unique())
        self.assertTrue("Loxodonta africana" in rest_df["species_name"].unique())
        self.assertFalse("Loxodonta africana" in bor_df["species_name"].unique())
        self.assertTrue(len(bor_df) == 51)
        self.assertTrue(len(rest_df) == 49)

        ucsc_tax_table = fastaUtility.ucsc_tax_table
        # with pd.option_context('display.max_rows',None):
        #     display(ucsc_tax_table.loc[:,['common_name','tax_name','NCBI_taxid']])

    def test_poorly_formatted_orthologs(self):
        from Bio import SeqIO
        test_fpath = "tests/test_data/NCBI_orthologs/2350.fasta"
        fastas = SeqIO.parse(test_fpath,'fasta')
        fasta_ids = [fasta.id for fasta in fastas]
        self.assertTrue(len(fasta_ids)==0)

    def test_parse_NCBI_df_unformatted_fasta(self):
        from utility.directoryUtility import taxid_dict
        test_fpath = "tests/test_data/best_raw/ENST00000236850.4.fasta"
        ncbi_df = fastaUtility.load_NCBI_fasta_df(test_fpath,taxid_dict)
        with pd.option_context('display.max_columns',None):
            display(ncbi_df)


if __name__ == "__main__":
    unittest.main()