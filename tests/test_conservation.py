import unittest
from utility import NCBIfilter, fastaUtility as fautil, directoryUtility, UCSCerrors
import pandas as pd
import os
from IPython.display import display
from conservation import analysis_record_filter as ar_filt
from conservation import conservation_summary as ar_filt
from Bio import SeqIO

class ConservationTest(unittest.TestCase):

    def test_analysis_subset(self):
        #Tests proper analysis set filtering for AGS (and other NCBI species for which UCSC data isn't available)
        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        test_records_fpath = "tests/tmp/ENST00000002596.6_9999records.tsv"
        test_filtered_path = "tests/tmp/ENST00000002596_6_9999filtered.fasta"
        full_test_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
        # print(full_test_df.index)
        boreo_df, rest_df = fautil.partition_UCSC_by_clade(full_test_df,'boreoeutheria')
        ncbi_idx = full_test_df.loc[full_test_df["NCBI_taxid"]==9999].index
        for path in [test_records_fpath,test_filtered_path]:
            if os.path.exists(path):
                os.remove(path)
        self.assertFalse(os.path.exists(test_filtered_path))
        self.assertTrue(10181 in boreo_df["NCBI_taxid"].unique())

        test_analysis_taxids = {9999:'Urocitellus parryii',10181:'Heterocephalus glaber',29073:'Ursus maritimus',
                                9994:'Marmota marmota marmota'}
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_combined_fpath,test_records_fpath,
                                                        UCSC_analysis_subset=boreo_df.index,NCBI_record_subset=ncbi_idx,
                                                        filtered_outpath=test_filtered_path,taxid_dict=test_analysis_taxids)
        self.assertTrue(os.path.exists(filtered_outpath))
        self.assertTrue(29073 not in records_df["NCBI_taxid"].unique())
        self.assertTrue(10181 not in records_df["NCBI_taxid"].unique())
        self.assertTrue(9999 in records_df["NCBI_taxid"].unique())
        for tid in rest_df["NCBI_taxid"]:
            self.assertTrue(tid not in records_df["NCBI_taxid"].unique())
        fasta_ids = [fasta.id for fasta in SeqIO.parse(filtered_outpath,'fasta')]
        self.assertTrue(len(fasta_ids) == len(records_df))

    def test_analysis_subset2(self):

        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        test_records_fpath = "tests/tmp/ENST00000002596.6_10181records.tsv"
        test_filtered_path = "tests/tmp/ENST00000002596_6_10181filtered.fasta"
        full_test_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
        # print(full_test_df.index)
        boreo_df, rest_df = fautil.partition_UCSC_by_clade(full_test_df,'boreoeutheria')
        ncbi_idx = full_test_df.loc[full_test_df["NCBI_taxid"]==10181].index
        for path in [test_records_fpath,test_filtered_path]:
            if os.path.exists(path):
                os.remove(path)
        self.assertFalse(os.path.exists(test_filtered_path))
        self.assertTrue(10181 in boreo_df["NCBI_taxid"].unique())

        test_analysis_taxids = {9999:'Urocitellus parryii',10181:'Heterocephalus glaber',29073:'Ursus maritimus',
                                9994:'Marmota marmota marmota'}
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_combined_fpath,test_records_fpath,
                                                        UCSC_analysis_subset=boreo_df.index,NCBI_record_subset=ncbi_idx,
                                                        filtered_outpath=test_filtered_path,taxid_dict=test_analysis_taxids)
        self.assertTrue(os.path.exists(filtered_outpath))
        self.assertTrue(29073 not in records_df["NCBI_taxid"].unique())
        self.assertTrue(len(records_df.loc[records_df['NCBI_taxid']==10181,:]) == 1)
        for tid in rest_df["NCBI_taxid"]:
            self.assertTrue(tid not in records_df["NCBI_taxid"].unique())
        fasta_ids = [fasta.id for fasta in SeqIO.parse(filtered_outpath,'fasta')]
        self.assertTrue(len(fasta_ids) == len(records_df))

    def test_drop_redundant_UCSC(self):
        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        combined_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)

        boreo_df, rest_df = fautil.partition_UCSC_by_clade(combined_df, 'boreoeutheria')
        #Test AGS as NCBI index => drops 13LGS in redundant_filtered
        ncbi_idx = combined_df.loc[combined_df["NCBI_taxid"] == 9999].index
        single_NCBI_combined = fautil.load_UCSC_NCBI_df(test_combined_fpath,UCSC_subset=boreo_df.index,NCBI_subset=ncbi_idx)

        redundant_filtered = ar_filt.drop_redundant_UCSC_records(single_NCBI_combined,ncbi_idx)
        self.assertTrue(43179 in single_NCBI_combined['NCBI_taxid'].unique())
        self.assertFalse(43179 in redundant_filtered['NCBI_taxid'].unique())

        #Tests NMR as NCBI index => drops NMR UCSC data
        ncbi_idx = combined_df.loc[(combined_df["NCBI_taxid"] == 10181) & \
                                   (~combined_df.index.str.contains('ENST'))].index
        self.assertTrue(len(ncbi_idx)==1)
        single_NCBI_combined = fautil.load_UCSC_NCBI_df(test_combined_fpath, UCSC_subset=boreo_df.index,
                                                        NCBI_subset=ncbi_idx)
        redundant_filtered = ar_filt.drop_redundant_UCSC_records(single_NCBI_combined, ncbi_idx)
        self.assertTrue(len(single_NCBI_combined.loc[single_NCBI_combined['NCBI_taxid']==10181]) == 2)
        self.assertTrue(len(redundant_filtered.loc[redundant_filtered['NCBI_taxid'] == 10181]) == 1)
        self.assertTrue(ncbi_idx[0] in redundant_filtered.index)

        #Tests dropping UCSC records if NCBI index contains multiple records.
        ncbi_idx = pd.Index(['XP_026259682.1','XP_021115942.1'],name="record_id")
        double_NCBI_combined = fautil.load_UCSC_NCBI_df(test_combined_fpath, UCSC_subset=boreo_df.index,
                                                        NCBI_subset=ncbi_idx)
        redundant_filtered = ar_filt.drop_redundant_UCSC_records(double_NCBI_combined, ncbi_idx)
        self.assertTrue(len(double_NCBI_combined.loc[double_NCBI_combined['NCBI_taxid'] == 10181]) == 2)
        self.assertTrue(len(redundant_filtered.loc[redundant_filtered['NCBI_taxid'] == 10181]) == 1)
        self.assertTrue(43179 in double_NCBI_combined['NCBI_taxid'].unique())
        self.assertFalse(43179 in redundant_filtered['NCBI_taxid'].unique())
        for idx_val in ncbi_idx:
            self.assertTrue(idx_val in redundant_filtered.index)


if __name__ == "__main__":
    unittest.main()