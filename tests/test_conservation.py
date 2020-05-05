import unittest
from utility import NCBIfilter, fastaUtility as fautil, directoryUtility, UCSCerrors
import pandas as pd
import os
from IPython.display import display
from conservation import analysis_record_filter as ar_filt
from conservation import conservation_summary as cs
from Bio import SeqIO
from conservation.analysis_calc import blos_df

TEST_ANALYSIS_DICT = {9999: 'Urocitellus parryii', 10181: 'Heterocephalus glaber', 29073: 'Ursus maritimus',
                        9994: 'Marmota marmota marmota'}

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


        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_combined_fpath,test_records_fpath,
                                                        UCSC_analysis_subset=boreo_df.index,NCBI_record_subset=ncbi_idx,
                                                        filtered_outpath=test_filtered_path,taxid_dict=TEST_ANALYSIS_DICT)
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
                                                        filtered_outpath=test_filtered_path,taxid_dict=TEST_ANALYSIS_DICT)
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

    def test_gene_summary(self):

        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        test_out_records_fpath = "tests/tmp/ENST00000002596.6_9999records_gstest.tsv"
        test_out_summary_fpath = "tests/tmp/ENST00000002596.6_9999summary_gstest.tsv"
        combined_fasta_df = fautil.load_UCSC_NCBI_df(test_combined_fpath, TEST_ANALYSIS_DICT)
        ncbi_taxid = 9999
        ncbi_partition = (combined_fasta_df['NCBI_taxid'] == ncbi_taxid) & \
                         ~(combined_fasta_df.index.str.contains("ENST"))
        ncbi_idx = combined_fasta_df.loc[ncbi_partition, :].index

        boreo_df, rest_df = fautil.partition_UCSC_by_clade(combined_fasta_df, 'boreoeutheria')
        UCSC_test_subset = boreo_df.index
        records_df, filtered_aln_path = ar_filt.filter_analysis_subset(test_combined_fpath, test_out_records_fpath,
                                                                       UCSC_test_subset,
                                                                       NCBI_record_subset=ncbi_idx,
                                                                       taxid_dict=TEST_ANALYSIS_DICT, drop_redundant=True,
                                                                       drop_ncbi_from_ucsc=True)
        align_df = fautil.align_fasta_to_df(filtered_aln_path,ucsc_flag=True)
        summary_df = cs.gene_summary_table(align_df, ncbi_idx, test_idx=ncbi_idx, blos_df=blos_df,drop_NCBI=False,
                                        summary_table_outpath=test_out_summary_fpath,
                                        use_jsd_gap_penalty=True)


    def test_length_metrics_filter(self):
        import io
        import contextlib
        import numpy as np
        best_NCBI_parent = "combined_alignments/hg38AA_knownCanonical/4species1/bestNCBI"
        test_dir_vars = {'bestNCBI_parent':best_NCBI_parent}
        lm_fpath = "combined_alignments/hg38AA_knownCanonical/4species1/length_metrics.tsv"
        suppl_fpath = "combined_alignments/hg38AA_knownCanonical/4species1/length_metrics_suppl.tsv"
        out_buf = io.StringIO()
        with contextlib.redirect_stdout(out_buf):
            force_overwrite = False
            skip_ovewrwrite = os.path.exists(suppl_fpath) and not force_overwrite
            suppl_lm_fpath = ar_filt.length_metric_supplement(TEST_ANALYSIS_DICT,test_dir_vars,lm_fpath,
                                                          force_overwrite=force_overwrite,suppl_lm_df_fpath=suppl_fpath)
            if skip_ovewrwrite:
                check_output_msg = "Run with force_overwrite=True if recalculation desired."
                self.assertIn(check_output_msg,out_buf.getvalue())
                print("Successfully pre-checked file")
            else:
                self.assertTrue(os.path.exists(suppl_lm_fpath))

        suppl_lm_df = pd.read_csv(suppl_lm_fpath, sep='\t', index_col=0)
        ncbi_col_labels = ['9999_length', '10181_length', '29073_length', '9994_length']

        row_gen = suppl_lm_df.iterrows()
        for idx, row in row_gen:
            try:
                ncbi_cols = row[ncbi_col_labels]
                precalc_bn_mean, precalc_bn_med = row['bestncbi_mean'], row['bestncbi_median']
                recalc_bn_med = np.nanmedian(ncbi_cols)
                recalc_bn_med = np.nanmedian(ncbi_cols)
                self.assertTrue(recalc_bn_med == precalc_bn_med)
                self.assertTrue(precalc_bn_mean == precalc_bn_mean)
            except AssertionError as ae:
                display(row)
                display(ncbi_cols)
                raise ae

    def test_lm_filter(self):
        best_NCBI_parent = "combined_alignments/hg38AA_knownCanonical/4species1/bestNCBI"
        test_dir_vars = {'bestNCBI_parent': best_NCBI_parent}
        lm_fpath = "combined_alignments/hg38AA_knownCanonical/4species1/length_metrics.tsv"
        ar_filt.length_metric_filter(TEST_ANALYSIS_DICT,test_dir_vars,lm_fpath)

    def test_overall_summary(self):
        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"


if __name__ == "__main__":
    unittest.main()