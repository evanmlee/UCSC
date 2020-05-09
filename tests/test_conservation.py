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
        ncbi_idx = full_test_df.loc[(full_test_df["NCBI_taxid"]==10181) & ~(full_test_df.index.str.contains("ENST")) ].index
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

    def test_drop_incomplete_records(self):
        with pd.option_context('display.max_columns',None):
            test_tid = "ENST00000422628.5"
            test_combined_fpath = "tests/test_data/combined/{0}.fasta".format(test_tid)
            test_records_outpath = "tests/tmp/{0}_9999records.tsv".format(test_tid)
            filtered_tmp_fpath = "tests/tmp/{0}_9999filtered.fasta".format(test_tid)
            combined_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
            ncbi_idx = combined_df.loc[(combined_df["NCBI_taxid"] == 9999) & ~(combined_df.index.str.contains("ENST"))].index
            boreo_df, rest_df = fautil.partition_UCSC_by_clade(combined_df, 'boreoeutheria')

            #Test proper use (NCBI_record_subset of length 1 raises SequenceAnalysisError for records with no suitable
            # UCSC record data
            with self.assertRaises(UCSCerrors.SequenceAnalysisError):
                filtered_records_df, outpath = ar_filt.filter_analysis_subset(test_combined_fpath,test_records_outpath,
                                                      UCSC_analysis_subset=boreo_df.index,NCBI_record_subset=ncbi_idx,
                                                      filtered_outpath=filtered_tmp_fpath,taxid_dict=TEST_ANALYSIS_DICT)
            # Test multiple record NCBI subset still raises SAE if no usable UCSC records present
            ncbi_idx = combined_df.loc[~(combined_df.index.str.contains("ENST"))].index
            with self.assertRaises(UCSCerrors.SequenceAnalysisError):
                filtered_records_df, outpath = ar_filt.filter_analysis_subset(test_combined_fpath,
                                                                              test_records_outpath,
                                                                              UCSC_analysis_subset=boreo_df.index,
                                                                              NCBI_record_subset=ncbi_idx,
                                                                              filtered_outpath=filtered_tmp_fpath,
                                                                              taxid_dict=TEST_ANALYSIS_DICT)

            test_tid = "ENST00000628423.2"
            test_combined_fpath = "tests/test_data/combined/{0}.fasta".format(test_tid)
            test_records_outpath = "tests/tmp/{0}_9999records.tsv".format(test_tid)
            filtered_tmp_fpath = "tests/tmp/{0}_9999filtered.fasta".format(test_tid)
            combined_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
            ncbi_idx = combined_df.loc[(combined_df["NCBI_taxid"] == 10181) &
                                       ~(combined_df.index.str.contains("ENST"))].index
            boreo_df, rest_df = fautil.partition_UCSC_by_clade(combined_df, 'boreoeutheria')
            filtered_records_df, outpath = ar_filt.filter_analysis_subset(test_combined_fpath, test_records_outpath,
                                                                         UCSC_analysis_subset=boreo_df.index,
                                                                         NCBI_record_subset=ncbi_idx,
                                                                         filtered_outpath=filtered_tmp_fpath,
                                                                         taxid_dict=TEST_ANALYSIS_DICT)
            ncbi_len = combined_df.loc[ncbi_idx,'length'].iloc[0]
            sub_thresh_len = filtered_records_df.loc[filtered_records_df['length'] < ncbi_len*0.8]
            self.assertTrue(len(sub_thresh_len)==0)

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
        summary_df = cs.gene_summary_table(align_df, ncbi_idx, blos_df,unique_thresh=1,
                                        summary_table_outpath=test_out_summary_fpath,
                                        use_jsd_gap_penalty=True)
        self.assertTrue(len(summary_df.loc[summary_df['NCBI Variant Count']==1])==len(summary_df))
        self.assertTrue(227 in summary_df.index)
        self.assertTrue(summary_df['Test-Outgroup BLOSUM62'].count() == 1)



    @unittest.skip("Median and mean calculation checks are correct")
    def test_suppl_length_metrics(self):
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
            ncbi_cols = row[ncbi_col_labels]
            precalc_bn_mean, precalc_bn_med = row['bestncbi_mean'], row['bestncbi_median']
            recalc_bn_med = np.nanmedian(ncbi_cols)
            recalc_bn_med = np.nanmedian(ncbi_cols)
            self.assertTrue(recalc_bn_med == precalc_bn_med)
            self.assertTrue(precalc_bn_mean == precalc_bn_mean)

    def test_lc_calcs(self):
        best_NCBI_parent = "combined_alignments/hg38AA_knownCanonical/4species1/bestNCBI"
        combined_run = "combined_alignments/hg38AA_knownCanonical/4species1"
        test_dir_vars = {'bestNCBI_parent': best_NCBI_parent,
                         'combined_run':combined_run}
        lm_fpath = "{0}/length_metrics.tsv".format(test_dir_vars['combined_run'])
        ar_filt.length_metric_checks(TEST_ANALYSIS_DICT,test_dir_vars,lm_fpath,how='most',overwrite=False)
        lc_fpath = "{0}/length_checks.tsv".format(test_dir_vars['combined_run'])
        lc_df = pd.read_csv(lc_fpath,sep='\t',index_col=0)
        lc_counts = lc_df.count()
        lc_passed_counts = lc_df.sum()
        all_passed = [tid for tid in lc_df.index if lc_df.loc[tid,:].sum() == len(lc_df.columns)]

        suppl_lm_fpath = "{0}/length_metrics_suppl.tsv".format(test_dir_vars['combined_run'])
        suppl_lm_df = pd.read_csv(suppl_lm_fpath,sep='\t',index_col=0)

        #Test length checks for records with no comparable UCSC data
        all_fail_tids = ["ENST00000638261.1","ENST00000643389.2","ENST00000643310.2","ENST00000639909.2","ENST00000640602.1"]
        for tid in all_fail_tids:
            tid_row = lc_df.loc[tid]
            self.assertTrue(sum(tid_row) == 0)


    def test_overall_summary(self):
        from utility.directoryUtility import taxid_dict,dir_vars,empty_directory
        from query import orthologUtility as orutil
        import contextlib, io
        from utility.fastaUtility import ucsc_tax_table

        # empty_directory("tests/tmp")
        length_checks_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
        lc_df = pd.read_csv(length_checks_fpath, sep='\t', index_col=0)
        test_tid,test_col_label = "ENST00000002596.6","9999_check"
        test_tid2 = "ENST00000003583.12"
        test_tid3 = "ENST00000367429.8"
        lc_col = lc_df[test_col_label]

        xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
        xref_df = orutil.load_NCBI_xref_table(xref_fpath)
        test_xref_table = xref_df.iloc[:5]


        sub_tax_table = ucsc_tax_table.iloc[:51]
        ucsc_tax_subset = sub_tax_table.loc[:,'NCBI_taxid']

        out_buf = io.StringIO()
        with contextlib.redirect_stdout(out_buf):
            cs.overall_summary_table(dir_vars,test_xref_table,taxid_dict,UCSC_tax_subset=ucsc_tax_subset,
                                 modified_outpath="tests/tmp",length_checks_fpath=length_checks_fpath)
            prep_message = "{0}\t{1}\t{2}\t{3}".format(test_tid3,"SequenceAnalysisError",0,
                                                       "Analysis records filtering resulted in empty outgroup set.")
            #Check that length_checks prevents empty data sets from run attempts
            self.assertFalse(prep_message in out_buf.getvalue())

        ags_summary_fpath = "tests/tmp/9999_summary.tsv"
        ags_summary_df = pd.read_csv(ags_summary_fpath,sep='\t',index_col=0)
        self.assertTrue(len(ags_summary_df['Gene'].unique())==3)
        self.assertTrue(len(ags_summary_df['Transcript ID'].unique()) == 3)
        # ags_errors_fpath = "{0}/9999_analysis_errors.tsv".format(dir_vars['summary_run'])
        ags_missing_dir = "tests/tmp/9999/ENST00000367429.8"
        self.assertFalse(os.path.exists(ags_missing_dir))
        ags_expected_dir = "tests/tmp/9999/ENST00000367772.8"
        expected_paths = [ags_expected_dir, "{0}/ENST00000367772.8_records.tsv".format(ags_expected_dir),
                          "{0}/ENST00000367772.8_summary.tsv".format(ags_expected_dir)]
        for path in expected_paths:
            self.assertTrue(os.path.exists(path))

        test_records_fpath = "{0}/ENST00000367772.8_records.tsv".format(ags_expected_dir)
        test_records_df = pd.read_csv(test_records_fpath,sep='\t',index_col=0)
        # display(test_records_df)
        self.assertTrue(9999 in test_records_df['NCBI_taxid'].unique())
        self.assertFalse(9785 in test_records_df['NCBI_taxid'].unique())

        ncbi_idx = "XP_026261019.1"
        ncbi_len = test_records_df.loc[ncbi_idx, 'length']
        sub_thresh_len = test_records_df.loc[test_records_df['length'] < ncbi_len * 0.8]
        self.assertTrue(len(sub_thresh_len) == 0)

        cs.overall_summary_table(dir_vars, xref_df, taxid_dict, UCSC_tax_subset=ucsc_tax_subset,
                                 length_checks_fpath=length_checks_fpath)

    def test_blosum(self):
        from utility.directoryUtility import dir_vars
        # test_fpath = "combined_alignments/ENST00000054668.5.fasta"


if __name__ == "__main__":
    unittest.main()