import unittest
from utility import UCSCerrors
import pandas as pd
from IPython.display import display
import os

from utility.directoryUtility import ucsc_taxid_dict
from utility import fastaUtility as fautil
from conservation import analysis_record_filter as ar_filt
from Bio import SeqIO

class analysisRecordFilterTest(unittest.TestCase):

    def test_analysis_subset(self):
        # Tests proper analysis set filtering for AGS (and other NCBI species for which UCSC data isn't available)
        test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"
        test_records_fpath = "tests/tmp/ENST00000002596.6_9999records.tsv"
        test_filtered_path = "tests/tmp/ENST00000002596_6_9999filtered.fasta"
        full_test_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)

        boreo_df, rest_df = fautil.partition_UCSC_by_clade(full_test_df ,'boreoeutheria')
        ncbi_idx = full_test_df.loc[full_test_df["NCBI_taxid" ] ==9999].index
        for path in [test_records_fpath ,test_filtered_path]:
            if os.path.exists(path):
                os.remove(path)
        self.assertFalse(os.path.exists(test_filtered_path))
        self.assertTrue(10181 in boreo_df["NCBI_taxid"].unique())


        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_combined_fpath ,test_records_fpath ,9999,
                                                                      test_source="NCBI",
                                                                      UCSC_analysis_subset=boreo_df.index
                                                                      ,NCBI_record_subset=ncbi_idx,
                                                                      filtered_outpath=test_filtered_path)
        self.assertTrue(os.path.exists(filtered_outpath))
        self.assertTrue(29073 not in records_df["NCBI_taxid"].unique())
        self.assertTrue(10181 not in records_df["NCBI_taxid"].unique())
        self.assertTrue(9999 in records_df["NCBI_taxid"].unique())
        for tid in rest_df["NCBI_taxid"]:
            self.assertTrue(tid not in records_df["NCBI_taxid"].unique())

        from utility.directoryUtility import ucsc_taxid_dict
        for taxid in ucsc_taxid_dict:
            self.assertTrue(taxid not in records_df["NCBI_taxid"].unique())

        fasta_ids = [fasta.id for fasta in SeqIO.parse(filtered_outpath ,'fasta')]
        self.assertTrue(len(fasta_ids) == len(records_df))

    def test_analysis_subset2(self):
        # Test record filtering for UCSC analysis species (speTri)
        test_combined_fpath = "tests/test_data/ucsc_raw/boreoeutheria/ENST00000002596.6.fasta"
        test_records_fpath = "tests/tmp/ENST00000002596.6_43179records.tsv"
        test_filtered_path = "tests/tmp/ENST00000002596_6_43179filtered.fasta"
        full_test_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
        for path in [test_records_fpath ,test_filtered_path]:
            if os.path.exists(path):
                os.remove(path)
        self.assertFalse(os.path.exists(test_filtered_path))
        self.assertTrue(10181 in full_test_df["NCBI_taxid"].unique())

        test_taxid, test_source = 43179,"UCSC"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_combined_fpath ,test_records_fpath,test_taxid,
                                                                      test_source=test_source,
                                                                      filtered_outpath=test_filtered_path)
        self.assertTrue(os.path.exists(filtered_outpath))
        self.assertTrue(len(records_df.loc[records_df['NCBI_taxid' ] ==10181 ,:]) == 0)

        for taxid in ucsc_taxid_dict:
            if taxid != 43179:
                self.assertTrue(taxid not in records_df["NCBI_taxid"].unique())
        fasta_ids = [fasta.id for fasta in SeqIO.parse(filtered_outpath ,'fasta')]
        self.assertTrue(len(fasta_ids) == len(records_df))

        with pd.option_context('display.max_columns',None):
            display(records_df)

    def test_drop_incomplete_records(self):
        with pd.option_context('display.max_columns' ,None):
            test_tid = "ENST00000422628.5"
            test_combined_fpath = "tests/test_data/combined/{0}.fasta".format(test_tid)
            test_records_outpath = "tests/tmp/{0}_9999records.tsv".format(test_tid)
            filtered_tmp_fpath = "tests/tmp/{0}_9999filtered.fasta".format(test_tid)
            combined_df = fautil.load_UCSC_NCBI_df(test_combined_fpath)
            ncbi_idx = combined_df.loc[(combined_df["NCBI_taxid"] == 9999) & ~(combined_df.index.str.contains("ENST"))].index
            boreo_df, rest_df = fautil.partition_UCSC_by_clade(combined_df, 'boreoeutheria')

            # Test proper use (NCBI_record_subset of length 1 raises SequenceAnalysisError for records with no suitable
            # UCSC record data
            with self.assertRaises(UCSCerrors.SequenceAnalysisError):
                filtered_records_df, outpath = ar_filt.filter_analysis_subset(test_combined_fpath ,test_records_outpath,
                                                                              9999,test_source="NCBI",
                                                                              UCSC_analysis_subset=boreo_df.index
                                                                              ,NCBI_record_subset=ncbi_idx,
                                                                              filtered_outpath=filtered_tmp_fpath)
            # Test multiple record NCBI subset still raises SAE if no usable UCSC records present
            ncbi_idx = combined_df.loc[~(combined_df.index.str.contains("ENST"))].index
            with self.assertRaises(UCSCerrors.SequenceAnalysisError):
                filtered_records_df, outpath = ar_filt.filter_analysis_subset(test_combined_fpath,test_records_outpath,
                                                                              9999,test_source="NCBI",
                                                                              UCSC_analysis_subset=boreo_df.index,
                                                                              NCBI_record_subset=ncbi_idx,
                                                                              filtered_outpath=filtered_tmp_fpath)

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
                                                                          filtered_outpath=filtered_tmp_fpath)
            ncbi_len = combined_df.loc[ncbi_idx ,'length'].iloc[0]
            sub_thresh_len = filtered_records_df.loc[filtered_records_df['length'] < ncbi_len* 0.8]
            self.assertTrue(len(sub_thresh_len) == 0)

    def test_id_record_check(self):
        from utility.directoryUtility import taxid_dict, dir_vars
        from utility.NCBIfilter import CLOSEST_EVO_TAXIDS, SECONDARY_EVO_TAXIDS, length_metrics_df
        from numpy import isnan
        from query import orthologUtility as orutil

        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
        xref_df = orutil.load_NCBI_xref_table(xref_fpath)
        lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath, taxid_dict)

        test_taxid, test_tid = 9999, "ENST00000367772.8"
        test_exon_ins_taxid, test_exon_ins_tid = 10181, "ENST00000602312.2"
        test_dist_evos_taxid, test_dist_evos_tid = 9994, "ENST00000359326.9"
        test_dist_evos_taxid, test_dist_evos_tid = 9994, "ENST00000003583.12"
        # test_taxid,test_tid = test_exon_ins_taxid,test_exon_ins_tid
        # test_taxid, test_tid = test_missing_records_taxid,test_missing_records_tid
        test_taxid, test_tid = test_dist_evos_taxid, test_dist_evos_tid
        print("Test TaxID: {0}, Test Transcript ID: {1}".format(test_taxid, test_tid))
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids, quiet=True)
        self.assertTrue(sum(length_checks.values()) == 5)

        ###missing_records_test###
        test_missing_records_taxid, test_missing_records_tid = 29073, "ENST00000354619.10"
        print("Test TaxID: {0}, Test Transcript ID: {1}".format(test_taxid, test_tid))
        test_taxid, test_tid = test_missing_records_taxid, test_missing_records_tid
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath, taxid_dict)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids, quiet=True)
        # Missing records check
        self.assertTrue("9646_check" not in length_checks)
        self.assertTrue("9708_check" in length_checks)
        self.assertAlmostEqual(sum(length_checks.values()) / len(length_checks), 0.5)

        # Assignment test for missing records
        test_spec_lc_df = pd.DataFrame(
            columns=["clade_check", "bestncbi_check", "9646_check", "9708_check", "9713_check"])
        test_spec_lc_df.loc[test_tid, :] = length_checks
        self.assertFalse(
            test_spec_lc_df.loc[test_tid, ["clade_check", "bestncbi_check", "9708_check", "9713_check"]].isna().any())
        self.assertTrue(isnan(test_spec_lc_df.loc[test_tid, "9646_check"]))

        # bad_homology_test
        test_taxid, test_tid = 29073, "ENST00000377627.7"
        print("Test TaxID: {0}, Test Transcript ID: {1}".format(test_taxid, test_tid))
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids, quiet=True)
        self.assertTrue(sum(length_checks.values()) == 0)

        # test ATP5MC1, AGS
        test_taxid, test_tid = 9999, "ENST00000355938.9"
        print("Test TaxID: {0}, Test Transcript ID: {1}".format(test_taxid, test_tid))
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids, quiet=False)

        # Check failed length/id checks:
        lc_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
        lc_df = pd.read_csv(lc_fpath, sep='\t', index_col=0)
        failed_9999 = lc_df.loc[lc_df['9999_check'] == False, :]

        test_taxid = 9999
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]

        iter_idx = 0
        iter_limit = 20
        for tid in failed_9999.index:
            if iter_idx >= iter_limit:
                break
            print("==TID: {0}==".format(tid))
            symbol = xref_df.loc[tid, 'HGNC_symbol']
            print("Gene: {0}".format(symbol))
            test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, tid)
            lm_row = lm_df.loc[tid, :]
            length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids, quiet=False)
            print(length_checks)
            iter_idx += 1

            # self.assertTrue(sum(length_checks.values()) == 0)

    def test_spec_length_check_file_load(self):
        test_fpath = "combined_alignments/hg38AA_knownCanonical/4speciesv1/9999_length_checks.tsv"
        lc_df = pd.read_csv(test_fpath, sep='\t', index_col=0)
        lc_row = lc_df.loc["ENST00000368174.5", :]
        nonna = lc_row.dropna()
        pass_ratio = sum(nonna) / len(nonna)
        self.assertTrue(pass_ratio == 0.8)

        # ncbi_only_test
        lc_row = lc_df.loc["ENST00000381566.6", :]
        nonna = lc_row.dropna()
        ncbi_only_check = (nonna['bestncbi_check'] == True and sum(nonna) == 1)
        self.assertTrue(ncbi_only_check)

        lc_row = lc_df.loc["ENST00000594199.3", :]
        nonna = lc_row.dropna()
        ncbi_only_check = (nonna['bestncbi_check'] == True and sum(nonna) == 1)
        self.assertFalse(ncbi_only_check)

    @unittest.skip("Speed test comparing filtered/ unfiltered input file for identity matrix loading. Faster method "
                   "used in id_length_check")
    def test_id_length_speed(self):
        from utility.directoryUtility import taxid_dict, dir_vars
        import time
        from utility.NCBIfilter import CLOSEST_EVO_TAXIDS, SECONDARY_EVO_TAXIDS, length_metrics_df

        test_taxid, test_tid = 9999, "ENST00000367772.8"
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        # Test speed test for filtered/ unfiltered alignment.
        combined_fasta_fpath = test_comb_fpath
        combined_df = fautil.load_UCSC_NCBI_df(combined_fasta_fpath, taxid_dict)
        ucsc_mask = combined_df.index.str.contains("ENST")
        ucsc_df, ncbi_df = combined_df.loc[ucsc_mask], combined_df.loc[~ucsc_mask]

        filtered_ids = ucsc_df.loc[ucsc_df['NCBI_taxid'].isin(check_taxids)].index
        filtered_ids.append(ncbi_df.index)

        print("==Unfiltered ID DM calcs.==")
        for i in range(10):
            t0 = time.process_time()
            id_dm, align_srs = fautil.construct_id_dm(combined_df.index, combined_fasta_fpath, filtered=True,
                                                      aligned=True)
            assert ((align_srs.index == combined_df.index).all())
            t1 = time.process_time()
            print(t1 - t0)
        print("==Filtered ID DM calcs.==")
        dm_pos = [align_srs.index.get_loc(record_id) for record_id in filtered_ids]
        filtered_dm = id_dm[dm_pos, :]
        filtered_dm = filtered_dm[:, dm_pos]
        print(filtered_dm)
        for i in range(10):
            t0 = time.process_time()
            filtered_ids = ucsc_df.loc[ucsc_df['NCBI_taxid'].isin(check_taxids)].index
            filtered_ids.append(ncbi_df.index)
            tmp_fpath = "tmp/lm_checks.fasta"
            fautil.filter_fasta_infile(filtered_ids, combined_fasta_fpath, tmp_fpath)
            id_dm, align_srs = fautil.construct_id_dm(filtered_ids, tmp_fpath, filtered=True, aligned=True)
            assert ((align_srs.index == filtered_ids).all())
            t1 = time.process_time()
            print(t1 - t0)
        print(id_dm)

    def test_lc_calcs(self):
        from utility.directoryUtility import taxid_dict, dir_vars
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        ar_filt.length_metric_checks(lm_fpath, pass_how='most')
        lc_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
        lc_df = pd.read_csv(lc_fpath, sep='\t', index_col=0)
        lc_counts = lc_df.count()
        lc_passed_counts = lc_df.sum()
        all_passed = [tid for tid in lc_df.index if lc_df.loc[tid, :].sum() == len(lc_df.columns)]
        print("Number of transcript IDs for which all orthologs passed: {0}".format(len(all_passed)))

        # Test length checks for records with no comparable UCSC data
        all_fail_tids = ["ENST00000638261.1", "ENST00000643389.2", "ENST00000643310.2", "ENST00000639909.2",
                         "ENST00000640602.1"]
        for tid in all_fail_tids:
            tid_row = lc_df.loc[tid]
            print(tid_row)