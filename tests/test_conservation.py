import unittest
from utility import NCBIfilter, fastaUtility as fautil, directoryUtility, UCSCerrors
import pandas as pd
import os
from IPython.display import display
from conservation import analysis_record_filter as ar_filt
from conservation import conservation_summary as cs
from conservation import analysis_calc as ac
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

    def test_id_record_check(self):
        from utility.directoryUtility import taxid_dict,dir_vars
        from utility.NCBIfilter import CLOSEST_EVO_TAXIDS,SECONDARY_EVO_TAXIDS,length_metrics_df
        from numpy import isnan

        test_taxid,test_tid = 9999,"ENST00000367772.8"
        test_exon_ins_taxid, test_exon_ins_tid = 10181, "ENST00000602312.2"
        test_dist_evos_taxid, test_dist_evos_tid = 9994, "ENST00000359326.9"
        test_dist_evos_taxid, test_dist_evos_tid = 9994, "ENST00000003583.12"
        # test_taxid,test_tid = test_exon_ins_taxid,test_exon_ins_tid
        # test_taxid, test_tid = test_missing_records_taxid,test_missing_records_tid
        test_taxid,test_tid = test_dist_evos_taxid, test_dist_evos_tid
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'],test_taxid,test_tid)
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath,taxid_dict)
        lm_row = lm_df.loc[test_tid,:]
        length_checks = ar_filt.id_length_check(test_comb_fpath,lm_row,check_taxids)
        self.assertTrue(sum(length_checks.values())==5)

        ###missing_records_test###
        test_missing_records_taxid, test_missing_records_tid = 29073, "ENST00000354619.10"
        test_taxid, test_tid = test_missing_records_taxid,test_missing_records_tid
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath, taxid_dict)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids)
        #Missing records check
        self.assertTrue("9646_check" not in length_checks)
        self.assertTrue("9708_check" in length_checks)
        self.assertAlmostEqual(sum(length_checks.values())/len(length_checks),0.5)

        #Assignment test for missing records
        test_spec_lc_df = pd.DataFrame(columns=["clade_check","bestncbi_check","9646_check","9708_check","9713_check"])
        test_spec_lc_df.loc[test_tid,:] = length_checks
        self.assertFalse(test_spec_lc_df.loc[test_tid,["clade_check","bestncbi_check","9708_check","9713_check"]].isna().any())
        self.assertTrue(isnan(test_spec_lc_df.loc[test_tid,"9646_check"]))

        #bad_homology_test
        test_taxid, test_tid = 29073, "ENST00000377627.7"
        check_taxids = [CLOSEST_EVO_TAXIDS[test_taxid]] + SECONDARY_EVO_TAXIDS[test_taxid]
        test_comb_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'], test_taxid, test_tid)
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath, taxid_dict)
        lm_row = lm_df.loc[test_tid, :]
        length_checks = ar_filt.id_length_check(test_comb_fpath, lm_row, check_taxids)
        self.assertTrue(sum(length_checks.values())==0)

    def test_spec_length_check_file_load(self):
        test_fpath = "combined_alignments/hg38AA_knownCanonical/4speciesv1/9999_length_checks.tsv"
        lc_df = pd.read_csv(test_fpath,sep='\t',index_col=0)
        lc_row = lc_df.loc["ENST00000368174.5",:]
        nonna = lc_row.dropna()
        pass_ratio = sum(nonna)/len(nonna)
        self.assertTrue(pass_ratio == 0.8)

        #ncbi_only_test
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
        #Test speed test for filtered/ unfiltered alignment.
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
        from utility.directoryUtility import taxid_dict,dir_vars
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        ar_filt.length_metric_checks(lm_fpath,pass_how='most')
        lc_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
        lc_df = pd.read_csv(lc_fpath,sep='\t',index_col=0)
        lc_counts = lc_df.count()
        lc_passed_counts = lc_df.sum()
        all_passed = [tid for tid in lc_df.index if lc_df.loc[tid,:].sum() == len(lc_df.columns)]
        print("Number of transcript IDs for which all orthologs passed: {0}".format(len(all_passed)))

        #Test length checks for records with no comparable UCSC data
        all_fail_tids = ["ENST00000638261.1","ENST00000643389.2","ENST00000643310.2","ENST00000639909.2","ENST00000640602.1"]
        for tid in all_fail_tids:
            tid_row = lc_df.loc[tid]
            print(tid_row)
            # self.assertTrue(sum(tid_row) == 1)

    def test_exon_diff_splitting(self):
        from utility.directoryUtility import taxid_dict,dir_vars
        bestncbi_parent = dir_vars['bestNCBI_parent']
        #Test tail insertion
        test_fpath = "{0}/combined/9994/ENST00000371274.8.fasta".format(bestncbi_parent)
        test_records_fpath = "tmp/ENST00000371274.8_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath,test_records_fpath)
        ncbi_record_id = records_df.loc[~records_df.index.str.contains("ENST"),:].index[0]
        align_df = fautil.align_fasta_to_df(filtered_outpath, ucsc_flag=True)
        unique_thresh = 1
        uniques = ac.find_uniques(align_df,unique_thresh,ncbi_record_id)
        indel_pos = ac.indel_postiions(uniques,ncbi_record_id,unique_thresh)
        self.assertTrue(346 in uniques.columns)
        self.assertTrue(346 not in indel_pos.columns)
        # display(indel_pos)
        exon_diffs = ac.id_exon_diffs(indel_pos)
        self.assertTrue(len(exon_diffs)==197)
        self.assertTrue(exon_diffs[-1] == 569)

        test_fpath = "{0}/combined/9994/ENST00000003583.12.fasta".format(bestncbi_parent)
        test_records_fpath = "tmp/ENST00000003583.12_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath, test_records_fpath)
        ncbi_record_id = records_df.loc[~records_df.index.str.contains("ENST"), :].index[0]
        align_df = fautil.align_fasta_to_df(filtered_outpath, ucsc_flag=True)
        unique_thresh = 1
        uniques = ac.find_uniques(align_df, unique_thresh, ncbi_record_id)
        indel_pos = ac.indel_postiions(uniques, ncbi_record_id, unique_thresh)
        exon_diffs = ac.id_exon_diffs(indel_pos)
        self.assertTrue(25 in exon_diffs)
        self.assertTrue(39 in exon_diffs)
        self.assertTrue(26 not in exon_diffs)
        self.assertTrue(38 not in exon_diffs)

    def test_filter_uniques(self):
        from utility.directoryUtility import taxid_dict, dir_vars
        bestncbi_parent = dir_vars['bestNCBI_parent']
        # Test tail insertion
        test_fpath = "{0}/combined/9994/ENST00000371274.8.fasta".format(bestncbi_parent)
        test_records_fpath = "tmp/ENST00000371274.8_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath, test_records_fpath)
        ncbi_record_id = records_df.loc[~records_df.index.str.contains("ENST"), :].index[0]
        align_df = fautil.align_fasta_to_df(filtered_outpath, ucsc_flag=True)
        unique_thresh = 1
        uniques = ac.find_uniques(align_df, unique_thresh, ncbi_record_id)
        filt_uniques, exon_diffs = ac.filter_exon_diffs(uniques,ncbi_record_id,unique_thresh)
        self.assertTrue(115 in filt_uniques.columns)
        self.assertTrue(205 in filt_uniques.columns)
        self.assertTrue(373 not in filt_uniques.columns)
        # with pd.option_context('display.max_columns',None):
        #     print("Filtered unique positions")
        #     display(filt_uniques)
        # print("Exon differences")
        # display(exon_diffs)


    def test_gene_summary(self):

        # test_combined_fpath = "tests/test_data/combined/ENST00000002596.6.fasta"

        test_out_records_fpath = "tests/tmp/ENST00000371274.8_9994records.tsv"
        # test_out_summary_fpath = "tests/tmp/ENST00000002596.6_9999summary_gstest.tsv"
        test_taxid = 9994
        tid = "ENST00000371274.8"
        test_combined_fpath = "tests/test_data/combined/{0}/{1}.fasta".format(test_taxid,tid)
        tid_summary_dict = "tests/tmp/{0}".format(tid)
        directoryUtility.create_directory(tid_summary_dict)
        combined_fasta_df = fautil.load_UCSC_NCBI_df(test_combined_fpath, TEST_ANALYSIS_DICT)
        ncbi_taxid = test_taxid
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
        summary_df = cs.gene_summary_table(align_df, ncbi_idx, blos_df,tid_summary_dict,tid,unique_thresh=1,
                                        use_jsd_gap_penalty=True)
        # self.assertTrue(len(summary_df.loc[summary_df['NCBI Variant Count']==1])==len(summary_df))
        # self.assertTrue(227 in summary_df.index)
        # self.assertTrue(summary_df['Test-Outgroup BLOSUM62'].count() == 1)

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
                                 length_checks_fpath=length_checks_fpath)#,skip_overall=False)

    def test_blosum(self):
        from utility.directoryUtility import dir_vars
        # test_fpath = "combined_alignments/ENST00000054668.5.fasta"


if __name__ == "__main__":
    unittest.main()