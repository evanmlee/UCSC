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

    def test_indel_exon_diff_splitting(self):
        from utility.directoryUtility import taxid_dict,dir_vars
        bestncbi_parent = dir_vars['bestNCBI_parent']
        #Test tail insertion
        test_fpath = "{0}/combined/9994/ENST00000371274.8.fasta".format(bestncbi_parent)
        test_records_fpath = "tests/tmp/ENST00000371274.8_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath,test_records_fpath,9994,test_source="NCBI")
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
        test_records_fpath = "tests/tmp/ENST00000003583.12_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath, test_records_fpath,9994,test_source="NCBI")
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

    def test_single_unique_identification(self):
        from utility.directoryUtility import taxid_dict,dir_vars
        bestncbi_parent = dir_vars['bestNCBI_parent']
        test_fpath = "{0}/combined/9994/ENST00000371274.8.fasta".format(bestncbi_parent)
        test_records_fpath = "tests/tmp/ENST00000371274.8_records.tsv"

        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath, test_records_fpath, 9994,
                                                                      test_source="NCBI")
        ncbi_record_id = records_df.loc[~records_df.index.str.contains("ENST"), :].index[0]
        align_df = fautil.align_fasta_to_df(filtered_outpath, ucsc_flag=True)
        unique_thresh = 1
        uniques = ac.find_uniques(align_df, unique_thresh, ncbi_record_id)

        single_unique_pos = ac.id_single_uniques(uniques)
        for test_pos in [206,341,345,373]:
            self.assertTrue(test_pos not in single_unique_pos)

        filt_uni, exon_diff = ac.filter_uniques(uniques,ncbi_record_id,unique_thresh,how='single_only')
        display(filt_uni)

    def test_filter_uniques(self):
        from utility.directoryUtility import taxid_dict, dir_vars
        bestncbi_parent = dir_vars['bestNCBI_parent']
        # Test tail insertion
        test_fpath = "{0}/combined/9994/ENST00000371274.8.fasta".format(bestncbi_parent)
        test_records_fpath = "tmp/ENST00000371274.8_records.tsv"
        records_df, filtered_outpath = ar_filt.filter_analysis_subset(test_fpath, test_records_fpath,test_taxid=9994,test_source="NCBI")
        ncbi_record_id = records_df.loc[~records_df.index.str.contains("ENST"), :].index[0]
        align_df = fautil.align_fasta_to_df(filtered_outpath, ucsc_flag=True)
        unique_thresh = 1
        uniques = ac.find_uniques(align_df, unique_thresh, ncbi_record_id)
        filt_uniques, exon_diffs = ac.filter_exon_diffs(uniques,ncbi_record_id,unique_thresh,how='indel')
        self.assertTrue(115 in filt_uniques.columns)
        self.assertTrue(205 in filt_uniques.columns)
        self.assertTrue(373 not in filt_uniques.columns)
        # with pd.option_context('display.max_columns',None):
        #     print("Filtered unique positions")
        #     display(filt_uniques)
        print("Exon differences")
        display(exon_diffs)

        filt_uniques, exon_diffs = ac.filter_exon_diffs(uniques, ncbi_record_id, unique_thresh, how='all')


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
        records_df, filtered_aln_path = ar_filt.filter_analysis_subset(test_combined_fpath, test_out_records_fpath,9994,
                                                                       test_source="NCBI",UCSC_analysis_subset=UCSC_test_subset,
                                                                       NCBI_record_subset=ncbi_idx)
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

        # cs.overall_summary_table(dir_vars, xref_df, taxid_dict, UCSC_tax_subset=ucsc_tax_subset,
        #                          length_checks_fpath=length_checks_fpath)#,skip_overall=False)

    def test_suppl_calcs1(self):
        from utility.directoryUtility import dir_vars,taxid_dict
        # cs.overall_suppl_calculations(9994,check_percentiles=[99])
        lc_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
        xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
        for taxid in taxid_dict:
            cs.write_background_gene_set(taxid,lc_fpath,xref_fpath,overwrite=True)
            cs.overall_suppl_calculations(taxid, check_percentiles=[95,97,99])

    def test_gene_sets(self):
        test_fpath = "{0}/conservation/gene_sets/9999_analysis_genes.txt".format(directoryUtility.dir_vars['summary_run'])
        print(os.path.exists(test_fpath))

        with open(test_fpath,'rt') as ag_f:
            ag_list = [fline.strip().upper() for fline in ag_f.readlines()]
        print(len(ag_list))

        ag_set = set(ag_list)
        self.assertTrue(len(ag_set)==len(ag_list))


if __name__ == "__main__":
    unittest.main()