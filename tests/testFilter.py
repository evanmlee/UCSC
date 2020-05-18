import unittest
from query import NCBI_homology, sql_orthologs as so, api_orthologs as ao, orthologUtility as orutil
from utility import NCBIfilter, fastaUtility, directoryUtility, UCSCerrors
import pandas as pd
from IPython.display import display
import os

class NCBIHomologyTest(unittest.TestCase):

    def setUP(self):
        pass
    def test_transcript_table(self):
        config, taxid_dict, dir_vars = directoryUtility.config_initialization()
        transcripts_fpath = "reorganized_data/xref_summary_knownCanonical/ensembl_xrefs.tsv"
        xref_tbl = NCBIfilter.load_transcript_table(transcripts_fpath)
        tid_index = xref_tbl.index
        t_tbl = so.SQL_transcript_table(dir_vars,tid_index)
        self.assertTrue(t_tbl["stable_transcript_id"].str.contains("ENST00000367772").any())
        self.assertTrue(t_tbl["version"].str.contains("2").any())
        self.assertTrue([1 in t_tbl.index])

    def test_compare_orthologs(self):
        test_idx = pd.Index([1,2,3,4])
        col1 = pd.Series(index=test_idx,data=["1","3","3",""])
        col2 = pd.Series(index=test_idx, data=["1", "2", "3", "4"])
        test1_df = pd.DataFrame(data={"col":col1})
        test1_fpath = "tmp/test1.tsv"
        test2_fpath = "tmp/test2.tsv"
        test1_df.to_csv(test1_fpath,sep='\t')

        test2_df = pd.DataFrame(data={"col":col2})
        test2_df.to_csv(test2_fpath, sep='\t')
        diff_df = NCBI_homology._compare_formatted_orthologs("tmp/test1.tsv","tmp/test2.tsv")
        # display(diff_df)
        self.assertTrue('2' in diff_df["col"].unique())
        #Check that rows returned are from correct_fpath (test_df2)
        self.assertFalse('3' in diff_df["col"].unique())
        self.assertTrue(4 in diff_df.index)
        os.remove(test1_fpath)
        os.remove(test2_fpath)

    def test_UCSC_NCBI_query(self):
        driver = ao.headless_driver()
        id_df = pd.DataFrame()
        hgid = "4494"
        tids = [9999, 29073, 10181, 9994]
        spec_names = ['Urocitellus parryii',"Ursus maritimus","Heterocephalus glaber","Marmota marmota marmota"]
        name_dict = dict(zip(tids,spec_names))
        colname_dict = dict(zip(tids,["{0}_gid".format(tid) for tid in tids]))
        with self.assertRaises(UCSCerrors.NCBIQueryError):
            ao.UCSC_NCBI_ID_query(driver,id_df,hgid,hgid,name_dict,
                                             colname_dict)
        n_load_test = 5
        hgid = "6656"
        for _ in range(n_load_test):
            id_df = pd.DataFrame(dtype=str)
            ao.UCSC_NCBI_ID_query(driver,id_df,hgid,hgid,name_dict,colname_dict,
                                             initial_timeout=0)
            try:
                self.assertTrue('9999_gid' in id_df.columns)
                self.assertTrue('113175396' in id_df['9999_gid'].unique())
                self.assertTrue(len(id_df.dropna(axis=1).columns) == 2)
            except AssertionError as e:
                display(id_df)
                raise e

    @unittest.skip("Works")
    def test_map_AGS_ids(self):
        tsv_in_full = pd.read_csv("reorganized_data/xref_summary_knownCanonical/NCBI_xrefs.tsv",sep='\t',dtype=str)
        tsv_lim = tsv_in_full.iloc[:5,:]
        lim_fpath = "tmp/xrefs_lim.tsv"
        tsv_lim.to_csv(lim_fpath,sep='\t')
        tids = [9999, 29073, 10181, 9994]
        spec_names = ['Urocitellus parryii', "Ursus maritimus", "Heterocephalus glaber", "Marmota marmota marmota"]
        name_dict = dict(zip(tids, spec_names))
        colname_dict = dict(zip(tids, ["{0}_gid".format(tid) for tid in tids]))

        test_outfpath = "tmp/test_orthologs.tsv"
        out_df = ao.map_AGS_geneIDs(lim_fpath,name_dict,test_outfpath,"tmp/errors.tsv")
        self.assertTrue(len(out_df.dropna(axis=0,how='any')) == 4)

        test_cache_fpath = 'tmp/test_cache.tsv'
        out_df.loc['2268','9994_gid'] = 'garbage'
        out_df.loc['3075', '9994_gid'] = 'garbage'
        out_df.to_csv(test_outfpath,sep='\t')
        out_df.to_csv(test_cache_fpath, sep='\t')
        out_df = ao.map_AGS_geneIDs(lim_fpath, name_dict, test_outfpath, "tmp/errors.tsv",overwrite_gid=['3075'])
        #Test overwrite_gid
        self.assertFalse(out_df.loc['3075', '9994_gid']=='garbage')
        self.assertTrue(out_df.loc['2268','9994_gid'] == 'garbage')
        test_fpaths = [lim_fpath,test_outfpath,test_cache_fpath]
        for fpath in test_fpaths:
            os.remove(fpath)


class NCBIfilterTest(unittest.TestCase):
    def test_NCBIdf(self):
        config, taxid_dict, dir_vars = directoryUtility.config_initialization()
        test_fpath = "{0}/NCBI_alignments/2629.fasta".format(dir_vars["allNCBI_parent"])
        test_df = NCBIfilter.NCBI_fasta_df(test_fpath,taxid_dict)
        self.assertTrue('ANJ16998.1' in test_df.index)
        self.assertFalse('Urocitellus parryii kodiacensis' in test_df["species_name"].unique())
        self.assertTrue('Urocitellus parryii' in test_df["species_name"].unique())
        self.assertTrue(len(test_df.loc[test_df["species_name"]=='Urocitellus parryii',:]) == 10)


    def test_length_metrics(self ):
        from utility.directoryUtility import config, taxid_dict,dir_vars
        from query import orthologUtility as orutil
        xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
        xref_df = orutil.load_NCBI_xref_table(xref_fpath)
        print(type(xref_df.loc["ENST00000367772.8",'NCBI_gid']))
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df = pd.read_csv(lm_fpath,sep='\t',index_col='UCSC_transcript_id')
        missing_gid = lm_df.loc[lm_df['NCBI_gid'].isnull()]

        display(missing_gid)

    def test_evo_rel_aln(self):
        from utility.directoryUtility import config, taxid_dict, dir_vars
        from utility import NCBIfilter as filt

        #Test Case with multiple orthologs, no evolutionary relative data
        test_taxid, test_tid, test_gid = 10181, "ENST00000625469.1", 338599
        with self.assertRaises(UCSCerrors.SequenceDataError):
            evo_record_idx = filt.evo_relative_alignment(dir_vars, test_taxid, test_tid, test_gid)

        # test_taxid,test_tid,test_gid =  10181, "ENST00000381223.9", 9189
        # test_taxid, test_tid, test_gid = 10181, "ENST00000401088.9", 26099
        # test_taxid, test_tid, test_gid = 10181, "ENST00000538487.6", 2914
        # test_taxid, test_tid, test_gid = 10181, "ENST00000279394.7", 9344
        # test_taxid, test_tid, test_gid = 10181, "ENST00000602312.2", 100288072
        # test_taxid, test_tid, test_gid = 10181, "ENST00000438661.2", 3749
        # test_taxid, test_tid, test_gid = 10181, "ENST00000367485.4", 10216
        # test_taxid, test_tid, test_gid = 10181, "ENST00000373235.4", 339488
        # test_taxid, test_tid, test_gid = 10181, "ENST00000526589.5", 3964
        # test_taxid, test_tid, test_gid = 10181, "ENST00000372827.8", 656
        # test_taxid, test_tid, test_gid = 10181, "ENST00000375105.8", 26279
        # test_taxid, test_tid, test_gid = 10181, "ENST00000322875.8", 4179
        # test_taxid, test_tid, test_gid = 10181, "ENST00000236051.2", 10969
        # test_taxid, test_tid, test_gid = 10181, "ENST00000244174.10", 3581
        # test_taxid, test_tid, test_gid = 10181, "ENST00000440428.5", 92002
        # test_taxid, test_tid, test_gid = 10181, "ENST00000373040.4", 29935
        # test_taxid, test_tid, test_gid = 10181, "ENST00000357033.8", 1756
        # test_taxid, test_tid, test_gid = 10181, "ENST00000372466.8", 54830
        test_taxid, test_tid, test_gid = 10181, "ENST00000374005.8", 2268
        evo_record_idx = filt.evo_relative_alignment(dir_vars, test_taxid, test_tid, test_gid)

        er_aln_outpath = "{0}/evo_relative_alignments/{1}/{2}.fasta".format(dir_vars['allNCBI_parent'],test_taxid,test_tid)
        check_paths = ["tmp/evo_relatives_subset.fasta","tmp/all_orthologs.fasta",er_aln_outpath]
        for path in check_paths:
            self.assertTrue(os.path.exists(path))
        filt.write_bestNCBI_record(dir_vars,test_taxid,test_tid,test_gid,evo_record_idx,ucsc_clade='boreoeutheria')

    def test_new_length_metrics(self):

        from utility import NCBIfilter as filt
        from utility.directoryUtility import dir_vars,taxid_dict, empty_directory
        # test_taxid, test_tid, test_gid = 10181, "ENST00000357033.8", 1756
        test_taxid,test_tid,test_gid =  10181, "ENST00000381223.9", 9189

        empty_directory("tmp")
        test_lm_fpath = "tmp/length_metrics.tsv"
        self.assertFalse(os.path.exists(test_lm_fpath))
        # lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        lm_df, ordered_ucsc_taxids = filt.length_metrics_df(test_lm_fpath,taxid_dict)
        self.assertTrue(len(lm_df) == 0)
        filt.length_metrics_row(dir_vars, taxid_dict, lm_df,ordered_ucsc_taxids,test_tid,test_gid,
                                 ucsc_clade='boreoeutheria')
        self.assertTrue(len(lm_df) == 1)
        display(lm_df)
        self.assertTrue(lm_df.loc[test_tid,'43179_length'] == 0)
        self.assertTrue(lm_df.loc[test_tid, '43179_length'] == 0)
        self.assertTrue(type(lm_df.loc[test_tid,'NCBI_hgid']) == str)

    def test_ucsc_partitioning(self):
        from utility import fastaUtility as fautil
        from Bio import SeqIO
        test_fpath = "tests/test_data/ucsc_raw/ENST00000381223.9.fasta"
        ucsc_df = fautil.UCSC_fasta_df(test_fpath)
        incl,rest = fautil.partition_UCSC_by_clade(ucsc_df,'all')
        self.assertTrue(len(rest)==0)

        be_records = fautil.UCSC_subtax_generator(test_fpath,clade='boreoeutheria')
        test_outpath = "tests/tmp/ENST00000381223.9_boreoeutheria.fasta"
        SeqIO.write(be_records,test_outpath,'fasta')

        self.assertTrue(os.path.exists(test_outpath))
        filt_df = fautil.UCSC_fasta_df(test_outpath)
        self.assertTrue(len(filt_df)==51)
        self.assertFalse(9785 in filt_df['NCBI_taxid'].unique())
        self.assertTrue(143302 in filt_df['NCBI_taxid'].unique())



if __name__ == "__main__":
    unittest.main()