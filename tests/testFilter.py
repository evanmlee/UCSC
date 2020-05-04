import unittest
from query import NCBI_homology
from utility import UCSCfilter, fastaUtility, directoryUtility, UCSCerrors
import pandas as pd
from IPython.display import display
import os

class NCBIHomologyTest(unittest.TestCase):

    def setUP(self):
        pass
    def test_transcript_table(self):
        config, taxid_dict, dir_vars = directoryUtility.config_initialization()
        transcripts_fpath = "reorganized_data/xref_summary_knownCanonical/ensembl_xrefs.tsv"
        xref_tbl = UCSCfilter.load_transcript_table(transcripts_fpath)
        tid_index = xref_tbl.index
        t_tbl = NCBI_homology.SQL_transcript_table(dir_vars,tid_index)
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
        driver = NCBI_homology.headless_driver()
        id_df = pd.DataFrame()
        hgid = "4494"
        tids = ['9999','29073','10181','9994']
        spec_names = ['Urocitellus parryii',"Ursus maritimus","Heterocephalus glaber","Marmota marmota marmota"]
        name_dict = dict(zip(tids,spec_names))
        colname_dict = dict(zip(tids,["{0}_gid".format(tid) for tid in tids]))
        with self.assertRaises(UCSCerrors.NCBIQueryError):
            NCBI_homology.UCSC_NCBI_ID_query(driver,id_df,hgid,hgid,name_dict,
                                             colname_dict)
        n_load_test = 5
        hgid = "6656"
        for _ in range(n_load_test):
            id_df = pd.DataFrame(dtype=str)
            NCBI_homology.UCSC_NCBI_ID_query(driver,id_df,hgid,hgid,name_dict,colname_dict,
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
        tids = ['9999', '29073', '10181', '9994']
        spec_names = ['Urocitellus parryii', "Ursus maritimus", "Heterocephalus glaber", "Marmota marmota marmota"]
        name_dict = dict(zip(tids, spec_names))
        colname_dict = dict(zip(tids, ["{0}_gid".format(tid) for tid in tids]))

        test_outfpath = "tmp/test_orthologs.tsv"
        out_df = NCBI_homology.map_AGS_geneIDs(lim_fpath,name_dict,test_outfpath,"tmp/errors.tsv")
        self.assertTrue(len(out_df.dropna(axis=0,how='any')) == 4)

        test_cache_fpath = 'tmp/test_cache.tsv'
        out_df.loc['2268','9994_gid'] = 'garbage'
        out_df.loc['3075', '9994_gid'] = 'garbage'
        out_df.to_csv(test_outfpath,sep='\t')
        out_df.to_csv(test_cache_fpath, sep='\t')
        out_df = NCBI_homology.map_AGS_geneIDs(lim_fpath, name_dict, test_outfpath, "tmp/errors.tsv",overwrite_gid=['3075'])
        #Test overwrite_gid
        self.assertFalse(out_df.loc['3075', '9994_gid']=='garbage')
        self.assertTrue(out_df.loc['2268','9994_gid'] == 'garbage')
        test_fpaths = [lim_fpath,test_outfpath,test_cache_fpath]
        for fpath in test_fpaths:
            os.remove(fpath)

    def test_xref_load(self):
        #TODO: Check dtypes from loading
        xrefs = NCBI_homology.load_NCBI_xref_table("reorganized_data/xref_summary_knownCanonical/NCBI_xrefs.tsv",gid_dtype='int')
        assert(xrefs["NCBI_gid"].dtype == int)

class UCSCfilterTest(unittest.TestCase):
    def test_NCBIdf(self):
        config, taxid_dict, dir_vars = directoryUtility.config_initialization()
        test_fpath = "{0}/NCBI_alignments/2629.fasta".format(dir_vars["allNCBI_parent"])
        test_df = UCSCfilter.NCBI_fasta_df(test_fpath,taxid_dict)
        self.assertTrue('ANJ16998.1' in test_df.index)
        self.assertFalse('Urocitellus parryii kodiacensis' in test_df["species_name"].unique())
        self.assertTrue('Urocitellus parryii' in test_df["species_name"].unique())
        self.assertTrue(len(test_df.loc[test_df["species_name"]=='Urocitellus parryii',:]) == 10)



if __name__ == "__main__":
    unittest.main()