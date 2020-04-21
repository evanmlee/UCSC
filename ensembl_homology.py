# homology_dirpath = "ensembl_homology_data"
# create_directory(homology_dirpath)
# homology_API_errors_fpath = "{0}/homology_API_errors.tsv".format(homology_dirpath)

#TODO: Clean up, debug
def fetch_ensembl_homology_data(species_str, gid_list, API_errors_fpath, homology_dirpath="ensembl_homology_data"):
    """Downloads ensembl homology data using the REST API. Data acquisition works but homology data itself is
     less preferable than NCBI homology (see NCBI_homology)"""
    homology_spec_dirpath = "{0}/{1}".format(homology_dirpath, species_str)
    create_directory(homology_spec_dirpath)
    SERVER = "https://rest.ensembl.org"
    for ens_gid in gid_list:
        ext = "/homology/id/{0}?content-type=text/xml;compara=vertebrates;aligned=0;target_species={1}".format(ens_gid,
                                                                                                               species_str)
        homology_xml_fpath = "{0}/{1}.txt".format(homology_spec_dirpath, ens_gid)
        if not os.path.exists(homology_xml_fpath):
            homology_r = ensembl_xml_request(ext, "homology/id", xml_fpath=homology_xml_fpath,
                                             errors_fpath=API_errors_fpath)
