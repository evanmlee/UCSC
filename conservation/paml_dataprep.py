from utility.directoryUtility import dir_vars,config, taxid_dict, ucsc_taxid_dict
import os
import re


def process_nh_tree(in_nh_fpath,out_nh_fpath,add_ncbi_species=False,ncbi_subset=[],remove_branch_lengths=True):
    from utility.NCBIfilter import CLOSEST_EVO_TAXIDS
    from utility.fastaUtility import ucsc_tax_table

    if add_ncbi_species:
        if ncbi_subset:
            ce_taxid_values = [CLOSEST_EVO_TAXIDS[k] for k in ncbi_subset]
        else:
            ce_taxid_values = [CLOSEST_EVO_TAXIDS[k] for k in CLOSEST_EVO_TAXIDS]
        ce_taxid_spec_names = [ucsc_tax_table.loc[ucsc_tax_table['NCBI_taxid']==taxid,'tax_name'].iloc[0].replace(" ","_") \
                               for taxid in ce_taxid_values]
        ce_spec_dict = dict(zip(ce_taxid_spec_names,ce_taxid_values))
        ce_tree_modifications = {}
        for ce_spec in ce_spec_dict:
            ce_taxid = ce_spec_dict[ce_spec]
            mod_str = ""
            for taxid in taxid_dict:
                formatted_spec_name = taxid_dict[taxid].replace(" ","_")
                if CLOSEST_EVO_TAXIDS[taxid] == ce_taxid:
                    # if mod_str:
                    mod_str = mod_str + ","+formatted_spec_name
                    # else:
                        # mod_str = formatted_spec_name
            ce_tree_modifications[ce_spec] = mod_str
    with open(in_nh_fpath) as nh_f:
        in_nh_flines = nh_f.readlines()
    with open(out_nh_fpath,'wt') as nh_f:
        indent_level = 0
        for fline in in_nh_flines:
            if remove_branch_lengths:
                fline = re.sub(r":\d\.\d+","",fline)
            indent_level = indent_level + fline.count('(') - fline.count(')')
            # print(indent_level)
            if add_ncbi_species:
                for spec_name in ce_spec_dict:
                    if re.search(spec_name,fline):
                        print(fline)
                        indent_str = ',\n'+(' '*(indent_level-1))
                        indented_clade = ce_tree_modifications[spec_name].replace(',',indent_str)
                        mod_clade = "({0}{1})".format(spec_name,indented_clade)#ce_tree_modifications[spec_name])
                        fline = re.sub(spec_name,mod_clade,fline)
            nh_f.write(fline)
def main():
    # in_nh_fpath = "config/4species1_species_tree_manual_edit.nh"
    in_nh_fpath = "config/hg38.100way.scientificNames.nh"
    out_nh_fpath = "config/4species1_species_tree_updated.nh"
    process_nh_tree(in_nh_fpath,out_nh_fpath,add_ncbi_species=True)

if __name__ == "__main__":
    main()