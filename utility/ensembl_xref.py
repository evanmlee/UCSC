# Given list of UCSC transcript IDs (ENSEMBL IDs with version numbers), make corresponding GET lookup/id requests
# Store in file specified at ens_tsv_fpath
import pandas as pd
import requests, sys
from requests.exceptions import HTTPError
import time
import os

def write_ensembl_ids(ts_id, write_format="at", ens_tsv_df=pd.DataFrame(),
                      ens_tsv_fpath="reorganized_data/ensembl_stable_ids.tsv"):
    """Given a UCSC transcript id (ts_id), write ENSEMBL stable transcript and parent gene IDs to ens_tsv_fpath.
    ts_id: Should contain version number suffix i.e. ENST############.#
    write_format: If set to "wt", will erase the file at ens_tsv_fpath and start a new file.

    """
    SERVER = "https://rest.ensembl.org"
    if not os.path.exists(ens_tsv_fpath) or write_format == "wt":
        ens_tsv_f = open(ens_tsv_fpath, 'wt')
        ens_tsv_f.write("UCSC_transcript_id\tstable_transcript_id\tstable_gene_id\n")
    else:
        if not ens_tsv_df.empty and ts_id in ens_tsv_df["UCSC_transcript_id"].unique():
            return
        ens_tsv_f = open(ens_tsv_fpath, 'at')
    stable_ts_id_re = re.search("(ENST\d+)\.\d+", ts_id)
    if not stable_ts_id_re:
        print("Poorly formatted transcript ID: {0}".format(ts_id))
        print("transcript ID should be ENSEMBL transcript ID with version number suffix (ENST############.#)")
    else:
        stable_ts_id = stable_ts_id_re.groups()[0]
        ext = "/lookup/id/{0}?".format(stable_ts_id)
        r = requests.get(SERVER + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            status_code = r.status_code
            try:
                r.raise_for_status()
            except HTTPError as e:
                if status_code >= 400 and status_code < 500:
                    headers = r.headers
                    if headers['X-RateLimit-Remaining'] == 0:
                        time.sleep(headers['X-RateLimit-Reset'] + 5)
                        r = requests.get(SERVER + ext, headers={"Content-Type": "application/json"})
                        if not r.ok:
                            r.raise_for_status()
                            sys.exit()
                    else:
                        write_ensembl_errors(ts_id, "UCSC_transcript_id", "lookup/id", status_code, str(e))
                else:
                    print("Status Code: {0}".format(status_code))
                    print("Server error, will retry every 60s")
                    while not r.ok:
                        time.sleep(60)
                        r = requests.get(SERVER + ext, headers={"Content-Type": "application/json"})

        else:
            decoded = r.json()
            repr_dict = dict(decoded)
            parent_gid = repr_dict['Parent']
            ens_tsv_f.write("{0}\t{1}\t{2}\n".format(ts_id, stable_ts_id, parent_gid))
            ens_tsv_f.close()

def ensembl_xml_request(ext,request_type,xml_fpath="",errors_fpath="reorganized_data/ensembl_API_errors.tsv"):
    """Wrapper function for Ensembl API xml type requests. Attempts request, returns request object or writes errors.
    ext: Ensembl REST API endpoint function (ie homology/id/{id} or xrefs/id/{id})
    request_type: string used for logging errors
    xml_fpath: If provided, the results will be saved to a txt file at the file path
    errors_fpath: File for script API error log
    """
    SERVER = "https://rest.ensembl.org"
    xml_r = requests.get(SERVER+ext)
    if not xml_r.ok:
        status_code = xml_r.status_code
        if status_code >= 500:
            time.sleep(60)
            xml_r = requests.get(SERVER+ext)
            if not xml_r.ok:
                xml_r.raise_for_status
                sys.exit()
        elif status_code >= 400:
            try:
                xml_r.raise_for_status()
            except HTTPError as e:
                status_code = xml_r.status_code
                write_ensembl_errors(ens_gid,"Ensembl_gene_id",request_type,status_code,str(e),errors_fpath=errors_fpath)
    if xml_fpath:
        xml_text = xml_r.text
        xml_f = open(xml_fpath,'w')
        xml_f.write(xml_text)
        xml_f.close()
    return xml_r

def get_ensembl_xrefs(gid_list,xml_data_dir):
    for ens_gid in gid_list:
        xml_fpath = "{0}/{1}.txt".format(xml_data_dir,ens_gid)
        if not os.path.exists(xml_fpath):
            ext = "/xrefs/id/{0}?content-type=text/xml".format(ens_gid)
            xref_r = ensembl_xml_request(ext,"xrefs/id",xml_fpath=xml_fpath)


def write_xref_data(ens_tsv_fpath, xref_dirpath, xref_tsv_fpath):
    ens_tsv_df = pd.read_csv(ens_tsv_fpath, delimiter='\t', index_col='UCSC_transcript_id')
    #     gid_list = ens_tsv_df["stable_gene_id"]
    #     tsid_list = ens_tsv_df["UCSC_transcript_id"]

    #     xref_df = pd.DataFrame(columns=['UCSC_transcript_id','stable_transcript_id',\
    #                                 'NCBI_gid','HGNC_symbol','primary_UniProt_ID','other_UniProt_IDs'])
    xref_df = pd.DataFrame(columns=['stable_transcript_id', 'stable_gene_id', \
                                    'NCBI_gid', 'HGNC_symbol', 'primary_UniProt_ID', 'other_UniProt_IDs'])
    #     for ens_gid in gid_list[:]:
    for UCSC_ts_id, ens_tsv_row in ens_tsv_df.iterrows():
        #         ens_tsv_row = ens_tsv_df.loc[UCSC_ts_id,:]

        #         UCSC_ts_id,stable_ts_id = ens_tsv_row["UCSC_transcript_id"],ens_tsv_row["stable_transcript_id"]
        #         xref_row = {'UCSC_transcript_id':UCSC_ts_id,'stable_transcript_id':stable_ts_id}
        stable_ts_id, ens_gid = ens_tsv_row["stable_transcript_id"], ens_tsv_row["stable_gene_id"]
        xref_row = {'stable_transcript_id': stable_ts_id, 'stable_gene_id': ens_gid}
        xml_fpath = "{0}/{1}.txt".format(xref_dirpath, ens_gid)
        xml_f = open(xml_fpath, 'rt')
        xml_text = xml_f.read()
        xml_f.close()
        root = ET.fromstring(xml_text)
        NCBI_elems = root.findall(".//data[@db_display_name='NCBI gene']")
        HGNC_elems = root.findall(".//data[@db_display_name='HGNC Symbol']")
        NCBI_gid = ""
        if len(HGNC_elems) == 1:
            HGNC_elem = HGNC_elems[0]
            HGNC_symbol = HGNC_elem.attrib['display_id']
            xref_row['HGNC_symbol'] = HGNC_symbol
            if len(NCBI_elems) > 1:
                HGNC_matched = [elem for elem in NCBI_elems if elem.attrib['display_id'] == HGNC_symbol]
                assert (len(HGNC_matched) == 1)
                NCBI_gid = HGNC_matched[0].attrib['primary_id']
            elif len(NCBI_elems) == 1:
                NCBI_gid = NCBI_elems[0].attrib['primary_id']
        else:
            if len(NCBI_elems) > 1:
                #                 print("No HGNC Symbol: {0}".format(ens_gid))
                SERVER = "https://rest.ensembl.org"
                ext = "/lookup/id/{0}?content-type=application/json".format(ens_gid)
                lookup_r = requests.get(SERVER + ext)
                ens_display_name = dict(lookup_r.json())['display_name']
                ens_dn_matched = [elem for elem in NCBI_elems if elem.attrib['display_id'] == ens_display_name]
                if len(ens_dn_matched) == 1:
                    NCBI_gid = ens_dn_matched[0].attrib['primary_id']
            elif len(NCBI_elems) == 1:
                NCBI_gid = NCBI_elems[0].attrib['primary_id']
        if NCBI_gid:
            xref_row['NCBI_gid'] = NCBI_gid
        UP_elems = root.findall(".//data[@db_display_name='UniProtKB Gene Name']")
        if len(UP_elems) == 0:
            #         print(ens_gid)
            pass
        elif len(UP_elems) == 1:
            UP_pid = UP_elems[0].attrib['primary_id']
            xref_row['primary_UniProt_ID'] = UP_pid
        else:
            UP_pids = [elem.attrib['primary_id'] for elem in UP_elems]
            UP_dids = [elem.attrib['display_id'] for elem in UP_elems]
            UP_pid = UP_pids[0]
            other_UP_pids = ','.join(UP_pids[1:])
            xref_row['primary_UniProt_ID'] = UP_pid
            xref_row['other_UniProt_IDs'] = other_UP_pids
        xref_row = pd.Series(xref_row, name=UCSC_ts_id)
        xref_df = xref_df.append(xref_row)
    xref_df.index.name = "UCSC_transcript_id"
    xref_df.to_csv(xref_tsv_fpath, sep='\t')
    return xref_df
