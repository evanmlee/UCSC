import utility.directoryUtility as dirutil
import utility.fastaUtility as fautil
config, dir_vars = dirutil.config,dirutil.dir_vars

import re
import os

def split_raw_UCSC_fasta(seq_format):
    dataset_config = config['DATASET']
    hg_version,dataset_name = dataset_config['HGVersion'],dataset_config['DatasetName']
    reorganized_dir_path = dir_vars['reorg_parent']
    reorganized_seq_dir_path = "{0}/{1}{2}_{3}/all".format(reorganized_dir_path,hg_version,seq_format,dataset_name)
    FASTA_PATH = "hg38_multiz100way/{0}.exon{1}.fa".format(dataset_name, seq_format)
    TRANSCRIPT_LIST_PATH = "{0}/hg38{1}_{2}_transcripts.txt".format(reorganized_dir_path, seq_format, dataset_name)
    transcript_list_f = open(TRANSCRIPT_LIST_PATH, 'at')

    curr_ts_id = "@@@"
    # Skip this block

    # EOF = True
    # header_dict,seq_dict = {}, {}
    block_idx = 0

    with open(FASTA_PATH, 'rt') as fasta_f:
        for i in range(1000000):
            if EOF:
                break
                # Iterate over blocks of 201 lines (200 lines corresponding to alternating UCSC Fasta headers and sequence data)
            # Followed by one blank line
            for j in range(200):
                fline = fasta_f.readline()
                if j == 0:
                    # Set ts_id to be the corresponding Ensembl transcript ID for this fasta block.
                    try:
                        ts_id = re.search(">(ENST\d+\.\d+)_", fline).groups()[0].strip()
                    except:
                        # Either end of transcript (1+ blank lines), or end of file; attempt to read 5 lines
                        for n_attempts in range(5):
                            fline = fasta_f.readline()
                            attempted_match = re.search(">(ENST\d+\.\d+)_", fline)
                            if attempted_match:
                                break
                        if attempted_match:
                            ts_id = re.search(">(ENST\d+\.\d+)_", fline).groups()[0].strip()
                        # Fails re.search if end of file.
                        else:
                            EOF = True
                            ts_id = curr_ts_id
                    out_fasta_fpath = "{0}/{1}.fasta".format(reorganized_seq_dir_path, ts_id)
                    if os.path.exists(out_fasta_fpath):
                        for _ in range(199):
                            fasta_f.readline()
                        block_idx = 0
                        break
                    # Start of new transcript block
                    # 1: append transcript_id to transcript list file
                    # 2: reinitialize all variables for merging each block of UCSC multiz into one line
                    # 3: Open new fasta file
                    if not ts_id == curr_ts_id or EOF:
                        # Below condition ignores first time initializing new transcript block
                        # Format and write output at this point (single line seq records, eliminate
                        # extraneous information in UCSC headers and keep only transcript/genome name and coordinates)
                        if not curr_ts_id == "@@@":
                            #                         if curr_ts_id not in ts_id_set:
                            #                             ts_id_set.add(curr_ts_id)
                            transcript_list_f.write(curr_ts_id + "\n")
                            out_fasta_fpath = "{0}/{1}.fasta".format(reorganized_seq_dir_path, curr_ts_id)
                            out_fasta_f = open(out_fasta_fpath, 'wt')
                            #                         print(curr_ts_id)
                            #                         print(out_fasta_fpath)
                            for record_tag in header_dict:
                                header_line = header_dict[record_tag]
                                seq_line = seq_dict[record_tag] + "\n"
                                if record_tag in f_coords_dict:
                                    f_coords = f_coords_dict[record_tag]
                                    l_coords = l_coords_dict[record_tag]
                                    pos_re = "([\w\.]+:)(\d+)\-(\d+)([\-\+])"
                                    try:
                                        chr_, f1, f2, strand = re.search(pos_re, f_coords).groups()
                                        _, l1, l2, _ = re.search(pos_re, l_coords).groups()
                                        coords_str = "{0}{1}-{2}{3}".format(chr_, f1, l2, strand)
                                        formatted_header = "{0} {1}\n".format(record_tag, coords_str)
                                    except:
                                        print("Reg exp error")
                                        print(f_coords)
                                        print(l_coords)
                                else:
                                    # No coords info/ no sequence record. Use record_tag as formatted_header
                                    formatted_header = "{0}\n".format(record_tag)
                                out_fasta_f.write(formatted_header)
                                out_fasta_f.write(seq_line)
                            out_fasta_f.close()
                        if EOF:
                            break
                        curr_ts_id = ts_id
                        block_idx = 0
                        header_dict, seq_dict, f_coords_dict, l_coords_dict = {}, {}, {}, {}
                # First block: extract number of blocks from header lines (n_blocks), extract record_tag
                # (species/ genome) information for each sequence and store 1) first block genome coordinates
                # (f_coords) 2) header line information (other info encoded in UCSC fasta headers) and 3) first
                # block of sequence data in seq_dict
                if block_idx == 0:
                    if j % 2 == 0:
                        if j == 0:
                            schema_re = re.search(">ENST\d+\.\d+_[a-zA-Z0-9]+_\d+_(\d+) ", fline)
                            n_blocks = schema_re.groups()[0]
                        record_tag = re.search("(>ENST\d+\.\d+_[a-zA-Z0-9]+)_", fline).groups()[0]
                        schema_re = re.search(">ENST\d+\.\d+_[a-zA-Z0-9]+_\d+_\d+ \d+ \d \d ([\w\.:\+\-]+)", fline)
                        if schema_re:
                            # If coordinates provided = UCSC record for that species, store in f_coords_dict
                            f_coords = schema_re.groups()[0]
                            f_coords_dict[record_tag] = f_coords
                            l_coords_dict[record_tag] = f_coords
                        header_dict[record_tag] = fline
                    elif j % 2 == 1:
                        # even lines = sequence data, store in seq_dict
                        seq_dict[record_tag] = fline.strip()
                else:
                    if j % 2 == 0:
                        # If this block has sequence data for record_tag, update last coordinates
                        # Ensures that l_coords_dict has last possible coordinate positions corresponding to record_tag
                        record_tag = re.search("(>ENST\d+\.\d+_[a-zA-Z0-9]+)_", fline).groups()[0]
                        schema_re = re.search(">ENST\d+\.\d+_[a-zA-Z0-9]+_\d+_\d+ \d+ \d \d ([\w\.:\+\-]+)", fline)
                        if schema_re:
                            l_coords = schema_re.groups()[0]
                            l_coords_dict[record_tag] = l_coords

                    elif j % 2 == 1:
                        # Append block sequence data into entry in seq_dict
                        prev_line = seq_dict[record_tag]
                        seq_dict[record_tag] = prev_line + fline.strip()
                        #             out_fasta_f.write(fline)
            # End of block. Skip blank line delimiting blocks, increment block_idx
            fasta_f.readline()
            block_idx += 1
    transcript_list_f.close()

