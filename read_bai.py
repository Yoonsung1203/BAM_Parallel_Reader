#%%
import struct


def read_bai(path_bai):
    with open(path_bai, "rb") as handler:
        contents = handler.read()
    
    ind_now = 0
    
    assert contents[ind_now:ind_now+4] == b"BAI\x01", "File is not BAI file"
    ind_now += 4
    assert len((struct.unpack("<I", contents[ind_now:ind_now+4]))) == 1
    n_ref = struct.unpack("<I", contents[ind_now:ind_now+4])[0]
    ind_now += 4
    
    dict_bin_per_ref = dict()
    for i_ref in range(n_ref):
        dict_bin_per_ref[i_ref] = {"bins":{}, "intervals":{}}
        assert len((struct.unpack("<I", contents[ind_now:ind_now+4]))) == 1
        n_bin = struct.unpack("<I", contents[ind_now:ind_now+4])[0]
        ind_now += 4
        
        for i_bin in range(n_bin):
            assert len((struct.unpack("<I", contents[ind_now:ind_now+4]))) == 1
            bin_name = struct.unpack("<I", contents[ind_now:ind_now+4])[0]
            ind_now += 4
            dict_bin_per_ref[i_ref]["bins"][i_bin] = {
                "bin":bin_name,
                "chunks":{}
            }
            assert len((struct.unpack("<I", contents[ind_now:ind_now+4]))) == 1
            n_chunk = struct.unpack("<I", contents[ind_now:ind_now+4])[0]
            ind_now += 4
            
            for i_chunk in range(n_chunk):
                assert len((struct.unpack("<Q", contents[ind_now:ind_now+8]))) == 1
                chunk_begin = struct.unpack("<Q", contents[ind_now:ind_now+8])[0]
                ind_now += 8
                assert len((struct.unpack("<Q", contents[ind_now:ind_now+8]))) == 1
                chunk_end = struct.unpack("<Q", contents[ind_now:ind_now+8])[0]
                ind_now += 8
                
                dict_bin_per_ref[i_ref]["bins"][i_bin]["chunks"][i_chunk] = {
                    "begin":chunk_begin,
                    "end":chunk_end
                }
        
        assert len((struct.unpack("<I", contents[ind_now:ind_now+4]))) == 1
        n_interval = struct.unpack("<I", contents[ind_now:ind_now+4])[0]
        ind_now += 4
        
        for i_interval in range(n_interval):
            assert len((struct.unpack("<Q", contents[ind_now:ind_now+8]))) == 1
            ioffset = struct.unpack("<Q", contents[ind_now:ind_now+8])[0]
            ind_now += 8
            
            dict_bin_per_ref[i_ref]["intervals"][i_interval] = ioffset
    
    n_no_coor = contents[ind_now:]
    if len(n_no_coor) > 0:
        assert len(n_no_coor) == 8
        assert len((struct.unpack("<Q", contents[ind_now:ind_now+8]))) == 1
        n_no_coor = struct.unpack("<Q", contents[ind_now:ind_now+8])[0]
        print("Unplaced unmapped reads: ", n_no_coor, sep = '')
    return dict_bin_per_ref
# %%
dict_bai = read_bai("/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.sorted.bam.bai")

from Bio import bgzf
# %%
file_handler = bgzf.BgzfReader("/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.sorted.bam")
# %%
