#%%
from pathlib import Path
import os, sys, subprocess
from collections import OrderedDict
from time import time

from joblib import Parallel, delayed

sys.path.append(str(Path(__file__).parents[0]))
from bam_util import *

_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

class BamPairSorter():
    def __init__(self, path_file, parallel = 1):
        self.path = path_file
        self.parallel = parallel
        
        self.file_handler = None
        
        self.__data_header = None
        self.__offset_after_header = None
        
        self.__dict_refID_to_readpair_pos_and_offset = dict()
        self.__list_readpair_offset_diffchr = list()
        
        self.__dict_refID_to_sorted_readpair_offset = dict()
        
    def run_sorting(self):
        self.__set_file_reader()
        print("Checking BAM Header...", flush = True)
        self.__skip_header()
        print("Checking Read pair coordinates...", flush = True)
        self.__get_read_pair_block_position_and_coordinates()
        print("Sorting Read pairs...", flush = True)
        self.__sort_readpairs_by_coordinate()
        
    def __set_file_reader(self):
        # Open BGzipped BAM file
        self.file_handler = open(self.path, "rb")
        
    def __skip_header(self):
        # Record the annotations and header on the top of BAM file.
        # Also, check the offset where header ends ( = where variant information starts)
        with open(self.path, "rb") as handle:
            bsize, header_data = load_bgzf_block(handle)
            self.__offset_after_header = handle.tell()
        self.__header_size = bsize
        self.__data_header = header_data
        plain_header_text, dict_ref = extract_data_from_binary_bam_header(header_data)
        self.plain_header_text = plain_header_text
        self.dict_refID = dict_ref

    def __get_read_pair_block_position_and_coordinates(self):
        self.file_handler.seek(self.__offset_after_header)
        
        bef_refid = None
        bef_coord = None
        bef_curroffset = None
        bef_blockcoffset = None
        bef_tlen = None
        bef_readname = None
        while 1:
            curr_offset = self.file_handler.tell()
            bsize, bdata = load_bgzf_block(self.file_handler)
            if bdata:
                list_reads, list_read_start_bytes = split_bgzf_block_into_reads(bdata, True)
                
                for read_data, start_bytes_in_bdata in zip(list_reads, list_read_start_bytes):
                    refID = get_refID_on_read_binary_data(read_data)
                    coords = get_pos_on_read_binary_data(read_data)
                    readname = get_readname_on_read_binary_data(read_data)
                    tlen = get_tlen_on_read_binary_data(read_data)
                    
                    if bef_refid == None:
                        bef_refid = refID
                        bef_coord = coords
                        bef_curroffset = curr_offset
                        bef_blockcoffset = start_bytes_in_bdata
                        bef_tlen = tlen
                        bef_readname = readname
                    else:
                        is_samechr = bef_refid == refID
                        
                        assert bef_readname == readname, "Read pair of bam file is not preserved"
                        front_coord = bef_coord if bef_coord<=coords else coords
                        
                        if is_samechr:
                            if self.__dict_refID_to_readpair_pos_and_offset.get(refID) == None:
                                self.__dict_refID_to_readpair_pos_and_offset[refID] = list()
                            self.__dict_refID_to_readpair_pos_and_offset[refID].append(
                                [front_coord, bef_curroffset, bef_blockcoffset, curr_offset, start_bytes_in_bdata]
                            )
                            assert bef_tlen+tlen == 0, "Read pair of bam file is not preserved"
                        else:
                            self.__list_readpair_offset_diffchr.append([bef_curroffset, bef_blockcoffset, curr_offset, start_bytes_in_bdata])
                        
                        bef_refid = None
            else:
                break
        assert bef_refid == None, "Process ended with leftover read"
        
    def __sort_readpairs_by_coordinate(self):
        for refID in self.__dict_refID_to_readpair_pos_and_offset.keys():
            list_readpair_pos_and_offset = self.__dict_refID_to_readpair_pos_and_offset[refID]
            sorted_readpair_by_pos = sorted(list_readpair_pos_and_offset, key = lambda val: val[0])
            sorted_readpair_offsets = list(map(lambda pair_info: pair_info[1:], sorted_readpair_by_pos))
            self.__dict_refID_to_sorted_readpair_offset[refID] = sorted_readpair_offsets
            
    def save_sorted_reads(self, path_save, path_save_diffchr = None):
        print("Save sorted read pairs...")
        path_save_header = f"{path_save}.__tmp.header"
        write_bam_header(path_save_header, self.__data_header)
        
        list_files_concat = [path_save_header]
        for refID in sorted(self.__dict_refID_to_sorted_readpair_offset.keys()):
            list_readpairs_for_writing = self.__dict_refID_to_sorted_readpair_offset[refID]
            threads_for_ref = min(len(list_readpairs_for_writing), self.parallel)
            list_splitted_readpairs_for_threads = self.__split_list_into_n_lists(list_readpairs_for_writing, threads_for_ref)
            list_temp_files_for_threads = list(map(lambda ind: f"{path_save}.__tmp.refID{refID}.{ind}", range(threads_for_ref)))
            with Parallel(n_jobs = threads_for_ref) as parallel:
                list_time_check_results = parallel(delayed(write_part_of_sorted_bam_with_offsets)(
                    self.path,
                    list_temp_files_for_threads[ind],
                    list_splitted_readpairs_for_threads[ind]
                )for ind in range(threads_for_ref))
            path_save_refID = f"{path_save}.__tmp.refID{refID}"
            subprocess.run(f"cat {' '.join(list_temp_files_for_threads)} > {path_save_refID}", shell = True)
            [os.remove(path_tmp) for path_tmp in list_temp_files_for_threads]
            list_files_concat.append(path_save_refID)
        subprocess.run(f"cat {' '.join(list_files_concat)} > {path_save}", shell = True)
        [os.remove(path_tmp) for path_tmp in list_files_concat]
        write_eof(path_save)
        return list_time_check_results
        
    def __split_list_into_n_lists(self, data, n_split):
        n_data_per_split = len(data) // n_split
        n_leftover_data = len(data) % n_split
        
        list_splitted_lists = list()
        ind_start = 0
        for ind in range(n_split):
            n_data_this_ind = n_data_per_split
            if ind < n_leftover_data:
                n_data_this_ind += 1
            list_splitted_lists.append(data[ind_start:ind_start+n_data_this_ind])
            ind_start += n_data_this_ind
        return list_splitted_lists
    
    def close_reader(self):
        self.file_handler.close()
        
    def __del__(self):
        self.close_reader()
        
def get_readname_on_read_binary_data(data):
    byte_start_l_read_name = 8    
    byte_start_read_name = 32
    
    l_read_name = struct.unpack("<B", get_part_of_binary_string(data, byte_start_l_read_name, len_data = 1))[0]
    read_name = get_part_of_binary_string(data, byte_start_read_name, len_data = l_read_name)
    
    return read_name

def get_tlen_on_read_binary_data(data):
    byte_start_tlen = 28  
    tlen = struct.unpack("<i", get_part_of_binary_string(data, byte_start_tlen, len_data = 4))[0]
    
    return tlen
      
def get_refID_on_read_binary_data(data):
    byte_start_refID = 0
    refID = struct.unpack("<i", get_part_of_binary_string(data, byte_start_refID, len_data = 4))[0]
    
    return refID
      
def get_pos_on_read_binary_data(data):
    byte_start_pos = 4
    pos = struct.unpack("<i", get_part_of_binary_string(data, byte_start_pos, len_data = 4))[0]
    
    return pos

def write_part_of_sorted_bam_with_offsets(path_bam, path_save, list_read_offsets):    
    file_reader = open(path_bam, "rb")
    file_writer = open(path_save, "wb")
    
    cache = OrderedDict()
    buffer = b''
    
    def read_block_data_from_offset_cached(file_reader, offset, maxcache = 1000):
        if offset in cache:
            cache.move_to_end(offset)
            return cache[offset]
        elif len(cache) == maxcache:
            cache.popitem(last = False)
        cache[offset] = read_block_data_from_offset(file_reader, offset, ignore_checking = True)
        return cache[offset]
    
    for pair_offsets in list_read_offsets:
        read1_offset, read1_startbytes, read2_offset, read2_startbytes = pair_offsets
        read1_block_data = read_block_data_from_offset_cached(file_reader, read1_offset)
        read1_data = get_single_read_data_from_block_data(read1_block_data, read1_startbytes)
        
        read2_block_data = read_block_data_from_offset_cached(file_reader, read2_offset)
        read2_data = get_single_read_data_from_block_data(read2_block_data, read2_startbytes)
        
        readpair_data = read1_data+read2_data
        if len(buffer) + len(readpair_data) >= 65536:
            write_block(file_writer, buffer)
            buffer = readpair_data
        else:
            buffer += readpair_data
        
    if len(buffer) > 0:
        write_block(file_writer, buffer)
    
    file_reader.close()
    file_writer.flush()
    file_writer.close()

#%%
if __name__ == "__main__":
    path_bam = "/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.small_test.bam"
    
    path_save_sorted = "/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.small_test.sorted_by_pair.parallel1.bam"
    bps = BamPairSorter(path_bam, 1)
    bps.run_sorting()
    list_time_check_results = bps.save_sorted_reads(path_save_sorted)
# %%
