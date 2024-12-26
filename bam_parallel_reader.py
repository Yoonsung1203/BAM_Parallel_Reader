#%%
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parents[0]))
from bam_util import *


_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

def search_nearest_bgzip_block(handler, coffset_start_search, wsize_search = 5000, step_search = 4900):
    handler.seek(coffset_start_search)
    curr_coffset = coffset_start_search
    
    while 1:
        handler.seek(curr_coffset)
        dat_search = handler.read(wsize_search)
        try:
            ind_bgzf_magic = dat_search.index(_bgzf_magic)
            curr_coffset += ind_bgzf_magic
            handler.seek(curr_coffset)
            dat_search = handler.read(wsize_search)
            if validate_bgzip_block_header(dat_search):
                break
            else:
                curr_coffset += 4
                continue
        except Exception:
            curr_coffset += step_search
            continue
    return curr_coffset

def validate_bgzip_block_header(data):
    check_gid1 = data[0] == 31
    check_gid2 = data[1] == 139
    check_gcm = data[2] == 8
    check_gflag = data[3] == 4
    check_si1 = data[12] == 66
    check_si2 = data[13] == 67
    check_slen = struct.unpack("<H", data[14:16])[0] == 2
    
    return check_gid1 + check_gid2 + check_gcm + check_gflag + check_si1 + check_si2 + check_slen == 7 

class BAMParallelReader():
    # Split BGzipped BAM for parallelization
    def __init__(self, path_file, parallel = 1, path_gzi = None):
        self.path = path_file
        self.path_gzi = path_gzi
        self.parallel = parallel
        
        self.file_handler = None
        
        self.__header_size = None
        self.plain_header_text = list()
        self.dict_refID = dict()
        
        self.__offset_eof = None
        
        self.__offsets_for_parallelizing = list()
        
        self.list_splitted_bam_reader = list()
        self.__list_reader_offset_start = list()
        self.__list_reader_offset_end = list()
    
    def split_bgzip_bam_into_multiple_readers(self):
        # Search for offsets to split the BGzipped BAM into multiple portions
        self.__set_file_reader()
        self.__skip_header()
        self.__check_file_start_end_offset()
        self.__split_bam_offset_for_parallelization()
        self.__generate_blockgzipbamreaders_for_parallelizing()
        print(self.__offsets_for_parallelizing)
        
    def __set_file_reader(self):
        # Open BGzipped BAM file
        self.file_handler = open(self.path, "rb")
        
    def __skip_header(self):
        # Record the annotations and header on the top of BAM file.
        # Also, check the offset where header ends ( = where variant information starts)
        with open(self.path, "rb") as handle:
            bsize, header_data = load_bgzf_block(handle)
        self.__header_size = bsize
        plain_header_text, dict_ref = extract_data_from_binary_bam_header(header_data)
        self.plain_header_text = plain_header_text
        self.dict_refID = dict_ref
    
    def __check_file_start_end_offset(self):
        # Find the offsets to split BGzipped BAM file into multiple portions
        
        # Check "End of File (EOF)" context
        # BGzipped file must ends with EOF context. If is not, it means this file is incomplete.
        file_object = open(self.path, "rb")
        filesize = file_object.seek(0, 2)
        filesize_except_eof = filesize - len(_bgzf_eof)
        file_object.seek(filesize_except_eof)
        eof = file_object.read(len(_bgzf_eof))
        file_object.close()
        assert eof == _bgzf_eof, "Block gzip file does not ends with 'end of file' context"
        self.__offset_eof = filesize_except_eof
        
    def __split_bam_offset_for_parallelization(self):
        if self.path_gzi:
            self.__split_bam_offset_for_parallelization_with_gzi()
        else:
            self.__split_bam_offset_for_parallelization_wo_gzi()
    
    def __split_bam_offset_for_parallelization_with_gzi(self):
        list_block_offsets, _ = read_bgzip_index(self.path_gzi)
        list_ind_block_offset_for_parallelizing = list(map(int, np.linspace(1, len(list_block_offsets)-1, self.parallel+1)))
        
        assert len(set(list_ind_block_offset_for_parallelizing)) == len(list_ind_block_offset_for_parallelizing), "Too many cores for small bam file"
        
        self.__offsets_for_parallelizing = list(map(lambda ind: list_block_offsets[ind], list_ind_block_offset_for_parallelizing))       
        
    def __split_bam_offset_for_parallelization_wo_gzi(self):
        # The offset of EOF equals to the end of file
        # Also, BGzipped file needs different offset because of it's "Blocked" nature.
        # Generally, offset is the size of file. (Unit: Bytes)
        # This code splits the size of file into the number of "parallel"
        list___offsets_for_parallelizing = list(map(int, np.linspace(self.__header_size, self.__offset_eof, self.parallel+1)[:-1]))
        list_block_start___offsets_for_parallelizing = list(map(lambda offset_start_search: search_nearest_bgzip_block(self.file_handler, offset_start_search), list___offsets_for_parallelizing))
        
        # This code does not optimize the offset for parallelization
        # If this assertion error occurs, please reduce the number of jobs or just do no parallelize it
        assert len(set(list_block_start___offsets_for_parallelizing)) == len(list_block_start___offsets_for_parallelizing), "Too many cores for small bam file"
        self.__offsets_for_parallelizing = list_block_start___offsets_for_parallelizing
        
    def reset_bgzip_bam_readers(self):
        # Reset the cursor of each split BamPartReader object 
        for bgvr in self.list_splitted_bam_reader:
            del bgvr
        self.list_splitted_bam_reader = list()
        self.__generate_blockgzipbamreaders_for_parallelizing()
    
    def __generate_blockgzipbamreaders_for_parallelizing(self):
        # Generate "parallel" number of BamPartReader object for parallelizing
        # Each object starts reading BAM from each split offsets and ends reading BAM until the file ends, or meet the next start offset
        self.__list_reader_offset_start = self.__offsets_for_parallelizing
        self.__list_reader_offset_end = self.__offsets_for_parallelizing[1:] + [self.__offset_eof]
        self.list_splitted_bam_reader = list(map(lambda bstart, bend: BamPartReader(self.path, bstart, bend), self.__list_reader_offset_start, self.__list_reader_offset_end))        
    
    def __close_reader(self): 
        # Close the BgzfReader for checking Block Information. If it's already closed, don't do anything
        # Also, it parallely closes the BamPartReader objects under this object
        if hasattr(self.file_handler, "close"):
            self.file_handler.close()
        for bam_reader in self.list_splitted_bam_reader:
            del bam_reader
    
    def __del__(self):
        # Automatically close the all file objects while this object is deleted. Necessary for saving resources
        self.__close_reader()
        del self.file_handler
        
class BamPartReader():
    def __init__(self, path_file, block_start = None, block_end = None):
        self.path = path_file
        self.bstart = block_start
        self.bend = block_end
        assert self.bstart <= self.bend, "File reading start offset must be smaller than end offset"
        
        self.file_handler = None
        
        self.__curr_block = None
        self.__curr_reads = list()
        
        self.generator_reads = None
        
    def set_file_handler(self):
        self.file_handler = open(self.path, "rb")
        self.file_handler.seek(self.bstart)
        self.generator_reads = self.__generate_nextread_binary()
        
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.get_nextread_binary()
    
    def __generate_nextread_binary(self):
        while 1:
            self.__read_block()
            if self.__curr_block:
                self.__curr_reads = split_bgzf_block_into_reads(self.__curr_block)
                for single_binaryread in self.__curr_reads:
                    yield single_binaryread
            else:
                break
    
    def get_nextread_binary(self):
        return next(self.generator_reads)
    
    def __read_block(self):
        curr_offset = self.file_handler.tell()
        if curr_offset >= self.bend:
            self.__curr_block = None
        else:
            bsize, self.__curr_block = load_bgzf_block(self.file_handler)
    
    def __close_reader(self):
        if hasattr(self.file_handler, "close"):
            self.file_handler.close()
    
    def __del__(self):
        self.__close_reader()
        del self.file_handler
        
        