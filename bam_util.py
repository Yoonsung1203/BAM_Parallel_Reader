#%%
import numpy as np
import struct,math,re,sys,zlib
from time import time
from functools import lru_cache


_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bytes_BC = b"BC"

_bam_read_binary_format_order = [
    "refID",
    "pos",
    "l_read_name",
    "mapq",
    "bin",
    "n_cigar_op",
    "flag",
    "l_seq",
    "next_refID",
    "next_pos",
    "tlen"
]
_bam_read_binary_format = {
    "refID" : {"fmt":"<i", "byte":4},
    "pos" : {"fmt":"<i", "byte":4},
    "l_read_name" : {"fmt":"<B", "byte":1},
    "mapq" : {"fmt":"<B", "byte":1},
    "bin" : {"fmt":"<H", "byte":2},
    "n_cigar_op" : {"fmt":"<H", "byte":2},
    "flag" : {"fmt":"<H", "byte":2},
    "l_seq" : {"fmt":"<I", "byte":4},
    "next_refID" : {"fmt":"<i", "byte":4},
    "next_pos" : {"fmt":"<i", "byte":4},
    "tlen" : {"fmt":"<i", "byte":4},
}
_bam_tag_type = {
    'A':{"fmt":'<c', 'byte_len':1}, # char
    'c':{"fmt":'<b', 'byte_len':1}, # signed char, int8_t
    'C':{"fmt":'<B', 'byte_len':1}, # unsigned char, uint8_t
    's':{"fmt":'<h', 'byte_len':2}, # signed short, int16_t
    'S':{"fmt":'<H', 'byte_len':2}, # unsigned short, uint16_t
    'i':{"fmt":'<i', 'byte_len':4}, # signed int, int32_t
    'I':{"fmt":'<I', 'byte_len':4}, # unsigned int, uint32_t
    'q':{"fmt":'<q', 'byte_len':8}, # signed long long, int64_t
    'Q':{"fmt":'<Q', 'byte_len':8}, # unsigned long long, uint64_t
    'f':{"fmt":'<f', 'byte_len':4}, # float, float
    'd':{"fmt":'<d', 'byte_len':8}  # double, double
}
_cigar_op_str = dict(zip(list(range(9)), "MIDNSHP=X"))

def read_bgzip_index(path_gzi):
    list_coffset_of_block_start = [0]
    list_ucoffset_of_block_start = [0]
    with open(path_gzi, "rb") as handle:
        n_entries = struct.unpack("<Q", handle.read(8))[0]
        
        for _ in range(n_entries):
            coffset = struct.unpack("<Q", handle.read(8))[0]
            ucoffset = struct.unpack("<Q", handle.read(8))[0]
            
            list_coffset_of_block_start.append(coffset)
            list_ucoffset_of_block_start.append(ucoffset)
            
    return list_coffset_of_block_start, list_ucoffset_of_block_start

def make_virtual_offset_from_bytes(bsize):
    return bsize << 16

def split_virtual_offset(virtual_offset):
    bsize = virtual_offset >> 16
    within_block = virtual_offset ^ (bsize << 16)
    return bsize, within_block    

def load_bgzf_block(handle, text_mode=False):
    """Load the next BGZF block of compressed data (PRIVATE).

    Returns a tuple (block size and data), or at end of file
    will raise StopIteration.
    ** Copied code from Biopython Github **
    """
    magic = handle.read(4)
    if not magic:
        # End of file - should we signal this differently now?
        # See https://www.python.org/dev/peps/pep-0479/
        raise StopIteration
    if magic != _bgzf_magic:
        raise ValueError(
            r"A BGZF (e.g. a BAM file) block should start with "
            r"%r, not %r; handle.tell() now says %r"
            % (_bgzf_magic, magic, handle.tell())
        )
    gzip_mod_time, gzip_extra_flags, gzip_os, extra_len = struct.unpack(
        "<LBBH", handle.read(8)
    )

    block_size = None
    x_len = 0
    while x_len < extra_len:
        subfield_id = handle.read(2)
        subfield_len = struct.unpack("<H", handle.read(2))[0]  # uint16_t
        subfield_data = handle.read(subfield_len)
        x_len += subfield_len + 4
        if subfield_id == _bytes_BC:
            if subfield_len != 2:
                raise ValueError("Wrong BC payload length")
            if block_size is not None:
                raise ValueError("Two BC subfields?")
            block_size = struct.unpack("<H", subfield_data)[0] + 1  # uint16_t
    if x_len != extra_len:
        raise ValueError(f"x_len and extra_len differ {x_len}, {extra_len}")
    if block_size is None:
        raise ValueError("Missing BC, this isn't a BGZF file!")
    # Now comes the compressed data, CRC, and length of uncompressed data.
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15)  # Negative window size means no headers
    data = d.decompress(handle.read(deflate_size)) + d.flush()
    expected_crc = handle.read(4)
    expected_size = struct.unpack("<I", handle.read(4))[0]
    if expected_size != len(data):
        raise RuntimeError("Decompressed to %i, not %i" % (len(data), expected_size))
    # Should cope with a mix of Python platforms...
    crc = zlib.crc32(data)
    if crc < 0:
        crc = struct.pack("<i", crc)
    else:
        crc = struct.pack("<I", crc)
    if expected_crc != crc:
        raise RuntimeError(f"CRC is {crc}, not {expected_crc}")
    if text_mode:
        # Note ISO-8859-1 aka Latin-1 preserves first 256 chars
        # (i.e. ASCII), but critically is a single byte encoding
        return block_size, data.decode("latin-1")
    else:
        return block_size, data
    
def load_bgzf_block_compact(handle, text_mode=False):
    """Load the next BGZF block of compressed data (PRIVATE).

    Returns a tuple (block size and data), or at end of file
    will raise StopIteration.
    ** Copied code from Biopython Github **
    """
    cblock = handle.read(18)
    extra_len = struct.unpack("<H", cblock[10:12])[0]

    block_size = struct.unpack("<H", cblock[16:18])[0] + 1  # uint16_t
    
    # Now comes the compressed data, CRC, and length of uncompressed data.
    deflate_size = block_size - 1 - extra_len - 19
    
    data_with_sizeinfo = handle.read(deflate_size+8)
    d = zlib.decompressobj(-15)  # Negative window size means no headers
    data = d.decompress(data_with_sizeinfo[:deflate_size]) + d.flush()
    
    return data

def split_bgzf_block_into_reads(block_data, return_read_start_bytes = False):
    ind_check = 0
    list_reads_data = list()
    list_read_start_bytes = list()
    while 1:
        list_read_start_bytes.append(ind_check)
        block_size = struct.unpack("<I", get_part_of_binary_string(block_data, ind_check, len_data = 4))[0]
        read_data = get_part_of_binary_string(block_data, ind_check+4, len_data = block_size)
        list_reads_data.append(read_data)
        ind_check += (4+block_size)
        if ind_check >= len(block_data):
            break
    if return_read_start_bytes:
        return list_reads_data, list_read_start_bytes
    else:
        return list_reads_data

def extract_readid_from_binary_read(read_data):
    dict_data = dict()
    ind_check = 0
    for bin_key in _bam_read_binary_format_order:
        binary_format = _bam_read_binary_format[bin_key]
        binary_read_data_part = get_part_of_binary_string(read_data, ind_check, len_data = binary_format["byte"])
        dict_data[bin_key] = struct.unpack(binary_format["fmt"], binary_read_data_part)[0]
        ind_check += binary_format["byte"]
    
    dict_data["read_name"] = get_part_of_binary_string(read_data, ind_check, len_data = dict_data["l_read_name"])
    return dict_data["read_name"]

def extract_data_from_binary_read(read_data):
    dict_data = dict()
    
    ind_check = 0
    for bin_key in _bam_read_binary_format_order:
        binary_format = _bam_read_binary_format[bin_key]
        binary_read_data_part = get_part_of_binary_string(read_data, ind_check, len_data = binary_format["byte"])
        dict_data[bin_key] = struct.unpack(binary_format["fmt"], binary_read_data_part)[0]
        ind_check += binary_format["byte"]
    
    dict_data["read_name"] = get_part_of_binary_string(read_data, ind_check, len_data = dict_data["l_read_name"])
    ind_check += dict_data["l_read_name"]
    
    dict_data["cigar"] = []
    for _ in range(dict_data["n_cigar_op"]):
        binary_cigar = get_part_of_binary_string(read_data, ind_check, len_data = 4)
        val_cigar_thisop = struct.unpack("<I", binary_cigar)[0]
        dict_data["cigar"].append(val_cigar_thisop)
        ind_check += 4
    
    bytes_seq = math.floor((dict_data["l_seq"]+1)/2)
    dict_data["seq"] = get_part_of_binary_string(read_data, ind_check, len_data=bytes_seq)
    ind_check += bytes_seq
    
    dict_data["qual"] = get_part_of_binary_string(read_data, ind_check, len_data=dict_data["l_seq"])
    ind_check += dict_data["l_seq"]
    
    rest_data = read_data[ind_check:]
    while 1:
        if len(rest_data) == 0: break
        btr_data = bytearray(rest_data)
        n_for_this_tag = 3
        
        tag = get_part_of_binary_string(btr_data, 0, len_data = 2).decode()
        val_type = get_part_of_binary_string(btr_data, 2).decode()
        
        if val_type in _bam_tag_type:
            fmt_value = _bam_tag_type[val_type]["fmt"]
            byte_len_value = _bam_tag_type[val_type]["byte_len"]
            value_binary = get_part_of_binary_string(btr_data, 3, len_data = byte_len_value)
            value = struct.unpack(fmt_value, value_binary)[0]
            n_for_this_tag += byte_len_value
        elif val_type == 'Z':
            # Special case: list of char
            ind_end, value = get_parsed_data_and_end_index_of_binary_characters_until_null(btr_data, 3)
            n_for_this_tag += (ind_end-2)
        elif val_type == 'H':
            # Special case: list of hexhex
            ind_end, value = get_parsed_data_and_end_index_of_binary_characters_until_null(btr_data, 3)
            value = [get_part_of_binary_string(value, ind, len_data = 2) for ind in range(0, len(value), 2)]            
        elif val_type == 'B':
            # Special case: list of typed data
            val_type = get_part_of_binary_string(btr_data, 3).decode()
            val_length = struct.unpack("<I", get_part_of_binary_string(btr_data, 4, len_data=4))
            
            fmt_value = _bam_tag_type[val_type]["fmt"]
            byte_len_value = _bam_tag_type[val_type]["byte_len"]
            
            total_bytes_of_value = byte_len_value * val_length
            value_binary = get_part_of_binary_string(btr_data, 8, total_bytes_of_value)
            value = [struct.unpack(fmt_value, get_part_of_binary_string(value_binary, ind, len_data = byte_len_value)) for ind in range(0, len(value_binary, byte_len_value))]
            n_for_this_tag += 5 + total_bytes_of_value # 5: byte of val_type (1) + bytes of val_length (4)
        else:
            raise Exception(f"Wrong Datatype: {val_type}") 
        dict_data[tag] = value
        rest_data = rest_data[n_for_this_tag:]
    return dict_data

def convert_binary_to_seq(data, len_data):
    dict_val_to_seq = dict(zip(range(16), "=ACMGRSVTWYHKDBN"))
    hex_val = list(map(int,data.hex()))
    seq = ''.join(list(map(dict_val_to_seq.__getitem__, hex_val))[:len_data])
    return seq

def convert_binary_to_phred_qual(data):
    qual = list(data)
    ascii_qual = list(map(lambda val: val + 33, qual))
    str_ascii_qual = ''.join(list(map(chr, ascii_qual)))
    return str_ascii_qual

def convert_cigar_list_to_cigarstring(list_cigar):
    list_str_cigar = list()
    for cigar_op in list_cigar:
        op_len = cigar_op >> 4
        op = cigar_op ^ (op_len << 4)
        op_str = _cigar_op_str[op]
        list_str_cigar.append(str(op_len))
        list_str_cigar.append(op_str)
    return ''.join(list_str_cigar)

def extract_data_from_binary_bam_header(data):
    bam_magic = data[:4]
    assert bam_magic == b"BAM\x01", "BAM file does not start with BAM magic string. Maybe not a BAM file."
    l_text = struct.unpack("<I", get_part_of_binary_string(data, 4, len_data=4))[0]
    text = get_part_of_binary_string(data, 8, len_data=l_text)
    n_ref = struct.unpack("<I", get_part_of_binary_string(data, l_text+8, len_data=4))[0]
    
    dict_refID = dict()
    
    curr_ind = l_text+12
    for ind_ref in range(n_ref):
        l_name = struct.unpack("<I", get_part_of_binary_string(data, curr_ind, len_data=4))[0]
        curr_ind += 4
        name = get_part_of_binary_string(data, curr_ind, len_data=l_name).decode().strip('\x00')
        curr_ind += l_name
        l_ref = struct.unpack("<I", get_part_of_binary_string(data, curr_ind, len_data=4))[0]
        curr_ind += 4
        
        dict_refID[ind_ref] = {"name":name, "l_ref":l_ref}
    return text, dict_refID
        
def get_parsed_data_and_end_index_of_binary_characters_until_null(data, ind_start):
    curr_ind = ind_start
    while 1:
        char_next = get_part_of_binary_string(data, curr_ind)
        if char_next == b'\x00':
            break
        else:
            curr_ind += 1
            continue
    parsed_data = get_part_of_binary_string(data, ind_start, ind_end = curr_ind-1)
    return curr_ind, parsed_data

def get_part_of_binary_string(data, ind_start, ind_end = None, len_data = None):
    '''
    ind_start: Start Index for parsing data
    ind_end: End Index for parsing data (kwargs priority: 1)
    len_data: Length for parsing data (kwargs priority: 2)
    '''
    if ind_end == None or ind_start == ind_end:
        if len_data != None and len_data > 1:
            return data[ind_start:ind_start+len_data]
        else:
            return bytes(chr(data[ind_start]), "latin-1")
    else:
        assert ind_start <= ind_end, "End Index must be higher or equal to Start Index!"
        return data[ind_start:ind_end+1]  
    
def get_compressed_block_of_bam_data(data, compresslevel = 6):
    if len(data) > 65536:
        raise ValueError(f"{len(data)} Block length > 65536")
    compressor = zlib.compressobj(
        compresslevel, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
    )
    compressed_data = compressor.compress(data) + compressor.flush()
    del compressor
    
    bsize = struct.pack("<H", len(compressed_data) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(data) & 0xFFFFFFFF)
    uncompressed_length = struct.pack("<I", len(data))
    
    compressed_block = _bgzf_header + bsize + compressed_data + crc + uncompressed_length
    return compressed_block

def write_bam_header(path_save, header_data):
    file_handler = open(path_save, "wb")
    write_block(file_handler, header_data)
    file_handler.flush()
    file_handler.close()
            
def read_block_data_from_offset(file_handler, coffset, ignore_checking = False):
    file_handler.seek(coffset)
    if ignore_checking:
        block_data = load_bgzf_block_compact(file_handler)
    else:
        bsize, block_data = load_bgzf_block(file_handler)
    return block_data

def get_single_read_data_from_block_data(block_data, start_bytes):
    block_size = struct.unpack("<I", get_part_of_binary_string(block_data, start_bytes, len_data = 4))[0]
    read_data = get_part_of_binary_string(block_data, ind_start=start_bytes, len_data = 4+block_size)
    return read_data
    
def write_block(file_handler, data):
    if isinstance(file_handler, str):
        file_handler = open(file_handler, "ab")
    compressed_data = get_compressed_block_of_bam_data(data)
    file_handler.write(compressed_data)
    
def write_eof(file_handler):
    if isinstance(file_handler, str):
        file_handler = open(file_handler, "ab")
    file_handler.write(_bgzf_eof)
    file_handler.flush()
    file_handler.close()
# %%
