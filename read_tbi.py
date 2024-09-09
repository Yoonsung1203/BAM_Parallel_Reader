#%%
import struct


def read_tbi(path_tbi):
    file_handler = open(path_tbi, "rb")
    
    magic = file_handler.read(4)
    assert magic == b"TBI\x01", f"The file {path_tbi} seems not a TBI format"
    
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack(
        "<iiiiiiii", file_handler.read(32)
    )
    names = file_handler.read(l_nm)
    
    for ind_ref in range(n_ref):
        n_bin = struct.unpack("<i", file_handler.read(4))[0]
        for ind_bin in range(n_bin):
            bin_number, n_chunk = struct.unpack("<Ii", file_handler.read(8))
            for ind_chunk in range(n_chunk):
                cnk_beg, cnk_end = struct.unpack("<QQ", file_handler.read(16))
        n_intv = struct.unpack("<i", file_handler.read(4))[0]
        for ind_intv in range(n_intv):
            ioff = struct.unpack("<Q", file_handler.read(8))[0]
    n_no_coor = struct.unpack("<Q", file_handler.read(8))[0]
            