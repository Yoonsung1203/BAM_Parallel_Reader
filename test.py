#%%
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parents[0]))
from bam_util import *
from sort_bam_with_pairedread import *
from bam_parallel_reader import *

# %%
path_bam = "/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.bam"
open_handle = open(path_bam, "rb")

#%%
header_bsize, header_data = load_bgzf_block(open_handle)
next_bsize, next_data = load_bgzf_block(open_handle)
list_reads = split_bgzf_block_into_reads(next_data)
dict_data = extract_data_from_binary_read(list_reads[94])
# %%
# CIGAR:
# op = int(bin({value})[-4:], 2)
# op_len = int(bin({value})[:-4], 2)
# %%
path_gzi = "/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.bam.gzi"
list_coffset, list_ucoffset = read_bgzip_index(path_gzi)
# %%
bpr = BAMParallelReader("/BiO/Access/yoonsung/Research/Test_bam_parallelize/U10K-00751_L01_R1.trimmed_bismark_bt2_pe.deduplicated.sorted_by_pair.parallel10.chr1.bam", 10)
bpr.split_bgzip_bam_into_multiple_readers()
bpr.reset_bgzip_bam_readers()

def count_n_reads(bam_reader):
    bam_reader.set_file_handler()
    cnt_read = 0
    for read in bam_reader:
        cnt_read += 1
    return cnt_read

from joblib import Parallel, delayed

with Parallel(5) as parallel:
    list_n_read = parallel(delayed(count_n_reads)(bam_reader) for bam_reader in bpr.list_splitted_bam_reader)
print(sum(list_n_read))
# %%
