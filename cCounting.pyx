from cpython.string cimport PyString_AsString
from libc.stdlib cimport malloc, free
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import numpy as np

def c_count(ref, read, buckets, winSize, k):
    cdef int error = 0

    #convert python variables into C variables
    cdef char *read_seq = PyString_AsString(read)
    cdef char *ref_seq = PyString_AsString(ref)
    cdef int read_len = len(read)
    cdef int ref_len = len(ref)

    cdef int num_rows = buckets[0]
    cdef int num_cols = buckets[1]
    cdef int window_size = winSize
    cdef int k_size = k
    
    cdef int **counts = <int **> PyMem_Malloc(num_rows * sizeof(int*))
    if not counts:
        raise MemoryError()
    
    cdef int row_index, col_index
    
    for row_index in range(num_rows):
        counts[row_index] = <int *> PyMem_Malloc(num_cols * sizeof(int))
        if not counts[row_index]:
            raise MemoryError
        for col_index in range(num_cols):
            counts[row_index][col_index] = 0
    
    #count
    cdef int same, k_index
    for row_index in range(ref_len - k_size + 1):
        for col_index in range(read_len - k_size + 1):
            same = 1
            k_index = 0
            while same and k_index < k_size:
                if ref_seq[row_index+k_index] != read_seq[col_index+k_index]:
                    same = 0
                    break
                k_index += 1
            if same:
                counts[row_index / window_size][col_index / window_size] += 1
    
#~     cdef int offset, index
#~     cdef int same, k_index
#~     if ref_len > read_len:
#~         for offset in range(-ref_len, 0):
#~             k_index = 0
#~             index = 0
#~             while index < (min(ref_len + offset, read_len) - k_size + 1):
#~                 same = 1
#~                 while same and k_index < k_size:
#~                     if ref_seq[index+k_index-offset] != read_seq[index+k_index]:
#~                         same = 0
#~                         index += k_index  
#~                         k_index = 0
#~                         break
#~                     k_index += 1
#~                 if same:
#~                     counts[(index+k_index-offset) / window_size][(index+k_index) / window_size] += 1
#~                     k_index = k_size - 1
#~                 index += 1           
#~         for offset in range(read_len):
#~             k_index = 0
#~             index = 0
#~             while index < (read_len - offset - k_size + 1):
#~                 same = 1
#~                 while same and k_index < k_size:
#~                     if ref_seq[index+k_index] != read_seq[index+k_index+offset]:
#~                         same = 0
#~                         index += k_index
#~                         k_index = 0
#~                         break
#~                     k_index += 1
#~                 if same:
#~                     counts[(index+k_index) / window_size][(index+k_index+offset) / window_size] += 1
#~                     k_index = k_size - 1                
#~                 index += 1            
#~     else:
#~         for offset in range(-ref_len, 0):
#~             k_index = 0
#~             index = 0
#~             while index < (ref_len + offset - k_size + 1):
#~                 same = 1
#~                 while same and k_index < k_size:
#~                     if ref_seq[index+k_index-offset] != read_seq[index+k_index]:
#~                         same = 0
#~                         index += k_index  
#~                         k_index = 0
#~                         break
#~                     k_index += 1
#~                 if same:
#~                     counts[(index+k_index-offset) / window_size][(index+k_index) / window_size] += 1
#~                     k_index = k_size - 1
#~                 index += 1
#~         for offset in range(read_len):
#~             k_index = 0
#~             index = 0
#~             while index < (min(read_len - offset, ref_len) - k_size + 1):
#~                 same = 1
#~                 while same and k_index < k_size:
#~                     if ref_seq[index+k_index] != read_seq[index+k_index+offset]:
#~                         same = 0
#~                         index += k_index  
#~                         k_index = 0
#~                         break
#~                     k_index += 1
#~                 if same:
#~                     counts[(index+k_index) / window_size][(index+k_index+offset) / window_size] += 1
#~                     k_index = k_size - 1
#~                 index += 1

    #convert back to python return type
    np_counts = np.empty((num_rows, num_cols), dtype=int)
    for i in range(num_rows):
        for j in range(num_cols):
            np_counts[i,j] = counts[i][j]
        free(counts[i])
    free(counts)
    return np_counts
