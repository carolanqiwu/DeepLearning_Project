#!/usr/bin/env python

'''Script for processing scATAC fragment files.

Usage example: 
screen

python chromafold/scripts/scATAC_preparation.py \
--cell_type_prefix cell_type_prefix \
--fragment_file /data/merged_fragments.tsv.gz \
--barcode_file /data/archr_data/archr_filtered_barcode.csv \
--lsi_file /data/archr_data/archr_filtered_lsi.csv \
--genome_assembly "mm10" \
--save_path /data/atac_data 

'''

import random
import sys
import numpy as np
import pysam
import pandas as pd
import os, argparse
import itertools
import scipy
from scipy.sparse import csr_matrix
from scipy import sparse
import pickle
from scipy.sparse import coo_matrix, csr_matrix, find
from sklearn.neighbors import NearestNeighbors

parser = argparse.ArgumentParser(description="Set-up data preparations",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--cell_type_prefix", help="Cell type and its prefix for the file names.",
                    type=str, default="")
parser.add_argument("--fragment_file", help="Path and name of the input fragment file.",
                    type=str, default="")
parser.add_argument("--genome_assembly", help="Genome assembly of the scATAC files.",
                    type=str, default="")
parser.add_argument("--lsi_file", help="Path and name of the LSI file from ArchR running.",
                    type=str, default="")
parser.add_argument("--save_path", help="Path to save all the intermediate and final files.",
                    type=str, default="")

args = parser.parse_args()
config = vars(args)
print(config)
CELL_TYPE = config['cell_type_prefix']
FRAG_FILE = config['fragment_file']
GENOME_ASSEMBLY = config['genome_assembly']
LSI_PATH = config['lsi_file']
SAVE_PATH = config['save_path']

#########################
#     Load functions    #
#########################

def atac_processing(lsi_path, frag_path, genome_assembly, 
                    save_path, cell_type_prefix):
    print("Loading genome sizes")
    chrom_size = pd.read_csv(f'/data1/lesliec/carolw/genome/hg38/{genome_assembly}.chrom.sizes', sep = '\t', header = None, index_col = 0)
    # add chr prefix
    chrom_size.index = ["chr" + str(c) for c in chrom_size.index]
    valid_chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5','chr6','chr7',
                   'chr8', 'chr9','chr10','chr11','chr12','chr13','chr14',
                   'chr15','chr16','chr17','chr18','chr19','chr20',
                   'chr21','chr22','chrX']
    valid_chrom = [i for i in valid_chrom if i in chrom_size.index.to_list()]
    
    print("Loading valid barcodes...")
    lsi = pd.read_csv(lsi_path, index_col=0)
    lsi_barcodes = [x.split("#")[1] for x in lsi.index]

    barcode_list = lsi_barcodes
    valid_barcode = set(barcode_list)
    barcode_list = sorted(lsi_barcodes)
    barcode_to_id = {b: i for i, b in enumerate(barcode_list)}
    n_cells = len(barcode_list)
    print(f"{n_cells} barcodes loaded")

    # Initialize empty CSR matrices per chromosome
    tile500_dict = {}
    tile50_dict = {}
    pbulk50 = {}

    for chrom in valid_chrom:
        n_tiles500 = chrom_size.loc[chrom, 1] // 500 + 1
        n_tiles50 = chrom_size.loc[chrom, 1] // 50 + 1
        tile500_dict[chrom] = csr_matrix((n_cells, n_tiles500), dtype=np.float32)
        tile50_dict[chrom] = csr_matrix((n_cells, n_tiles50), dtype=np.float32)
        pbulk50[chrom] = np.zeros((n_tiles50, 1), dtype=np.int32)

    dtypes = {0:"category", 1:"int32", 2:"int32", 3:"category", 4:"int16"}

    reader = pd.read_csv(frag_path, sep="\t", header=None, dtype=dtypes,
                         compression="gzip", chunksize=500_000)

    print("Streaming fragments and updating matrices...")
    for i, chunk in enumerate(reader):
        chunk = chunk[chunk[3].isin(valid_barcode)]
        if chunk.empty:
            continue

        cell_ids = chunk[3].map(barcode_to_id).values
        chroms = chunk[0].astype(str)

        for pos_col in (1, 2):  # start and end
            pos500 = chunk[pos_col].values // 500
            pos50 = chunk[pos_col].values // 50

            for chrom in np.unique(chroms):
                if chrom not in valid_chrom:
                    continue

                mask = chroms == chrom
                ids_masked = cell_ids[mask]

                # 500bp per-cell
                rows500 = ids_masked.astype(np.int32)
                cols500 = pos500[mask].astype(np.int32)
                data500 = np.ones_like(rows500, dtype=np.float32)
                tile500_dict[chrom] += csr_matrix((data500, (rows500, cols500)),
                                                  shape=tile500_dict[chrom].shape)

                # 50bp per-cell
                rows50 = ids_masked.astype(np.int32)
                cols50 = pos50[mask].astype(np.int32)
                data50 = np.ones_like(rows50, dtype=np.float32)
                tile50_dict[chrom] += csr_matrix((data50, (rows50, cols50)),
                                                 shape=tile50_dict[chrom].shape)

                # pseudo-bulk 50bp
                np.add.at(pbulk50[chrom][:, 0], pos50[mask], 1)

        if i % 10 == 0:
            print(f"Processed chunk {i}")
        
    
    print("Saving results...")
    # 500bp per-cell tiles
    pickle.dump(tile500_dict, open(f"{save_path}/atac/{cell_type_prefix}_tile_500bp_dict.p", "wb"))
    np.save(f"{save_path}/atac/{cell_type_prefix}_tile_500bp_barcode.npy", np.array(barcode_list))

    # 50bp per-cell tiles
    pickle.dump(tile50_dict, open(f"{save_path}/atac/{cell_type_prefix}_tile_50bp_dict.p", "wb"))
    np.save(f"{save_path}/atac/{cell_type_prefix}_tile_50bp_barcode.npy", np.array(barcode_list))

    # 50bp pseudo-bulk
    pickle.dump(pbulk50, open(f"{save_path}/atac/{cell_type_prefix}_tile_pbulk_50bp_dict.p", "wb"))

    print("ATAC preprocessing complete.")
    
    return None
  

def generate_cicero_metacell_sparse(graph, max_overlap, sampled_id=[0]):
    order = np.arange(graph.shape[0])
    np.random.seed(10)
    np.random.shuffle(order)
    selected = np.zeros(graph.shape[0], dtype=bool)
    selected[sampled_id] = True

    for idx in order:
        candidate = graph[idx]  # CSR row
        if selected.any():
            selected_cells = graph[selected]  # subset of rows
            overlap = candidate.dot(selected_cells.T).max()
        else:
            overlap = 0
        if overlap < max_overlap:
            selected[idx] = True
    return selected
    
def metacell_computing(lsi_path, n_neighbors = 100, max_overlap = 33,
                      save_path = SAVE_PATH, cell_type_prefix = CELL_TYPE):
    lsi = pd.read_csv(lsi_path, index_col = 0)
    lsi.index = [x.split('#')[1] for x in lsi.index]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean').fit(lsi.values) 
    graph = nbrs.kneighbors_graph(lsi.values)
    selected = generate_cicero_metacell_sparse(graph, max_overlap = max_overlap)
    
    metacell_assignment = graph[selected, :]  # (n_metacells × n_cells)

    pd.DataFrame(
    metacell_assignment.toarray()
).to_csv(f"{save_path}/atac/{cell_type_prefix}_metacell_mask.csv")

    return None

#########################
#     Script running    #
#########################

def main():
    args = parser.parse_args()
    config = vars(args)
    print(config)
    CELL_TYPE = config['cell_type_prefix']
    FRAG_FILE = config['fragment_file']
    GENOME_ASSEMBLY = config['genome_assembly']
    LSI_PATH = config['lsi_file']
    SAVE_PATH = config['save_path']

    atac_processing(lsi_path = LSI_PATH,
                    frag_path = FRAG_FILE, 
                    genome_assembly=GENOME_ASSEMBLY,
                    save_path = SAVE_PATH, 
                    cell_type_prefix = CELL_TYPE)
    metacell_computing(lsi_path = LSI_PATH, 
                       save_path = SAVE_PATH, 
                       cell_type_prefix = CELL_TYPE)
    
if __name__ == '__main__':
    main()
