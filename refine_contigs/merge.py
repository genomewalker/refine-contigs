"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""


import logging
from refine_contigs.utils import (
    fasta_to_dataframe,
    fast_flatten,
    df_to_seq,
    dereplicate_fragments,
    get_components_clique,
    get_components_par,
    get_graph,
    process_minimus2,
    concat_df,
    clean_up
)
import gzip
import networkx as nx
import pathlib, os, sys
from Bio import SeqIO
import uuid
import pandas as pd
import pyfaidx

log = logging.getLogger("my_logger")

sys.setrecursionlimit(10**6)

def merge_contigs(args):
    
    prefix = args.prefix
    dname = str(uuid.uuid4())
    tmp_dir = pathlib.Path(args.tmp_dir, dname).absolute()
    min_id = float(100.0 * args.min_id)

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    # Read contigs file
    logging.info("Processing contig file")
    contigs = fasta_to_dataframe(args.contigs)
    contigs["name_original"] = contigs["name"]
    contigs["name"] = contigs.reset_index().index
    contigs["name"] = contigs["name"].apply(lambda x: f"{prefix}_{x:012d}", 1)
    contigs["length"] = contigs.sequence.map(len)

    logging.info(f"Read and processed {len(contigs.index)} contigs")
    contigs_tmp = pathlib.Path(tmp_dir, dname).with_suffix(".fasta")
    seq_records = df_to_seq(contigs[["name", "sequence"]])
    with open(contigs_tmp, "w") as handle:
        SeqIO.write(seq_records, handle, "fasta")

    pyfaidx.Faidx(str(contigs_tmp))

    max_seq_len = contigs.length.max()
    assm_len = sum(contigs.length)
    contigs_len = contigs[["name", "length"]].set_index("name").T.to_dict()

    G, results = get_graph(
        contigs=str(contigs_tmp),
        tmp_dir=tmp_dir,
        max_seq_len=max_seq_len,
        threads=args.threads,
        min_id=min_id,
        min_cov=args.min_cov,
    )

    if nx.number_connected_components(G) > 0:
        # Get components
        G_components = get_components_clique(G)
        if G_components.count(None) == len(G_components):
            logging.info("Couldn't find any component")
            exit(0)

        d = {
            name: f"comp-{k}"
            for k, comp in enumerate(list(G_components))
            for name in comp
        }
        comps = (
            pd.DataFrame.from_dict(d, orient="index", columns=["component"])
            .rename_axis("Chromosome")
            .reset_index()
        )
        ids_overlaps = fast_flatten([list(n.nodes()) for n in G_components])
        # component = G_components[0]
        # For each component extrac aligned and non-aligned regions
        logging.info(
            f"Trying to merge fragments with Minimus [id:{args.minimus2_minid}%; ovl:{args.minimus2_overlap}]"
        )

        parms = {
            "contigs": str(contigs_tmp),
            "overlap": args.minimus2_overlap,
            "minid": args.minimus2_minid,
            "maxtrim": args.minimus2_maxtrim,
            "threads": args.minimus2_threads,
            "tmp": tmp_dir,
            "conserr": args.minimus2_conserr,
        }

        mn2_res = get_components_par(
            parms=parms,
            components=G_components,
            threads=args.threads,
            func=process_minimus2,
        )
        # Write sequence ids that have been merged
        #  contigs["name"] = contigs.reset_index().index
        # contigs["name"] = contigs["name"].apply(lambda x: f"{prefix}_{x:012d}", 1)
        # contigs["length"] = contigs.sequence.map(len)
        # Write table with merges
        # Combine and cluster
        # rename
        # save
        logging.info("Combining merged and non-merged contigs")
        ids_components = fast_flatten([list(n.nodes()) for n in G_components])
        to_include = contigs[~contigs["name"].isin(ids_components)].copy()
        to_include["m_type"] = "added"
        
        mn2_df = concat_df([mn2_res, to_include])
        mn2_df["old_name"] = mn2_df["name"]
        mn2_df["name"] = mn2_df.index
        mn2_df["name"] = mn2_df["name"].apply(lambda x: f"{prefix}_mn_{x:012d}", 1)
        mn2_df["length"] = mn2_df.sequence.map(len)

        logging.info(
            f"Global clustering [id:{args.global_cls_id*100}%; cov:{args.global_cls_cov*100}]"
        )

        derep_global, derep_tsv_global = dereplicate_fragments(
            frags=mn2_df,
            threads=args.threads,
            tmp_dir=tmp_dir,
            cls_id=args.global_cls_id,
            cls_cov=args.global_cls_cov,
            cls_step="global",
        )

        dfs = fasta_to_dataframe(derep_global)
        #dfs["old_name"] = dfs["name"]
        #dfs["name"] = dfs.index
        #dfs["name"] = dfs["name"].apply(lambda x: f"{prefix}_mn_{x:012d}", 1)

        seq_records = df_to_seq(dfs)

        fname = f"{args.output}.merged.fasta.gz"
        logging.info(f"Saving contigs to {fname} file")
        with gzip.open(fname, "wt") as handle:
            SeqIO.write(seq_records, handle, "fasta")

        cls_global = pd.read_csv(
            derep_tsv_global, sep="\t", names=["rep_global", "name"]
        )

        contigs_names = contigs[["name", "name_original"]]
        contigs_names.columns = ["old_name", "name_original"]
        
        mappings = (
            pd.merge(right=contigs_names,
            left=mn2_df[["name", "old_name", "m_type", "length"]],
            how="left")
            .merge(
                right=cls_global,
                how="left",
            )
        )

        mappings.columns = [
            "contig_name_merged",
            "contig_name_minimus2",
            "contig_minimus2_type",
            "contig_length",
            "contig_name_original",
            "cls_rep_global",
        ]

        mapping_fname = f"{args.output}.merged.mapping.tsv.gz"
        logging.info(f"Saving name mappings to {mapping_fname} file")
        mappings.to_csv(mapping_fname, sep="\t", compression="gzip", index=False)

        comp_fname = f"{args.output}.merged.all-vs-all.tsv.gz"
        logging.info(f"Saving all-vs-all comparison to {comp_fname} file")
        results.to_csv(comp_fname, sep="\t", compression="gzip", index=False)

        g_fname = f"{args.output}.merged.graph-edgelist.tsv.gz"
        logging.info(f"Saving graph edgelist to {g_fname} file")
        nx.to_pandas_edgelist(G).to_csv(
            g_fname, sep="\t", compression="gzip", index=False
        )
        clean_up(keep=args.keep_files, temp_dir=str(tmp_dir))

    else:
        logging.info("Couldn't find any overlaps")
