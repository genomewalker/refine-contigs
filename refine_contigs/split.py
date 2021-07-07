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
    get_components,
    fasta_to_dataframe,
    fast_flatten,
    df_to_seq,
    dereplicate_fragments,
    combine_fragment_files,
    get_components_par,
    get_graph,
    process_alns_reg,
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

def split_contigs(args):

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
        G_components = get_components(G)
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
        logging.info("Finding overlaps between contigs")
        parms = {
            "contigs": str(contigs_tmp),
            "results": results,
            "min_id": min_id,
            "min_cov": args.min_cov,
            "contigs_len": contigs_len,
            "threads": args.threads,
        }

        aln_reg = get_components_par(
            parms=parms,
            components=G_components,
            threads=args.threads,
            func=process_alns_reg,
        )
        miss_contigs = aln_reg["Chromosome"].unique()
        miss_contigs_ovl = aln_reg["Class"].tolist().count("overlap")
        miss_contigs_novl = aln_reg["Class"].tolist().count("non-overlap")
        miss_contigs_ovl_nt = sum(
            aln_reg[aln_reg["Class"] == "overlap"]["length"].tolist()
        )
        miss_contigs_nt = sum(
            [contigs_len[k]["length"] for k in miss_contigs if k in contigs_len]
        )
        miss_contigs_ovl_nt_prop = 100 * (miss_contigs_ovl_nt / assm_len)
        logging.info(
            f"Found {miss_contigs_ovl} overlaps in {len(miss_contigs)} contigs ({miss_contigs_ovl_nt_prop:.2f}% of the assembly)"
        )

        aln_reg_fname = f"{args.output}.split.overlaps.tsv.gz"
        logging.info(f"Saving overlaps to {aln_reg_fname}")

        contigs_names = contigs[["name", "name_original"]]
        contigs_names.columns = ["Chromosome", "name_original"]
        pd.merge(
            left=aln_reg[["Chromosome", "Start", "End", "Class", "length"]],
            right=contigs_names,
            how="inner",
        ).to_csv(aln_reg_fname, sep="\t", compression="gzip", index=False)

        logging.info(
            f"Clustering fragments longer than {args.frag_min_len} NTs [id:{args.frag_cls_id*100}%; cov:{args.frag_cls_cov*100}]"
        )

        aln_reg["frag"] = aln_reg.groupby(["Chromosome"]).cumcount()
        aln_reg["name"] = aln_reg["Chromosome"] + str("-") + aln_reg["frag"].astype(str)

        derep_frag, derep_tsv_frag = dereplicate_fragments(
            frags=aln_reg[(aln_reg["length"] > args.frag_min_len)],
            threads=args.threads,
            tmp_dir=tmp_dir,
            cls_id=args.frag_cls_id,
            cls_cov=args.frag_cls_cov,
            cls_step="fragment",
        )

        logging.info("Combining merged and non-merged contigs")
        dfs = combine_fragment_files(df1=derep_frag, df2=contigs, ids=ids_overlaps)

        logging.info(
            f"Global clustering [id:{args.global_cls_id*100}%; cov:{args.global_cls_cov*100}]"
        )

        derep_global, derep_tsv_global = dereplicate_fragments(
            frags=dfs,
            threads=args.threads,
            tmp_dir=tmp_dir,
            cls_id=args.global_cls_id,
            cls_cov=args.global_cls_cov,
            cls_step="global",
        )

        dfs = fasta_to_dataframe(derep_global)
        dfs["old_name"] = dfs["name"]
        dfs["name"] = dfs.index
        dfs["name"] = dfs["name"].apply(lambda x: f"{prefix}_sp_{x:012d}", 1)
        seq_records = df_to_seq(dfs)

        fname = f"{args.output}.split.fasta.gz"
        logging.info(f"Saving contigs to {fname} file")
        with gzip.open(fname, "wt") as handle:
            SeqIO.write(seq_records, handle, "fasta")

        cls_frag = pd.read_csv(derep_tsv_frag, sep="\t", names=["rep_fragment", "name"])
        cls_global = pd.read_csv(
            derep_tsv_global, sep="\t", names=["rep_global", "name"]
        )
        mappings = (
            pd.merge(
                left=contigs_names,
                right=aln_reg[["Chromosome", "name", "length"]],
                how="left",
            )
            .merge(
                right=cls_frag,
                how="left",
            )
            .merge(
                right=cls_global,
                how="left",
            )
            .merge(right=comps, how="left")
        )
        mappings["old_name"] = mappings.apply(
            lambda x: x["Chromosome"]
            if pd.isnull(x["rep_global"])
            else x["rep_global"],
            axis=1,
        )

        dfs = dfs[["name", "old_name"]]
        dfs.columns = ["final_name", "old_name"]
        mappings = pd.merge(left=mappings, right=dfs, how="left",)[
            [
                "name_original",
                "Chromosome",
                "name",
                "component",
                "rep_fragment",
                "rep_global",
                "final_name",
            ]
        ]
        mappings.columns = [
            "contig_name_original",
            "contig_renamed",
            "contig_name_fragment",
            "component",
            "cls_rep_frag",
            "cls_rep_global",
            "contig_name_refined",
        ]

        mapping_fname = f"{args.output}.split.mapping.tsv.gz"
        logging.info(f"Saving name mappings to {mapping_fname} file")
        mappings.to_csv(mapping_fname, sep="\t", compression="gzip", index=False)

        comp_fname = f"{args.output}.split.all-vs-all.tsv.gz"
        logging.info(f"Saving all-vs-all comparison to {comp_fname} file")
        results.to_csv(comp_fname, sep="\t", compression="gzip", index=False)

        g_fname = f"{args.output}.split.graph-edgelist.tsv.gz"
        logging.info(f"Saving graph edgelist to {g_fname} file")
        nx.to_pandas_edgelist(G).to_csv(
            g_fname, sep="\t", compression="gzip", index=False
        )
        clean_up(keep=args.keep_files, temp_dir=str(tmp_dir))
    else:
        logging.info("Couldn't find any overlaps")
