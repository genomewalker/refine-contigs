import argparse
import sys
import shutil
import gzip
import os
import pathlib

import logging
import pandas as pd
from multiprocessing import Pool
from functools import partial
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
import tqdm
from refine_contigs import __version__
import time
import networkx as nx
from Bio import SeqIO, Seq, SeqRecord
import pathlib
import uuid
import subprocess
from itertools import chain
from networkx.algorithms.approximation import clique
from statistics import mean
import pyranges as pr

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return value


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        if os.path.isfile(arg):
            parser.error("argument %s: The file %s does not exist!" % (var, arg))
        else:
            parser.error("argument %s: The directory %s does not exist!" % (var, arg))
    else:
        if os.path.isfile(arg):
            return get_open_func(arg)(arg, "rt")  # return an open file handle
        else:
            return arg


help_msg = {
    "search_results": "MMseqs2 search results",
    "contigs": "Contig file to check for misassemblies",
    "min_id": "Minimun id to use for the overlap",
    "min_cov": "Minimun percentage of the coverage for the overlap",
    "frag_min_len": "Minimum fragment length to keep",
    "frag_cls_id": "Minimum identity to cluster the fragments",
    "frag_cls_cov": "Minimum coverage to cluster the fragments",
    "global_cls_id": "Minimum identity to cluster the refined dataset",
    "global_cls_cov": "Minimum coverage to cluster the refined dataset",
    "prefix": "Prefix for contigs name",
    "output": "Fasta file name to save the merged contigs",
    "threads": "Number of threads",
    "debug": "Print debug messages",
    "version": "Print program version",
    "tmp": "Temporary directory",
    "keep_files": "Keep temporary data",
}


def get_arguments(argv=None):

    parser = argparse.ArgumentParser(
        description="Finds misassemblies in ancient data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    optional = parser._action_groups.pop()  # Edited this line
    required = parser.add_argument_group("required arguments")
    ovl_args = parser.add_argument_group("Overlap identification arguments")
    frag_args = parser.add_argument_group("Fragment refinement arguments")
    cls_args = parser.add_argument_group("Final clustering arguments")
    parser._action_groups.append(optional)  # added this line

    # required.add_argument(
    #     "--search-results",
    #     dest="search_results",
    #     default=argparse.SUPPRESS,
    #     required=True,
    #     metavar="FILE",
    #     type=lambda x: is_valid_file(parser, x, "--in-dir"),
    #     help=help_msg["search_results"],
    # )
    required.add_argument(
        "--contigs",
        required=True,
        metavar="FILE",
        default=argparse.SUPPRESS,
        type=lambda x: is_valid_file(parser, x, "--contigs"),
        dest="contigs",
        help=help_msg["contigs"],
    )
    optional.add_argument(
        "--tmp",
        type=str,
        default="./tmp",
        metavar="DIR",
        dest="tmp_dir",
        help=help_msg["tmp"],
    )
    optional.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        dest="threads",
        default=16,
        help=help_msg["threads"],
    )
    optional.add_argument(
        "--prefix",
        type=str,
        default="contig",
        metavar="PREFIX",
        dest="prefix",
        help=help_msg["prefix"],
    )
    optional.add_argument(
        "--output",
        type=str,
        default="contigs-merged",
        metavar="OUT",
        dest="output",
        help=help_msg["output"],
    )
    ovl_args.add_argument(
        "--min-id",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=100, parser=parser, var="--min-id"
        ),
        dest="min_id",
        default=95.0,
        help=help_msg["min_id"],
    )
    ovl_args.add_argument(
        "--min-cov",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--min-cov"
        ),
        default=0.25,
        help=help_msg["min_cov"],
        dest="min_cov",
    )
    frag_args.add_argument(
        "--frag-min-len",
        metavar="INT",
        type=lambda x: check_values(
            x, minval=0, maxval=1e12, parser=parser, var="--frag-min-len"
        ),
        default=1000,
        help=help_msg["frag_min_len"],
        dest="frag_min_len",
    )
    frag_args.add_argument(
        "--frag-cls-id",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--frag-cls-id"
        ),
        default=0.95,
        help=help_msg["frag_cls_id"],
        dest="frag_cls_id",
    )
    frag_args.add_argument(
        "--frag-cls-cov",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--frag-cls-cov"
        ),
        default=0.6,
        help=help_msg["frag_cls_cov"],
        dest="frag_cls_cov",
    )
    cls_args.add_argument(
        "--glob-cls-id",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--global-cls-id"
        ),
        default=0.90,
        help=help_msg["global_cls_id"],
        dest="global_cls_id",
    )
    cls_args.add_argument(
        "--glob-cls-cov",
        metavar="FLOAT",
        type=lambda x: check_values(
            x, minval=0, maxval=1, parser=parser, var="--global-cls-cov"
        ),
        default=0.9,
        help=help_msg["global_cls_cov"],
        dest="global_cls_cov",
    )
    optional.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    optional.add_argument(
        "--keep-files",
        dest="keep_files",
        action="store_false",
        help=help_msg["keep_files"],
    )
    optional.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )

    args = parser.parse_args()

    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def applyParallel(dfGrouped, func, threads, parms):
    p = Pool(threads)
    func = partial(func, parms=parms)
    ret_list = tqdm.tqdm(
        p.map(func, [group for name, group in dfGrouped]),
        total=len([group for name, group in dfGrouped]),
    )
    return pd.concat(ret_list)


def create_graph(results, min_id, min_cov):
    log.debug("Filtering graph with min-id {} and min-cov {}".format(min_id, min_cov))
    results_filt = results[results["qcov"] >= min_cov][
        ["source", "target", "pident", "qcov"]
    ]
    results_filt.columns = ["source", "target", "weight", "qcov"]
    # print(results_filt)
    # results_filt["weight"] = 1
    G = nx.Graph()
    M = nx.from_pandas_edgelist(
        results_filt,
        edge_attr=["weight"],
        create_using=nx.MultiGraph(),
    )
    for u, v, data in M.edges(data=True):
        if not G.has_edge(u, v):
            # set weight to 1 if no weight is given for edge in M
            weight_raw = mean(
                d.get("weight_raw", 0) for d in M.get_edge_data(u, v).values()
            )
            weight = mean(d.get("weight", 0) for d in M.get_edge_data(u, v).values())
            G.add_edge(u, v, weight=weight, weight_raw=weight_raw)

    G.remove_edges_from(list(nx.selfloop_edges(G)))
    # Identify isolated nodes
    edges2rm = [
        (u, v) for u, v, d in G.edges(data=True) if float(d["weight"]) < float(min_id)
    ]
    G.remove_edges_from(edges2rm)
    isolated = list(nx.isolates(G))
    G.remove_nodes_from(isolated)
    return G.to_undirected()


def get_components_clique(G):
    components = []
    if G.number_of_edges() > 2:
        # TODO: Check how many components do we have
        n_comp = nx.number_connected_components(G)
        if n_comp >= 1:
            log.debug("Graph with {} component(s".format(n_comp))
            for component in sorted(nx.connected_components(G), key=len, reverse=True):
                component_sg = G.subgraph(component)
                if component_sg.number_of_nodes() > 1:
                    clq = list(clique.max_clique(component_sg))
                    component_sg_clq = component_sg.subgraph(clq)
                    components.append(component_sg_clq.copy())
        else:
            log.debug("Skipping getting nodes in component")
            components = [None]
    elif G.number_of_edges() == 1:
        components.append(G.copy())
    else:
        log.debug("Skipping getting nodes in component")
        components = [None]
    # partition = community_louvain.best_partition(G_filt, resolution=1.0)
    return components


def get_components(G):
    components = []
    if G.number_of_edges() > 2:
        # TODO: Check how many components do we have
        n_comp = nx.number_connected_components(G)
        if n_comp >= 1:
            log.debug("Graph with {} component(s".format(n_comp))
            for component in sorted(nx.connected_components(G), key=len, reverse=True):
                component_sg = G.subgraph(component)
                if component_sg.number_of_nodes() > 1:
                    components.append(component_sg.copy())
        else:
            log.debug("Skipping getting nodes in component")
            components = [None]
    elif G.number_of_edges() == 1:
        components.append(G.copy())
    else:
        log.debug("Skipping getting nodes in component")
        components = [None]
    # partition = community_louvain.best_partition(G_filt, resolution=1.0)
    return components


def get_strand(start, end):
    if start < end:
        return str("+")
    else:
        return str("-")


def get_starts(start, end):
    if start < end:
        Start = start
        End = end
        Strand = get_strand(start, end)
    else:
        Start = end
        End = start
        Strand = get_strand(start, end)
    return (Start, End, Strand)


def get_coords(x, ref):
    if x["source"] == ref:
        x["Chromosome"] = x["source"]
        x["Start"], x["End"], x["Strand"] = get_starts(x["qstart"], x["qend"])
    else:
        x["Chromosome"] = x["target"]
        x["Start"], x["End"], x["Strand"] = get_starts(x["tstart"], x["tend"])
    return x


def get_best_aln(alns, ref, contigs, contigs_len):
    alns = alns[(alns["source"] == ref) | (alns["target"] == ref)][
        ["source", "target", "qstart", "qend", "tstart", "tend"]
    ].copy()
    aln = alns.apply(get_coords, ref=ref, axis=1)[
        ["Chromosome", "Start", "End", "Strand"]
    ]

    aln = pr.PyRanges(aln).merge(count=False, strand=False)

    ref_len = contigs_len[ref]["length"]
    df = pr.from_dict({"Chromosome": ref, "Start": [1], "End": [ref_len]})

    aln1 = df.subtract(aln)
    aln1.Class = "non-overlap"
    aln.Class = "overlap"
    aln2 = pr.concat([aln, aln1])
    aln2.length = aln2.lengths()
    seqs = pr.get_fasta(aln2, contigs)
    aln2.sequence = seqs
    return aln2.df


def process_alns(node, parms):
    res = get_best_aln(
        ref=node,
        alns=parms["alns"],
        contigs_len=parms["contigs_len"],
        contigs=parms["contigs"],
    )
    return res


def aligned_regions(
    contigs, results, component, min_cov, min_id, contigs_len, threads=1
):
    ids = list(component.nodes())
    results_filt = results[
        (results["source"].isin(ids))
        & (results["target"].isin(ids))
        & (results["qcov"] >= min_cov)
        & (results["pident"] >= min_id)
        & (results["source"] != results["target"])
    ][
        [
            "source",
            "target",
            "pident",
            "alnlen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "qlen",
            "tlen",
        ]
    ].copy()

    parms = {"alns": results_filt, "contigs": contigs, "contigs_len": contigs_len}

    dfs = list(map(partial(process_alns, parms=parms), ids))

    dfs = concat_df(dfs)
    return dfs


def fasta_to_dataframe(infile, header_sep=None, key="name", seqkey="sequence"):
    """Get fasta proteins into dataframe"""
    recs = SeqIO.parse(infile, "fasta")
    keys = [key, seqkey]
    data = [(r.name, str(r.seq)) for r in recs]
    df = pd.DataFrame(data, columns=(keys))
    # fix bad names
    if header_sep not in ["", None]:
        df[key] = df[key].apply(lambda x: x.split(header_sep)[0], 1)
    # df[key] = df[key].str.replace('|','_')
    return df


def df_to_seq(df):
    seq_records = []
    for i, row in df.iterrows():
        record = SeqRecord.SeqRecord(
            seq=Seq.Seq(row["sequence"]), id=row["name"], description=""
        )
        seq_records.append(record)
    return seq_records


def refine_contigs(
    component,
    contigs,
    keep_files,
    overlap,
    minid,
    maxtrim,
    threads,
    conserr,
    tmp,
):
    ids = component.nodes()
    contigs = fasta_to_dataframe(contigs)
    ctg = contigs[contigs["name"].isin(ids)]
    seq_records = df_to_seq(ctg)

    out_suffix = ".fasta"
    fname = str(uuid.uuid4())
    outfile = pathlib.Path(tmp, fname).with_suffix(out_suffix)

    with open(outfile, "w") as handle:
        SeqIO.write(seq_records, handle, "fasta")

    seq, singletons = run_minimus2(
        fname=fname,
        overlap=overlap,
        minid=minid,
        maxtrim=maxtrim,
        tmp_dir=tmp,
        threads=threads,
        conserr=conserr,
    )
    seqs = fasta_to_dataframe(seq)
    seqs["m_type"] = str("merged")
    sglt = fasta_to_dataframe(singletons)
    sglt["m_type"] = str("singleton")

    if keep_files:
        for p in pathlib.Path(tmp).glob("{}*".format(fname)):
            if os.path.isfile(p) or os.path.islink(p):
                os.remove(p)  # remove the file
            elif os.path.isdir(p):
                shutil.rmtree(p)  # remove dir and all contains
            else:
                raise ValueError("file {} is not a file or dir.".format(p))

    return concat_df([seqs, sglt])


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def run_minimus2(fname, overlap, minid, maxtrim, tmp_dir, threads, conserr):
    seq_file = pathlib.Path(tmp_dir, fname).with_suffix(".fasta")
    amos_file = pathlib.Path(tmp_dir, fname).with_suffix(".amos")
    amos_fasta = amos_file.with_suffix(".amos.fasta")
    amos_sglt = amos_file.with_suffix(".amos.singletons.seq")
    amos_cmd = [
        "toAmos",
        "-s",
        seq_file,
        "-o",
        amos_file.with_suffix(".amos.afg"),
    ]
    subprocess.run(amos_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    minimus2_cmd = [
        "minimus2_mod",
        amos_file,
        "-D OVERLAP=" + str(int(overlap)),
        "-D MINID=" + str(float(minid)),
        "-D MAXTRIM=" + str(int(maxtrim)),
        "-D THREADS=" + str(int(threads)),
        "-D CONSERR=" + str(float(conserr)),
    ]
    subprocess.run(minimus2_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return amos_fasta, amos_sglt


def run_mmseqs2(contigs, tmp_dir, max_seq_len, threads):

    mmseqs_tmp = pathlib.Path(tmp_dir, "mmseqs-search/tmp/")
    if not os.path.isdir(mmseqs_tmp):
        os.makedirs(mmseqs_tmp, exist_ok=True)
    mmseqs_db = pathlib.Path(tmp_dir, "mmseqs-search/contigs").with_suffix(".db")
    mmseqs_res = pathlib.Path(tmp_dir, "mmseqs-search/contigs-results").with_suffix(
        ".tsv"
    )
    logging.info("Creating MMseqs2 DB")
    mmseqs_createdb_cmd = [
        "mmseqs",
        "createdb",
        contigs,
        mmseqs_db,
    ]

    subprocess.run(
        mmseqs_createdb_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    logging.info(f"Computing all-vs-all contig search [Max seq len: {max_seq_len}]")
    mmseqs2_search_cmd = [
        "mmseqs",
        "easy-search",
        contigs,
        str(mmseqs_db),
        str(mmseqs_res),
        str(mmseqs_tmp),
        "--format-output",
        "query,theader,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov",
        "--max-seq-len",
        str(max_seq_len),
        "--search-type",
        "3",
        "--threads",
        str(threads),
    ]
    proc = subprocess.run(
        mmseqs2_search_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return mmseqs_res


def dereplicate_fragments(frags, threads, tmp_dir, cls_step, cls_id=0.9, cls_cov=0.6):

    mmseqs_tmp = pathlib.Path(tmp_dir, "mmseqs-cluster/tmp/")
    if not os.path.isdir(mmseqs_tmp):
        os.makedirs(mmseqs_tmp, exist_ok=True)
    mmseqs_fasta = pathlib.Path(tmp_dir, f"mmseqs-cluster/{cls_step}").with_suffix(
        ".fasta.gz"
    )
    mmseqs_fasta_derep = pathlib.Path(
        tmp_dir, f"mmseqs-cluster/{cls_step}-derep_rep_seq"
    ).with_suffix(".fasta")
    mmseqs_tsv_derep = pathlib.Path(
        tmp_dir, f"mmseqs-cluster/{cls_step}-derep_cluster"
    ).with_suffix(".tsv")

    mmseqs_res = pathlib.Path(tmp_dir, f"mmseqs-cluster/{cls_step}-derep")

    seq_records = df_to_seq(frags)
    with gzip.open(mmseqs_fasta, "wt") as handle:
        SeqIO.write(seq_records, handle, "fasta")

    mmseqs2_cluster_cmd = [
        "mmseqs",
        "easy-cluster",
        str(mmseqs_fasta),
        str(mmseqs_res),
        str(mmseqs_tmp),
        "--diag-score",
        "1",
        "--cov-mode",
        "1",
        "--cluster-mode",
        "2",
        "-c",
        str(cls_cov),
        "--min-seq-id",
        str(cls_id),
        "--threads",
        str(threads),
    ]

    proc = subprocess.run(
        mmseqs2_cluster_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return mmseqs_fasta_derep, mmseqs_tsv_derep


def process_components(component, parms):
    res = refine_contigs(
        component=component,
        contigs=parms["contigs"],
        overlap=parms["overlap"],
        minid=parms["minid"],
        maxtrim=parms["maxtrim"],
        threads=parms["threads"],
        tmp=parms["tmp"],
        keep_files=parms["keep_files"],
        conserr=parms["conserr"],
    )
    return res


def process_alns_reg(component, parms):
    res = aligned_regions(
        component=component,
        contigs=parms["contigs"],
        results=parms["results"],
        min_id=parms["min_id"],
        min_cov=parms["min_cov"],
        contigs_len=parms["contigs_len"],
        threads=parms["threads"],
    )
    return res


def get_graph(contigs, tmp_dir, max_seq_len, threads, min_id, min_cov):
    """
    Run mmseqs2
    """
    mmseqs2_res = run_mmseqs2(
        contigs=contigs,
        tmp_dir=tmp_dir,
        max_seq_len=max_seq_len,
        threads=threads,
    )

    """
    Read file with mmseqs2 results
    """
    results = pd.read_csv(
        mmseqs2_res,
        sep="\t",
    )
    results.columns = [
        "source",
        "target",
        "pident",
        "alnlen",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits",
        "qlen",
        "tlen",
        "qcov",
        "tcov",
    ]

    # aln_seqs = get_alignments(results=results, min_id=args.min_id, min_cov=args.min_cov, contigs=contigs)

    G = create_graph(results=results, min_id=min_id, min_cov=min_cov)
    log.info(
        "Graph properties: nodes={} edges={} density={} components={}".format(
            G.number_of_nodes(),
            G.number_of_edges(),
            f"{nx.density(G):.3f}",
            nx.number_connected_components(G),
        )
    )
    return G, results


def get_components_par(parms, components, threads):
    if is_debug():
        dfs = list(map(partial(process_alns_reg, parms=parms), components))
    else:
        p = Pool(threads)
        if len(components) > 1000:
            c_size = int((len(components) / 100))
        else:
            c_size = int((len(components) / threads))
            if c_size == 0:
                c_size = 1
        dfs = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(process_alns_reg, parms=parms),
                    components,
                    chunksize=c_size,
                ),
                total=len(components),
                leave=False,
                ncols=80,
            )
        )
    return concat_df(dfs)


def combine_fragment_files(df1, df2, ids):
    dfs = fasta_to_dataframe(df1)
    # to_include = derep[~derep["name"].isin(ids_components)].copy()
    # to_include["m_type"] = "added"
    # dfs = concat_df([dfs, to_include])
    to_include = df2[~df2["name"].isin(ids)].copy()
    to_include["m_type"] = "added"
    dfs = concat_df([dfs, to_include])
    # dfs["name"] = dfs.index
    # dfs["name"] = dfs["name"].apply(lambda x: f"{name_str}_{x:012d}", 1)
    return dfs
