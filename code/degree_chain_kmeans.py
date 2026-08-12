import argparse
import os
import re
import warnings

import h5py
import networkx as nx
import numpy as np
import scipy.io as sio
from sklearn.cluster import KMeans
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import davies_bouldin_score
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler

try:
    import community as community_louvain
    HAS_LOUVAIN = True
except ImportError:
    community_louvain = None
    HAS_LOUVAIN = False

GLOBAL_SEED = 915
KMAX = 10
DATA_ROOT = r"F:\YAN1\PDENB-Copy\Motif_subspace\TSSKEA-main\Data"
OUTPUT_ROOT = r"F:\YAN1\Motif"


DATASETS = [
    {
        "cancer": "BRCA",
        "network": "PGIN",
        "input_dir": "BRCA_node_net",
        "input_pattern": "DE_BRCA_sample_result_*.mat",
        "matrix_key": "subnetwork_adjacency",
        "output_dir": "degree_chain_kmeans_cluster1_PGIN_BRCA_result",
        "output_prefix": "DE_BRCA_node_biomarker",
    },
    {
        "cancer": "LUSC",
        "network": "PGIN",
        "input_dir": "LUSC_node_net",
        "input_pattern": "DE_LUNG_sample_result_*.mat",
        "matrix_key": "subnetwork_adjacency",
        "output_dir": "degree_chain_kmeans_cluster1_PGIN_LUSC_result",
        "output_prefix": "DE_LUNG_node_biomarker",
    },
    {
        "cancer": "LUAD",
        "network": "PGIN",
        "input_dir": "LUAD_node_net",
        "input_pattern": "DE_LUNG_sample_result_*.mat",
        "matrix_key": "subnetwork_adjacency",
        "output_dir": "degree_chain_kmeans_cluster1_PGIN_LUAD_result",
        "output_prefix": "DE_LUNG_node_biomarker",
    },
    {
        "cancer": "BRCA",
        "network": "PEN",
        "input_dir": "BRCA_edge_net",
        "input_pattern": "DE_BRCA_edge_network_*.mat",
        "matrix_key": "edge_net",
        "output_dir": "degree_chain_kmeans_cluster1_PEN_BRCA_result",
        "output_prefix": "DE_BRCA_edge_biomarker",
    },
    {
        "cancer": "LUSC",
        "network": "PEN",
        "input_dir": "LUSC_edge_net",
        "input_pattern": "DE_LUSC_edge_network_*.mat",
        "matrix_key": "edge_net",
        "output_dir": "degree_chain_kmeans_cluster1_PEN_LUSC_result",
        "output_prefix": "DE_LUSC_edge_biomarker",
    },
    {
        "cancer": "LUAD",
        "network": "PEN",
        "input_dir": "LUAD_edge_net",
        "input_pattern": "DE_LUAD_edge_network_*.mat",
        "matrix_key": "edge_net",
        "output_dir": "degree_chain_kmeans_cluster1_PEN_LUAD_result",
        "output_prefix": "DE_LUAD_edge_biomarker",
    },
]


def sample_id_from_name(filename):
    match = re.search(r"_(\d+)\.mat$", filename)
    if not match:
        raise ValueError(f"Cannot parse sample id from {filename}")
    return int(match.group(1))


def clean_adjacency(A):
    A_abs = np.abs(np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0))
    A_abs[A_abs < 1e-12] = 0
    return A_abs


def load_matrix(mat_path, matrix_key):
    try:
        data = sio.loadmat(mat_path)
        if matrix_key not in data:
            raise KeyError(matrix_key)
        return np.array(data[matrix_key], dtype=float)
    except NotImplementedError:
        with h5py.File(mat_path, "r") as handle:
            if matrix_key not in handle:
                raise KeyError(matrix_key)
            return np.array(handle[matrix_key], dtype=float).T


def degree_chain_features(A):
    A_abs = clean_adjacency(A)
    degree_vec = np.sum(A_abs, axis=1).reshape(-1, 1)
    diag_a2 = np.einsum("ij,ji->i", A_abs, A_abs).reshape(-1, 1)
    chain_vec = np.abs(np.square(degree_vec) - diag_a2)
    features = np.hstack([degree_vec, chain_vec])
    return np.nan_to_num(features), degree_vec, chain_vec


def compact_labels(labels):
    labels = np.asarray(labels).reshape(-1)
    unique_labels = np.unique(labels)
    label_map = {old: new for new, old in enumerate(unique_labels)}
    return np.array([label_map[label] for label in labels], dtype=np.int32)


def choose_kmeans_labels(features, kmax=KMAX):
    n_samples = features.shape[0]
    if n_samples < 2 or np.allclose(features, features[0]):
        return np.zeros(n_samples, dtype=np.int32), 1, np.nan

    scaled = StandardScaler().fit_transform(features)
    search_k = min(kmax, n_samples - 1)
    best_k = 1
    best_score = np.inf
    best_labels = None

    for k in range(2, search_k + 1):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            labels = KMeans(
                n_clusters=k,
                random_state=GLOBAL_SEED,
                n_init=5,
            ).fit_predict(scaled)

        if len(np.unique(labels)) < 2:
            continue

        score = davies_bouldin_score(scaled, labels)
        if score < best_score:
            best_k = k
            best_score = score
            best_labels = labels

    if best_labels is None:
        return np.zeros(n_samples, dtype=np.int32), 1, np.nan

    final_labels = KMeans(
        n_clusters=best_k,
        random_state=GLOBAL_SEED,
        n_init=10,
    ).fit_predict(scaled)
    return compact_labels(final_labels), best_k, best_score


def build_cluster_adjacency(A, labels):
    A_abs = clean_adjacency(A)
    labels = compact_labels(labels)
    cluster_count = len(np.unique(labels))
    CAG = np.zeros((cluster_count, cluster_count), dtype=float)
    rows, cols = np.nonzero(A_abs)

    for i, j in zip(rows, cols):
        if i >= j:
            continue

        ci, cj = labels[i], labels[j]
        if ci == cj:
            continue

        CAG[ci, cj] += A_abs[i, j]
        CAG[cj, ci] += A_abs[i, j]

    return CAG


def build_super_adjacency(CAG, labels):
    labels = compact_labels(labels)
    super_count = len(np.unique(labels))
    super_cag = np.zeros((super_count, super_count), dtype=float)

    for i in range(CAG.shape[0]):
        for j in range(i + 1, CAG.shape[1]):
            ci, cj = labels[i], labels[j]
            if ci == cj:
                continue
            super_cag[ci, cj] += CAG[i, j]
            super_cag[cj, ci] += CAG[i, j]

    return super_cag


def super_cluster(CAG, use_louvain=True):
    cluster_count = CAG.shape[0]
    if cluster_count <= 3:
        return np.arange(cluster_count, dtype=np.int32)

    if use_louvain and HAS_LOUVAIN:
        graph = nx.from_numpy_array(CAG)
        partition = community_louvain.best_partition(
            graph,
            weight="weight",
            random_state=GLOBAL_SEED,
        )
        current_labels = np.zeros(cluster_count, dtype=np.int32)
        for node, label in partition.items():
            current_labels[node] = label
        current_labels = compact_labels(current_labels)

        if len(np.unique(current_labels)) <= 3:
            return current_labels

        while len(np.unique(current_labels)) > 3:
            unique_labels, sizes = np.unique(current_labels, return_counts=True)
            smallest = unique_labels[np.argmin(sizes)]
            super_cag = build_super_adjacency(CAG, current_labels)

            connections = super_cag[int(smallest), :].copy()
            connections[int(smallest)] = -1
            target = int(np.argmax(connections))

            if connections[target] <= 0:
                largest = unique_labels[np.argmax(sizes)]
                target = int(largest)
                if target == int(smallest):
                    target = int(unique_labels[unique_labels != smallest][0])

            current_labels[current_labels == smallest] = target
            current_labels = compact_labels(current_labels)

        return current_labels

    if np.count_nonzero(CAG) == 0:
        labels = np.arange(cluster_count, dtype=np.int32)
        labels[3:] = 2
        return labels

    spectral = SpectralClustering(
        n_clusters=3,
        affinity="precomputed",
        random_state=GLOBAL_SEED,
    )
    return compact_labels(spectral.fit_predict(CAG))


def process_file(dataset, mat_path):
    sample_id = sample_id_from_name(os.path.basename(mat_path))
    A = load_matrix(mat_path, dataset["matrix_key"])
    features, degree_vec, chain_vec = degree_chain_features(A)
    labels, best_k, best_db_score = choose_kmeans_labels(features)
    CAG = build_cluster_adjacency(A, labels)
    super_label = super_cluster(CAG)
    labels_super = compact_labels(super_label[labels])
    final_k = len(np.unique(labels_super))

    out_dir = os.path.join(OUTPUT_ROOT, dataset["output_dir"])
    os.makedirs(out_dir, exist_ok=True)
    out_name = f"{dataset['output_prefix']}_{sample_id}_Motifcluster_chain.mat"
    out_path = os.path.join(out_dir, out_name)

    sio.savemat(
        out_path,
        {
            "motif_label": labels,
            "motif_label_super": labels_super,
            "super_label": super_label,
            "CAG": CAG,
            "chain_vec": chain_vec,
            "degree_vec": degree_vec,
            "feature_matrix": features,
            "best_k": np.array([[best_k]], dtype=np.int32),
            "final_k": np.array([[final_k]], dtype=np.int32),
            "best_db_score": np.array([[best_db_score]], dtype=float),
            "feature_names": np.array(["degree", "chain_motif"], dtype=object),
        },
    )

    print(
        f"[{dataset['network']} {dataset['cancer']}] sample {sample_id:03d} | "
        f"nodes={A.shape[0]} | K={best_k} -> {final_k}"
    )
    return out_path


def run_dataset(dataset):
    input_dir = os.path.join(DATA_ROOT, dataset["input_dir"])
    files = [
        os.path.join(input_dir, name)
        for name in os.listdir(input_dir)
        if re.fullmatch(dataset["input_pattern"].replace("*", r"\d+"), name)
    ]
    files.sort(key=lambda path: sample_id_from_name(os.path.basename(path)))

    completed = 0
    for mat_path in files:
        process_file(dataset, mat_path)
        completed += 1

    print(
        f"[{dataset['network']} {dataset['cancer']}] completed {completed} files -> "
        f"{os.path.join(OUTPUT_ROOT, dataset['output_dir'])}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="KMeans clustering using degree and chain-motif features."
    )
    parser.add_argument(
        "--network",
        choices=["all", "PGIN", "PEN"],
        default="all",
        help="Network type to run.",
    )
    parser.add_argument(
        "--cancer",
        choices=["all", "BRCA", "LUSC", "LUAD"],
        default="all",
        help="Cancer dataset to run.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    np.random.seed(GLOBAL_SEED)
    args = parse_args()
    selected = [
        dataset
        for dataset in DATASETS
        if (args.network == "all" or dataset["network"] == args.network)
        and (args.cancer == "all" or dataset["cancer"] == args.cancer)
    ]
    for dataset in selected:
        run_dataset(dataset)
