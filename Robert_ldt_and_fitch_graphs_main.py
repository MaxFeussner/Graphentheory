from pathlib import Path
from random import sample

import asymmetree.hgt as hgt
import asymmetree.treeevolve as te
import networkx as nx
import pandas as pd
from asymmetree.cograph import Cotree
from asymmetree.datastructures import PhyloTree

import asymmetree.tools.GraphTools as gt

from cluster_deletion import calculate_cl_fitch
from local_stuff import SERIALIZE_PATH

SERIALIZE_PATH = Path(SERIALIZE_PATH)
TREE_PATH = SERIALIZE_PATH / "tree_pickles"
SUBGRAPH_FRACTIONS = [1.0, 0.8, 0.6, 0.4, 0.2]

def calculate_recall(TP, FN):
    if TP + FN != 0:
        RECALL = (TP / (TP + FN))
        return RECALL
    else:
        return 0

def calculate_precision(TP, FP):
    if TP + FP != 0:
        PRECISION = (TP / (TP + FP))
        return PRECISION
    else:
        return 0

def calculate_accuracy(TP, TN, FP, FN):
    nenner = (TP + TN + FP + FN)
    if nenner != 0:
        ACCURACY = ((TP + TN) / nenner)
        return ACCURACY
    else:
        return 0

def compare_fitches(true_graph, other_graph):
    # 'true_graph ' and ' other_graph ' are NetworkX Graph
    # other_graph ~ cl-Fitch oder rs-Fitch
    ct = gt.contingency_table(true_graph, other_graph)
    # Ausgabe Anzahl von True-Positives, True-Negatives, False-Positives und False-Negatives
    tfpn = (ct['tp'], ct['tn'], ct['fp'], ct['fn'])
    return tfpn


def sample_nodes(G, fraction):
    nodes = G.nodes
    nodes = sample(nodes, round(len(nodes) * fraction))
    return nodes


def calculate_rs_fitch(LDT):
    # 'ldt_subgraph' is an LDT graph extracted from a species / gene tree pair
    constructor = hgt.RsScenarioConstructor(LDT)
    if constructor.run():  # construct and check success
        S2 = constructor.S  # species tree
        T2 = constructor.T  # gene tree ( does not contain losses )

        transfer_edges2 = hgt.rs_transfer_edges(T2, S2)
        recognition_fitch = hgt.undirected_fitch(T2, transfer_edges2)
        return recognition_fitch


def get_genetree_triples_from_ldt(ldt, real_triples):
    """
    Ein Tripel ist so definiert, dass jeder induzierter Teilgraph aus 3 Knoten x,y,z ein Tripel xy|z ist, gdw.
    xy eine Kante ist und weder xz noch yz Kanten sind. xy|z == yx|z.
    :param ldt:
    :param real_triples: Sorted real triples
    :return:
    """
    if not ldt.edges:
        return []
    edges = ldt.edges
    nodes = ldt.nodes
    triples = []
    for edge in edges:
        e1, e2 = edge[0], edge[1]
        other_nodes = [node for node in nodes if node not in (e1, e2)]
        for x in other_nodes:
            if any((ldt.has_edge(e1, x), ldt.has_edge(e2, x))):
                continue
            triple = (*sorted(edge), x)
            assert triple in real_triples
            triples.append(triple)
    return triples if triples else []


def colors_not_unique(node_color_tuples, *args):
    """
    Funktioniert, weil set((1,1,1)) -> {1}. Es werden also nur die unique Werte im Set behalten.
    :param node_color_tuples:
    :param args: any number of nodes from same graph as node_color_tuples
    :return: Boolean: True if colors are not unique.
    """
    colors = []
    for node in args:
        colors.append(node_color_tuples[node])
    return not len(set(colors)) == len(colors)


def get_speciestree_triples_from_ldt(ldt, real_triples):
    """
    Ein Tripel ab|c == ba|c, mit a = color(x), b = color(y) und c = color(z) wird erstellt, gdw.
        - x,y,z sind Knoten von ldt, a,b,c sind paarweise verschieden
        - Es gibt Kanten xz und yz aber nicht xy
    :param ldt:
    :param real_triples: Sorted triples
    :return:
    """
    edges = ldt.edges
    node_color_tuples = ldt.nodes(data='color')
    triples = []
    #  foreach edge find a node that is connected to exaclty one of the edge-nodes
    for edge in edges:
        e1, e2 = edge[0], edge[1]
        if colors_not_unique(node_color_tuples, e1, e2):
            continue
        other_nodes = [node[0] for node in node_color_tuples if node[0] not in (e1, e2)]
        for x in other_nodes:
            if colors_not_unique(node_color_tuples, e1, e2, x):
                continue
            if ldt.has_edge(e1, x) ^ ldt.has_edge(e2, x):  # bool ^ bool entspricht xor
                triple = ()
                if ldt.has_edge(e1, x) and not ldt.has_edge(e2, x):
                    triple = (*sorted((node_color_tuples[e2], node_color_tuples[x])), node_color_tuples[e1])
                elif ldt.has_edge(e2, x) and not ldt.has_edge(e1, x):
                    triple = (*sorted((node_color_tuples[e1], node_color_tuples[x])), node_color_tuples[e2])
                assert triple in real_triples, f"{triple, e1, e2, x}"
                triples.append(triple)
    return triples


def sort_triples(triple_list):
    return [(*sorted((t[0], t[1])), t[2]) for t in triple_list]


def informative_triples(ldt, OGT, S):
    """
    Vergleich von den wahren Tripel Mengen des Simulationsszenarios und den Tripeln, die aus dem LDT Graphen
    extrahiert werden (diese sollen Teilmengen der wahren Tripel sein!).
    Gen und Speziesbäume haben verschiedene Tripelmengen.
    :return: Meaningful quantification??
    """
    ldt = ldt  # G=(L,E)
    real_species_triples = sort_triples(S.get_triples(id_only=True))
    real_gene_triples = sort_triples(OGT.get_triples(id_only=True))
    species_triples_from_ldt = get_speciestree_triples_from_ldt(ldt, real_species_triples)
    gene_triples_from_ldt = get_genetree_triples_from_ldt(ldt, real_gene_triples)

    s = len(species_triples_from_ldt) / len(real_species_triples) if len(real_gene_triples) != 0 else None
    g = len(gene_triples_from_ldt) / len(real_gene_triples) if len(real_gene_triples) != 0 else None
    return {"s_triples_frac": s, "g_triples_frac": g}


def assert_nodes_count(ldt, real_fitch):
    len_ldt = ldt.number_of_nodes()
    len_fitch = real_fitch.number_of_nodes()
    assert len_ldt == len_fitch, f"ldt: {len_ldt}, fitch: {len_fitch}"


if __name__ == '__main__':
    #  iterate through serialized trees
    metadata = pd.read_csv(SERIALIZE_PATH / "metadata.csv")
    row_list = []  # List of dictionaries to create DataFrame at the end
    counter = 0
    for index, row in metadata.iterrows():
        if index % 100 == 0:
            print(f"{index} Dateien verarbeitet")
        s_filemane = row["S_filename"]
        tgt_filename = row["TGT_filename"]
        S = PhyloTree.load(TREE_PATH / s_filemane)
        TGT = PhyloTree.load(TREE_PATH / tgt_filename)
        OGT = te.observable_tree(TGT)
        """
        Later Divergence Time Graph G< where nodes are genes and edges denote that lca(x,y) is later than lca(s(x), s(y)),
        i.e. the two genes x,y diverged later than the species mapped to them.
        """
        ldt = hgt.ldt_graph(OGT, S)
        transfer_edges = hgt.rs_transfer_edges(OGT, S)
        real_fitch = hgt.undirected_fitch(OGT, transfer_edges)

        #  Prüfe ob ldt und fitch gleiche knoten haben
        assert_nodes_count(ldt, real_fitch)
        #  Prüfe ob graph leer ist, wenn ja wird er nicht weiter verarbeitet, das führt zu fehlenden Zeilen in der csv
        if nx.is_empty(ldt):
            continue
        #  Prüfe ob ldt Cograph ist
        assert Cotree.cotree(ldt), "ldt ist kein Cograph"

        spotted_triple_fractions = informative_triples(ldt, OGT, S)


        for subgraph_frac in SUBGRAPH_FRACTIONS:
            """
            Bilde subgraphen und vergleiche die verschiedenen Methoden zur Bestimmung der Xenologe.
            In dieser (inneren) Schleife wir jeweils ein dict erstellt. Gespeichert wird das zunächst in der
            "row_list", das ist eine Liste von dicts, die kann von pandas zu eine DataFrame verwandelt werden, das
            dann als csv gespeichert wird.
            """
            cols_from_metadata = ["S_filename", "s_n_leafs", "gt_dupl_rate", "gt_loss_rate", "gt_hgt_rate",
                                  "s_non_binary_prop"]
            new_data = row[cols_from_metadata].to_dict()

            nodes_for_induced_subgraph = sample_nodes(ldt, subgraph_frac)

            ldt_subgraph = ldt.subgraph(nodes_for_induced_subgraph)
            real_fitch_subgraph = real_fitch.subgraph(nodes_for_induced_subgraph)

            rs_fitch = calculate_rs_fitch(ldt)
            cl_fitch = calculate_cl_fitch(ldt)
            rs_fitch_subgraph = rs_fitch.subgraph(nodes_for_induced_subgraph)
            cl_fitch_subgraph = cl_fitch.subgraph(nodes_for_induced_subgraph)

            # Vergleiche Kanten des wahren Fitch-Graph (und Subgraphen) mit den jeweiligen cl- und rs-Fitches
            tfpn_cl_fitch = compare_fitches(real_fitch_subgraph, cl_fitch_subgraph)
            tfpn_rs_fitch = compare_fitches(real_fitch_subgraph, rs_fitch_subgraph)

            # ldt und real fitch graph
            new_data["subgraph_frac"] = subgraph_frac
            new_data["nodes_ldt"] = ldt_subgraph.number_of_nodes()
            new_data["edges_ldt_subgraph"] = ldt_subgraph.number_of_edges()
            new_data["edges_real_fitch"] = real_fitch_subgraph.number_of_edges()

            #TODO Auch für Subgraphen implementieren
            new_data.update(spotted_triple_fractions)

            new_data["edges_cl_fitch"] = cl_fitch_subgraph.number_of_edges()
            new_data["edges_rs_fitch"] = rs_fitch_subgraph.number_of_edges()

            # results from compare edges of real-fitch vs. cl-fitch
            # (ct['tp'], ct['tn'], ct['fp'], ct['fn'])
            new_data["cl_fitch_tp"] = tfpn_cl_fitch[0]
            new_data["cl_fitch_tn"] = tfpn_cl_fitch[1]
            new_data["cl_fitch_fp"] = tfpn_cl_fitch[2]
            new_data["cl_fitch_fn"] = tfpn_cl_fitch[3]

            cl_RECALL = calculate_recall(tfpn_cl_fitch[0], tfpn_cl_fitch[3]) # (tfpn_cl_fitch[0]/(tfpn_cl_fitch[0] + tfpn_cl_fitch[3]))
            cl_PRECISION = calculate_precision(tfpn_cl_fitch[0], tfpn_cl_fitch[2]) #(tfpn_cl_fitch[0]/(tfpn_cl_fitch[0] + tfpn_cl_fitch[2]))
            cl_ACCURACY = calculate_accuracy(tfpn_cl_fitch[0],tfpn_cl_fitch[1],tfpn_cl_fitch[2],tfpn_cl_fitch[3]) #((tfpn_cl_fitch[0] + tfpn_cl_fitch[1])/(tfpn_cl_fitch[0] + tfpn_cl_fitch[1] + tfpn_cl_fitch[2] + tfpn_cl_fitch[3]))
            new_data["cl_RECALL"] = cl_RECALL
            new_data["cl_PRECISION"] = cl_PRECISION
            new_data["cl_ACCURACY"] = cl_ACCURACY

            # results from compare edges of real-fitch vs. rs-fitch
            # (ct['tp'], ct['tn'], ct['fp'], ct['fn'])
            new_data["rs_fitch_tp"] = tfpn_rs_fitch[0]
            new_data["rs_fitch_tn"] = tfpn_rs_fitch[1]
            new_data["rs_fitch_fp"] = tfpn_rs_fitch[2]
            new_data["rs_fitch_fn"] = tfpn_rs_fitch[3]

            rs_RECALL = calculate_recall(tfpn_rs_fitch[0], tfpn_rs_fitch[3]) # (tfpn_rs_fitch[0]/(tfpn_rs_fitch[0] + tfpn_rs_fitch[3]))
            rs_PRECISION = calculate_precision(tfpn_rs_fitch[0], tfpn_rs_fitch[2]) #(tfpn_rs_fitch[0]/(tfpn_rs_fitch[0] + tfpn_rs_fitch[2]))
            rs_ACCURACY = calculate_accuracy(tfpn_rs_fitch[0],tfpn_rs_fitch[1],tfpn_rs_fitch[2],tfpn_rs_fitch[3]) #((tfpn_rs_fitch[0] + tfpn_rs_fitch[1])/(tfpn_rs_fitch[0] + tfpn_rs_fitch[1] + tfpn_rs_fitch[2] + tfpn_rs_fitch[3]))
            new_data["rs_RECALL"] = rs_RECALL
            new_data["rs_PRECISION"] = rs_PRECISION
            new_data["rs_ACCURACY"] = rs_ACCURACY

            row_list.append(new_data)

    results = pd.DataFrame(row_list)
    results.to_csv(SERIALIZE_PATH / "results.csv")
