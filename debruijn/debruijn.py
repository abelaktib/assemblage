#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import os
import sys
import statistics
from random import randint
import argparse
import pickle
from operator import itemgetter
import random
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
random.seed(9001)
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as file:
        line = file.readlines()
        for i in range(1, len(line), 4):
            yield line[i].strip()


def cut_kmer(read, kmer_size):
    for i in range(0, len(read)-2, 1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmers_dict = {}
    seq = read_fastq(fastq_file)
    for kmer in seq:
        kmers_cut = cut_kmer(kmer, kmer_size)
        for i in kmers_cut:
            if i in kmers_dict.keys():
                kmers_dict[i] += 1
            else:
                kmers_dict[i] = 1
    return kmers_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for k, i in kmer_dict.items():
        graph.add_edge(k[:-1], k[1:], weight=i)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    start = 1
    i = -1
    if delete_entry_node is True:
        start = 0
    if delete_sink_node is True:
        i = None
    for path in path_list:
        graph.remove_nodes_from(path[start:i]) 
    return graph


def std(data):
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    weight_list = []
    for i, j in enumerate(weight_avg_list):
        if j == max(weight_avg_list) :
            weight_list.append(i)
    if len(weight_list) > 1:
        l = []
        for i in weight_list:
            if path_length[i] == max(path_length):
                l.append(i)
        if len(l) > 1:
            best = random.randint[l]
        else:
            best = l[0]
    else:
        best = weight_list[0]

    tmp=[]
    for path in path_list:
        tmp.append(path)
    tmp.pop(best)

    return remove_paths(graph, tmp, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    moy = 0
    for i in graph.subgraph(path).edges(data=True) :
        moy += i[2]["weight"]
    moy = moy / (len(path)-1)
    return moy  


def solve_bubble(graph, ancestor_node, descendant_node):
    weight_avg_list = []
    path_lengt = []

    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    for path in path_list:
        path_lengt += [len(path)]
        weight_avg_list += [path_average_weight(graph, path)]

    graph = select_best_path(graph, path_list, path_lengt, weight_avg_list) 

    return graph


def simplify_bubbles(graph):
    flag = False

    for node in graph:
        if graph.in_degree(node) >= 2:
            node_pred = []

            for node in graph.predecessors(node):
                node_pred += [node]

            if len(node_pred) <= 2:
                flag = True
                break
        if flag:
            break
    if flag is False:
        ancestor = nx.lowest_common_ancestor(graph, node_pred[0], node_pred[1])
        graph = solve_bubble(graph, ancestor, node)

    return graph
    



def solve_entry_tips(graph, starting_nodes):
    n_list = []
    if len(starting_nodes) > 2:
        for n_i in starting_nodes:
            for n_j in starting_nodes:
                n_list.append( [(n_i, n_j)])
    else:
        n_list = [(starting_nodes)]

    a_list = []

    for node in n_list:
        a_list.append( [nx.lowest_common_ancestor(graph.reverse(),
                                             node[0], node[1])])

    path_list = []
    weight_avg_list = []
    path_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, node[0], a_list[i]))[0]
        path_2 = list(nx.all_simple_paths(graph, node[1], a_list[i]))[0]

        path_list += [path_1, path_2]

        path_len += [len(path_1), len(path_2)]
        weight_avg_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, path_list, path_len, weight_avg_list,
                             delete_entry_node=True)

    return graph




def solve_out_tips(graph, ending_nodes):
    n_list = []

    if len(ending_nodes) > 2:
        for n_i in ending_nodes:
            for n_j in ending_nodes:
                n_list += [(n_i, n_j)]
    else:
        n_list = [tuple(ending_nodes)]

    a_list = []

    for node in n_list:
        a_list += [nx.lowest_common_ancestor(graph, node[0], node[1])]

    path_list = []
    weight_avg_list = []
    path_len = []

    for i, node in enumerate(n_list):
        path_1 = list(nx.all_simple_paths(graph, a_list[i], node[0]))[0]
        path_2 = list(nx.all_simple_paths(graph, a_list[i], node[1]))[0]

        path_list += [path_1, path_2]

        path_len += [len(path_1), len(path_2)]
        weight_avg_list += [path_average_weight(graph, path_1),
                   path_average_weight(graph, path_2)]

    graph = select_best_path(graph, path_list, path_len, weight_avg_list,
                             delete_sink_node=True)

    return graph


def get_starting_nodes(graph):
    nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            nodes.append(node)
    return nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for enter in starting_nodes:
        for sortie in ending_nodes:
            if nx.has_path(graph, enter, sortie):
                paths = nx.all_simple_paths(graph, enter, sortie)
                for x in paths:
                    contig = "".join(x[::len(x[0])])
                    contigs.append((contig, len(contig)))
    return contigs


def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as file:
        for num, long in enumerate(contigs_list):
            file.write(f">contig_{num} len={long[1]}\n")
            file.write(fill(long[0]))
            file.write("\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)



# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer)
    graph = simplify_bubbles(graph)
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)

    graph = solve_entry_tips(graph,start)

    graph = solve_out_tips(graph, end)

    # Getting contigs.
    contigs_list = get_contigs(graph, start, end)
    save_contigs(contigs_list, args.output_file)

    draw_graph(graph, args["graphimg_file"])

    save_graph(graph,args["output_file"])



if __name__ == '__main__':
    main()
