def parse_reads_file(reads_fn):

    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads...")
            all_reads = []
            for line in rFile:
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None
    
def get_reads(input):

    reads = []
    for i in range(len(input)):
        if input[i][0][0] == '>':
            continue
        reads.append(input[i][0])
    return reads

def get_de_bruijn_graph(kmers):

    adj_list = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix in adj_list:
            adj_list[prefix].append(suffix)
        else:
            adj_list[prefix] = [suffix]

    return adj_list

def get_node_degrees(adjacency_list):
    node_degrees = {}
    for outgoing_node, incoming_nodes in adjacency_list.items():
        # increment the out-degree for outgoing_node by the number of incoming nodes
        if outgoing_node in node_degrees:
            node_degrees[outgoing_node][1] += len(incoming_nodes)
        else:
            node_degrees[outgoing_node] = [0, len(incoming_nodes)]

        # increment the in-degree for each incoming node by 1
        for incoming_node in incoming_nodes:
            if incoming_node in node_degrees:
                node_degrees[incoming_node][0] += 1
            else:
                node_degrees[incoming_node] = [1, 0]

    return node_degrees


def find_paths(adjacency_list):

    # find in-degrees and out-degrees for each node
    node_degrees = get_node_degrees(adjacency_list)

    # find maximal non-branching paths

    # Paths ← empty list
    paths = []

    # for each node v in graph
    for v in list(node_degrees.keys()):
        # if v is not a 1-in-1-out node
        if node_degrees[v] != [1,1]:
            # if out(v) > 0
            if node_degrees[v][1] > 0:
                # for each outgoing edge (v, w) from v
                for w in adjacency_list[v]:
                    # NonBranchingPath ← the path consisting of single edge (v, w)
                    non_branching_path = [v,w]
                    # while w is a 1-in-1-out node
                    while node_degrees[w] == [1,1]:
                        # extend NonBranchingPath by the edge (w, u)
                        u = adjacency_list[w][0]
                        non_branching_path.append(u)
                        # w ← u
                        w = u
                    # add NonBranchingPath to the set Paths
                    paths.append(non_branching_path)

    # for each isolated cycle Cycle in Graph
    #   add Cycle to Paths

    # isolated cycles won't contain any of the nodes traversed in paths, so remove those nodes from our graph
    for path in paths:
        for node in path:
            if node in adjacency_list:
                del adjacency_list[node]

    # remove any nodes that aren't 1-in-1-out nodes
    for node in list(adjacency_list.keys()):
        # if node isn't a 1-in-1-out node
        if node_degrees[node] != [1,1]:
            del adjacency_list[node]

    # adjacency_list now consists of only 1-in-1-out nodes that weren't traversed in any of our other paths
    # while adjacency_list is not empty
    while adjacency_list:
        # select any node as the starting node for our cycle
        start_node = list(adjacency_list.keys())[0]
        
        curr_node = start_node
        next_node = adjacency_list[start_node][0]

        cycle = [start_node]

        first_time_visiting_start_node = True
        while curr_node != start_node or first_time_visiting_start_node:
            first_time_visiting_start_node = False

            # remove curr_node from adjacency_list
            del adjacency_list[curr_node]

            # add next_node to cycle
            cycle.append(next_node)

            # update curr_node and next_node
            curr_node = next_node

            # make sure we haven't reached a dead end
            if next_node not in adjacency_list:
                continue
            next_node = adjacency_list[next_node][0]
        
        # add our cycle to paths
        paths.append(cycle)

    return paths

def open_file(name, outputs):

    with open(name, 'w') as f:
        for output in outputs:
            f.write(output)
            if output != outputs[-1]:
                f.write('\n')
        f.close()


def get_kmer_frequencies(kmers):
    kmer_to_frequency = {}

    for kmer in kmers:
        if kmer in kmer_to_frequency:
            kmer_to_frequency[kmer] += 1
        else:
            kmer_to_frequency[kmer] = 1

    return kmer_to_frequency



if __name__ == "__main__":

    input = parse_reads_file('project3a/project3a_10000_spectrum.fasta')
    reads = get_reads(input)

    k = get_kmer_frequencies(reads)

    k = {kmer: freq for kmer, freq in k.items() if freq > 0}

    l = get_de_bruijn_graph(k)

    m = find_paths(l)
    #print(m)
    for val in m:
        print(len(val))
    # cycle = m[3]
    # cycle2 = m[1]
    j = 0
    for x in m:
        print(j, x[0], x[-1])
        j += 1

    order = [0, 6, 4, 3, 1, 12, 5, 10, 1, 11, 0, 7, 5, 9, 4, 2]
    final = []
    for val in order:
        cycle = m[val]
        for i in range(len(cycle) - 1):
            final.append(reads.index(cycle[i] + cycle[i+1][-1]))
    final_output = []
    for val in final:
        final_output.append('>read_' + str(val))
    print(final_output)

    open_file('read.txt', final_output)
    

