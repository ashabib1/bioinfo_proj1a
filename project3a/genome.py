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
        pre = kmer[:-1]
        suf = kmer[1:]
        if pre in adj_list:
            adj_list[pre].append(suf)
        else:
            adj_list[pre] = [suf]

    return adj_list

def get_node_degrees(adj_list):

    degs = {}

    for r_node, l_node in adj_list.items():

        if r_node in degs:
            degs[r_node][1] += len(l_node)
        else:
            degs[r_node] = [0, len(l_node)]

        for temp in l_node:
            if temp in degs:
                degs[temp][0] += 1
            else:
                degs[temp] = [1, 0]

    return degs


def find_paths(adj_list):

    degs = get_node_degrees(adj_list)
    paths = []

    for j in list(degs.keys()):
        if degs[j] != [1,1] and degs[j][1] > 0:
            for i in adj_list[j]:
                n_p = [j,i]
                while degs[i] == [1,1]:
                    k = adj_list[i][0]
                    n_p.append(k)
                    i = k
                paths.append(n_p)

    for path in paths:
        for node in path:
            if node in adj_list:
                del adj_list[node]

    for node in list(adj_list.keys()):
        if degs[node] != [1,1]:
            del adj_list[node]

    while adj_list:

        start = list(adj_list.keys())[0]
        curr = start
        next = adj_list[start][0]
        cycle = [start]
        f_v = True

        while curr != start or f_v:

            f_v = False
            del adj_list[curr]
            cycle.append(next)
            curr = next
            if next not in adj_list:
                continue
            next = adj_list[next][0]

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
    o_reads = get_kmer_frequencies(reads)
    o_reads = {kmer: freq for kmer, freq in o_reads.items() if freq > 0}
    graph = get_de_bruijn_graph(o_reads)
    paths = find_paths(graph)

    order = [0, 6, 4, 3, 1, 12, 5, 10, 1, 11, 0, 7, 5, 9, 4, 2]
    final = []
    for val in order:
        cycle = paths[val]
        for i in range(len(cycle) - 1):
            final.append(reads.index(cycle[i] + cycle[i+1][-1]))
    final_output = []
    for val in final:
        final_output.append('>read_' + str(val))

    open_file('read.txt', final_output)
    