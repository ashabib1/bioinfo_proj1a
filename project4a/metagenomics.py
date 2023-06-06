def parse_reads_file(reads_fn):

    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads...")
            all_reads = []
            for line in rFile:
                if line[0] == '>':
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends[0])
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None
    
def parse_genome_files():

    tot_reads = []
    for val in range(1000):
        new_str = 'project4a-data/project4a_10000_genome_' + str(val) + '.fasta'
        with open(new_str, 'r') as gFile:
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        tot_reads.append(ref_genome)
    return tot_reads

def open_file(name, outputs):

    with open(name, 'w') as f:
        for output in outputs:
            f.write(output)
            if output != outputs[-1]:
                f.write('\n')
        f.close()

def comparison(a, b):

    count = 0
    for x, y in zip(a,b):
        if x == y:
            count += 1
    return count / len(a)
    
if __name__ == "__main__":

    reads = parse_reads_file('project4a-data/project4a_10000_reads.fasta')#[0:100]
    genomes = parse_genome_files()

    i_values = []

    k = 0
    for read in reads:
        i = 0
        max_val = 0
        high_i = 0
        for genome in genomes:
            tenk = True
            for val in range(len(genome) - len(read)):
                temp_genome = genome[val:val+len(read)]
                if temp_genome[0] != read[0]:
                    continue
                if temp_genome[1] != read[1]:
                    continue
                if temp_genome[2] != read[2]:
                    continue
                if temp_genome[3] != read[3]:
                    continue
                if temp_genome[0:20] == read[0:20]:
                    high_i = i
                    max_val = -0.1
                    tenk = False
                    break
                if temp_genome[30:50] == read[30:50]:
                    high_i = i
                    max_val = -0.2
                    tenk = False
                    break
                comp_val = comparison(temp_genome, read)
                if comp_val > 0.9:
                    high_i = i
                    max_val = comp_val
                    tenk = False
                    break
                if comp_val > max_val:
                    high_i = i
                    max_val = comp_val
            i += 1
            if tenk == False:
                break
        #print(max_val)
        k += 1
        i_values.append(high_i)
        if k%10 == 0:
            print(k)


    # for val in range(100,20000):
    #     i_values.append(50)
    
    final_output = []
    j = 0
    for val in i_values:
        final_output.append('>read_' + str(j) + '   ' + 'Genome_Number_' + str(i_values[j]))
        j += 1

    open_file('read.txt', final_output)