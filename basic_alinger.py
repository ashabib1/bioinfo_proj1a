
def parse_reads_file(reads_fn):

    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing...")
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None
    
def parse_ref_file(ref_fn):

    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing...")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None
    
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
    return count


def find_error_index(a, b):

    i = 0
    indices = []
    for x, y in zip(a, b):
        if x != y:
            indices.append(i)
        i += 1
    return indices

    
if __name__ ==  "__main__":

    THRESHOLD = 0.97

    reads = parse_reads_file('project1a_10000_with_error_paired_reads.fasta')
    reference = parse_ref_file('project1a_10000_reference_genome.fasta')
    max_match_list = []
    outputs = []

    for j in range(len(reads)):
        if reads[j][0][0] == '>':
            continue
        current_read = reads[j][0]
        match_list = []
        for i in range(0, len(reference)-len(current_read)):
            score = comparison(current_read, reference[i:i+len(current_read)])
            match_list.append(score)
            if score / len(current_read) > THRESHOLD:
                err_index = find_error_index(current_read, reference[i:i+len(current_read)])
                try:
                    if err_index[0] > len(current_read) - 5:
                        continue
                except:
                    pass
                try:
                    if err_index[len(err_index) - 1] < 5:
                        continue
                except:
                    pass
                for index in err_index:
                    potential_output = '>S' + str(match_list.index(max(match_list)) + index) + ' ' + str(reference[i + index]) + ' ' + str (current_read[index])
                    if potential_output not in outputs:
                        outputs.append(potential_output)
                    #print(potential_output)
        max_match_list.append(max(match_list))
    outputs.sort(key=lambda x:int(x.split(' ')[0][2:]))
    
    open_file('readme2.txt', outputs)