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
    
def parse_ref_file(ref_fn):

    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Reference...")
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

def remove(string, i):

    a = string[:i]
    b = string[i+1:]
    return a + b



    
if __name__ ==  "__main__":

    THRESHOLD = 0.95
    DELETION_THRESHOLD = 0.6

    reads = parse_reads_file('project1a_10000_with_error_paired_reads.fasta')
    reference = parse_ref_file('project1a_10000_reference_genome.fasta')
    max_match_list = []
    outputs = []
    initial = 0
    insertion = False

    for j in range(len(reads)):
        if reads[j][0][0] == '>':
            continue
        if j%4==3 and insertion==True:
            initial = saved_index
        else:
            initial = 0
        insertion = False
        current_read = reads[j][0]
        match_list = []
        for i in range(initial, len(reference)-len(current_read)):
            score = comparison(current_read, reference[i:i+len(current_read)])
            match_list.append(score)
            if score / len(current_read) > THRESHOLD:
                err_index = find_error_index(current_read, reference[i:i+len(current_read)])
                try:
                    if err_index[0] > len(current_read) - 5 or err_index[len(err_index) - 1] < 5:
                        continue
                except:
                     pass
                saved_index = i + len(current_read)
                insertion=True
                for index in err_index:
                    potential_output = '>S' + str(match_list.index(max(match_list)) + index) + ' ' + str(reference[i + index]) + ' ' + str (current_read[index])
                    #if potential_output not in outputs:
                    outputs.append(potential_output)
            # elif score / len(current_read) > DELETION_THRESHOLD:
            #     current_reference = reference[i:i+len(current_read)]
            #     for j in range(len(current_reference)):
            #         new_reference = remove(current_reference, j)
            #         new_score = comparison(current_read[:-1], new_reference)
            #         if new_score / (len(current_read) - 1) > THRESHOLD:
            #             # err_index = find_error_index(current_read, new_reference[j:j+len(current_read)])
            #             # for index in err_index:
            #             #     potential_output = '>D' + str(match_list.index(max(match_list)) + j) + ' ' + str(new_reference[j + index])
            #             #     if potential_output not in outputs:
            #             #         outputs.append(potential_output)
            #             if j < 5 or j > 45:
            #                 continue
            #             potential_output = '>D' + str(match_list.index(max(match_list)) + j - 1) + ' ' + str(reference[i + j])
            #             if potential_output not in outputs:
            #                 outputs.append(potential_output)
        max_match_list.append(max(match_list))
    outputs.sort(key=lambda x:int(x.split(' ')[0][2:]))
    final_outputs = []
    for output in outputs:
        if outputs.count(output) > 1 and output not in final_outputs:
            final_outputs.append(output)
    
    open_file('readme2.txt', final_outputs)