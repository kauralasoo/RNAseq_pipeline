import pandas as pd
import yaml


def extract_samples_metadata(filepath):

    samples = []
    metadata = []

    with open(filepath) as f:
        for i, line in enumerate(f):
            split_line = line.split('\t')

            samples.append(split_line[1])

            data_pieces = split_line[2].rstrip().split(';')

            split_pieces = [x.split('=') for x in data_pieces]

            split_pieces_concat = []

            for j in range(len(split_pieces) - 1):
                curr = split_pieces[j]

                if len(curr) == 1:
                    continue

                counter = 1
                next_item = split_pieces[j + counter]
                while len(next_item) == 1:
                    val = next_item[0]

                    if len(val) > 0:
                        curr.append(val)
                    counter += 1
                    try:
                        next_item = split_pieces[j + counter]
                    except IndexError:
                        next_item = ''

                split_pieces_concat.append(curr)
            metadata.append(split_pieces_concat)

    return samples, metadata


t_cells_samples, t_cells_metadata = extract_samples_metadata('Tcell_analysis_sample_meta_info.map')
monocytes_samples, monocytes_metadata = extract_samples_metadata('monocyte_analysis_sample_meta_info.map')
neutrophils_samples, neutrophils_metadata = extract_samples_metadata('neutrophil_analysis_sample_meta_info.map')

total_raw_samples = t_cells_samples + monocytes_samples + neutrophils_samples
total_raw_metadata = t_cells_metadata + monocytes_metadata + neutrophils_metadata

total_raw_metadata_dicts = []
for metadata in total_raw_metadata:
    metadata_dict = dict()

    for elem in metadata:
        metadata_dict[elem[0]] = ','.join(elem[1:])  # If there are more than one elements for keyword, comma separate them

    total_raw_metadata_dicts.append(metadata_dict)

# Acquire all of the possible columns

datacols = set()
for row in total_raw_metadata:
    colnames = [x[0] for x in row]
    datacols |= set(colnames)

data = dict()

for datacol in datacols:
    data[datacol] = []

featurecounts = pd.read_csv('./gene_expression_featureCounts.txt', sep='\t')

samples = featurecounts.columns.tolist()[2:]
del featurecounts

for sample in samples:
    ind = total_raw_samples.index(sample)
    sample_metadata_dict = total_raw_metadata_dicts[ind]

    for key in data:
        #print(key, sample_metadata_dict.get(key, 'NA'))
        data[key].append(sample_metadata_dict.get(key, 'NA'))

df = pd.DataFrame(data=data, index=samples)

df.to_csv('./metadata.csv')