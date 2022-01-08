from scripts import *
from reconstruction import *
from sequencing import *

file = open('test', 'r')
out = open('output.txt', 'w')


def testComp():
    k = int(file.readline())
    text = file.readline()
    for elem in composition(text, k):
        print(elem, end='\n')


def testOverlap():
    p = list()
    for line in file:
        p.append(line.strip('\n'))
    print(overlap(p))


def testOverlapGraph():
    lines = list()
    for line in file:
        lines.append(line.strip())
    overlap_graph(lines)


def test_deBruijin():
    k = int(file.readline())
    text = file.readline()
    comp = composition(text, k)
    out.write('\n'.join(deBruijin_graph(comp)))


def test_deBruijin2():
    text = list()
    for line in file:
        text.append(line[:-1])
    out.write('\n'.join(deBruijin_graph(text)))


def test_cycle():
    edges = {}
    for line in file.readlines():
        line = line.strip().split('->')
        edges[int(line[0])] = []
        if ',' in line[1]:
            for elem in line[1].split(','):
                edges[int(line[0])].append(int(elem))
        else:
            edges[int(line[0])].append(int(line[1]))
    out.write('->'.join(map(str, eulerian_cycle(edges))))


def test_path():
    edges = {}
    for line in file.readlines():
        line = line.strip().split('->')
        edges[int(line[0])] = []
        if ',' in line[1]:
            for elem in line[1].split(','):
                edges[int(line[0])].append(int(elem))
        else:
            edges[int(line[0])].append(int(line[1]))
    out.write('->'.join(map(str, eulerian_path(edges))))


def test_reconstruction():
    k = file.readline().strip()
    patterns = []
    for line in file:
        patterns.append(line.strip())
    out.write(string_reconstruction(patterns, k))


def test_universal_string():
    k = int(file.readline())
    out.write(universal_string(k))


def test_pairReconstruction():
    temp = file.readline().strip().split(' ')
    k = int(temp[0])
    d = int(temp[1])
    pairs = list()
    for line in file.readlines():
        pairs.append(tuple(line.strip().split('|')))
    out.write(pair_reconstruction(k, d, pairs))


def test_gapped():
    temp = file.readline().strip().split(' ')
    k = int(temp[0])
    d = int(temp[1])
    pairs = list()
    for line in file.readlines():
        pairs.append(tuple(line.strip().split('|')))
    out.write(string_from_gapped_patterns(pairs, k, d))


def test_translate():
    rna = file.readline()
    out.write(translate_RNA(rna))


def test_peptide_encoding():
    dna = file.readline().strip()
    peptide = file.readline()
    for elem in peptide_encoding(dna, peptide):
        out.write(elem + '\n')


def get_aa_mass():
    mass_dict = {}
    with open('integer_mass_table.txt', 'r') as mass_table:
        for line in mass_table.readlines():
            mass_dict[line.split()[0]] = line.strip().split()[1]
    return mass_dict


def test_cyclospectrum():
    peptide = file.readline()
    spectrum = str(cyclospectrum(peptide, get_aa_mass())).replace(',', '')
    out.write(spectrum)


def test_lienarspectrum():
    peptide = file.readline().strip()
    spectrum = str(linear_spectrum(peptide)).replace(',', '')
    out.write(spectrum)


def test_count():
    out.write(str(count_mass(int(file.readline()), {}, get_aa_mass())[0]))


def test_cyclopeptide():
    spectrum = []
    for x in file.readline().strip().split():
        spectrum.append(int(x))
    output = cyclopeptide_sequencing(spectrum)
    for elem in output[0]:
        out.write('-'.join(elem) + ' ')
    out.write('\n' + str(output[1]))


def test_score():
    peptide = file.readline().strip()
    spectrum = [int(x) for x in file.readline().split()]
    out.write(str(linear_score(peptide, spectrum)))


def test_trim():
    board = [peptide for peptide in file.readline().strip().split()]
    spectrum = [int(x) for x in file.readline().strip().split()]
    n = int(file.readline())
    out.write(str(trim(board, spectrum, n)).replace(',', '').replace("'", ''))


def test_leaderboard():
    n = int(file.readline().strip())
    spectrum = [int(x) for x in file.readline().split()]
    peptide = leaderboard_sequencing(spectrum, n)
    out.write(peptide[0] + '\n' + '-'.join(peptide[1]))


def test_convolution():
    spectrum = [int(mass) for mass in file.readline().split()]
    out.write(''.join(str(convolution(spectrum))))


def test_convolution_sequencing():
    m = file.readline().strip()
    n = file.readline().strip()
    spectrum = [int(mass) for mass in file.readline().split()]
    out.write(convolution_sequencing(m, n, spectrum))


# testComp()
# testOverlap()
# testOverlapGraph()
# test_deBruijin()
# test_deBruijin2()
# test_cycle()
# test_path()
test_reconstruction()
# test_universal_string()
# test_pairReconstruction()
# test_gapped()
# test_translate()
# test_peptide_encoding()
# test_cyclospectrum()
# test_lienarspectrum()
# test_count()
# test_cyclopeptide()
# test_score()
# test_trim()
# test_leaderboard()
# test_convolution()

file.close()
