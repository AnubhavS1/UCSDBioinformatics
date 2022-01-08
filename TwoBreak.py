# function to convert chromosome to cycle
def chromosome_cycle(chromosome):
    nodes = [0] * len(chromosome) * 2
    for i in range(len(chromosome)):
        j = chromosome[i]
        if j > 0:
            nodes[2*i] = 2 * j - 1
            nodes[2*i+1] = 2 * j
        else:
            nodes[2*i] = -2 * j
            nodes[2*i+1] = -2 * j - 1
    return nodes


# find edges connecting synteny blocks to each other
def colored_edges(p):
    edges = {}
    for chromosome in p:
        nodes = chromosome_cycle(chromosome)
        for i in range(1, len(chromosome)):
            edges[nodes[2*i-1]] = nodes[2*i]
        edges[nodes[len(nodes)-1]] = nodes[0]
    return edges


# count number of cycles in a breakpoint graph and return a set of them
def allCycles(bp_graph, available):
    # construct cycles while counting them
    num_cycles = 0
    cycles = []
    while len(available) > 0:
        num_cycles += 1
        cycles.append([])
        cycle = [available[0]]
        cycles[num_cycles-1].append(cycle[0])
        # as long as a cycle doesn't (i.e. it doesn't reach a value its already been to) continue
        while cycle:
            keys = list(bp_graph.keys())
            values = list(bp_graph.values())
            current = cycle.pop(0)
            if current not in available:
                break
            if current in values:
                cycle.append(keys[values.index(current)])
                cycles[num_cycles-1].append(cycle[0])
                bp_graph.pop(keys[values.index(current)])
            elif current in keys:
                if type(bp_graph[current]) == list and len(list(bp_graph[current])) > 0:
                    cycle.append(bp_graph[current][0])
                    cycles[num_cycles-1].append(cycle[0])
                    bp_graph[current].remove(bp_graph[current][0])
                elif type(bp_graph[current]) == int:
                    cycle.append(bp_graph[current])
                    cycles[num_cycles-1].append(cycle[0])
                    bp_graph.pop(current)
            else:
                for i in range(len(bp_graph.keys())):
                    if len(cycle) > 0:
                        break
                    if type(bp_graph[keys[i]]) == list:
                        for num in bp_graph[keys[i]]:
                            if num == current:
                                cycle.append(keys[i])
                                cycles[num_cycles-1].append(cycle[0])
                                bp_graph[keys[i]].remove(num)
                                break
            available.remove(current)
    return num_cycles, cycles


# calculate 2 break distance between two genomes P and Q
def two_break_dist(p, q):
    blocks = 0
    for elem in p:
        for num in elem:
            blocks += 1
    # construct break point graphs of p and q
    bp_graph = {}
    vals = []
    colp = colored_edges(p)
    for edge in colp:
        bp_graph[edge] = [colp[edge]]
        vals.append(colp[edge])
    colq = colored_edges(q)
    for edge in colq:
        if edge in bp_graph.keys():
            bp_graph[edge].append(colq[edge])
        else:
            bp_graph[edge] = [colq[edge]]
        vals.append(colq[edge])
    for key in bp_graph:
        if len(bp_graph[key]) == 1:
            bp_graph[key] = bp_graph[key][0]
    # count number of cycles in graph
    available = list(bp_graph.keys())
    for val in vals:
        if val not in available:
            available.append(val)
    return blocks - allCycles(bp_graph, available)[0]


# get the genome from a colored edges on genome graph
def black_edges(col: dict):
    edges = []
    for pair in list(col.items()):
        tmp = []
        for num in pair:
            if num % 2 == 0:
                tmp.append(num/2)
            else:
                tmp.append(-(num+1)/2)
        edges.append(tmp)
    return edges


# convert cycle in graph to a genome form
def cycletochromosome(nodes):
    chromosome = []
    for i in range(int(len(nodes)/2)):
        if nodes[2*i-1] < nodes[2*i]:
            chromosome.append(int(nodes[2*i] / 2))
        else:
            chromosome.append(int(-nodes[2*i-1] / 2))
    return chromosome


def graphtogenome(edges: list):
    cycles = []
    for num in edges:
        if edges[0] % 2 == 0:
            if num == edges[0] - 1:
                cycles.append(edges[0:edges.index(num)+1])
                edges = edges[edges.index(num)+1:]
        else:
            if num == edges[0] + 1:
                cycles.append(edges[0:edges.index(num) + 1])
                edges = edges[edges.index(num)+1:]
    p = []
    for cycle in cycles:
        p.append(cycletochromosome(cycle))
    return p


def two_breakgenomegraph(edges, i, j, k, l):
    if [i, j] in edges:
        edges.remove([i, j])
    else:
        edges.remove([j, i])
    if [k, l] in edges:
        edges.remove([k, l])
    else:
        edges.remove([l, k])
    edges.append([i, k])
    edges.append([j, l])
    return edges


def two_breakgenome(p, i, j, k, l):
    return two_breakgenomegraph(colored_edges(p), i, j, k, l)
