import copy


# calculate the indegree and outdegree of each node
def inout(graph):
    in_dict = {}
    out_dict = {}
    for key in graph.keys():
        out_dict[key] = len(graph[key])
    for value in graph.values():
        for v in value.keys():
            if v in in_dict:
                in_dict[v] += 1
            else:
                in_dict[v] = 1

    return in_dict, out_dict


# trim all nodes that have 0 indegree or 0 outdegree that aren't start/end nodes
def trim(graph, degree, start, end):
    trimmed = copy.deepcopy(graph)
    for key in graph.keys():
        if key in degree[1] and key not in degree[0] and key != start:
            del trimmed[key]
    graph = copy.deepcopy(trimmed)
    for key in graph.keys():
        for k in graph[key].keys():
            if k in degree[0] and k not in degree[1] and k != end:
                del trimmed[key][k]
    return trimmed


def topological_ordering(graph: dict, start, degree):
    ordered = {start: graph[start]}
    in_sorted = {k: v for k, v in sorted(degree[0].items(), key=lambda item: item[1])}
    for key in in_sorted:
        try:
            ordered[key] = graph[key]
        except KeyError:
            continue
    return ordered


# global variable for longest path from start to sink
path = []


def score(node, graph, start, degree):
    if node == start:
        return 0
    # ind = degree[0][node]
    if degree[0][node] == 1:
        for n in graph:
            if node in graph[n]:
                path.append(str(node))
                return score(n, graph, start, degree) + graph[n][node]
    else:
        path.append(str(node))
        predecessors = [key for key in graph if node in graph[key]]
        pred_scores = {score(n, graph, start, degree) + graph[n][node]: n for n in predecessors}
        for value in pred_scores.values():
            check = pred_scores[max(pred_scores)]
            if value != check and str(value) in path:
                path.remove(str(value))
        return max(pred_scores)


def longest_path(graph: dict, start, end):
    trim_graph = None
    degree = inout(graph)
    for i in range(100):
        trim_graph = topological_ordering(trim(graph, degree, start, end), start, degree)
        degree = inout(trim_graph)
        trim_graph = topological_ordering(trim(graph, degree, start, end), start, degree)
    path.reverse()
    return score(end, trim_graph, start, degree), (str(start) + '->' + '->'.join(path))
