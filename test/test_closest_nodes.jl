using VascularTrees.TreesProcessing: find_closest_nodes!

coords = permutedims(hcat(
    [[0.0, 0.0, 0.0],  [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # tree 1: reference tree
    [[0.0, 10.1, 0.0], [0.1, 0.0, 0.0],  [0.0, 0.0, 10.1], [10.1, 0.0, 0.0]],  # tree 2
    [[0.0, 0.0, 9.9],  [9.9, 0.0, 0.0],  [0.0, 0.1, 0.0],  [0.0, 9.9, 0.0]],   # tree 3
    [[10.0, 0.2, 0.0], [0.0, 0.0, 10.2], [0.0, 9.8, 0.0],  [0.0, 0.0, 0.2]],   # tree 4
))
k = 3
compute_bound = true

find_closest_nodes!(coords, k, compute_bound)