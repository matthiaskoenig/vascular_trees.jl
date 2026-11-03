using VascularTrees.TreesProcessing: find_closest_nodes!

coords = permutedims(hcat(
    [[0.0, 0.0, 0.0],  [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # tree 1: reference tree
    [[0.0, 10.1, 0.0], [0.1, 0.0, 0.0],  [0.0, 0.0, 10.1], [10.1, 0.0, 0.0]],  # tree 2
    [[0.0, 0.0, 9.9],  [9.9, 0.0, 0.0],  [0.0, 0.1, 0.0],  [0.0, 9.9, 0.0]],   # tree 3
    [[10.0, 0.2, 0.0], [0.0, 0.0, 10.2], [0.0, 9.8, 0.0],  [0.0, 0.0, 0.2]],   # tree 4
))
# greedy step + alternating step (auction or Hungarian) for the closest terminal nodes
m_options = MatchingOptions(
    k_greedy = 3,
    max_iter = 10,
    tol = 1e-8,
    algorithm = AuctionAssignment(),
    compute_bound = false,
)

find_closest_nodes!(coords, m_options)