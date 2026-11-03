"""
    assign_tree(algorithm, workspace, centroids, trees, groups, t, cost) -> assignment

Give the nodes of tree `t` to the groups: `assignment[g]` is the node of tree `t` for group `g`,
every node is used once and the sum of `d(c_g, p_i)` is the smallest.
`centroids[g]` is the centroid of group `g` without tree `t`, `cost` the current total cost.
"""
function assign_tree(::HungarianAssignment, workspace, centroids, trees, groups, t, cost)
    (; cost_matrix) = workspace
    # cost_matrix[g, i] = d(c_g, p_i), p_i - coordinates of node i of tree t
    for i in axes(trees, 2), g in axes(groups, 2)          # column-major: g inner
        cost_matrix[g, i] = sq_dist(centroids[g], trees[t, i])
    end
    assignment, _ = hungarian(cost_matrix)
    return assignment
end