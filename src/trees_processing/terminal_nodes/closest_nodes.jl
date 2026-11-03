# d(p, q) = (p_x - q_x)² + (p_y - q_y)² + (p_z - q_z)²
"""
    sq_dist(p, q)

Squared distance between two nodes, given by their coordinates `p = (x, y, z)` and `q`.
"""
sq_dist(p, q) = sum((p .- q) .^ 2)

"""
    calculate_group_cost(trees, groups, group_idx)

Cost of one group (column `group_idx` of `groups`, it holds one node index per tree): the
centroid of the coordinates of the group's nodes is found, then the squared distances from every
node of the group to this centroid are summed.
Smaller cost - the nodes of the group are closer to each other.
"""
function calculate_group_cost(trees, groups, group_idx)
    # p_t - coordinates of the node of tree t in this group, m - number of trees
    # centroid: μ = (p_1 + p_2 + ... + p_m) / m
    μ = mean(trees[t, groups[t, group_idx]] for t in axes(groups, 1))

    # cost = d(p_1, μ) + d(p_2, μ) + ... + d(p_m, μ)
    group_cost = sum(sq_dist(trees[t, groups[t, group_idx]], μ) for t in axes(groups, 1))
    return group_cost
end

"""
    calculate_total_cost(trees, groups)

Sum of the costs of all groups (see [`calculate_group_cost`](@ref)). This is the value we want
to make as small as possible.
"""
function calculate_total_cost(trees, groups)
    # total cost = cost(group 1) + cost(group 2) + ... + cost(group n)
    return sum(group_idx -> calculate_group_cost(trees, groups, group_idx), axes(groups, 2))
end

"""
    calculate_lower_bound(n_trees, trees)

The smallest total cost that is possible in theory. For every pair of trees the best pairing of
their nodes (smallest sum of squared distances between paired nodes) is found with the Hungarian
algorithm; the costs of all pairs are summed and divided by `n_trees`. No grouping can have a
smaller cost, so it shows how far our solution can still be improved.

Why it works: the cost of one group can also be written with the distances between its nodes,
`d(p_1, μ) + ... + d(p_m, μ) = (sum of d(p_s, p_t) over all pairs of trees s < t) / m`.
"""
function calculate_lower_bound(n_trees, trees)
    # calculate the cost that can be reached when using hungarian algorithm to change douplets
    # bound = (best_pairing(1, 2) + best_pairing(1, 3) + ... + best_pairing(m-1, m)) / m
    bound = 0.0
    for (X, Y) in combinations(1:n_trees, 2)
        # n×n: every node of tree X against every node of tree Y
        # cost_matrix[i, j] = d(node i of tree X, node j of tree Y)
        cost_matrix = [sq_dist(x, y) for x in view(trees, X, :), y in view(trees, Y, :)]
        # best_pairing(X, Y) = smallest sum of cost_matrix[i, j] where every i and every j is used once
        _, pair_cost = hungarian(cost_matrix)
        bound += pair_cost
    end
    return bound / n_trees
end

"""
    find_neighbours!(neighbours, trees, reference_tree, k)

For every node of the reference tree (usually A) find the `k` closest nodes (by coordinates) in
each of the other trees (knn). `neighbours[t]` holds these nodes for tree `t + 1`.
"""
function find_neighbours!(neighbours, trees, reference_tree, k)
    # use knn algorithm to find closest nodes to the reference tree (usually, it is A) from other trees
    # only tree 1 is the reference; the other trees are searched for its neighbours
    for t in 1:size(trees, 1) - 1
        tree = view(trees, t + 1, :)
        # neighbours[t] holds the nodes of tree t + 1 closest to each node of the reference tree
        neighbours[t], _ = knn(KDTree(tree), reference_tree, min(k, length(tree)))
    end
end

"""
    select_greedy!(groups, order, candidates, taken, n_groups) -> n_groups

Go through the candidate groups from the best (smallest cost) to the worst and take a candidate
only if none of its nodes is taken yet. Taken candidates are written to `groups`, their nodes
are marked in `taken`. Returns the new number of filled groups.
"""
function select_greedy!(groups, order, candidates, taken, n_groups)
    # select pairs - each node can be taken in the par only once
    for c in order
        candidate = view(candidates, :, c)
        any(taken[t, candidate[t]] for t in eachindex(candidate)) && continue
        n_groups += 1
        groups[:, n_groups] .= candidate
        for t in eachindex(candidate)
            taken[t, candidate[t]] = true # mark nodes as taken
        end
    end
    return n_groups
end

"""
    pair_free_nodes!(groups, taken, trees, n_groups, idx) -> n_groups

Pair the nodes of the trees `idx[1]` and `idx[2]` that are still free after
[`select_greedy!`](@ref): all combinations are scored by their cost and taken greedily
(see [`select_greedy!`](@ref)). Returns the new number of filled groups.
"""
function pair_free_nodes!(groups, taken, trees, n_groups, idx)
    # indices of the nodes that are not taken yet (taken holds true/false flags)
    free_1 = findall(!, taken[idx[1], :])
    free_2 = findall(!, taken[idx[2], :])

    # every combination (a, b): a - free node of tree idx[1], b - free node of tree idx[2]
    pairs = Matrix{Int}(undef, 2, length(free_1) * length(free_2))
    combinations = Iterators.product(free_1, free_2)
    for (c, (a, b)) in enumerate(combinations)
        pairs[:, c] .= (a, b)
    end
    # cost of a pair = d(p_a, μ) + d(p_b, μ), μ = (p_a + p_b) / 2  (= d(p_a, p_b) / 2)
    costs = [calculate_group_cost(trees, pairs, c) for c in axes(pairs, 2)]
    order = sortperm(costs)

    return select_greedy!(view(groups, 1:2, :), order, pairs, taken, n_groups)
end

"""
    assign_free_nodes!(groups, taken, trees)

Fill the groups that still miss a node from tree 3, 4, ...: for each such tree, the free
nodes are given to the open groups, the node closest to the group's centroid first. Each node
and each group is used only once.
"""
function assign_free_nodes!(groups, taken, trees)
    # get the pair to the nodes that were freed before (pair_free_nodes!) - their pair with other nodes was worse that these other nodes with another node from this tree
    for t in 3:size(trees, 1)
        # indices of the nodes of tree t that are not taken yet
        free_nodes = findall(!, taken[t, :])
        # groups that still have no node from tree t
        open_groups = findall(==(0), groups[t, :])
        # centroid of each open group from its already assigned trees
        # c_g = (sum of p_s over trees s that already have a node in group g) / (number of such trees)
        centroids = [mean(trees[t_other, groups[t_other, group_idx]] for t_other in axes(groups, 1) if groups[t_other, group_idx] != 0) for group_idx in open_groups]

        # costs[g, i] = d(c_g, p_i), p_i - coordinates of free node i of tree t
        costs = [sq_dist(centroids[g], trees[t, free_nodes[i]])
            for g in eachindex(open_groups), i in eachindex(free_nodes)]
        # smallest distance first
        order = sortperm(vec(costs))
        idx = CartesianIndices(costs)

        group_filled = falses(length(open_groups))
        node_taken = falses(size(trees, 2))
        # lower sq distance - better combination, no duplications, principle as in select_greedy!
        for o in order
            g, i = Tuple(idx[o])
            node_idx = free_nodes[i]
            (group_filled[g] || node_taken[node_idx]) && continue
            groups[t, open_groups[g]] = node_idx
            group_filled[g] = true                # mark the group and the node as taken
            node_taken[node_idx] = true
        end
    end
end

"""
    get_greedy_solution!(groups, reference_tree, trees, k) -> greedy_cost

First (greedy) solution:
1. for every node of the reference tree all combinations with its `k` neighbours (closest by
coordinates) from the other trees are built and scored (see [`calculate_group_cost`](@ref));
2. the best combinations are taken, each node only once ([`select_greedy!`](@ref));
3. the nodes that are left are paired ([`pair_free_nodes!`](@ref), trees 1 and 2) and added to
the groups ([`assign_free_nodes!`](@ref), trees 3, ...).
`groups` is filled in place. Returns the total cost of this solution.
"""
function get_greedy_solution!(groups, reference_tree, trees, k)
    n_trees = size(trees, 1)
    n_other = n_trees - 1
    n_neighbours = min(k, size(trees, 2))
    n_reference = length(reference_tree)
    # every reference node × k choices in each of the other (m - 1) trees:
    # n_candidates = n_reference * k^(m - 1)
    n_candidates = n_reference * n_neighbours^n_other

    neighbours = Vector{Vector{Vector{Int}}}(undef, n_other)
    candidates = Matrix{Int}(undef, n_trees, n_candidates)
    taken = falses(n_trees, size(trees, 2)) # one flag vector per tree
    costs = Vector{Float64}(undef, n_candidates)

    find_neighbours!(neighbours, trees, reference_tree, k)

    c = 0
    for node_idx in eachindex(reference_tree)
        # for this node collect all neighbour nodes from all trees
        neighbour_lists = [neighbours[t][node_idx] for t in 1:n_other]
        # get all possible combinations of this node and neighbours to check them
        for neighbour_combo in Iterators.product(neighbour_lists...)
            c += 1
            candidates[:, c] .= (node_idx, neighbour_combo...)
            costs[c] = calculate_group_cost(trees, candidates, c)
        end
    end

    order = sortperm(costs)

    n_groups = 0 # number of filled groups (columns)
    n_groups = select_greedy!(groups, order, candidates, taken, n_groups)
    n_groups = pair_free_nodes!(groups, taken, trees, n_groups, [1, 2])

    # --- attach free nodes of trees 3.. to the pairs of trees 1–2 ---
    assign_free_nodes!(groups, taken, trees)

    greedy_cost = calculate_total_cost(trees, groups)

    return greedy_cost
end

"""
    alternating_assignment!(groups, trees, greedy_cost, max_iter, tol) -> cost

Improve the groups step by step. For every tree: take the centroids of the groups without this
tree and give the nodes of the tree to the groups again in the best way (Hungarian algorithm).
A tree is solved again only if another tree's assignment changed since its last solve.
Stops when no tree has to be solved again, the cost improves by less than `tol` (relative), or
`max_iter` rounds are done. Returns the new total cost.
"""
function alternating_assignment!(groups, trees, greedy_cost, max_iter, tol)
    cost = greedy_cost
    n = size(groups, 2)
    m = size(trees, 1)

    needs_solve = trues(size(trees, 1))   # needs_solve[t]: tree t has to be solved (again)
    cost_matrix = Matrix{Float64}(undef, n, n)
    centroids = Vector{SVector{3,Float64}}(undef, n)
    # S[g] = p_1 + p_2 + ... + p_m, sum of the nodes of group g
    S = [sum(trees[s, groups[s, g]] for s in 1:m) for g in axes(groups, 2)]

    for _ in 1:max_iter
        for t in axes(trees, 1)
            needs_solve[t] || continue
            # centroid of group g without tree t = (sum of all m nodes - node of tree t) / (m - 1)
            # c_g = (S[g] - p_t) / (m - 1)
            for g in axes(groups, 2)
                centroids[g] = (S[g] - trees[t, groups[t, g]]) / (m - 1)
            end
            # cost_matrix[g, i] = d(c_g, p_i), p_i - coordinates of node i of tree t
            for i in axes(trees, 2), g in 1:n          # column-major: g inner
                cost_matrix[g, i] = sq_dist(centroids[g], trees[t, i])
            end
            # assignment[g] = node of tree t for group g; every node used once, sum of costs smallest
            assignment, _ = hungarian(cost_matrix)
            # update S
            for g in axes(groups, 2)                                   
                S[g] += trees[t, assignment[g]] - trees[t, groups[t, g]]
            end
            if assignment != view(groups, t, :)
                groups[t, :] .= assignment
                needs_solve .= true        # every other tree's subproblem changed
            end
            needs_solve[t] = false      # t was just solved
        end
        any(needs_solve) || break          # converged: nothing left to re-solve
        new_cost = calculate_total_cost(trees, groups)
        # rel_change = (old cost - new cost) / old cost
        rel_change = (cost - new_cost) / cost
        cost = new_cost
        rel_change <= tol && break
    end
    return cost
end

"""
    find_closest_nodes!(coords, k, compute_bound) -> groups

Group the nodes of several trees so that each group has one node from every tree and the nodes
in a group are as close as possible. `coords[t, i]` holds the coordinates `(x, y, z)` of node
`i` of tree `t`; `k` is the number of neighbours checked in the greedy step. If
`compute_bound` is `true`, the lower bound ([`calculate_lower_bound`](@ref)) is also computed
and it is logged how much the result can still be improved at most; this is slow
(m(m - 1)/2 Hungarian runs).

`coords` is reordered in place: afterwards column `g` holds the coordinates of the nodes of
group `g`, one node per tree. Returns `groups`: `groups[t, g]` is the original index of the node
of tree `t` in group `g`.
"""
function find_closest_nodes!(coords::Matrix, k::Int, compute_bound::Bool)
    trees = SVector{3,Float64}.(coords)
    n_trees = size(coords, 1)
    n_nodes = size(coords, 2)
    groups = zeros(Int, n_trees, n_nodes)

    reference_tree = @view trees[1, :]

    max_iter = 10
    tol = 10e-09

    greedy_cost = get_greedy_solution!(groups, reference_tree, trees, k)
    cost = alternating_assignment!(groups, trees, greedy_cost, max_iter, tol)

    if compute_bound
        bound = calculate_lower_bound(n_trees, trees)
        # gap = (cost - bound) / bound * 100 [%]
        gap = (cost - bound) / bound * 100
        @info "At most $(gap)% improvement is left"
    end

    # new coords[t, g] = old coords[t, groups[t, g]]: column g holds the nodes of group g
    # (the right side is copied first, so the row is not overwritten while it is read)
    for t in axes(coords, 1)
        coords[t, :] .= coords[t, groups[t, :]]
    end
    return groups
end

# coords = permutedims(hcat(
#     [[0.0, 0.0, 0.0],  [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # tree 1: reference tree
#     [[0.0, 10.1, 0.0], [0.1, 0.0, 0.0],  [0.0, 0.0, 10.1], [10.1, 0.0, 0.0]],  # tree 2
#     [[0.0, 0.0, 9.9],  [9.9, 0.0, 0.0],  [0.0, 0.1, 0.0],  [0.0, 9.9, 0.0]],   # tree 3
#     [[10.0, 0.2, 0.0], [0.0, 0.0, 10.2], [0.0, 9.8, 0.0],  [0.0, 0.0, 0.2]],   # tree 4
# ))
# find_closest_nodes!(coords, 1)