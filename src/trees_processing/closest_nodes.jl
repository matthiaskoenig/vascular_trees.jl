using NearestNeighbors,
    StaticArrays,
    Statistics,
    Hungarian

sq_dist(p, q) = sum((p .- q) .^ 2)

function calculate_group_cost(trees, groups, group_idx)
    μ = mean(trees[t, groups[t, group_idx]] for t in axes(groups, 1))

    group_cost = sum(sq_dist(trees[t, groups[t, group_idx]], μ) for t in axes(groups, 1))
    return group_cost
end

function calculate_total_cost(trees, groups)
    return sum(group_idx -> calculate_group_cost(trees, groups, group_idx), axes(groups, 2))
end

# neighbours[t] holds the nodes of tree t + 1 closest to each node of the reference tree
function find_neighbours!(neighbours, trees, reference_tree, k)
    for (t, tree) in enumerate(eachrow(trees[2:end, :]))
        neighbours[t], _ = knn(KDTree(tree), reference_tree, min(k, length(tree)))
    end
end

function select_greedy!(groups, order, candidates, taken, n_groups)
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

function assign_free_nodes!(groups, taken, trees)
    for t in 3:size(trees, 1)
        free_nodes = taken[t, :]
        # groups that still have no node from tree t
        open_groups = findall(==(0), groups[t, :])
        # centroid of each open group from its already assigned trees
        centroids = [mean(trees[t_other, groups[t_other, group_idx]] for t_other in axes(groups, 1) if groups[t_other, group_idx] != 0) for group_idx in open_groups]

        # all (open group index, node index) combinations, scored by distance from the node to the group's centroid
        options = collect(Iterators.product(eachindex(open_groups), free_nodes))
        costs = [sq_dist(centroids[group_idx], trees[t, node_idx]) for (group_idx, node_idx) in options]
        # smallest distance first
        order = sortperm(vec(costs))

        group_filled = falses(length(open_groups))
        node_taken = falses(size(trees, 2))
        for o in order
            group_idx, node_idx = options[o]
            (group_filled[group_idx] || node_taken[node_idx]) && continue
            groups[t, open_groups[group_idx]] = node_idx
            group_filled[group_idx] = true                # mark the group and the node as taken
            node_taken[node_idx] = true
        end
    end
end

function pair_free_nodes!(groups, taken, trees, n_groups, idx)
    free_1 = taken[idx[1], :]
    free_2 = taken[idx[2], :]

    pairs = Matrix{Int}(undef, 2, length(free_1) * length(free_2))
    combinations = Iterators.product(free_1, free_2)
    for (c, (a, b)) in enumerate(combinations)
        pairs[:, c] .= (a, b)
    end
    costs = [calculate_group_cost(trees, pairs, c) for c in axes(pairs, 2)]
    order = sortperm(costs)

    return select_greedy!(view(groups, 1:2, :), order, pairs, taken, n_groups)
end

function get_greedy_solution!(groups, reference_tree, trees, k)
    n_trees = size(trees, 1)
    n_other = n_trees - 1
    n_neighbours = min(k, size(trees, 2))
    n_reference = length(reference_tree)
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

function alternating_assignment(groups, trees, greedy_cost, max_iter, tol)
    rel_change = 1
    cost = greedy_cost
    for _ in 1:max_iter
        for t in axes(trees, 1)

            # centroid of each group without tree t
            other_trees = [t_other for t_other in axes(groups, 1) if t_other != t]
            centroids = [mean(trees[t_other, groups[t_other, group_idx]] for t_other in other_trees) for group_idx in axes(groups, 2)]

            cost_matrix = [sq_dist(centroids[group_idx], trees[t, node_idx])
                           for group_idx in axes(groups, 2), node_idx in axes(trees, 2)]
            assignment, _ = hungarian(cost_matrix)
            groups[t, :] .= assignment
        end

        new_cost = calculate_total_cost(trees, groups)
        rel_change = (cost - new_cost) / cost
        cost = new_cost
        rel_change <= tol && break
    end
    return groups
end

function find_closest_nodes(coords::Matrix, k::Int)
    trees = SVector{3,Float64}.(coords)
    n_trees = size(coords, 1)
    n_nodes = size(coords, 2)
    groups = zeros(Int, n_trees, n_nodes)

    reference_tree = @view trees[1, :]

    max_iter = 10
    tol = 10e-09

    greedy_cost = get_greedy_solution!(groups, reference_tree, trees, k)
    alternating_assignment(groups, trees, greedy_cost, max_iter, tol)
end

coords = permutedims(hcat(
    [[0.0, 0.0, 0.0],  [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],  # tree 1: reference tree
    [[0.0, 10.1, 0.0], [0.1, 0.0, 0.0],  [0.0, 0.0, 10.1], [10.1, 0.0, 0.0]],  # tree 2
    [[0.0, 0.0, 9.9],  [9.9, 0.0, 0.0],  [0.0, 0.1, 0.0],  [0.0, 9.9, 0.0]],   # tree 3
    [[10.0, 0.2, 0.0], [0.0, 0.0, 10.2], [0.0, 9.8, 0.0],  [0.0, 0.0, 0.2]],   # tree 4
))
result = find_closest_nodes(coords, 1)
