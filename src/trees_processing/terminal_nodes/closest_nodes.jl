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
    μ = zero(eltype(trees))
    for t in axes(groups, 1)
        μ += trees[t, groups[t, group_idx]]
    end
    μ /= size(groups, 1)

    # cost = d(p_1, μ) + d(p_2, μ) + ... + d(p_m, μ)
    group_cost = zero(eltype(eltype(trees)))
    for t in axes(groups, 1)
        group_cost += sq_dist(trees[t, groups[t, group_idx]], μ)
    end
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
        # n*n: every node of tree X against every node of tree Y
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
    kdtrees = Vector{KDTree}(undef, size(trees, 1) - 1)
    # use knn algorithm to find closest nodes to the reference tree (usually, it is A) from other trees
    # only tree 1 is the reference; the other trees are searched for its neighbours
    for t in 1:size(trees, 1) - 1
        tree = view(trees, t + 1, :)
        kdtrees[t] = KDTree(tree)
        # neighbours[t] holds the nodes of tree t + 1 closest to each node of the reference tree
        neighbours[t], _ = knn(kdtrees[t], reference_tree, min(k, length(tree)))
    end
    return kdtrees
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
    pair_free_nodes!(groups, taken, trees, n_groups, idx, k) -> n_groups

Pair the nodes of the trees `idx[1]` and `idx[2]` that are still free after
[`select_greedy!`](@ref): for every free node of tree `idx[1]` only its `k` closest free nodes
of tree `idx[2]` are scored by their cost and taken greedily (see [`select_greedy!`](@ref)).
Nodes that are still free after that are paired in index order; [`alternating_assignment!`](@ref)
improves them later. Returns the new number of filled groups.
"""
function pair_free_nodes!(groups, taken, trees, n_groups, idx, k)
    # indices of the nodes that are not taken yet (taken holds true/false flags)
    free_1 = findall(!, view(taken, idx[1], :))
    free_2 = findall(!, view(taken, idx[2], :))
    isempty(free_1) && return n_groups

    # neighbours[i] - positions (in free_2) of the free nodes of tree idx[2] closest to free_1[i]
    n_neighbours = min(k, length(free_2))
    neighbours, _ = knn(KDTree(trees[idx[2], free_2]), trees[idx[1], free_1], n_neighbours)

    # combinations (a, b): a - free node of tree idx[1], b - one of its closest free nodes of tree idx[2]
    pairs = Matrix{Int}(undef, 2, length(free_1) * n_neighbours)
    c = 0
    for (a, neighbours_a) in zip(free_1, neighbours), j in neighbours_a
        c += 1
        pairs[1, c] = a
        pairs[2, c] = free_2[j]
    end
    # cost of a pair = d(p_a, μ) + d(p_b, μ), μ = (p_a + p_b) / 2  (= d(p_a, p_b) / 2)
    costs = [calculate_group_cost(trees, pairs, c) for c in axes(pairs, 2)]
    order = sortperm(costs)

    n_groups = select_greedy!(view(groups, 1:2, :), order, pairs, taken, n_groups)

    # nodes whose closest nodes were all taken: pair the rest in index order
    for (a, b) in zip(findall(!, view(taken, idx[1], :)), findall(!, view(taken, idx[2], :)))
        n_groups += 1
        groups[1:2, n_groups] .= (a, b)
        taken[idx[1], a] = true
        taken[idx[2], b] = true
    end
    return n_groups
end

"""
    assign_free_nodes!(groups, taken, trees, k)

Fill the groups that still miss a node from tree 3, 4, ...: for each such tree, every open group
gets the `k` free nodes closest to its centroid as candidates, and the candidates are taken
the closest first. Each node and each group is used only once. Groups that are still open after
that get the remaining nodes in index order; [`alternating_assignment!`](@ref) improves them later.
"""
function assign_free_nodes!(groups, taken, trees, k)
    # get the pair to the nodes that were freed before (pair_free_nodes!) - their pair with other nodes was worse that these other nodes with another node from this tree
    for t in 3:size(trees, 1)
        # indices of the nodes of tree t that are not taken yet
        free_nodes = findall(!, view(taken, t, :))
        # groups that still have no node from tree t
        open_groups = findall(==(0), view(groups, t, :))
        isempty(open_groups) && continue
        # centroid of each open group from its already assigned trees
        # c_g = (sum of p_s over trees s that already have a node in group g) / (number of such trees)
        centroids = [mean(trees[t_other, groups[t_other, group_idx]] for t_other in axes(groups, 1) if groups[t_other, group_idx] != 0) for group_idx in open_groups]

        # neighbours[g] - positions (in free_nodes) of the free nodes closest to c_g, dists[g] - their
        # distances (not squared, but the order is the same as for d(c_g, p_i))
        n_neighbours = min(k, length(free_nodes))
        neighbours, dists = knn(KDTree(trees[t, free_nodes]), centroids, n_neighbours)

        # candidate c: group cand_group[c] with free node free_nodes[cand_node[c]]
        n_candidates = length(open_groups) * n_neighbours
        cand_group = Vector{Int}(undef, n_candidates)
        cand_node = Vector{Int}(undef, n_candidates)
        cand_dist = Vector{Float64}(undef, n_candidates)
        c = 0
        for g in eachindex(open_groups), j in eachindex(neighbours[g])
            c += 1
            cand_group[c] = g
            cand_node[c] = neighbours[g][j]
            cand_dist[c] = dists[g][j]
        end
        # smallest distance first
        order = sortperm(cand_dist)

        group_filled = falses(length(open_groups))
        node_taken = falses(length(free_nodes))
        # lower distance - better combination, no duplications, principle as in select_greedy!
        for o in order
            g, i = cand_group[o], cand_node[o]
            (group_filled[g] || node_taken[i]) && continue
            groups[t, open_groups[g]] = free_nodes[i]
            group_filled[g] = true                # mark the group and the node as taken
            node_taken[i] = true
        end

        # groups whose closest nodes were all taken: give them the rest in index order
        for (g, i) in zip(findall(!, group_filled), findall(!, node_taken))
            groups[t, open_groups[g]] = free_nodes[i]
        end
    end
end

"""
    fill_candidates!(candidates, costs, neighbours, combos, trees)

Write every candidate group to `candidates` and its cost to `costs`: each node of the reference
tree combined with every choice (`combos`) of its neighbours in the other trees. A separate
function so that it is compiled for the concrete number of trees in `combos`.
"""
function fill_candidates!(candidates, costs, neighbours, combos, trees)
    c = 0
    for node_idx in eachindex(neighbours[1]), combo in combos
        c += 1
        candidates[1, c] = node_idx
        # combo[t] - which of the neighbours of tree t + 1 is used
        for t in eachindex(neighbours)
            candidates[t + 1, c] = neighbours[t][node_idx][combo[t]]
        end
        costs[c] = calculate_group_cost(trees, candidates, c)
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
    # every reference node * k choices in each of the other (m - 1) trees:
    # n_candidates = n_reference * k^(m - 1)
    n_candidates = n_reference * n_neighbours^n_other

    neighbours = Vector{Vector{Vector{Int}}}(undef, n_other)
    candidates = Matrix{Int}(undef, n_trees, n_candidates)
    taken = falses(n_trees, size(trees, 2)) # one flag vector per tree
    costs = Vector{Float64}(undef, n_candidates)

    kdtrees = find_neighbours!(neighbours, trees, reference_tree, k)

    # all possible choices of one neighbour in each of the other trees
    combos = CartesianIndices(ntuple(_ -> n_neighbours, n_other))
    fill_candidates!(candidates, costs, neighbours, combos, trees)

    order = sortperm(costs)

    n_groups = 0 # number of filled groups (columns)
    n_groups = select_greedy!(groups, order, candidates, taken, n_groups)
    n_groups = pair_free_nodes!(groups, taken, trees, n_groups, [1, 2], k)

    # --- attach free nodes of trees 3.. to the pairs of trees 1–2 ---
    assign_free_nodes!(groups, taken, trees, k)

    greedy_cost = calculate_total_cost(trees, groups)

    return kdtrees, greedy_cost
end

"""
    make_workspace(algorithm, trees, n) -> workspace

Arrays that [`assign_tree`](@ref) needs, allocated once before the alternating step.
- [`HungarianAssignment`](@ref): the dense `n × n` cost matrix.
- [`AuctionAssignment`](@ref): one kd-tree per tree (node coordinates never change, so they are
built once), `k` of every group of every tree, the candidates of all groups and the state of
the auction.
"""
make_workspace(::HungarianAssignment, trees, n) = (; cost_matrix = Matrix{Float64}(undef, n, n))

function make_workspace(algorithm::AuctionAssignment, trees, n)
    # k can not be larger than the number of nodes
    k_min = min(algorithm.k_min, n)
    k_max = min(algorithm.k_max, n)
    kdtrees = [KDTree(view(trees, t, :)) for t in axes(trees, 1)]
    # k_groups[t, g] - number of candidates of group g in tree t, kept between the rounds
    k_groups = fill(clamp(algorithm.k_start, k_min, k_max), size(trees, 1), n)
    # candidates of group g: entries first_entry[g]:(first_entry[g + 1] - 1),
    # entry e: node entry_nodes[e] of tree t with reward entry_rewards[e] = -d(c_g, p_i)
    # at most k_max neighbours + the current node per group
    max_entries = n * (k_max + 1)
    first_entry = Vector{Int}(undef, n + 1)
    entry_nodes = Vector{Int}(undef, max_entries)
    entry_rewards = Vector{Float64}(undef, max_entries)
    # auction state: prices[i] - price of node i, owners[i] - group that has node i (0 - none),
    # assignment[g] - node of group g (0 - none)
    prices = zeros(n)
    owners = zeros(Int, n)
    assignment = zeros(Int, n)
    free_groups = sizehint!(Int[], n)
    return (; kdtrees, k_groups, k_min, k_max, first_entry, entry_nodes, entry_rewards,
            prices, owners, assignment, free_groups)
end

"""
    fill_auction_candidates!(workspace, neighbours, centroids, trees, groups, t) -> n_entries

Write the candidates of every group `g` of tree `t` into the workspace: its `k_groups[t, g]`
closest nodes (`neighbours[g]`, sorted by distance) and its current node, with
reward = -d(c_g, p_i). Entry `r` of group `g` is its `r`-th closest node; entry `k + 1` (if
present) is its current node, which is not among its `k` closest. Returns the number of entries.
"""
function fill_auction_candidates!(workspace, neighbours, centroids, trees, groups, t)
    (; k_groups, first_entry, entry_nodes, entry_rewards) = workspace
    n_entries = 0
    for g in axes(groups, 2)
        first_entry[g] = n_entries + 1
        k = k_groups[t, g]
        candidates = view(neighbours[g], 1:k)
        for i in candidates
            n_entries += 1
            entry_nodes[n_entries] = i
            entry_rewards[n_entries] = -sq_dist(centroids[g], trees[t, i])
        end
        # the current node is always a candidate, so a full assignment always exists
        current = groups[t, g]
        if current ∉ candidates
            n_entries += 1
            entry_nodes[n_entries] = current
            entry_rewards[n_entries] = -sq_dist(centroids[g], trees[t, current])
        end
    end
    first_entry[end] = n_entries + 1
    return n_entries
end

"""
    auction_assignment!(workspace, n, ε_final, θ, max_bids) -> converged

Forward auction with ε-scaling on the candidates in the workspace (rewards are maximised).
Every free group bids for its best node (value = reward - price): the price of this node rises
by best value - second best value + ε, its previous owner becomes free. When all groups have a
node, ε is divided by `θ` and the next phase starts (prices are kept), until `ε_final` is
reached; then the result is at most `n * ε_final` worse than the best assignment of the
candidates. The result is in `workspace.assignment`. Returns `false` if more than `max_bids`
bids were needed.
"""
function auction_assignment!(workspace, n, ε_final, θ, max_bids)
    (; first_entry, entry_nodes, entry_rewards, prices, owners, assignment, free_groups) = workspace
    largest_reward = maximum(abs, view(entry_rewards, 1:(first_entry[end] - 1)))
    ε = max(largest_reward / θ, ε_final)
    n_bids = 0
    while true                                       # ε-scaling phases
        fill!(owners, 0)
        fill!(assignment, 0)
        empty!(free_groups)
        append!(free_groups, 1:n)
        while !isempty(free_groups)
            g = pop!(free_groups)
            best, second, best_node = -Inf, -Inf, 0
            @inbounds for e in first_entry[g]:(first_entry[g + 1] - 1)
                value = entry_rewards[e] - prices[entry_nodes[e]]
                if value > best
                    second, best, best_node = best, value, entry_nodes[e]
                elseif value > second
                    second = value
                end
            end
            prices[best_node] += best - second + ε
            previous_owner = owners[best_node]
            if previous_owner != 0
                assignment[previous_owner] = 0
                push!(free_groups, previous_owner)
            end
            owners[best_node] = g
            assignment[g] = best_node
            (n_bids += 1) > max_bids && return false
        end
        ε <= ε_final && return true
        ε = max(ε / θ, ε_final)
    end
end

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

# k of every group of tree t is doubled (at most k_max)
function grow_k!(k_groups, t, k_max)
    for g in axes(k_groups, 2)
        k_groups[t, g] = min(2k_groups[t, g], k_max)
    end
end

# rank of the node that group g got: r - its r-th closest node, k + 1 - its current node,
# which is not among its k closest
function assigned_rank(workspace, g)
    (; first_entry, entry_nodes, assignment) = workspace
    for e in first_entry[g]:(first_entry[g + 1] - 1)
        entry_nodes[e] == assignment[g] && return e - first_entry[g] + 1
    end
end

function assign_tree(algorithm::AuctionAssignment, workspace, centroids, trees, groups, t, cost)
    (; kdtrees, k_groups, k_min, k_max, prices, assignment) = workspace
    n = size(groups, 2)
    # one group, or total cost 0 (every group is one point): nothing can be improved
    (n == 1 || iszero(cost)) && return groups[t, :]
    # final ε = auction_tol * cost / n -> error of the solve <= n * ε = auction_tol * cost
    ε_final = algorithm.auction_tol * cost / n
    fill!(prices, 0.0)

    while true
        # neighbours[g] - the nodes of tree t closest to the centroid of group g, sorted
        neighbours, _ = knn(kdtrees[t], centroids, maximum(view(k_groups, t, :)), true)
        n_entries = fill_auction_candidates!(workspace, neighbours, centroids, trees, groups, t)
        max_bids = algorithm.bids_per_entry * n_entries

        # too many bids: groups compete for too few candidates -> more candidates for all groups
        if !auction_assignment!(workspace, n, ε_final, algorithm.θ, max_bids)
            if all(==(k_max), view(k_groups, t, :))
                @warn "Auction for tree $t needed too many bids with k_max = $k_max candidates, its assignment is kept" maxlog = 1
                return groups[t, :]
            end
            grow_k!(k_groups, t, k_max)
            fill!(prices, 0.0)
            continue
        end

        # group g got its farthest candidate (rank k) or its current node outside its k closest
        # (rank k + 1): its list was probably too short, a farther node could be better ->
        # k of g is doubled and the problem is solved again (prices are kept)
        any_grown = false
        any_at_k_max = false
        for g in axes(groups, 2)
            if assigned_rank(workspace, g) >= k_groups[t, g]
                if k_groups[t, g] < k_max
                    k_groups[t, g] = min(2k_groups[t, g], k_max)
                    any_grown = true
                else
                    any_at_k_max = true
                end
            end
        end
        any_grown && continue
        # the result is still the best of the candidates, so it is never worse than the
        # current assignment (it is a candidate); with k_max = n all nodes are candidates
        any_at_k_max && k_max < n && @warn "k_max = $k_max candidates are not enough for some groups of tree $t, the result can be not optimal" maxlog = 1

        # group g got one of its k ÷ 4 closest nodes -> fewer candidates in the next round
        for g in axes(groups, 2)
            if k_groups[t, g] > k_min && assigned_rank(workspace, g) <= k_groups[t, g] ÷ 4
                k_groups[t, g] = max(k_min, k_groups[t, g] ÷ 2)
            end
        end
        return assignment
    end
end

"""
    alternating_assignment!(groups, trees, greedy_cost, options) -> cost

Improve the groups step by step. For every tree: take the centroids of the groups without this
tree and give the nodes of the tree to the groups again in the best way ([`assign_tree`](@ref)
with `options.algorithm`). A tree is solved again only if another tree's assignment changed
since its last solve. Stops when no tree has to be solved again, the cost improves by less than
`options.tol` (relative), or `options.max_iter` rounds are done. Returns the new total cost.
"""
function alternating_assignment!(groups, trees, greedy_cost, options::MatchingOptions)
    (; max_iter, tol, algorithm) = options
    cost = greedy_cost
    n = size(groups, 2)
    m = size(trees, 1)

    workspace = make_workspace(algorithm, trees, n)
    needs_solve = trues(size(trees, 1))   # needs_solve[t]: tree t has to be solved (again)
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
            # assignment[g] = node of tree t for group g
            assignment = assign_tree(algorithm, workspace, centroids, trees, groups, t, cost)
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
    find_closest_nodes!(coords, options::MatchingOptions) -> groups

Group the nodes of several trees so that each group has one node from every tree and the nodes
in a group are as close as possible. `coords[t, i]` holds the coordinates `(x, y, z)` of node
`i` of tree `t`. `options` (see [`MatchingOptions`](@ref)) sets the number of neighbours
checked in the greedy step, the algorithm and stopping rule of the alternating step and whether
the lower bound ([`calculate_lower_bound`](@ref)) is computed and it is logged how much the
result can still be improved at most; this is slow (m(m - 1)/2 Hungarian runs).

`coords` is reordered in place: afterwards column `g` holds the coordinates of the nodes of
group `g`, one node per tree. Returns `groups`: `groups[t, g]` is the original index of the node
of tree `t` in group `g`.
"""
function find_closest_nodes!(coords::Matrix, options::MatchingOptions)
    trees = SVector{3,Float64}.(coords)
    n_trees = size(coords, 1)
    n_nodes = size(coords, 2)
    groups = zeros(Int, n_trees, n_nodes)

    reference_tree = @view trees[1, :]

    _, greedy_cost = get_greedy_solution!(groups, reference_tree, trees, options.k_greedy)
    cost = alternating_assignment!(groups, trees, greedy_cost, options)

    if options.compute_bound
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

