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