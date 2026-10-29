using NearestNeighbors,
    StaticArrays

function calculate_squared_distance(coordinates)
    x_coord = [c[1] for c in coordinates]
    y_coord = [c[2] for c in coordinates]
    z_coord = [c[3] for c in coordinates]

    x_centroid = sum(x_coord) / length(x_coord)
    y_centroid = sum(y_coord) / length(y_coord)
    z_centroid = sum(z_coord) / length(z_coord)

    squared_dist = []
    for c in coordinates
        x, y, z = c
        distance = (x - x_centroid)^2 + (y - y_centroid)^2 + (z - z_centroid)^2
        push!(squared_dist, distance)
    end
    return sum(squared_dist)
end

function find_closest_nodes(nodes::Matrix, k::Int)
    fixed_nodes = @view nodes[1, :]
    fixed = reduce(hcat, collect.(nodes[1, :]))

    n_trees = size(nodes, 1) - 1
    kdtrees = Vector{KDTree}(undef, n_trees)
    knn_idx = Vector{Vector{Vector{Int}}}(undef, n_trees)
    knn_dist = Vector{Vector{Vector{Float64}}}(undef, n_trees)
    for (kr, row) in enumerate(eachrow(@view nodes[2:end, :]))
        data = reduce(hcat, collect.(row))
        kdtrees[kr] = KDTree(data)
        knn_idx[kr], knn_dist[kr] = knn(kdtrees[kr], fixed, min(k, size(data, 2)))
    end

    sq_distances = []
    candidates = []
    for (kn, c) in enumerate(fixed_nodes)
        # for this node collect all neighbour nodes from all trees 
        neighbour_lists = [knn_idx[kr][kn] for kr in 1:n_trees]
        
        # get all possible combinations of this node and neighbours to check them
        combinations = Iterators.product(neighbour_lists...)
        for combo in combinations
            coords = [c]
            for (kr, idx) in enumerate(combo)
                push!(coords, nodes[kr + 1, idx])
            end
            push!(candidates, (kn, combo...))
            push!(sq_distances, calculate_squared_distance(coords))
        end
    end

    p = sortperm(sq_distances, rev=false)
    used = [falses(size(nodes, 2)) for _ in 1:n_trees + 1] # one flag vector per row
    triplets = []
    for idx in p
        cand = candidates[idx]
        any(used[r][cand[r]] for r in eachindex(cand)) && continue

        push!(triplets, cand)
        for r in eachindex(cand)
            used[r][cand[r]] = true # consume the points
        end
    end
    
    leftovers_fixed = findall(.!used[1])
    leftovers_b = findall(.!used[2])
    leftover = [leftovers_fixed, leftovers_b]
    sq_distances_l = []
    combinations = Iterators.product(leftover...)
    for (kc, combo) in enumerate(combinations)
        a = combo[1]
        b = combo[2]
        coords = [nodes[1, a], nodes[2, b]]
        push!(sq_distances_l, calculate_squared_distance(coords))
    end
    combinations = collect(Iterators.product(leftover...)) 
    p = sortperm(sq_distances_l, rev=false)
    used_l = [falses(size(nodes, 2)) for _ in 1:2 + 1]
    triplets_l = []
    for idx in p
        cand = combinations[idx]
        any(used_l[r][cand[r]] for r in eachindex(cand)) && continue

        push!(triplets_l, cand)
        for r in eachindex(cand)
            used_l[r][cand[r]] = true # consume the points
        end
    end

    # --- attach leftover C points to the leftover A–B pairs ---
    midpoint(p, q) = (p .+ q) ./ 2
    sq_dist(p, q) = sum((p .- q) .^ 2)

    leftovers_c = findall(.!used[3])
    midpoints = [midpoint(nodes[1, a], nodes[2, b]) for (a, b) in triplets_l]

    # all (pair index, c index) combinations, scored by distance from c to the pair's midpoint
    pair_c = collect(Iterators.product(eachindex(triplets_l), leftovers_c))
    sq_distances_c = [sq_dist(midpoints[ip], nodes[3, c]) for (ip, c) in pair_c]

    p = sortperm(vec(sq_distances_c))        # smallest distance first
    used_pairs = falses(length(triplets_l))
    used_c = falses(size(nodes, 2))
    triplets_abc = []
    for idx in p
        ip, c = pair_c[idx]
        (used_pairs[ip] || used_c[c]) && continue

        a, b = triplets_l[ip]
        push!(triplets_abc, (a, b, c))
        used_pairs[ip] = true                # consume the pair and the point
        used_c[c] = true
    end

end

"""
End of Step 1 (leftovers only). After the greedy pass, a few points remain unused. Only those are matched: the remaining A points with the remaining B points, then the remaining C points attached to the midpoints of those pairs. The greedy triplets stay as they are. The purpose is simply to make the solution complete.

Step 2 (everything). Starting from that complete solution, each step takes all n pairs and reassigns all n points of the third set at once. For example, it fixes all n A–B pairs and finds the best assignment of all n C points to their midpoints. Every triplet can change, including the ones greedy picked as excellent.

Why it has to be global: the leftover triplets are usually poor because their points are far apart. They can’t be fixed using only each other, since no good partners exist among the leftovers. The fix is for nearby good triplets to give up a point and take a slightly worse one in return. A good triplet gets a little looser, while a bad one gets much tighter, and the total drops. That’s only possible if Step 2 is allowed to change the good triplets too.
"""