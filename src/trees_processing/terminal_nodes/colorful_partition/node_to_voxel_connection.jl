function connect_nodes_to_voxels(grid, P, inside)

    owner = assign_owners(inside, grid, P)

    outside = Int[]
    for i in axes(P, 2)
        v = voxel_of(grid, view(P, :, i))
        if !checkbounds(Bool, inside, v...) || !inside[v...]
            push!(outside, i)
        end
    end
    !isempty(outside) && @warn "$(length(outside)) data points lie outside V (or within one voxel of its surface)"

    merge_stray_pieces!(owner, grid, P)
    svolume = small_cell_volumes(owner, grid.h, size(P, 2))

    pairs = find_neighbor_pairs(owner)
    orphans = attach_orphans!(pairs, svolume, P)  
    offsets, nbrs = create_adjacency_vector(pairs, size(P, 2))

    return (owner = owner, outside = outside, svolume = svolume,
        orphans = orphans, offsets = offsets, nbrs = nbrs)
end


find_voxel_center(g::VoxelGrid, i, j, k) =
    (g.origin[1] + (i - 0.5) * g.h,
     g.origin[2] + (j - 0.5) * g.h,
     g.origin[3] + (k - 0.5) * g.h)

"""
    assign_owners(inside, g, P) -> owner

`owner[i, j, k]` is the index of the data point nearest to the center of voxel
(i, j, k), or 0 if the voxel is outside V. Works slice by slice along z, one
slice per thread, so start Julia with several threads (`julia -t auto`).
"""
function assign_owners(inside::BitArray{3}, g::VoxelGrid, P::AbstractMatrix{Float64})
    nx, ny, nz = g.dims
    owner = zeros(Int32, nx, ny, nz)
    tree = KDTree(P)

    Threads.@threads for k in 1:nz
        n_in = count(@view inside[:, :, k])
        n_in == 0 && continue
        centers = Matrix{Float64}(undef, 3, n_in)
        where_ = Vector{CartesianIndex{2}}(undef, n_in)
        q = 0
        for j in 1:ny, i in 1:nx
            inside[i, j, k] || continue
            q += 1
            x, y, z = find_voxel_center(g, i, j, k)
            centers[1, q] = x
            centers[2, q] = y
            centers[3, q] = z
            where_[q] = CartesianIndex(i, j)
        end
        idxs, _ = nn(tree, centers)
        for q in 1:n_in
            owner[where_[q][1], where_[q][2], k] = Int32(idxs[q])
        end
    end
    return owner
end

"""
    voxel_of(g, p) -> (i, j, k)

Voxel that contains the location `p` (may lie outside the grid).
"""
voxel_of(g::VoxelGrid, p) =
    ntuple(d -> floor(Int, (p[d] - g.origin[d]) / g.h) + 1, 3)


"""
    small_cell_volumes(owner, h, N) -> vol

vol[i] = number of voxels owned by point i, times the volume of one voxel.
"""
function small_cell_volumes(owner::Array{Int32,3}, h, N)
    counts = zeros(Int, N)
    for o in owner
        o == 0 || (counts[o] += 1)     # 0 = voxel outside V, skip it
    end
    return counts .* h^3
end

"""
    label_pieces(owner) -> (comp, comp_owner, comp_size)

comp[v] is the number of the piece that voxel v belongs to (0 = outside V).
A piece is a group of voxels with the same owner, connected through shared faces.
"""
function label_pieces(owner::Array{Int32,3})
    comp = zeros(Int32, size(owner))
    comp_owner = Int32[]
    comp_size = Int[]
    stack = CartesianIndex{3}[]
    for v in CartesianIndices(owner)
        (owner[v] == 0 || comp[v] != 0) && continue
        push!(comp_owner, owner[v])
        push!(comp_size, 0)
        id = Int32(length(comp_owner))
        comp[v] = id
        push!(stack, v)
        while !isempty(stack)
            u = pop!(stack)
            comp_size[id] += 1
            for d in NB6
                w = u + d
                if checkbounds(Bool, owner, w) && comp[w] == 0 && owner[w] == owner[u]
                    comp[w] = id
                    push!(stack, w)
                end
            end
        end
    end
    return comp, comp_owner, comp_size
end

"""
    merge_stray_pieces!(owner, g, P; max_passes = 10)

Makes every point's small cell a single piece: each point keeps the piece it sits
in, and every other piece goes to the neighboring point whose home piece it shares
the most faces with. Changes `owner` in place.
"""
function merge_stray_pieces!(owner::Array{Int32,3}, g::VoxelGrid, P::AbstractMatrix;
                             max_passes = 10)
    N = size(P, 2)
    for _ in 1:max_passes
        comp, comp_owner, comp_size = label_pieces(owner)

        # Home piece: the piece containing the point's own voxel...
        home = zeros(Int32, N)
        for i in 1:N
            v = CartesianIndex(voxel_of(g, view(P, :, i)))
            if checkbounds(Bool, owner, v) && owner[v] == i
                home[i] = comp[v]
            end
        end
        # ...or, if another point owns that voxel, the point's largest piece.
        found = home .!= 0            # true: home piece came from the point's own voxel
        for c in eachindex(comp_owner)
            i = comp_owner[c]
            if !found[i] && (home[i] == 0 || comp_size[c] > comp_size[home[i]])
                home[i] = c
            end
        end

        is_home = falses(length(comp_owner))
        for i in 1:N
            home[i] != 0 && (is_home[home[i]] = true)
        end
        n_stray = count(!, is_home)
        n_stray == 0 && return owner

        # Shared faces between each stray piece and the home pieces next to it.
        contacts = Dict{Int32,Dict{Int32,Int}}()
        for v in CartesianIndices(owner)
            c = comp[v]
            (c == 0 || is_home[c]) && continue
            for d in NB6
                w = v + d
                checkbounds(Bool, owner, w) || continue
                cw = comp[w]
                (cw == 0 || cw == c || !is_home[cw]) && continue
                tally = get!(contacts, c, Dict{Int32,Int}())
                tally[owner[w]] = get(tally, owner[w], 0) + 1
            end
        end

        # Give each stray piece to the point it touches most.
        new_owner = Dict(c => argmax(tally) for (c, tally) in contacts)
        isempty(new_owner) && (@warn "$n_stray pieces touch no other cell (parts of V with no data point)"; return owner)
        for v in CartesianIndices(owner)
            c = comp[v]
            c != 0 && haskey(new_owner, c) && (owner[v] = new_owner[c])
        end
    end
    @warn "stray pieces remain after $max_passes passes"
end

"""
    find_neighbor_pairs(owner) -> pairs

All pairs (a, b), a < b, of points whose voxels share at least one face.
"""
function find_neighbor_pairs(owner::Array{Int32,3})
    nx, ny, nz = size(owner)
    pairs = Tuple{Int32,Int32}[]
    for k in 1:nz, j in 1:ny, i in 1:nx
        a = owner[i, j, k]
        a == 0 && continue
        if i < nx
            b = owner[i+1, j, k]
            (b != 0 && b != a) && push!(pairs, minmax(a, b))
        end
        if j < ny
            b = owner[i, j+1, k]
            (b != 0 && b != a) && push!(pairs, minmax(a, b))
        end
        if k < nz
            b = owner[i, j, k+1]
            (b != 0 && b != a) && push!(pairs, minmax(a, b))
        end
    end
    sort!(pairs)
    unique!(pairs)
    return pairs
end

"""
    create_adjacency_vector(pairs, N) -> (offsets, nbrs)

The neighbors of point i are nbrs[offsets[i] : offsets[i+1]-1].
"""
function create_adjacency_vector(pairs, N)
    deg = zeros(Int, N)
    for (a, b) in pairs
        deg[a] += 1
        deg[b] += 1
    end
    offsets = Vector{Int}(undef, N + 1)
    offsets[1] = 1
    for i in 1:N
        offsets[i+1] = offsets[i] + deg[i]
    end
    nbrs = Vector{Int32}(undef, offsets[end] - 1)
    pos = offsets[1:N]
    for (a, b) in pairs
        nbrs[pos[a]] = b; pos[a] += 1
        nbrs[pos[b]] = a; pos[b] += 1
    end
    return offsets, nbrs
end

neighbors(offsets, nbrs, i) = view(nbrs, offsets[i]:offsets[i+1]-1)

edge_length(P, i, j) = sqrt((P[1,i]-P[1,j])^2 + (P[2,i]-P[2,j])^2 + (P[3,i]-P[3,j])^2)

"""
    attach_orphans!(pairs, vol, P; k = 16) -> orphans

Every point with no voxels (vol == 0) gets one edge to the nearest point that has
some. Adds those edges to `pairs` and returns the list of such points.
"""
function attach_orphans!(pairs, vol, P; k = 16)
    tree = KDTree(P)
    orphans = findall(==(0), vol)
    for i in orphans
        idxs, _ = knn(tree, P[:, i], k, true)      # k nearest points, closest first
        q = findfirst(j -> j != i && vol[j] > 0, idxs)
        if q === nothing
            @warn "point $i: none of its $k nearest points has any voxels"
            continue
        end
        push!(pairs, minmax(Int32(i), Int32(idxs[q])))
    end
    sort!(pairs)
    unique!(pairs)
    return orphans
end