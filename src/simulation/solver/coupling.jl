"""
    get_synchronization_indices(graph_id, species_ODE_groups, perm)

Indices of the species of the tree `graph_id` that connect it with the terminal part:
preterminal species for inflow trees, terminal species for outflow trees.
`species_ODE_groups` holds the ODE group of every species of the tree.

`perm` is the row of `TerminalParameters.synch_indices` that belongs to this tree. The indices
are reordered with it, so that element k is the species connected to terminal k (column k of
the terminal block).
"""
function get_synchronization_indices(graph_id, species_ODE_groups, perm)
    group = graph_id in flow_direction.inflow_trees ? groups.preterminal : groups.terminal
    return get_indices(species_ODE_groups, [group])[perm]
end


"""
    get_indices(collection, targets)

Indices of the elements of `collection` that are in `targets`.
"""
function get_indices(collection, targets)
    return findall(in(targets), collection)
end
