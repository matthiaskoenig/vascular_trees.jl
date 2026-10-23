"""
    get_synchronization_indices(graph_id, species_ODE_groups)

Indices of the species of the tree `graph_id` that connect it with the terminal part:
preterminal species for inflow trees, terminal species for outflow trees.
`species_ODE_groups` holds the ODE group of every species of the tree.
"""
function get_synchronization_indices(graph_id, species_ODE_groups)
    if graph_id in flow_direction.inflow_trees
        return get_indices(species_ODE_groups, [groups.preterminal])
    else
        return get_indices(species_ODE_groups, [groups.terminal])
    end
end


"""
    get_indices(collection, targets)

Indices of the elements of `collection` that are in `targets`.
"""
function get_indices(collection, targets)
    return findall(in(targets), collection)
end
