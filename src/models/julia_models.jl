module Julia_models

    export jf_dxdt!

    function jf_dxdt!(dx,
                        x, 
                        p, 
                        t)

        is_inflow = p[2]
        edges = p[3]
        terminals = p[4]
        #preterminalssrc/models/julia_models.jl
        flows=p[5]
        volumes=p[6]
        groups=p[7]
        pre_element = p[8]
        post_element = p[9]
        
        if is_inflow
            jf_inflow!(dx, x, edges, terminals, flows, volumes, groups, pre_element, post_element, t)
        else
            jf_outflow!(dx, x, edges, terminals, flows, volumes, groups, t)
        end
        # println()
        # println(size(dx))
        # show(dx)
        # println()
    end

    function jf_inflow!(dx,
                        x,
                        edges,
                        terminals,
                        flows,
                        volumes,
                        groups,
                        pre_element,
                        post_element,
                        t)
                        # (dx::Vector{Float64},
                        # x::Vector{Float64},
                        # edges::Vector{Tuple{Int32, Int32}},
                        # terminals::Vector{Tuple{Int32, Int32}},
                        # flows::Vector{Float64},
                        # volumes::Vector{Float64},
                        # groups::Vector{Int16},
                        # t::Union{Tuple{Float64, Float64}, Float64})

        #x[length(x)] = f_intervention(t)

        for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            # is_preterminal = in(element, preterminals)
            # is_terminal = in(element, terminals)
            group = groups[ke]
            pre_elements = pre_element[ke]
            post_elements = post_element[ke]

            # species before the element
            # pre_elements = findfirst(x -> x[2]==source_id, edges) 
            # species after the element
            #post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal (input) element
            if group == 0 # source_id == 0
                # constant
                dx[ke] = 0.0
            # in -> inflow & inflow element not connected to terminal
            elseif group == 1 #(!is_preterminal .&& !is_terminal .&& source_id != 0) 
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                    
                # species before
                #for pre ∈ pre_elements
                dx[ke] = sum(element_flow .* x[pre_elements])
                #end
                # species after
                #for post ∈ post_elements
                dx[ke] = dx[ke] - sum(flows[post_elements] .* x[ke])
                #end

            # inflow element connected to terminal
            elseif group == 2 # (is_preterminal) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                #for pre ∈ pre_elements
                dx[ke] = sum(element_flow .* x[pre_elements])
                #end
                # species after
                #for _ ∈ post_elements
                dx[ke] -= element_flow * x[ke] * length(post_elements)
                #end

            # terminal element
            elseif group == 3 # (is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = QA * A / VT - QA * T / VT
    
                # species before and output
                #for pre ∈ pre_elements
                dx[ke] = sum(flows[pre_elements] .* x[pre_elements]) - sum(flows[pre_elements] .* x[ke])
                #end
            end

            if group != 0
                dx[ke] = dx[ke] / element_volume
            end
        end

    end

    function jf_outflow!(dx::Vector{Float64},
                        x::Vector{Float64},
                        edges::Vector{Tuple{Int32, Int32}},
                        terminals::Vector{Tuple{Int32, Int32}},
                        flows::Vector{Float64},
                        volumes::Vector{Float64},
                        groups::Vector{Int16},
                        t::Union{Tuple{Float64, Float64}, Float64})

        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            # is_preterminal = in(element, preterminals)
            # is_terminal = in(element, terminals)
            # is_start = in(element, start)
            group = groups[ke]

            dx[ke] = 0.0

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal (outflow) element
            if group == 0 #target_id == 0
                # species before
                #for pre_element ∈ pre_elements
                dx[ke] = sum((flows[pre_elements] .* x[pre_elements] .- flows[pre_elements] .* x[pre_elements]) ./ volumes[pre_elements])
                #end
            # outflow -> out & outflow element not connected to terminal
            elseif group == 1 #(!is_preterminal .&& !is_terminal .&& target_id != 0) 
                # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
                    
                # species before
                #for pre_element ∈ pre_elements
                dx[ke] += sum(flows[pre_elements] .* x[pre_elements]) #flows[pre_element] * x[pre_element]
                #end

                # species after
                #for _ ∈ post_elements
                dx[ke] -= element_flow .* x[ke] .* length(post_elements)
                #end

            # outflow element connected to terminal
            elseif group == 2 #(is_preterminal) 
                # dV/dt = QV * T / VV - QV * V / VV;
                
                # species before
                #for pre_element ∈ pre_elements
                dx[ke] += sum(element_flow .* x[pre_elements])
                #end
                # species after
                #for _ ∈ post_elements
                dx[ke] -= element_flow .* x[ke] .* length(post_elements)
                #end

            # terminal element
            elseif group == 3 #(is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = 0
                # species before and output
                for pre_element ∈ pre_elements
                    (!in(edges[pre_element], terminals)) && (dx[ke] = 0.0)
                end
            end

            if target_id != 0
                dx[ke] = dx[ke] / element_volume
            end
        end
    end

end