module Julia_models

    export pyf_dxdt!, jf_dxdt!, str_jf_dxdt!
    

    function pyf_dxdt!(dx::Vector{Float64}, x::Vector{Float64}, p::Matrix{Float64}, t::Float64)
        """
        """
        # create views for convenience
        # p contains edges / elements (source id, target id), volumes, flows, 
        # is the edge from inflow system or not,
        # is it source for terminal node, ist it target of terminal node
        elements = @view p[:, 1:2]
        volumes = @view p[:, 3]
        flows = @view p[:, 4]
        is_inflow = @view p[:, 5]
        source_for_terminal = @view p[:, 6]
        target_for_terminal = @view p[:, 7]

        @inbounds for (ke, element) ∈ enumerate(eachrow(elements))
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_is_inflow = is_inflow[ke]
            element_source_for_terminal = source_for_terminal[ke]
            element_target_for_terminal = target_for_terminal[ke]
        
            # species before the element
            pre_elements = findall(x -> x==source_id, @view elements[:, 2]) 
            # species after the element
            post_elements = findall(x -> x==target_id, @view elements[:, 1]) 
        
            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (element_source_for_terminal == 0.0 .&& element_is_inflow == 1.0) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                if length(pre_elements) != 0
                    for pre_element ∈ pre_elements
                        dx[ke] = flows[ke] * x[pre_element] / element_volume
                    end
                else
                    dx[ke] = flows[ke] * x[ke] / element_volume
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flows[post_element] * x[ke] / element_volume
                end
    
            # inflow element connected to terminal
            elseif (element_source_for_terminal == 1.0 .&& element_is_inflow == 1.0 .&& length(pre_elements) != 0) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] = flows[ke] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flows[ke] * x[ke] / element_volume
                end
            
            # terminal element
            elseif (element_is_inflow == 0.0 .&& source_id == target_id)
                # dT/dt = QA * A / VT - QV * T / VT
    
                # species before
                for pre_element ∈ pre_elements
                    if elements[pre_element][1] != elements[pre_element][2]
                        dx[ke] = flows[pre_element] * x[pre_element] / element_volume
                    end
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flows[post_element] * x[ke] / element_volume
                end
    
            # outflow element connected to terminal
            elseif (element_target_for_terminal == 0.0 .&& element_is_inflow == 0.0 .&& source_id != target_id)
                # dV/dt = QV * T / VT - QV * V / VV;
    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] = flows[ke] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flows[ke] * x[ke] / element_volume
                end
    
            # outflow element & outflow -> out not connected to terminal
            elseif (element_target_for_terminal == 0.0 .&& element_is_inflow == 0.0 .&& source_id != target_id)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] = flows[pre_element] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flows[ke] * x[ke] / element_volume
                end
    
            end
        end
        #@show dx
    end

    function jf_dxdt!(dx::Vector{Float64}, 
                        x::Vector{Float64}, 
                        p::Tuple{Bool, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Tuple{Int32, Int32}}, Vector{Float64}, Vector{Float64}}, 
                        t::Union{Tuple{Float64, Float64}, Float64})

        is_inflow = p[1]
        edges = p[2]
        terminals = p[3]
        preterminals = p[4]
        flows = p[5]
        volumes = p[6]
        
        if is_inflow
            jf_inflow!(dx, x, edges, terminals, preterminals, flows, volumes, t)
        else
            jf_outflow!(dx, x, edges, terminals, preterminals, flows, volumes, t)
        end
        # println()
        # println(size(dx))
        # show(dx)
        # println()
    end

    function jf_inflow!(dx::Vector{Float64},
                        x::Vector{Float64},
                        edges::Vector{Tuple{Int32, Int32}},
                        terminals::Vector{Tuple{Int32, Int32}},
                        preterminals::Vector{Tuple{Int32, Int32}},
                        flows::Vector{Float64},
                        volumes::Vector{Float64},
                        t::Union{Tuple{Float64, Float64}, Float64})

        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)

            dx[ke] = 0.0

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal (input) element
            if source_id == 0
                # constant
                dx[ke] = 0.0
            # in -> inflow & inflow element not connected to terminal
            elseif (!is_preterminal .&& !is_terminal .&& source_id != 0) 
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += element_flow * x[pre_element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flows[post_element] * x[ke]
                end

            # inflow element connected to terminal
            elseif (is_preterminal) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += element_flow * x[pre_element]
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= element_flow * x[ke]
                end

            # terminal element
            elseif (is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = QA * A / VT - QA * T / VT
    
                # species before and output
                for pre_element ∈ pre_elements
                    (!in(edges[pre_element], terminals)) && (dx[ke] = flows[pre_element] * x[pre_element] - flows[pre_element] * x[ke])
                end
            end

            if source_id != 0
                dx[ke] = dx[ke] / element_volume
            end
        end

    end

    function jf_outflow!(dx::Vector{Float64},
                        x::Vector{Float64},
                        edges::Vector{Tuple{Int32, Int32}},
                        terminals::Vector{Tuple{Int32, Int32}},
                        preterminals::Vector{Tuple{Int32, Int32}},
                        flows::Vector{Float64},
                        volumes::Vector{Float64},
                        t::Union{Tuple{Float64, Float64}, Float64})

        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)
            # is_start = in(element, start)

            dx[ke] = 0.0

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal (outflow) element
            if target_id == 0
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flows[pre_element] * x[pre_element]
                end
            # outflow -> out & outflow element not connected to terminal
            elseif (!is_preterminal .&& !is_terminal .&& target_id != 0) 
                # dV/dt = QV_pre * V_pre / VV - Q * V / VV;
                    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flows[pre_element] * x[pre_element]
                end

                # species after
                for _ ∈ post_elements
                    dx[ke] -= element_flow * x[ke]
                end

            # outflow element connected to terminal
            elseif (is_preterminal) 
                # dV/dt = QV * T / VV - QV * V / VV;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += element_flow * x[pre_element]
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= element_flow * x[ke]
                end

            # terminal element
            elseif (is_terminal)
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

    function str_jf_dxdt(dx_str,
                        dx_vstr, 
                        x, 
                        p, 
                        t)

        edges = p[3]
        terminals = p[4]
        start = p[5]
        preterminals = p[6]

        flows = p[7]
        volumes = p[8]

        flow_ids = p[9]
        volume_ids = p[10]


        @inbounds for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            element_volume_id = volume_ids[ke]
            element_flow_id = flow_ids[ke]
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)
            is_start = in(element, start)

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal (input) element
            if source_id == 0
                # species after
                # species after
                # dx[ke] = 0
                dx_str[ke] *= "0.0" 
                dx_vstr[ke] *= "0.0"
            # in -> inflow & inflow element not connected to terminal
            elseif (!is_preterminal .&& !is_terminal .&& source_id != 0) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;

                # species before
                for pre_element ∈ pre_elements
                    # dx[ke] += element_flow * x[pre_element]
                    dx_str[ke] *= " + $element_flow_id*$(x[pre_element])" 
                    dx_vstr[ke] *= " + $element_flow*$(x[pre_element])"
                end

                # species after
                for post_element ∈ post_elements
                    # dx[ke] -= flows[post_element] * x[ke]
                    dx_str[ke] *= " - $(flow_ids[post_element])*$(x[ke])"
                    dx_vstr[ke] *= " - $(flows[post_element])*$(x[ke])"
                end

            # inflow element connected to terminal
            elseif (is_preterminal) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    # dx[ke] += flows[ke] * x[pre_element]
                    dx_str[ke] *= " + $(flow_ids[ke])*$(x[pre_element])"
                    dx_vstr[ke] *= " + $(flows[ke])*$(x[pre_element])"
                end
                # species after
                for _ ∈ post_elements
                    # dx[ke] -= flows[ke] * x[ke]
                    dx_str[ke] *= " - $(flow_ids[ke])*$(x[ke])"
                    dx_vstr[ke] *= " - $(flows[ke])*$(x[ke])"
                end

            # terminal element
            elseif (is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = QA * A / VT - QA * T / VT
    
                # species before and output
                for pre_element ∈ pre_elements
                    # (!in(edges[pre_element], terminals)) && (dx[ke] += flows[pre_element] * x[pre_element])
                    if !in(edges[pre_element], terminals)
                        dx_str[ke] *= " + $(flow_ids[pre_element])*$(x[pre_element])" # - $(flow_ids[pre_element])*$(x[ke])/$element_volume_id")
                        dx_vstr[ke] *= " + $(flows[pre_element])*$(x[pre_element])"
                    end
                end
            end

            # marginal (input) element
            if source_id != 0
                dx_str[ke] = "(" * dx_str[ke] * " ) / $element_volume_id"
                dx_vstr[ke] = "(" * dx_vstr[ke] * " ) / $element_volume"
            end
        end
        return dx_str, dx_vstr
    end

end