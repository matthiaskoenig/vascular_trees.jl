module Julia_models

    export pyf_dxdt!, jf_dxdt!, str_jf_dxdt!

    function pyf_dxdt!(dx::Vector{Float64}, x::Vector{Float64}, p::Matrix{Float64}, t::Float64)
        """
        check ODEs! - double loops for terminal edges!
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
                        dx[ke] += flows[ke] * x[pre_element] / element_volume
                    end
                else
                    dx[ke] += flows[ke] * x[ke] / element_volume
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
                    dx[ke] += flows[ke] * x[pre_element] / element_volume
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
                    dx[ke] += flows[pre_element] * x[pre_element] / element_volume
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
                    dx[ke] += flows[ke] * x[pre_element] / element_volume
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
                    dx[ke] += flows[pre_element] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flows[ke] * x[ke] / element_volume
                end
    
            end
        end
        #@show dx
    end

    function jf_dxdt!(dx, 
                        x, 
                        p, 
                        t)

        edges = p[1]
        terminals = p[2]
        start = p[3]
        preterminals = p[4]
        flows = p[5]
        volumes = p[6]

        for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volumes[ke]
            element_flow = flows[ke]
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (!is_preterminal .&& !is_terminal) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                if length(pre_elements) != 0
                    for pre_element ∈ pre_elements
                        dx[ke] += element_flow * x[pre_element] / element_volume  
                    end
                else
                    dx[ke] += element_flow * x[ke] / element_volume  
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flows[post_element] * x[ke] / element_volume
                end

            # inflow element connected to terminal
            elseif (is_preterminal) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flows[ke] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flows[ke] * x[ke] / element_volume
                end

            # terminal element
            elseif (is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = QA * A / VT - QA * T / VT
    
                # species before and output
                for pre_element ∈ pre_elements
                    (!in(edges[pre_element], terminals)) && (dx[ke] += flows[pre_element] * x[pre_element] / element_volume - flows[pre_element] * x[ke] / element_volume)
                end
            end
        end
        #return dx
    end

    function str_jf_dxdt(dx, 
                        x, 
                        p, 
                        t)

        edges = p[1]
        terminals = p[2]
        start = p[3]
        preterminals = p[4]
        flow_ids = p[7]
        volume_ids = p[8]

        for (ke, element) in enumerate(edges)
            # retrieve information for element
            source_id, target_id = element
            element_volume = volume_ids[ke]
            element_flow = flow_ids[ke]
            is_preterminal = in(element, preterminals)
            is_terminal = in(element, terminals)

            # species before the element
            pre_elements = findall(x -> x[2]==source_id, edges) 
            # species after the element
            post_elements = findall(x -> x[1]==target_id, edges)

            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (!is_preterminal .&& !is_terminal) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                if length(pre_elements) != 0
                    for pre_element ∈ pre_elements
                        dx[ke] *= "+$element_flow * $(x[pre_element]) / $element_volume"  
                    end
                else
                    dx[ke] *= "+$element_flow * $(x[ke]) / $element_volume"
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] *= "-$(flow_ids[post_element]) * $(x[ke]) / $element_volume"
                end

            # inflow element connected to terminal
            elseif (is_preterminal) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] *= "+$(flow_ids[ke]) * $(x[pre_element]) / $element_volume"
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] *= "-$(flow_ids[ke]) * $(x[ke]) / $element_volume"
                end

            # terminal element
            elseif (is_terminal)
                # THIS IS DIFFERENT THAN IN OTHER MODELS
                # dT/dt = QA * A / VT - QA * T / VT
    
                # species before and output
                for pre_element ∈ pre_elements
                    (!in(edges[pre_element], terminals)) && (dx[ke] *= "+$(flow_ids[pre_element]) * $(x[pre_element]) / $element_volume - $(flow_ids[pre_element]) * $(x[pre_element]) / $element_volume")
                end
            end
        end
        return dx
    end

end