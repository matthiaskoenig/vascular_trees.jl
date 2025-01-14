module Julia_models
    """
    FIXME: preterminals are wrong!!!!
    """
    function f_dxdt!(dx::Vector{Float64}, x::Vector{Float64}, p::Matrix{Float64}, t::Float64)

        # create views for convenience
        # p contains edges / elements (source id, target id), volumes, flows, 
        # is the edge from inflow system or not,
        # is it source for terminal node, ist it target of terminal node
        elements = @view p[:, 1:2]
        volume_values = @view p[:, 3]
        flow_values = @view p[:, 4]
        is_inflow = @view p[:, 5]
        source_for_terminal = @view p[:, 6]
        target_for_terminal = @view p[:, 7]

        @inbounds for (ke, element) ∈ enumerate(eachrow(elements))
            # retrieve information for element
            source_id, target_id = element
            element_volume = volume_values[ke]
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
                        dx[ke] += flow_values[ke] * x[pre_element] / element_volume
                    end
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flow_values[post_element] * x[ke] / element_volume
                end
    
            # inflow element connected to terminal
            elseif (element_source_for_terminal == 1.0 .&& element_is_inflow == 1.0 .&& length(pre_elements) != 0) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[ke] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / element_volume
                end
            
            # terminal element
            elseif (element_is_inflow == 0.0 .&& source_id == target_id)
                # dT/dt = QA * A / VT - QV * T / VT
    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[pre_element] * x[pre_element] / element_volume
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flow_values[post_element] * x[ke] / element_volume
                end
    
            # outflow element connected to terminal
            elseif (element_target_for_terminal == 0.0 .&& element_is_inflow == 0.0 .&& source_id != target_id)
                # dV/dt = QV * T / VT - QV * V / VV;
    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[ke] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / element_volume
                end
    
            # outflow element & outflow -> out not connected to terminal
            elseif (element_target_for_terminal == 0.0 .&& element_is_inflow == 0.0 .&& source_id != target_id)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[pre_element] * x[pre_element] / element_volume
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / element_volume
                end
    
            end
        end
        #@show dx
    end

    function f_dxdt_threaded!(dx::Vector{Float64}, x::Vector{Float64}, p::Matrix{Float64}, t::Float64)
        # create views for convenience
        # p contains edges / elements (source id, target id), volumes, flows, 
        # is the edge from inflow system or not,
        # is it source for terminal node, ist it target of terminal node
        elements = @view p[:, 1:2]
        volume_values = @view p[:, 3]
        flow_values = @view p[:, 4]
        is_inflow = @view p[:, 5]
        source_for_terminal = @view p[:, 6]
        target_for_terminal = @view p[:, 7]

        Threads.@threads for (ke, element) ∈ collect(enumerate(eachrow(elements)))
            # retrieve information for element
            # source_id, target_id = element
            # element_volume = volume_values[ke]
            # element_is_inflow = is_inflow[ke]
            # element_source_for_terminal = source_for_terminal[ke]
            # element_target_for_terminal = target_for_terminal[ke]

            # species before the element
            pre_elements[threadid()] = findall(x -> x==element[1], @view elements[:, 2]) 
            # species after the element
            post_elements[threadid()] = findall(x -> x==element[2], @view elements[:, 1]) 
    
            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (source_for_terminal[ke] == 0.0 .&& is_inflow[ke] == 1.0) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[ke] * x[pre_element] / volume_values[ke]
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flow_values[post_element] * x[ke] / volume_values[ke]
                end
    
            # inflow element connected to terminal
            elseif (source_for_terminal[ke] == 1.0 .&& is_inflow[ke] == 1.0 .&& length(pre_elements) != 0) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[ke] * x[pre_element] / volume_values[ke]
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / volume_values[ke]
                end
            
            # terminal element
            elseif (is_inflow[ke] == 0.0 .&& element[1] == element[2])
                # dT/dt = QA * A / VT - QV * T / VT
    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[pre_element] * x[pre_element] / volume_values[ke]
                end
                # species after
                for post_element ∈ post_elements
                    dx[ke] -= flow_values[post_element] * x[ke] / volume_values[ke]
                end
    
            # outflow element connected to terminal
            elseif (target_for_terminal[ke] == 0.0 .&& is_inflow[ke] == 0.0 .&& element[1] != element[2])
                # dV/dt = QV * T / VT - QV * V / VV;
    
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[ke] * x[pre_element] / volume_values[ke]
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / volume_values[ke]
                end
    
            # outflow element & outflow -> out not connected to terminal
            elseif (target_for_terminal[ke] == 0.0 .&& is_inflow[ke] == 0.0 .&& element[1] != element[2])
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[ke] += flow_values[pre_element] * x[pre_element] / volume_values[ke]
                end
                # species after
                for _ ∈ post_elements
                    dx[ke] -= flow_values[ke] * x[ke] / volume_values[ke]
                end
    
            end
        end
    end
end