module Julia_models
    function f_dxdt!(dx::Vector{Float64}, x::Vector{Float64}, p::Matrix{Float64}, t::Float64)

        # create views for convenience
        # p contains edges / elements (source id, target id), volumes, flows, 
        # is the edge from inflow system or not
        elements = @view p[:, 1:2]
        volume_values = @view p[:, 3]
        flow_values = @view p[:, 4]
        is_inflow = @view p[:, 5]

        @inbounds for (ke, element) ∈ enumerate(elements)
            # retrieve information for element
            source_id, target_id = element
            element_index = source_id + (target_id - 1) * size(A, 1)
            #
            #
            #
            # species before the slement
            pre_elements = findall(x -> x==, @view A[:, source_id]) 
            # species after the slement
            post_elements = findall(!iszero, @view A[target_id, :])

            # distinguish elements connected to terminal nodes
            terminal_source::Bool = false
            if (length(post_elements) == 1) && (target_id == post_elements[1])
                terminal_source = true
            end
            terminal_target::Bool = false
            if (length(pre_elements) == 1) && (source_id == pre_elements[1])
                terminal_target = true
            end

            # write equations
            # marginal inflow element & in -> inflow & inflow element not connected to terminal
            if (terminal_source == false .&& is_inflow[element] == 1.0) 
                # dA_marginal/dt = -QAmarginal * A_marginal / VAmarginal;
                # or
                # dA/dt = QA * A_pre / VA - Q_post * A / VA;
                # or
                # dA/dt = QA * A_marginal / VA - Q_post * A / VA; (QA and VA here are equal to QAmarginal and VAmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end

            # inflow element connected to terminal
            elseif (terminal_source == true .&& is_inflow[element] == 1.0 .&& length(pre_elements) != 0) 
                # dA/dt = QA * A_pre / VA - QA * A / VA;
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end
            
            # terminal element
            elseif (is_inflow[element] == 0.0 .&& source_id == target_id)
                # dT/dt = QA * A / VT - QV * T / VT

                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for post_element ∈ post_elements
                    dx[element_index] -= flow_values[target_id, post_element] * x[element] / volume_values[element]
                end

            # outflow element connected to terminal
            elseif (terminal_target == true .&& is_inflow[element] == 0.0 .&& source_id != target_id)
                # dV/dt = QV * T / VT - QV * V / VV;

                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[element] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            # outflow element & outflow -> out not connected to terminal
            elseif (terminal_target == false .&& is_inflow[element] == 0.0 .&& source_id != target_id)
                # dV/dt = QV_pre * V_pre / VV - QV * V / VV;
                # or
                # dV_marginal = QV_pre * V_pre / VVmarginal - QV_marginal * V_marginal / VV_marginal; (QV and VV here are equal to QVmarginal and VVmarginal - done in processing julia graph)
                
                # species before
                for pre_element ∈ pre_elements
                    dx[element_index] += flow_values[pre_element, source_id] * x[pre_element, source_id] / volume_values[element]
                end
                # species after
                for _ ∈ post_elements
                    dx[element_index] -= flow_values[element] * x[element] / volume_values[element]
                end

            end

        end
        #@show dx
    end
end