module Interventions

export f_intervention

const dose = 1.0
const hill = 1.0
const stop_time = 7.0

function f_intervention(t)
    # return f_saturation(t)
    return f_soft_step(t)
end

f_saturation(t) = dose * (t^hill / (0.1^hill + t^hill)) # amount per time = mmole / min
f_decline(t)= dose * (-(t-stop_time)^hill / (0.1^hill + (t-stop_time)^hill) + 1.0) 

function f_soft_step(t)

    if t < stop_time
        f = t -> f_saturation(t)
    elseif t >= stop_time
        f = t -> f_decline(t)
    end
    return f(t)
end

end