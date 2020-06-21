using Random
using Plots
using DelimitedFiles


# Parametres
equib = 100_000
n = 1_000_000
dir = "ising-mod-data/"
Tc = 2 / (log(1 + sqrt(2)))

# number of different tempratures to simulate
temps = 81
Ts = LinRange(0.4*Tc, 1.4*Tc, temps)
# Ts[31] = Tc
# The different sizes of the grid to simutale
Ns = [2, 4, 8, 12, 16, 20]
sizes, = size(Ns)

names_obs = ["energy",  "abs_mag", "energy_squared", "abs_mag_squared", "fraction"]
names_func = ["tension", "heat_cap", "susc"]
num_quant, = size(names_obs)


# Utillities
# Small functions needed throught the program

function get_s(N)
    return rand!(Array{Int}(undef, (N, N)), [-1, 1])
end

function indxShft(i, n, N)
    """ Shift indexes with n steps, subject to peridic boundaries """
    return mod(i-1+n, N) + 1    # Fucking 1-indexation should be illegal
end

function sum_neigh(s, i, j, N)
    return (
        s[indxShft(i, -1, N), j] 
        + s[indxShft(i, +1, N), j] 
        + s[i, indxShft(j, -1, N)] 
        + s[i, indxShft(j, +1, N)]
    )
end


# Physics

function get_delta_H(s, i, j, N)
    return 2*s[i, j] * sum_neigh(s, i, j, N)
end

function energy(s, N)
    E = 0
    for i in 1:N
        for j in 1:N
            E += -1/2 * s[i, j] * sum_neigh(s, i, j, N)
        end
    end
    return E
end


 
function abs_mag(s, N)
    return abs(sum(s))
end


function fraction(s, N, T)
    diff_H = 2*sum(s[end, :].*s[1, :])
    return exp(-diff_H/T)
end

function tension(frac, N, T)
    return -T/N * log(frac)
end

function heat_cap(E, E2, N, T)
    return (E2 - E^2) / T^2 * N^2
end

function susc(M, M2, N, T)
    return (M2 - M^2) / T * N^2
end

function analytic_mag(Ts)
    return (1 - sinh(2/ Ts)^(-4))^(1 / 8)
end

# Monte Carlo simulation

function MC_sweep!(s, N, T)
    for _ in 1:N^2
        i, j = rand(1:N, 2)
        delta_H = get_delta_H(s, i, j, N)
        if delta_H < 0
            s[i, j] *= -1
            continue
        end
        if exp(-delta_H/T) > rand()
            s[i, j] *= -1
        end
    end
end

function get_sample_from_stat(s, N, T)
    observables = Array{Float64}(undef, (size(names_obs)[1]))
    observables[1] = energy(s, N)/ (N^2)
    observables[2] = abs_mag(s, N) / N^2
    observables[3] = observables[1]^2
    observables[4] = observables[2]^2
    observables[5] = fraction(s, N, T)

    return observables
end


function get_samples(N, T, n, equib)
    s = get_s(N)
    observables = zeros(size(names_obs)[1])
    for i in 1:equib
        MC_sweep!(s, N, T)
    end
    for i in 1:n
        MC_sweep!(s, N, T)
        observables += get_sample_from_stat(s, N, T)
    end
    return observables / n
end


function write_data(data, num_quant, names_obs, dir)
    for i in 1:num_quant
        path = dir * names_obs[i] * ".csv"
        file = open(path, "w")
        writedlm(file, data[i, :, :], ",")
        close(file)
    end
end

function read_data(num_quant, sizes, temps, names_obs, dir)
    data  = Array{Float64}(undef, (num_quant, sizes, temps))
    for i in 1:num_quant
        path = dir * names_obs[i] * ".csv"
        data[i, :, :] = readdlm(path, ',', Float64, '\n')
    end
    return data
end

function plot_equibliration()
    N = 100
    num_sweeps = 5000
    s = get_s(N)
    frames = 50
    heatmap(1:N, 1:N, s)
    heatmap!(show=true)

    T = 0.01
    for i in 1:num_sweeps
        MC_sweep!(s, N, T)
        if i%frames==1
            heatmap(1:N, 1:N, s)
            heatmap!(show=true)
        end
    end
end


function sample_observables()
    data  = Array{Float64}(undef, (num_quant, sizes, temps))
    steps = sizes * temps
    for i in eachindex(Ns)
        for j in eachindex(Ts)
            N = Ns[i]; T = Ts[j]
            data[:, i, j] = get_samples(N, T, n, equib)

            print(ceil(((i-1)*temps+j)/steps * 100), "%\n")
        end
    end
    write_data(data, num_quant, names_obs, dir)
end

function plot_observables()
    data = read_data(num_quant, sizes, temps, names_obs, dir)

    for i in eachindex(names_obs)
        p = plot()
        for j in eachindex(Ns)
            plot!(Ts, data[i, j, :], 
            marker=:circle, label=string(Ns[j]))
        end

        if i == 2
            T = LinRange(Ts[1], Tc, 500) 
            plot!(T, analytic_mag.(T),
            linewidth=3, linestyle=:dash, linecolor=:black) 
        end

        lim = [minimum(data[i, :, :]), maximum(data[i, :, :])]
        plot!([Tc, Tc], lim, label="\$T_c\$")
        plot!(title=names_obs[i])
        savefig(dir * names_obs[i] * ".png")
    end
end

function plot_tension()
    data = read_data(num_quant, sizes, temps, names_obs, dir)
    p = plot()
    for j in eachindex(Ns)
        tau = tension.(data[5, j, :], Ns[j], Ts)
        plot!(Ts, tau, marker=:dot, label=string(Ns[j]))
    end

    plot!([Tc, Tc], [0, 2])
    plot!(title=names_func[1])
    savefig(dir * names_func[1] * ".png")

end

function plot_heat_cap()
    data = read_data(num_quant, sizes, temps, names_obs, dir)
    p = plot()
    for j in eachindex(Ns)
        tau = heat_cap.(data[1, j, :], data[3, j, :], Ns[j], Ts)
        plot!(Ts, tau, marker=:dot, label=string(Ns[j]))
    end
    
    
    plot!([Tc, Tc], [0, 2], label="\$T_c\$")
    plot!(title=names_func[2])
    savefig(dir * names_func[2] * ".png")
end

function plot_susc()
    data = read_data(num_quant, sizes, temps, names_obs, dir)
    p = plot()
    for j in eachindex(Ns)
        tau = heat_cap.(data[2, j, :], data[4, j, :], Ns[j], Ts)
        plot!(Ts, tau, marker=:dot, label=string(Ns[j]))
    end
    plot!([Tc, Tc], [0, 4], label="\$T_c\$")
    plot!(title=names_func[3])
    savefig(dir * names_func[3] * ".png")
end

function more_plots()
    data = read_data(num_quant, sizes, temps, names_obs, dir)
    p = plot(size=(1000, 500))

    ts = (Tc .- Ts)./Ts
    m = 49
    for i in 1:sizes
        tau = tension.(data[5, i, m], Ns[i], Ts[m])
        plot!([Ns[i],], [tau,], marker=:dot, label=string(Ns[i]))
    end
    x = LinRange(log(2), log(20), 100)
    plot!(exp.(x), exp.(-x .+ 0.75),  xaxis=:log, yaxis=:log)
    plot!(title=string(Ts[m]))
    savefig(dir * "p1" * ".png")
    closeall()

    p = plot()
    indx = 40:48
    for i in indx
        t = ts[i]
        tau = tension.(data[5, :, i], Ns, Ts[i])
        plot!( 1 ./ (Ns .* t), tau ./ t, label=string(t))
    end
    
    x = LinRange(0, 40, 100)
    plot!(x, 2 .*x, linestyle=:dash, linecolor=:black)
    plot!(legend=:bottomright)
    savefig(dir * "p2" * ".png")
end

sample_observables()
plot_observables()
plot_tension()
plot_heat_cap()
plot_susc()
more_plots()
