"""
    wetted!(area_size, drop_pos, maxheight, height, t; hthresh = 0.055)

Measures the wetted area and maximal height of the film at time step `t`. 
"""
function wetted!(area_size, maxheight, height, t; hthresh = 0.055)
    area_size[t] = length(findall(height .> hthresh))
    maxheight[t] = maximum(height)

    return nothing
end

function wetted!(area_size, state::LBM_state; hthresh = 0.055)
    push!(area_size, length(findall(state.height .> hthresh)))

    return nothing
end
"""
    fluid_dry!(fluid, dummy, height, t; hthresh = 0.055)

Tracks the location of the thin film as a boolean field.
"""
function fluid_dry!(fluid, dummy, height, t; hthresh = 0.055)
    dummy .= false
    dummy[height .> hthresh] .= true
    fluid[t, :] .= vec(dummy)
    
    return nothing
end


"""
    contact_line(state::StateMultiLayer_1D; th1=0.05, th2=0.052, j=1, side=1)

Detects the contact line of a droplet. `side=1`` is the left contact line, `side=-1` the right one. 
"""
function contact_line(state::StateMultiLayer_1D; th1=0.05, th2=0.052, j=1, side=1)
    input=state.height
	start=side==1 ? 1 : length(input[:,j])
	target=side==1 ? length(input[:,j]) : 1
	height=input[:,j]
	for i in start:side:target
		if height[i]<th1
			start=i
			break
		end
	end
	contact_line=1
	for i in start:side:target
		if height[i]>th2
			contact_line=i
			return contact_line
			break
		end
	end
	return NaN
end			

"""
	function dewetted!(state:: LBM_state, sys::SysConstActive_1D)

Dewtermines the area of the substrate tha is dewetted
"""
function dewetted!(state:: LBM_state, sys::SysConstActive_1D)
    hthresh = sys.hmin-sys.hcrit + 0.005
    state.precursor .= 0
    state.precursor[state.height .< hthresh] .=1
end

"""
	function rupture_time!(rupt, state, t; htresh=0.06)

overwrites `rupt` with `t` if the film is ruptured somewhere
"""
function rupture_time!(rupt, state, t, sys)
    hthresh = sys.hmin-sys.hcrit + 0.005
	if minimum(state.height)<hthresh && isnan(rupt[1])
        rupt.=[t,true]
    end
    return nothing
end


"""
    t0(;hᵦ=0.07, γ=0.01, μ=1/6, θ=1/6)

Computes a characteristic time scale for an spinodally dewetting film.

# Arguments

- `hᵦ :: Float64`: height at which the disjoining pressure vanishes
- `γ :: Float64`: surface tension value
- `μ :: Float64`: kinematic viscosity, same as dynamic for ρ=1
- `θ :: Float64`: highest contact angle given as radiant, e.g. θ=π/9 for 20 degrees
e
# Mathematics



# Examples
```jldoctest
julia> using Swalbe, Statistics, Test

julia> x = ones(50,50); y = ones(50,50); h = ones(50,50);

julia> Swalbe.thermal!(x, y, h, 0.1, 1/6, 1)

julia> @test mean(x) ≈ 0.0 atol=1e-2
Test Passed

julia> @test mean(y) ≈ 0.0 atol=1e-2
Test Passed

julia> @test var(x) ≈ 2*0.1/11 atol=(2*0.1/11)/10 # var = 2kbt*6*μ/slip
Test Passed

```

# References

- [Mecke, Rauscher](https://iopscience.iop.org/article/10.1088/0953-8984/17/45/042/meta)
"""
function t0(;hᵦ=0.07, γ=0.01, μ=1/6, θ=1/6)
    qsq = hᵦ * (1 - cospi(θ)) * (2 - 3 * hᵦ) 
    charT = 3 * μ / (γ * qsq^2)

    return charT
end

"""
    snapshot!(snap, field, t; dumping)

Makes a copy of the current state of an input array `in` and saves the vectorized values as a column to `out`.

Function that fills a preallocated array out with a time series of system snap shots of e.g. the height field `h`.

# Arguments

- `snap :: Array{Number,2}`: Array that stores the snap shots as columns
- `field :: Array{Number,2}`: Input argument, e.g. ``h(\\mathbf{x},t)``
- `t :: Int`: The current time step
- `dumping :: Int`: Sampling frequency

# Examples

```jldoctest
julia> using Swalbe, Test

julia> h1 = reshape(collect(1:25),5,5); h2 = reshape(collect(5:5:125),5,5);

julia> snapshot = zeros(2, 25);

julia> Swalbe.snapshot!(snapshot,h1,10,dumping=10)

julia> Swalbe.snapshot!(snapshot,h2,20,dumping=10)

julia> @test all(h1 .== reshape(snapshot[1,:],5,5))
Test Passed

julia> @test all(h2 .== reshape(snapshot[2,:],5,5))
Test Passed

```
# References

See also: The [scripts folder](https://github.com/Zitzeronion/Swalbe.jl/tree/master/scripts) 
"""
function snapshot!(snap, field, t; dumping = 1000)
    if t % dumping == 0
        snap[t÷dumping, :] .= vec(Array(field))
    end
    
    return nothing
end

"""
    surfacearea!(area_lv, red_energy, height, θ, ∇hx, ∇hy, dgrad, surface,)

Measures the surface area of the liquid vapor interface and the reduced surface energy.

# Arguments

- `area_lv :: Vector{Float64}`: array to store the computed liquid vapor area
- `red_energy :: Vector{Float64}`: array the stores the computed reduced surface energy 
- `height :: Matrix{Float64}`: current height configuration
- `θ :: Matrix{Float64}`: contact angle field distribution
- `∇hx :: Matrix{Float64}`: height gradient x-component 
- `∇hy :: Matrix{Float64}`: height gradient y-component
- `dgrad :: Array{Float64,3}`: dummy array to store derivatives 
- `surface :: Matrix{Float64}`: array that computes locally the liquid vapor surface area 
- `t :: Int`: current time step
- `hthresh :: Float64`: height threshold below which the substrate is considered *dry*

"""
function surfacearea!(area_lv, red_energy, height, θ::Float64, ∇hx, ∇hy, dgrad, surface, t; htresh = 0.055)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surf = 0.0
    surface .= sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surf = sum(surface[height .> htresh])
    area_lv[t] = surf
    red_energy[t] = surf - length(height[height .> htresh]) * cospi(θ) 

    return nothing
end

function surfacearea!(area_lv, red_energy, height, θ::Matrix, ∇hx, ∇hy, dgrad, surface, t; htresh = 0.055)
    ∇f_simple!(∇hx, ∇hy, height, dgrad)
    surf = 0.0
    surface .= sqrt.(∇hx.^2 .+ ∇hy.^2 .+ 1)
    surf = sum(surface[height .> htresh])
    area_lv[t] = surf
    red_energy[t] = surf - sum(cospi.(θ[height .> htresh]))

    return nothing
end

"""
    ∇f_simple!(outputx, outputy, f, dgrad) 

Simple gradient calculation for the differential surface area.
"""
function ∇f_simple!(outputx, outputy, f, dgrad)
    fip, fjp, fim, fjm, fipjp, fimjp, fimjm, fipjm = Swalbe.viewneighbors(dgrad)
    # Straight elements  j+1, i+1, i-1, j-1
    circshift!(fip, f, (1,0))
    circshift!(fjp, f, (0,1))
    circshift!(fim, f, (-1,0))
    circshift!(fjm, f, (0,-1))
    # Diagonal elements  
    circshift!(fipjp, f, (1,1))
    circshift!(fimjp, f, (-1,1))
    circshift!(fimjm, f, (-1,-1))
    circshift!(fipjm, f, (1,-1))
    # In the end it is just a weighted sum...
    outputx .= -1/3 .* (fip .- fim) .- 1/12 .* (fipjp .- fimjp .- fimjm .+ fipjm)
    outputy .= -1/3 .* (fjp .- fjm) .- 1/12 .* (fipjp .+ fimjp .- fimjm .- fipjm)

    return nothing
end

"""
    center_of_mass(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D)

	function center_of_mass(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D)

Calculates the center of mass for a 1D height field for a droplet

The center of mass is computed as the weighted average position of the height field, where the weight is given by the height at each position.

# Arguments
- `state :: Swalbe.State_1D`: A state object containing the height field `state.height` representing the heights at each position along the 1D domain.
- `sys :: Swalbe.SysConst_1D`: A system constant object that provides the length `sys.L` of the domain.

# Returns
- `output :: Float64`: The position of the center of mass as a scalar value along the 1D domain.

# Example
```jldoctest
julia> center_of_mass(state, sys)
```
"""
function center_of_mass(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D)
    mass = sum(state.height)
    output = 0
    for i in 1:sys.L
        output += i*state.height[i]
    end
    output = output/mass
    return output
end

"""
    center_of_mass_corrected(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D; fac=1.1)

Calculate the corrected center of mass in 1D, excluding regions where the height is below a threshold.

# Arguments
- `state::Swalbe.State_1D`: The state of the system, containing the height field `state.height`.
- `sys::Swalbe.SysConst_1D`: The system constants, including `sys.hmin`, the minimum height for consideration.
- `fac::Real`: A factor used to adjust the threshold for the minimum height (`hmin*fac`). Default is 1.1.

# Returns
- A scalar value representing the center of mass, considering only the region with height above the threshold.
"""
function center_of_mass_corrected(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D; fac=1.1)
    a=findfirst(h->h>sys.hmin*fac, state.height)
    b=findlast(h->h>sys.hmin*fac, state.height)
    mass = sum(state.height[a:b])
    output = 0
    for i in a:b
        output += i*state.height[i]
    end
    output = output/mass
    return output
end


"""
    center_of_volume(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D; fac=1.1)

Calculate the center of volume in 1D, considering regions where the height exceeds a threshold.

# Arguments
- `state::Swalbe.State_1D`: The state of the system, containing the height field `state.height`.
- `sys::Swalbe.SysConst_1D`: The system constants, including `sys.hmin`.
- `fac::Real`: A factor used to adjust the threshold for the minimum height (`hmin*fac`). Default is 1.1.

# Returns
- A scalar value representing the center of volume in the 1D domain, within the region where the height is above the threshold.
"""
function center_of_volume(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D; fac=1.1)
        a=findfirst(h->h>sys.hmin * fac, state.height)
        b=findlast(h->h>sys.hmin * fac, state.height)
        return a+ (b-a)/2
end


"""
    center_of_mass(state::Swalbe.State, sys::Swalbe.SysConst)

Calculate the center of mass for a 2D system based on the height field.

# Arguments
- `state::Swalbe.State`: The state of the system, which includes the height field `state.height`.
- `sys::Swalbe.SysConst`: The system constants, including `sys.Lx` and `sys.Ly`, the dimensions of the system.

# Returns
- A 2D vector representing the center of mass in the x and y directions.
"""
function center_of_mass(state::Swalbe.State, sys::Swalbe.SysConst)
    mass = sum(state.height)
    output = [0.0, 0.0]
    for i in 1:sys.Lx
        for j in 1:sys.Ly
            output .= output .+ state.height[i,j] .* [i,j]
        end
    end
    output .= output./mass
    return output
end


"""
    center_of_volume_x(state::Swalbe.State, sys::Swalbe.SysConst; fac=1.1)

Calculate the center of volume in the x-direction for a 2D system, considering regions where the height exceeds a threshold.

# Arguments
- `state::Swalbe.State`: The state of the system, containing the height field `state.height`.
- `sys::Swalbe.SysConst`: The system constants, including `sys.Ly`.
- `fac::Real`: A factor used to adjust the threshold for the minimum height (`hmin*fac`). Default is 1.1.

# Returns
- A scalar value representing the center of volume along the x-direction in the region where the height is above the threshold.
"""
function center_of_volume_x(state::Swalbe.State, sys::Swalbe.SysConst; fac=1.1)
        a=findfirst(h->h>sys.hmin * fac, state.height[:,Int(floor(sys.Ly/2))])
        b=findlast(h->h>sys.hmin * fac, state.height[:,Int(floor(sys.Ly/2))])
        try
                return a+ (b-a)/2
        catch
                return NaN
        end
end

"""
    center_of_mass(height_in; fac=1.1)

Calculate the center of mass for a 2D height matrix by considering only regions above a height threshold.

# Arguments
- `height_in::Matrix{T}`: A 2D matrix representing the height field.
- `fac::Real`: A factor used to adjust the threshold for the minimum height (`min(height_in)*fac`). Default is 1.1.

# Returns
- A 2D vector representing the center of mass in the x and y directions.
"""
function center_of_mass(height_in; fac=1.1)
        Lx, Ly = size(height_in)
        height = zeros(Lx, Ly)
        height .= height_in
        prec = minimum(height)
        height[height .< prec * fac] .= 0.0
        mass = sum(height)
        output = [0.0, 0.0]
        for i in 1:Lx
                for j in 1:Ly
                        output .= output .+ height[i,j] .* [i,j]
                end
        end
        output .= output ./ mass
        return output
end


"""
    center_of_area(height_in; fac=1.1)

Calculate the center of area for a 2D height matrix, considering only regions above a height threshold.

# Arguments
- `height_in::Matrix{T}`: A 2D matrix representing the height field.
- `fac::Real`: A factor used to adjust the threshold for the minimum height (`min(height_in)*fac`). Default is 1.1.

# Returns
- A 2D vector representing the center of area in the x and y directions, where regions with height above the threshold are considered.
"""
function center_of_area(height_in; fac=1.1)
        Lx, Ly = size(height_in)
        height = zeros(Lx, Ly)
        height .= height_in
        prec = minimum(height)
        precursor = zeros(Lx, Ly)
        precursor[height .> fac * prec] .= 1.0
        mass = sum(precursor)
        output = [0.0, 0.0]
        for i in 1:Lx
                for j in 1:Ly
                        output .= output .+ precursor[i,j] .* [i,j]
                end
        end
        output .= output ./ mass
        return output
end

"""
    free_surface(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D)

Compute the free surface area in the 1D system by measuring the vertical changes in height across the domain.

# Arguments
- `state::Swalbe.State_1D`: The state of the system, containing the height field `state.height`.
- `sys::Swalbe.SysConst_1D`: The system constants, including `sys.hcrit`, the critical height for considering the free surface.

# Returns
- A scalar value representing the total vertical change in height across the free surface in the domain.
"""
function free_surface(state::Swalbe.State_1D, sys::Swalbe.SysConst_1D)
    output = 0
    for i in 22:sys.L
        if state.height[i] > sys.hcrit + 0.01
            output += abs(state.height[i] - state.height[i-1])
        end
    end
    return output
end


"""
	function count_connected_clusters(vector)

Counts the numbers of concected clusters in a 1D array containing 0 and 1 as values
"""
function count_connected_clusters(vector)
# Find the indices of ones in the vector
one_indices = findall(x -> x == 1, vector)
count=0
for i in 1:length(vector)
	if vector[Swalbe.modified_modulo(i+1,length(vector))]!=vector[i]
		count+=1
	end
end
	count /= 2
	return count
end

"""
	function bridge_height(state::Swalbe.StateMiscible_1D, sys::Swalbe.SysConstMiscible_1D;corrector=20,verbose=true)

Measures the bridge height of two coalescing droplets. 
"""
function bridge_height(state::Swalbe.StateMiscible_1D, sys::Swalbe.SysConstMiscible_1D;corrector=20,verbose=true, threshb=0.001, shift=6, prec_fac=2.0)
	a=findmax(state.height[:,1])[2]
	# b=findmax(state.height[:,2])[2]
    	hip, him = viewneighbors_1D(state.dgrad[:,:,1])
    	# Straight elements j+1, i+1, i-1, j-1
	circshift!(hip, state.height[:,2], shift)
	circshift!(him, state.height[:,2], -shift)
	b=findlast(x->x, (hip .+ threshb .< state.height[:,2]) .&& (him .+ threshb .< state.height[:,2]) .&& state.height[:,2].> prec_fac*minimum(state.height[:,2]) )[1]
	b=max(a,b)
	mm=findmin((state.height[:,1].+state.height[:,2])[a:b])
	nm=mm[2]+a-1==b ? findmax(state.height[:,1] .+ state.height[:,2]) : [mm[1],mm[2]+a-1]
	if verbose
		println("a=$a,b=$b,bridge_height=$(mm[1]),bridge_position=$(nm[2])")
	end
	return a, b, nm
end

"""
    contact_line(state::LBM_state_1D)

Will find the first three phase contact line coming from the left, returns 1 if none found
"""
function contact_line(state::LBM_state_1D)
    start=1
	for i in start:length(state.height)
		if state.height[i]<0.051
			start=i
			break
		end
	end
    for x in start:length(state.height)	
        if state.height[x]>0.052
            return x
            break
        end
    end
    return 1
end

"""
	function countdroplets(state, sys; thresh_extra=0.05)

Counts the number of droplets in a 1D height field array
"""
function countdroplets(state, sys; thresh_extra=0.05)
	thresh = sys.hmin -sys.hcrit .+ thresh_extra
	y=zeros(sys.L)
	y[state.height.<=thresh].= 0
	y[state.height.>thresh].= 1
	return count_connected_clusters(y)
end

"""
    function arc_contact_angle_theta(str; theta=pi/9,th=0.06, j=1)

helper function for `contact_angle`, Fits a arc of angle theta onto the data of `height[:,j]`
"""
function arc_contact_angle_theta(height; theta=pi/9,th=0.06, j=1, exponent=2)
	cl=findfirst(h-> h>th, height[:,j])[1]
	cr=findlast(h-> h>th, height[:,j])[1]
	h_m, x_m=findmax(height[:,j])
	hmin=minimum(height[:,j])
	r=h_m/(1-cos(theta))
	extra_h=j==2 ? height[x_m,1] : 0
	arc=hmin*ones(size(height[:,j]))
	for x in cl:cr
		arc[x]=max(sqrt(abs.(r^2-(x-x_m)^2))-r+h_m + extra_h, hmin)
	end
	l2norm=sum(abs.(height[cl:cr,j] .- arc[cl:cr]) .^exponent)
	return l2norm
end


"""
    function contact_angle(state, sys;  j=1,  step=0.01,min_angle=0.9*pi/18, max_angle=2.1*pi/9, extra_prec=0.1)

meassures the contact angle of a droplet by fitting an arc to it. This is done by iterating over `min_angle:step:max_angle` return the one with the smallest difference in l_2 norm. `extra_prec` defines how much beyond precursor heigth to define the contact line
"""
function contact_angle(state, sys;  j=1,  step=0.0001,min_angle=pi/64, max_angle=pi/2, extra_prec=0.1, exponent=2)
	dev=[0.0 for theta in min_angle:step:max_angle]
	for (i,theta) in enumerate(min_angle:step:max_angle)
		dev[i]=arc_contact_angle_theta(state.height, theta=theta, th=sys.hmin - sys.hcrit + extra_prec, j=j)
	end
	return min_angle+(findmin(dev)[2]-1)*step
end

"""
    function contact_angle_error(state, sys, theta_0;  j=1,  step=0.01,min_angle=0.9*pi/18, max_angle=2.1*pi/9, extra_prec=0.1)

meassures the contact angle of a droplet by fitting an arc to it. This is done by iterating over `min_angle:step:max_angle` return the one with the smallest difference in l_2 norm. `extra_prec` defines how much beyond precursor heigth to define the contact line. Return the arrow with respect to an expected contact angle
"""
function contact_angle_error(state, sys, theta_0;  j=1,   extra_prec=0.1, exponent=1)
	th=sys.hmin - sys.hcrit + extra_prec
	cl=findfirst(h-> h>th, state.height[:,j])[1]
	cr=findlast(h-> h>th, state.height[:,j])[1]
	return arc_contact_angle_theta(state.height, theta=theta_0, th=th, j=j, exponent=exponent) / sum(state.height[cl:cr,j].^exponent)
end


"""
    contact_line_right(state::LBM_state_1D)

Will find the first three phase contact line coming from the right, returns length(state.height) if none found
"""
function contact_line_right(state::LBM_state_1D)
    start=length(state.height)
	for i in length(state.height):-1:1
		if state.height[i]<0.051
			start=i
			break
		end
	end
    for x in start:-1:1
        if state.height[x]>0.052
            return x
            break
        end
    end
    return length(state.height)
end


"""
    contact_angle(state::LBM_state_1D; n=1)

Will calculate the apparent contact angle at the first contact line from the left, is just calculating an angle at x=1 if nothing is dewetted 

# Arguments

- `n:: Int`: size of the triangle used to calculate the contact angle
"""
function contact_angle(state::LBM_state_1D; n=5)
    x=Swalbe.contact_line(state)
    return asin((state.height[x+n]-state.height[x])/n)
end


"""
    drop_height(state)

Return maximum(state.height)
"""
function drop_height(state::LBM_state_1D)
    return maximum(state.height)
end


"""
    drop_width(state)

Gives width of a droplet, if there are multiple droplets it will return length of the leftest droplet to the rightest 
"""
function drop_width(state::LBM_state_1D)
    return -Swalbe.contact_line(state)+Swalbe.contact_line_right(state)
end

"""
    mymean(m:: Vector)

returns location of the mean of the Input vector interpreted as a probability density
"""
function mymean(m:: Vector)
	μ=0
	for i in 1:length(m)
		μ+= m[i]*i
	end
	return μ/sum(m)
end


"""
    mystd(m:: Vector)

returns standard deviation of the Input vector interpreted as a probability density
"""
function mystd(m:: Vector)
	σ=0
	μ=mymean(m)
	for i in 1:length(m)
		σ += m[i]*((i-μ)^2)
	end
	return sqrt(σ/length(m))
end

"""
    mykurtosis(m:: Vector)

returns kurtosis of the Input vector interpreted as a probability density
"""
function mykurtosis(m:: Vector)
	σ=mystd(m)
	μ=mymean(m)
	k=0
	for i in 1:length(m)
		k += m[i]*(((i-μ)/σ)^4)
	end
	return k/length(m)-3
end

"""
	function bridge_height_surf(state::Swalbe.StateActive_1D, sys::Swalbe.SysConstActive_1D;corrector=20,verbose=true)

Measures the bridge height of two coalescing droplets. 
"""
function bridge_height(state::Swalbe.StateActive_1D, sys::Swalbe.SysConstActive_1D;corrector=20,verbose=true, threshb=0.001, shift=6)
    	hip, him = viewneighbors_1D(state.dgrad[:,:,1])
    	# Straight elements j+1, i+1, i-1, j-1
	circshift!(hip, state.height, shift)
	circshift!(him, state.height, -shift)
	a=findfirst(x->x, (hip .+ threshb .< state.height) .&& (him .+ threshb .< state.height))[1]
	b=findlast(x->x, (hip .+ threshb .< state.height) .&& (him .+ threshb .< state.height))[1]
	mm=findmin((state.height)[a:b])
	nm=mm[2]+a-1==b ? findmax(state.height) : [mm[1],mm[2]+a-1]
	if verbose
		println("a=$a,b=$b,bridge_height=$(mm[1]),bridge_position=$(nm[2])")
	end
	return a, b, nm
end

"""
    function catalytic_yield(state::Swalbe.StateActive_1D, sys::Swalbe.SysConstActive_1D)

meassures the amount of product that is given to the atmosphere at one time step, thus it has units of `[rho_A]m^(d-1)/s`, i.e. either `1/s` or `kg/s` depending on whether you consider `rho_A` a number or mass density
""" 
function catalytic_yield(state::Swalbe.StateActive_1D, sys::Swalbe.SysConstActive_1D)
    @. state.precursor = ifelse(state.height <= sys.hmin - sys.hcrit + 0.01, 0,1)
    yield = sum((sys.sigma_A_up .* state.rho_A )./ (state.height ) .* state.precursor)
    return yield
end
    
"""
    function arc_neuman_angle_theta_up(height; theta=pi/9,th=0.06)

Creates a upper arc between the two contact lines of a given simulation snapshot with Neuman angle `theta`, compare to `arc_contact_angle_theta`
"""
function arc_neuman_angle_theta_up(height; theta=pi/9,th=0.06, exponent=2)
	h_low=maximum(height[:,1])
	cl=findfirst(h-> h>th, height[:,2])[1]
	cr=findlast(h-> h>th, height[:,2])[1]
	h_m, x_m=findmax(height[:,2] .+ height[:,1] .- h_low)
	r=h_m/(1-cos(theta))
	extra_h=h_low
	arc=h_low*ones(size(height[:,2]))
	for x in cl:cr
		arc[x]=max(sqrt(abs.(r^2-(x-x_m)^2))-r+h_m + extra_h, h_low)
	end
	l2norm=sum(abs.(height[cl:cr,1] .+ height[cl:cr,2] .- arc[cl:cr]) .^exponent)
	return l2norm
end
"""
    function arc_neuman_angle_theta_low(height; theta=pi/9,th=0.06)

Creates a lower arc between the two contact lines of a given simulation snapshot with Neuman angle `theta`, compare to `arc_contact_angle_theta`
"""
function arc_neuman_angle_theta_low(height; theta=pi/9,th=0.06, exponent=2)
	h_low=maximum(height[:,1])
	cl=findfirst(h-> h>th, height[:,2])[1]
	cr=findlast(h-> h>th, height[:,2])[1]
	h_m, x_m=findmax(-height[:,1] .+ h_low)
	r=h_m/(1-cos(theta))
	extra_h=h_low
	arc=h_low*ones(size(height[:,2]))
	for x in cl:cr
		arc[x]=min(r-sqrt(abs.(r^2-(x-x_m)^2))-h_m + extra_h, h_low)
	end
	l2norm=sum(abs.(height[cl:cr,1]  .- arc[cl:cr]) .^exponent)
	return l2norm
end
"""
    function neuman_angle(state,sys; step=0.01,min_angle=0.9*pi/18, max_angle=2.1*pi/9, extra_prec=0.1)

Measures the Neuman angle of a simulation snapshot by fitting an arc with the interface, compare to `contact_angle`
"""
function neuman_angle(state,sys; step=0.0001,min_angle=pi/128, max_angle=pi/2, extra_prec=0.1, exponent=2)
	dev_up=[0.0 for theta in min_angle:step:max_angle]
    th=sys.hmin - sys.hcrit + extra_prec
	for (i,theta) in enumerate(min_angle:step:max_angle)
		dev_up[i]=arc_neuman_angle_theta_up(state.height, theta=theta, th=th, exponent=exponent)
	end
	dev_low=[0.0 for theta in min_angle:step:max_angle]
	for (i,theta) in enumerate(min_angle:step:max_angle)
		dev_low[i]=arc_neuman_angle_theta_low(state.height, theta=theta, th=th, exponent=exponent)
	end
	return min_angle+(findmin(dev_up)[2]-1)*step, min_angle+(findmin(dev_low)[2]-1)*step
end

"""
    function neuman_angle_error(state,sys, theta_0; step=0.01,min_angle=0.9*pi/18, max_angle=2.1*pi/9, extra_prec=0.1)

Measures the Neuman angle of a simulation snapshot by fitting an arc with the interface, compare to `contact_angle`
"""
function neuman_angle_error(state,sys,theta_01, theta_02;  extra_prec=0.1, exponent=1)
    th=sys.hmin - sys.hcrit + extra_prec
	cl=findfirst(h-> h>th, state.height[:,2])[1]
	cr=findlast(h-> h>th, state.height[:,2])[1]
	e1=arc_neuman_angle_theta_up(state.height, theta=theta_01, th=th, exponent=exponent)
	e2=arc_neuman_angle_theta_low(state.height, theta=theta_02, th=th, exponent=exponent)
	return (e1+e2)/sum(state.height[cl:cr,2].^exponent)
end
