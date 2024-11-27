using Swalbe

"""
    function_rho_LBM!(state,sys)
does a LBM time step for ``\\rho``

## Numerics
Evolve on the same lattice as ``h`` resp ``f``
Let 
```
\\begin{align}
\\sum g_i &= \\rho\\
V &= \\rho v
\\end{align}
```
do
```
g_i(x+ \\Delta t c_i, t+ \\Delta t) = g_i(x,t) + \\frac{\\Delta t}{\\tau}(g_i^{eq}-f_i) + \\Delta t \\frac{w_i}{c_s^2} c_i \\cdot F_{tot}
```
where 
```
\\begin{align}
     g_i^{eq} = \\begin{cases}
          \\rho  - \\frac{2\\rho v^2}{3c_s^2}& i = 0 \\
 \\frac{\\rho c_i\\cdot \\mathbf{v}}{3 c_s^2} + \\frac{\\rho (c_i\\cdot \\mathbf{v})^2}{2 c_s^4}-\\frac{\\rho v^2}{6c_s^2} &  i=1,3,5,7 \\
 \\frac{\\rho c_i\\cdot\\mathbf{v}}{12 c_s^2} + \\frac{\\rho (c_i\\cdot\\mathbf{v})^2}{8c_s^4}-\\frac{\\rho v^2}{24c_s^2} & i=2,4,6,8
     \\end{cases}
 \\end{align}
 ```
and 
```
V = -\\frac{D\\nabla\\rho(x,t)}{\\rho}+M\\Gamma\\rho(x,t)\\nabla\\rho(x,t)-
```
then ``\\rho ,v`` solves by assumption ``\\frac{\\rho*}{\\lambda*}<<1`` and ``Re<<1`` at leading order
 ```
 \\begin{align}
 \\begin{split}
     \\partial_t \\rho + \\text{div} (\\rho v) &= 0\\
     0 \\overset{\\text{@equil}}{=} \\partial_t (\\rho v) &=-u-  F(\\rho)
     \\end{split}
 \\end{align}
 ```
so finally
```
\\partial_t \\rho - \\text{div}(D\\nabla\\rho(x,t)+M\\Gamma\\rho(x,t)\\nabla\\rho(x,t)-D \\rho(x,t)\\nabla \\ln(h(x,t)))=0
```
as desired. 

In termes of code we therefore do the following algorithm
```
for t in 1:Tmax
    do LBM for h
    forcing!(F,rho,h)
    equilibrium!(g^eq,g)
    BGKandStream!(g,g^eq,F)
    moments!((rho,u),g)
end 
```

Not as easy to implement as an explicit method but at the end of the day it is coded straight forward and may work. 
At least LBM is known to be super well at mass conservation with is indeed our problem. 
Also that way we don't have to much trouble coupling the time steps. 

"""
function update_rho_LBM!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_pressure!(state, sys)
    Swalbe.rho_grad_disj_p!(state)
    state.rho_F.= state.rho∇p
    Swalbe.rho_equilibrium_quadratic!(state)
    Swalbe.rho_BGKandStream!(state,sys)
    if sys.alpha != 0
        Swalbe.nabla_F!(state,sys)
        Swalbe.V_p!(state,sys)
        Swalbe.V_gamma!(state,sys)
    end
    Swalbe.rho_moments!(state,sys)
end


"""
	function update_rho_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)

LBM routine for the update of the chemical compounds in the active/catalysis model of [Richter et all.](https://arxiv.org/abs/2402.14635), see the publication and ´update_rho_LBM!´. 

The difference here as compared to ´update_rho_LBM!´ is that the precursor film is handeled seperately by a if - else condition to avoid unwanted effects in the precursor layer. This difference only manifests in the forcing. See ´nabla_F_precursor!´. 
"""
function update_rho_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_pressure!(state, sys)
    Swalbe.rho_grad_disj_p!(state)
    state.rho_F.= state.rho∇p
    Swalbe.rho_equilibrium_quadratic!(state)
    Swalbe.rho_BGKandStream!(state,sys)
    if sys.alpha != 0
        Swalbe.nabla_F_precursor!(state,sys)
        Swalbe.V_p!(state,sys)
        Swalbe.V_gamma!(state,sys)
    end
    Swalbe.rho_moments!(state,sys)
end


"""
	function update_rho_A_LBM!(state::Active_1D, sys:: SysConstActive_1D)


updates the rho_A acording to the LBM algorithm described in `update_rho_LBM!` see also [Richter et all.](https://arxiv.org/abs/2402.14635)
"""
function update_rho_A_LBM!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_A_pressure!(state, sys)
    Swalbe.rho_A_grad_disj_p!(state)
    state.rho_A_F.= state.rho_A∇p
    Swalbe.rho_A_equilibrium_quadratic!(state)
    Swalbe.rho_A_BGKandStream!(state,sys)
    Swalbe.rho_A_moments!(state,sys)
    Swalbe.react_rho_A!(state,sys)
end

"""
	function update_rho_B_LBM!(state::Active_1D, sys:: SysConstActive_1D)


updates the rho_B acording to the LBM algorithm described in `update_rho_LBM!` see also [Richter et all.](https://arxiv.org/abs/2402.14635)
"""

function update_rho_B_LBM!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_B_pressure!(state, sys)
    Swalbe.rho_B_grad_disj_p!(state)
    state.rho_B_F.= state.rho_B∇p
    Swalbe.rho_B_equilibrium_quadratic!(state)
    Swalbe.rho_B_BGKandStream!(state,sys)
    Swalbe.rho_B_moments!(state,sys)
    Swalbe.react_rho_B!(state,sys)
end

"""
	function update_rho_A_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)


updates the rho_A acording to the LBM algorithm described in `update_rho_LBM!` but handels the reaction in the precursor differently, see `react_rho_A_precurso`, see also [Richter et all.](https://arxiv.org/abs/2402.14635)
"""
function update_rho_A_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_A_pressure!(state, sys)
    Swalbe.rho_A_grad_disj_p!(state)
    state.rho_A_F.= state.rho_A∇p
    Swalbe.rho_A_equilibrium_quadratic!(state)
    Swalbe.rho_A_BGKandStream!(state,sys)
    Swalbe.rho_A_moments!(state,sys)
    Swalbe.react_rho_A_precursor!(state,sys)
end



"""
	function update_rho_A_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)


updates the rho_A acording to the LBM algorithm described in `update_rho_LBM!` but handels the reaction in the precursor differently, see `react_rho_A_precurso`, see also [Richter et all.](https://arxiv.org/abs/2402.14635)
"""

function update_rho_B_LBM_precursor!(state::Active_1D, sys:: SysConstActive_1D)
    Swalbe.rho_B_pressure!(state, sys)
    Swalbe.rho_B_grad_disj_p!(state)
    state.rho_B_F.= state.rho_B∇p
    Swalbe.rho_B_equilibrium_quadratic!(state)
    Swalbe.rho_B_BGKandStream!(state,sys)
    Swalbe.rho_B_moments!(state,sys)
    Swalbe.react_rho_B_precursor!(state,sys)
end
# Surfactant
"""
	function surface_tension_surfactant!(state::Active_1D, sys::SysConstActive_1D)

calculates the surface tension of a thin film with a surfactant where `state.height` is the film heigt, and `state.rho` is the surfactant concentration according to
```
\\gamma=\\gamma_0 - \\Gamma \\rho
```
"""
function surface_tension_surfactant!(state::Active_1D, sys::SysConstActive_1D)
    state.γ .= sys.γ_0 .- sys.Γ .* state.rho
end

# Solutes
"""
    function surface_tension_solute!(state::Active_1D, sys::SysConstActive_1D)

Calculates the surface tension of a thin film with a dissolved solute. Here, `state.height` is the film height, and `state.rho` is the solute concentration. The surface tension is calculated as

γ = γ_0 - Γ * ρ / h

where `γ_0` is the base surface tension and `Γ` is the solute parameter.
"""
function surface_tension_solute!(state::Active_1D, sys::SysConstActive_1D)
    state.γ .= sys.γ_0 .- sys.Γ .* state.rho ./ state.height
end

"""
    function surface_tension_fourier!(state::Active_1D, sys::SysConstActive_1D)

Calculates the surface tension of a thin film in Fourier space, allowing for more efficient convolution with a kernel. The surface tension is given by

γ = γ_0 - Γ * irfft(rfft(ρ) * fourier_kernel)

where `fourier_kernel` is a filter applied in Fourier space, and `irfft` is the inverse Fourier transform.
"""
function surface_tension_fourier!(state::Active_1D, sys::SysConstActive_1D)
    state.γ .= sys.γ_0 .- sys.Γ .*  irfft( rfft(state.rho) .* state.fourier_kernel, sys.L)
end


# catalysis as done in the catalysis in thin liuqid films paper
"""
    function surface_tension_4_fields!(state::Active_1D, sys::SysConstActive_1D)

Calculates the surface tension of a thin film with two additional concentration fields, the product `state.rho_A` and the reactant `state.rho_B`, each contributing to surface tension based on parameters `Γ` and `GammaB` as

γ = γ_0 - (Γ * ρ_A + GammaB * ρ_B) / height

This function models surface tension considering dual surfactants or solutes, with each term affecting the overall surface tension.

see [Richter et all.](https://arxiv.org/abs/2402.14635)

"""
function surface_tension_4_fields!(state::Active_1D, sys::SysConstActive_1D)
    state.γ .= sys.γ_0 .- (sys.Γ .* state.rho_A  .+ sys.GammaB .* state.rho_B) ./state.height
end





"""
	function window_integration!(state::Active_1D, sys::SysConstActive_1D)

smoothes out the concentration fields by a moving averadge, in some cases this helps to control simulations that suffer from low frequency numerical instabilities. 
"""
function window_integration!(state::Active_1D, sys::SysConstActive_1D)
    rip, rim = viewneighbors_1D(state.dgrad)
    circshift!(rip, state.rho, 1)
    circshift!(rim, state.rho, -1)
    state.rho .= state.rho  .* sys.weight_center .+ (1-sys.weight_center)/2 .* rip .+ (1-sys.weight_center)/2 .* rim 
    circshift!(rip, state.rho_A, 1)
    circshift!(rim, state.rho_A, -1)
    state.rho_A .= state.rho_A  .* sys.weight_center .+ (1-sys.weight_center)/2 .* rip .+ (1-sys.weight_center)/2 .* rim 
    circshift!(rip, state.rho_B, 1)
    circshift!(rim, state.rho_B, -1)
    state.rho_B .= state.rho_B  .* sys.weight_center .+ (1-sys.weight_center)/2 .* rip .+ (1-sys.weight_center)/2 .* rim 
    # println(sys.weight_center)
end


"""

    function react_rho_A!(state::Active_1D, sys::SysConstActive_1D)

performs reaction as 
`\\rho_A(t+1)=\\rho_A(t)+\\omega \\rho(t) - \\sigma_1 \\rho_A(t)-\\sigma_2\\rho_A(t)`` 
"""
function react_rho_A!(state::Active_1D, sys::SysConstActive_1D)
	state.rho_A .+= (sys.production_rate .* state.rho .* state.rho_B - sys.sigma_A_up .* state.rho_A )./ (state.height )
end


# Corrected for a factor 1/h in the production term of the chemical reaction
"""
    function react_rho_A_precursor!(state::Active_1D, sys::SysConstActive_1D)

Performs a reaction update for `rho_A` that includes a precursor factor dependent on film height, such that the reaction is handeled differently in the precursor and in the bulk film. Remember that the precursor is to be considered as a dewetted part of the film. The reaction is scaled by the precursor factor:

ρ_A(t+1) = ρ_A(t) + (ω * ρ * ρ_B - σ₁ * ρ_A) / height * precursor

where `precursor` is a binary field (1 or 0) determined by a height threshold condition, scaling the reaction terms based on a minimal height criterion.

Compare `react_rho_A!
"""
function react_rho_A_precursor!(state::Active_1D, sys::SysConstActive_1D)
    @. state.precursor = ifelse(state.height <= sys.hmin - sys.hcrit + 0.01, 0,1)
    state.rho_A .+= (sys.production_rate .* state.rho .* state.rho_B - sys.sigma_A_up .* state.rho_A )./ (state.height ) .* state.precursor
end

"""
    function react_rho_B_precursor!(state::Active_1D, sys::SysConstActive_1D)

Updates `rho_B` by applying a reaction with a precursor factor dependent on film height. The update follows:

ρ_B(t+1) = ρ_B(t) + [(-ω * ρ * ρ_B - σ₂ * ρ_B) / height + ρ_BRes * σ₂_down] * precursor

where `ω` is the reaction rate, `σ₂` is a decay or removal term, `height` is a scaling factor, and `precursor` is determined based on a height criterion.

Compare `react_rho_A!`, `react_rho_B!`, and `react_rho_A_precursor!`
"""
function react_rho_B_precursor!(state::Active_1D, sys::SysConstActive_1D)
    @. state.precursor = ifelse(state.height <= sys.hmin - sys.hcrit + 0.01, 0,1)
    state.rho_B .+= ((- sys.production_rate .* state.rho .* state.rho_B - sys.sigma_B_up .* state.rho_B )./ (state.height ) .+ sys.rho_BRes_sigma_B_down) .* state.precursor
end

function react_rho_B!(state::Active_1D, sys::SysConstActive_1D)
	state.rho_B .+= (- sys.production_rate .* state.rho .* state.rho_B - sys.sigma_B_up .* state.rho_B )./ (state.height ) .+ sys.rho_BRes_sigma_B_down
end







# linear stability analysis coefficients
"""
	function f(h_0, sys::SysConstActive_1D)

Computes the disjoining pressure for analysis purposes. Not used in simulation but only in in-script postprocesing. 

``f(h)=\\frac{(n-1)(m-1)}{(n-m)h^*}\\left[\\left(\\frac{h^*}{h_0}\\right)^n-\\frac{h^*}{h_0}\\right)^m\\right]``
"""
function f(h_0, sys::SysConstActive_1D)
    fh = (sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)*((sys.hmin/h_0)^sys.n-(sys.hmin/h_0)^sys.m)
    return fh
end


"""
	function fprime(h_0, sys::SysConstActive_1D)

Computes the disjoining pressure for analysis purposes. Not used in simulation but only in in-script postprocesing. 

``f(h)=\\frac{(n-1)(m-1)}{(n-m)h^{*}}\\left[\\left(-n\\frac{h^*}{h_0}\\right)^n+m\\frac{h^*}{h_0}\\right)^{m}\\right]\\frac{1}{h_0}``
"""
function fprime(h_0, sys::SysConstActive_1D)
	fp = (sys.n-1)*(sys.m-1)/((sys.n-sys.m)*sys.hmin)*(-sys.n*(sys.hmin/h_0)^sys.n+sys.m*(sys.hmin/h_0)^sys.m)/h_0
	return fp
end

"""
	function Mob(h_0, sys::SysConstActive_1D)

The mobility
``\\frac{2h_0^2+6\\delta h_0 + 3 \\delta^2}{6}``
Not used in simulation but only in in-script postprocesing. 
"""
function Mob(h_0, sys::SysConstActive_1D)
	return(2*h_0^2+6*sys.δ* h_0+ 3*sys.δ^2 )/(6) 
end


"""
    function A_surfactant(q, sys::SysConstActive_1D; rho_0=1.0, h_0=1.0)

Computes the LSA matrix `A` for a system with a surfactant layer. This matrix is part of the stability analysis.

The function calculates elements based on system constants, including mobility `Mob(h_0)`, surface tension, and disjoining pressure.

Not used in simulation but only in in-script postprocessing.
"""
function A_surfactant(q, sys::SysConstActive_1D; rho_0=1.0, h_0=1.0)
    A=zeros(2,2)
    M=Mob(h_0,sys)
    Sigma=q^2*(sys.γ_0-sys.Γ*rho_0/h_0)- (sys.γ_0 - sys.Γ* rho_0/h_0 - sys.γ_0 * cospi(sys.θ_0))*fprime(h_0,sys)
    foo=f(h_0,sys)
    N=(sys.δ^2 + 2 *h_0*sys.δ + h_0^2)/2
    A[1,1]=-q^2 * M * h_0   / sys.μ * Sigma
    A[1,2]=-q^2 * N * rho_0 / sys.μ * Sigma
    A[2,1]=             - q^2 * M * h_0   / sys.μ * sys.Γ * foo - q^2 * h_0   * (h_0 + 2*sys.δ) / (2*sys.μ) * sys.Γ
    A[2,2]=-q^2 * sys.D - q^2 * N * rho_0 / sys.μ * sys.Γ * foo - q^2 * rho_0 * (sys.δ + h_0 )  / (sys.μ)   * sys.Γ 
    return A
end

"""
    function A_solute(q, sys::SysConstActive_1D; rho_0=1.0, h_0=1.0)

Computes the LSA matrix `A` for a system with solutes. This matrix is used for stability analysis with a solute, including mobility, surface tension, and disjoining pressure.

Not used in simulation but only in in-script postprocessing.
"""
function A_solute(q, sys::SysConstActive_1D; rho_0=1.0, h_0=1.0)
    A=zeros(2,2)
    M=Mob(h_0,sys)
    Sigma=q^2*(sys.γ_0-sys.Γ*rho_0)- (sys.γ_0 - sys.Γ* rho_0 - sys.γ_0 * cospi(sys.θ_0))*fprime(h_0,sys)
    foo=f(h_0,sys)
    N=(h_0 + 2* sys.δ)/2
    A[1,1]=-q^2 * M * h_0   / sys.μ * Sigma + q^2 * M * h_0     / sys.μ * sys.Γ * rho_0 / h_0^2 * foo + q^2 * h_0    * N / sys.μ * sys.Γ * rho_0 / h_0^2
    A[1,2]=-q^2 * M * rho_0 / sys.μ * Sigma + q^2 * M * rho_0   / sys.μ * sys.Γ * rho_0 / h_0^2 * foo + q^2 * rho_0  * N / sys.μ * sys.Γ * rho_0 / h_0^2 + q^2 * sys.D * rho_0 / h_0 
    A[2,1]=             - q^2 * M * h_0   / sys.μ * sys.Γ / h_0 * foo - q^2 * h_0   * N / sys.μ * sys.Γ / h_0
    A[2,2]=-q^2 * sys.D - q^2 * M * rho_0 / sys.μ * sys.Γ / h_0 * foo - q^2 * rho_0 * N / sys.μ * sys.Γ / h_0
    return A
end



"""
    function Sigma(q, sys::SysConstActive_1D; h_0=1.0, rho_0A=10, rho_0B=10)

Computes the dispersion relation `Sigma` for the stability analysis of a system with multiple fields, incorporating surface tensions and densities.

Not used in simulation but only in in-script postprocessing.
"""
function Sigma(q, sys::SysConstActive_1D;  h_0=1.0, rho_0A=10, rho_0B=10)
    return q^2 * ( sys.γ_0 - cospi(sys.θ_0 )* sys.γ_0 -( sys.Γ* rho_0A .+ sys.GammaB * rho_0B) / h_0)*fprime(h_0,sys)- q^4* ( sys.γ_0 - (sys.Γ * rho_0A + sys.GammaB * rho_0B) / h_0)
end


# Corrected for a factor 1/h in the production term of the chemical reaction
"""
	function A(q, sys::SysConstActive_1D; rho_0=1.0, h_0=1.0)

Calculates the LSA matrix. Formula not reportet here. See [Richter et all.](https://arxiv.org/abs/2402.14635)
Not used in simulation but only in in-script postprocesing. 
"""

function A(q, sys::SysConstActive_1D; h_0=1.0, rho_0=1.0, rho_0A=10, rho_0B=10)
    Sig = Sigma(q,sys,h_0=h_0, rho_0A=rho_0A, rho_0B=rho_0B) 
    V_gamma_h=(h_0 + 2* sys.δ)/(2)
    A=zeros(4,4)
    A[1,1]=                                              q^2 * V_gamma_h          * h_0    /sys.μ *( sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 + q^2 * Mob(h_0,sys)*     h_0    /sys.μ * (sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 * f(h_0,sys) + Mob(h_0,sys)*    h_0    /sys.μ * Sig  
    A[2,1]= -q^2 * sys.D         * rho_0  * ∇F(h_0,sys)+ q^2 * V_gamma(h_0,sys)   * rho_0  /sys.μ *( sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 + q^2 * V_p(h_0,sys)*     rho_0  /sys.μ * (sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 * f(h_0,sys) + V_p(h_0,sys)*    rho_0  /sys.μ * Sig
    A[3,1]= -q^2 * sys.D_Product * rho_0A * ∇F_0(h_0)  + q^2 * V_gamma_0(h_0,sys) * rho_0A /sys.μ *( sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 + q^2 * V_p_0(h_0, sys)*  rho_0A /sys.μ * (sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 * f(h_0,sys) + V_p_0(h_0, sys)* rho_0A /sys.μ * Sig +sys.sigma_A_up*rho_0A/h_0^2 - sys.production_rate*rho_0B*rho_0/h_0^2
    A[4,1]= -q^2 * sys.D_B       * rho_0B * ∇F_0(h_0)  + q^2 * V_gamma_0(h_0,sys) * rho_0B /sys.μ *( sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 + q^2 * V_p_0(h_0, sys)*  rho_0B /sys.μ * (sys.Γ * rho_0A +sys.GammaB*rho_0B)/h_0^2 * f(h_0,sys) + V_p_0(h_0, sys)* rho_0B /sys.μ * Sig +sys.sigma_B_up*rho_0B/h_0^2 + sys.production_rate*rho_0B*rho_0/h_0^2
    A[1,2]= 0
    A[2,2]= -sys.D*q^2
    A[3,2]=  sys.production_rate*rho_0B/h_0
    A[4,2]= -sys.production_rate*rho_0B/h_0
    A[1,3]=                      -q^2 * V_gamma_h          *h_0    / sys.μ   *sys.Γ/h_0      - q^2 * Mob(h_0,sys)    * h_0    /sys.μ* sys.Γ      /h_0 * f(h_0,sys)
    A[2,3]=                      -q^2 * V_gamma(h_0,sys)   *rho_0  / sys.μ   *sys.Γ/h_0      - q^2 * V_p(h_0,sys )   * rho_0  /sys.μ* sys.Γ      /h_0 * f(h_0,sys)
    A[3,3]= -q^2 * sys.D_Product -q^2 * V_gamma_0(h_0, sys)*rho_0A / sys.μ   *sys.Γ/h_0      - q^2 * V_p_0(h_0, sys) * rho_0A /sys.μ* sys.Γ      /h_0 * f(h_0,sys) -sys.sigma_A_up / h_0
    A[4,3]=                      -q^2 * V_gamma_0(h_0, sys)*rho_0B / sys.μ   *sys.Γ/h_0      - q^2 * V_p_0(h_0, sys) * rho_0B /sys.μ* sys.Γ      /h_0 * f(h_0,sys) 
    A[1,4]=                      -q^2 * V_gamma_h          *h_0    / sys.μ   *sys.GammaB/h_0 - q^2 * Mob(h_0,sys)    * h_0    /sys.μ* sys.GammaB /h_0 * f(h_0,sys)
    A[2,4]=                      -q^2 * V_gamma(h_0,sys)   *rho_0  / sys.μ   *sys.GammaB/h_0 - q^2 * V_p(h_0,sys )   * rho_0  /sys.μ* sys.GammaB /h_0 * f(h_0,sys)
    A[3,4]=                      -q^2 * V_gamma_0(h_0, sys)*rho_0A / sys.μ   *sys.GammaB/h_0 - q^2 * V_p_0(h_0, sys) * rho_0A /sys.μ* sys.GammaB /h_0 * f(h_0,sys) + sys.production_rate*rho_0/h_0                            
    A[4,4]= -q^2 * sys.D_Product -q^2 * V_gamma_0(h_0, sys)*rho_0B / sys.μ   *sys.GammaB/h_0 - q^2 * V_p_0(h_0, sys) * rho_0B /sys.μ* sys.GammaB /h_0 * f(h_0,sys) - sys.production_rate*rho_0/h_0 - sys.sigma_B_up / h_0
    return A
end


"""
    function V_gamma(h_0, sys::SysConstActive_1D)

Computes the function `V_gamma` for post-processing analysis based on linear stability analysis. This piecewise function varies according to `alpha` in `sys`, covering cases such as `alpha == -Inf`, `alpha < 0`, `alpha == 0`, `alpha > 0`, and `alpha == Inf`.

Returns:
- A calculated value based on the system parameters (`sys`) and height `h_0`.
"""
function V_gamma(h_0, sys::SysConstActive_1D)
    if sys.alpha == -Inf
        return sys.δ + h_0
    elseif sys.alpha < 0
        return (sys.δ * (exp(sys.alpha * h_0) - 1) + (exp(sys.alpha * h_0) - 1) / sys.alpha - h_0) / (exp(sys.alpha * h_0) - 1)
    elseif sys.alpha == 0
        return (2 * sys.δ + h_0) / 2
    elseif sys.alpha > 0 && sys.alpha < Inf
        return (-sys.δ * exp(-sys.alpha * h_0) + (1 - exp(-sys.alpha * h_0) * (sys.alpha * h_0 + 1)) / sys.alpha + sys.δ) / (1 - exp(-sys.alpha * h_0))
    else
        return sys.δ
    end
end

"""
    function V_gamma_0(h_0, sys::SysConstActive_1D)

Calculates `V_gamma` in the simplified case where `alpha = 0`, as:
V_gamma_0 = (2 * δ + h_0) / 2.

Returns:
- A baseline or default value used for comparison, based on `h_0` and `δ`.
"""
function V_gamma_0(h_0, sys::SysConstActive_1D)
    return (2 * sys.δ + h_0) / 2
end

"""
    function V_p(h_0, sys::SysConstActive_1D)

Computes the function `V_p`, defined as:
V_p = - (2 - 2αh - α^2 δ (δ + 2h) + exp(-αh) (α^2 (δ + h)^2 - 2)) / (2α^2(1 - exp(-αh)))

This is a post-processing function used in linear stability analysis. Similar to `V_gamma`, this function covers cases such as `alpha == -Inf`, `alpha < 0`, `alpha == 0`, `alpha > 0`, and `alpha == Inf`.

Returns:
- Computed `V_p` value based on `sys` parameters and height `h_0`.
"""
function V_p(h_0, sys::SysConstActive_1D)
    if sys.alpha == -Inf
        return (sys.δ + h_0)^2 / 2
    elseif sys.alpha < 0
        return -(-2 + sys.alpha^2 * (sys.δ + h_0)^2 - exp(sys.alpha * h_0) * (-2 + 2 * sys.alpha * h_0 + sys.alpha^2 * sys.δ * (sys.δ + 2 * h_0))) / (2 * sys.alpha^2 * (exp(sys.alpha * h_0) - 1))
    elseif sys.alpha == 0
        return Mob(h_0, sys)
    elseif sys.alpha > 0 && sys.alpha < Inf
        return -(2 - 2 * sys.alpha * h_0 - sys.alpha^2 * sys.δ * (sys.δ + 2 * h_0) + exp(-sys.alpha * h_0) * (sys.alpha^2 * (sys.δ + h_0)^2 - 2)) / (2 * sys.alpha^2 * (1 - exp(-sys.alpha * h_0)))
    else
        return sys.δ^2 / 2 + sys.δ * h_0
    end
end

"""
    function V_p_0(h_0, sys::SysConstActive_1D)

Provides the baseline or default `V_p` calculation when `alpha = 0`.

Returns:
- Computed `V_p` value using the simplified case where `alpha == 0`.
"""
function V_p_0(h_0, sys::SysConstActive_1D)
    return Mob(h_0, sys)
end

"""
    function ∇F(h_0, sys::SysConstActive_1D)

Computes the gradient function `∇F`, defined as:
∇F = (α * e^(αh) * ∂x h) / (e^(αh) - 1)

This function is used for post-processing and varies according to the value of `alpha` in `sys`. It considers cases such as `alpha == -Inf`, `alpha < 0`, `alpha == 0`, `alpha > 0`, and `alpha == Inf`.

Returns:
- Calculated gradient `∇F` based on system parameters and `h_0`.
"""
function ∇F(h_0, sys::SysConstActive_1D)
    if sys.alpha == -Inf
        return 0
    elseif sys.alpha < 0
        return -sys.alpha * exp(sys.alpha * h_0) / (exp(sys.alpha * h_0) - 1)
    elseif sys.alpha == 0
        return -1 / h_0
    elseif sys.alpha > 0 && sys.alpha < Inf
        return -sys.alpha / (exp(sys.alpha * h_0) - 1)
    else
        return 0
    end
end

"""
    function ∇F_0(h_0)

Computes the gradient function `∇F` in the case where `alpha = 0`.

Returns:
- Simplified gradient value, `-1 / h_0`.
"""
function ∇F_0(h_0)
    return -1 / h_0
end

"""
    function maxq(h_0, rho_0, rho_0A, rho_0B, sys::SysConstActive_1D)

Calculates the wavenumber where the largest eigenvalue of matrix `A` reaches its maximum. Iterates over values of `q` from 0 to 0.5 with a step of 0.0001, and finds the `q` that maximizes the eigenvalue of interest.

Returns:
- Wavenumber `q` where the maximum eigenvalue occurs.
"""
function maxq(h_0, rho_0, rho_0A, rho_0B, sys::SysConstActive_1D)
    lambdas = []
    step = 0.0001
    for q in 0:step:0.5
        lambda = eigvals(A(q, sys, rho_0=rho_0, h_0=h_0, rho_0A=rho_0A, rho_0B=rho_0B))[4]
        append!(lambdas, real(lambda))
    end
    return step * first(findall(x -> x == maximum(lambdas), lambdas)) - step
end

