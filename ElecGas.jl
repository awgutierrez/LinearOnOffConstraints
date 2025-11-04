PWA = PiecewiseAffineApprox
optimizer = HiGHS.Optimizer


@info "Building the model"
Emin = [0 0]
Emax = [45 75]

Gmin = [17.5 35.7]
Gmax = [180 110]

m = Model(optimizer)

@variable(m, Emin[i] ≤ e[i=1:2] ≤ Emax[i])
@variable(m, Gmin[i] ≤ g[i=1:2] ≤ Gmax[i])

P̲₁ , P̅₁ = 50, 900
P̲₂ , P̅₂ = 50, 800

@variable(m, P̲₁ ≤ p₁ ≤ P̅₁)
@variable(m, P̲₂ ≤ p₂ ≤ P̅₂)

@variable(m, b, Bin)
Γ₁₂ = 1.2
Γ₂₁ = 1.5

@constraint(m, p₂ - Γ₂₁*p₁ ≤ 0)
@constraint(m, p₂ - Γ₁₂*p₁ ≥ 0)

# Change of variable
@variable(m, 0 ≤ u ≤ max(P̅₁ - P̲₂, P̅₂ - P̲₁))
@variable(m, 0 ≤ v ≤ P̅₁ + P̅₂)
@constraint(m, v == p₁ + p₂)
@constraint(m, b --> { u == p₁ - p₂ })
@constraint(m, !b --> { u == p₂ - p₁ })

@info "Approximating the nonlinear term..."
# Approximate the positive part and negative part of f
planes = 8
pwa1 = PWA.approx(
    x -> sqrt(x[1]*x[2]),
    [(0, max(P̅₁ - P̲₂, P̅₂ - P̲₁)) , (0 , P̅₁ + P̅₂)], 
    Concave(),
    Cluster(; optimizer = optimizer, planes = planes, metric = :l1),
)
@info "Approximation done"

@info "Adding the approximation to the model"
@variable(m, f)
@variable(m, 0 ≤ f⁺ ≤ PWA.evaluate(pwa1, [P̅₁ - P̲₂ P̅₁ + P̅₂]))
@variable(m, 0 ≤ f⁻ ≤ PWA.evaluate(pwa1, [P̅₂ - P̲₁ P̅₁ + P̅₂]))
@constraint(m, f == f⁺ - f⁻)

#######################################################
# Add the piecewise linear approximation to the model #
for p in pwa1.planes
	@constraint(m, b --> { f⁺ <= -sum(p.α .* [u v]) - p.β })
	@constraint(m, !b --> { f⁻ >= -sum(p.α .* [u v]) - p.β })
end
#######################################################
@info "Approximated linear constraint done"

# Gas flow as the average of the inlet and outlet flows
@variable(m, f⁺ᴵ ≥ 0)
@variable(m, f⁺ᴼ ≥ 0)
@variable(m, f⁻ᴵ ≥ 0)
@variable(m, f⁻ᴼ ≥ 0)
@constraint(m, 2*f⁺ == f⁺ᴵ + f⁺ᴼ)
@constraint(m, 2*f⁻ == f⁻ᴵ + f⁻ᴼ)

# Linepack
@variable(m, ℓ ≥ 0)
γ = 2.1
ℓ₀ = 10 
@constraint(m, 2*ℓ == γ*(p₁+p₂))
@constraint(m, ℓ == ℓ₀ + f⁺ᴵ - f⁺ᴼ + f⁻ᴵ - f⁻ᴼ)

# Balance constraint for the coupling (e,g)
Dᵍ = 800
δ = [0.5 0.3]
@constraint(m, sum(g) - sum(δ.*e) - f⁺ᴵ + f⁻ᴼ == Dᵍ)

# Balance for the variable e
Dᵉ = 100
@variable(m, -10 ≤ h ≤ 10)
@constraint(m, sum(e) - h == Dᵉ)

# objective
Cᵉ = [1.4 3.7]
Cᵍ = [1.2 1.4]
@objective(m, Min, sum(Cᵉ.*e) + sum(Cᵍ.*g))

@info "Model ready"

println(m)

@info "Solving the model ..."
optimize!(m);
@info "Model solved"

@info "Showing optimal values ..."
println("e = ", value(e))
println("g = ", value(g))
println("f = ", value(f))
println("p = ", value([p₁ p₂]))
println("ℓ = ", value(ℓ))
