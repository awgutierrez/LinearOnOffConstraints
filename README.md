# LinearOnOffConstraints
An example that illustrates the use of [JuMP Indicator Constraints](https://jump.dev/JuMP.jl/stable/manual/constraints#Indicator-constraints) 
and the Julia package [PiecewiseAffineApprox.jl](https://github.com/sintefore/PiecewiseAffineApprox.jl) to transform a nonlinear constraint to a set of linear on/off constraints.

```julia
using PiecewiseAffineApprox, JuMP, HiGHS

include("ElecGas.jl")
```

# Dealing with a nonlinear constraint in a coupling (e,g) model

## :triangular_flag_on_post: The problem
Transform a constraint of the form

$$f |f| = p_1^2 - p_2^2$$

to a linear on/off constraint.

## :bulb: A possible solution
Use a combination of 
- Julia package [PiecewiseAffineApprox.jl](https://github.com/sintefore/PiecewiseAffineApprox.jl)
- [JuMP Indicator Constraints](https://jump.dev/JuMP.jl/stable/manual/constraints#Indicator-constraints)

## :beginner: A simple example
Minimize the linear function

$$C^e \cdot e + C^g \cdot g$$

over a set of constraints containing

$$\boxed{\quad f |f| = p_1^2 - p_2^2 \quad}$$

$$2f = f^I+f^O$$

$$1\cdot g - \delta \cdot e - f^I + f^O = D^g$$

$$\ell = \gamma (p_1 + p_2)/2 = \ell_0 + f^I - f^O$$

$$1\cdot e - h = D^e$$

$$p_2 \leq \Gamma p_1$$ 

$$1\cdot ------ \cdot  2$$

## :feet: Change of variable
$$f |f| = p_1^2 - p_2^2 = (p_1 - p_2)(p_1 + p_2)$$

$$u:=p_1 - p_2 \quad\quad v:=p_1 + p_2$$

$$f \geq 0 \iff u\geq 0\quad\mbox{ : in this case } \quad f = \sqrt{uv}$$

## :feet: Approximation

Use Julia package [PiecewiseAffineApprox.jl](https://github.com/sintefore/PiecewiseAffineApprox.jl) to approximate the function

$$F(u,v) = \sqrt{uv}$$

in a suitable domain. 

## :feet: Reformulation

Decompose the variable $f$ into its positive and negative part

$$f = f^+ - f^-$$

Use 
[JuMP Indicator Constraints](https://jump.dev/JuMP.jl/stable/manual/constraints#Indicator-constraints)
to formulate the following

$$u\geq 0 \iff f^+ = F(u,v) = \sqrt{uv}$$

$$u\leq 0 \iff f^- = F(-u,v) = \sqrt{-uv}$$


## :feet: Implementation in Julia

```julia
PWA = PiecewiseAffineApprox
optimizer = HiGHS.Optimizer
```
```julia
@info "Approximating the nonlinear term..."
# Approximate the positive part and negative part of f
planes = 8
pwa1 = PWA.approx(
    x -> sqrt(x[1]*x[2]),
    [(0, max(P̅₁ - P̲₂, P̅₂ - P̲₁)) , (0 , P̅₁ + P̅₂)], 
    Concave(),
    Cluster(optimizer=optimizer, planes=planes, metric= :l1))
@info "Approximation done"
```
```julia
@info "Adding the approximation to the model"
@variable(m, f)
@variable(m, 0 ≤ f⁺ ≤ PWA.evaluate(pwa1, [P̅₁ - P̲₂ P̅₁ + P̅₂]))
@variable(m, 0 ≤ f⁻ ≤ PWA.evaluate(pwa1, [P̅₂ - P̲₁ P̅₁ + P̅₂]))
@constraint(m, f == f⁺ - f⁻)

for p in pwa1.planes
	@constraint(m, b --> { f⁺ <= -sum(p.α .* [u v]) - p.β })
	@constraint(m, !b --> { f⁻ >= -sum(p.α .* [u v]) - p.β })
end
@info "Approximated linear constraint done"
```
