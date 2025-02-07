### Author: Christoph Berke (2023-2025)
### Distributed under MIT license.

# Interactive simulation of N = 200 magnetic pendula.
# The number/strength/position of the magnets, the friction
# and the visual appearance is adjustable during runtime.

using OrdinaryDiffEq
using GLMakie
using DataStructures: CircularBuffer
using MutableNamedTuples
using ColorSchemes

fontsize_theme = Theme(fontsize = 40)
set_theme!(fontsize_theme)

# Implementation of the equations of motion of N (independent) magnetic pendula
# with gravitation and friction.
function Nmagneticpendulum(du, u, p, t)

	kg, gamma, h = p.params 
	Nmags = size(p.magposition,1)

	# Gravitation and friction
	for i in 1:p.Npendula
		du[4*(i-1)+1] = u[4*(i-1)+3]
		du[4*(i-1)+2] = u[4*(i-1)+4]
		du[4*(i-1)+3] = -kg*u[4*(i-1)+1] - gamma*u[4*(i-1)+3]
		du[4*(i-1)+4] = -kg*u[4*(i-1)+2] - gamma*u[4*(i-1)+4]

		# Magnetic forces.
		for k in 1:Nmags 
			dx = (u[4*(i-1)+1]-p.magposition[k,1])
			dy = (u[4*(i-1)+2]-p.magposition[k,2])
			d3 = sqrt(dx^2+dy^2+h^2)^3
			du[4*(i-1)+3] += -p.magstrength[k] * dx/d3
			du[4*(i-1)+4] += -p.magstrength[k] * dy/d3
		end 
	end 

	return 
end 


TOGGLEWIDTH = true

# Initial colors.
bgcol = Observable(:white)
fgcol = Observable(:black)

# Initial set of parameters. 
kg, gamma, h = 0.2, 0.3, 0.2

# Positions and strengths of magnets.
magpos = zeros(5,2)
magstr = 20 * ones(5)

for i in 1:5
	magpos[i,:] = [sin(i*2*pi/5), cos(i*2*pi/5)]
end

# Nearby initial conditions for Npendula magnetic pendula.
Npendula = 200 

u0N = zeros(4*Npendula)
for i in 1:Npendula
	u0N[4*(i-1)+1] = 1.5 + i*0.0002
	u0N[4*(i-1)+2] = 1.5 
end 


# Wrap all parameter in named tuple. 
paramsN = MutableNamedTuple(params = [kg,gamma,h], magstrength = magstr, 
	magposition = magpos, Npendula = Npendula)

# Integrate system by one time step. Update positions (balls) and
# earlier positions (trajs).
function animstepN!(integ, balls, trajs)
	step!(integ, 0.005, true)
	for i in 1:integ.p.Npendula 
		balls[][i] = Point2f(integ[4*(i-1)+1],integ[4*(i-1)+2])
		push!(trajs[i][], Point2f(integ[4*(i-1)+1],integ[4*(i-1)+2]))
	end
	notify.(trajs)
	balls[] = balls[]
end 

# Define ODE Problem.
tspan = (0.,10.)
probN = ODEProblem(Nmagneticpendulum, u0N, tspan, paramsN)
integ = init(probN, Tsit5())

# Arrays to current position (ballN) and trajectories. i.e., 'tail' earlier 
# positions (trajN).
ballN = Array{Point2f,1}(undef,Npendula)
trajN = []
tail = 200
for i in 1:Npendula
	ballN[i] =  Point2f(u0N[4*(i-1)+1],u0N[4*(i-1)+2])
	push!(trajN, CircularBuffer{Point2f}(tail))
end 

for i in 1:Npendula
	# Initial trajectory is only a point.
	fill!(trajN[i], ballN[i])
	# Promote arrays to observables.
	trajN[i] = Observable(trajN[i])
end 

# Make ballN Observable.
ballN = Observable(ballN)

# Generate main figure and axis.
fig = Figure()
display(fig) 
ax = Axis(fig[1:7,5:11], 
	backgroundcolor = bgcol,
	bottomspinecolor = fgcol,
	rightspinecolor = fgcol,
	topspinecolor = fgcol,
	leftspinecolor = fgcol,
	halign = :right
	)


# Observable for Colorlist for the individual pendula.
ColMap = Observable(:curl)
c = lift(ColMap) do cmp 
	colorrange = range(0, 1, length = Npendula)
	# Get color list from 'colorschemes[cmp]'' and numerical data 'colorrange'.
	colors = [get(colorschemes[cmp], colorrange[i]) for i in 1:Npendula]
	to_color.(colors)
end 

# Colors for trajectory: Basic color of each pendula taken from colorlist 'c'.
# Add fadeout in tail via alpha channel.
tailcol = [Observable([RGBAf(c[][j].r, c[][j].g, c[][j].b, 
	(i/tail)^2) for i in 1:tail]) for j in 1:Npendula]

# Adjust colors from trajectory when colorlist 'c' changes.
on(c) do c 
	for i in 1:Npendula
		tailcol[i][] = [RGBAf(c[i].r, c[i].g, c[i].b, (j/tail)^2) for j in 1:tail]
	end 
end 

# lineplot for trajectories.
for i in 1:Npendula
	lines!(ax, trajN[i]; linewidth = 3, color = tailcol[i])
end

# scatterplot for current position.
scatter!(ax, ballN; marker = :circle, strokewidth = 2, 
	  strokecolor = c,
	  color = c, markersize = 12
)

# Layout tweaking.
limits!(ax,-2,2,-2,2)
hidedecorations!(ax)

# Position of magnets.
corners = [Point2f(sin(phi), cos(phi)) for phi in range(0, 2*pi, 6)][1:end-1]
corners = Observable(corners)

Makie.deactivate_interaction!(ax, :rectanglezoom)
spoint = select_point(ax.scene, marker = :circle)

# Start / stop button.
run = Button(fig[7,2]; 
		label = "Start / Stop", 
		valign = :bottom,
		strokecolor = fgcol,
        strokewidth = 10,
        padding = (20,20,20,20),
        buttoncolor = bgcol,
        cornerradius = 10,
        buttoncolor_active = bgcol,
        labelcolor = fgcol,
        tellwidth = false 
		)

# Boolean that determines whether simulation is running.
isrunning = Observable(false)
# Change value of isrunning via click on run button.
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
# If isrunning = true, integrate system, update data.
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animstepN!(integ, ballN, trajN)
        sleep(0.001)
    end
end

# Pick new initial condition with mouse.
on(spoint) do z

	# New initial conditions
    x, y = z
    r = sqrt(x^2+y^2)
    (y > 0) ? (phi = acos(x/r)) : (phi = -acos(x/r))
    # phi = acos(x/r)
    phi_range = range(phi-0.005 *2*pi, phi+0.005 *2*pi, length = Npendula)

    u = zeros(4*Npendula)
	for i in 1:Npendula 
		u[4*(i-1)+1] = r * cos(phi_range[i])
		u[4*(i-1)+2] = r * sin(phi_range[i])
	end 

	# Reinitialize integrator with new initial condition.
    reinit!(integ, u)

    # Update arrays that hold current positions and trajectories.
    for i in 1:Npendula
		ballN[][i] =  Point2f(u[4*(i-1)+1],u[4*(i-1)+2])
	end 

	for i in 1:Npendula 
		for j in 1:tail
			trajN[i][][j] = Point2f(u[4*(i-1)+1],u[4*(i-1)+2])
		end 
	end 

	notify.(trajN)

	ballN[] = ballN[]

end

# Show magnets.
sccorn = scatter!(ax, 
		corners, 
		markersize = 100, color = (:black, 0), 
		strokecolor = fgcol, 
		strokewidth = 10)

# Toggle to display / hide magnets. 
toggle1 = Toggle(fig[4,3], active = true, tellwidth = true)
connect!(sccorn.visible, toggle1.active)

# Slider grid to change parameters.
sg = SliderGrid(fig[1, 1:4],
    (label = "Reibungskraft", range = 0:0.1:2, startvalue = 0.3, format = "{:.1f}", linewidth = 30),
    (label = "Magnetstärke", range = 0:1:40, startvalue = 20, format = "{:.1f}", linewidth = 30),
    (label = "Abstand der Magnete", range = 0:0.1:1, startvalue = 1, format = "{:.1f}", linewidth = 30),
    (label = "Anzahl der Magnete", range = 1:1:10, startvalue = 5, format = "{:.1f}", linewidth = 30)
)

# Change parameters of integrator when associated quantity is adjusted via the 
# slider

# Friction
on(sg.sliders[1].value) do fric 
    integ.p.params[2] = fric 
end

# Strength of magnets.
on(sg.sliders[2].value) do magstr 
    integ.p.magstrength[:] = magsign[] * magstr * ones(sg.sliders[4].value[]) 
end

# Position of magnets.
on(sg.sliders[4].value) do nmgs 

	magpos = zeros(nmgs,2)

	for i in 1:nmgs
		magpos[i,:] = [sin(i*2*pi/nmgs), cos(i*2*pi/nmgs)]
	end
    
    integ.p.magstrength = magsign[] * sg.sliders[2].value[] * ones(nmgs)
    integ.p.magposition = sg.sliders[3].value[]*magpos

    corners[] = [Point2f(sg.sliders[3].value[]*magpos[i,1], sg.sliders[3].value[]*magpos[i,2]) for i in 1:nmgs]
end

on(sg.sliders[3].value) do rmgs 
	nmgs = sg.sliders[4].value[]
	magpos = zeros(nmgs,2)
	for i in 1:nmgs
		magpos[i,:] = rmgs * [sin(i*2*pi/nmgs), cos(i*2*pi/nmgs)]
	end
    integ.p.magposition[:] = magpos
    corners[] = [Point2f(magpos[i,1], magpos[i,2]) for i in 1:nmgs]
end


# Button to reset axis limits.
home = Button(fig[7, 4], label = "Home", halign = :left, valign = :bottom, tellwidth = false,
	strokecolor = fgcol,
        strokewidth = 10,
        padding = (20,20,20,20),
        buttoncolor = bgcol,
        cornerradius = 10,
        buttoncolor_active = bgcol,
        labelcolor = fgcol)

on(home.clicks) do _
	limits!(ax,-2,2,-2,2)
end

# Toggle to switch between attractive and repulsive force.
toggle3 = Toggle(fig[3,3], active = true, tellheight = false, tellwidth = true )

on(toggle3.active) do _ 
	integ.p.magstrength[:] .*= -1 #magstr
end 

magsign = lift(toggle3.active) do z
	2 * Int(z) - 1
end 

# Toggle to switch between dark and bright layout.
toggle2 = Toggle(fig[5,3], active = true, tellheight = false)

on(toggle2.active) do tog 
	if tog 
		bgcol[] = :white 
		fgcol[] = :black
		ColMap[] = :curl
		fig.scene.backgroundcolor[] = colorant"white"
		for i in 1:4
			sg.labels[i].color[] = colorant"black"
			sg.valuelabels[i].color[] = colorant"black"
		end 
	else 
		bgcol[] = :black
		fgcol[] = :white 
		ColMap[] = :coolwarm
		fig.scene.backgroundcolor[] = colorant"black"
		for i in 1:4
			sg.labels[i].color[] = colorant"white"
			sg.valuelabels[i].color[] = colorant"white"
		end 
	end 
end 


# Add labels.
Label(fig[3,1],  text = "Magnetkraft:", color = fgcol, tellheight = true, tellwidth = false, halign = :left)
Label(fig[3,2],  text = "abstoßend", color = fgcol, tellheight = false, tellwidth = false, halign = :right)
Label(fig[3,4],  text = "anziehend", color = fgcol, tellheight = false, tellwidth = false, halign = :left)

Label(fig[4,1],  text = "Magnete:", color = fgcol, tellheight = true, tellwidth = false, halign = :left)
Label(fig[4,2],  text = "ausblenden", color = fgcol, tellheight = false, tellwidth = true, halign = :right)
Label(fig[4,4],  text = "anzeigen", color = fgcol, tellheight = false, tellwidth = false, halign = :left)

Label(fig[5,1],  text = "Lichtschalter:", color = fgcol, tellheight = true, tellwidth = false, halign = :left)
Label(fig[5,2],  text = "aus", color = fgcol, tellheight = false, tellwidth = false, halign = :right)
Label(fig[5,4],  text = "an", color = fgcol, tellheight = false, tellwidth = false, halign = :left)

# Layout tweaking.
for i in 1:7
    colsize!(fig.layout, i+4, Aspect(i, 1.0))
end