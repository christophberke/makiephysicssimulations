### Author: Christoph Berke (2023-2025)
### Distributed under MIT license.

# Simple demo of a single magnetic pendulum moving 
# above three magnets.

using DifferentialEquations
using GLMakie
using DataStructures: CircularBuffer

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

# Perform one integration step of the ODE and update the Makie observables 
# that visualize the current position of the pendulum (balls) and the past
# trajectory in form of the outfading tail (trajs)
function animstep!(integ, balls, trajs)
	
	# Progress the integrator integ by one time step.
	step!(integ, 0.005, true)
	# Update the Makie observables balls (for the current pendula positions) ...
	balls[] = Point2f(integ[1],integ[2])
	# ... and push the current position to the pendulum "tail".
	push!(trajs[], Point2f(integ[1],integ[2]))
	# Notify the observable trajs that its value has been updated.
	notify(trajs)

end 

# Set fontsize for Makie labels and buttons.
fontsize_theme = Theme(fontsize = 40)
set_theme!(fontsize_theme)

# Generate three magnets 
magpos = zeros(3,2) 	# Array for positions of magnets
magstr = 0. * ones(3) 	# Strength of magnets 

# Fill positions: Magnets are placed on the corners of a regular triangle. 
for i in 1:3
	magpos[i,:] = [sin(i*2*pi/3), cos(i*2*pi/3)]
end

kg, gamma, h = 1, 0.3, 0.2	# gravitation constant, friction and height
u0N = [1.5,1.5,0,0]			# Initial conditions
Npendula = 1 				# Only a single pendulum
tspan = (0.,10.) 			# Evolution time span

# Wrap all parameters in a named tuple.
paramsN = (params = [kg,gamma,h], magstrength = magstr, 
	magposition = magpos, Npendula = Npendula)

# Define ODE problem and ODE integrator.
probN = ODEProblem(Nmagneticpendulum, u0N, tspan, paramsN)
integ = init(probN, Tsit5())


tail = 1000							# Length of the trajectory tail
ball = Point2f(u0N[1],u0N[2])		# Current position of pendulum
traj = CircularBuffer{Point2f}(tail)# Array for trajectory
fill!(traj, ball)

# Promote current position and trajectory to observables.
traj = Observable(traj)
ball = Observable(ball)

# Create and display figure and main axis environment.
fig = Figure(size = (800,800))
display(fig) 
ax = Axis(fig[1,1], aspect = DataAspect())

# Create RGBA array with decreasing alpha (opacity) value for fadeout of traj. 
c = to_color(:royalblue3)
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]

# Plot current position and trajectory. 
lines!(ax, traj; linewidth = 3, color = tailcol)
scatter!(ax, ball; marker = :circle, strokewidth = 2, 
	strokecolor = c, color = c, markersize = 25)

limits!(ax,-2,2,-2,2)
hidedecorations!(ax)

# Observable for positions of magnets. 
corners = [Point2f(sin(phi), cos(phi)) for phi in range(0, 2*pi, 4)][1:end-1]
corners = Observable(corners)

# GridLayout for buttons and labens below main axis.
gl = GridLayout(fig[2,1], tellwidth = false)

# Start / stop button for animation.
run = Button(gl[1,1]; label = "Start / Stop", 
	valign = :bottom, halign = :left, strokecolor = :black, strokewidth = 10, 
	padding = (20,20,20,20), tellwidth = false)

isrunning = Observable(false)	# Observable to start / stop animation.
# Change value of isrunning by clicking on run button.
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end

# Run animation if isrunning == true
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break 		# stop computation if figure was closed.
        animstep!(integ, ball, traj) 	# perform one integration step.
        sleep(0.001) 
    end
end

# Observable to select initial conditions by clicking.
spoint = select_point(ax.scene, marker = :circle)

# Connect clicking with restart of ODE solver with new initial conditions.
on(spoint) do z

	# Set current system positions to selected point and momenta to zero.
    x, y = z
    u = [x,y,0,0]
    reinit!(integ, u) 				# Restart ODE solver 
	ball[] =  Point2f(u[1],u[2]) 	# Reset obersable for current position

	# Reset observable for trajectory
	for j in 1:tail
		traj[][j] = Point2f(u[1],u[2])
	end 
	
	# Notify observables of their value update.	
	notify(traj)
	ball[] = ball[]

end

# Switch off zoom
Makie.deactivate_interaction!(ax, :rectanglezoom)
scatter_magnets = scatter!(ax, corners, markersize = 100, color = :goldenrod2)

# Toggle to switch magnets on / off.
Label(gl[1,2], "Magnets on / off", halign = :right)
magnets_switch = Toggle(gl[1,3], active = false, tellwidth = false)
connect!(scatter_magnets.visible, magnets_switch.active)

# Switch magnetes on and off by updating p.magstrength 
# (and also the gravitation constant, just for visualization purposes.)
on(magnets_switch.active) do tog 
	if tog 
		integ.p.magstrength[:] = 20 * ones(3)
		integ.p.params[1] = 0.2
	else 
		integ.p.magstrength[:] = zeros(3)
		integ.p.params[1] = 1
	end 
end 

# Layout tweaking.
colsize!(fig.layout, 1, Aspect(1, 1))
colgap!(gl, 2, Relative(.05))
colsize!(gl, 1, Auto())
colsize!(gl, 3, Auto())
colsize!(gl, 2, Relative(1))