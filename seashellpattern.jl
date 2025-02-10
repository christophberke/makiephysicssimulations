### Author: Christoph Berke (2023-2025)
### Distributed under MIT license.

### This simulation enables one to explore one-dimensional cellular automata
### and how the simple update rules contribute to the formation of the complex 
### patterns found on some exotic sea shells.

# Load packages
using GLMakie
using FileIO
using Colors 


#################################
# A number of useful functions. #
#################################

# Initial conditions for the first row of the sea shell pattern. init_xxx(N) 
# returns a BitArrays of length N. A '1' ('0') encodes a dark (bright) cell.

# A single dark cell in the middle of the array
function init_single(N)
    x = BitArray(zeros(N))
    x[round(Int, N/2)+1] = 1
    return x 
end

# Two dark cells 
function init_two(N)
    x = BitArray(zeros(N))
    x[round(Int, N/3)] = 1
    x[round(Int, 2*N/3)] = 1
    return x 
end

# Three dark cells
function init_three(N)
    x = BitArray(zeros(N))
    x[round(Int, N/4)] = 1
    x[round(Int, N/2)+1] = 1
    x[round(Int, 3*N/4)] = 1
    return x
end

# Random distribution of dark and bright cells
# (where bright and dark cells occur with equal probability)
function init_rand(N)
    x = BitArray([rand([true,false]) for _ in 1:N])
    return x 
end


# ruleX(a, X) updates the array a according to rule X and returns the result.
# The implementation assumes fixed boundary conditions (fbc), for periodic 
# boundary conditions (pbc) see the comments in the code below.
function ruleX(a::BitArray, X::Integer)

    anew = similar(a)   # Array to store the Array a after one update step.
    # a = CircularArray(a)  # Convert a to circular array for pbc

    # Translate integer X to update rule (see e.g. sec. "The numbering system" 
    # in https://en.wikipedia.org/wiki/Elementary_cellular_automaton)
    newstate = BitArray(digits(X, base=2, pad = 8))

    # Update all inner entries (fbc)
    for i in 2:length(a)-1 # change to 1:length(a)-1 for pbc.
        x = a[i-1:i+1]  # The three neighboring cells of the previous row "a"
                        # that determine the color of the next row "anew".
        
        # Loop over the eight possible color configurations of the three cells x
        for j in 0:7
            # Find the correct color configuration and store updated cell color.
            if x == digits(j, base = 2, pad = 3)
                anew[i] = newstate[j+1]
                break
            end
        end
    end

    # Boundaries (fbc, comment for pbc)
    anew[1], anew[end] =  a[1], a[end]

    # Return updated cell color pattern.
    return anew

end



##########################
# Visualization Settings #
##########################

# Primary color (orange to mimic the actual seashell color)
coltheme = "Oranges"
cmap = colormap(coltheme)

fontsize_theme = Theme(fontsize = 40)

# Settings for all buttons.
button_theme = Theme(
    Button = (
        strokecolor = cmap[end],
        strokewidth = 10,
        padding = (20,20,20,20),
        buttoncolor = cmap[10],
        cornerradius = 10,
        buttoncolor_active = cmap[end],
        buttoncolor_hover = cmap[30]
    )
)

# Setting for dropdown menu.
menu_theme = Theme(
    Menu = (
        dropdown_arrow_size = 40, 
        textpadding = (15,15,15,15), 
        fontsize = 36,
        cell_color_active = cmap[end-10],
        cell_color_hover = cmap[30],
        cell_color_inactive_even = cmap[10],
        cell_color_inactive_odd = cmap[10],
        selection_cell_color_inactive = cmap[10]
    )
)

final_theme = merge(button_theme, fontsize_theme, menu_theme)
set_theme!(final_theme)

# Initialize a figure.
fig = Figure(backgroundcolor = cmap[1], size = (1600, 900))

# Add example images of sea shell patterns and a title.
isdir("makiephysicssimulations") ? path="makiephysicssimulations/" : path=""
imshell = rotr90(load(path*"seashells.png"))
aximag = Axis(fig[1,1:4], tellheight = true, 
    halign = :right, aspect = DataAspect())
image!(aximag, imshell)
hidedecorations!(aximag)
hidespines!(aximag)
Label(fig[2,1:4], "Wie kommt das Muster auf die Schnecke?", 
    tellheight = false, halign = :center, fontsize = 56)


#################################################################
# Main axis environment for visualization of sea shell pattern. #
#################################################################

# Optional: Show a grid for the smallest array size to facilitate understanding
# of the update rules, akin to what is shown in presentation.pdf
showgrid = Observable(false)
grid_x, grid_y = 1.5:1:19.5, 1.5:1:20.5

# Axis environment for 2D visualization of sea shell / 1D cellular automaton.
ax = Axis(fig[1:9,5:13], 
    aspect = 1, 
    yreversed = true, 
    yminorgridvisible = showgrid, 
    yminorticks = grid_y,
    xminorgridvisible = showgrid, 
    xminorticks = grid_x,
    yminorgridcolor = :black, 
    xminorgridcolor = :black,
    halign = :right
)

hidedecorations!(ax) # Hide decorations, minor grid can later be adjusted.


###############################################
# Sliders and button to adjust sea shell size #
###############################################

# Slider to adjust number of cells in x- and y-direction.
sg = SliderGrid(fig[3, 1:4],
    (label = "Breite", range = 21:10:1021, startvalue = 100, 
        format = x -> string(x, " Zellen"), linewidth = 30, 
        color_active = cmap[end], 
        color_active_dimmed = cmap[20]),
    (label = "Höhe", range = 21:10:1001, startvalue = 100, 
        format = x -> string(x-2, " Zellen"), linewidth = 30, 
        color_active = cmap[end], 
        color_active_dimmed = cmap[20])
)

# Adjust x limits to slider value
on(sg.sliders[1].value) do v 
    xlims!(ax,[0.5,v-1.5])
end

# Adjust y limits to slider value
on(sg.sliders[2].value) do v 
    ylims!(ax,[v+0.5,0.5])
end

# Button to reset x- and y-limits (e.g., after zooming in).
home = Button(fig[7, 4], label = "Home", 
    tellheight = false, tellwidth = false, halign = :right)

# Connect click on button with reset of axis limits to slider values 
on(home.clicks) do _
    xlims!(ax,[0.5,sg.sliders[1].value[]-1.5])
    ylims!(ax,[sg.sliders[2].value[]+0.5,0.5])
end


###########################################################################
# Axis environment and settings for interactive update of automaton rule. #
###########################################################################

# Idea: Display a 2D bit array that describes the update rule in the form
# 000 001 010 011 100 101 110 111  <---- current pattern. 
#  0   0   1   0   1   0   1   0   <---- new pattern.
# (0 = bright, 1 = dark)
# First row: The eight possible patterns of three adjacent cells.
# Second row: define update rule. A rule can be selected by clicking in the 
# array cells of the second row. Each click changes a 0 to 1 and vice versa. 


# Array for the visualization of the update rules in the above described form
# The NaNs add empty spaces where necessary. 
# (Note that each triple blocks in the first line is mirrored compared to the 
# standard representation (e.g., 100 instead of 001). This has no deeper meaning
# and is only due to the way, the arrays are displayed (reversed ydir, etc.) 
updaterule = NaN * zeros(2,31)
updaterule[1,:] = [
    0 0 0 NaN 1 0 0 NaN 0 1 0 NaN 1 1 0 NaN 0 0 1 NaN 1 0 1 NaN 0 1 1 NaN 1 1 1]
for i in 2:4:30
    updaterule[2,i] = 0
end 


# Observable for selecting and reading the update rule.
Oupdaterule = Observable(updaterule')

# Extract rule from Oupdaterule whenever 
# an entry of Oupdaterule changes its value.
ruleno = lift(Oupdaterule) do O 
    # Binary representation of rule == digits in second row.
    digs = Oupdaterule[][2:4:end,2]
    # Convert binary representation to integer.
    x = sum(digs[k]*2^(k-1) for k=1:length(digs))
    Int(x)
end

# Axis environment to display the array Oupdaterule
ax2 = Axis(fig[9,1:4], aspect = DataAspect(), 
    yreversed = true, backgroundcolor = :gray80)

# Show array
hm = heatmap!(ax2, Oupdaterule, colormap = coltheme)

# Add thick black border to improve readibility.
for i in 0:4:28
    lines!(ax2, i .+ [0.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,0.5,0.5], 
        [0.5,1.5,1.5,2.5,2.5,1.5,1.5,0.5,0.5,1.5], 
        color = :black, linewidth = 5)
    lines!(ax2, i .+ [1.5,1.5,2.5,2.5], [0.5,1.5,1.5,0.5], 
        color = :black, linewidth = 5)
end

# Beautification of axis environment.
limits!(ax2, 0.25,31.75,2.75,0.25)
hidespines!(ax2)
hidedecorations!(ax2)

# Switch off zooming with mouse.
Makie.deactivate_interaction!(ax2, :rectanglezoom)

# Select a point by clicking left mouse button
spoint = select_point(ax2.scene, marker = :circle)

# Connect selection of point with update of Oupdaterule.
on(spoint) do z
    # Round selected point to get the cell number that contains the point.
    x, y = round.(Int, z)
    # Update of cell entry only if 
    # (i) the cell is in the second row 
    # (ii) the cell entry is not NaN.
    if (y == 2) && (mod(x+2,4)==0)
        # Update cell entry: 0 <--> 1
        Oupdaterule[][x,y] = mod(Oupdaterule[][x,y]+1,2)
    end
    # Notify Oupdaterule of its update.
    Oupdaterule[] = Oupdaterule[]
end


##############################################################################
# Dropdown menu to choose between a selection of predefined automaton rules. #
##############################################################################

# Description of dropdown menu.
Label(fig[7,1], "Update Regel", tellheight = true, halign = :left)

# Selection of predefined rules
rules = [0,5,18,22,30,45,57,60,73,90,99,110,126]

# Labels for dropdown menu
rulelabels = "Regel " .* string.(rules)

# Actual menu
menu = Menu(fig[7,2], tellwidth = true, halign = :right, 
    options = zip(rulelabels, rules))

# Connect selection of a predefined rule with an update of Oupdaterule.
on(menu.selection) do s 
    # ruleno[] = s     # Regel anpassen (Brauch man diese Zeile?)

    # Get binary representation of update rule.
    newstate = digits(s, base=2, pad = 8)
    # Update the corresponding entries (i.e., second row, every fourth column,
    # starting with column 2)
    for i in 1:8
        Oupdaterule[][4*(i-1)+2,2] = newstate[i]
    end
    # Notify Oupdaterule of its update.
    Oupdaterule[] = Oupdaterule[]
end


########################################################
# Settings for the first row of the sea shell pattern. #
########################################################

# Buttons to choose between different number of dark cells in first row.
fig[5, 2:4] = buttongrid = GridLayout(tellwidth = true, halign = :right)
buttonlabels = ["1", "2", "3", "zufällig"]
buttons = buttongrid[1, 1:4] = [Button(fig, label = l) for l in buttonlabels]

# Add description to buttongrid.
Label(fig[5,1:2], "Dunkle Zellen in 1. Reihe:", 
    tellheight = true, tellwidth = true, halign =:left)

# The array inits contains the above defined functions init_xxx that implement
# the desired number (1,2,3,random) of dark cells in the first row.
inits = [init_single,init_two,init_three,init_rand]

# The observable func holds the currently selected first-row condition.
func = Observable{Any}(inits[1]) 

# Connect the click on one of the buttons with an update of the observable func.
for i in 1:4
    on(buttons[i].clicks) do _
        func[] = inits[i]
    end
end

# Observable for the first row of the sea shell pattern that is updated when
# (i) the width of the sea shell changes by tuning the corresponding slider.
# (ii) the initial condition / number of dark cells changes.
# The observable is not updated when the automaton rule changes (as desired).
firstrow = @lift( $func($(sg.sliders[1].value)) )


#########################
# The actual sea shell. #
#########################

# A is the two-dimensional BitArray that is updated whenever 
# (i) The height or width is adjusted via the slider.
# (ii) The observable firstrow (number of dark cells in row 1) is adjusted.
# (iii) The update rule ruleno changes its value.
A = @lift begin 
    # BitArray with dimension given by slider values for width and height.
    x = BitArray(zeros($(sg.sliders[1].value),$(sg.sliders[2].value)))

    # First row determined by obersable first row.
    x[:,1] = $firstrow 

    # Update all other rows. The lower rows i>1 contain the result of the 
    # application of rule ruleno to the previous row i-1.
    for i in 2:$(sg.sliders[2].value)
        x[:,i] = ruleX(x[:,i-1], $ruleno)
    end 

    x[2:end-1,:] # Discard first and last column (never change due to fixed bc).
end

# Show sea shell pattern as heatmap.
hm1 = heatmap!(ax, A, colormap = coltheme)

# Optional: Translate heatmap in negative z-direction to make grid visible.
translate!(hm1, 0, 0, -15) 



##############################
# Final layout adjustements. #
##############################

# Adjust the width of columns 5 to 13 to the height of the rows 1 to 9.
# As the main axis is contained in fig[1:9,5:13], this results in a square
# image for the sea shell.
for i in 1:9
    colsize!(fig.layout, i+4, Aspect(i, 1.0))
end

# Thin colored boxes as separation elements between buttons / sliders / menu
# in row 4, 6, and 8
for i in 4:2:8
    rowsize!(fig.layout, i, 8)  # Fix height of row.
    # Add colored box
    Box(fig[i,1:4], color = cmap[end], strokewidth = 0)
end

display(fig)
