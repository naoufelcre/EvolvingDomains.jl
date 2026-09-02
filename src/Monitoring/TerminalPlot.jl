module TerminalPlot

using ..Geometric: EvolvingDiscreteGeometry, grid_info

export plot

"""
    plot(geom::EvolvingDiscreteGeometry; io=stdout, size=nothing, label=nothing)

Render the negative level-set region with Unicode block characters. Repeated calls in a
terminal replace the previous frame and adapt to the current terminal size. `label`
adds a line above the plot, for example `label="Iteration 100"`.
"""
function plot(geom::EvolvingDiscreteGeometry; io::IO=stdout,
    size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing)
    rows, columns = isnothing(size) ? displaysize(io) : size
    label = isnothing(label) ? nothing : first(split(string(label), '\n'))[1:min(end, columns)]
    width = max(columns - 2, 1)
    height = max(rows - 2 - !isnothing(label), 1)
    nx, ny = grid_info(geom.grid).dims
    nx > 1 && ny > 1 || throw(ArgumentError("plot requires at least 2 nodes in each dimension"))

    # Most terminal cells are roughly twice as tall as they are wide.
    width = min(width, max(1, round(Int, 2 * height * (nx - 1) / (ny - 1))))
    height = min(height, max(1, round(Int, width * (ny - 1) / (2 * (nx - 1)))))
    phi = reshape(geom.levelset, nx, ny)

    inside(column, halfrow) = any(phi[i, j] <= 0
        for i in fld((column - 1) * nx, width) + 1:fld(column * nx - 1, width) + 1,
            j in fld((halfrow - 1) * ny, 2height) + 1:fld(halfrow * ny - 1, 2height) + 1)

    frame = IOBuffer()
    tty = io isa Base.TTY || get(io, :terminal, false)
    clearline = tty ? "\e[K" : ""
    isnothing(label) || println(frame, label, clearline)
    println(frame, '+', repeat("-", width), '+', clearline)
    for row in height:-1:1
        print(frame, '|')
        for column in 1:width
            upper = inside(column, 2row)
            lower = inside(column, 2row - 1)
            print(frame, upper ? (lower ? '█' : '▀') : (lower ? '▄' : ' '))
        end
        println(frame, '|', clearline)
    end
    print(frame, '+', repeat("-", width), '+', clearline, tty ? "\e[J" : "")
    tty && print(io, "\e[H")
    write(io, take!(frame))
    flush(io)
    return nothing
end

end # module
