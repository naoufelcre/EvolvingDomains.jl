module TerminalPlot

using ..Geometric: EvolvingDiscreteGeometry, grid_info

export plot

# ANSI approximations of Julia's purple, blue, green, and red logo colors.
const FIELD_COLORS = (54, 55, 56, 57, 63, 69, 75, 81, 80, 79, 78, 77, 76,
    112, 148, 184, 220, 214, 208, 202, 196)

"""
    plot(geom::EvolvingDiscreteGeometry; io=stdout, size=nothing, label=nothing,
         field=nothing, colorrange=nothing)

Render the negative level-set region with Unicode block characters. Repeated calls in a
terminal replace the previous frame and adapt to the current terminal size. `label`
adds a line above the plot, for example `label="Iteration 100"`.
Pass a nodal scalar `field` to color the interior in ANSI terminals.
"""
function plot(geom::EvolvingDiscreteGeometry; io::IO=stdout,
    size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing,
    field::Union{Nothing,AbstractVector{<:Real}}=nothing, colorrange=nothing)
    rows, columns = isnothing(size) ? displaysize(io) : size
    tty = io isa Base.TTY || get(io, :terminal, false)
    label = isnothing(label) ? nothing : first(split(string(label), '\n'))[1:min(end, columns)]
    width = max(columns - 2 - (tty && !isnothing(field) ? 12 : 0), 1)
    height = max(rows - 2 - !isnothing(label), 1)
    nx, ny = grid_info(geom.grid).dims
    nx > 1 && ny > 1 || throw(ArgumentError("plot requires at least 2 nodes in each dimension"))

    # Most terminal cells are roughly twice as tall as they are wide.
    width = min(width, max(1, round(Int, 2 * height * (nx - 1) / (ny - 1))))
    height = min(height, max(1, round(Int, width * (ny - 1) / (2 * (nx - 1)))))
    phi = reshape(geom.levelset, nx, ny)
    isnothing(field) || length(field) == length(phi) ||
        throw(DimensionMismatch("field has $(length(field)) values, expected $(length(phi))"))
    values = isnothing(field) ? nothing : reshape(field, nx, ny)
    finite = isnothing(field) ? nothing : filter(isfinite, field)
    isnothing(finite) || !isempty(finite) || throw(ArgumentError("field has no finite values"))
    lo, hi = isnothing(values) ? (0.0, 1.0) :
        (isnothing(colorrange) ? extrema(finite) : colorrange)

    inside(column, halfrow) = any(phi[i, j] <= 0
        for i in fld((column - 1) * nx, width) + 1:fld(column * nx - 1, width) + 1,
            j in fld((halfrow - 1) * ny, 2height) + 1:fld(halfrow * ny - 1, 2height) + 1)
    function value(column, halfrow)
        cell = (values[i, j] for
            i in fld((column - 1) * nx, width) + 1:fld(column * nx - 1, width) + 1,
            j in fld((halfrow - 1) * ny, 2height) + 1:fld(halfrow * ny - 1, 2height) + 1
            if phi[i, j] <= 0 && isfinite(values[i, j]))
        total, count = 0.0, 0
        for v in cell
            total += v
            count += 1
        end
        return count == 0 ? lo : total / count
    end
    color(v) = FIELD_COLORS[1 + round(Int, (length(FIELD_COLORS) - 1) *
        clamp((v - lo) / max(hi - lo, eps()), 0, 1))]
    number(v) = string(round(v, sigdigits=4))

    frame = IOBuffer()
    clearline = tty ? "\e[K" : ""
    isnothing(label) || println(frame, label, clearline)
    println(frame, '+', repeat("-", width), '+', clearline)
    for row in height:-1:1
        print(frame, '|')
        for column in 1:width
            upper = inside(column, 2row)
            lower = inside(column, 2row - 1)
            if tty && !isnothing(values) && (upper || lower)
                upper && print(frame, "\e[38;5;$(color(value(column, 2row)))m")
                lower && print(frame, "\e[48;5;$(color(value(column, 2row - 1)))m")
                print(frame, upper ? '▀' : ' ')
                print(frame, "\e[0m")
            else
                print(frame, upper ? (lower ? '█' : '▀') : (lower ? '▄' : ' '))
            end
        end
        print(frame, '|')
        if tty && !isnothing(values)
            barvalue = lo + (hi - lo) * (row - 0.5) / height
            text = row == height ? number(hi) : row == cld(height, 2) ? number((lo + hi) / 2) : row == 1 ? number(lo) : ""
            print(frame, "  \e[48;5;$(color(barvalue))m \e[0m ", text)
        end
        println(frame, clearline)
    end
    print(frame, '+', repeat("-", width), '+', clearline, tty ? "\e[J" : "")
    tty && print(io, "\e[H")
    write(io, take!(frame))
    flush(io)
    return nothing
end

end # module
