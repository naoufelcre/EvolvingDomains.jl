module TerminalPlot

using Printf: @sprintf
using ..Geometric: EvolvingDiscreteGeometry, grid_info

export plot, plot_geometry, plot_curves

# ANSI approximations of Julia's purple, blue, green, and red logo colors.
const FIELD_COLORS = (54, 55, 56, 57, 63, 69, 75, 81, 80, 79, 78, 77, 76,
    112, 148, 184, 220, 214, 208, 202, 196)

# Curve palette sampled from the field ramp, so a field and its summary curves
# read as one visual system. Order: min, avg, max, then cycling.
const CURVE_COLORS = (FIELD_COLORS[4], FIELD_COLORS[11], FIELD_COLORS[19])
const CURVE_GLYPHS = ('*', '+', 'x', 'o', '#', '@')
const MIN_CURVE_WIDTH = 30
const MAX_CURVE_WIDTH = 60
const COLORBAR_WIDTH = 10

_istty(io::IO) = io isa Base.TTY || get(io, :terminal, false)

function _clean_label(label, columns::Int)
    isnothing(label) && return nothing
    return join(Iterators.take(first(split(string(label), '\n')), columns))
end

# Left-aligned title padded/truncated to the panel's visible width.
_fit(text, width::Int) = rpad(join(Iterators.take(string(text), width)), width)

_curve_title(curvetitle, ylabel) =
    isnothing(curvetitle) ? (isnothing(ylabel) ? "trace" : "$(ylabel)(t)") : string(curvetitle)

_number(v::Real) = lpad(@sprintf("%.2f", v), 6)

function _emit(io::IO, tty::Bool, frame::IOBuffer)
    tty && print(io, "\e[H")
    write(io, take!(frame))
    flush(io)
    return nothing
end

# --- Geometry panel -----------------------------------------------------------

function _geometry_rows(geom::EvolvingDiscreteGeometry, width::Int, height::Int,
    field, colorrange, tty::Bool; colorbar::Bool)
    nx, ny = grid_info(geom.grid).dims
    nx > 1 && ny > 1 || throw(ArgumentError("plot requires at least 2 nodes in each dimension"))
    phi = reshape(geom.levelset, nx, ny)
    if !isnothing(field)
        length(field) == length(phi) ||
            throw(DimensionMismatch("field has $(length(field)) values, expected $(length(phi))"))
    end
    values = isnothing(field) ? nothing : reshape(field, nx, ny)
    finite = isnothing(field) ? nothing : filter(isfinite, field)
    if !isnothing(finite) && isempty(finite)
        throw(ArgumentError("field has no finite values"))
    end
    lo, hi = isnothing(values) ? (0.0, 1.0) :
        (isnothing(colorrange) ? extrema(finite) : colorrange)
    color(v) = FIELD_COLORS[1 + round(Int, (length(FIELD_COLORS) - 1) *
        clamp((v - lo) / max(hi - lo, eps()), 0, 1))]

    inside(column, halfrow) = any(phi[i, j] <= 0
        for i in fld((column - 1) * nx, width) + 1:fld(column * nx - 1, width) + 1,
            j in fld((halfrow - 1) * ny, 2height) + 1:fld(halfrow * ny - 1, 2height) + 1)
    function cellvalue(column, halfrow)
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

    rows = String[]
    push!(rows, '+' * repeat("-", width) * '+')
    for row in height:-1:1
        buf = IOBuffer()
        print(buf, '|')
        for column in 1:width
            upper = inside(column, 2row)
            lower = inside(column, 2row - 1)
            if tty && !isnothing(values) && (upper || lower)
                upper && print(buf, "\e[38;5;$(color(cellvalue(column, 2row)))m")
                lower && print(buf, "\e[48;5;$(color(cellvalue(column, 2row - 1)))m")
                print(buf, upper ? '▀' : ' ')
                print(buf, "\e[0m")
            else
                print(buf, upper ? (lower ? '█' : '▀') : (lower ? '▄' : ' '))
            end
        end
        print(buf, '|')
        if colorbar && tty && !isnothing(values)
            # Gradient swatch on every row keeps all rows the same visual width,
            # so a side-by-side curve panel stays aligned; numbers only top/mid/bottom.
            barvalue = lo + (hi - lo) * (row - 0.5) / height
            # Visual row from the top, matching the curve panel's tick convention
            # (row counts from the bottom here); equal for odd and even heights.
            visual = height - row + 1
            text = visual == 1 ? _number(hi) : visual == cld(height, 2) ? _number((lo + hi) / 2) : visual == height ? _number(lo) : "      "
            print(buf, "  \e[48;5;$(color(barvalue))m \e[0m ", text)
        end
        push!(rows, String(take!(buf)))
    end
    push!(rows, '+' * repeat("-", width) * '+')
    return rows
end

"""
    plot_geometry(geom::EvolvingDiscreteGeometry; io=stdout, size=nothing, label=nothing,
                  field=nothing, colorrange=nothing)

Render the negative level-set region with Unicode block characters. Repeated calls in a
terminal replace the previous frame and adapt to the current terminal size. `label`
adds a line above the plot, for example `label="Iteration 100"`.
Pass a nodal scalar `field` to color the interior in ANSI terminals.
"""
function plot_geometry(geom::EvolvingDiscreteGeometry; io::IO=stdout,
    size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing,
    field::Union{Nothing,AbstractVector{<:Real}}=nothing, colorrange=nothing)
    rows, columns = isnothing(size) ? displaysize(io) : size
    tty = _istty(io)
    lab = _clean_label(label, columns)
    nx, ny = grid_info(geom.grid).dims
    nx > 1 && ny > 1 || throw(ArgumentError("plot requires at least 2 nodes in each dimension"))
    width = max(columns - 2 - (tty && !isnothing(field) ? COLORBAR_WIDTH : 0), 1)
    height = max(rows - 2 - !isnothing(lab), 1)

    # Most terminal cells are roughly twice as tall as they are wide.
    width = min(width, max(1, round(Int, 2 * height * (nx - 1) / (ny - 1))))
    height = min(height, max(1, round(Int, width * (ny - 1) / (2 * (nx - 1)))))
    grows = _geometry_rows(geom, width, height, field, colorrange, tty; colorbar=true)

    frame = IOBuffer()
    clearline = tty ? "\e[K" : ""
    isnothing(lab) || println(frame, lab, clearline)
    for r in grows[1:end-1]
        println(frame, r, clearline)
    end
    println(frame, grows[end], clearline)
    tty && print(frame, "\e[J")
    _emit(io, tty, frame)
    return nothing
end

# --- Curve panel --------------------------------------------------------------

function _normalize_series(t::AbstractVector, y::AbstractVector)
    if eltype(y) <: AbstractVector
        series = [collect(Float64, s) for s in y]
        all(s -> length(s) == length(t), series) ||
            throw(DimensionMismatch("series lengths $(map(length, series)) do not match t length $(length(t))"))
        return collect(Float64, t), series
    else
        length(y) == length(t) ||
            throw(DimensionMismatch("y length $(length(y)) does not match t length $(length(t))"))
        return collect(Float64, t), [collect(Float64, y)]
    end
end

function _normalize_series(t::AbstractVector, Y::AbstractMatrix)
    size(Y, 1) == length(t) ||
        throw(DimensionMismatch("Y has $(size(Y, 1)) rows, expected $(length(t))"))
    return collect(Float64, t), [collect(Float64, Y[:, k]) for k in axes(Y, 2)]
end

function _curve_limits(t, series, xrange, yrange)
    fx = filter(isfinite, t)
    fy = [v for s in series for v in s if isfinite(v)]
    x0, x1 = isnothing(xrange) ? (isempty(fx) ? (0.0, 1.0) : extrema(fx)) : xrange
    y0, y1 = isnothing(yrange) ? (isempty(fy) ? (0.0, 1.0) : extrema(fy)) : yrange
    x1 > x0 || (x0 -= 0.5; x1 += 0.5)
    y1 > y0 || (y0 -= 0.5; y1 += 0.5)
    return x0, x1, y0, y1
end

function _curve_rows(t, series, width::Int, height::Int, xrange, yrange, tty::Bool)
    x0, x1, y0, y1 = _curve_limits(t, series, xrange, yrange)
    cells = fill(' ', height, width)
    owners = zeros(Int, height, width)
    colof(x) = clamp(1 + round(Int, (x - x0) / (x1 - x0) * (width - 1)), 1, width)
    rowof(y) = clamp(1 + round(Int, (y1 - y) / (y1 - y0) * (height - 1)), 1, height)
    for (s, ys) in enumerate(series)
        glyph = CURVE_GLYPHS[mod1(s, length(CURVE_GLYPHS))]
        prev = nothing
        for (x, y) in zip(t, ys)
            if !(isfinite(x) && isfinite(y))
                prev = nothing
                continue
            end
            c, r = colof(x), rowof(y)
            cells[r, c] = glyph
            owners[r, c] = s
            if !isnothing(prev)
                pc, pr = prev
                for rr in min(r, pr):max(r, pr)
                    if cells[rr, c] == ' '
                        cells[rr, c] = '·'
                        owners[rr, c] = s
                    end
                end
            end
            prev = (c, r)
        end
    end
    rows = String[]
    push!(rows, '+' * repeat("-", width) * '+')
    for r in 1:height
        buf = IOBuffer()
        print(buf, '|')
        for c in 1:width
            if tty && owners[r, c] != 0
                print(buf, "\e[38;5;$(CURVE_COLORS[mod1(owners[r, c], length(CURVE_COLORS))])m",
                    cells[r, c], "\e[0m")
            else
                print(buf, cells[r, c])
            end
        end
        print(buf, '|')
        v = y1 - (r - 0.5) / height * (y1 - y0)
        text = r == 1 ? _number(y1) : r == cld(height, 2) ? _number((y0 + y1) / 2) : r == height ? _number(y0) : ""
        isempty(text) || print(buf, ' ', text)
        push!(rows, String(take!(buf)))
    end
    push!(rows, '+' * repeat("-", width) * '+')
    return rows
end

function _legend(t, series, xrange, yrange, ylabel, labels, tty::Bool)
    names = isnothing(labels) ? ["series $i" for i in eachindex(series)] : collect(labels)
    length(names) == length(series) ||
        throw(DimensionMismatch("got $(length(names)) labels for $(length(series)) series"))
    x0, x1, _, _ = _curve_limits(t, series, xrange, yrange)
    parts = String[]
    for (s, name) in enumerate(names)
        glyph = CURVE_GLYPHS[mod1(s, length(CURVE_GLYPHS))]
        if tty
            push!(parts, "\e[38;5;$(CURVE_COLORS[mod1(s, length(CURVE_COLORS))])m$glyph\e[0m $name")
        else
            push!(parts, "$glyph $name")
        end
    end
    head = isnothing(ylabel) ? "" : string(ylabel, ", ")
    return head * "t ∈ [$(_number(x0)), $(_number(x1))]: " * join(parts, "   ")
end

function _single_curves(io::IO, size, label, t, series, xrange, yrange, labels, ylabel, curvetitle)
    rows, columns = isnothing(size) ? displaysize(io) : size
    tty = _istty(io)
    lab = _clean_label(label, columns)
    legend = _legend(t, series, xrange, yrange, ylabel, labels, tty)
    width = min(max(columns - 2, 1), MAX_CURVE_WIDTH)
    height = max(rows - 2 - !isnothing(lab) - 2, 1)
    crows = _curve_rows(t, series, width, height, xrange, yrange, tty)
    frame = IOBuffer()
    clearline = tty ? "\e[K" : ""
    isnothing(lab) || println(frame, lab, clearline)
    println(frame, _fit(_curve_title(curvetitle, ylabel), width + 2), clearline)
    for r in crows
        println(frame, r, clearline)
    end
    println(frame, legend, clearline)
    tty && print(frame, "\e[J")
    _emit(io, tty, frame)
    return nothing
end

"""
    plot_curves(t, y; io=stdout, size=nothing, label=nothing, curvetitle=nothing,
                xrange=nothing, yrange=nothing, labels=nothing, ylabel=nothing)

GNUplot-style ASCII time series. `y` is one vector or a vector of vectors / matrix
columns for multi-series. Pass fixed `xrange`/`yrange` in loops to stop axis flashing.
"""
function plot_curves(t::AbstractVector, y::AbstractVector; io::IO=stdout,
    size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing, curvetitle=nothing,
    xrange=nothing, yrange=nothing, labels=nothing, ylabel=nothing)
    tt, series = _normalize_series(t, y)
    _single_curves(io, size, label, tt, series, xrange, yrange, labels, ylabel, curvetitle)
    return nothing
end

function plot_curves(t::AbstractVector, Y::AbstractMatrix; io::IO=stdout,
    size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing, curvetitle=nothing,
    xrange=nothing, yrange=nothing, labels=nothing, ylabel=nothing)
    tt, series = _normalize_series(t, Y)
    _single_curves(io, size, label, tt, series, xrange, yrange, labels, ylabel, curvetitle)
    return nothing
end

# --- Dual layout --------------------------------------------------------------

function _dual(io::IO, size, label, geom::EvolvingDiscreteGeometry, t, series,
    field, colorrange, xrange, yrange, labels, ylabel, geotitle, curvetitle)
    rows, columns = isnothing(size) ? displaysize(io) : size
    tty = _istty(io)
    lab = _clean_label(label, columns)
    legend = _legend(t, series, xrange, yrange, ylabel, labels, tty)
    nx, ny = grid_info(geom.grid).dims
    nx > 1 && ny > 1 || throw(ArgumentError("plot requires at least 2 nodes in each dimension"))
    showbar = tty && !isnothing(field)
    avail_h = max(rows - 2 - !isnothing(lab) - 2, 1)
    avail_w = max(columns - 4 - 2 - (showbar ? COLORBAR_WIDTH : 0), 1)
    wg_fit = max(1, round(Int, 2 * avail_h * (nx - 1) / (ny - 1)))
    wg = min(wg_fit, max(avail_w - MIN_CURVE_WIDTH, 1))
    wc = min(max(avail_w - wg, 1), MAX_CURVE_WIDTH)
    hg = min(avail_h, max(1, round(Int, wg * (ny - 1) / (2 * (nx - 1)))))
    grows = _geometry_rows(geom, wg, hg, field, colorrange, tty; colorbar=true)
    crows = _curve_rows(t, series, wc, hg, xrange, yrange, tty)
    # Legend starts under the curve panel, past geometry + colorbar + separator.
    indent = " " ^ (wg + 2 + (showbar ? COLORBAR_WIDTH : 0) + 2)
    frame = IOBuffer()
    clearline = tty ? "\e[K" : ""
    isnothing(lab) || println(frame, lab, clearline)
    println(frame, _fit(geotitle, wg + 2), "  ", _fit(_curve_title(curvetitle, ylabel), wc + 2), clearline)
    for i in eachindex(grows)
        println(frame, grows[i], "  ", crows[i], clearline)
    end
    println(frame, indent, legend, clearline)
    tty && print(frame, "\e[J")
    _emit(io, tty, frame)
    return nothing
end

"""
    plot(geom::EvolvingDiscreteGeometry; kwargs...)

Dispatch entry point: `plot_geometry` for a lone geometry, `plot_curves` for
`plot(t, y)`, dual geometry+curves panel for `plot(geom, t, y)`.
"""
plot(geom::EvolvingDiscreteGeometry; kwargs...) = plot_geometry(geom; kwargs...)
plot(t::AbstractVector, y::AbstractVector; kwargs...) = plot_curves(t, y; kwargs...)
plot(t::AbstractVector, Y::AbstractMatrix; kwargs...) = plot_curves(t, Y; kwargs...)

function plot(geom::EvolvingDiscreteGeometry, t::AbstractVector, y::AbstractVector;
    io::IO=stdout, size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing,
    field::Union{Nothing,AbstractVector{<:Real}}=nothing, colorrange=nothing,
    xrange=nothing, yrange=nothing, labels=nothing, ylabel=nothing,
    geotitle="φ ≤ 0", curvetitle=nothing)
    tt, series = _normalize_series(t, y)
    _dual(io, size, label, geom, tt, series, field, colorrange, xrange, yrange,
        labels, ylabel, geotitle, curvetitle)
    return nothing
end

function plot(geom::EvolvingDiscreteGeometry, t::AbstractVector, Y::AbstractMatrix;
    io::IO=stdout, size::Union{Nothing,Tuple{Int,Int}}=nothing, label=nothing,
    field::Union{Nothing,AbstractVector{<:Real}}=nothing, colorrange=nothing,
    xrange=nothing, yrange=nothing, labels=nothing, ylabel=nothing,
    geotitle="φ ≤ 0", curvetitle=nothing)
    tt, series = _normalize_series(t, Y)
    _dual(io, size, label, geom, tt, series, field, colorrange, xrange, yrange,
        labels, ylabel, geotitle, curvetitle)
    return nothing
end

end # module
