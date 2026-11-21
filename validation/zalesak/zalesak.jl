# The Zalesak problem: a slotted disk or sphere carried by a rigid rotation.
#
# Only the analytic problem lives here — initial level set, exact solution, velocity, exact
# measure — so that every variant script reads like the mathematical setup and nothing else.

# --- 2D slotted disk in the unit box; v = 20π(-(y-½), x-½); one turn per T = 0.1 ----

const OMEGA_2D = 20π
const PERIOD_2D = 2π/OMEGA_2D
const EXACT_AREA_2D = π*0.15^2 - 0.025*sqrt(0.15^2 - 0.025^2) - 0.15^2*asin(0.025/0.15)

# Signed distance to the slot: x ∈ [0.475, 0.525], y ∈ [0.60, 0.75].
function slot_2d(x, y)
    dx = max(0.475 - x, x - 0.525)
    dy = max(0.60 - y, y - 0.75)
    (dx <= 0 && dy <= 0) && return max(dx, dy)
    dx <= 0 && return dy
    dy <= 0 && return dx
    hypot(dx, dy)
end

"Disk of radius 0.15 centered at (0.5, 0.75) minus the slot."
disk(x, y) = max(hypot(x - 0.5, y - 0.75) - 0.15, -slot_2d(x, y))

"Exact level set at time `t`: the initial one rotated about the box center."
function disk_rotated(x, y, t)
    iszero(t) && return disk(x, y)
    s, c = sincos(OMEGA_2D*t)
    disk(0.5 + c*(x-0.5) + s*(y-0.5), 0.5 - s*(x-0.5) + c*(y-0.5))
end

disk_velocity(x, y) = SVector(-OMEGA_2D*(y - 0.5), OMEGA_2D*(x - 0.5))

# --- 3D slotted sphere in [-1,1]³; v = π(y,-x,0); one revolution per T = 2 ----

const CENTER_3D = (-1/3, 0.0, 0.0)
const RADIUS_3D = 0.5
const SLOT_3D = 1/6
const PERIOD_3D = 2.0
const EXACT_VOLUME_3D = 4π*RADIUS_3D^3/3 - π*SLOT_3D*(RADIUS_3D^2 - SLOT_3D^2/3)

"Sphere of radius 0.5 at (-1/3,0,0) minus a slot through all z."
function sphere(x, y, z)
    cx, cy, cz = CENTER_3D
    ball = hypot(x - cx, y - cy, z - cz) - RADIUS_3D
    slot = max(abs(x - cx) - SLOT_3D, abs(y - (cy + RADIUS_3D)) - RADIUS_3D)
    max(ball, -slot)
end

function sphere_rotated(x, y, z, t)
    iszero(t) && return sphere(x, y, z)
    s, c = sincos(π*t)
    sphere(c*x - s*y, s*x + c*y, z)
end

sphere_velocity(x, y, z) = SVector(π*y, -π*x, 0.0)

# --- 3D stress variant: small sphere at (0.5,0.75,0.5), tilted axis, ten turns ----

const OMEGA_TILT = 20π
const PERIOD_TILT = 0.1
const AXIS_TILT = (1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
const CENTER_TILT = (0.5, 0.75, 0.5)
const RADIUS_TILT = 0.15
const HALF_WIDTH_TILT = 0.025
const EXACT_VOLUME_TILT = 4π*RADIUS_TILT^3/3 -
    π*HALF_WIDTH_TILT*(RADIUS_TILT^2 - HALF_WIDTH_TILT^2/3)

function tilted_sphere(x, y, z)
    cx, cy, cz = CENTER_TILT
    ball = hypot(x - cx, y - cy, z - cz) - RADIUS_TILT
    slot = max(abs(x - cx) - HALF_WIDTH_TILT, abs(y - 0.675) - 0.075)
    max(ball, -slot)
end

"Exact level set at `t`: Rodrigues rotation of the initial data about the box center."
function tilted_rotated(x, y, z, t)
    iszero(t) && return tilted_sphere(x, y, z)
    s, c = sincos(OMEGA_TILT*t)
    ax, ay, az = AXIS_TILT
    rx, ry, rz = 0.5, 0.5, 0.5
    vx, vy, vz = x - rx, y - ry, z - rz
    wx = ay*vz - az*vy
    wy = az*vx - ax*vz
    wz = ax*vy - ay*vx
    d = ax*vx + ay*vy + az*vz
    tilted_sphere(rx + c*vx - s*wx + ax*d*(1-c),
                  ry + c*vy - s*wy + ay*d*(1-c),
                  rz + c*vz - s*wz + az*d*(1-c))
end

function tilted_velocity(x, y, z)
    ax, ay, az = AXIS_TILT
    SVector(OMEGA_TILT*(ay*(z-0.5) - az*(y-0.5)),
            OMEGA_TILT*(az*(x-0.5) - ax*(z-0.5)),
            OMEGA_TILT*(ax*(y-0.5) - ay*(x-0.5)))
end
