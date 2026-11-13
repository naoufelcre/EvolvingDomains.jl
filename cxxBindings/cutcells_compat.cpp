#define EVOLVINGDOMAINS_CUTCELLS_BUILD
#include "cutcells_compat.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cut_quadrilateral.h"
#include "cut_triangle.h"
#include "cut_tetrahedron.h"

struct edcut_result {
    std::array<cutcells::cell::CutCell<double>, 3> regions;
};

namespace {

void set_error(char *buffer, uint64_t capacity, const char *message) noexcept
{
    if (buffer == nullptr || capacity == 0) return;
    const size_t n = static_cast<size_t>(capacity - 1);
    std::strncpy(buffer, message, n);
    buffer[n] = '\0';
}

void clear_error(char *buffer, uint64_t capacity) noexcept
{
    if (buffer != nullptr && capacity > 0) buffer[0] = '\0';
}

const cutcells::cell::CutCell<double> *get_region(
    const edcut_result *result, int32_t region) noexcept
{
    if (result == nullptr || region < EDCUT_NEGATIVE || region > EDCUT_POSITIVE)
        return nullptr;
    return &result->regions[static_cast<size_t>(region)];
}

bool valid_output_count(size_t count) noexcept
{
    return count <= static_cast<size_t>(std::numeric_limits<int32_t>::max());
}

void append_region(cutcells::cell::CutCell<double> &destination,
                   const cutcells::cell::CutCell<double> &source)
{
    const int base = static_cast<int>(destination._vertex_coords.size() / 2);
    destination._vertex_coords.insert(destination._vertex_coords.end(),
                                      source._vertex_coords.begin(), source._vertex_coords.end());
    for (size_t cell = 0; cell < source._types.size(); ++cell) {
        const int begin = source._offset[cell];
        const int end = source._offset[cell + 1];
        for (int i = begin; i < end; ++i)
            destination._connectivity.push_back(base + source._connectivity[i]);
        destination._offset.push_back(static_cast<int>(destination._connectivity.size()));
        destination._types.push_back(source._types[cell]);
    }
}

// Upstream's sign tolerance collapses to zero when the first sample is zero.
void rotate_zero_first(const double *coordinates, const double *values,
                        std::array<double, 8> &rotated_coordinates,
                        std::array<double, 4> &rotated_values)
{
    size_t rotation = 0;
    while (rotation < 4 && values[rotation] == 0.0) ++rotation;
    if (rotation == 4) throw std::invalid_argument("all level-set values are zero");
    if (rotation == 0) return;
    for (size_t i = 0; i < 4; ++i) {
        const size_t source = (i + rotation) % 4;
        rotated_coordinates[2 * i] = coordinates[2 * source];
        rotated_coordinates[2 * i + 1] = coordinates[2 * source + 1];
        rotated_values[i] = values[source];
    }
}

// Upstream emits cyclic quads but its generic triangulator assumes Basix ordering.
void triangulate_cyclic_quads(cutcells::cell::CutCell<double> &region)
{
    std::vector<int> connectivity;
    std::vector<int> offsets{0};
    std::vector<cutcells::cell::type> types;
    for (size_t cell = 0; cell < region._types.size(); ++cell) {
        const int begin = region._offset[cell];
        const int end = region._offset[cell + 1];
        if (begin < 0 || end < begin || static_cast<size_t>(end) > region._connectivity.size())
            throw std::runtime_error("CutCells returned invalid quadrilateral topology");
        const int count = end - begin;
        if (region._types[cell] == cutcells::cell::type::quadrilateral) {
            if (count != 4)
                throw std::runtime_error("CutCells returned a non-four-vertex quadrilateral");
            const int a = region._connectivity[begin];
            const int b = region._connectivity[begin + 1];
            const int c = region._connectivity[begin + 2];
            const int d = region._connectivity[begin + 3];
            connectivity.insert(connectivity.end(), {a, b, c, a, c, d});
            types.insert(types.end(), 2, cutcells::cell::type::triangle);
            offsets.push_back(static_cast<int>(connectivity.size()) - 3);
            offsets.push_back(static_cast<int>(connectivity.size()));
        } else {
            connectivity.insert(connectivity.end(), region._connectivity.begin() + begin,
                                region._connectivity.begin() + end);
            types.push_back(region._types[cell]);
            offsets.push_back(static_cast<int>(connectivity.size()));
        }
    }
    region._connectivity = std::move(connectivity);
    region._offset = std::move(offsets);
    region._types = std::move(types);
}

// Copy upstream tetrahedra verbatim, rejecting any other cell type.
void append_tetrahedra(cutcells::cell::CutCell<double> &destination,
                       const cutcells::cell::CutCell<double> &source)
{
    for (size_t cell=0; cell<source._types.size(); ++cell) {
        const int begin=source._offset[cell], end=source._offset[cell+1];
        if(source._types[cell]!=cutcells::cell::type::tetrahedron || end-begin!=4)
            throw std::runtime_error("CutCells returned a non-tetrahedron in a hexahedron phase");
        const int base=static_cast<int>(destination._vertex_coords.size()/3);
        for(int i=begin;i<end;++i) {
            const int id=source._connectivity[i];
            if(id<0 || static_cast<size_t>(3*id+2)>=source._vertex_coords.size())
                throw std::runtime_error("CutCells returned invalid 3-D connectivity");
            destination._vertex_coords.insert(destination._vertex_coords.end(),
                source._vertex_coords.begin()+3*id, source._vertex_coords.begin()+3*id+3);
        }
        destination._connectivity.insert(destination._connectivity.end(),{base,base+1,base+2,base+3});
        destination._types.push_back(cutcells::cell::type::tetrahedron);
        destination._offset.push_back(static_cast<int>(destination._connectivity.size()));
    }
}

// Append one tetrahedron given its own coordinates (used for uncut tets).
void append_tet(cutcells::cell::CutCell<double> &destination,
                const std::array<double,12> &tet_coordinates)
{
    const int base=static_cast<int>(destination._vertex_coords.size()/3);
    destination._vertex_coords.insert(destination._vertex_coords.end(),
                                      tet_coordinates.begin(), tet_coordinates.end());
    destination._connectivity.insert(destination._connectivity.end(),{base,base+1,base+2,base+3});
    destination._types.push_back(cutcells::cell::type::tetrahedron);
    destination._offset.push_back(static_cast<int>(destination._connectivity.size()));
}

// Constant level-set gradient of a tetrahedron, zero when it is degenerate.
std::array<double,3> tet_gradient(const std::array<double,12> &tet_coordinates,
                                  const std::array<double,4> &tet_values)
{
    const auto cross=[](const std::array<double,3> &a,const std::array<double,3> &b) {
        return std::array<double,3>{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
    };
    const auto dot=[](const std::array<double,3> &a,const std::array<double,3> &b) {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    };
    const std::array<double,3> e1{tet_coordinates[3]-tet_coordinates[0],
                                  tet_coordinates[4]-tet_coordinates[1],
                                  tet_coordinates[5]-tet_coordinates[2]};
    const std::array<double,3> e2{tet_coordinates[6]-tet_coordinates[0],
                                  tet_coordinates[7]-tet_coordinates[1],
                                  tet_coordinates[8]-tet_coordinates[2]};
    const std::array<double,3> e3{tet_coordinates[9]-tet_coordinates[0],
                                  tet_coordinates[10]-tet_coordinates[1],
                                  tet_coordinates[11]-tet_coordinates[2]};
    const auto c23=cross(e2,e3), c31=cross(e3,e1), c12=cross(e1,e2);
    const double det=dot(e1,c23);
    if(det==0.0) return {0.0,0.0,0.0};
    const double r1=tet_values[1]-tet_values[0];
    const double r2=tet_values[2]-tet_values[0];
    const double r3=tet_values[3]-tet_values[0];
    return {(r1*c23[0]+r2*c31[0]+r3*c12[0])/det,
            (r1*c23[1]+r2*c31[1]+r3*c12[1])/det,
            (r1*c23[2]+r2*c31[2]+r3*c12[2])/det};
}

// CutCells v0.4.0's quadrilateral lookup overcounts alternating-sign cases.
void cut_ambiguous_quad(const double *coordinates, const double *values,
                        std::array<cutcells::cell::CutCell<double>, 3> &regions)
{
    for (auto &region : regions) {
        region._gdim = 2;
        region._tdim = 2;
        region._offset.push_back(0);
    }
    const double scale = std::max({std::abs(values[0]), std::abs(values[1]),
                                   std::abs(values[2]), std::abs(values[3])});
    const double decider = (values[0] / scale) * (values[2] / scale) -
                           (values[1] / scale) * (values[3] / scale);
    const std::array<std::array<int, 3>, 2> triangles = decider > 0.0
        ? std::array<std::array<int, 3>, 2>{{{{0, 1, 2}}, {{0, 2, 3}}}}
        : std::array<std::array<int, 3>, 2>{{{{0, 1, 3}}, {{1, 2, 3}}}};
    const std::vector<std::string> names = {"phi<0", "phi=0", "phi>0"};
    for (const auto &triangle : triangles) {
        std::array<double, 6> triangle_coordinates;
        std::array<double, 3> triangle_values;
        for (int i = 0; i < 3; ++i) {
            const int vertex = triangle[i];
            triangle_coordinates[2 * i] = coordinates[2 * vertex];
            triangle_coordinates[2 * i + 1] = coordinates[2 * vertex + 1];
            triangle_values[i] = values[vertex];
        }
        std::vector<cutcells::cell::CutCell<double>> pieces;
        cutcells::cell::triangle::cut<double>(
            std::span<const double>(triangle_coordinates), 2,
            std::span<const double>(triangle_values), names, pieces, true);
        if (pieces.size() != regions.size())
            throw std::runtime_error("CutCells returned an unexpected triangle region count");
        for (size_t i = 0; i < regions.size(); ++i)
            append_region(regions[i], pieces[i]);
    }
}

} // namespace

extern "C" int32_t edcut_quad_create(
    const double *coordinates, uint64_t coordinate_count,
    const double *values, uint64_t value_count,
    edcut_result **out, char *error, uint64_t error_capacity)
{
    clear_error(error, error_capacity);
    if (out == nullptr) {
        set_error(error, error_capacity, "out must not be null");
        return 1;
    }
    *out = nullptr;
    if (coordinates == nullptr || coordinate_count != 8 || values == nullptr || value_count != 4) {
        set_error(error, error_capacity, "a quadrilateral requires 8 coordinate scalars and 4 level-set values");
        return 1;
    }
    for (uint64_t i = 0; i < coordinate_count; ++i) {
        if (!std::isfinite(coordinates[i])) {
            set_error(error, error_capacity, "coordinates must be finite");
            return 1;
        }
    }
    for (uint64_t i = 0; i < value_count; ++i) {
        if (!std::isfinite(values[i])) {
            set_error(error, error_capacity, "level-set values must be finite");
            return 1;
        }
    }

    try {
        std::array<double, 8> rotated_coordinates{};
        std::array<double, 4> rotated_values{};
        rotate_zero_first(coordinates, values, rotated_coordinates, rotated_values);
        const double *cut_coordinates = values[0] == 0.0 ? rotated_coordinates.data() : coordinates;
        const double *cut_values = values[0] == 0.0 ? rotated_values.data() : values;
        auto result = std::make_unique<edcut_result>();
        const int sign_mask = (values[0] < 0.0 ? 1 : 0) |
                              (values[1] < 0.0 ? 2 : 0) |
                              (values[2] < 0.0 ? 4 : 0) |
                              (values[3] < 0.0 ? 8 : 0);
        if ((sign_mask == 5 || sign_mask == 10) &&
            std::all_of(cut_values, cut_values + 4, [](double value) { return value != 0.0; })) {
            cut_ambiguous_quad(cut_coordinates, cut_values, result->regions);
        } else {
            const std::array<std::string, 3> regions = {"phi<0", "phi=0", "phi>0"};
            std::vector<cutcells::cell::CutCell<double>> cut_parts;
            cutcells::cell::quadrilateral::cut<double>(
                std::span<const double>(cut_coordinates, 8), 2,
                std::span<const double>(cut_values, 4),
                std::vector<std::string>(regions.begin(), regions.end()),
                cut_parts, false);
            if (cut_parts.size() != result->regions.size()) {
                set_error(error, error_capacity, "CutCells returned an unexpected region count");
                return 2;
            }
            for (size_t i = 0; i < cut_parts.size(); ++i) {
                triangulate_cyclic_quads(cut_parts[i]);
                result->regions[i] = std::move(cut_parts[i]);
            }
        }
        for (const auto &region : result->regions) {
            if (region._vertex_coords.size() % 2 != 0 ||
                !valid_output_count(region._vertex_coords.size()) ||
                !valid_output_count(region._connectivity.size()) ||
                !valid_output_count(region._offset.size()) ||
                !valid_output_count(region._types.size())) {
                set_error(error, error_capacity, "CutCells returned output too large for the C ABI");
                return 2;
            }
        }
        *out = result.release();
        return 0;
    } catch (const std::exception &e) {
        set_error(error, error_capacity, e.what());
        return 2;
    } catch (...) {
        set_error(error, error_capacity, "unknown C++ exception in CutCells");
        return 2;
    }
}

extern "C" int32_t edcut_hex_create(
    const double *coordinates, uint64_t coordinate_count,
    const double *values, uint64_t value_count,
    edcut_result **out, char *error, uint64_t error_capacity)
{
    clear_error(error, error_capacity);
    if (out == nullptr) { set_error(error, error_capacity, "out must not be null"); return 1; }
    *out = nullptr;
    if (coordinates == nullptr || coordinate_count != 24 || values == nullptr || value_count != 8) {
        set_error(error, error_capacity, "a hexahedron requires 24 coordinate scalars and 8 level-set values"); return 1;
    }
    for (uint64_t i=0;i<coordinate_count;++i) if (!std::isfinite(coordinates[i])) {
        set_error(error,error_capacity,"coordinates must be finite"); return 1;
    }
    for (uint64_t i=0;i<value_count;++i) if (!std::isfinite(values[i])) {
        set_error(error,error_capacity,"level-set values must be finite"); return 1;
    }
    if (std::all_of(values,values+8,[](double x){return x==0.0;})) {
        set_error(error,error_capacity,"all level-set values are zero"); return 1;
    }
    try {
        auto result=std::make_unique<edcut_result>();
        for (auto &r:result->regions) { r._gdim=3; r._tdim=3; r._offset.push_back(0); }
        const bool has_negative=std::any_of(values,values+8,[](double x){return x<0;});
        const bool has_positive=std::any_of(values,values+8,[](double x){return x>0;});
        constexpr std::array<std::array<int,4>,6> tets={{{{0,1,2,6}},{{0,2,3,6}},{{0,3,7,6}},{{0,7,4,6}},{{0,4,5,6}},{{0,5,1,6}}}};
        auto full_region=[&](int region) {
            auto &dst=result->regions[static_cast<size_t>(region)];
            for (const auto &tet:tets) {
                const int base=static_cast<int>(dst._vertex_coords.size()/3);
                for(int id:tet) dst._vertex_coords.insert(dst._vertex_coords.end(),coordinates+3*id,coordinates+3*id+3);
                dst._connectivity.insert(dst._connectivity.end(),{base,base+1,base+2,base+3});
                dst._types.push_back(cutcells::cell::type::tetrahedron);
                dst._offset.push_back(static_cast<int>(dst._connectivity.size()));
            }
        };
        if (!(has_negative && has_positive)) {
            full_region(has_negative ? EDCUT_NEGATIVE : EDCUT_POSITIVE);
        } else {
            std::array<double,8> cut_values{};
            std::copy(values, values + 8, cut_values.begin());
            // Upstream needs strictly nonzero samples and counts phi==0 as phi<0,
            // so nudge exact zeros onto that side.
            double value_scale=0;
            for(double x:cut_values) value_scale=std::max(value_scale,std::abs(x));
            for(double &x:cut_values) if(x==0.0) x=-1e-12*value_scale;

            // Direction toward phi>0 used only when a tetrahedron gradient vanishes.
            std::array<double,3> fallback_direction{};
            {
                std::array<double,3> negative_mean{}, positive_mean{};
                int negative_count=0, positive_count=0;
                for(int i=0;i<8;++i) {
                    if(values[i]<0) { for(int d=0;d<3;++d) negative_mean[d]+=coordinates[3*i+d]; ++negative_count; }
                    else if(values[i]>0) { for(int d=0;d<3;++d) positive_mean[d]+=coordinates[3*i+d]; ++positive_count; }
                }
                if(negative_count>0 && positive_count>0)
                    for(int d=0;d<3;++d)
                        fallback_direction[d]=positive_mean[d]/positive_count-negative_mean[d]/negative_count;
            }

            // Cut the cube as six tetrahedra: unlike the generated hexahedron tables,
            // tetrahedron pieces are always convex and tile the cube exactly.
            const std::array<std::string,3> names={"phi<0","phi=0","phi>0"};
            for(const auto &tet:tets) {
                std::array<double,12> tet_coordinates{};
                std::array<double,4> tet_values{};
                for(int i=0;i<4;++i) {
                    for(int d=0;d<3;++d) tet_coordinates[3*i+d]=coordinates[3*tet[i]+d];
                    tet_values[i]=cut_values[tet[i]];
                }
                const bool tet_negative=std::any_of(tet_values.begin(),tet_values.end(),[](double x){return x<0;});
                const bool tet_positive=std::any_of(tet_values.begin(),tet_values.end(),[](double x){return x>0;});
                if(!(tet_negative && tet_positive)) {
                    append_tet(result->regions[tet_negative ? EDCUT_NEGATIVE : EDCUT_POSITIVE],tet_coordinates);
                    continue;
                }
                std::array<cutcells::cell::CutCell<double>,3> pieces;
                for(int region=0;region<3;++region)
                    cutcells::cell::tetrahedron::cut<double>(std::span<const double>(tet_coordinates),3,
                        std::span<const double>(tet_values),names[static_cast<size_t>(region)],pieces[static_cast<size_t>(region)],true);
                append_tetrahedra(result->regions[EDCUT_NEGATIVE],pieces[EDCUT_NEGATIVE]);
                append_tetrahedra(result->regions[EDCUT_POSITIVE],pieces[EDCUT_POSITIVE]);
                auto gradient=tet_gradient(tet_coordinates,tet_values);
                if(gradient[0]==0.0 && gradient[1]==0.0 && gradient[2]==0.0)
                    gradient=fallback_direction;
                auto &surface=result->regions[EDCUT_INTERFACE];
                const auto &part=pieces[EDCUT_INTERFACE];
                for(size_t cell=0;cell<part._types.size();++cell) {
                    const int begin=part._offset[cell], end=part._offset[cell+1];
                    if(part._types[cell]!=cutcells::cell::type::triangle || end-begin!=3)
                        throw std::runtime_error("CutCells returned a non-triangular interface cell");
                    std::array<std::array<double,3>,3> tri{};
                    for(int i=0;i<3;++i) {
                        const int id=part._connectivity[begin+i];
                        if(id<0 || static_cast<size_t>(3*id+2)>=part._vertex_coords.size())
                            throw std::runtime_error("CutCells returned invalid interface connectivity");
                        for(int d=0;d<3;++d) tri[i][d]=part._vertex_coords[3*id+d];
                    }
                    // Orient the triangle toward increasing phi, i.e. toward phi>0.
                    const std::array<double,3> n{(tri[1][1]-tri[0][1])*(tri[2][2]-tri[0][2])-(tri[1][2]-tri[0][2])*(tri[2][1]-tri[0][1]),
                        (tri[1][2]-tri[0][2])*(tri[2][0]-tri[0][0])-(tri[1][0]-tri[0][0])*(tri[2][2]-tri[0][2]),
                        (tri[1][0]-tri[0][0])*(tri[2][1]-tri[0][1])-(tri[1][1]-tri[0][1])*(tri[2][0]-tri[0][0])};
                    if(n[0]*gradient[0]+n[1]*gradient[1]+n[2]*gradient[2]<0) std::swap(tri[1],tri[2]);
                    const int base=static_cast<int>(surface._vertex_coords.size()/3);
                    for(const auto &v:tri) surface._vertex_coords.insert(surface._vertex_coords.end(),v.begin(),v.end());
                    surface._connectivity.insert(surface._connectivity.end(),{base,base+1,base+2});
                    surface._types.push_back(cutcells::cell::type::triangle);
                    surface._offset.push_back(static_cast<int>(surface._connectivity.size()));
                }
            }
        }
        for (const auto &region:result->regions) {
            if(region._vertex_coords.size()%3 || !valid_output_count(region._vertex_coords.size()) ||
               !valid_output_count(region._offset.size()) || !valid_output_count(region._connectivity.size()) ||
               !valid_output_count(region._types.size())) {
                set_error(error,error_capacity,"CutCells returned output too large for the C ABI"); return 2;
            }
        }
        *out=result.release(); return 0;
    } catch(const std::exception &e) { set_error(error,error_capacity,e.what()); return 2; }
      catch(...) { set_error(error,error_capacity,"unknown C++ exception in CutCells"); return 2; }
}

extern "C" int32_t edcut_region_sizes(
    const edcut_result *result, int32_t region,
    uint64_t *points_count, uint64_t *offsets_count,
    uint64_t *connectivity_count, uint64_t *types_count,
    char *error, uint64_t error_capacity)
{
    clear_error(error, error_capacity);
    const auto *cell = get_region(result, region);
    if (cell == nullptr) {
        set_error(error, error_capacity, "invalid result, region");
        return 1;
    }
    if (points_count == nullptr || offsets_count == nullptr ||
        connectivity_count == nullptr || types_count == nullptr) {
        set_error(error, error_capacity, "output counts must not be null");
        return 1;
    }
    *points_count = static_cast<uint64_t>(cell->_vertex_coords.size());
    *offsets_count = static_cast<uint64_t>(cell->_offset.size());
    *connectivity_count = static_cast<uint64_t>(cell->_connectivity.size());
    *types_count = static_cast<uint64_t>(cell->_types.size());
    return 0;
}

extern "C" int32_t edcut_region_copy(
    const edcut_result *result, int32_t region,
    double *points, uint64_t points_capacity,
    int32_t *offsets, uint64_t offsets_capacity,
    int32_t *connectivity, uint64_t connectivity_capacity,
    int32_t *types, uint64_t types_capacity,
    char *error, uint64_t error_capacity)
{
    clear_error(error, error_capacity);
    const auto *cell = get_region(result, region);
    if (cell == nullptr) {
        set_error(error, error_capacity, "invalid result or region");
        return 1;
    }
    if (points_capacity < cell->_vertex_coords.size() ||
        offsets_capacity < cell->_offset.size() ||
        connectivity_capacity < cell->_connectivity.size() ||
        types_capacity < cell->_types.size() ||
        (!cell->_vertex_coords.empty() && points == nullptr) ||
        (!cell->_offset.empty() && offsets == nullptr) ||
        (!cell->_connectivity.empty() && connectivity == nullptr) ||
        (!cell->_types.empty() && types == nullptr)) {
        set_error(error, error_capacity, "output capacity is too small or a required output pointer is null");
        return 1;
    }
    for (size_t i = 0; i < cell->_offset.size(); ++i) {
        if (!valid_output_count(static_cast<size_t>(cell->_offset[i]))) {
            set_error(error, error_capacity, "CutCells offset is outside the C ABI range");
            return 2;
        }
    }
    for (size_t i = 0; i < cell->_connectivity.size(); ++i) {
        if (!valid_output_count(static_cast<size_t>(cell->_connectivity[i]))) {
            set_error(error, error_capacity, "CutCells connectivity is outside the C ABI range");
            return 2;
        }
    }
    if (!cell->_vertex_coords.empty())
        std::copy(cell->_vertex_coords.begin(), cell->_vertex_coords.end(), points);
    if (!cell->_offset.empty())
        std::copy(cell->_offset.begin(), cell->_offset.end(), offsets);
    if (!cell->_connectivity.empty())
        std::transform(cell->_connectivity.begin(), cell->_connectivity.end(), connectivity,
                       [](int value) { return static_cast<int32_t>(value); });
    if (!cell->_types.empty())
        std::transform(cell->_types.begin(), cell->_types.end(), types,
            [](cutcells::cell::type value) {
                switch (value) {
                    case cutcells::cell::type::interval: return static_cast<int32_t>(EDCUT_INTERVAL);
                    case cutcells::cell::type::triangle: return static_cast<int32_t>(EDCUT_TRIANGLE);
                    case cutcells::cell::type::tetrahedron: return static_cast<int32_t>(EDCUT_TETRAHEDRON);
                    default: return static_cast<int32_t>(EDCUT_UNSUPPORTED);
                }
            });
    return 0;
}

extern "C" void edcut_result_free(edcut_result *result)
{
    delete result;
}

extern "C" const char *edcut_upstream_version(void)
{
    return "0.4.0";
}
