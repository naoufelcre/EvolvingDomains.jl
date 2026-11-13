#ifndef EVOLVINGDOMAINS_CUTCELLS_COMPAT_H
#define EVOLVINGDOMAINS_CUTCELLS_COMPAT_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
  #if defined(EVOLVINGDOMAINS_CUTCELLS_BUILD)
    #define EDCUT_API __declspec(dllexport)
  #else
    #define EDCUT_API __declspec(dllimport)
  #endif
#else
  #define EDCUT_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct edcut_result edcut_result;

enum edcut_region {
    EDCUT_NEGATIVE = 0,
    EDCUT_INTERFACE = 1,
    EDCUT_POSITIVE = 2
};

enum edcut_cell_type {
    EDCUT_UNSUPPORTED = 0,
    EDCUT_INTERVAL = 1,
    EDCUT_TRIANGLE = 2,
    EDCUT_TETRAHEDRON = 3
};

/* Build a quadrilateral's phi<0, phi=0 and phi>0 pieces.
   Coordinates are four (x,y) pairs in VTK cyclic order: bottom-left,
   bottom-right, top-right, top-left;
   values has one signed level-set sample per vertex. All sizes are checked.
   On failure, returns a nonzero status and writes a bounded message to error. */
EDCUT_API int32_t edcut_quad_create(
    const double *coordinates, uint64_t coordinate_count,
    const double *values, uint64_t value_count,
    edcut_result **out, char *error, uint64_t error_capacity);

/* Build a hexahedron's phi<0, phi=0 and phi>0 pieces. Coordinates are eight
   (x,y,z) triples in VTK hexahedron order: (000,100,110,010,001,101,111,011).
   Values has one signed level-set sample per vertex. Volume pieces are
   tetrahedra and interface pieces are oriented triangles (normal toward phi>0).
   The cut is piecewise linear over the hexahedron's six-tetrahedron subdivision,
   so the phase pieces always tile the cell exactly; it is not an exact trilinear
   level-set reconstruction. Samples equal to zero are treated as phi<0.
   All sizes are checked; failure returns nonzero and writes a bounded message. */
EDCUT_API int32_t edcut_hex_create(
    const double *coordinates, uint64_t coordinate_count,
    const double *values, uint64_t value_count,
    edcut_result **out, char *error, uint64_t error_capacity);

/* Query flattened output sizes. points_count counts scalar coordinates;
   offsets_count is number of cells plus one; connectivity_count is the
   number of vertex indices; types_count is number of output cells. */
EDCUT_API int32_t edcut_region_sizes(
    const edcut_result *result, int32_t region,
    uint64_t *points_count, uint64_t *offsets_count,
    uint64_t *connectivity_count, uint64_t *types_count,
    char *error, uint64_t error_capacity);

/* Copy output arrays into caller-owned storage. All capacities are checked
   before any destination is written. Indices and offsets are zero-based. */
EDCUT_API int32_t edcut_region_copy(
    const edcut_result *result, int32_t region,
    double *points, uint64_t points_capacity,
    int32_t *offsets, uint64_t offsets_capacity,
    int32_t *connectivity, uint64_t connectivity_capacity,
    int32_t *types, uint64_t types_capacity,
    char *error, uint64_t error_capacity);

EDCUT_API void edcut_result_free(edcut_result *result);
EDCUT_API const char *edcut_upstream_version(void);

#ifdef __cplusplus
}
#endif

#endif
