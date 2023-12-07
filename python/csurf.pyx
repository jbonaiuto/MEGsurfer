# csurf.pyx
import numpy as np
cimport numpy as cnp
from scipy.sparse import coo_matrix

cdef extern from "float.h":
    double DBL_MAX

cdef void fill_adjacency_matrix(cnp.ndarray[int, ndim=2] faces, int num_vertices, cnp.ndarray[int, ndim=1] rows, cnp.ndarray[int, ndim=1] cols):
    cdef int idx = 0
    for i in range(faces.shape[0]):
        for j in range(3):
            rows[idx] = faces[i, j]
            cols[idx] = faces[i, (j + 1) % 3]
            idx += 1
            rows[idx] = faces[i, j]
            cols[idx] = faces[i, (j + 2) % 3]
            idx += 1

cdef void fill_distance_matrix(cnp.ndarray[double, ndim=2] vertices, cnp.ndarray[int, ndim=2] faces, cnp.ndarray[double, ndim=1] distances, cnp.ndarray[int, ndim=1] rows, cnp.ndarray[int, ndim=1] cols):
    cdef int idx = 0
    for i in range(faces.shape[0]):
        for j in range(3):
            v1 = vertices[faces[i, j]]
            v2 = vertices[faces[i, (j + 1) % 3]]
            rows[idx] = faces[i, j]
            cols[idx] = faces[i, (j + 1) % 3]
            distances[idx] = np.sqrt(np.sum((v1 - v2) ** 2))
            idx += 1

cdef void relax_edges(int u, cnp.ndarray[double, ndim=1] dist, cnp.ndarray[int, ndim=1] indices, cnp.ndarray[int, ndim=1] indptr, cnp.ndarray[double, ndim=1] data, int[:] remaining):
    cdef int v
    cdef double alt
    for v in indices[indptr[u]:indptr[u + 1]]:
        alt = dist[u] + data[indptr[u]:indptr[u + 1]][indices[indptr[u]:indptr[u + 1]] == v]
        if alt < dist[v]:
            dist[v] = alt

cdef cnp.ndarray[double, ndim=1] compute_geodesic_distances_internal(cnp.ndarray[double, ndim=2] vertices, cnp.ndarray[int, ndim=2] faces, cnp.ndarray[int, ndim=1] source_indices, double max_dist):
    cdef int num_vertices = vertices.shape[0]
    cdef cnp.ndarray[double, ndim=1] dist = np.full(num_vertices, np.inf)
    cdef cnp.ndarray[int, ndim=1] D_indices
    cdef cnp.ndarray[int, ndim=1] D_indptr
    cdef cnp.ndarray[double, ndim=1] D_data

    # Get data, indices, and indptr from the CSR matrix
    D_data = compute_mesh_distances(vertices, faces).data
    D_indices = compute_mesh_distances(vertices, faces).indices
    D_indptr = compute_mesh_distances(vertices, faces).indptr

    cdef int[:] remaining = np.arange(num_vertices, dtype=np.int32)
    cdef int remaining_size = remaining.size
    cdef int[:] new_remaining
    cdef int u, v, idx

    for src in source_indices:
        dist[src] = 0.0

    while remaining_size > 0:
        u = remaining[np.argmin(dist[remaining])]
        if dist[u] > max_dist:
            break
        relax_edges(u, dist, D_indices, D_indptr, D_data, remaining)

        new_remaining = np.empty(remaining_size - 1, dtype=np.int32)
        idx = 0
        for v in remaining:
            if v != u:
                new_remaining[idx] = v
                idx += 1
        remaining = new_remaining
        remaining_size -= 1

    return dist

def compute_mesh_adjacency(faces):
    num_vertices = np.max(faces) + 1
    rows = np.empty(6 * faces.shape[0], dtype=np.int32)
    cols = np.empty(6 * faces.shape[0], dtype=np.int32)
    fill_adjacency_matrix(faces, num_vertices, rows, cols)
    data = np.ones(len(rows), dtype=np.int32)
    return coo_matrix((data, (rows, cols)), shape=(num_vertices, num_vertices)).tocsr()

def compute_mesh_distances(vertices, faces):
    num_edges = 3 * faces.shape[0]
    distances = np.empty(num_edges, dtype=np.float64)
    rows = np.empty(num_edges, dtype=np.int32)
    cols = np.empty(num_edges, dtype=np.int32)
    fill_distance_matrix(vertices, faces, distances, rows, cols)
    return coo_matrix((distances, (rows, cols)), shape=(vertices.shape[0], vertices.shape[0])).tocsr()

def compute_geodesic_distances(vertices, faces, source_indices, max_dist=np.inf):
    return compute_geodesic_distances_internal(vertices.astype(np.float64), faces.astype(np.int32), np.array(source_indices, dtype=np.int32), max_dist)
