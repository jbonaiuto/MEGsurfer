import os
import sys
import subprocess
import numpy as np
import nibabel as nib
import trimesh
from scipy.spatial import KDTree

from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
from vtkmodules.vtkFiltersCore import vtkDecimatePro
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def create_surf_gifti(vertices, faces, normals=None):
    """
    Create a Gifti image object from surface mesh data.

    This function creates a GiftiImage object from the provided vertices, faces, and optional normals.
    The vertices and faces are required, while normals are optional. If normals are provided, they are
    added to the Gifti image. The function returns the GiftiImage object.

    Parameters:
    vertices (numpy.ndarray): Array of vertices. Each row represents a vertex with its x, y, z coordinates.
    faces (numpy.ndarray): Array of faces. Each row represents a face with three integers corresponding to vertex indices.
    normals (numpy.ndarray, optional): Array of vertex normals. Each row represents a normal vector corresponding to a vertex.

    Returns:
    nibabel.gifti.GiftiImage: The GiftiImage object created from the provided mesh data.

    Notes:
    - Vertex, face, and normal arrays should be NumPy arrays.
    - Vertices and normals should be in float32 format, and faces should be in int32 format.

    Example:
    >>> import numpy as np
    >>> import nibabel as nib
    >>> vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    >>> faces = np.array([[0, 1, 2], [0, 2, 3]])
    >>> normals = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]])
    >>> gifti_img = create_surf_gifti(vertices, faces, normals)
    """

    # Create new gifti object
    new_gifti = nib.gifti.GiftiImage()

    # Cast vertices and faces to the appropriate data types
    vertices = vertices.astype(np.float32)
    faces = faces.astype(np.int32)

    # Add the vertices and faces to the gifti object
    new_gifti.add_gifti_data_array(
        nib.gifti.GiftiDataArray(data=vertices, intent=nib.nifti1.intent_codes['NIFTI_INTENT_POINTSET']))
    new_gifti.add_gifti_data_array(
        nib.gifti.GiftiDataArray(data=faces, intent=nib.nifti1.intent_codes['NIFTI_INTENT_TRIANGLE']))

    # If normals are provided and not empty, cast them to float32 and add them to the Gifti image
    if normals is not None and len(normals) > 0:
        normals = normals.astype(np.float32)
        new_gifti.add_gifti_data_array(
            nib.gifti.GiftiDataArray(data=normals, intent=nib.nifti1.intent_codes['NIFTI_INTENT_VECTOR']))

    return new_gifti


def remove_vertices(gifti_surf, vertices_to_remove):
    """
    Remove specified vertices from a Gifti surface and update the faces accordingly.

    This function modifies a Gifti surface by removing the specified vertices. It also updates
    the faces of the surface so that they only reference the remaining vertices. If normals
    are present in the surface, they are also updated to correspond to the new set of vertices.

    Parameters:
    gifti_surf (nibabel.gifti.GiftiImage): The Gifti surface object from which vertices will be removed.
    vertices_to_remove (array_like): An array of vertex indices to be removed from the surface.

    Returns:
    nibabel.gifti.GiftiImage: A new GiftiImage object with the specified vertices removed and faces updated.

    Notes:
    - The function assumes that the GiftiImage object contains at least two data arrays: one for vertices
      and one for faces. If normals are present, they are also updated.
    - Vertex indices in `vertices_to_remove` should be zero-based (following Python's indexing convention).
    - The returned GiftiImage object is a new object; the original `gifti_surf` object is not modified in place.

    Example:
    >>> import nibabel as nib
    >>> gifti_surf = nib.load('path_to_gifti_file.gii')
    >>> vertices_to_remove = np.array([0, 2, 5])  # Indices of vertices to remove
    >>> new_gifti_surf = remove_vertices(gifti_surf, vertices_to_remove)
    """

    # Extract vertices and faces from the gifti object
    vertices_data = [da for da in gifti_surf.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_POINTSET']][0]
    faces_data = [da for da in gifti_surf.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_TRIANGLE']][0]

    vertices = vertices_data.data
    faces = faces_data.data

    # Determine vertices to keep
    vertices_to_keep = np.setdiff1d(np.arange(vertices.shape[0]), vertices_to_remove)

    # Create new array of vertices
    new_vertices = vertices[vertices_to_keep, :]

    # Find which faces to keep - ones that point to kept vertices
    face_x = np.isin(faces[:, 0], vertices_to_keep)
    face_y = np.isin(faces[:, 1], vertices_to_keep)
    face_z = np.isin(faces[:, 2], vertices_to_keep)
    faces_to_keep = np.where(face_x & face_y & face_z)[0]

    # Re-index faces
    x_faces = faces[faces_to_keep, :].reshape(-1)
    idxs = np.searchsorted(vertices_to_keep, x_faces)
    new_faces = idxs.reshape(-1, 3)

    # Create new gifti object
    normals=None
    normals_data = [da for da in gifti_surf.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_VECTOR']]
    if normals_data:
        normals = normals_data[0].data[vertices_to_keep, :]
    new_gifti = create_surf_gifti(new_vertices, new_faces, normals=normals)

    return new_gifti



def downsample_mesh_vtk(vertices, faces, reduction_ratio=0.1):
    """
    Downsample a mesh using the VTK library.

    This function takes a mesh defined by its vertices and faces, and downsamples it using
    VTK's vtkDecimatePro algorithm. The reduction ratio determines the degree of downsampling.
    The function returns the vertices and faces of the downsampled mesh.

    Parameters:
    vertices (numpy.ndarray): Array of vertices. Each row represents a vertex with its x, y, z coordinates.
    faces (numpy.ndarray): Array of faces. Each row represents a face with indices of the vertices forming the face.
    reduction_ratio (float): The proportion of the mesh to remove. For example, a reduction ratio of 0.1
                             retains 90% of the original mesh.

    Returns:
    tuple: A tuple containing two numpy arrays:
           - reduced_vertices: The vertices of the downsampled mesh.
           - reduced_faces: The faces of the downsampled mesh.

    Notes:
    - The input faces array should be triangulated, i.e., each face should consist of exactly three vertex indices.
    - The VTK library is used for mesh decimation, which must be installed and properly configured.

    Example:
    >>> import numpy as np
    >>> vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    >>> faces = np.array([[0, 1, 2], [0, 2, 3]])
    >>> reduced_vertices, reduced_faces = downsample_mesh_vtk(vertices, faces, 0.1)
    """

    # Convert vertices and faces to a VTK PolyData object
    points = vtkPoints()
    for point in vertices:
        points.InsertNextPoint(point)

    cells = vtkCellArray()
    for face in faces:
        cells.InsertNextCell(len(face))
        for vertex in face:
            cells.InsertCellPoint(vertex)

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    # Apply vtkDecimatePro for decimation
    decimate = vtkDecimatePro()
    decimate.SetInputData(polydata)
    decimate.SetTargetReduction(1 - reduction_ratio)
    decimate.Update()

    # Extract the decimated mesh
    decimated_polydata = decimate.GetOutput()

    # Convert back to numpy arrays
    reduced_vertices = vtk_to_numpy(decimated_polydata.GetPoints().GetData())

    # Extract and reshape the face data
    face_data = vtk_to_numpy(decimated_polydata.GetPolys().GetData())
    # Assuming the mesh is triangulated, every fourth item is the size (3), followed by three vertex indices
    reduced_faces = face_data.reshape(-1, 4)[:, 1:4]

    return reduced_vertices, reduced_faces


def downsample_multiple_surfaces(in_surfs, ratio):
    """
    Downampled multiple surface meshes using the VTK decimation algorithm.

    This function takes a list of input surface meshes (in Gifti format) and applies a dowsampling
    process to each surface. The downsampling is performed using VTK's vtkDecimatePro algorithm. The
    first surface in the list is downsampled, and its vertex mapping is then applied to all other
    surfaces in the list. The function returns a list of downsampled surface meshes.

    Parameters:
    in_surfs (list of nibabel.gifti.GiftiImage): Input Gifti surface meshes to be downsampled.
    ratio (float): The reduction ratio for the downsampling process. For example, a ratio of 0.1
                   implies that the mesh will be reduced to 90% of its original size.

    Returns:
    list of nibabel.gifti.GiftiImage: List of downsampled Gifti surface meshes.

    Notes:
    - The function prints the percentage of vertices retained in the first surface after downsampling.
    - If normals are present in the input surfaces, they are also downsampled and mapped to the new surfaces.
    - The resulting surfaces maintain the original topology and are suitable for visualization and further processing.

    Example:
    >>> import nibabel as nib
    >>> in_surfs = [nib.load('path/to/input_surf1.gii'), nib.load('path/to/input_surf2.gii')]
    >>> ratio = 0.1
    >>> out_surfs = downsample_multiple_surfaces(in_surfs, ratio)
    >>> for i, ds_surf in enumerate(out_surfs):
    ...     nib.save(ds_surf, f'path/to/output_surf{i+1}.gii')
    """
    out_surfs=[]

    primary_surf = in_surfs[0]
    reduced_vertices, reduced_faces = downsample_mesh_vtk(primary_surf.darrays[0].data,
                                                          primary_surf.darrays[1].data,
                                                          reduction_ratio=ratio)

    # Find the original vertices closest to the downsampled vertices
    kdtree = KDTree(primary_surf.darrays[0].data)
    # Calculate the percentage of vertices retained
    decim_orig_dist, orig_vert_idx = kdtree.query(reduced_vertices, k=1)
    print(f"{(1 - np.mean(decim_orig_dist > 0)) * 100}% of the vertices in the decimated first surface belong to the initial first surface vertices.")

    reduced_normals = None
    if len(primary_surf.darrays) > 2 and primary_surf.darrays[2].intent == nib.nifti1.intent_codes['NIFTI_INTENT_VECTOR']:
        reduced_normals = primary_surf.darrays[2].data[orig_vert_idx]

    # Save the downsampled primary surface with normals
    ds_surf = create_surf_gifti(reduced_vertices, reduced_faces, normals=reduced_normals)
    out_surfs.append(ds_surf)

    # Process other surfaces
    for i in range(1, len(in_surfs)):
        surf = in_surfs[i]

        reduced_normals = None
        if len(surf.darrays) > 2 and surf.darrays[2].intent == nib.nifti1.intent_codes['NIFTI_INTENT_VECTOR']:
            reduced_normals = surf.darrays[2].data[orig_vert_idx]

        ds_surf = create_surf_gifti(surf.darrays[0].data[orig_vert_idx,:], reduced_faces, normals=reduced_normals)
        out_surfs.append(ds_surf)
    return out_surfs


def combine_surfaces(surfaces):
    """
    Combine multiple surface meshes into a single surface mesh.

    This function takes a list of surface mesh files and combines them into a single surface mesh.
    It concatenates the vertices, faces, and normals (if present) from each surface. The faces are
    re-indexed appropriately to maintain the correct references to the combined vertex array.

    Parameters:
    surfaces (list of str): Paths to the surface mesh files (Gifti format) to be combined.
    out_fname (str): Path to the output file where the combined surface mesh will be saved.

    Notes:
    - The input surface mesh files should be in Gifti format.
    - The vertices, faces, and normals (if present) from each surface are concatenated.
    - The faces are re-indexed to reference the correct vertices in the combined vertex array.
    - The combined surface mesh is saved in Gifti format to the specified output file.

    Raises:
    ValueError: If the vertex or face arrays do not have the expected dimensions.

    Example:
    >>> surfaces = ['path/to/surface1.gii', 'path/to/surface2.gii']
    >>> out_fname = 'path/to/combined_surface.gii'
    >>> combine_surfaces(surfaces, out_fname)
    """

    combined_vertices = []
    combined_faces = []
    combined_normals = []

    face_offset = 0

    for mesh in surfaces:

        # Extract vertices and faces
        vertices = np.concatenate(
            [da.data for da in mesh.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_POINTSET']])
        faces = np.concatenate(
            [da.data for da in mesh.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_TRIANGLE']])

        # Check for normals
        normal_arrays = [da.data for da in mesh.darrays if da.intent == nib.nifti1.intent_codes['NIFTI_INTENT_VECTOR']]
        normals = np.concatenate(normal_arrays) if normal_arrays else np.array([])

        if vertices.ndim != 2 or vertices.shape[1] != 3:
            raise ValueError("Vertices array should have shape [n, 3]")
        if faces.ndim != 2 or faces.shape[1] != 3:
            raise ValueError("Faces array should have shape [n, 3]")

        combined_vertices.append(vertices)
        combined_faces.append(faces + face_offset)
        if normals.size:
            combined_normals.append(normals)

        face_offset += vertices.shape[0]

    # Combine the arrays
    combined_vertices = np.vstack(combined_vertices).astype(np.float32)
    combined_faces = np.vstack(combined_faces).astype(np.int32)
    if combined_normals:
        combined_normals = np.vstack(combined_normals).astype(np.float32)

    combined_surf = create_surf_gifti(combined_vertices, combined_faces, normals=combined_normals)
    return combined_surf
