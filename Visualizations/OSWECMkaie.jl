using GeometryBasics
using FileIO
using MeshIO
using Makie
using GLMakie

# Load the following .stl files and mesh them together 
mesh1 = load("Models/geometryFiles/OSWEC/flap.stl")
vertices1 = GeometryBasics.coordinates(mesh1)
faces1 = GeometryBasics.faces(mesh1)


mesh2 = load("Models/geometryFiles/OSWEC/base.stl")
vertices2 = GeometryBasics.coordinates(mesh2)
faces2 = GeometryBasics.faces(mesh2)

# Calculate the offset to place the first mesh on top of the second mesh
max_height2 = maximum(v[3] for v in vertices2)
min_height1 = minimum(v[3] for v in vertices1)
vertical_offset = max_height2 - min_height1 + 1
vertices1 = [Point3f0(v[1], v[2], v[3] + vertical_offset) for v in vertices1] # -10 part was hardcoded to better match the visualization compared to MATLAB

# Combine vertices and faces
vertices_combined = vcat(vertices1, vertices2)
faces_combined = vcat(faces1, [NgonFace(f.data .+ length(vertices1)) for f in faces2])

# Function to plot the mesh
function plot_mesh(vertices, faces)
    set_theme!(theme_black())
    fig = Figure(resolution = (400, 400))
    ax = LScene(fig[1, 1], scenekw = (raw = true,))
    
    # Plot each face of the mesh
    for face in faces
        verts = [vertices[i] for i in face.data]
        mesh!(ax, verts, color = :brown)
    end
    
    display(fig)
    
    # Update camera position for better view
    cam = cameracontrols(ax.scene)
    update_cam!(cam, Vec3f0(120, 120, 120), Vec3f0(0, 0, 0), Vec3f0(0, 1, 0))
end

# Plot the combined mesh
plot_mesh(vertices_combined, faces_combined)