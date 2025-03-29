# Load necessary PyMOL objects
loadall *_fixed.pqr
loadall *.dx
[cmd.set_name(obj, "apbs_pot") for obj in cmd.get_names("all") if obj.endswith("-pot")]

# Create an isosurface colored per APBS-calculated electrostatics 
isosurface pos_surf, apbs_pot, level=1
color red, pos_surf
isosurface neg_surf, apbs_pot, level=-1
color blue, neg_surf

# Adjust scene
bg grey90
turn x, 90
center apbs_pot

# Save PyMOL scene
save apbs_pot_isosurface.pse

# Set viewport
viewport 1920, 1080

# Set PyMOL memory usage
set hash_max, 400

# Take screenshot 1
png apbs_pot_isosurface_sideview.png

# Take screenshot 2
turn x, 90
png apbs_pot_isosurface_topview.png