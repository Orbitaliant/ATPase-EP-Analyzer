import glob
import csv
from mmcif.io.IoAdapterCore import IoAdapterCore

io = IoAdapterCore()

cif_files = glob.glob("./1_input_pdb/*.cif")
cif_c_rings_dict = {}
for cif in cif_files:
  list_data_container = io.readFile(cif)
  data_container = list_data_container[0]
  for i in range(len(data_container.getObj('entity_poly'))):
    entity = data_container.getObj('entity')[i]
    entity_poly = data_container.getObj('entity_poly')[i]
    if int(entity[5]) >= 8:
      c_ring_helices_ids = entity_poly[6]
  cif_c_rings_dict[cif[-8:-4]] = c_ring_helices_ids

print(cif_c_rings_dict)

with open('pdb_data.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['pdb_id', 'c_ring_helices_ids'])
    for key, value in cif_c_rings_dict.items():
        writer.writerow([key, value])

print('pdb_data.csv created successfully.')