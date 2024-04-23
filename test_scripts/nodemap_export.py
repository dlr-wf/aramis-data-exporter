# -*- coding: utf-8 -*-
"""
This script exports the Zeiss GOM Aramis data in a neutral data format. It can only be used with the Python interpreter
of Aramis Professional or GOM Inspect.

Output:

     - nodemap.txt: Facet data with coordinates, displacement vector and strain tensor.
     - Connections.txt: Connection information of the mesh.
     - .vtk: Export of the data in vtk format, which can be used directly in ParaView.
"""

import gom
from aramis_exporter.utils import AramisExporter

# Define the project, specimen, and wmp number, the aramis_project belongs to
PROJECT_NAME = "DELETE"
SPECIMEN_NAME = "SP"
WMP_NUMBER = "123"

# Definition of the name of the output file
export_file_name = f"{PROJECT_NAME}_{SPECIMEN_NAME}_{SPECIMEN_NAME}"

# Initialize gom_project with tom functionalities
exporter = AramisExporter(gom.app, gom.script,
                          export_file_name,
                          project_name=PROJECT_NAME,
                          specimen_name=SPECIMEN_NAME,
                          experiment_name=WMP_NUMBER,
                          )
# export data of given stages or all stages
exporter.export_data(stage_indxs='all')
exporter.export_data_to_vtk(stage_indxs='all')

print(f"Done.")
