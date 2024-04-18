import os
import time
import numpy as np


class AramisExporter:
    """We use this class as a HELPER module to read and export the data of an OPEN GOM Aramis Professional project
       as a txt file in a specific format.
       You can only use this module IF YOU WORK ON YOUR GOM ARAMIS SYSTEM. And only, directly from the
       GOM Aramis scripting editor.

    Methods:
        * show_last_stage - shows last stage of project and returns the stage object
        * show_stage - shows specific stage and returns the aramis stage object
        * get_result_dict - gets all results (facet coordinates, displacements, strain) as np.arrays and save to a dictionary
        * gather_process_data - gathers all process data from the aramis project
        * gather_stage_process_data - gathers process data which may change by stage
        * expot_data - export all stages to txt files or specific stages to txt files
        * export_data_to_vtk - export all stages to vtk files or specific stages to vtk files

    """

    def __init__(self, gom_app, gom_script, project_name: str, specimen_name: str, experiment_name: str):
        """Initializes class arguments. This is only possible to initialize if called from the gom python instance.

        Args:
            gom_app: (gom.app object) read from the open GOM Aramis project
            gom_script: (gom.script object) read from the open GOM Aramis project
            project_name: name of associated project, might be the gom project or free given string
            specimen_name: name of the specimen analyzed, can be any string
            experiment_name: name of the experiment, can be any string
        """
        print("Reading data from the GOM Aramis Professional project in CrackPy...")

        self.gom_app = gom_app
        self.project = gom_app.project
        self.script = gom_script
        self.project_name = project_name
        self.specimen_name = specimen_name
        self.experiment_name = experiment_name

        # auto functions
        self.ref_stage = self._get_ref_stage()
        self.num_stages = self._get_num_stages()
        self.result_types, self.disp_directions = self._init_result_type_dict()
        self._check_surface_components()
        self.current_surface_comp = self._get_current_component()
        self._check_rbmc()
        self.process_data, self.calibration_data = self.gather_process_data()

    def _get_num_stages(self) -> int:
        """
        Returns:
            total number of stages
        """
        return len(self.project.stages)

    def _get_ref_stage(self) -> int:
        """
        Returns:
            reference stage index
        """
        return int(self.project.get('reference_stage.index') - 1)

    def _check_export_folder_exists(self, export_folder_name: str) -> os.path:
        """Check if the export folder exists. If not creates.
        Args:
            export_folder_name: name of export folder relative to project directory

        Returns:
            path object for export directory
        """
        if not os.path.isdir(os.path.join(os.path.dirname(self.project.project_file), export_folder_name)):
            os.mkdir(os.path.join(os.path.dirname(self.project.project_file), export_folder_name))
            print(f"Creating a folder '{export_folder_name}'.")
        print(f"Data will be saved in '{export_folder_name}'.")
        return os.path.join(os.path.dirname(self.project.project_file), export_folder_name)

    def _check_all_stages_active(self):
        """Checks if all stages are active.
            Raises Error if not."""
        for stage in self.project.stages:
            if not stage.is_active:
                raise ValueError("Exports only work if ALL stages are set active.")
            else:
                pass

    @staticmethod
    def _init_result_type_dict():
        """Initialization of a result type dictionary mapping dic result types to general result types.

        Returns:
            (dicts) result types, displacement directions

        """
        rtd = {".epsX": "epsilon_x",
               ".epsY": "epsilon_y",
               ".epsXY": "epsilon_xy",
               ".phiM": "mises_strain",
               ".dX": "displacement",
               ".dY": "displacement",
               ".dZ": "displacement"
               }
        disp_directions = {".dX": "x",
                           ".dY": "y",
                           ".dZ": "z"
                           }
        return rtd, disp_directions

    def _get_current_component(self, index=0):
        """Gets the surface component of index as current_component.

        Returns:
            index of current surface component

        """
        return self.project.actual_elements.filter("type", "surface_component")[index]

    def _get_current_stage_indx(self) -> int:
        """Gets the current stage index. Starting from 0!

        Returns:
            index of current stage

        """
        return int(self._get_current_component().get('stage.index') - 1)

    def _check_rbmc(self):
        """Checks if rigid body motion compensation was conducted and is active."""
        rbmc_applied = False
        gom_elements = self.project.alignments.filter('type', 'transformation_object')
        for gom_elem in gom_elements:
            if gom_elem.get('object_family') == 'alignment_rbmc':
                rbmc_applied = True

        if not rbmc_applied:
            _ = self.script.alignment.create_rbmc_by_component(
                alignment_stage_creation_policy='separate_alignment_for_each_stage',
                component=self.project.actual_elements[self.current_surface_comp.name],
                name_expression='SKBK über $creation_sequence_args[\'component\'] != Unknown ? creation_sequence_args[\'component\'].name : \'?\'$')

        gom_elements = self.project.alignments.filter('type', 'transformation_object')
        for gom_elem in gom_elements:
            if gom_elem.get('object_family') == 'alignment_rbmc':
                rbmc_object_name = gom_elem.get('name')

        self.script.manage_alignment.set_alignment_active(
            movement_correction=self.project.alignments[rbmc_object_name])
        self.script.sys.recalculate_project(with_reports=False)

    def _check_surface_components(self):
        """Checks if all necessary data are calculated."""
        gom_surface_component_elements = self.project.actual_elements.filter('type', 'surface_component')
        surface_component_element_name = gom_surface_component_elements[0].get('name')

        actual_surf_elements = []
        for gom_element in self.project.inspection:
            actual_surf_elements.append(gom_element.get('name'))

        for result_string in self.result_types.keys():
            if surface_component_element_name + result_string not in actual_surf_elements:
                if self.result_types[result_string] == "displacement":
                    distance_restriction = self.disp_directions[result_string]
                    _ = self.script.inspection.inspect_dimension(
                        elements=[self.project.actual_elements[surface_component_element_name]],
                        distance_restriction=distance_restriction,
                        nominal_value=0.0,
                        nominal_value_source='fixed_value',
                        type=self.result_types[result_string])
                else:
                    _ = self.script.inspection.inspect_dimension(
                        elements=[self.project.actual_elements[surface_component_element_name]],
                        nominal_value=0.0,
                        nominal_value_source='fixed_value',
                        type=self.result_types[result_string])
                print(f"Creating surface element '{surface_component_element_name + result_string}' "
                      f"against nominal value = 0.0.")
        print("Recalculating...")
        self.script.sys.recalculate_project(with_reports=False)
        print("...done.")

    def show_last_stage(self):
        """Shows last stage of project and returns the stage object.

        Returns:
            aramis stage object
        """
        stage = self.project.stages[-1]
        self.script.sys.show_stage(stage=stage)
        print(f"Showing stage with name {self.project.stages[-1].get('name')} "
              f"and index {self.project.stages[-1].get('index')}...")
        return stage

    def show_stage(self, index: int):
        """Shows specific stage and returns the aramis stage object.

        Args:
            index: index used to specify stage. Starting with 0 for first stage

        Returns:
            aramis stage object

        """
        stage = self.project.stages[index]
        self.script.sys.show_stage(stage=stage)
        print(
            f"Showing stage with name {self.project.stages[index].get('name')} "
            f"and index {self.project.stages[index].get('index')}...")
        return stage

    def get_result_dict(self, current_stage_index: int) -> dict:
        """Gets all results as np.arrays and save to a dictionary.
        Results will be saved for facet coordinates (x, y, z), displacements(x, y), strain (x, y, xy)

        Args:
            current_stage_index: index of current stage starting from 0

        Returns:
            dictionary of the structure {ResultType : np.array(Result)}

        """
        print(f"Getting result dictionary for stage with name {self.project.stages[current_stage_index].get('name')}"
              f" and index {self.project.stages[current_stage_index].get('index')}...")
        current_surface_comp = self._get_current_component(index=0)
        facet_coordinates = np.array(
            self.project.actual_elements[current_surface_comp.name].data.coordinate[
                current_stage_index])  # check for s.th. like 'data.ref_coordinate
        facet_coordinates = np.where(np.abs(facet_coordinates) > 1e-30, facet_coordinates, 0)
        disp_x = np.array(
            self.project.inspection[current_surface_comp.name + '.dX'].data.result_dimension.deviation[
                current_stage_index])[
            0].flatten()
        disp_y = np.array(
            self.project.inspection[current_surface_comp.name + '.dY'].data.result_dimension.deviation[
                current_stage_index])[
            0].flatten()
        disp_z = np.array(
            self.project.inspection[current_surface_comp.name + '.dZ'].data.result_dimension.deviation[
                current_stage_index])[
            0].flatten()
        eps_x = \
            np.array(self.project.inspection[current_surface_comp.name + '.epsX'].data.result_dimension.deviation[
                         current_stage_index])[0].flatten()
        eps_y = \
            np.array(self.project.inspection[current_surface_comp.name + '.epsY'].data.result_dimension.deviation[
                         current_stage_index])[0].flatten()
        eps_xy = \
            np.array(self.project.inspection[current_surface_comp.name + '.epsXY'].data.result_dimension.deviation[
                         current_stage_index])[0].flatten()
        eps_eqv = \
            np.array(self.project.inspection[current_surface_comp.name + '.phiM'].data.result_dimension.deviation[
                         current_stage_index])[0].flatten()

        x_undef = facet_coordinates[:, 0] - disp_x
        y_undef = facet_coordinates[:, 1] - disp_y
        z_undef = facet_coordinates[:, 2] - disp_z

        return {"facet_coordinates": facet_coordinates,
                "x_undef": x_undef,
                "y_undef": y_undef,
                "z_undef": z_undef,
                "disp_x": disp_x,
                "disp_y": disp_y,
                "disp_z": disp_z,
                "eps_x": eps_x,
                "eps_y": eps_y,
                "eps_xy": eps_xy,
                "eps_eqv": eps_eqv}

    def gather_process_data(self):
        """Gathers all process data from the open Aramis dic project.

        Returns:
             (dict, dict) process data, calibration data

        """
        process_data = {"project_name": self.project_name, "specimen": self.specimen_name,
                        "experiment_number": self.experiment_name, "gom_project_file": self.project.project_file,
                        "project_creation_time": self.project.project_creation_time,
                        "sensor_name": self.gom_app.get("sys_sensor_configuration.name"),
                        "camera_type": self.gom_app.get("sys_sensor_configuration.camera_type"),
                        "camera_focal_length":
                        self.current_surface_comp.deformation_measurement_information.calibration.camera_focal_length,
                        "measuring_distance": self.gom_app.get(
                            "sys_sensor_configuration.scan_measuring_volume.measuring_distance"),
                        "camera_angle": self.gom_app.sys_calibration_camera_angle,
                        "camera_angle_degrees": self.gom_app.sys_calibration_camera_angle * 180.0 / np.pi,

                        }
        calibration_data = {"calibration_date": self.gom_app.sys_calibration_date,
                            "calibration_object": self.gom_app.sys_calibration_object_name,
                            "calibration_volume_width": self.gom_app.sys_calibration_volume_width,
                            "calibration_volume_length": self.gom_app.sys_calibration_volume_length,
                            "calibration_volume_depth": self.gom_app.sys_calibration_volume_depth,
                            "calibration_deviation":
                                self.current_surface_comp.deformation_measurement_information.calibration.deviation
                            }

        return process_data, calibration_data

    def gather_stage_process_data(self):
        """Gathers process data which may change by stage.

        Returns:
            (dicts) process data of actual stage, rigid body motion data for actual stage

        """
        stage_process_data = {"facet_size": self.current_surface_comp.facet_size,
                              "facet_distance": self.current_surface_comp.point_distance,
                              "exposure_time":
                                self.project.measurement_series['Deformation 1'].measurements['D1'].
                                    get('acquisition_parameters.exposure_time'),
                              "current_stage_index": self.current_surface_comp.get('stage.index') - 1,
                              "current_stage_name": self.current_surface_comp.get('stage.name'),
                              "current_stage_date": self.current_surface_comp.get('stage.absolute_time_stamp'),
                              "current_stage_relative_date": self.current_surface_comp.get('stage.relative_time'),
                              "reference_stage_index": self.current_surface_comp.get('reference_stage.index') - 1,
                              "reference_stage_name": self.current_surface_comp.get('reference_stage.display_name'),
                              "reference_stage_date": self.current_surface_comp.get(
                                  'reference_stage.absolute_time_stamp')
                              }

        gom_elements = self.project.alignments.filter('type', 'transformation_object')
        for gom_element in gom_elements:
            if gom_element.get('object_family') == 'alignment_rbmc':
                rbmc_object_name = gom_element.get('name')
        rmbc_data = {"alignment_rotation_x": self.project.alignments[rbmc_object_name].alignment.rotation.x,
                     "alignment_rotation_y": self.project.alignments[rbmc_object_name].alignment.rotation.y,
                     "alignment_rotation_z": self.project.alignments[rbmc_object_name].alignment.rotation.z,
                     "alignment_translation_x": self.project.alignments[rbmc_object_name].alignment.translation.x,
                     "alignment_translation_y": self.project.alignments[rbmc_object_name].alignment.translation.y,
                     "alignment_translation_z": self.project.alignments[rbmc_object_name].alignment.translation.z,
                     "alignment_deviation": self.project.alignments[rbmc_object_name].alignment.deviation,
                     }

        return stage_process_data, rmbc_data


    def export_data(self, stage_indxs: list or str = "all",
                    export_folder_name: str = "aramis_to_txt"):
        """Can be called to export all stages to txt files or specific stages to txt files.

        Args:
            stage_indxs: list of stage indexes or "all" or "last"
            export_folder_name: name of export folder. Will be created if not exists. Relative to aramis project path

        """
        self._check_all_stages_active()

        if stage_indxs == "all":
            stages = self.project.stages

        elif stage_indxs == "last":
            stages = [self.project.stages[-1]]
        else:
            stages = []
            for index in stage_indxs:
                for stage in self.project.stages:
                    if stage.get("index") == index:
                        stages.append(stage)
        export_directory = self._check_export_folder_exists(export_folder_name)
        print(f'Number of stages {len(stages)}')
        for current_stage in stages:
            current_stage_index = int(current_stage.get('index'))
            self.script.sys.show_stage(stage=current_stage)
            current_surface_comp = self._get_current_component(index=0).get('name')

            if self.project.inspection[current_surface_comp + '.dX'].computation_status == "computed":
                self._export_stage_to_txt(export_directory=export_directory,
                                          current_stage_index=current_stage_index)

    def export_data_to_vtk(self, stage_indxs: list or str = "all", export_folder_name: str = "aramis_to_vtk",
                           results: list or str = "all"):
        """
        Can be called to export all stages to vtk files or specific stages to vtk files.

        Args:
            stage_indxs: list of stage indexes or "all" or "last"
            export_folder_name: name of export folder. Will be created if not exists. Relative to aramis project path
            results: list of result types or 'all'

        """
        if results == "All" or results == "all":
            results = ["disp_x", "disp_y", "disp_z", "eps_x", "eps_y", "eps_xy", "eps_eqv"]

        self._check_all_stages_active()

        if stage_indxs == "all":
            stages = self.project.stages

        else:
            stages = []
            for index in stage_indxs:
                for stage in self.project.stages:
                    if stage.get("index") == index:
                        stages.append(stage)

        export_directory = self._check_export_folder_exists(export_folder_name)

        for current_stage in stages:
            current_stage_index = int(current_stage.get('index'))

            self.script.sys.show_stage(stage=current_stage)
            self._export_stage_to_vtk(export_directory=export_directory,
                                      results=results,
                                      current_stage_index=current_stage_index)

    def _export_stage_to_txt(self, export_directory: os.path, current_stage_index: int):
        """Exports exactly one stage to txt.

        Args:
            export_directory: relative path to export directory
            current_stage_index: index of stage to be exported

        """
        start = time.time()
        print(f"Exporting stage with index {current_stage_index} in...")
        # getting the current stage as object

        out_file_name = f"{self.project_name}_{self.specimen_name}_{self.experiment_name}" \
                        f"_{self.project.project_name}_dic_results_{self.ref_stage}_{current_stage_index}"

        out_file = open(os.path.join(export_directory, out_file_name + '.txt'), 'w')
        self._write_header(out_file)
        self._write_data(out_file, current_stage=current_stage_index)
        out_file.close()

        connection_file_name = out_file_name + "_connections"
        connection_file = open(os.path.join(export_directory, connection_file_name + '.txt'), 'w')
        self._write_connections(connection_file, current_stage=current_stage_index)

        connection_file.close()

        time_taken = time.time() - start
        print(f"{time_taken} seconds.")

    def _write_header(self, out_file):
        """Adds a header of metadata to an output file.

        Args:
            out_file: (open file object) Output file which should get the header

        """
        out_file.write("# Process data:\n")
        for each_key in self.process_data.keys():
            out_file.write(f"# {each_key.ljust(30)}: {str(self.process_data[each_key])}\n")
        out_file.write("#" * 100 + "\n")
        out_file.write("# Calibration data:\n")
        for each_key in self.calibration_data.keys():
            out_file.write(f"# {each_key.ljust(30)}: {str(self.calibration_data[each_key])}\n")
        out_file.write("#" * 100 + "\n")
        out_file.write("# stage data:\n")
        stage_data, rmbc_data = self.gather_stage_process_data()
        for each_key in stage_data:
            out_file.write(f"# {each_key.ljust(30)}: {str(stage_data[each_key])}\n")
        out_file.write("#" * 100 + "\n")

        out_file.write("# rigid body motion compensation:\n")
        for each_key in rmbc_data.keys():
            out_file.write(f"# {each_key.ljust(30)}: {str(rmbc_data[each_key])}\n")
        out_file.write("#" * 100 + "\n")

        # get input signals
        out_file.write("# SIGNALS:\n")
        gom_value_elements = self.project.inspection.filter('type', 'inspection_value_element')
        for gom_value_element in gom_value_elements:
            out_file.write(
                f"# {str(gom_value_element.get('name')).ljust(30)}: "
                f"{str(gom_value_element.get('type')).ljust(20)}: "
                f"{str(gom_value_element.get('input_value'))}\n")
        out_file.write('#\n')
        gom_value_elements = self.project.actual_elements.filter('type', 'value_element')
        for gom_value_element in gom_value_elements:
            out_file.write(
                f"# {str(gom_value_element.get('name')).ljust(30)}: "
                f"{str(gom_value_element.get('type')).ljust(20)}: "
                f"{str(gom_value_element.get('input_value'))}\n")
        out_file.write('#\n')
        gom_value_elements = self.project.actual_elements.filter('type', 'analog_input')
        for gom_value_element in gom_value_elements:
            out_file.write(
                f"# {str(gom_value_element.get('name')).ljust(30)}: "
                f"{str(gom_value_element.get('type')).ljust(20)}: "
                f"{str(gom_value_element.get('dimension'))}\n")
        out_file.write("#" * 100 + "\n")

    def _write_data(self, out_file, current_stage):
        """Adds data to an output file.

        Args:
            out_file: (open file object) Output file which should get the header

        """

        result_dic = self.get_result_dict(current_stage)
        out_file.write(
            '#{:>10}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}; {:>20}\n'.format(
                'ID', 'x_undef [mm]', 'y_undef [mm]', 'z_undef [mm]', 'u [mm]', 'v [mm]', 'w [mm]', 'epsx [%]',
                'epsy [%]', 'epsxy [1]', 'epseqv [%]'))

        for facet_index in range(len(result_dic["facet_coordinates"][:, 0])):
            out_file.write(
                '{:10.0f}; {:20.10f}; {:20.10f}; {:20.10f}; '
                '{:20.15f}; {:20.15f}; {:20.15f}; {:20.15f}; {:20.15f}; {:20.15f}; {:20.15f}\n'.format(
                    facet_index + 1,
                    result_dic["x_undef"][facet_index],
                    result_dic["y_undef"][facet_index],
                    result_dic["z_undef"][facet_index],
                    result_dic["disp_x"][facet_index],
                    result_dic["disp_y"][facet_index],
                    result_dic["disp_z"][facet_index],
                    result_dic["eps_x"][facet_index],
                    result_dic["eps_y"][facet_index],
                    result_dic["eps_xy"][facet_index],
                    result_dic["eps_eqv"][facet_index])
            )

    def _write_connections(self, connection_file, current_stage: int):
        """Internal routine to write a connection file, i.e. triangular connection of facets' center points.

        Args:
            connection_file: (open file object) Output file which should get the header
            current_stage: number of current stage

        """
        connection_file.write(
            f'{"Type":>10}; {"Element #":>10}; {"Node 1":>10}; {"Node 2":>10}; {"Node 3":>10}\n')

        connection_array = np.array(
            self.project.actual_elements[self.current_surface_comp.name].data.triangle[current_stage])
        for elem in range(len(connection_array[:, 0])):
            connection_file.write(
                f' {3:>10}; {elem + 1:>10}; {connection_array[elem, 0]:>10}; '
                f'{connection_array[elem, 1]:>10}; {connection_array[elem, 2]:>10}\n'
                )

    def _export_stage_to_vtk(self, export_directory: os.path, results: list or str, current_stage_index: int):
        """Exports exactly one stage to vtk.

        Args:
            export_directory: relative path to export directory
            results: list of result types or 'all'
            current_stage_index: index of stage to be exported

        """
        start = time.time()
        print(f"Exporting to vtk stage with index {current_stage_index} in...")
        # getting the current stage as object
        current_stage = self.project.stages[current_stage_index]
        self.script.sys.show_stage(stage=current_stage)
        out_file_name = f"{self.project_name}_{self.specimen_name}_{self.experiment_name}" \
                        f"_{self.project.project_name}_{self.ref_stage}_{current_stage_index}"

        with open(os.path.join(export_directory, out_file_name + '.vtk'), 'w') as out_file:
            out_file.write("# vtk DataFile Version 2.0\n"
                           "3D unstructured mesh of FE model with tetra elements\n"
                           "ASCII\n"
                           "\n"
                           "DATASET UNSTRUCTURED_GRID\n")
            self._write_data_to_vtk(out_file, results=results, current_stage=current_stage_index)

        time_taken = time.time() - start
        print(f"{time_taken} seconds.")

    def _write_data_to_vtk(self, out_file, results: list, current_stage: int):
        """Internal routine to fill an open .vtk file.

        Args:
            out_file: (open file object) Output .vtk file
            results: list of result types
            current_stage: number of current stage

        """
        result_dict = self.get_result_dict(current_stage)
        triangle_connections = np.array(
            self.project.actual_elements[self.current_surface_comp.name].data.triangle[current_stage])

        out_file.write(f'POINTS {len(result_dict["x_undef"][:])} float\n')
        for point_index in range(len(result_dict["x_undef"][:])):
            if np.isnan(result_dict["x_undef"][point_index]):
                out_file.write('0 0 0\n')
            else:
                out_file.write(
                    f'{result_dict["x_undef"][point_index]} '
                    f'{result_dict["y_undef"][point_index]} '
                    f'{result_dict["z_undef"][point_index]}\n')

        # get missing triangles
        no_of_missing_cells = 0
        for element_index in range(len(triangle_connections[:, 0])):
            if triangle_connections[element_index, 0] == -1:
                no_of_missing_cells += 1
            else:
                pass

        out_file.write('\n')
        out_file.write(
            f'CELLS {len(triangle_connections[:, 0]) - no_of_missing_cells} '
            f'{4 * (len(triangle_connections[:, 0]) - no_of_missing_cells)}\n')

        for elem in range(len(triangle_connections[:, 0])):
            if triangle_connections[elem, 0] == -1:
                pass
            else:
                out_file.write(
                    f'3 {triangle_connections[elem, 0]} '
                    f'{triangle_connections[elem, 1]} '
                    f'{triangle_connections[elem, 2]}\n')

        out_file.write('\n')
        out_file.write(f'CELL_TYPES {len(triangle_connections[:, 0]) - no_of_missing_cells}\n')
        for _ in range(len(triangle_connections[:, 0]) - no_of_missing_cells):
            out_file.write('5\n')
        del _
        # point data
        out_file.write(f'POINT_DATA {len(result_dict["disp_x"])}\n')

        for key in results:
            out_file.write(f'SCALARS {key} float\n')
            out_file.write('LOOKUP_TABLE default\n')
            for point_index in range(len(result_dict[key][:])):
                if np.isnan(result_dict[key][point_index]):
                    out_file.write('0.0\n')
                else:
                    out_file.write(f'{result_dict[key][point_index]}\n')

