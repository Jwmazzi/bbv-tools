# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0613,W0511,R0903,E0401,R0902,R0913

"""
Name:         BBV.pyt
Author:       Jeffrey Scarmazzi 2023
Requirements: ArcGIS Pro 3.X 
Description:  This Python Toolbox takes in a portion of a Utility Network 
              and generates "segments" based on the connectivity of the network. There are two tools
              in this toolbox: SegmentBuilder and SegmentTrace.
"""

from collections import namedtuple
from pathlib import Path
from uuid import uuid4

import arcpy

# Generic representation of the objects traversed in SegmentBuilder
Node = namedtuple("Node", ["global_id", "geometry"])


class Toolbox:
    """ Entry Point for Python Toolbox """
    def __init__(self):
        # Define Name and Alias (Shown in Tool Execution) for the toolbox
        self.label = "BBV"
        self.alias = "Bounded by Valve"

        # List of tool classes associated with this toolbox
        self.tools = [SegmentBuilder, SegmentTrace]


class SegmentBuilder:
    """Generate a 'Segment Model' from a Utility Network."""

    def __init__(self):
        self.label = "Segment Builder"
        self.description = ""
        self.canRunInBackground = False

        # Keeping tracke of UN path for simplicity in function signatures
        self.un_dataset: Path = None

    def getParameterInfo(self):
        """Define parameter definitions"""

        param_0 = arcpy.Parameter(
            displayName="Input Utility Network",
            name="network",
            datatype="DEUtilityNetwork",
            parameterType="Required",
            direction="Input",
        )

        param_1 = arcpy.Parameter(
            displayName="Overwrite Existing Segment Tables",
            name="overwrite",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input",
        )
        param_1.value = True  # Default Value

        param_2 = arcpy.Parameter(
            displayName="Run Subnetwork",
            name="subnetwork",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            multiValue=False,
        )
        param_2.filter.type = "ValueList"
        param_2.filter.list = ["All Features"]

        return [param_0, param_1, param_2]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # If the Utility Network Parameter Has a Value
        if parameters[0].value:
            self.un_input = parameters[0].value

            # Get the Feature Dataset Path for UN and
            # Return Subnetwork Names from Water Line Feature Class
            un_path = str(Path(arcpy.Describe(parameters[0].value).catalogPath).parent)
            subnetworks = self.get_subnetwork_names(un_path)

            # Push Network List + All Features to the Subnetwork Parameter
            parameters[2].filter.list = sorted(subnetworks) + ["All Features"]

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    @staticmethod
    def get_subnetwork_names(un_dataset: str) -> list[str]:
        """Return unique subnetwork names from the WaterLine Feature Class

        Args:
            un_dataset (str): Path to the Utility Network Dataset

        Raises:
            ValueError: If the WaterLine Feature Class is not part of the Utility Network

        Returns:
            list[str]: Unique Subnetwork Names
        """

        arcpy.env.workspace = un_dataset

        targets = arcpy.ListFeatureClasses("WaterLine", feature_type="Polyline")

        if len(targets) != 1:
            raise ValueError(
                "Searching for WaterLine Was Ambiguous in target Utility Nework"
            )

        # Return Unique Subnetwork Names
        return list(
            set(
                [
                    row[0]
                    for row in arcpy.da.SearchCursor(
                        targets[0], ["PressureSubnetworkName"]
                    )
                ]
            )
        )

    def get_line_table(self, overwrite: bool) -> str:
        """Return the path to the Segment_Lines table and creates it if it does not exist.

        Args:
            overwrite (bool): Whether or not to truncate the existing table.

        Returns:
            str: Path to the Segment_Lines table.
        """

        # Unpack pathing information
        table_path = Path(self.un_dataset.parent)
        table_name = "Segment_Lines"
        full_path  = str(Path(table_path, table_name))

        # Create Empty Schema if it Does Not Exist
        if not arcpy.Exists(full_path):
            new_table = arcpy.management.CreateTable(str(table_path), table_name)
            arcpy.management.AddField(new_table, "segment_id", "TEXT")
            arcpy.management.AddField(new_table, "un_line_guid", "GUID")
            return new_table

        if overwrite:
            arcpy.AddMessage(f"Truncating Existing Table: {table_name}")
            arcpy.management.TruncateTable(full_path)

        return full_path

    def get_point_table(self, overwrite: bool) -> str:
        """Return the path to the Segment_Points table and creates it if it does not exist.

        Args:
            overwrite (bool): Whether or not to truncate the existing table.

        Returns:
            str: Path to the Segment_Points table.
        """

        # Unpack pathing information
        table_path = Path(self.un_dataset.parent)
        table_name = "Segment_Points"
        full_path  = str(Path(table_path, table_name))

        # Create Empty Schema if it Does Not Exist
        if not arcpy.Exists(full_path):
            new_table = arcpy.management.CreateTable(str(table_path), table_name)
            arcpy.management.AddField(new_table, "segment_id", "TEXT")
            arcpy.management.AddField(new_table, "un_point_guid", "GUID")
            return new_table

        if overwrite:
            arcpy.AddMessage(f"Truncating Existing Table: {table_name}")
            arcpy.management.TruncateTable(full_path)

        return full_path

    def stage_segment_tables(self) -> None:
        """Helper function to set class attributes for the Segment Tables."""

        arcpy.AddMessage("Staging Segment Tables")

        # Establish pathing the Segment Tables for later use
        self.point_table = self.get_point_table(self.overwrite)
        self.line_table  = self.get_line_table(self.overwrite)

    def get_lines(self, subnetwork: str = None) -> list[Node]:
        """Return lines used in segment building process.

        Args:
            subnetwork (str, optional): Name of subnetwork to trace. Defaults to None.

        Raises:
            ValueError: If the WaterLine Feature Class is not part of the Utility Network

        Returns:
            list[Node]: Lines in the network.
        """

        arcpy.env.workspace = str(self.un_dataset)

        targets = arcpy.ListFeatureClasses("WaterLine", feature_type="Polyline")

        if len(targets) != 1:
            raise ValueError(
                "Search for WaterLine Was Empty or Ambiguous in target Utility Nework"
            )

        # Only return list of Nodes below for the specified subnetwork
        if subnetwork and subnetwork != "All Features":
            query = f"PressureSubnetworkName = '{subnetwork}'"
        else:
            query = None

        return [
            Node(row[0], row[1])
            for row in arcpy.da.SearchCursor(targets[0], ["GLOBALID", "SHAPE@"], query)
        ]

    def get_points(self, subnetwork: str = None) -> list[Node]:
        """Return points used in segment building process.

        Args:
            subnetwork (str, optional): Name of the subnetwork to trace. Defaults to None.

        Raises:
            ValueError: If the WaterDevice Feature Class is not part of the Utility Network

        Returns:
            list[Node]: Points in the network.
        """

        arcpy.env.workspace = str(self.un_dataset)

        targets = arcpy.ListFeatureClasses("WaterDevice", feature_type="Point")

        if len(targets) != 1:
            raise ValueError(
                "Search for WaterDevice Was Empty or Ambiguous in target Utility Nework"
            )

        # Only return list of Nodes below for the specified subnetwork
        if subnetwork and subnetwork != "All Features":
            query = f"PressureSubnetworkName = '{subnetwork}'"
        else:
            query = None

        return [
            Node(row[0], row[1])
            for row in arcpy.da.SearchCursor(targets[0], ["GLOBALID", "SHAPE@"], query)
        ]

    @staticmethod
    def get_con_points(line: Node, points: [Node]) -> list[Node]:
        """Return any points that touch the input line. At each step of the traversal,
        we want to know which valves are connected to the current line.

        Args:
            line (Node): A line in the network.
            points (Node]): Points (e.g. Valves) in the network.

        Returns:
            list[Node]: All of the points that touch the input line.
        """

        connected_points = []

        # Simple check for spatial coincidence
        # Touches may not be as good as other methods, but works on the baseline
        for point in points:
            if line.geometry.touches(point.geometry):
                connected_points.append(point)

        return connected_points

    @staticmethod
    def get_con_lines(line: Node, lines: [Node], points: [Node]) -> list[Node]:
        """Return any lines that touch the input line. In addition to checking if the
        line is touching, we need to ensure that the line is not touching a valve
        that has been seen at the current location in the traversal.

        Args:
            line (Node): A line in the network.
            lines (Node]): Lines in the network.
            points (Node]): Points (e.g. Valves) in the network.

        Returns:
            list[Node]: All of the valid lines that touch the input line.
        """

        connected_lines = []

        for neighbor_line in lines:
            if line.geometry.touches(neighbor_line.geometry):

                # Check if the line touches a valve that has already been seen
                # If so, we don't want to add it to the connected lines because it
                # would lead to inaccurate segment results
                # E.G. We stop at valves and without this check we would return lines
                # on the wrong side of a given valve
                touches_valve = False
                for point in points:
                    if neighbor_line.geometry.touches(point.geometry):
                        touches_valve = True

                if not touches_valve:
                    connected_lines.append(neighbor_line)

        return connected_lines

    @staticmethod
    def traverse_segment(line, visited, points, lines, valves, segments):
        """Depth-first search of the network to find all connected lines and valves.

        Args:
            line (_type_): Starting point for a given segment.
            visited (_type_): Visited nodes of the traversal.
            points (_type_): Points (e.g. Valves) in the network.
            lines (_type_): Lines in the network.
            valves (_type_): Valves that are part of the segment.
            segments (_type_): Lines that are part of the segment.
        """

        visited.append(line.global_id)

        # Any points toching the current line should be part of the segment
        connected_points = SegmentBuilder.get_con_points(line, points)

        for point in connected_points:
            valves.append(point)

        # Fetch any lines that are touching the current line
        connected_lines = SegmentBuilder.get_con_lines(line, lines, connected_points)

        # Remove line to prevent duplicate traversals & improve speed of traversal
        for conn_line in connected_lines:
            try:
                lines.remove(conn_line)
            except ValueError:
                pass  # TODO - Find These Missing Lines

            segments.append(conn_line)

            # Continue through the network to find more connected lines
            SegmentBuilder.traverse_segment(
                conn_line, visited, points, lines, valves, segments
            )

    @staticmethod
    def generate_segments(points: list[Node], lines: list[Node]) -> dict[str, dict]:
        """Driver logic for segment generation. Generally, we take a line, build a segment from it,
           and then continue to build segments from any remaining lines until there are no more lines.

        Args:
            points (list[Node]): Points in the subnetwork that will be traversed.
            lines (list[Node]): Lines in the subnetwork that will be traversed.

        Returns:
            dict[str, dict]: A dictionary of segments, where each segment is a dictionary of lines and valves.
        """        

        arcpy.AddMessage("Generating Segments")

        all_segments = {}

        while lines:

            # Get a line from the list of lines & unique ID for the new segment
            start_line = lines.pop()
            segment_id = str(uuid4())

            # Create a new segment based off the start line
            all_segments[segment_id] = {"lines": [start_line], "valves": []}

            # Each segment is made up of valves and lines. In the future, valves might be too
            # restrictive of a term, so the plan is to use points wherever possible.
            visited  = []
            valves   = []
            segments = []

            # Execute DFS to find all connected lines and valves
            SegmentBuilder.traverse_segment(
                start_line, visited, points, lines, valves, segments
            )

            # Add the traversed lines and valves to the segment
            all_segments[segment_id]["valves"] = list(valves)
            all_segments[segment_id]["lines"] += list(segments)

        return all_segments

    def write_segments(self, segments: dict):
        """Write segments to the Segment_Lines and Segment_Points tables.

        Args:
            segments (dict): Segment information.
        """

        # For each Segment, write the Segment ID and the Global ID of each Line
        # Traversing the segments Dictionary twice is likely not ideal, but for
        # simplicity sake, we will do it for now.
        with arcpy.da.InsertCursor(
            self.line_table, ["segment_id", "un_line_guid"]
        ) as cursor:
            for segment_id, segment_info in segments.items():
                for line in segment_info["lines"]:
                    cursor.insertRow([segment_id, line.global_id])

        # For each Segment, write the Segment ID and the Global ID of each Point
        with arcpy.da.InsertCursor(
            self.point_table, ["segment_id", "un_point_guid"]
        ) as cursor:
            for segment_id, segment_info in segments.items():
                for valve in segment_info["valves"]:
                    cursor.insertRow([segment_id, valve.global_id])

    def execute(self, parameters, messages):
        """Primary Business Logic for the Segment Builder Tool"""

        # Unpack UN Dataset Parent & Overwrite Preference & Stage Segment Tables
        self.un_dataset = Path(arcpy.Describe(parameters[0]).catalogPath).parent
        self.overwrite  = parameters[1]
        self.stage_segment_tables()

        # Unpack Subnetwork Parameter & Get Lines & Points for Traversal
        subnetwork = parameters[2].value
        lines      = self.get_lines(subnetwork)
        points     = self.get_points(subnetwork)
        arcpy.AddMessage(
            f"Found {len(lines)} Lines & {len(points)} Points for Processing"
        )

        # Generate Segments from Lines & Points & Get Dictionary of Results Back
        segments = self.generate_segments(points, lines)
        arcpy.AddMessage(f"Generated {len(segments)} Segments")

        # Write Segments to Tables
        self.write_segments(segments)

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return


class SegmentTrace:
    """ Conduct simple 'trace' to validate segment building process. """
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Segment Trace"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        param_1 = arcpy.Parameter(
            displayName="Input Utility Network",
            name="network",
            datatype="DEUtilityNetwork",
            parameterType="Required",
            direction="Input",
        )

        param_2 = arcpy.Parameter(
            displayName="Feature Selection",
            name="tracestart",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
        )

        return [param_1, param_2]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    @staticmethod
    def fetch_lines(table: str, segment_id: str) -> list[str]:
        """Get lines for input segment id

        Args:
            table (str): Segment table to query
            segment_id (str): Segment ID to query

        Returns:
            list[str]: Global IDs of lines in the segment
        """
        un_guids = []

        # Hard coded table names are likely not required and could
        # be fixed with some better OOP inheritance
        with arcpy.da.SearchCursor(
            table, ["un_line_guid"], f"segment_id = '{segment_id}'"
        ) as cursor:
            for row in cursor:
                un_guids.append(row[0])

        return un_guids

    @staticmethod
    def fetch_points(table: str, segment_id: str) -> list[str]:
        """Get points for input segment id

        Args:
            table (str): Segment table to query
            segment_id (str): Segment ID to query

        Returns:
            list[str]: Global IDs of points in the segment
        """
        un_guids = []

        # Hard coded table names are likely not required and could
        # be fixed with some better OOP inheritance
        with arcpy.da.SearchCursor(
            table, ["un_point_guid"], f"segment_id = '{segment_id}'"
        ) as cursor:
            for row in cursor:
                un_guids.append(row[0])

        return un_guids

    @staticmethod
    def get_globalid(layer: 'Layer', oid: int) -> str:
        """Return the global id for the input oid on ArcPy Mapping Layer

        Args:
            layer (Layer): Input Layer from user
            oid (int): Object ID to match for global id

        Returns:
            str: Global ID for selected feature
        """

        where_clause = f"OBJECTID = {oid}"

        # We only expect 1 feature, sp the first is returned
        # Raising an execption is the case no records are found is likely
        # better than returning None
        with arcpy.da.SearchCursor(
            layer.dataSource, ["GLOBALID"], where_clause
        ) as cursor:
            for row in cursor:
                return row[0]

    @staticmethod
    def get_segment_id(target_table: str, target_guid: str) -> str:
        """Return the segment id for the input global id

        Args:
            target_table (str): Table to query
            target_guid (str): Global ID to match for segment id

        Returns:
            str: Segment ID
        """

        if target_table.endswith("Segment_Lines"):
            field = "un_line_guid"
        else:
            field = "un_point_guid"

        # We only expect 1 feature, sp the first is returned
        # Raising an execption is the case no records are found is likely
        # better than returning None
        with arcpy.da.SearchCursor(
            target_table, ["segment_id"], f"{field} = '{target_guid}'"
        ) as cursor:
            for row in cursor:
                return row[0]

    def apply_selection(self, target_layer: str, guids: list[str]) -> None:
        """Apply selection to active Map Layer

        Args:
            target_layer (str): Name of the layer to select against
            guids (list[str]): GUIDs to select (must be converted to oids first)
        """
        arcpy.AddMessage(f"Applying Selecction to {target_layer}")

        current_map = arcpy.mp.ArcGISProject("CURRENT").activeMap

        target_layer = [
            layer for layer in current_map.listLayers() if layer.name == target_layer
        ]

        if len(target_layer) != 1:
            arcpy.AddWarning(f"Target Layer For Selection Not Found: {target_layer}")
            arcpy.AddWarning("No Feature Will Be Selected")
            return

        target_layer = target_layer[0]

        oids = []

        # An exeception should likely be raised if no oids are found
        with arcpy.da.SearchCursor(
            target_layer.dataSource, ["OBJECTID"], f"GLOBALID IN {tuple(guids)}"
        ) as cursor:
            for row in cursor:
                oids.append(row[0])

        target_layer.setSelectionSet(oids, "NEW")

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Unpack User Inputs
        un_root_path = str(
            Path(arcpy.Describe(parameters[0]).catalogPath).parent.parent
        )
        feature_layer = parameters[1].value

        point_table = str(Path(un_root_path, "Segment_Points"))
        line_table = str(Path(un_root_path, "Segment_Lines"))

        # Get Selected Feature from Inpout Map Layer
        selected_features = list(feature_layer.getSelectionSet())

        # In the future, we will support multiple features
        if len(selected_features) != 1:
            arcpy.AddWarning("Please Select a Single Feature to Trace")
            return

        # Collect GlobalId for Lookup in Segment Tables
        target_guid = self.get_globalid(feature_layer, selected_features[0])
        arcpy.AddMessage(f"Fetching Segment by GlobalId: {target_guid}")

        shape_type = arcpy.Describe(feature_layer).shapeType

        if shape_type == "Point":
            segment_id = self.get_segment_id(point_table, target_guid)
            arcpy.AddMessage(f"Target Segment ID: {segment_id}")

            guids = self.fetch_lines(line_table, segment_id)
            arcpy.AddMessage(f"Found {len(guids)} Lines in Segment")
            self.apply_selection("Water Line", guids)

        elif shape_type == "Polyline":
            segment_id = self.get_segment_id(line_table, target_guid)
            arcpy.AddMessage(f"Target Segment ID: {segment_id}")

            guids = self.fetch_points(point_table, segment_id)
            arcpy.AddMessage(f"Found {len(guids)} Points in Segment")
            self.apply_selection("Water Device", guids)

        else:
            arcpy.AddMessage(f"Please Start with Point or Polyline, Not {shape_type}")
            return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
