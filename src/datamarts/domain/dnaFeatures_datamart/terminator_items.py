import multigenomic_api


class TerminatorDNAFeatures:

    @property
    def objects(self):
        # Terminator Objects
        terminator_objects = multigenomic_api.terminators.get_all()
        for terminator in terminator_objects:
            dtt_datamart = TerminatorDNAFeatures.DTTDatamart(terminator)
            yield dtt_datamart
        del terminator_objects

    class DTTDatamart:

        def __init__(self, entity):
            self.entity = entity
            self.line_type = entity.citations

        @property
        def line_type(self):
            return self._line_type

        @line_type.setter
        def line_type(self, citations):
            self._line_type = 1
            if citations:
                for citation in citations:
                    if citation.evidences_id:
                        evidence = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        if evidence.type == "S":
                            self._line_type = 2
                            break
                        elif evidence.type == "C":
                            self._line_type = 3
                            break

        def to_dict(self):
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "leftEndPosition": self.entity.transcriptionTerminationSite.left_end_position,
                "lineRGBColor": "0,0,0",
                # TODO: this gonna be defined by evidence, if is weak or strong
                "lineType": self._line_type,
                "lineWidth": 1,
                "objectType": "terminator",
                "objectRGBColor": "0,0,0",
                "organism": {
                    "organism_id": self.entity.organisms_id
                },
                "rightEndPosition": self.entity.transcriptionTerminationSite.right_end_position,
                # TODO: ask what will be showed on this tooltip
                "tooltip": f"{self.entity.class_}; \n"
                           f"Genome Position: "
                           f"{self.entity.transcriptionTerminationSite.left_end_position}-"
                           f"{self.entity.transcriptionTerminationSite.right_end_position}; "
            }
            return dttDatamart