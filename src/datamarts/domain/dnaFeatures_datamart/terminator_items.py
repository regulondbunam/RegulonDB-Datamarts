import multigenomic_api


class TerminatorDNAFeatures:

    @property
    def objects(self):
        # Terminator Objects
        terminator_objects = multigenomic_api.terminators.get_all()
        for terminator in terminator_objects:
            dtt_datamart = TerminatorDNAFeatures.DTTDatamart(terminator)
            yield dtt_datamart
            print(terminator.id)
        del terminator_objects

    class DTTDatamart:

        def __init__(self, entity):
            self.entity = entity

        def to_dict(self):
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "leftEndPosition": self.entity.transcriptionTerminationSite.left_end_position,
                "lineRGBColor": "0,0,0",
                # TODO: this gonna be defined by evidence, if is weak or strong
                "lineType": 1,
                "lineWidth": 1,
                "objectType": "terminator",
                "objectRGBColor": "0,0,0",
                "organism": {
                    "organism_id": self.entity.organisms_id
                },
                "rightEndPosition": self.entity.transcriptionTerminationSite.right_end_position,
                # TODO: ask what will be showed on this tooltip
                # "tooltip": f"{self.entity.name}; "
                           # f"{self.positions['leftEndPosition']}->{self.positions['rightEndPosition']}; "
                           # f"{self.product}"
            }
            return dttDatamart