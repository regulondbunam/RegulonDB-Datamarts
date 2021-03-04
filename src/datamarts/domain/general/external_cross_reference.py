class ExternalCrossReference:

    def __init__(self, external_cross_reference, object_id):
        self.external_cross_reference = external_cross_reference
        self.object_id = object_id

    def to_json(self):
        external_cross_reference = {
            "externalCrossReferenceId": self.external_cross_reference.id,
            "externalCrossReferenceName": self.external_cross_reference.name,
            "objectId": self.object_id,
            "url": ExternalCrossReference.replace_url(self.external_cross_reference.url, self.object_id)
        }
        return external_cross_reference

    @staticmethod
    def replace_url(url, object_id):
        if "/A" in url[-2:]:
            url_parsed = url
        else:
            url_parsed = url.replace("~A", object_id)
        return url_parsed