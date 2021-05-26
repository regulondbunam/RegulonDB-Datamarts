class Summary():
    def __init__(self, regulated_objects, reg_ints):
        self.regulated_objects = regulated_objects
        self.reg_ints = reg_ints

    def get_summary(self):
        summary_object = {}
        for obj, regulated_list in self.regulated_objects.items():
            count = get_counts_of_regulated_object(regulated_list)
            summary_object.update({obj: count})
        summary_object = insert_ri_bs_counts(self.reg_ints, summary_object)
        return summary_object


def get_counts_of_regulated_object(regulated_list):
    count_object = {"repressed": 0, "activated": 0, "dual": 0, "unknown": 0, "total": len(regulated_list)}
    for regulated_item in regulated_list:
        if regulated_item["function"]:
            if regulated_item["function"] == "repressor":
                count_object["repressed"] = count_object.get("repressed", 0) + 1
            if regulated_item["function"] == "activator":
                count_object["activated"] = count_object.get("activated", 0) + 1
            if regulated_item["function"] == "dual":
                count_object["dual"] = count_object.get("dual", 0) + 1
        else:
            count_object["unknown"] = count_object.get("unknown", 0) + 1
    return count_object


def insert_ri_bs_counts(reg_ints, summary_object):
    binding_sites = []
    reg_int_count = get_counts_of_regulated_object(reg_ints)
    summary_object.update({"regulatoryInteractions": reg_int_count})
    for reg_int in reg_ints:
        if reg_int["regulatoryBindingSites"]:
            binding_sites_object = reg_int["regulatoryBindingSites"]
            binding_sites_object["function"] = reg_int["function"]
            if binding_sites_object not in binding_sites:
                binding_sites.append(binding_sites_object)
    binding_sites_count = get_counts_of_regulated_object(binding_sites)
    summary_object.update({"bindingSites": binding_sites_count})
    return summary_object
