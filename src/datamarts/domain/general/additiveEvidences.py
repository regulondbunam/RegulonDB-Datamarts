import json
import itertools


def load_rules():
    with open("./lib/evidences_rules.json", 'r') as add_json_rules:
        rules = json.load(add_json_rules)
        return rules


class AdditiveEvidences:
    def __init__(self, all_citations):
        self.additive_evidences = all_citations
        self.citations = all_citations

    @property
    def additive_evidences(self):
        return self._additive_evidences

    @additive_evidences.setter
    def additive_evidences(self, all_citations):
        all_rules = load_rules()
        evs_with_rules = []
        for citation in all_citations:
            if citation.get("evidence", None):
                evidence = citation.get("evidence", None)
                if evidence.get('additiveEvidenceCodeRule', None):
                    if evidence not in evs_with_rules:
                        evs_with_rules.append(evidence)
        evs_found_accepted_rules = get_all_valid_rules(evs_with_rules, all_rules)
        additive_evidences = build_additive_evidence(evs_found_accepted_rules, evs_with_rules)

        self._additive_evidences = additive_evidences

    def to_dict(self):
        return self.additive_evidences

    def get_confidence_level(self):
        if len(self.additive_evidences) > 0:
            for add_ev in self.additive_evidences:
                if add_ev.get('type', None) == "C":
                    return "C"
            for add_ev in self.additive_evidences:
                if add_ev.get('type', None) == "S":
                    return "S"
            for add_ev in self.additive_evidences:
                if add_ev.get('type', None) == "W":
                    return "W"
        else:
            for citation in self.citations:
                if citation.get("evidence", None):
                    ev = citation.get("evidence", None)
                    if ev.get('type', None) == "C":
                        return "C"
            for citation in self.citations:
                if citation.get("evidence", None):
                    ev = citation.get("evidence", None)
                    if ev.get('type', None) == "S":
                        return "S"
            for citation in self.citations:
                if citation.get("evidence", None):
                    ev = citation.get("evidence", None)
                    if ev.get('type', None) == "W":
                        return "W"
            return None


def get_rules_with_code(code_rule, rules):
    valid_rules = []
    for rule in rules:
        if code_rule in rule.get('rules_values'):
            valid_rules.append(rule)
        else:
            continue
    return valid_rules


def remove_rule_with_code(code_rule, rules):
    valid_rules = []
    for rule in rules:
        if code_rule not in rule.get('rules_values'):
            valid_rules.append(rule)
        else:
            continue
    return valid_rules


def get_all_valid_rules(evs_with_rules, rules):
    valid_rules = []
    found_codes = []
    codes = []
    for evidence in evs_with_rules:
        code_rule = evidence.get('additiveEvidenceCodeRule', None)
        if code_rule not in found_codes:
            found_codes.append(code_rule)
        left_rules = get_rules_with_code(code_rule, rules)
        for rule in left_rules:
            if rule not in valid_rules:
                valid_rules.append(rule)
    for rule in valid_rules:
        for code in rule.get("rules_values", None):
            if code not in codes:
                codes.append(code)
    for additional_code in codes:
        if additional_code not in found_codes:
            valid_rules = remove_rule_with_code(additional_code, valid_rules)
    return valid_rules


def get_all_rules_combination(rule, evidences):
    rule_dict = {}
    rule_list = []
    for code in rule.get("rules_values", None):
        for ev in evidences:
            if ev.get('additiveEvidenceCodeRule', None) == code:
                rule_dict.setdefault(f"{code}", []).append(ev)
    for code in rule.get("rules_values", None):
        rule_list.append(rule_dict.get(f"{code}"))
    return list(itertools.product(*rule_list))


def find(lst, key, value):
    for i, dic in enumerate(lst):
        if dic.get(key, None) == value:
            return dic
    return -1


def build_additive_evidence(rules, evidences):
    additive_rules = []
    for rule in rules:
        rules_comb = get_all_rules_combination(rule, evidences)
        for i in rules_comb:
            code = ""
            codes = []
            for v in i:
                code = code + v.get("code", None) + "/"
                codes.append(v.get('additiveEvidenceCodeRule', None))
            rule = find(rules, "rules_values", codes)
            additive_rules.append({
                "code": f"AE({code[:-1]})",
                "type": rule.get("type", None),
                "category": rule.get("evidenceCategory", None)
            })
    return additive_rules
