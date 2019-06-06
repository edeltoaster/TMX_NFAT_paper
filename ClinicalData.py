class ClinicalData:
    def __init__(self, path, patients = None):
        self.header = []
        self.patient_data = {}
        # header:
        # case_id	submitter_id	project_id	gender	year_of_birth	race	ethnicity	year_of_death	classification_of_tumor	last_known_disease_status	primary_diagnosis	tumor_stage	age_at_diagnosis	vital_status	morphology	days_to_death	days_to_last_known_disease_status	days_to_recurrence	tumor_grade	tissue_or_organ_of_origin	days_to_birth	progression_or_recurrence	prior_malignancy	site_of_resection_or_biopsy	days_to_last_follow_up	therapeutic_agents	days_to_treatment	treatment_intent_type	treatment_or_therapy
        for line in open(path):
            line = line.strip().split("\t")
            
            if line[0] == "case_id":
                # read header
                self.header = line
            else:
                # read data
                self.patient_id = line[1]
                if patients == None or self.patient_id in patients: # add only subset
                    self.patient_data[self.patient_id] = line
                    
    def getAnnotations(self):
        return self.header
        
    def getAnnotation(self, annotation):
        data = {}
        
        # find annotation_index
        ann_ind = [i for i,x in enumerate(self.header) if x == annotation]
        if len(ann_ind) == 0:
            return data
        ann_ind = ann_ind[0]
        
        # set data
        for patient in self.patient_data:
            data[patient] = self.patient_data[patient][ann_ind]
            
        return data
    
    def getSurvivalData(self, patients = None):
        df = self.getAnnotation("days_to_last_follow_up")
        dd = self.getAnnotation("days_to_death")
        alive = self.getAnnotation("vital_status")
        data = {}
        for p in df:
            if patients != None and not p in patients:
                continue
            if alive[p] == "alive":
                data[p] = (df[p], "False")
            else:
                data[p] = (dd[p], "True")
                
        return data
    
    def writeSurvCSV(self, path, patient_group_map):
        surv = self.getSurvivalData(patient_group_map.keys())
        open(path, "w").writelines(["patient,time,status,group\n"] + [p + "," + surv[p][0] + "," + surv[p][1] + "," + patient_group_map[p] + "\n" for p in surv])
        
    def getGender(self):
        return self.getAnnotation("gender")
    
    def getAge(self):
        age = self.getAnnotation("age_at_diagnosis")
        for p in age:
            if age[p] == "--":
                age[p] = float("nan")
            else:
                a = float(age[p])
                age[p] = a / 365.25
        return age
        