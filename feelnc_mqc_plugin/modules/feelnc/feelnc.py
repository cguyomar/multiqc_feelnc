from multiqc.modules.base_module import BaseMultiqcModule
import logging
from multiqc import config
from multiqc.plots import table, linegraph, bargraph
from copy import deepcopy
from collections import OrderedDict
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Detection of long non coding transcripts (FEELnc)', anchor='feelnc',
        href="https://github.com/tderrien/FEELnc",
        info="results for the transcripts of the novel annotation. FEELnc is a classifier that identifies and characterizes long non coding RNAs in a set of annotated transcripts.")

        target="FEELnc"
        self.mnameshort = '<a href="{}" target="_blank">{}</a>'.format(self.href, target)
        self.intro = "<p>{} {}</p>{}".format(self.mnameshort, self.info, self.extra)

        self.feelnc_classes_data = dict()
        self.feelnc_classes_counts_all = dict()
        self.feelnc_classes_counts_sense = dict()
        self.feelnc_classes_counts_antisense = dict()
        self.feelnc_roc_curves = dict()
        self.feelnc_classification_summary = dict()
        self.feelnc_filter_counts = dict()

        # Parse logs

        # Read filter log for unclassified transcripts
        for f in self.find_log_files("feelnc/filter_log", filehandles=True):
            self.feelnc_filter_counts["Transcripts of the novel annotation"] = self.parse_filter_log(f)
            break

        # Read classification summary file
        for f in self.find_log_files("feelnc/classification_summary", filehandles=False):
            self.feelnc_classification_summary["Transcripts of the novel annotation"] = self.parse_asis(f)
            break

        # Add unclassified transciprts to the summary
        classification_counts = self.feelnc_classification_summary["Transcripts of the novel annotation"]
        filter_counts = self.feelnc_filter_counts["Transcripts of the novel annotation"]
        classification_counts["Not evaluated by FEELnc (overlapping coding transcript in sense)"] = 0
        classification_counts["Not evaluated by FEELnc (other reason)"] = 0

        if "overlap" in filter_counts.keys():
            classification_counts["Not evaluated by FEELnc (coding transcripts)"] = filter_counts["overlap"]
        for val in ["monoexonic","biexonic","Size"]:
            if val in filter_counts.keys():
                classification_counts["Not evaluated by FEELnc (other reason)"] += filter_counts[val]

        # Parse ROC curve
        for f in self.find_log_files("feelnc/roc", filehandles=True):
           self.feelnc_roc_curves[f['s_name']] = self.parse_roc(f)

        # Parse interactions classes
        for f in self.find_log_files("feelnc/lnc_classes", filehandles=True):
            self.feelnc_classes_data["lncRNA class"]= self.parse_classes(f)
            self.feelnc_classes_counts_all["lncRNA class"],self.feelnc_classes_counts_sense["lncRNA class"],self.feelnc_classes_counts_antisense[f["s_name"]] = self.count_classes(f)
            break

        # Make report
        config = {'id':'classification_summary',
        'title':'Summary of FEELnc classification'
        }
        config_table = {"id": "feelnc_classification", "namespace": "feelnc","sortRows": False,"title": "Classification of transcripts by FEELnc"}
        config = {"id": "feelnc_classification", "namespace": "feelnc","sortRows": False}
        self.add_section(
            name = 'Transcript classification',
            description = "For all potential lnc transcripts (excluding known coding and short transcripts), FEELnc uses a random forest classifier to determine the coding potential of transcripts",
            anchor = 'feelnc_classification',
            plot =  bargraph.plot(self.feelnc_classification_summary,pconfig=config_table)
        )

        if 'exons_RF_TGROC' in self.feelnc_roc_curves.keys():
            self.add_section(
                name = 'ROC curve',
                description = "ROC curve of the FEELnc classifier",
                anchor = 'feelnc_roc',
                content = self.feelnc_roc_curves['exons_RF_TGROC']
            )

        cats = OrderedDict()
        cats['Genic - Exonic'] = {
            'color': '#b2d1ff'
        }
        cats['Genic - Intronic'] = {
            'color': '#2660c1'
        }
        cats['Intergenic - Downstream'] = {
            'color': '#b32a2a'
        }
        cats['Intergenic - Upstream'] = {
            'color': '#fe6e62'
        }
        config = {"data_labels": [{'name': 'All'}, {'name': 'Sense'}, {'name': 'Antisense'}]}

        self.add_section(
            name = 'LncRNA position with respect to the closest mRNA',
            description = "The relative position of each lncRNA transcript is \
            assessed with respect to the closest coding transcript. \
            Coding transcripts are either transcripts from the reference whose\
            'transcript biotype' is 'protein_coding' or transcripts evaluated as \
            'mRNA' by FEELnc",
            anchor = 'feelnc_classification',
            plot = bargraph.plot([self.feelnc_classes_counts_all,self.feelnc_classes_counts_sense,self.feelnc_classes_counts_antisense],[cats,cats,cats],config),
        )

    def parse_filter_log(self,f):
        counts = dict()
        for line in f['f']:
            if line.split(" ")[0] == "Filter":
                reason = line.split(" ")[1]
                if reason in counts.keys():
                    counts[reason] += 1
                else:
                    counts[reason] = 1
        return counts

    def parse_asis(self,f):
        parsed_data = {}
        raw_data = f['f']
        raw_data = raw_data.split("\n")
        for l in raw_data:
            l = l.split("\t")
            if len(l) >= 2:
                parsed_data[l[0]] = l[1]
        return(parsed_data)

    def parse_roc(self,f):
        '<div class="mqc-custom-content-image"><img src="{}" /></div>'.format(f['fn'])

    def parse_classes(self,f):
        data = {}
        line_number = 2
        col_names = []
        for l in f['f']:
            if col_names == []:
                col_names = l.split("\t")
                col_names = [s.strip() for s in col_names]
                log.info(col_names)
            else:
                row = dict(zip(col_names,l.split("\t")))
                row = {k:v.strip() for k, v in row.items()}
                row["isBest"] = bool(int(row["isBest"]))
                row["distance"] = int(row["distance"])

                data[str(line_number)] = row
                line_number += 1
        return(data)

    def count_classes(self,f):
        template = {'Genic - Exonic':0,
                    'Genic - Intronic':0,
                    'Intergenic - Upstream':0,
                    'Intergenic - Downstream':0}
        classes_counts = dict(template)
        classes_counts_sense = dict(template)
        classes_counts_antisense = dict(template)
        subclasses_counts = OrderedDict()

        for r in self.feelnc_classes_data["lncRNA class"].values():
            if r["isBest"]:
                subclass = ' '.join([r["type"],r["direction"],r["subtype"],r["location"]])
                if subclass in subclasses_counts.keys():
                    subclasses_counts[subclass]['count'] += 1
                else:
                    subclasses_counts[subclass] = {'direction': r["direction"],
                                                    'type': r["type"],
                                                    'subtype': r["subtype"],
                                                    'location': r["location"],
                                                    'count': 0}
                classes_counts[r["type"].title() + " - " + r["location"].title()]+=1
                if r["direction"]=="sense":
                    classes_counts_sense[r["type"].title() + " - " + r["location"].title()]+=1
                else:
                    classes_counts_antisense[r["type"].title() + " - " + r["location"].title()]+=1

        # Rename dict keys for readability
        new_keys = range(len(subclasses_counts.keys()))
        subclasses_counts_renamed = dict(zip(new_keys,subclasses_counts.values()))
        log.info(classes_counts)
        return([classes_counts,classes_counts_sense,classes_counts_antisense])
