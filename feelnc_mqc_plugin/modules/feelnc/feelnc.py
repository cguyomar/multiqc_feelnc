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
        super(MultiqcModule, self).__init__(name='FeelNC', anchor='feelnc',
        href="https://github.com/tderrien/FEELnc",
        info="FEELnc allows to annotate long non coding RNAs (lncRNAs) based on reconstructed transcripts from RNA-seq data ")

        self.feelnc_roc_curves = dict()
        self.feelnc_classes_data = dict()
        self.feelnc_classes_counts_all = dict()
        self.feelnc_classes_counts_sense = dict()
        self.feelnc_classes_counts_antisense = dict()
        self.feelnc_classifier_stats = dict()
        self.feelnc_transcripts_counts = OrderedDict()

        # Parse logs
        for f in self.find_log_files("feelnc/roc", filehandles=True):
           self.feelnc_roc_curves[f['s_name']] = self.parse_roc(f)

        for f in self.find_log_files("feelnc/lnc_classes_log", filehandles=False):
            self.feelnc_classifier_stats[f["s_name"]] = self.parse_log(f)

        for f in self.find_log_files("feelnc/lnc_classes", filehandles=True):
            self.feelnc_classes_data[f["s_name"]]= self.parse_classes(f)
            self.feelnc_classes_counts_all[f["s_name"]],self.feelnc_classes_counts_sense[f["s_name"]],self.feelnc_classes_counts_antisense[f["s_name"]] = self.count_classes(f)
        log.info(self.feelnc_classes_counts_all)

        for f in self.find_log_files("feelnc/rf_summary", filehandles=False):
            self.feelnc_transcripts_counts[f["s_name"]] = self.parse_log(f)

        # Make report
        config = {'col1_header': 'FeelNC run'}
        self.add_section(
            name = 'Coding potential prediction',
            description = 'Coding potential among transcripts filtered by feelnc_filter',
            plot = table.plot(self.feelnc_transcripts_counts,{},config)
        )

        self.add_section(
            name = 'ROC curve',
            anchor = 'feelnc_roc',
            content = self.feelnc_roc_curves['exons_RF_TGROC']
        )

        table_config = {
            'namespace': '',                         # Name for grouping. Prepends desc and is in Config Columns modal
            'id': '<random string>',                 # ID used for the table
            'table_title': '<table id>',             # Title of the table. Used in the column config modal
            'save_file': False,                      # Whether to save the table data to a file
            'sortRows': True,                         # Whether to sort rows alphabetically
            'only_defined_headers': True,             # Only show columns that are defined in the headers config
            'col1_header': 'FeelNC run'     ,        # The header used for the first column
            'no_beeswarm': False    # Force a table to always be plotted (beeswarm by default if many rows)
        }

        self.add_section(
            name = 'lncRNA classification summary',
            anchor = 'feelnc_classification',
#            plot = table.plot(self.feelnc_classes_counts['lncRNA_classes'])
            plot = table.plot(self.feelnc_classifier_stats,{},config),
        )
        cats = OrderedDict()
        cats['genic - exonic'] = {
            'color': '#b2d1ff'
        }
        cats['genic - intronic'] = {
            'color': '#2660c1'
        }
        cats['intergenic - downstream'] = {
            'color': '#b32a2a'
        }
        cats['intergenic - upstream'] = {
            'color': '#fe6e62'
        }

        config_table = {"id": "feelnc_classes_table", "namespace": "feelnc","sortRows": False}
        config = {"data_labels": [{'name': 'All'}, {'name': 'Sense'}, {'name': 'Antisense'}]}
        self.add_section(
            name = 'lncRNA classes',
            description = "Distribution of lncRNA classes, only for bestHits",
            anchor = 'feelnc_classes',
            plot = bargraph.plot([self.feelnc_classes_counts_all,self.feelnc_classes_counts_sense,self.feelnc_classes_counts_antisense],[cats,cats,cats],config)
        )

        # Whole classes table
        # self.add_section(
        #     name = 'LNC classes',
        #     anchor = 'feelnc_classes',
        #     plot = table.plot(self.feelnc_classes_data["lncRNA_classes"])
        # )

    def parse_log(self,f):
        regexes = {
            "Number of lncRNAs": r"#Number of lncRNA : (\d+)",
            "Number of mRNAs": r"#Number of mRNA : (\d+)",
            "Number of interactions": r"#Number of interaction : (\d+)",
            "Number of lncRNAs without interaction": r"#Number of lncRNA without interaction : (\d+)",
            "Number of lncRNA transcripts": r"-Nb_lncRNAs:\t(\d+)",
            "Number of mRNA transcripts": r"-Nb_mRNAs:\t(\d+)",
        }

        raw_data = f['f']
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
        log.info(parsed_data)
        return(parsed_data)

    def parse_roc(self,f):
        img_html = ('<div class="mqc-custom-content-image"><img src="{}" /></div>'.format(f['fn']))
        return(img_html)

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
                row = { k:v.strip() for k, v in row.items()}
                row["isBest"] = bool(int(row["isBest"]))
                row["distance"] = int(row["distance"])

                data[str(line_number)] = row
                line_number += 1

        return(data)

    def count_classes(self,f):
        #classes_counts = {'nb_lincrna':0,'nb_antisense':0,'nb_sense_intronic':0,'nb_sense_exonic':0}
        classes_counts = {'genic - exonic':0,'genic - intronic':0,'intergenic - upstream':0,'intergenic - downstream':0}
        classes_counts_sense = {'genic - exonic':0,'genic - intronic':0,'intergenic - upstream':0,'intergenic - downstream':0}
        classes_counts_antisense = {'genic - exonic':0,'genic - intronic':0,'intergenic - upstream':0,'intergenic - downstream':0}
        subclasses_counts = OrderedDict()
        subclasses = []

        for r in self.feelnc_classes_data[f["s_name"]].values():
            if r["isBest"]:
                subclass = r["type"] + ' ' + r["direction"] + ' ' + r["subtype"] + ' ' + r["location"]
                if subclass in subclasses_counts.keys():
                    subclasses_counts[subclass]['count'] += 1
                else:
                    subclasses_counts[subclass] = {'direction': r["direction"], 'type': r["type"], 'subtype': r["subtype"], 'location': r["location"], 'count': 0}
                classes_counts[r["type"] + " - " + r["location"]]+=1
                if r["direction"]=="sense":
                    classes_counts_sense[r["type"] + " - " + r["location"]]+=1
                else:
                    classes_counts_antisense[r["type"] + " - " + r["location"]]+=1

        # Rename dict keys for readability
        new_keys = list(range(len(subclasses_counts.keys())))
        subclasses_counts_renamed = dict(zip(new_keys,list(subclasses_counts.values())))
        log.info(classes_counts)
        #return({"all":classes_counts,"sense":classes_counts_sense,"antisense":classes_counts_antisense})
        return([classes_counts,classes_counts_sense,classes_counts_antisense])


    # def parse_gtf(self,f):
    #     for l in f['f']:
    #         row = l.split("\t")
    #         if row[0]==1
    #
