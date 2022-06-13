import re
from typing import Dict

import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt

from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

from Scripts.core import Model
from Scripts.utils import ReconstructionTool
from Scripts.xrefs_converters import ReactionsConverter


class ReactionsAssessor:

    def __init__(self, reference_model: Model, conversion_method: str, reference_model_format: str):
        self.reference_model = reference_model
        self.reference_model_format = reference_model_format
        self.conversion_method = conversion_method

        if self.conversion_method == "modelseed":
            self.reactions_converter = ReactionsConverter("Xrefs files/ModelSEED-reactions.csv")
        else:
            self.reactions_converter = ReactionsConverter("Xrefs files/MetaNetX-reactions.csv")

        self.reference_reaction_sets = self.convert_reference_model_with_converter()
        self.general_reference_reaction_set = []

        for reactions in list(self.reference_reaction_sets.values()):
            self.general_reference_reaction_set.extend([re.upper() for re in reactions])
        self.general_reference_reaction_set = set(self.general_reference_reaction_set)

        self.y_true = None
        self.reference_model_reactions = None

    @staticmethod
    def get_reactions_to_convert(model):
        reactions_not_to_convert = []
        all_reactions = model.model.reactions
        reactions_to_convert = []

        transport_reactions = {}

        exchanges = model.model.exchanges
        demands = model.model.demands
        sinks = model.model.sinks

        for group in model.model.groups:
            reactions_list = []
            if "transport" in str(group.name).lower() or "drain" in str(group.name).lower():
                for member in group.members:
                    reactions_list.append(member.id)
                transport_reactions[str(group.name)] = reactions_list

        for reaction in all_reactions:
            reaction_id = reaction.id

            add = True

            # exclude exchange reactions from the comparison
            if len(exchanges) > 0 and add:
                if reaction in exchanges:
                    add = False

            if "drain" in str(reaction.name).lower() and add:
                add = False

            if "exchange" in str(reaction.name).lower() and add:
                add = False

            if str(reaction_id).upper().startswith("EX_") and add:
                add = False

            # exclude demand reactions from the comparison
            if len(demands) > 0 and add:
                if reaction in demands:
                    add = False

            # exclude transport reactions from the comparison
            if len(sinks) > 0 and add:
                if reaction in sinks:
                    add = False

            # exclude transport reactions from the comparison
            if "transport" in str(reaction.name).lower() and add:
                add = False

            if str(reaction_id).upper().startswith("TRANS-RXN") and add:
                add = False

            if "TRANS-RXN" in str(reaction.name).upper() and add:
                add = False

            if len(transport_reactions) > 0 and add:
                for group in transport_reactions.values():
                    if reaction_id in group:
                        add = False
                        break

            if add:
                reactions_to_convert.append(reaction_id)
            else:
                reactions_not_to_convert.append(reaction_id)
        return reactions_to_convert

    def convert_reference_model_with_converter(self):
        reaction_sets = {}
        reactions_to_convert = ReactionsAssessor.get_reactions_to_convert(self.reference_model)

        self.reference_model.reaction_converter = self.reactions_converter
        ModelSEED_report = self.reference_model.get_reactions_other_version(database=self.reference_model_format,
                                                                            reactions=reactions_to_convert,
                                                                            preprocess_ids=True)
        ModelSEED_convertable_reactions = ModelSEED_report.convertable

        print("Tool: " + str(self.reference_model.reconstruction_tool))
        print("Total number of reactions: " + str(len(self.reference_model.model.reactions)))
        print("Removed reactions: " + str(len(self.reference_model.model.reactions) - len(reactions_to_convert)))

        print("Reactions converted: " + str(len(ModelSEED_convertable_reactions.keys())))
        print("Reactions not converted: " + str(len(ModelSEED_report.non_convertable)))

        for convertable_reaction in ModelSEED_convertable_reactions:
            converted_reactions = ModelSEED_convertable_reactions[convertable_reaction]

            for converted_reaction in converted_reactions:
                if convertable_reaction in reaction_sets:
                    reaction_sets[convertable_reaction].append(converted_reaction.upper())

                else:
                    reaction_sets[convertable_reaction] = [converted_reaction.upper()]

        return reaction_sets

    def convert_model_reactions_with_converter(self, model):
        reaction_set = []
        model_info = {}

        reactions_to_convert = ReactionsAssessor.get_reactions_to_convert(model)
        # reactions_to_convert = [reaction.id for reaction in model.model.reactions]

        model.reaction_converter = self.reactions_converter
        ModelSEED_report = model.get_reactions_other_version(database=self.reference_model_format,
                                                             reactions=reactions_to_convert,
                                                             preprocess_ids=True)
        ModelSEED_convertable_reactions = ModelSEED_report.convertable

        print("Tool: " + str(model.reconstruction_tool))
        print("Total number of reactions: " + str(len(model.model.reactions)))

        model_info["reactions number"] = len(model.model.reactions)
        model_info["removed reactions"] = len(model.model.reactions) - len(reactions_to_convert)
        model_info["converted reactions"] = len(ModelSEED_convertable_reactions.keys())
        model_info["non-converted reactions"] = len(ModelSEED_report.non_convertable)

        print("Reactions converted: " + str(len(ModelSEED_convertable_reactions.keys())))
        print("Reactions not converted: " + str(len(ModelSEED_report.non_convertable)))

        for convertable_reaction in ModelSEED_convertable_reactions:
            converted_reactions = ModelSEED_convertable_reactions[convertable_reaction]

            if len(converted_reactions) > 1:
                found_reaction = False
                for converted_reaction in converted_reactions:

                    if converted_reaction in self.general_reference_reaction_set:
                        found_reaction = True
                        reaction_set.append(converted_reaction)
                        break

                if not found_reaction:
                    reaction_set.append(converted_reactions[0])

            else:
                for converted_reaction in converted_reactions:
                    if converted_reaction not in reaction_set:
                        reaction_set.append(str(converted_reaction).upper())

        duplicated = len(ModelSEED_convertable_reactions.keys()) - len(set(reaction_set))
        model_info["duplicated reactions"] = duplicated
        return set(reaction_set), model_info

    def convert_reactions(self, models):
        reaction_sets = {}
        models_info = {}

        for model_name in models.keys():
            reaction_set, model_info = self.convert_model_reactions_with_converter(models[model_name])
            models_info[model_name] = model_info
            reaction_sets[model_name] = reaction_set

        return reaction_sets, models_info

    def get_confusion_matrix_cells_for_reactions(self, assessed_model_reactions):
        """

        :param assessed_model_reactions:
        :return: TP, FP, FN
        """

        true_positives = 0
        false_positives = 0
        false_negatives = 0
        for reaction in assessed_model_reactions:
            if reaction in self.general_reference_reaction_set:
                true_positives += 1
            elif reaction not in self.general_reference_reaction_set:
                false_positives += 1

        for reaction in self.reference_reaction_sets:
            found = False
            for converted_reaction in self.reference_reaction_sets[reaction]:
                if converted_reaction in assessed_model_reactions:
                    found = True
                    break

            if not found:
                false_negatives += 1

        return true_positives, false_positives, false_negatives

    def get_confusion_matrix_cells_for_genes(self, assessed_model_genes):

        if self.reference_model.reconstruction_tool == ReconstructionTool.T_GONDII_CURATED.value:
            reference_model_genes = [gene.lower() for gene in self.reference_model.model.genes]
        else:
            reference_model_genes = [gene.id.lower() for gene in self.reference_model.model.genes]
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        for gene in assessed_model_genes:
            if gene in reference_model_genes:
                true_positives += 1
            else:
                false_positives += 1

        for gene in reference_model_genes:

            if gene not in assessed_model_genes:
                false_negatives += 1

        return true_positives, false_positives, false_negatives

    @staticmethod
    def get_recall(true_positive, false_negative):
        return true_positive / (true_positive + false_negative)

    @staticmethod
    def get_precision(true_positive, false_positive):
        return true_positive / (true_positive + false_positive)

    @staticmethod
    def get_f1_score(recall, precision):
        return 2 * (recall * precision) / (recall + precision)

    @staticmethod
    def get_coverage_additional_information_ratio(true_positive, false_positive):
        return true_positive / false_positive

    @staticmethod
    def get_jaccard_distance(true_positive, false_positive, false_negative):
        return 1 - (true_positive / (true_positive + false_positive + false_negative))


class ResultsReport:

    def __init__(self, reference_model: Model, models_to_be_assessed: Dict[str, Model],
                 reference_model_format: str,
                 conversion_method: str):

        self.reference_model = reference_model

        if self.reference_model.reconstruction_tool == ReconstructionTool.T_GONDII_CURATED.value:
            reactions = [reaction.id for reaction in reference_model.model.reactions]
            genes = []
            for reaction in reactions:
                genes_in_reaction = [x.group() for x in re.finditer("TGME[0-9]+_[0-9]+", reaction)]
                genes.extend(genes_in_reaction)

            self.reference_model.model.genes = genes

        self.model_assessor = ReactionsAssessor(reference_model,
                                                reference_model_format=reference_model_format,
                                                conversion_method=conversion_method)
        self.models_to_be_assessed = models_to_be_assessed

    def generate_reactions_report(self, file_path):

        model_reactions, models_info = self.model_assessor.convert_reactions(self.models_to_be_assessed)
        report_df = DataFrame(columns=["model", "reactions number", "removed reactions", "converted reactions",
                                       "non-converted reactions", "recall", "precision", "f1"])

        for i, model in enumerate(model_reactions):
            true_positives, false_positives, false_negatives = \
                self.model_assessor.get_confusion_matrix_cells_for_reactions(model_reactions[model])

            precision = self.model_assessor.get_precision(true_positives, false_positives)
            recall = self.model_assessor.get_recall(true_positives, false_negatives)
            f1_score = self.model_assessor.get_f1_score(recall, precision)
            ratio = self.model_assessor.get_coverage_additional_information_ratio(true_positives, false_positives)
            jaccard_distance = self.model_assessor.get_jaccard_distance(true_positives, false_positives,
                                                                        false_negatives)

            report_df.at[i, "model"] = model
            report_df.at[i, "reactions number"] = models_info[model]["reactions number"]
            report_df.at[i, "removed reactions"] = models_info[model]["removed reactions"]
            report_df.at[i, "converted reactions"] = models_info[model]["converted reactions"]
            report_df.at[i, "non-converted reactions"] = models_info[model]["non-converted reactions"]
            report_df.at[i, "duplicated reactions"] = models_info[model]["duplicated reactions"]
            report_df.at[i, "true positives"] = true_positives
            report_df.at[i, "false positives"] = false_positives
            report_df.at[i, "false negatives"] = false_negatives
            report_df.at[i, "recall"] = recall
            report_df.at[i, "precision"] = precision
            report_df.at[i, "f1"] = f1_score
            report_df.at[i, "ratio"] = ratio
            report_df.at[i, "jaccard distance"] = jaccard_distance

        report_df.to_csv(file_path, index=False)

    def generate_genes_report(self, file_path):

        report_df = DataFrame(columns=["model", "genes number", "recall", "precision", "f1"])

        for i, model in enumerate(self.models_to_be_assessed):
            genes = [gene.id.lower() for gene in self.models_to_be_assessed[model].model.genes]
            true_positives, false_positives, false_negatives = \
                self.model_assessor.get_confusion_matrix_cells_for_genes(genes)

            precision = self.model_assessor.get_precision(true_positives, false_positives)
            recall = self.model_assessor.get_recall(true_positives, false_negatives)
            f1_score = self.model_assessor.get_f1_score(recall, precision)
            ratio = self.model_assessor.get_coverage_additional_information_ratio(true_positives, false_positives)
            jaccard_distance = self.model_assessor.get_jaccard_distance(true_positives, false_positives,
                                                                        false_negatives)

            report_df.at[i, "model"] = model
            report_df.at[i, "genes number"] = len(genes)
            report_df.at[i, "true positives"] = true_positives
            report_df.at[i, "false positives"] = false_positives
            report_df.at[i, "false negatives"] = false_negatives
            report_df.at[i, "recall"] = recall
            report_df.at[i, "precision"] = precision
            report_df.at[i, "f1"] = f1_score
            report_df.at[i, "ratio"] = ratio
            report_df.at[i, "jaccard distance"] = jaccard_distance

        report_df.to_csv(file_path, index=False)

    def radar_factory(self, num_vars, frame='circle'):
        """Create a radar chart with `num_vars` axes.

        This function creates a RadarAxes projection and registers it.

        Parameters
        ----------
        num_vars : int
            Number of variables for radar chart.
        frame : {'circle' | 'polygon'}
            Shape of frame surrounding axes.

        """
        # calculate evenly-spaced axis angles
        theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)

        class RadarAxes(PolarAxes):

            name = 'radar'

            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                # rotate plot such that the first axis is at the top
                self.set_theta_zero_location('N')

            def fill(self, *args, closed=True, **kwargs):
                """Override fill so that line is closed by default"""
                return super().fill(closed=closed, *args, **kwargs)

            def plot(self, *args, **kwargs):
                """Override plot so that line is closed by default"""
                lines = super().plot(*args, **kwargs)
                for line in lines:
                    self._close_line(line)

                return lines

            def _close_line(self, line):
                x, y = line.get_data()
                # FIXME: markers at x[0], y[0] get doubled-up
                if x[0] != x[-1]:
                    x = np.concatenate((x, [x[0]]))
                    y = np.concatenate((y, [y[0]]))
                    line.set_data(x, y)

            def set_varlabels(self, labels):
                self.set_thetagrids(np.degrees(theta), labels)

            def _gen_axes_patch(self):
                # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
                # in axes coordinates.
                if frame == 'circle':
                    return Circle((0.5, 0.5), 0.5)
                elif frame == 'polygon':
                    return RegularPolygon((0.5, 0.5), num_vars,
                                          radius=.5, edgecolor="k")
                else:
                    raise ValueError("unknown value for 'frame': %s" % frame)

            def draw(self, renderer):
                """ Draw. If frame is polygon, make gridlines polygon-shaped """
                if frame == 'polygon':
                    gridlines = self.yaxis.get_gridlines()
                    for gl in gridlines:
                        gl.get_path()._interpolation_steps = num_vars
                super().draw(renderer)

            def _gen_axes_spines(self):
                if frame == 'circle':
                    return super()._gen_axes_spines()
                elif frame == 'polygon':
                    # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                    spine = Spine(axes=self,
                                  spine_type='circle',
                                  path=Path.unit_regular_polygon(num_vars))
                    # unit_regular_polygon gives a polygon of radius 1 centered at
                    # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                    # 0.5) in axes coordinates.
                    spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                        + self.transAxes)

                    return {'polar': spine}
                else:
                    raise ValueError("unknown value for 'frame': %s" % frame)

        register_projection(RadarAxes)
        return theta

    def generate_radar_graph(self, metrics, metrics_names, output_file):
        results = []
        for i, model in enumerate(self.models_to_be_assessed):
            metrics_results = self.model_assessor.compare_model(model, metrics)
            metrics_results_list = []
            for metric in metrics_results:
                metrics_results_list.append(metrics_results[metric])

            results.append(metrics_results_list)

        data = [metrics_names,
                ('Approaches performance', results)]

        N = len(data[0])
        theta = self.radar_factory(N, frame='circle')

        spoke_labels = data.pop(0)
        title, case_data = data[0]

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='radar'))
        fig.subplots_adjust(top=0.85, bottom=0.05)

        ax.set_rgrids([0.2, 0.4, 0.6, 0.8])
        ax.set_title(title, position=(0.5, 1.1), ha='center')
        lax = []
        for i, d in enumerate(case_data):
            line = ax.plot(theta, d)
            ax.fill(theta, d, alpha=0.25)
            lax.append(line[0])

        index = [model.id for model in self.models_to_be_assessed]
        ax.legend(handles=lax, labels=index, loc='lower right')
        ax.set_varlabels(spoke_labels)

        plt.savefig(output_file)
