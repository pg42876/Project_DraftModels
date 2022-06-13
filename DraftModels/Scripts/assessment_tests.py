from unittest import TestCase

from Scripts.core import read_sbml_into_cobra_model
from Scripts.utils import ReconstructionTool
from Scripts.assessment import ResultsReport


class TestAssessment(TestCase):

    def test_merlin_reaction_ids_conversion(self):
        model2_model_seed = read_sbml_into_cobra_model(file_path="../Models/ModelSEED/Bpertussis.xml",
                                                       database_version="modelseed",
                                                       reconstruction_tool=ReconstructionTool.MODELSEED.value)

        model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/Bpertussis.xml",
                                           database_version="bigg", reconstruction_tool="curated")

        models_to_assess = {"model_seed": model2_model_seed}
        report = ResultsReport(model, models_to_assess, conversion_method="metanetx", reference_model_format="bigg")
        report.generate_reactions_report("coiso.csv")



    def test_genes_analysis(self):
        merlin_model = read_sbml_into_cobra_model(file_path="../Models/Merlin/Bpertussis.xml",
                                                  database_version="kegg",
                                                  reconstruction_tool=ReconstructionTool.MERLIN.value)

        models_to_assess = {"merlin": merlin_model}

        model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/Bpertussis.xml",
                                           database_version="bigg", reconstruction_tool="curated")

        report = ResultsReport(model, models_to_assess)
        report.generate_genes_report("/results/genes.csv")

    def test_genes_l_plantarum(self):
        from Scripts.core import read_sbml_into_cobra_model
        from Scripts.utils import ReconstructionTool

        organism = "Lplantarum"
        merlin_model = read_sbml_into_cobra_model(file_path="../Models/Merlin/" + str(organism) + ".xml",
                                                  database_version="kegg",
                                                  reconstruction_tool=ReconstructionTool.MERLIN.value)
        autokeggrec_model = read_sbml_into_cobra_model(file_path="../Models/AutoKEGGRec/" + str(organism) + ".xml",
                                                       database_version="kegg",
                                                       reconstruction_tool=ReconstructionTool.AUTOKEGGREC.value)

        carveme_model = read_sbml_into_cobra_model(file_path="../Models/CarveMe/" + str(organism) + ".xml",
                                                   database_version="bigg",
                                                   reconstruction_tool=ReconstructionTool.CARVEME.value)

        modelseed_model = read_sbml_into_cobra_model(file_path="../Models/ModelSEED/" + str(organism) + ".xml",
                                                     database_version="modelseed",
                                                     reconstruction_tool=ReconstructionTool.MODELSEED.value)

        pathwaytools_model = read_sbml_into_cobra_model(file_path="../Models/PathwayTools/" + str(organism) + ".xml",
                                                        database_version="metacyc",
                                                        reconstruction_tool=ReconstructionTool.PATHWAYTOOLS.value)

        raven_model = read_sbml_into_cobra_model(file_path="../Models/RAVEN/" + str(organism) + ".xml",
                                                 database_version="kegg",
                                                 reconstruction_tool=ReconstructionTool.RAVEN.value)

        merlin_bit_model = read_sbml_into_cobra_model(file_path="../Models/Merlin-BIT/" + str(organism) + ".xml",
                                                      database_version="bigg",
                                                      reconstruction_tool=ReconstructionTool.MERLIN.value)
        aureme_model = read_sbml_into_cobra_model(file_path="../Models/AuReMe/" + str(organism) + ".xml",
                                                  database_version="bigg",
                                                  reconstruction_tool=ReconstructionTool.AUREME.value)

        models_to_assess = {"merlin": merlin_model,
                            "aureme": aureme_model,
                            "autokeggrec model": autokeggrec_model,
                            "carveme model": carveme_model,
                            "modelseed model": modelseed_model,
                            "pathwaytools model": pathwaytools_model,
                            "raven model": raven_model,
                            "merlin bit model": merlin_bit_model,
                            }

        model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/Lplantarum.xml",
                                           database_version="bigg", reconstruction_tool="curated")
        report = ResultsReport(model, models_to_assess)
        report.generate_genes_report("results/Lplantarum_genes_results.csv")

    def test_tgondii_genes(self):
        organism = "Tgondii"
        merlin_model = read_sbml_into_cobra_model(file_path="../Models/Merlin/" + str(organism) + ".xml",
                                                  database_version="kegg",
                                                  reconstruction_tool=ReconstructionTool.MERLIN.value)
        autokeggrec_model = read_sbml_into_cobra_model(file_path="../Models/AutoKEGGRec/" + str(organism) + ".xml",
                                                       database_version="kegg",
                                                       reconstruction_tool=ReconstructionTool.AUTOKEGGREC.value)

        modelseed_model = read_sbml_into_cobra_model(file_path="../Models/ModelSEED/" + str(organism) + ".xml",
                                                     database_version="modelseed",
                                                     reconstruction_tool=ReconstructionTool.MODELSEED.value)

        pathwaytools_model = read_sbml_into_cobra_model(file_path="../Models/PathwayTools/" + str(organism) + ".xml",
                                                        database_version="metacyc",
                                                        reconstruction_tool=ReconstructionTool.PATHWAYTOOLS.value)

        raven_model = read_sbml_into_cobra_model(file_path="../Models/RAVEN/" + str(organism) + ".xml",
                                                 database_version="kegg",
                                                 reconstruction_tool=ReconstructionTool.RAVEN.value)

        merlin_bit_model = read_sbml_into_cobra_model(file_path="../Models/Merlin-BIT/" + str(organism) + ".xml",
                                                      database_version="bigg",
                                                      reconstruction_tool=ReconstructionTool.MERLIN.value)
        aureme_model = read_sbml_into_cobra_model(file_path="../Models/AuReMe/" + str(organism) + ".xml",
                                                  database_version="bigg",
                                                  reconstruction_tool=ReconstructionTool.AUREME.value)

        models_to_assess = {"merlin": merlin_model,
                            "aureme": aureme_model,
                            "autokeggrec model": autokeggrec_model,
                            "modelseed model": modelseed_model,
                            "pathwaytools model": pathwaytools_model,
                            "raven model": raven_model,
                            "merlin bit model": merlin_bit_model,
                            }

        model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/Tgondii.xml",
                                           database_version="kegg",
                                           reconstruction_tool=ReconstructionTool.T_GONDII_CURATED.value)

        report = ResultsReport(model, models_to_assess)

        report.generate_genes_report("results/Tgondii_genes_results.csv")

    def test_tgondii_reactions(self):
        model = read_sbml_into_cobra_model(file_path="../Models/Manually_curated/Tgondii.xml",
                                           database_version="kegg",
                                           reconstruction_tool=ReconstructionTool.T_GONDII_CURATED.value)

        model2_model_seed = read_sbml_into_cobra_model(file_path="../Models/ModelSEED/Tgondii.xml",
                                                       database_version="modelseed",
                                                       reconstruction_tool=ReconstructionTool.MODELSEED.value)


        models_to_assess = {"model_seed": model2_model_seed}

        print("Number of genes: ", len(model.model.genes))
        report_metanetx_all = ResultsReport(model, models_to_assess, conversion_method="metanetx",
                                            reference_model_format="kegg")
