import pandas as pd
import argparse
import sys

try:
    parser = argparse.ArgumentParser(
        prog="Annotate Variants from VEP and openCravat",
        description="change format from vcf to tsv and filter varaints.",
    )

    parser.add_argument("--output", type=str, required=True, help="Ouput path.")
    parser.add_argument("--sample", type=str, required=True, help="Sample.")
    parser.add_argument("--vcf", type=str, required=True, help="VCF from VEP.")

    args = parser.parse_args()

    output = args.output
    sample = args.sample
    vcf = args.vcf

    def normalize_headers(h1, h2):
        # propagate h1 values forward
        last = ""
        h1_norm = []
        for x in h1:
            if x != "":
                last = x
            h1_norm.append(last)

        # combine
        cols = [f"{a}_{b}".replace(" ", "_") for a, b in zip(h1_norm, h2)]
        return cols

    def load_tsv_until_next_hash(path):
        with open(path, "r") as f:
            # Skip # metadata line(s)
            for line in f:
                if not line.startswith("#"):
                    header1 = line.rstrip("\n").split("\t")
                    break

            header2 = next(f).rstrip("\n").split("\t")

            # normalize headers so columns match correctly
            columns = normalize_headers(header1, header2)

            expected = len(columns)

            rows = []
            app = rows.append

            for line in f:
                if line.startswith("#"):
                    break
                fields = line.rstrip("\n").split("\t")

                if len(fields) != expected:
                    # Skip malformed row
                    print("malformed line: " + fields)

                app(fields)

        return pd.DataFrame(rows, columns=columns, dtype=str)

    df = load_tsv_until_next_hash(vcf + ".tsv")

    # Select columns of interest
    cols = [
        "Variant_Annotation_Chrom",
        "Variant_Annotation_Position",
        "Variant_Annotation_End_Position",
        "Variant_Annotation_Ref_Base",
        "Variant_Annotation_Alt_Base",
        "VCF_Info_Phred",
        "VCF_Info_Alternate_reads",
        "VCF_Info_Total_reads",
        "VCF_Info_Zygosity",
        "Extra_VCF_INFO_Annotations_CSQ_STRAND",
        "Extra_VCF_INFO_Annotations_CSQ_SYMBOL",
        "Extra_VCF_INFO_Annotations_CSQ_Gene",
        "Extra_VCF_INFO_Annotations_CSQ_Feature",
        "Extra_VCF_INFO_Annotations_CSQ_BIOTYPE",
        "Extra_VCF_INFO_Annotations_CSQ_VARIANT_CLASS",
        "Extra_VCF_INFO_Annotations_CSQ_Consequence",
        "Extra_VCF_INFO_Annotations_CSQ_IMPACT",
        "dbSNP_rsID",
        "1000_Genomes_AF",
        "ALFA:_Allele_Frequency_Aggregator_Global_AF",
        "gnomAD4_Global_AF",
        "gnomAD4_Global_AN",
        "gnomAD4_Global_AC",
        "gnomAD4_Non-Fin_Eur_AF",
        "gnomAD4_Non-Fin_Eur_AN",
        "gnomAD4_Non-Fin_Eur_AC",
        "CADD_Phred",
        "PolyPhen-2_HDIV_Rank_Score",
        "PolyPhen-2_HVAR_Rank_Score",
        "REVEL_Rank_Score",
        "SIFT_Rank_Score",
        "SIFT_Prediction",
        "SIFT_Confidence",
        "AlphaMissense_Protein",
        "AlphaMissense_Score",
        "AlphaMissense_Class",
        "GERP++_RS_Score",
        "Phast_Cons_Vert_Ranked_Score",
        "PhyloP_Vert_Ranked_Score",
        "BayesDel_Prediction_(AF)",
        "Mutation_Assessor_Functional_Impact",
        "MutationTaster_Prediction",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_UniProt_Accession_Number",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_Prediction",
        "ALoFT_Classification",
        "ALoFT_Confidence",
        "GeneHancer_GeneHancer_Type",
        "GeneHancer_Gene_Targets",
        "MutPred-Indel_Rank_score",
        "MutPred-Indel_Property",
        "ClinGen_Gene_Disease",
        "ClinVar_Clinical_Significance",
        "GTEx_eQTLs_Target_Gene",
        "GWAS_Catalog_Disease/Trait",
        "GWAS_Catalog_Odds_Ratio/Beta_Coeff",
        "GWAS_Catalog_P-value",
        "GWAS_Catalog_PMID",
        "OMIM_Entry_ID",
        "Pangolin_Splicing_Gain_Δscore",
        "Pangolin_Splicing_Loss_Δscore",
        "SpliceAI_Acceptor_Gain_Score",
        "SpliceAI_Acceptor_Loss_Score",
        "SpliceAI_Donor_Gain_Score",
        "SpliceAI_Donor_Loss_Score",
        "MITOMAP_Disease",
        "MITOMAP_MitoTip_Score",
        "MITOMAP_PubMed_ID",
        "Variant_Annotation_Samples",
    ]
    df_cols = df[cols]

    column_mapping = {
        "Variant_Annotation_Chrom": "Chromosome",
        "Variant_Annotation_Position": "Start",
        "Variant_Annotation_End_Position": "End",
        "Variant_Annotation_Ref_Base": "ref",
        "Variant_Annotation_Alt_Base": "alt",
        "VCF_Info_Phred": "phred",
        "VCF_Info_Alternate_reads": "alt_reads",
        "VCF_Info_Total_reads": "total_reads",
        "VCF_Info_Zygosity": "zygosity",
        "Extra_VCF_INFO_Annotations_CSQ_STRAND": "strand",
        "Extra_VCF_INFO_Annotations_CSQ_SYMBOL": "gene_symbol",
        "Extra_VCF_INFO_Annotations_CSQ_Gene": "gene_id",
        "Extra_VCF_INFO_Annotations_CSQ_Feature": "transcript",
        "Extra_VCF_INFO_Annotations_CSQ_BIOTYPE": "biotype",
        "Extra_VCF_INFO_Annotations_CSQ_VARIANT_CLASS": "class",
        "Extra_VCF_INFO_Annotations_CSQ_Consequence": "consequence",
        "Extra_VCF_INFO_Annotations_CSQ_IMPACT": "impact",
        "dbSNP_rsID": "rsID",
        "1000_Genomes_AF": "1000_Genomes_AF",
        "ALFA:_Allele_Frequency_Aggregator_Global_AF": "ALFA_AF",
        "gnomAD4_Global_AF": "gnomAD4_global_AF",
        "gnomAD4_Global_AC": "gnomAD4_global_AC",
        "gnomAD4_Global_AN": "gnomAD4_global_AN",
        "gnomAD4_Non-Fin_Eur_AF": "gnomAD4_NFE_AF",
        "gnomAD4_Non-Fin_Eur_AC": "gnomAD4_NFE_AC",
        "gnomAD4_Non-Fin_Eur_AN": "gnomAD4_NFE_AN",
        "CADD_Phred": "CADD_score",
        "PolyPhen-2_HDIV_Rank_Score": "polyPhen2_HDIV_score",
        "PolyPhen-2_HVAR_Rank_Score": "polyPhen2_HVAR_score",
        "REVEL_Rank_Score": "REVEL_score",
        "SIFT_Rank_Score": "SIFT_score",
        "SIFT_Prediction": "SIFT_prediction",
        "SIFT_Confidence": "SIFT_confidence",
        "AlphaMissense_Protein": "AlphaMissense_protein",
        "AlphaMissense_Score": "AlphaMissense_score",
        "AlphaMissense_Class": "AlphaMissense_classification",
        "GERP++_RS_Score": "GERP++_score",
        "Phast_Cons_Vert_Ranked_Score": "PhastCons_score",
        "PhyloP_Vert_Ranked_Score": "PhyloP_score",
        "BayesDel_Prediction_(AF)": "BayesDel_prediction",
        "Mutation_Assessor_Functional_Impact": "MutationAssessor_functional_impact",
        "MutationTaster_Prediction": "MutationTaster_prediction",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_UniProt_Accession_Number": "PROVEAN_uniprot",
        "PROVEAN:_Protein_Variant_Effect_Analyzer_Prediction": "PROVEAN_prediction",
        "ALoFT_Classification": "ALoFT_classification",
        "ALoFT_Confidence": "ALoFT_confidence",
        "GeneHancer_GeneHancer_Type": "GeneHancer_type",
        "GeneHancer_Gene_Targets": "GeneHancer_gene_targets",
        "MutPred-Indel_Rank_score": "MutPred-Indel_score",
        "MutPred-Indel_Property": "MutPred-Indel_proprety",
        "ClinGen_Gene_Disease": "ClinGen_gene_disease",
        "ClinVar_Clinical_Significance": "ClinVar_clinical_significance",
        "GTEx_eQTLs_Target_Gene": "GTEx_eQTLs_target_genes",
        "GWAS_Catalog_Disease/Trait": "GWAS_Catalog_trait",
        "GWAS_Catalog_Odds_Ratio/Beta_Coeff": "GWAS_Catalog_odds_ratio",
        "GWAS_Catalog_P-value": "GWAS_Catalog_pvalue",
        "GWAS_Catalog_PMID": "GWAS_Catalog_PMID",
        "OMIM_Entry_ID": "OMIM",
        "Pangolin_Splicing_Gain_Δscore": "Pangolin_splicing_gain_deltascore",
        "Pangolin_Splicing_Loss_Δscore": "Pangolin_splicing_loss_deltascore",
        "SpliceAI_Acceptor_Gain_Score": "SpliceAI_acceptor_gain_score",
        "SpliceAI_Acceptor_Loss_Score": "SpliceAI_acceptor_loss_score",
        "SpliceAI_Donor_Gain_Score": "SpliceAI_donor_gain_score",
        "SpliceAI_Donor_Loss_Score": "SpliceAI_donor_loss_score",
        "MITOMAP_Disease": "MITOMAP_disease",
        "MITOMAP_MitoTip_Score": "MITOMAP_MitoTip_score",
        "MITOMAP_PubMed_ID": "MITOMAP_PubMed_id",
        "Variant_Annotation_Samples": "sample",
    }

    df_small = df_cols.rename(columns=column_mapping)

    # All
    df_small.to_csv(
        output + sample + "_all.tsv", sep="\t", index=False, lineterminator="\n"
    )

    # Change to numeric
    cols_to_fix = [
        "alt_reads",
        "total_reads",
        "phred",
        "1000_Genomes_AF",
        "ALFA_AF",
        "gnomAD4_global_AF",
        "CADD_score",
        "GERP++_score",
        "PhastCons_score",
        "PhyloP_score",
        "polyPhen2_HDIV_score",
        "polyPhen2_HVAR_score",
        "SpliceAI_acceptor_gain_score",
        "SpliceAI_acceptor_loss_score",
        "SpliceAI_donor_gain_score",
        "SpliceAI_donor_loss_score",
    ]
    df_small[cols_to_fix] = df_small[cols_to_fix].apply(pd.to_numeric, errors="coerce")

    # Rare variants
    df_rare = df_small[
        (df_small["alt_reads"] >= 2)
        & (df_small["total_reads"] >= 5)
        & (df_small["phred"] >= 20)
        & (
            (df_small["1000_Genomes_AF"] <= 0.05)
            | (df_small["ALFA_AF"] <= 0.05)
            | (df_small["gnomAD4_global_AF"] <= 0.05)
        )
    ]

    print('>>> Filters "Rare":')
    print("\t>>> Quality Phred >= 20")
    print("\t>>> Number of alt reads >= 2")
    print("\t>>> Allele Frequency <= 5%")

    df_rare.to_csv(
        output + sample + "_rare.tsv", sep="\t", index=False, lineterminator="\n"
    )

    # Likely Pathogenic variants
    df_patho = df_rare[
        (
            ((df_rare["impact"] == "MODERATE") | (df_rare["impact"] == "HIGH"))
            | (df_rare["CADD_score"] >= 10)
            | (df_rare["GERP++_score"] >= 2)
            | (df_rare["PhastCons_score"] >= 0.5)
            | (df_rare["PhyloP_score"] >= 2.27)
            | (
                (df_rare["MutationTaster_prediction"] == "Automatic Disease Causing")
                | (df_rare["MutationTaster_prediction"] == "Damaging")
            )
            | (
                (df_rare["polyPhen2_HDIV_score"] >= 0.5)
                | (df_rare["polyPhen2_HVAR_score"] >= 0.5)
            )
            | (
                (df_rare["SIFT_score"] == "Damaging")
                & (df_rare["SIFT_confidence"] == "High")
            )
            | (
                (df_rare["SpliceAI_acceptor_gain_score"] >= 0.8)
                | (df_rare["SpliceAI_acceptor_loss_score"] >= 0.8)
                | (df_rare["SpliceAI_donor_gain_score"] >= 0.8)
                | (df_rare["SpliceAI_donor_loss_score"] >= 0.8)
            )
        )
    ]

    print('\n>>> Filters "Pathogenic":')
    print('\t>>> Same as "Rare" and :')
    print('\t>>> Variant impact is "Moderate" or "High"')
    print("\t>>> or CADD (pathogenecity) >= 10")
    print("\t>>> or GERP++ (conservation) >= 2")
    print("\t>>> or PhastCons (conservation) >= 0.5")
    print("\t>>> or PhyloP (conservation) >= 2.27")
    print(
        '\t>>> or MutationTaster  (function prediction) is "Damaging" or "Disease Causing"'
    )
    print("\t>>> or PolyPhen-2 (function prediction) HDIV or HVAR >= 0.5")
    print('\t>>> or SIFT (function prediction) is "Damaging" with high confidence')
    print("\t>>> or SpliceAI (splicing) one of the four scores >= 0.8")

    if len(df_patho) > 0:
        df_patho.to_csv(
            output + sample + "_pathogenic.tsv",
            sep="\t",
            index=False,
            lineterminator="\n",
        )

    # Get number of variants
    all_var = len(df_small.index)
    rare_var = len(df_rare.index)
    patho_var = len(df_patho.index)

    print(
        "\n>>> Number of variants in all: ",
        all_var,
        "\n>>> Number of variants in rare: ",
        rare_var,
        "\n>>> Number of variants in pathogenic: ",
        patho_var,
    )

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
