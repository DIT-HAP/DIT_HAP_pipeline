"""
Enrichment functions.
"""

# ================================= Imports =================================

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Literal
import time
import requests
from requests.exceptions import RequestException
from io import StringIO
import altair as alt
from datetime import date
from dataclasses import dataclass

from .utils import read_file

# == GOATools imports ==
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt

# ================================= Constants =================================

STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = Literal["tsv", "tsv-no-header", "json", "xml"]
STRING_METHOD = Literal["get_string_ids", "enrichment"]
STRING_SPECIES_ID = "4896"
STRING_CALLER_IDENTITY = "dit-hap.streamlit.app"
STRING_SEPARATOR = "%0d"
MAX_RETRIES = 5
RETRY_DELAY = 10


# ================================= Data Configuration Classes =================================
@dataclass
class OntologyData:
    """Container for ontology data."""

    ontology_obo_path: Path
    ontology_association_file: Path
    slim_term_dataframe: pd.DataFrame


@dataclass
class OntologyDataConfig:
    """Configuration for ontology data file paths."""

    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms_table: list[Path]

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.ontology_obo,
            self.ontology_association_gaf,
            *self.slim_terms_table,
        ]

        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene ontology file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    def load_data(self) -> OntologyData:
        """Load ontology data from files and return as OntologyData dataclass."""
        self.validate_paths()

        slim_dfs = [
            read_file(path, **{"header": None, "names": ["Term", "Description"]})
            for path in self.slim_terms_table
        ]
        slim_df = pd.concat(slim_dfs, ignore_index=True)

        return OntologyData(
            ontology_obo_path=self.ontology_obo,
            ontology_association_file=self.ontology_association_gaf,
            slim_term_dataframe=slim_df,
        )


@dataclass
class GeneMetaData:
    """Container for gene metadata data."""

    gene_info_with_essentiality: pd.DataFrame
    id2name: dict


@dataclass
class GeneMetaConfig:
    """Configuration for gene metadata file paths."""

    gene_IDs_names_products: Path
    deletion_library_essentiality: Path

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.gene_IDs_names_products,
            self.deletion_library_essentiality,
        ]

        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene metadata file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    def load_data(self) -> GeneMetaData:
        """Load gene metadata data from configuration."""

        gene_IDs_names_products = read_file(self.gene_IDs_names_products)
        deletion_library_essentiality = read_file(self.deletion_library_essentiality)

        # Fill NA in gene_name with gene_systematic_id
        gene_IDs_names_products["gene_name"] = gene_IDs_names_products[
            "gene_name"
        ].fillna(gene_IDs_names_products["gene_systematic_id"])
        id2name = dict(
            zip(
                list(gene_IDs_names_products["gene_systematic_id"]),
                list(gene_IDs_names_products["gene_name"]),
            )
        )

        # Merge gene_IDs_names_products and deletion_library_essentiality
        gene_info_with_essentiality = gene_IDs_names_products.merge(
            deletion_library_essentiality[
                [
                    "Updated_Systematic_ID",
                    "Gene dispensability. This study",
                    "Deletion mutant phenotype description",
                    "Phenotypic classification used for analysis",
                    "Category",
                ]
            ],
            how="left",
            left_on="gene_systematic_id",
            right_on="Updated_Systematic_ID",
        ).drop(columns=["Updated_Systematic_ID"])

        return GeneMetaData(
            gene_info_with_essentiality=gene_info_with_essentiality, id2name=id2name
        )


# ================================= Helper Functions =================================


def assign_term_name(term_ID, term_dag):
    if term_ID in term_dag:
        return term_dag[term_ID].name
    else:
        return "No record for {}".format(term_ID)


def format_phaf_file(fypo_obo_file: Path, phaf_file: Path) -> Path:
    """Format the phaf file to the go style gaf file."""
    phaf_dag = GODag(str(fypo_obo_file))
    phaf = pd.read_csv(phaf_file, sep="\t").query(
        "(`Allele type` == 'deletion' or `Allele type` == 'disruption') and Condition.str.contains('FYECO:0000005')"
    )
    phaf["DB"] = "PomBase"
    phaf["DB_Object_ID"] = phaf["Gene systematic ID"]
    phaf["DB_Object_Symbol"] = phaf["Gene symbol"]
    phaf["Qualifier"] = ""
    phaf["GO_ID"] = phaf["FYPO ID"]
    phaf["DB:Reference"] = phaf["Reference"]
    phaf["Evidence"] = phaf["Evidence"]
    phaf["With"] = ""
    phaf["Aspect"] = "FYPO"
    phaf["DB_Object_Name"] = phaf["FYPO ID"].apply(assign_term_name, term_dag=phaf_dag)
    phaf["Synonym"] = ""
    phaf["DB_Object_Type"] = "protein"
    phaf["Taxon"] = "taxon:4896"
    phaf["Date"] = phaf["Date"].str.replace("-", "")
    phaf["Assigned_By"] = phaf["#Database name"]
    phaf["Annotation_Extension"] = phaf["Extension"]
    phaf["Gene_Product_Form_ID"] = ""
    reformat_phaf = phaf[
        [
            "DB",
            "DB_Object_ID",
            "DB_Object_Symbol",
            "Qualifier",
            "GO_ID",
            "DB:Reference",
            "Evidence",
            "With",
            "Aspect",
            "DB_Object_Name",
            "Synonym",
            "DB_Object_Type",
            "Taxon",
            "Date",
            "Assigned_By",
            "Annotation_Extension",
            "Gene_Product_Form_ID",
        ]
    ].copy()

    with open(phaf_file.parent / "phaf_go_style.tsv", "w") as f:
        f.write(
            f"!gaf-version: 2.2\n!generated-by: Yusheng Yang\n!date-generated: {date.today().strftime('%Y-%m-%d')}\n!URL: https://www.pombase.org/monthly_releases/2025/pombase-2025-09-01/phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv\n!contact: yangyusheng@nibs.ac.cn\n"
        )

    reformat_phaf.to_csv(
        phaf_file.parent / "phaf_go_style.tsv",
        sep="\t",
        index=False,
        header=False,
        mode="a",
    )

    return phaf_file.parent / "phaf_go_style.tsv"


def format_mondo_gaf_file(mondo_obo_file: Path, mondo_gaf_file: Path) -> Path:
    """Format the mondo gaf file to the go style gaf file."""
    mondo_dag = GODag(str(mondo_obo_file))
    mondo = pd.read_csv(mondo_gaf_file, sep="\t")
    mondo["DB"] = "Pombase"
    mondo["DB_Object_ID"] = mondo["#gene_systematic_id"]
    mondo["DB_Object_Symbol"] = mondo["gene_name"]
    mondo["Qualifier"] = ""
    mondo["GO_ID"] = mondo["mondo_id"]
    mondo["DB:Reference"] = mondo["reference"]
    mondo["Evidence"] = ""
    mondo["With"] = ""
    mondo["Aspect"] = "MONDO"
    mondo["DB_Object_Name"] = mondo["mondo_id"].apply(
        assign_term_name, term_dag=mondo_dag
    )
    mondo["Synonym"] = ""
    mondo["DB_Object_Type"] = "protein"
    mondo["Taxon"] = "taxon:4896"
    mondo["Date"] = mondo["date"].fillna("2025-09-01").str.replace("-", "")
    mondo["Assigned_By"] = "PomBase"
    mondo["Annotation_Extension"] = ""
    mondo["Gene_Product_Form_ID"] = ""
    reformat_mondo = mondo[
        [
            "DB",
            "DB_Object_ID",
            "DB_Object_Symbol",
            "Qualifier",
            "GO_ID",
            "DB:Reference",
            "Evidence",
            "With",
            "Aspect",
            "DB_Object_Name",
            "Synonym",
            "DB_Object_Type",
            "Taxon",
            "Date",
            "Assigned_By",
            "Annotation_Extension",
            "Gene_Product_Form_ID",
        ]
    ].copy()

    with open(mondo_gaf_file.parent / "mondo_go_style.tsv", "w") as f:
        f.write(
            f"!gaf-version: 2.2\n!generated-by: Yusheng Yang\n!date-generated: {date.today().strftime('%Y-%m-%d')}\n!URL: https://www.pombase.org/monthly_releases/2025/pombase-2025-08-01/ontologies_and_associations/human_disease_association.tsv\n!contact: yangyusheng@nibs.ac.cn\n"
        )

    reformat_mondo.to_csv(
        mondo_gaf_file.parent / "mondo_go_style.tsv",
        sep="\t",
        index=False,
        header=False,
        mode="a",
    )

    return mondo_gaf_file.parent / "mondo_go_style.tsv"


def mapslim(term: str, dag: GODag, slim_dag: dict) -> tuple[set[str], set[str]]:
    """Maps a term (accession) to it's slim terms."""

    all_ancestors = set()
    covered_ancestors = set()

    # get all paths for the term in the dag
    paths = dag.paths_to_top(term)
    for path in paths:
        # the next loop needs to run bottom->up, i.e. from the term item to
        # the root, thus we need to reverse the list prior to iteration
        path.reverse()

        got_leaf = False
        for term in path:
            if term.id in slim_dag:
                all_ancestors.add(term.id)
                if got_leaf:
                    covered_ancestors.add(term.id)
                got_leaf = True

    # get the direct ancestors, i.e. those that are not covered by a earlier
    # ancestor of the slim in _any_ path (in bottom->top order)
    direct_ancestors = all_ancestors - covered_ancestors
    return direct_ancestors, all_ancestors


def get_slim_ns2assoc(ns2assoc: dict, dag: GODag, slim_dag: dict) -> dict:
    """Get the slim ns2assoc."""

    term2slim = {}
    for term in dag:
        term2slim[term] = {}
        term2slim[term]["direct_ancestors"], term2slim[term]["all_ancestors"] = mapslim(
            term, dag, slim_dag
        )

    ns2slim_assoc = {"direct_ancestors": {}, "all_ancestors": {}}
    for ns, gene2terms in ns2assoc.items():
        ns2slim_assoc["direct_ancestors"][ns] = {}
        ns2slim_assoc["all_ancestors"][ns] = {}
        for gene, terms in gene2terms.items():
            ns2slim_assoc["direct_ancestors"][ns][gene] = set()
            ns2slim_assoc["all_ancestors"][ns][gene] = set()
            for term in terms:
                ns2slim_assoc["direct_ancestors"][ns][gene].update(
                    term2slim[term]["direct_ancestors"]
                )
                ns2slim_assoc["all_ancestors"][ns][gene].update(
                    term2slim[term]["all_ancestors"]
                )
    return ns2slim_assoc


def create_enrichment_dataframe(
    oea_results_sig: list[GOEnrichmentRecord], **kwargs
) -> pd.DataFrame:
    """Create an enrichment dataframe from the GOEnrichmentRecord objects manually."""
    results_list = []
    for result in oea_results_sig:
        res_dict = {
            "GO": result.GO,
            "NS": result.NS,
            "name": result.name,
            "level": result.goterm.level,
            "depth": result.goterm.depth,
            "p_uncorrected": result.p_uncorrected,
            "p_fdr_bh": result.p_fdr_bh,
            "study_count": result.study_count,
            "study_n": result.study_n,
            "pop_count": result.pop_count,
            "pop_n": result.pop_n,
            "ratio_in_study": "/".join(map(str, result.ratio_in_study)),
            "ratio_in_pop": "/".join(map(str, result.ratio_in_pop)),
            "study_items": result.study_items,
            "pop_items": result.pop_items,
        }

        # get the definition, but handle the case where it is not available
        try:
            res_dict["defn"] = result.goterm.defn
        except Exception as e:
            res_dict["defn"] = ""

        # convert the study and pop items to names, but handle the case where the itemid2name is not available
        if kwargs["itemid2name"] is not None:
            study_items = sorted([
                kwargs["itemid2name"][item] for item in result.study_items
            ])
            pop_items = sorted([
                kwargs["itemid2name"][item] for item in result.pop_items
            ])
            res_dict["study_items"] = ", ".join(study_items)
            res_dict["pop_items"] = ", ".join(pop_items)
        else:
            res_dict["study_items"] = ", ".join(sorted(result.study_items))
            res_dict["pop_items"] = ", ".join(sorted(result.pop_items))

        # append the result to the list
        results_list.append(res_dict)

    oea_results_sig_prt = pd.DataFrame(results_list)
    return oea_results_sig_prt


# ================================= Main Functions =================================
def load_ontology_data(
    ontology_data: OntologyData, **kwargs
) -> tuple[GODag, GafReader, dict, dict, dict, dict, dict]:
    """Load ontology data from obo file and association file."""
    try:
        dag = GODag(
            str(ontology_data.ontology_obo_path),
            optional_attrs=["def", "relationship"],
            load_obsolete=False,
        )
    except KeyError:
        dag = GODag(
            str(ontology_data.ontology_obo_path),
            optional_attrs=["def"],
            load_obsolete=False,
        )

    objanno = GafReader(str(ontology_data.ontology_association_file), godag=dag)

    slim_terms = ontology_data.slim_term_dataframe["Term"].to_list()
    slim_dag = {term: dag[term] for term in slim_terms}

    ns2assoc = objanno.get_ns2assc(**kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    gene2go = objanno.get_id2gos_nss(**kwargs)
    go2genes = objanno.get_id2gos_nss(go2geneids=True, **kwargs)
    return dag, objanno, ns2assoc, gene2go, go2genes, slim_dag, ns2slim_assoc


def ontology_enrichment(
    query_genes: list[str], bg_genes: list[str], **kwargs
) -> tuple[GOEnrichmentStudyNS, list[GOEnrichmentRecord]]:
    """Perform ontology enrichment analysis."""
    oea_obj = GOEnrichmentStudyNS(bg_genes, **kwargs)

    # Run enrichment analysis
    oea_results = oea_obj.run_study(query_genes, **kwargs)
    oea_results_sig = [
        r
        for r in oea_results
        if (r.p_fdr_bh < kwargs["alpha"]) and (r.enrichment == "e")
    ]

    return oea_obj, oea_results_sig


def format_ontology_enrichment_results(
    label: str, oea_results_sig: list[GOEnrichmentRecord], **kwargs
) -> pd.DataFrame:
    """Format ontology enrichment results into a pandas DataFrame."""
    # itemid2name to covert itemid to name
    # transform the result to a dataframe
    try:
        oea_results_sig_prt = pd.DataFrame(get_goea_nts_prt(oea_results_sig, **kwargs))
    except Exception:
        oea_results_sig_prt = create_enrichment_dataframe(oea_results_sig, **kwargs)

    if oea_results_sig_prt.empty:
        return pd.DataFrame()
    else:
        # calculate the gene_ratio and term_coverage
        oea_results_sig_prt["gene_ratio"] = round(
            oea_results_sig_prt["study_count"] / oea_results_sig_prt["study_n"], 2
        )
        oea_results_sig_prt["term_coverage"] = round(
            oea_results_sig_prt["study_count"] / oea_results_sig_prt["pop_count"], 2
        )

        # sort gene items
        oea_results_sig_prt["study_items"] = oea_results_sig_prt["study_items"].apply(
            lambda x: ",".join(sorted(x.split(", ")))
        )
        oea_results_sig_prt["pop_items"] = oea_results_sig_prt["pop_items"].apply(
            lambda x: ",".join(sorted(x.split(", ")))
        )

        # keep the columns we need and reorder them
        # prt_columns = ["GO", "NS", "enrichment", "name", "ratio_in_study", "ratio_in_pop", "p_uncorrected", "depth", "study_count", "p_fdr_bh", "study_items", "pop_items", "study_n", "pop_count", "pop_n", "item_id", "namespace", "level", "is_obsolete", "alt_ids", "defn"]
        kept_columns = [
            "NS",
            "GO",
            "name",
            "level",
            "depth",
            "p_uncorrected",
            "p_fdr_bh",
            "gene_ratio",
            "term_coverage",
            "study_count",
            "study_n",
            "pop_count",
            "pop_n",
            "ratio_in_study",
            "ratio_in_pop",
            "study_items",
            "pop_items",
        ]
        if "defn" in oea_results_sig_prt.columns:
            kept_columns.append("defn")
        oea_results_sig_prt = (
            oea_results_sig_prt[kept_columns]
            .copy()
            .rename(
                columns={
                    "NS": "namespace",
                    "GO": "term_id",
                    "name": "term",
                    "p_fdr_bh": "p_fdr",
                }
            )
        )
        return oea_results_sig_prt


def ontology_enrichment_pipeline(
    ontology_data: OntologyData,
    query_genes: list[str],
    bg_genes: list[str],
    load_kwargs: dict = {},
    enrichment_kwargs: dict = {},
    format_kwargs: dict = {},
) -> tuple[pd.DataFrame, pd.DataFrame, GODag, GafReader]:
    """Perform ontology enrichment analysis pipeline."""

    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag, ns2slim_assoc = (
        load_ontology_data(ontology_data, **load_kwargs)
    )

    enrichment_kwargs["godag"] = dag
    enrichment_kwargs["ns2assoc"] = ns2assoc

    slim_enrichment_kwargs = enrichment_kwargs.copy()
    slim_enrichment_kwargs["godag"] = slim_dag
    slim_enrichment_kwargs["ns2assoc"] = ns2slim_assoc["all_ancestors"]
    slim_enrichment_kwargs["propagate_counts"] = False

    oea_obj, oea_results_sig = ontology_enrichment(
        query_genes, bg_genes, **enrichment_kwargs
    )
    oea_obj_slim, oea_results_sig_slim = ontology_enrichment(
        query_genes, bg_genes, **slim_enrichment_kwargs
    )

    oea_results_sig_prt = format_ontology_enrichment_results(
        "full", oea_results_sig, **format_kwargs
    )
    oea_results_sig_prt_slim = format_ontology_enrichment_results(
        "slim", oea_results_sig_slim, **format_kwargs
    )

    return oea_results_sig_prt, oea_results_sig_prt_slim, dag, objanno


def stringdb_api_functions(
    output_format: STRING_OUTPUT_FORMAT = "xml",
    method: STRING_METHOD = "get_string_ids",
    params: dict = {},
) -> pd.DataFrame:
    """Perform STRING API functions."""
    request_url = "/".join([STRING_API_URL, output_format, method])
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(request_url, data=params)
            response.raise_for_status()
            if "Something went wrong!" in response.text:
                raise RequestException(
                    "communication_error: Something went wrong! -- please contact STRING maintainers if the issue persists"
                )
            break
        except RequestException as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
            else:
                raise ValueError(
                    f"Error performing STRING API functions after {MAX_RETRIES} retries: {e}"
                )

    match output_format:
        case "xml":
            return pd.read_xml(StringIO(response.text))
        case "tsv":
            return pd.read_csv(StringIO(response.text), sep="\t")
        case "tsv-no-header":
            return pd.read_csv(StringIO(response.text), sep="\t", header=None)
        case "json":
            return pd.read_json(StringIO(response.text))
        case _:
            raise ValueError(f"Invalid output format: {output_format}")


def format_string_enrichment_results(
    enrichment_df: pd.DataFrame, query_genes: list[str], background_genes: list[str]
) -> pd.DataFrame:
    """Format STRING enrichment results into a pandas DataFrame."""

    renamed_column_names = {
        "term": "term_id",
        "category": "namespace",
        "description": "term",
        "p_value": "p_uncorrected",
        "fdr": "p_fdr",
        "gene_ratio": "gene_ratio",
        "term_coverage": "term_coverage",
        "number_of_genes": "study_count",
        "study_n": "study_n",
        "number_of_genes_in_background": "pop_count",
        "pop_n": "pop_n",
        "ratio_in_study": "ratio_in_study",
        "ratio_in_pop": "ratio_in_pop",
        "preferredNames": "study_items",
    }
    enrichment_df = enrichment_df.rename(columns=renamed_column_names)
    enrichment_df["study_n"] = len(query_genes)
    enrichment_df["pop_n"] = len(background_genes)
    enrichment_df["ratio_in_study"] = enrichment_df.apply(lambda row: f"{row['study_count']}/{row['study_n']}", axis=1)
    enrichment_df["ratio_in_pop"] = enrichment_df.apply(lambda row: f"{row['pop_count']}/{row['pop_n']}", axis=1)
    enrichment_df["gene_ratio"] = round(
        enrichment_df["study_count"] / enrichment_df["study_n"], 2
    )
    enrichment_df["term_coverage"] = round(
        enrichment_df["study_count"] / enrichment_df["pop_count"], 2
    )
    enrichment_df = enrichment_df[list(renamed_column_names.values())]

    namespace_description = {
        "Process": "Biological Process (Gene Ontology)",
        "Component": "Cellular Component (Gene Ontology)",
        "Function": "Molecular Function (Gene Ontology)",
        "PMID": "Reference Publications (PubMed)",
        "NetworkNeighborAL": "Local Network Cluster (STRING)",
        "KEGG": "KEGG Pathways",
        "RCTM": "Reactome Pathways",
        "COMPARTMENTS": "Subcellular Localization (COMPARTMENTS)",
        "Keyword": "Annotated Keywords (UniProt)",
        "InterPro": "Protein Domains and Features (InterPro)",
        "SMART": "Protein Domains and Features (SMART)",
    }
    enrichment_df["namespace"] = enrichment_df["namespace"].map(namespace_description)

    namespace_order = list(namespace_description.values())
    enrichment_df["namespace_order"] = enrichment_df["namespace"].map({
        v: i for i, v in enumerate(namespace_order)
    })
    enrichment_df.sort_values(by="namespace_order", inplace=True)
    enrichment_df.drop("namespace_order", axis=1, inplace=True)

    return enrichment_df


def stringdb_enrichment(query_genes, bg_genes):
    """Perform STRING enrichment analysis using the STRING API."""

    get_string_id_params = {
        "identifiers": STRING_SEPARATOR.join(bg_genes),
        "species": STRING_SPECIES_ID,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": STRING_CALLER_IDENTITY,
    }

    get_string_id_response = stringdb_api_functions(
        output_format="xml", method="get_string_ids", params=get_string_id_params
    )
    get_string_id_response = get_string_id_response["stringId"].tolist()

    enrichment_params = {
        "identifiers": STRING_SEPARATOR.join(query_genes),
        "species": STRING_SPECIES_ID,
        "background_string_identifiers": STRING_SEPARATOR.join(get_string_id_response),
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    # Parse results
    enrichment_df = stringdb_api_functions(
        output_format="xml", method="enrichment", params=enrichment_params
    )

    # format the results
    enrichment_df = format_string_enrichment_results(
        enrichment_df, query_genes, bg_genes
    )

    return enrichment_df


def display_enrichment_results(enrichment_results: pd.DataFrame) -> alt.VConcatChart:
    """Display enrichment results."""
    charts = []
    for ns, ns_results in enrichment_results.groupby("namespace", sort=False):
        chart = (
            alt.Chart(ns_results)
            .mark_circle()
            .encode(
                alt.X("gene_ratio:Q", axis=alt.Axis(grid=True), title="Gene ratio"),
                alt.Y(
                    "term:N",
                    axis=alt.Axis(
                        grid=True, labelLimit=500, title="Term", orient="right"
                    ),
                    sort=alt.EncodingSortField(field="gene_ratio", order="descending"),
                ),
                size=alt.Size("term_coverage:Q", title="Term coverage"),
                color=alt.Color("p_fdr:Q", title="FDR").scale(
                    scheme="yelloworangered", reverse=True
                ),
                tooltip=ns_results.columns.tolist(),
            )
        ).properties(
            title=f"Enrichment results for {ns}",
        )
        charts.append(chart)

    return alt.vconcat(*charts)

def create_customized_enrichment_plot(
    enrichment_df: pd.DataFrame,
    title: str,
    x_col: str = "gene_ratio",
    y_col: str = "name",
    color_col: str = "p_fdr",
    size_col: str = "term_coverage"
):
    """Create an interactive Altair plot for enrichment results."""

    # Create scatter plot
    scatter = (
        alt.Chart(enrichment_df)
        .mark_circle()
        .encode(
            x=alt.X(f"{x_col}:N", title=f"{x_col}", axis=alt.Axis(grid=True)),
            y=alt.Y(
                f"{y_col}:N",
                sort=alt.EncodingSortField(field=x_col, order="ascending"),
                title="Enriched Terms",
                axis=alt.Axis(grid=True),
            ),
            color=alt.Color(
                f"{color_col}:Q", title=f"{color_col}", scale=alt.Scale(scheme="yelloworangered", reverse=True)
            ),
            size=alt.Size(
                f"{size_col}:Q", title=f"{size_col}"
            ),
            tooltip=enrichment_df.columns.tolist()
        )
    ).properties(
        title=title,
    )

    return scatter


def revigo_analysis(enrich_df: pd.DataFrame, cut_off: float = 0.7) -> pd.DataFrame:
    """Perform REVIGO analysis on enrichment results."""
    GOs = enrich_df["GO"].tolist()
    padj = enrich_df["p_fdr_bh"].tolist()
    data = "\n".join([f"{GO}\t{padj}" for GO, padj in zip(GOs, padj)])
    payload = {'cutoff':f'{cut_off}', 'valueType':'pvalue', 'speciesTaxon':'284812', 'measure':'SIMREL', 'goList':data}
    r = requests.post("http://revigo.irb.hr/Revigo", data=payload)
    revigo_res = pd.read_html(r.text)[0]
    id2name = dict(zip(revigo_res["Term ID"], revigo_res["Name"]))
    revigo_res["Representative"] = revigo_res.apply(lambda row: row["Name"] if np.isnan(row["Representative"]) else id2name["GO:" + str(int(row["Representative"])).rjust(7, "0")], axis=1)
    return revigo_res