"""Canonical schema constants for the Kluyveromyces BED dataset."""
IDENTIFIER_COLUMNS = ["experiment_id","source_id","project_id","run_id","species","strain","medium","record_type"]
PROVENANCE_COLUMNS = ["project_name","source_reference","data_type","notes"]
RESPONSE_COLUMNS = [
    "biomass_g_L","OD600","CFU_mL","mu_max_1_h","doubling_time_h",
    "productivity_g_L_h","ethanol_g_L","SCP_g_L","EPS_mg_L",
    "protein_content_pct","residual_protein_g_L","Yxs_lactose_g_g","Cell_yield_g_g"
]
RULES = {
    "grain":"one row = one experimental run/condition",
    "missing":"missing is unknown/unreported, not zero",
    "raw":"never overwrite files under data/raw/"
}
