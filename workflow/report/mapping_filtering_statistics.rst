Mapping and Filtering Statistics Summary
========================================

This table provides comprehensive statistics on read mapping and filtering performance across all samples and timepoints.

**Metrics Included:**

* **Total Reads**: Raw read counts before any filtering
* **Mapped Reads**: Successfully aligned reads to the reference genome
* **Mapping Rate**: Percentage of reads that mapped successfully
* **Properly Paired**: Read pairs with correct orientation and insert size
* **Filtered Reads**: Reads removed during quality filtering steps
* **Final Retained**: High-quality reads retained for downstream analysis

**Quality Thresholds:**

* Mapping rates >80% indicate good library quality
* Proper pairing rates >90% suggest intact DNA fragments
* Filtering rates <20% indicate minimal quality issues

**Sample Comparison:**

This table enables identification of:

* Samples with poor mapping performance
* Timepoints with systematic quality differences
* Batch effects or technical issues
* Overall experiment success rates

**Troubleshooting:**

* Low mapping rates may indicate contamination or poor library preparation
* High filtering rates suggest degraded DNA or sequencing issues
* Inconsistent metrics across replicates indicate technical problems

Use this summary to assess data quality before proceeding with downstream transposon insertion analysis.