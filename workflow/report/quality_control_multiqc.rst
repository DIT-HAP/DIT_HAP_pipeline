Quality Control MultiQC Report
==============================

This comprehensive MultiQC report aggregates quality control metrics from multiple stages of the transposon insertion sequencing pipeline:

**Preprocessing Quality Metrics:**

* **FastP Reports**: Read quality, adapter trimming, and filtering statistics
* **Demultiplexing Reports**: Sample assignment accuracy and barcode quality
* **FastQC Reports**: Per-base quality scores, GC content, and sequence composition

**Mapping Quality Metrics:**

* **Samtools Statistics**: Alignment rates, insert sizes, and mapping quality distributions
* **Flagstat Reports**: Read pair statistics and alignment classifications
* **Index Statistics**: Per-chromosome mapping coverage

**Key Quality Indicators:**

* Read quality scores should be >Q30 for most bases
* Adapter contamination should be minimal (<5%)
* Mapping rates should be >80% for good libraries
* Insert size distributions should be consistent across samples

**Interpretation Guidelines:**

* Red flags: Low quality scores, high adapter content, poor mapping rates
* Yellow warnings: Moderate quality issues that may affect downstream analysis
* Green indicators: High-quality data suitable for analysis

This report enables rapid identification of quality issues and comparison of metrics across all samples and timepoints in the experiment.