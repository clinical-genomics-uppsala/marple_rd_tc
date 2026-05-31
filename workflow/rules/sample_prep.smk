__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2025, Jessika Nordin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule sample_prep:
    input:
        samples="samples.tsv",
    output:
        sex=temp("sample_sex.tsv"),
    log:
        temp("config/sample_prep.tsv.log"),
    threads: rule_resource("sample_prep", "threads")
    resources:
        mem_mb=rule_resource("sample_prep", "mem_mb"),
        mem_per_cpu=rule_resource("sample_prep", "mem_per_cpu"),
        partition=rule_resource("sample_prep", "partition"),
        threads=rule_resource("sample_prep", "threads"),
        time=rule_resource("sample_prep", "time"),
    container:
        rule_container("sample_prep")
    message:
        "{rule}: Create tsv with sample and column with sex (all female) to run deepmosaic"
    shell:
        """(awk 'BEGIN {{ FS="\t"; OFS=","}} {{ print $1, $2, "female"}}' {input.samples} > {output.sex}) &&
        (sed -i "1s/female/predicted_sex/" {output.sex}) && (sed -i "1s/sample/sample_id/" {output.sex})
        &> {log}"""
