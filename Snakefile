configfile: "config/config.yml"

include: "rules/preprocessing.smk"
include: "rules/star_alignment.smk"
include: "rules/salmon_quantification.smk"