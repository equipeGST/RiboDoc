configfile: 
    'config/config.yaml'

ribodoc_version = "0.9.2"


# Sets the number of threads for multi-threading steps
multi_threads_nbr = 3
mem_mb_resources = 10000
# mem_mb_resources = (workflow.cores/3)*3000
utr_threshold = "0.25"


# Sets paths for inside the container
local_path = "//data/Enora/Amibe_2024-09-30/"
ribodoc_tools = "/RiboDoc/RiboDoc/tools/"

stats_path = local_path + "stats/"
logs_path = local_path + "logs/"
snakemake_log_path = local_path + ".snakemake/log/"


# Wildcards definition
SAMPLES, = glob_wildcards(local_path + "fastq/{sample}.fastq.gz")
SAMPLES.sort()
LENGTHS = list(map(str,range(int(config['readsLength_min']),int(config['readsLength_max'])+1)))

REFERENCE_COND = str(config['reference_condition'])
OTHERS_COND = config['conditions'].split(",")

print(OTHERS_COND)
print(SAMPLES)