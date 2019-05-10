#!/usr/local/Cluster-Apps/python/2.7.10/bin/python

########################################################################
# Written by Olga Shamardina
#
# Script to find the number of high-quality het SNVs within
# samples' deletions boundaries
#
########################################################################

# Run:
# $ for chr in $(seq 1 22) X; do spark-submit --master yarn --deploy-mode client --num-executors 100 --executor-cores 4 --executor-memory 6G --driver-memory 6G --packages com.databricks:spark-avro_2.11:3.2.0 hadoop_query.py $chr; done

from pyspark.sql.functions import col, create_map, lit, udf, collect_set, broadcast, explode
from pyspark.sql.types import StringType, IntegerType, StructType, StructField, ArrayType, BooleanType
from pyspark.sql import SparkSession
from pyspark import StorageLevel
from itertools import chain
import json
import csv
import sys
import re
import os.path

def het_in_del(samples, het):
    rv = []
    for s in samples:
        if s in het:
            rv.append(s)
    return rv

def prefilter(samples, het):
    for s in samples:
        if s in het:
            return True
    return False

udf_het_in_del = udf(lambda samples, het: het_in_del(samples, het), ArrayType(IntegerType()))
udf_prefilter = udf(lambda samples, het: prefilter(samples, het), BooleanType())

conf = {}
with open("config.sh", "r") as f:
    for l in f:
        if l.startswith("#"):
            continue
        if "=" in l:
            line = l.rstrip().split("=")
            k = line[0].strip()
            v = re.sub("#.+$", "", line[1].strip()).strip().replace("\"", "")
            conf[k] = v

dels_file = "file://" + os.path.join(conf["HADOOPOUTPUT"], "input/dels.filt1.chr{}.samp.bed")
manifest = conf["MANIFEST"]
avro_file = conf["AVRO"]
studyConf_file = conf["STUDYCONF"]
output = conf["HADOOPOUTPUT"]
log = os.path.join(conf["HADOOPOUTPUT"], "log")

our_to_illumina = {}  # our ID to Illumina ID
with open(manifest, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for line in reader:
        our_to_illumina[line["BRIDGE_ID"]] = line["ILMN_ID"]

spark = SparkSession.builder.appName("Number of SNPs per deletion").getOrCreate()
sc = spark.sparkContext

study = json.loads(sc.textFile(studyConf_file).first())["sampleIds"]
hadoop_to_our = {study[our_to_illumina[s]]:s for s in our_to_illumina.keys() if study.has_key(our_to_illumina[s])}
our_to_hadoop = {v:k for k, v in hadoop_to_our.items()}

OPR = 0.99
chromosome = sys.argv[1]

df = spark.read.format("com.databricks.spark.avro").load(avro_file)
df = df.filter("chromosome = '{}' and minOPR >= {} and type = 'SNV'".format(chromosome, OPR)).select(col("start").alias("pos"), "genotypes.het")

with open(log, "w") as f:
    f.write("df read\n")

dels_df = spark.read.csv(dels_file.format(chromosome), sep="\t", schema = StructType([StructField("chromosome", StringType(), False), StructField("start", StringType(), False), StructField("end", StringType(), False), StructField("wgsid", StringType(), False)]))
dels_df = dels_df.withColumn("start", dels_df.start.cast(IntegerType()))
dels_df = dels_df.withColumn("end", dels_df.end.cast(IntegerType()))
mapping_expr = create_map([lit(x) for x in chain(*our_to_hadoop.items())])
dels_df = dels_df.withColumn("hid", mapping_expr.getItem(col("wgsid"))).select("start", "end", "hid")
dels_df = dels_df.groupBy("start", "end").agg(collect_set("hid").alias("samples"))

with open(log, "a") as f:
    f.write("dels_df read\n")

df = df.repartition(5000, "pos")
dels_df = dels_df.repartition("start")
df = df.persist(StorageLevel.MEMORY_AND_DISK)
df.count()
dels_df = dels_df.cache()
dels_df.count()

with open(log, "a") as f:
    f.write("pre join finished\n")

df = df.join(broadcast(dels_df), ((df.pos >= dels_df.start) & (df.pos <= dels_df.end)), how = "inner")
df = df.repartition(5000)
df = df.persist(StorageLevel.MEMORY_AND_DISK)
df.count()

with open(log, "a") as f:
    f.write("join finished\n")

df = df.filter(udf_prefilter("samples", "het")).withColumn("snp_in_del", udf_het_in_del("samples", "het"))
df = df.select("start", "end", explode("snp_in_del").alias("hid")).groupBy("hid", "start", "end").count()

headers = ["sample", "start", "end", "count"]
with open("{}/{}.csv".format(output, chromosome), "w") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(headers)
    for r in df.collect():
        writer.writerow([hadoop_to_our[r["hid"]], r["start"], r["end"], r["count"]])

print "done"
