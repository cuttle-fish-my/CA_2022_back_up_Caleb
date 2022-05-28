import sys
import re

from pyspark import SparkContext, SparkConf


def flatMapFunc(document):
    """
    document[0] is the document ID (distinct for each document)
    document[1] is a string of all text in that document

    You will need to modify this code.
    """
    documentID = document[0]
    words = re.findall(r"\w+", document[1])

    for i in range(len(words)):
        words[i] = [words[i] + " " + documentID, i]
    return words


def mapFunc(arg):
    """
    You may need to modify this code.
    """
    return str(arg[0]), str(arg[1])


def reduceFunc(arg1, arg2):
    """
    You may need to modify this code.
    """
    return str(arg1) + " " + str(arg2)


def createIndices(file_name, output="spark-wc-out-createIndices"):
    sc = SparkContext("local[8]", "CreateIndices", conf=SparkConf().set("spark.hadoop.validateOutputSpecs", "false"))
    file = sc.sequenceFile(file_name)
    indices = file.flatMap(flatMapFunc) \
        .map(mapFunc) \
        .reduceByKey(reduceFunc) \
        .sortByKey()

    indices.coalesce(1).saveAsTextFile(output)


""" Do not worry about this """
if __name__ == "__main__":
    argv = sys.argv
    if len(argv) == 2:
        createIndices(argv[1])
    else:
        createIndices(argv[1], argv[2])
