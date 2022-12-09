
query = "PRJNA854613"
from Bio import Entrez
Entrez.email = "theo@theo.io"


import os
import re
def get_metadata(samn):
    command = f"esearch -db biosample -query {samn} | efetch > temp.txt"
    os.system(command)
    isolate = ""
    antiviral = ""
    with open("temp.txt", "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if "/isolate=" in line:
                isolate = line.split("=")[1].replace('"', "")
            if "/antiviral" in line:
                antiviral = line.split("=")[1].replace('"', "")
    return isolate, antiviral



# get all SRS ids from bio project
import os , tqdm
os.system(f"esearch -db sra -query {query} | efetch -format runinfo > temp.txt")

handle = open("temp.txt", "r")
lines = handle.readlines()
header = lines[0].split(",")
for data in tqdm.tqdm(lines[1:]):
    data = data.split(",")
    as_dict = dict(zip(header, data))
    srr = as_dict["Run"]
    samn = as_dict["BioSample"]
    isolate, antiviral = get_metadata(samn)
    print(f"{srr},{samn},{isolate},{antiviral}")

        
