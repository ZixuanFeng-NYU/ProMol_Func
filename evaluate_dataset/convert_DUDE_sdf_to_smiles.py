import os

target_list=os.listdir("DUD-E/all")
for target in target_list:
    files=os.listdir("DUD-E/all/"+target)
    for file_ in files:
        if file_.endswith(".gz"):
            os.system(f"gzip -d DUD-E/all/{target}/{file_}")
        
for target in target_list:
    files=os.listdir("DUD-E/all/"+target)
    for file_ in files:
        if file_.endswith(".sdf"):
            file_name=file_.split(".")[0]
            os.system(f"obabel -isdf DUD-E/all/{target}/{file_name}.sdf -osmi -O DUD-E/all/{target}/{file_name}.smi")

