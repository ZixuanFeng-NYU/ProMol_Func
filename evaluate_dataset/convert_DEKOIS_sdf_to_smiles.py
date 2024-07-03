import os

target_list=os.listdir("DEKOIS2.0_library/decoys/")
for target in target_list:
    if target.endswith(".sdf"):
        os.system(f"obabel -isdf DEKOIS2.0_library/decoys/{target} -osmi -O DEKOIS2.0_library/decoys/{target}.smi")

