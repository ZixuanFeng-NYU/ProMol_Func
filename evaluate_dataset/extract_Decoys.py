import rarfile
with rarfile.RarFile('DEKOIS2.0_library.rar', 'r') as rar:
    rar.extractall()

file_list=rar.namelist()
for file_name in file_list:
    print(file_name)
import gzip
import shutil
import os

extracted_dir = '.'
for file_name in file_list:
    if file_name.endswith('.gz'):
        with gzip.open(os.path.join(extracted_dir, file_name), 'rb') as f_in:
            uncompressed_file_name = os.path.splitext(file_name)[0]
            with open(os.path.join(extracted_dir, uncompressed_file_name), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
print("Files have been successfully unzipped.")

