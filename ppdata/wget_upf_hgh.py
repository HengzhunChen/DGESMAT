import urllib.request
import requests
import re
import os


"""
A python script used to download HGH type UPF files from Quantum-Espresso.

Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
                        Fudan University
This file is distributed under the terms of the MIT License.
"""


#
# figure out directory to download upf files
#
dir_name = './hgh_pbe_upf/'
if os.path.exists(dir_name):
    pass
else:
    os.makedirs(dir_name, exist_ok=True)

#
# generate url for each chemical element
#
basic_url = 'http://pseudopotentials.quantum-espresso.org/'
hgh_upf_url = 'http://pseudopotentials.quantum-espresso.org/legacy_tables/hartwigesen-goedecker-hutter-pp'
hgh_table_text = requests.get(hgh_upf_url).text
regexp_table = '<a .* href="/(.*)">.*</a>'
target_url = re.findall(regexp_table, hgh_table_text)
element_urls = [basic_url + url for url in target_url]

#
# choose XC type for HGH pseudo-potential
#
xc_type = 'pbe'
# xc_type = 'blyp'
# xc_type = 'pz'

#
# download upf files
#
for element_url in element_urls:
    element_name = element_url.split('/')[-1]
    element_text = requests.get(element_url).text
    regexp_element = '<a .* href="/(.*)">.*</a>'
    upf_file_urls = re.findall(regexp_element, element_text)
    access_flag = False
    for upf_url in upf_file_urls:
        # print(upf_url)
        if re.search(xc_type, upf_url):
            url = basic_url + upf_url
            access_flag = True
            full_name = url.split('//')[-1]
            file_name = full_name.split('/')[-1]
            urllib.request.urlretrieve(url, dir_name + file_name)
            print(f"download {file_name}")
            break
    if not access_flag:
        print(f"{xc_type} type UPF file is not found for chemical element {element_name}")
