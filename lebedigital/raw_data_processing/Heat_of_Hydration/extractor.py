import json as js
import os
import numpy as np
with open("calo_example.csv", "r") as f:
    data = f.readlines()
    start_time = data[2].replace(";", "").replace("\n", "").split(",")[1]
    # print(start_time)
    end_time = data[3].replace(";", "").replace("\n", "").split(",")[1]
    # print(end_time)
    mix_names = data[8].replace(";", "").replace("\n", "").split(",")[1:]
    # print(mix_names)
    sample_mass = data[9].replace(";", "").replace("\n", "").split(",")[1:]
    sample_mass = [float(m.replace("g", "")) for m in sample_mass]
    # print("time;" in data[24].lower())
    #find start line of data
    for start_line, line in enumerate(data):
        if "time;" in line.lower():
            headers = line
            break
    #find end line of data 
    for end_line, line in enumerate(data):
        if line.startswith(";;") and end_line > start_line:
            break
    headers = headers.split(";")
    numeric_data = []
    print(headers, start_line, end_line)
    for line in data[start_line+1:end_line]:
        line = line.replace("\n", "").split(";")
        line = [float(l) for l in line]
        numeric_data.append(line)
    numeric_data = np.array(numeric_data)
#write processed data to json  and csv:
n = len(mix_names)
for i, mix in enumerate(mix_names):
    mix_headers = headers[:2]
    mix_headers.extend(headers[2+(4*i):(2+4*i)+4])
    with open(os.path.join("processed_data", mix+"_calo"+".json"), "w") as f:
        js.dump({"start_time": start_time,
                 "end_time": end_time,
                 "mix_name": mix,
                 "sample_mass": sample_mass[i],
                 "lab": "KIT",
                 "headers": mix_headers
                 }, f)
    colids = [0, 1]
    colids.extend(range(2+(4*i), (2+4*i)+4))
    # print(colids)
    np.savetxt(os.path.join("processed_data", mix+"_calo"+".csv"), numeric_data[:, colids], delimiter=";", header=";".join(mix_headers).replace(",", " "))
