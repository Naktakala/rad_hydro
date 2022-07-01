import sys
import numpy as np
import matplotlib.pyplot as plt

a_const = 0.013722354852

num_args = len(sys.argv)
if num_args == 1:
    sys.exit(0)

fields_to_plot = []
for arg in sys.argv[1:]:
    fields_to_plot.append(arg)

print("Fields of choice: ", fields_to_plot)

ofile = open("ZRawOutput.txt", "r")

lines = ofile.readlines()

# Process keys
field_keys = {}
num_fields = 0
field_length = len(lines)-1
for line in lines[0:1]:
    words = line.split()
    field_num = 0
    for word in words:
        field_keys[word] = field_num
        field_num += 1
    num_fields = field_num

# Process fields
fields = np.zeros([field_length, num_fields])
i = 0
for line in lines[1:len(lines)]:
    words = line.split()
    j = 0
    for word in words:
        fields[i, j] = float(word)
        j += 1
    i += 1

ofile.close()



plt.figure()
x = fields[:, field_keys["z"]]
for fieldname in fields_to_plot:
    if fieldname == "radT":
        y = fields[:, field_keys["radE"]]
        y = (y/a_const)**0.25
    else:
        y = fields[:, field_keys[fieldname]]
    plt.plot(x, y, linewidth=1, label=fieldname)

plt.legend()
plt.xlim([-0.02+0.25, 0.02+0.25])
plt.savefig("ZRawOutput.png")

print("Done")
